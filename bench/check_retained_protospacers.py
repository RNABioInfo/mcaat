#!/usr/bin/env python3
"""Verify retained IDs/sequences and every emitted seed/extension association."""
import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path

def fasta(path):
    sequences = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]
            assert name not in sequences, f'duplicate FASTA ID: {name}'
            sequences[name] = ''
        elif line.strip():
            assert name is not None, f'sequence before FASTA header: {path}'
            sequences[name] += line.strip().upper()
    assert all(sequences.values()), f'empty FASTA record: {path}'
    return sequences

def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

def check(baseline, output):
    original = fasta(baseline)
    assert original, 'empty baseline'
    retained = fasta(output / 'Protospacers.fasta')
    assert retained == original, 'supplied IDs/sequences were lost, changed or added'
    assert list(retained) == list(original), 'supplied record order changed'
    params = json.loads((output / 'protospacer_parameters.json').read_text())
    assert params['mode'] == 'extend_supplied_protospacers'
    extensions = table(output / 'protospacer_extensions.tsv')
    by_seed = {r['seed_id']: r for r in extensions}
    assert len(by_seed) == len(extensions) == len(original) and set(by_seed) == set(original)
    assert all(r['retention_status'] == 'retained_input' for r in extensions)
    assert all(r['extension_status'] in {'extended', 'extended_ambiguous', 'extended_variant', 'extended_variant_ambiguous', 'not_extended'} for r in extensions)
    contigs = fasta(output / 'ContiguousPaths.fasta')
    regions = table(output / 'protospacer_regions.tsv')
    assert len({r['region_id'] for r in regions}) == len(regions) == len(contigs)
    for region in regions:
        sequence = contigs[region['region_id']]
        assert sequence == region['sequence'] and len(sequence) == int(region['length'])
        assert len(sequence) <= params['max_length']
        assert len(sequence) == len(region['vertices'].split(',')) + params['graph_k'] - 1
        assert len(region['edges'].split(',')) + 1 == len(region['vertices'].split(','))
    by_region = {r['region_id']: r for r in regions}
    links = table(output / 'protospacer_associations.tsv')
    assert len({r['seed_id'] for r in links}) == len(links), 'multiple final extensions for one input'
    extended = {r['seed_id'] for r in extensions if r['extension_status'].startswith('extended')}
    assert {r['seed_id'] for r in links} == extended
    assert {r['region_id'] for r in links} == set(contigs), 'unlinked extension exported'
    for link in links:
        seed = original[link['seed_id']]
        row = by_seed[link['seed_id']]
        sequence = contigs[link['region_id']]
        begin, end = int(link['seed_begin']), int(link['seed_end'])
        observed = sequence[begin:end]
        assert 0 <= begin < end <= len(sequence)
        oriented = reverse_complement(observed) if link['strand'] == '-' else observed
        if params.get('path_quality_filter'):
            assert oriented == seed, 'quality-filtered output changed the original target'
            assert max(begin, len(sequence) - end) >= params['min_extension_flank']
        if row['extension_status'].startswith('extended_variant'):
            assert params['schema'] >= 10 and params['target_variant_fallback']
            assert oriented == link['observed_target'] and len(oriented) == len(seed)
            qb, qe = int(link['query_begin']), int(link['query_end'])
            assert 0 <= qb < qe <= len(seed)
            matches = sum(a == b for a, b in zip(seed[qb:qe], oriented[qb:qe]))
            assert matches * 100 >= (qe - qb) * 90 and (qe - qb) * 100 >= len(seed) * 85
            assert abs(float(link['identity']) - matches / (qe - qb)) < 1e-9
            assert abs(float(link['query_coverage']) - (qe - qb) / len(seed)) < 1e-9
        else:
            assert oriented == seed, 'exact placement silently changed the input'
        assert len(sequence) > len(seed), 'seed-only placeholder counted as an extension'
        left, right = begin, len(sequence) - end
        assert left <= params['max_flank_length'] and right <= params['max_flank_length']
        if link['strand'] == '-':
            left, right = right, left
        assert (left, right) == (int(row['left_bp']), int(row['right_bp']))
        assert row['region_id'] == link['region_id'] and int(row['length']) == len(sequence)
        if 'max_node_visits' in params:
            vertices = by_region[link['region_id']]['vertices'].split(',')
            original_visits = Counter(vertices[begin:end - params['graph_k'] + 1])
            for vertex, visits in Counter(vertices).items():
                assert visits <= max(1, original_visits[vertex]) + params['max_node_visits'] - 1, 'node visit cap violated'
    for row in extensions:
        budget = params['max_query_states'] if params.get('target_variant_fallback') else min(params['max_states_per_side'] * 2, params['max_query_states'])
        assert int(row['expanded_states']) <= budget
        if params.get('allow_source_array_crossing'):
            flags = [int(row[side + '_crosses_source_array']) for side in ('left', 'right')]
            assert all(flag in (0, 1) for flag in flags) and sum(flags) <= 1, 'multiple source-array passages in one output'
        if row['extension_status'] == 'not_extended':
            assert row['region_id'] == '-' and int(row['length']) == 0
    if params.get('report_alternatives'):
        alternatives = table(output / 'protospacer_alternatives.tsv')
        alternative_sequences = fasta(output / 'ExtensionAlternatives.fasta')
        assert len(alternatives) == len(alternative_sequences)
        assert {r['candidate_id'] for r in alternatives} == set(alternative_sequences)
        placements = {r['seed_id']: r for r in table(output / 'protospacer_target_placements.tsv')}
        counts = Counter()
        selected = Counter()
        for row in alternatives:
            seed_id = row['seed_id']
            assert seed_id in original
            counts[seed_id] += 1
            assert int(row['rank']) == counts[seed_id] <= (6 if params.get('path_quality_filter') else 2) * params['beam_width']
            sequence = alternative_sequences[row['candidate_id']]
            assert sequence == row['sequence'] and len(sequence) == int(row['length']) <= params['max_length']
            begin, end = int(row['seed_begin']), int(row['seed_end'])
            assert 0 <= begin < end <= len(sequence)
            assert max(begin, len(sequence) - end) <= params['max_flank_length']
            observed = placements[seed_id]['observed_target']
            assert sequence[begin:end] == (original[seed_id] if observed == '-' else observed)
            if params.get('allow_source_array_crossing'):
                flags = [int(row[side + '_crosses_source_array']) for side in ('left', 'right')]
                assert all(flag in (0, 1) for flag in flags) and sum(flags) <= 1
            assert row['selected'] in ('0', '1')
            if row['selected'] == '1':
                selected[seed_id] += 1
                assert selected[seed_id] == 1 and seed_id in extended
                emitted = contigs[by_seed[seed_id]['region_id']]
                assert sequence in (emitted, reverse_complement(emitted)), 'selected alternative differs from primary output'
    assert params['input_records'] == params['retained_records'] == len(original)
    if params.get('path_quality_filter'):
        quality = table(output / 'protospacer_path_quality.tsv')
        chosen = [r for r in quality if r['selected'] == '1']
        assert len(chosen) == len(extended) and {r['seed_id'] for r in chosen} == extended
        assert all(r['accepted'] == '1' and r['reasons'] == 'pass' for r in chosen)
        if params.get('path_quality_version', 1) >= 3:
            assert params['path_quality_ambiguity_margin'] == 'relative_loss'
            for row in quality:
                for side in ('left', 'right'):
                    observations=int(row[side+'_observations'])
                    entering=int(row[side+'_entrance_observations'])
                    downstream=int(row[side+'_downstream_observations'])
                    pairs=int(row[side+'_transition_pairs'])
                    assert observations == entering + downstream
                    assert 0 <= pairs <= max(0, observations-1)
                    energy=float(row[side+'_entrance_energy'])+float(row[side+'_downstream_energy'])
                    roughness=float(row[side+'_local_roughness']); branch=float(row[side+'_branch_loss'])
                    assert all(math.isfinite(v) and v >= 0 for v in (energy,roughness,branch))
                    assert math.isclose(roughness,energy/pairs if pairs else 0,rel_tol=1e-9,abs_tol=1e-12)
                    if row['accepted']=='1' and observations:assert downstream>=64
        if params.get('path_quality_version', 1) >= 2:
            unresolved = fasta(output / 'UnresolvedExtensions.fasta')
            unresolved_rows = table(output / 'protospacer_unresolved.tsv')
            expected = {r['seed_id'] for r in extensions if r['reason'] == 'competing_quality_paths'}
            assert set(unresolved) == expected == {r['seed_id'] for r in unresolved_rows}
            assert not expected & extended
            for row in unresolved_rows:
                seq = unresolved[row['seed_id']]
                begin = int(row['left_bp']); end = len(seq) - int(row['right_bp'])
                assert seq[begin:end] == original[row['seed_id']]
                assert max(begin, len(seq) - end) >= params['min_extension_flank']
    assert params['extended_records'] == len(links)
    assert params['not_extended_records'] == len(original) - len(links)
    return dict(input_records=len(original), retained_unchanged=len(retained),
                extended_records=len(links), unique_regions=len(contigs),
                extension_statuses=dict(Counter(r['extension_status'] for r in extensions)),
                sequence_and_coordinate_checks='passed',
                reference_accuracy='not assessed by this checker')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    print(json.dumps(check(args.baseline, args.output), indent=2))
