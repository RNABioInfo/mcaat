#!/usr/bin/env python3
"""Verify retained IDs/sequences and every emitted seed/extension association."""
import argparse
import csv
import json
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
    assert all(r['extension_status'] in {'extended', 'extended_ambiguous', 'not_extended'} for r in extensions)
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
        assert (reverse_complement(observed) if link['strand'] == '-' else observed) == seed
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
        assert int(row['expanded_states']) <= min(params['max_states_per_side'] * 2, params['max_query_states'])
        if row['extension_status'] == 'not_extended':
            assert row['region_id'] == '-' and int(row['length']) == 0
    assert params['input_records'] == params['retained_records'] == len(original)
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
