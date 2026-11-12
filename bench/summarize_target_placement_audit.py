#!/usr/bin/env python3
"""Join all original/reference target placements with graph-only predictions.

Existence means a demonstrated path, not a proved global maximum over all
possible target variants. Unmapped reference alleles and work limits stay
distinct from exhausted searches.
"""
import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from check_retained_protospacers import check, fasta, table


def write(path, rows):
    assert rows
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def verify_alignment(reference, row):
    sequence, operations = row['sequence_outward'], row['operations']
    ri = qi = matches = 0
    for op in operations:
        if op in 'MX':
            assert (reference[ri] == sequence[qi]) == (op == 'M')
            matches += op == 'M'
            ri += 1
            qi += 1
        elif op == 'D':
            ri += 1
        elif op == 'I':
            qi += 1
        else:
            raise AssertionError('invalid operation')
    assert ri == int(row['checkpoint_ref_bp']) and qi == len(sequence)
    assert matches * 100 >= len(operations) * 90
    assert len(row['vertices'].split(',')) == qi + 1
    assert len(row['counts'].split(',')) == qi
    return 100 * matches / len(operations)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('audit', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    root, output = args.audit, args.output
    output.mkdir(parents=True, exist_ok=True)
    targets = {r['seed_id']: r['sequence'] for r in table(root / 'adapter-free/targets.tsv')}
    assert len(targets) == 81
    cases = {r['case_id']: r for r in table(root / 'adapter-free/reference_cases.tsv')}
    assert len(cases) == 1041
    exact = table(root / 'reference-alleles-exact/exact_native_routes.tsv')
    assert len(exact) == 2082
    original = {(r['case_id'], r['side']): r for r in table(root / 'report/reference_location_certificates.tsv')}
    original_inputs = {r['seed_id']: r for r in table(root / 'report/all81_recovery.tsv')}
    routes = {(r['case_id'], r['side']): r for r in table(root / 'reference-alleles-approx/oracle_routes.tsv')}
    paths = {(r['case_id'], r['side'], r['checkpoint_ref_bp']): r for r in table(root / 'reference-alleles-approx/oracle_paths.tsv')}
    positive = {cases[r['case_id']]['seed_id'] for r in routes.values() if r['status'] == 'verified_graph_path'}
    reuse = []
    # A cap in a repeated search cannot erase a previously completed proof for
    # the identical input, reference window, side and equal/looser path policy.
    for key, row in routes.items():
        case = cases[key[0]]
        sid = case['seed_id']
        if row['status'] == 'search_limit_unresolved' and targets[sid] == case['reference_allele']:
            previous = original[key]
            if previous['status'] == 'no_path_at_90pct_anchored_identity':
                assert previous['limited'] == '0' and previous['exhausted'] == '1'
                replacement = dict(previous)
                replacement['seed_id'] = row['seed_id']
                replacement['evidence_stage'] = 'original_identical_allele:' + previous['evidence_stage']
                routes[key] = replacement
                reuse.append(dict(case_id=key[0], side=key[1], seed_id=sid, original_stage=previous['evidence_stage']))

    by_input = defaultdict(list)
    certificates = []
    for native in exact:
        key = native['case_id'], native['side']
        case = cases[key[0]]
        sid = case['seed_id']
        if native['seed_mapped'] == '0':
            status, stage = 'reference_allele_not_mapped', 'exact_reference_allele_lookup'
        elif int(native['known_reference_bp']) < 100:
            status, stage = 'insufficient_known_reference_flank', 'reference_length'
        else:
            row = routes[key]
            status = row['status']
            stage = row.get('evidence_stage', 'reference-alleles-approx')
            if sid not in positive:
                assert status == 'no_path_at_90pct_anchored_identity', (sid, key, status)
                assert row['limited'] == '0' and row['exhausted'] == '1'
        record = dict(case_id=key[0], seed_id=sid, side=key[1], reference=case['reference'],
                      reference_allele=case['reference_allele'], reference_allele_mapped=int(native['seed_mapped']),
                      known_reference_bp=int(native['known_reference_bp']), exact_flank_bp=int(native['exact_bp']),
                      status=status, evidence=stage, exact_stop_vertex=native['vertex'],
                      expected_next_base=native['expected_base'], available_edges=native['available_edges'])
        certificates.append(record)
        by_input[sid].append(record)

    witnesses = []
    for sid in sorted(positive):
        choices = [r for r in routes.values() if cases[r['case_id']]['seed_id'] == sid and r['status'] == 'verified_graph_path']
        row = max(choices, key=lambda r: (int(r['checkpoint_ref_bp']), -int(r['edit_cost'])))
        case = cases[row['case_id']]
        path = paths[(row['case_id'], row['side'], row['checkpoint_ref_bp'])]
        identity = verify_alignment(case[row['side'] + '_outward'], path)
        witnesses.append(dict(seed_id=sid, case_id=row['case_id'], reference=case['reference'],
                              reference_allele=case['reference_allele'], side=row['side'],
                              reference_bp=int(row['checkpoint_ref_bp']), graph_bp=len(path['sequence_outward']),
                              identity=identity, sequence_outward=path['sequence_outward'],
                              operations=path['operations'], vertices=path['vertices'], counts=path['counts']))

    beam = root / 'final-beams'
    comparisons = {(r['run'], r['seed_id']): r for r in table(beam / 'per_seed.tsv')}
    extension_rows = {n: {r['seed_id']: r for r in table(beam / f'N{n}/protospacer_extensions.tsv')} for n in (16, 32, 64)}
    placement_rows = {r['seed_id']: r for r in table(beam / 'N32/protospacer_target_placements.tsv')}
    validations = {}
    expected_hash = '9a7886832f257d1d7e4bc69f0cd64aab9354195dbdad47a456e5c6b1376cff12'
    for n in (16, 32, 64):
        directory = beam / f'N{n}'
        assert fasta(directory / 'Protospacers.fasta') == targets
        digest = hashlib.sha256((directory / 'Protospacers.fasta').read_bytes()).hexdigest()
        assert digest == expected_hash
        validations[f'N{n}'] = dict(sha256=digest, retained=81)
    rows = []
    for sid, query in targets.items():
        values = by_input[sid]
        assert values
        mapped = any(r['reference_allele_mapped'] for r in values)
        selected = comparisons[('N32', sid)]['match_90_50'].lower() in ('1', 'true')
        if selected:
            assert sid in positive
            status = 'selected_reference_match'
        elif sid in positive:
            status = 'matching_graph_path_exists_but_another_target_is_selected'
        elif mapped:
            status = 'mapped_reference_alleles_have_no_qualifying_flank_at_tested_locations'
        else:
            status = 'reference_alleles_absent_original_anchor_has_no_qualifying_flank'
        if sid not in positive:
            assert original_inputs[sid]['classification'] == 'no_qualifying_adapter_free_source_cycle_free_path'
        best_exact = max(values, key=lambda r: r['exact_flank_bp'])
        record = dict(seed_id=sid, input_sequence=query, input_length=len(query), outcome=status,
                      reference_locations=len(values)//2, reference_allele_mapped=mapped,
                      longest_exact_reference_allele_flank=best_exact['exact_flank_bp'],
                      exact_reference=best_exact['reference'], exact_side=best_exact['side'],
                      reference_consistent_path_demonstrated=sid in positive,
                      original_anchor_result=original_inputs[sid]['classification'],
                      N32_observed_target=placement_rows[sid]['observed_target'],
                      N32_placement_mode=placement_rows[sid]['placement_mode'],
                      N32_placement_limited=placement_rows[sid]['placement_limited'])
        for n in (16, 32, 64):
            record[f'N{n}_length'] = extension_rows[n][sid]['length']
            record[f'N{n}_match_90_identity_50_coverage'] = int(comparisons[(f'N{n}', sid)]['match_90_50'].lower() in ('1', 'true'))
            record[f'N{n}_extension_reason'] = extension_rows[n][sid]['reason']
        rows.append(record)
    write(output / 'all81_results.tsv', rows)
    write(output / 'reference_location_results.tsv', certificates)
    write(output / 'demonstrated_paths.tsv', witnesses)
    if reuse:
        write(output / 'reused_certificates.tsv', reuse)
    summary = dict(inputs=81, reference_locations=1041,
                   demonstrated_matching_graph_paths_inputs=len(positive),
                   selected_reference_matches_N32=sum(r['N32_match_90_identity_50_coverage'] for r in rows),
                   outcome_counts=dict(Counter(r['outcome'] for r in rows)),
                   inputs_with_an_exact_reference_allele_in_graph=sum(r['reference_allele_mapped'] for r in rows),
                   inputs_without_an_exact_reference_allele_in_graph=sum(not r['reference_allele_mapped'] for r in rows),
                   input_preservation=validations,
                   initial_repeat_filter_used_for_absence_certificates=False,
                   reference_used_by_production=False,
                   global_maximum_proved=False,
                   scope='All original inputs and all 1041 reference-allele locations; known flanks 100..750 bases, >=90% aligned identity, adapters and closed source-array walks excluded. Original input preservation is separate from observed target-allele placement.',
                   limitation='No upper bound is claimed over arbitrary target variants absent from the reference collection, unknown reference gaps, or every larger approximate walk. One original-anchor search for an already recovered variant input remains unresolved. Alternative locations of already positive inputs can also remain work-limited.')
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
