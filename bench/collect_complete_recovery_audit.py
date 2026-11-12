#!/usr/bin/env python3
"""Combine complete-placement searches without turning work limits into negatives.

Reference-only diagnostic. Never imported by the MCAAT prediction executable.
"""
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path, rows):
    if not rows:
        path.write_text('')
        return
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def known(sequence):
    return next((i for i, base in enumerate(sequence) if base not in 'ACGT'), len(sequence))


def roots(rows):
    result, cursors = defaultdict(dict), Counter()
    for row in rows:
        key = row['seed_id'], row['side']
        node, prefix = int(row['vertex']), row['flank_prefix']
        root = (2**64 - 1 - cursors[key] - len(prefix)) if prefix else node
        cursors[key] += len(prefix)
        result[key][root] = row
    return result


def verify_alignment(path, reference):
    q = path['sequence_outward']
    r = reference[:int(path['checkpoint_ref_bp'])]
    i = j = matches = 0
    for op in path['operations']:
        assert op in 'MXID'
        if op in 'MX':
            assert (q[i] == r[j]) == (op == 'M')
            matches += op == 'M'
            i += 1
            j += 1
        elif op == 'I':
            i += 1
        else:
            j += 1
    assert i == len(q) and j == len(r)
    assert len(r) >= 100 and matches * 10 >= 9 * len(path['operations'])
    assert len(path['vertices'].split(',')) == len(q) + 1
    return 100 * matches / len(path['operations'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--allow-incomplete', action='store_true')
    args = parser.parse_args()
    base, out = args.root, args.output
    out.mkdir(parents=True, exist_ok=True)
    inputs = {r['seed_id']: r['sequence'] for r in table(base / 'complete-placement-final/targets.tsv')}
    cases = {r['case_id']: r for r in table(base / 'adapter-free/reference_cases.tsv')}
    placements = table(base / 'complete-placement-final/complete_placements.tsv')
    manifest = json.loads((base / 'complete-placement-final/manifest.json').read_text())
    assert manifest['inputs'] == len(inputs)
    assert manifest['completed_input_sides'] == 2 * len(inputs)
    assert manifest['work_limit'] is None and manifest['candidate_limit'] is None
    all_endpoints = defaultdict(set)
    for p in placements:
        query, observed = inputs[p['seed_id']], p['observed_target']
        begin, end, errors = int(p['query_begin']), int(p['query_end']), int(p['mismatches'])
        assert len(query) == len(observed) and 0 <= begin < end <= len(query)
        assert errors == sum(a != b for a, b in zip(query[begin:end], observed[begin:end]))
        assert 100 * (end - begin) >= 85 * len(query) and errors * 10 <= end - begin
        all_endpoints[p['seed_id'], p['side']].add((p['vertex'], p['flank_prefix']))

    certificates, witnesses, exact, stage_manifest = {}, [], defaultdict(list), []
    for stage in sorted(base.glob('complete-roots-*')):
        route_file = stage / 'oracle_routes.tsv'
        if not route_file.exists():
            continue
        for row in table(stage / 'targets.tsv'):
            assert inputs[row['seed_id']] == row['sequence'], 'stage changed an input'
        for row in table(stage / 'reference_cases.tsv'):
            assert cases[row['case_id']] == row, 'stage used a different reference window'
        supplied = table(stage / 'complete_placements.tsv')
        stage_endpoints = defaultdict(set)
        for p in supplied:
            stage_endpoints[p['seed_id'], p['side']].add((p['vertex'], p['flank_prefix']))
        lookup = roots(supplied)
        relaxed = 'relaxed' in stage.name
        stage_manifest.append(dict(stage=stage.name, routes_sha256=sha(route_file),
                                   roots_sha256=sha(stage / 'complete_placements.tsv'),
                                   source_cycle_guard='relaxed' if relaxed else 'enabled'))
        paths = {(r['case_id'], r['side']): r for r in table(stage / 'oracle_paths.tsv')}
        for row in table(stage / 'exact_native_routes.tsv'):
            exact[row['seed_id']].append(row)
        for row in table(route_file):
            key = row['seed_id'], row['side']
            pair = row['case_id'], row['side']
            if row['status'] == 'verified_graph_path' and not relaxed:
                path = paths[pair]
                identity = verify_alignment(path, cases[row['case_id']][row['side'] + '_outward'])
                root_row = lookup[key][int(path['vertices'].split(',')[0])]
                assert (root_row['vertex'], root_row['flank_prefix']) in all_endpoints[key]
                witnesses.append(dict(seed_id=row['seed_id'], case_id=row['case_id'], side=row['side'],
                    reference=row['reference'], reference_bp=int(path['checkpoint_ref_bp']),
                    graph_bp=int(path['query_bp']), identity=identity,
                    observed_target=root_row['observed_target'], native_vertex=root_row['vertex'],
                    flank_prefix=root_row['flank_prefix'], stage=stage.name,
                    sequence_outward=path['sequence_outward'], operations=path['operations'],
                    vertices=path['vertices'], counts=path['counts']))
            if (row['status'] == 'no_path_at_90pct_anchored_identity'
                    and row['limited'] == '0' and row['exhausted'] == '1'
                    and stage_endpoints[key] == all_endpoints[key]):
                certificates[pair] = dict(**row, stage=stage.name, source_cycles_allowed=relaxed)

    by_input = defaultdict(list)
    for w in witnesses:
        by_input[w['seed_id']].append(w)
    result, unresolved = [], []
    for sid, sequence in inputs.items():
        locations = [r for r in cases.values() if r['seed_id'] == sid]
        required = [(r['case_id'], side) for r in locations for side in ('left', 'right')
                    if known(r[side + '_outward']) >= 100]
        missing = [pair for pair in required if pair not in certificates]
        proof = max(by_input[sid], key=lambda r: (r['reference_bp'], r['identity']), default=None)
        status = 'yes' if proof else 'no' if required and not missing else 'unresolved'
        if status == 'unresolved':
            unresolved.extend(dict(seed_id=sid, case_id=c, side=s) for c, s in missing)
        result.append(dict(seed_id=sid, sequence=sequence, recoverable_100bp_90pct=status,
            reference_locations=len(locations), eligible_reference_sides=len(required),
            completed_negative_sides=sum(pair in certificates for pair in required),
            left_placement_endpoints=len(all_endpoints[sid, 'left']),
            right_placement_endpoints=len(all_endpoints[sid, 'right']),
            longest_exact_walk_bp=max((int(r['exact_bp']) for r in exact[sid]), default=0),
            witness_reference=proof['reference'] if proof else '',
            witness_side=proof['side'] if proof else '',
            witness_reference_bp=proof['reference_bp'] if proof else '',
            witness_identity=proof['identity'] if proof else '',
            witness_stage=proof['stage'] if proof else ''))
    counts = dict(Counter(r['recoverable_100bp_90pct'] for r in result))
    summary = dict(inputs=len(inputs), counts=counts, reference_locations=len(cases),
        placement_endpoints=len(placements), target_identity_percent=90, target_coverage_percent=85,
        target_alignment='ungapped; all permitted substitutions and clipped ends',
        continuation='a connected anchored flank aligning 100..750 reference bases at >=90% identity',
        continuation_alignment='substitutions and indels; actual graph edges',
        source_array_closed_walks='excluded', external_revisits='allowed in existence search',
        adapters='excluded', negative_requires='all admissible endpoints and all catalogued reference locations; exhausted search',
        exact_walk_note='longest_exact_walk_bp does not enforce the source-cycle history guard',
        positive_note='existence proof; does not require exhausting alternative locations',
        scope='this existing graph and the supplied reference-location catalogue; not an assertion about unknown phages',
        reference_used_during_production_prediction=False, stages=stage_manifest)
    write(out / 'all81_recoverability.tsv', result)
    write(out / 'connected_path_witnesses.tsv', witnesses)
    write(out / 'completed_negative_searches.tsv', list(certificates.values()))
    write(out / 'unresolved.tsv', unresolved)
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(dict(inputs=len(inputs), counts=counts, unresolved_sides=len(unresolved)), indent=2))
    if not args.allow_incomplete and counts.get('unresolved', 0):
        raise SystemExit('Incomplete: no final absence claim is permitted.')


if __name__ == '__main__':
    main()
