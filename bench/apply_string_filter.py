#!/usr/bin/env python3
"""Apply the explicit motif/conservation rules from a verified MCAAT string audit.

The audit is bound to the original report files by SHA-256 hashes. This command
rejects flagged candidates, runs production curation on the remaining evidence,
and writes accepted reports plus every rejected observation to a NEW directory.
No graph rebuild, cycle search, or new alignment is needed for an existing audit.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import time

INPUT_FILES = ('CRISPR_graph_summary.tsv', 'CRISPR_graph_systems.tsv',
               'CRISPR_graph_paths.tsv', 'CRISPR_spacer_variants.tsv')


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(chunk)
    return result.hexdigest()


def verify_input(source, audit):
    fingerprints = audit['input_fingerprints']
    for name in INPUT_FILES:
        expected = fingerprints[name]
        path = source / name
        if path.stat().st_size != expected['bytes'] or digest(path) != expected['sha256']:
            raise ValueError(f'{name} differs from the audited input; regenerate the audit')


def rejection_rows(scores):
    result = []
    for row in scores:
        values = [row[name] for name in ('short_motif_rule', 'conserved_spacer_rule', 'proposed_review')]
        if any(value not in ('0', '1') for value in values):
            raise ValueError('rule flags must be 0 or 1')
        motif, conserved, union = map(int, values)
        if union != int(bool(motif or conserved)):
            raise ValueError('inconsistent rule union: ' + row['array_id'])
        reasons = []
        if motif:
            reasons.append('short_motif_across_unit')
        if conserved:
            reasons.append('conservation_across_spacer')
        if reasons:
            result.append(dict(array_id=row['array_id'], reason='+'.join(reasons)))
    return result


def canonical(sequence):
    return min(sequence, sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1])


def apply_filter(source, audit_dir, destination, binary):
    source, audit_dir, destination, binary = map(Path, (source, audit_dir, destination, binary))
    audit = json.loads((audit_dir / 'summary.json').read_text())
    verify_input(source, audit)
    scores = rows(audit_dir / 'candidate_scores.tsv')
    systems = rows(source / 'CRISPR_graph_systems.tsv')
    ids = {row['array_id'] for row in systems}
    if len(ids) != len(systems) or len(scores) != len(systems) or {r['array_id'] for r in scores} != ids:
        raise ValueError('audit candidate IDs do not exactly match the input')
    rejects = rejection_rows(scores)
    if not binary.is_file():
        raise ValueError('build replay_array_candidates first; missing binary: ' + str(binary))
    if destination.resolve() in (source.resolve(), audit_dir.resolve()):
        raise ValueError('output must differ from input and audit directories')
    if destination.exists() and any(destination.iterdir()):
        raise ValueError('output directory must be new or empty')
    metadata = {r['metric']: r['value'] for r in rows(source / 'CRISPR_graph_summary.tsv')}
    minimum, maximum = int(metadata['cycle_min_length']), int(metadata['cycle_max_length'])
    if not 1 <= minimum <= maximum:
        raise ValueError('invalid input cycle bounds')
    destination.mkdir(parents=True, exist_ok=True)
    manifest = destination / 'string_filter_rejections.tsv'
    with manifest.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['array_id', 'reason'], delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rejects)
    command = [str(binary.resolve()), str(source.resolve()), str(destination.resolve()),
               '--array-filter-mode', 'conservative', '--cycle-min-length', str(minimum),
               '--string-filter', 'off',
               '--cycle-max-length', str(maximum), '--reject-candidates', str(manifest.resolve())]
    started = time.monotonic()
    with (destination / 'string_filter_run.log').open('w') as stream:
        subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, check=True)
    verify_input(source, audit)
    output_metrics = {r['metric']: r['value'] for r in rows(destination / 'CRISPR_graph_summary.tsv')}
    decisions = rows(destination / 'CRISPR_graph_decisions.tsv')
    id_by_candidate = {str(i+1): row['array_id'] for i, row in enumerate(systems)}
    applied = {id_by_candidate[r['candidate_id']]: r['status'][len('string_filter:'):]
               for r in decisions if r['status'].startswith('string_filter:')}
    if applied != {r['array_id']: r['reason'] for r in rejects}:
        raise RuntimeError('applied decisions differ from the string-filter manifest')
    if int(output_metrics['assigned_observations']) + int(output_metrics['review_observations']) != audit['observations']:
        raise RuntimeError('accepted and alternative observation counts do not preserve the input')
    result = dict(input=str(source.resolve()), audit=str(audit_dir.resolve()),
                  audit_scores_sha256=digest(audit_dir / 'candidate_scores.tsv'),
                  audit_summary_sha256=digest(audit_dir / 'summary.json'),
                  input_fingerprints=audit['input_fingerprints'], rule_parameters=audit['parameters'],
                  original_candidates=len(systems), string_rejected_candidates=len(rejects),
                  redundant_candidates_removed=int(output_metrics['redundant_candidates_removed']),
                  retained_candidates=int(output_metrics['arrays']),
                  original_spacer_entries=audit['observations'],
                  retained_spacer_entries=int(output_metrics['spacers']),
                  retained_unique_spacer_sequences=int(output_metrics['unique_spacer_sequences_across_arrays']),
                  alternative_observations=int(output_metrics['review_observations']),
                  seconds=time.monotonic()-started,
                  interpretation='Explicit computational string criteria, not validated biological labels or full benchmark precision/recall.')
    ref_file = audit_dir / 'reference_sets.json'
    if ref_file.exists():
        known = {canonical(s) for r in json.loads(ref_file.read_text()) for s in r['spacers']}
        before = {canonical(r['variant']) for r in rows(source / 'CRISPR_spacer_variants.tsv')}
        after = {canonical(r['variant']) for r in rows(destination / 'CRISPR_spacer_variants.tsv')}
        result['local_reference_matches_before'] = len(known & before)
        result['local_reference_matches_after'] = len(known & after)
        result['lost_local_reference_spacers'] = sorted((known & before) - after)
    (destination / 'string_filter_summary.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('audit', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--replay-binary', required=True, type=Path)
    args = parser.parse_args()
    try:
        result = apply_filter(args.input, args.audit, args.output, args.replay_binary)
    except (ValueError, KeyError, OSError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, f'Error: {error}\n')
    print(json.dumps({k: v for k, v in result.items() if k not in
                     ('input_fingerprints', 'rule_parameters', 'lost_local_reference_spacers')}, indent=2))


if __name__ == '__main__':
    main()
