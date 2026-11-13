#!/usr/bin/env python3
"""Compare graph-only beam widths on immutable inputs; BLAST only after prediction."""
import argparse
import csv
import hashlib
import json
import re
import shlex
import statistics
import subprocess
import time
from collections import Counter
from pathlib import Path
from check_retained_protospacers import check, fasta, table

FIELDS = 'qseqid sseqid qlen slen qstart qend sstart send length nident pident mismatch gapopen gaps evalue bitscore score sstrand qseq sseq'.split()

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def write_table(path, rows):
    if not rows:
        path.write_text('')
        return
    with path.open('w', newline='') as f:
        fields = list(dict.fromkeys(key for row in rows for key in row))
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('binary', 'arrays', 'graph', 'protospacers', 'baseline', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--database', required=True)
    p.add_argument('--blastn', default='blastn')
    p.add_argument('--blastdbcmd', default='blastdbcmd')
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--widths', type=int, nargs='+', default=[16, 32, 64])
    p.add_argument('--max-node-visits', type=int, default=3)
    p.add_argument('--max-states', type=int, default=50000)
    p.add_argument('--report-alternatives', action='store_true')
    p.add_argument('--analyze-only', action='store_true')
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    graph_files = sorted(a.graph.parent.glob(a.graph.name + '.sdbg*'))
    assert graph_files, 'graph parts missing'
    graph_state = {str(f): [f.stat().st_size, f.stat().st_mtime_ns] for f in graph_files}
    meta = dict(arguments={k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()},
                binary_sha256=sha(a.binary), input_sha256=sha(a.protospacers),
                arrays_sha256=sha(a.arrays / 'CRISPR_Arrays_1.txt'),
                graph_metadata_sha256=sha(Path(str(a.graph) + '.sdbg_info')),
                graph_files=graph_state, reference_used_during_search=False,
                identity_percent=90, coverage_percent=50,
                coverage_definition='single HSP query-coordinate span / full predicted extension length',
                baseline_parameters=json.loads((a.baseline / 'protospacer_parameters.json').read_text()))
    if not a.analyze_only:
        assert not (a.output / 'manifest.json').exists(), 'Use a fresh directory or --analyze-only'
        (a.output / 'manifest.json').write_text(json.dumps(meta, indent=2) + '\n')
        for n in a.widths:
            out = a.output / f'N{n}'
            cmd = [str(a.binary), 'phage-curate', '--arrays', str(a.arrays), '--graph', str(a.graph),
                   '--protospacers', str(a.protospacers), '--output', str(out), '--threads', str(a.threads),
                   '--beam-width', str(n), '--max-node-visits', str(a.max_node_visits),
                   '--max-length', '1550', '--max-flank-length', '750', '--max-states', str(a.max_states)]
            if a.report_alternatives:cmd += ['--report-alternatives','true']
            (a.output / f'N{n}.command.txt').write_text(shlex.join(cmd) + '\n')
            print(f'Starting N={n}', flush=True)
            start = time.monotonic()
            with (a.output / f'N{n}.log').open('w') as log:
                subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
            elapsed = time.monotonic() - start
            validation = check(a.protospacers, out)
            validation.update(wall_seconds=elapsed)
            (a.output / f'N{n}.validation.json').write_text(json.dumps(validation, indent=2) + '\n')
            print(f'N={n}: {validation}', flush=True)
    assert graph_state == {str(f): [f.stat().st_size, f.stat().st_mtime_ns] for f in graph_files}, 'graph modified'
    old_meta = json.loads((a.output / 'manifest.json').read_text())
    assert old_meta['binary_sha256'] == sha(a.binary), 'binary changed since benchmark'
    assert old_meta['input_sha256'] == sha(a.protospacers), 'input changed since benchmark'
    info = subprocess.check_output([a.blastdbcmd, '-db', a.database, '-info'], text=True)
    subjects = int(re.search(r'([\d,]+) sequences', info).group(1).replace(',', ''))
    (a.output / 'reference_info.txt').write_text(info)
    (a.output / 'blast_version.txt').write_text(subprocess.check_output([a.blastn, '-version'], text=True))
    # One shared BLAST input for identical baseline/new regions. Preserve all run
    # associations; matching regions and matching input seeds are separate counts.
    sequences = {}
    runs = [('baseline', a.baseline)] + [(f'N{n}', a.output / f'N{n}') for n in a.widths]
    for _, directory in runs:
        check(a.protospacers, directory)
        for rid, seq in fasta(directory / 'ContiguousPaths.fasta').items():
            assert rid not in sequences or sequences[rid] == seq
            sequences[rid] = seq
    hit_files = []
    for label, word, low, high in [('short', 7, 0, 100), ('long', 11, 100, 100000)]:
        query = a.output / (label + '.fasta')
        query.write_text(''.join(f'>{rid}\n{seq}\n' for rid, seq in sequences.items() if low < len(seq) <= high))
        if not query.stat().st_size:
            continue
        raw = a.output / (label + '.blast.tsv')
        cmd = [a.blastn, '-task', 'blastn', '-query', str(query), '-db', a.database, '-strand', 'both',
               '-word_size', str(word), '-dust', 'no', '-soft_masking', 'false', '-evalue', '1000',
               '-perc_identity', '90', '-qcov_hsp_perc', '50', '-max_target_seqs', str(subjects),
               '-num_threads', str(a.threads), '-outfmt', '6 ' + ' '.join(FIELDS), '-out', str(raw)]
        (a.output / (label + '.blast.command.txt')).write_text(shlex.join(cmd) + '\n')
        print('BLAST evaluation:', label, flush=True)
        with (a.output / (label + '.blast.log')).open('w') as log:
            subprocess.run(cmd, check=True, stderr=log)
        hit_files.append(raw)
    best = {}
    for f in hit_files:
        with f.open() as handle:
            for h in csv.DictReader(handle, fieldnames=FIELDS, delimiter='\t'):
                span = abs(int(h['qend']) - int(h['qstart'])) + 1
                if int(h['nident']) * 100 < 90 * int(h['length']) or span * 100 < 50 * int(h['qlen']):
                    continue
                h['query_coverage_percent'] = span * 100 / int(h['qlen'])
                if h['qseqid'] not in best or float(h['bitscore']) > float(best[h['qseqid']]['bitscore']):
                    best[h['qseqid']] = h
    write_table(a.output / 'matching_alignments.tsv', list(best.values()))
    summaries, comparisons = [], []
    baseline_matches = set()
    for name, directory in runs:
        rows = table(directory / 'protospacer_extensions.tsv')
        params = json.loads((directory / 'protospacer_parameters.json').read_text())
        matching = {r['seed_id'] for r in rows if r['region_id'] in best}
        if name == 'baseline':
            baseline_matches = matching
        lengths = [int(r['length']) for r in rows if r['region_id'] != '-']
        visits = []
        for r in table(directory / 'protospacer_regions.tsv'):
            visits.append(max(Counter(r['vertices'].split(',')).values()))
        summary = dict(run=name, retained=len(rows), extended_records=len(lengths),
                       unique_regions=params['regions'], matching_input_records=len(matching),
                       matching_unique_regions=len({r['region_id'] for r in rows if r['region_id'] in best}),
                       baseline_matches_kept=len(matching & baseline_matches),
                       baseline_matches_lost=len(baseline_matches - matching), new_matches=len(matching - baseline_matches),
                       median_extension_length=statistics.median(lengths) if lengths else 0,
                       min_extension_length=min(lengths, default=0), max_extension_length=max(lengths, default=0),
                       length_at_least_500=sum(n >= 500 for n in lengths), ambiguous=params['ambiguous_records'],
                       search_limited=params['search_limited_records'], expanded_states=sum(int(r['expanded_states']) for r in rows),
                       work_limited_records=sum('work_limit' in r['observed_stops'] for r in rows),
                       regions_with_revisited_nodes=sum(n > 1 for n in visits),
                       wall_seconds=json.loads((a.output / (name + '.validation.json')).read_text())['wall_seconds'] if name != 'baseline' else 'not re-timed')
        summaries.append(summary)
        for r in rows:
            h = best.get(r['region_id'])
            comparisons.append(dict(run=name, **r, match_90_50=bool(h), reference=h['sseqid'] if h else '',
                                    identity=h['pident'] if h else '', coverage=h['query_coverage_percent'] if h else ''))
        print(summary, flush=True)
    write_table(a.output / 'summary.tsv', summaries)
    write_table(a.output / 'per_seed.tsv', comparisons)
    (a.output / 'summary.json').write_text(json.dumps(summaries, indent=2) + '\n')
    (a.output / 'integrity.json').write_text(json.dumps(dict(all_input_records_and_coordinates_verified=True,
        graph_sizes_and_mtimes_unchanged=True, reference_used_during_search=False), indent=2) + '\n')

if __name__ == '__main__':
    main()
