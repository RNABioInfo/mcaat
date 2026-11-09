#!/usr/bin/env python3
"""Real-SDBG regression for adapter-mediated target extension.

The first spacer/adapter boundary is copied from the failing SRR14017553
prediction; all surrounding DNA and the other controls are synthetic. Builds
one graph from 150-base reads, runs both guarded and explicitly unguarded walks,
and verifies coordinates, one prediction per spacer, and graph immutability.
An optional pre-fix binary provides an independent before/after comparison.
"""
import argparse
import csv
import hashlib
import json
import random
import subprocess
import time
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('binary', type=Path)
p.add_argument('output', type=Path)
p.add_argument('--before-binary', type=Path)
a = p.parse_args()
binary, out = a.binary.resolve(), a.output.resolve()
out.mkdir(parents=True, exist_ok=True)
if any(out.iterdir()):
    raise SystemExit('Use an empty output directory')
rng = random.Random(902814)
def dna(n): return ''.join(rng.choice('ACGT') for _ in range(n))
def rc(s): return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
def table(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))
def run(name, arguments, exe=binary, success=True):
    started = time.monotonic()
    with (out / (name + '.log')).open('w') as log:
        r = subprocess.run([str(exe), *map(str, arguments)], stdout=log,
                           stderr=subprocess.STDOUT, timeout=180)
    assert (r.returncode == 0) == success, f'{name}: exit {r.returncode}; see log'
    return time.monotonic() - started
def hashes(prefix):
    return {f.name: hashlib.sha256(f.read_bytes()).hexdigest()
            for f in prefix.parent.glob(prefix.name + '.*')}

adapter = 'AGATCGGAAGAGC'
motifs = [adapter, rc(adapter), 'CTGTCTCTTATACACATCT', rc('CTGTCTCTTATACACATCT')]
repeat = dna(36)
s1 = 'AAACCATCGTTAGTTGTGATGCTACGAAT'
tail = 'AGATCGGAAGAGCACAAGTCTGAACACCAGTCACAGATCAGTATCTCGGATGCCGTCTTA'
s2, s3, s4 = dna(35), dna(34), dna(33)
bad = dna(24) + 'AGATCGGAAGAGCACACGT' + dna(25)
truth = [dna(1500) + s1, s2 + dna(1500), dna(800) + s3 + dna(800), dna(1500) + s4]
targets = [truth[0] + tail + dna(1500),
           dna(1500) + rc(tail) + truth[1],
           truth[2],
           truth[2][:800 + len(s3) + 250] + tail + dna(1500),
           truth[3] + 'CTGTCTCTTATACACATCT' + dna(1500),
           dna(1000) + bad + dna(1000)]
spacers = [s1, s2, s3, s4, bad]
hosts = [dna(200) + repeat + s + repeat + dna(200) for s in spacers]
incomplete = dna(35)
# Only five source-repeat-compatible bases exist on the left. An adapter on
# the right must not turn that unresolved left context into a completed negative.
targets.append(repeat[-5:] + incomplete + tail + dna(1500))
spacers.append(incomplete)
reads = out / 'reads.fasta'
with reads.open('w') as f:
    for index, sequence in enumerate(targets + hosts):
        for start in range(len(sequence) - 149):
            for copy in range(4):
                f.write(f'>{index}_{start}_{copy}\n{sequence[start:start+150]}\n')
arrays = out / 'arrays'
arrays.mkdir()
report = ''.join(f'>Array_{i}\n{repeat}\n        ' + len(repeat) * '-' + '\t' + s + '\n'
                 for i, s in enumerate(spacers))
report += f'>Array_RC\n{rc(repeat)}\n        ' + len(repeat) * '-' + '\t' + rc(s1) + '\n'
(arrays / 'CRISPR_Arrays_1.txt').write_text(report)
(out / 'truth.json').write_text(json.dumps(dict(clean_sequences=truth, targets=targets, spacers=spacers), indent=2) + '\n')
graph_run = out / 'graph_build'
run('build', ['--input-files', reads, '--output-folder', graph_run, '--threads', '2',
              '--ram', '2G', '--threshold-multiplicity', '1000000', '--autoclean', 'false'])
graph = graph_run / 'graph/graph'
original_hashes = hashes(graph)
common = ['phage-curate', '--arrays', arrays, '--graph', graph, '--max-length', '1200']
results = {}
cases = [('guarded', [], binary), ('unguarded', ['--adapter-sequences', 'none'], binary)]
if a.before_binary:
    cases.append(('before', [], a.before_binary.resolve()))
for name, flags, exe in cases:
    result = out / name
    elapsed = run(name, [*common, '--output', result, *flags], exe)
    regions = table(result / 'protospacer_regions.tsv')
    links = table(result / 'protospacer_associations.tsv')
    seeds = table(result / 'protospacer_candidates.tsv')
    queries = table(result / 'protospacer_queries.tsv')
    params = json.loads((result / 'protospacer_parameters.json').read_text())
    contaminated = sum(any(m in r['sequence'] for m in motifs) for r in regions)
    if name == 'guarded':
        assert contaminated == 0, 'adapter sequence retained'
        assert len(seeds) == len({s['query_id'] for s in seeds}) == 4, 'clean spacer lost or multiple final predictions'
        assert any(q['spacer'] == bad and q['mapping_status'] == 'adapter_in_spacer' for q in queries)
        assert any(q['spacer'] == incomplete and q['mapping_status'] == 'unresolved' for q in queries), 'opposite incomplete flank incorrectly treated as a completed negative'
        assert all(any(r['sequence'] in t or rc(r['sequence']) in t for t in truth) for r in regions), 'incorrect flank retained'
        assert all(int(r['length']) >= 250 for r in regions), 'clean flank erased with contaminating side'
        by_region = {r['region_id']: r for r in regions}
        by_source = {q['source_id']: q for q in queries}
        for link in links:
            seq = by_region[link['region_id']]['sequence']
            expected = by_source[link['source_id']]['spacer'][int(link['query_begin']):int(link['query_end'])]
            actual = seq[int(link['seed_begin']):int(link['seed_end'])]
            assert (rc(actual) if link['strand'] == '-' else actual) == expected, 'seed coordinate changed during adapter rollback'
        for r in regions:
            assert len(r['sequence']) == len(r['vertices'].split(',')) + 22
            assert len(r['edges'].split(',')) + 1 == len(r['vertices'].split(','))
        stops = table(result / 'protospacer_adapter_stops.tsv')
        assert stops and params['adapter_blocked_branches'] >= len(stops)
        assert {'left', 'right'} <= {s['side'] for s in stops}
        assert any(s['motif'] in ('CTGTCTCTTATACACATCT', rc('CTGTCTCTTATACACATCT')) for s in stops)
        assert params['adapter_contaminated_queries'] == 1
        parallel = out / 'parallel'
        run('parallel', [*common, '--output', parallel, '--threads', '2'])
        for f in result.iterdir():
            if f.name != 'protospacer_parameters.json':
                assert f.read_bytes() == (parallel / f.name).read_bytes(), f'thread-dependent {f.name}'
    else:
        assert contaminated > 0, 'fixture did not reproduce contaminated extensions'
    results[name] = dict(regions=len(regions), contaminated_regions=contaminated, seconds=elapsed)
assert hashes(graph) == original_hashes, 'extension modified the graph'
run('invalid_adapter', [*common, '--output', out / 'invalid', '--adapter-sequences', 'ACGT'], success=False)
(out / 'validation.json').write_text(json.dumps(results, indent=2) + '\n')
print(json.dumps(results, indent=2))
