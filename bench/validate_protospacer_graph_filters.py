#!/usr/bin/env python3
"""Graph-only acceptance controls; no read evidence is supplied to phage-curate.

Builds a real SDBG from controlled sequences. Checks informative flanks,
incompatible abundance contexts and equally scored extension alternatives.
The construction reads are used ONLY to build the graph, not for validation.
"""
import argparse
import csv
import hashlib
import json
import random
import subprocess
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('binary', type=Path)
p.add_argument('output', type=Path)
a = p.parse_args()
binary, out = a.binary.resolve(), a.output.resolve()
out.mkdir(parents=True, exist_ok=True)
if any(out.iterdir()):
    raise SystemExit('Use an empty output directory')
rng = random.Random(110028)
def dna(n): return ''.join(rng.choice('ACGT') for _ in range(n))
def rc(s): return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
def table(p):
    with p.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))
def run(name, arguments):
    with (out / (name + '.log')).open('w') as log:
        subprocess.run([str(binary), *map(str, arguments)], stdout=log,
                       stderr=subprocess.STDOUT, check=True, timeout=180)

repeat = 'CGGTTCATCCCCGCTTGCGCGGGGAACTC'
spacers = {name: dna(35) for name in ('true', 'short', 'depth_conflict', 'branch_tie', 'seed_tie')}
truth = dna(650) + spacers['true'] + dna(650)
short = 'AAAAAAAAAA' + spacers['short'] + 'TTTTTTTTTT'
depth_left, depth_right = dna(650), dna(650)
tie_left, tie_common = dna(650), dna(100)
tie_a, tie_b = 'A' + dna(549), 'C' + dna(549)
unlinked = [(dna(649) + 'A', 'A' + dna(649)), (dna(649) + 'C', 'C' + dna(649))]
genomes = [(truth, 8), (short, 120),
           (depth_left + spacers['depth_conflict'], 4),
           (spacers['depth_conflict'] + depth_right, 200),
           (tie_left + spacers['branch_tie'] + tie_common + tie_a, 4),
           (tie_left + spacers['branch_tie'] + tie_common + tie_b, 4)]
genomes += [(l + spacers['seed_tie'] + r, 8) for l, r in unlinked]
genomes += [(dna(250) + repeat + s + repeat + dna(250), 12) for s in spacers.values()]
reads = out / 'construction_reads.fasta'
with reads.open('w') as f:
    for index, (sequence, copies) in enumerate(genomes):
        length = min(150, len(sequence))
        for start in range(len(sequence) - length + 1):
            for copy in range(copies):
                f.write(f'>{index}_{start}_{copy}\n{sequence[start:start+length]}\n')
arrays = out / 'arrays'
arrays.mkdir()
(arrays / 'CRISPR_Arrays_1.txt').write_text(''.join(
    f'>Array_{name}\n{repeat}\n        ' + '-' * len(repeat) + '\t' + s + '\n'
    for name, s in spacers.items()))
graph_run = out / 'graph_build'
run('build', ['--input-files', reads, '--output-folder', graph_run, '--threads', '2',
              '--ram', '2G', '--threshold-multiplicity', '1000000', '--autoclean', 'false'])
graph = graph_run / 'graph/graph'
graph_hash = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in graph.parent.glob('graph.*')}
common = ['phage-curate', '--arrays', arrays, '--graph', graph, '--max-length', '1000']
results = {}
legacy_walk = ['--branch-guard', 'false', '--junction-check', 'false', '--prediction-score-margin', '0', '--max-flank-length', '2500']
for name, flags in [('default', []), ('exact_prediction_ties', ['--prediction-score-margin', '0']), ('post_selection_only', legacy_walk),
                    ('legacy_gates', [*legacy_walk, '--min-external-flank', '1', '--joint-flank-check', 'false', '--resolve-graph-ties', 'false'])]:
    result = out / name
    run(name, [*common, '--output', result, *flags])
    regions = {r['region_id']: r['sequence'] for r in table(result / 'protospacer_regions.tsv')}
    links = table(result / 'protospacer_associations.tsv')
    queries = {q['array_id']: q for q in table(result / 'protospacer_queries.tsv')}
    parameters = json.loads((result / 'protospacer_parameters.json').read_text())
    assert parameters['schema'] == 5 and parameters['reads1'] == ''
    for row in table(result / 'protospacer_regions.tsv'):
        assert len(row['sequence']) == int(row['length']) == len(row['vertices'].split(',')) + 22
        assert len(row['edges'].split(',')) + 1 == len(row['vertices'].split(','))
    for link in links:
        expected = queries[link['array_id']]['spacer'][int(link['query_begin']):int(link['query_end'])]
        observed = regions[link['region_id']][int(link['seed_begin']):int(link['seed_end'])]
        assert (rc(observed) if link['strand'] == '-' else observed) == expected, 'trimmed seed coordinates disagree'
    def predicted(case): return [regions[x['region_id']] for x in links if x['array_id'] == 'Array_' + case]
    def inside(sequence, target): return sequence in target or rc(sequence) in target
    assert all(len(predicted(case)) <= 1 for case in spacers), 'multiple predictions per input spacer'
    assert predicted('true') and all(inside(s, truth) for s in predicted('true')), 'unambiguous true target lost'
    assert all(len(s) >= 500 for s in predicted('true')), 'true target unnecessarily reduced'
    if name != 'legacy_gates':
        assert not predicted('short') and queries['Array_short']['mapping_status'] == 'unresolved', 'short tails treated as independent extension'
        assert all(inside(s, depth_left + spacers['depth_conflict']) or inside(s, spacers['depth_conflict'] + depth_right)
                   for s in predicted('depth_conflict')), 'incompatible-depth halves still joined'
        decisions = table(result / 'protospacer_decisions.tsv')
        assert any(d['reason'] == 'flank_multiplicity_conflict' for d in decisions)
        if name == 'default':
            assert not predicted('depth_conflict'), 'near-best incompatible choices forced into a prediction'
        if name == 'exact_prediction_ties':
            assert predicted('depth_conflict'), 'zero margin did not restore the slightly better one-sided candidate'
            qid=queries['Array_depth_conflict']['query_id']
            choice=next(s for s in table(result / 'protospacer_selection.tsv') if s['query_id']==qid)
            assert choice['score_tie']=='0' and 0<float(choice['score_gap'])<.1, 'fixture must exercise a nonzero near-score gap'
        assert predicted('branch_tie') and all(inside(s, tie_left + spacers['branch_tie'] + tie_common) for s in predicted('branch_tie')), 'arbitrary tied branch retained'
        assert not predicted('seed_tie') and queries['Array_seed_tie']['mapping_status'] == 'unresolved', 'indistinguishable seed flank combinations forced into a locus'
        if name == 'post_selection_only':
            assert any(s['ambiguity_trimmed'] == '1' for s in table(result / 'protospacer_selection.tsv'))
        else:
            assert any(s['decision'] == 'ambiguous_branch' for s in table(result / 'protospacer_branch_checks.tsv'))
        run('threads_'+name, [*common, '--output', out / ('threads_'+name), '--threads', '2', *flags])
        for f in result.iterdir():
            if f.name != 'protospacer_parameters.json':
                assert f.read_bytes() == (out / ('threads_'+name) / f.name).read_bytes(), f'thread-dependent output: {f.name}'
    else:
        assert predicted('short'), 'short-tail control did not reproduce previous acceptance'
        # The old stability ranking can already prefer a one-sided path here.
        # This case exercises the explicit conflict gate above; it does not
        # establish that the old final prediction was a crossed-depth join.
        assert not any(d['reason'] == 'flank_multiplicity_conflict'
                       for d in table(result / 'protospacer_decisions.tsv'))
        assert any(not inside(s, tie_left + spacers['branch_tie'] + tie_common) for s in predicted('branch_tie')), 'equal-score arbitrary extension was not reproduced'
        assert predicted('seed_tie'), 'forced seed-flank pairing was not reproduced'
    results[name] = {case: dict(predictions=len(predicted(case)), lengths=[len(s) for s in predicted(case)],
                               status=queries['Array_' + case]['mapping_status']) for case in spacers}
assert graph_hash == {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in graph.parent.glob('graph.*')}
(out / 'validation.json').write_text(json.dumps(results, indent=2) + '\n')
print(json.dumps(results, indent=2))
