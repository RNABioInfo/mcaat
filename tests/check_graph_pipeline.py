#!/usr/bin/env python3
"""Exercise the real CLI, saved-cycle replay and parameter validation on a small graph."""
import argparse
import json
from pathlib import Path
import subprocess
import time

p = argparse.ArgumentParser()
p.add_argument('binary', type=Path)
p.add_argument('graph', type=Path)
p.add_argument('output', type=Path)
a = p.parse_args()
a.output.mkdir(parents=True, exist_ok=True)
checks = []

def run(name, args, success=True):
    start = time.monotonic()
    result = subprocess.run([str(a.binary.resolve()), *map(str, args)], text=True, capture_output=True, timeout=60)
    (a.output / (name + '.log')).write_text(result.stdout + result.stderr)
    if (result.returncode == 0) != success:
        raise AssertionError(f'{name}: exit {result.returncode}; see log')
    checks.append({'name': name, 'seconds': time.monotonic() - start, 'exit': result.returncode})
    return result

for i, value in enumerate(('0', '-0.1', '1.1', 'nan', 'inf', '85%', '0.6junk', '85')):
    out = a.output / f'invalid-{i}'
    result = run(f'invalid-{i}', ['--graph', a.graph, '--output-folder', out, '--repeat-support', value], False)
    assert 'repeat-support' in result.stderr and not out.exists()
result = run('missing-support', ['--repeat-support'], False)
assert 'repeat-support' in result.stderr

for flag, values in [('--spacer-variant-similarity', ('nan','0','1.1')), ('--max-spacer-similarity', ('inf','-1','60%')), ('--max-repeat-spacer-similarity', ('bad','1.01')), ('--max-cycles', ('2','-3','500.5','99999999999999999999999'))]:
    for i,value in enumerate(values):
        for prefix in ([],['postprocess-arrays']):
            name=('post-' if prefix else '')+flag[2:]+str(i)
            result=run(name,[*prefix,flag,value],False)
            assert flag[2:] in result.stderr

common = ['--graph', a.graph, '--threads', '1', '--ram', '2G', '--autoclean', 'false']
default = a.output / 'default85'
run('default85', [*common, '--output-folder', default])
assert json.loads((default / 'parameters.json').read_text())['dna_sequence']['repeat_support'] == .85
replay = a.output / 'replay85'
run('replay85', ['postprocess-arrays', '--graph', a.graph, '--cycles', default / 'cycles.txt', '--output-folder', replay])
for filename in ('CRISPR_Arrays_1.txt', 'CRISPR_graph_paths.tsv', 'CRISPR_graph_systems.tsv', 'CRISPR_graph_candidates.tsv', 'CRISPR_graph_summary.tsv', 'CRISPR_graph_decisions.tsv', 'CRISPR_graph_alternatives.tsv', 'CRISPR_graph_review_cycles.tsv', 'CRISPR_spacer_variants.tsv'):
    assert (default / filename).read_bytes() == (replay / filename).read_bytes(), filename

config = a.output / 'support60.settings'
config.write_text('repeat-support=0.60\nspacer-variant-similarity=0.95\nmax-spacer-similarity=0.7\nmax-repeat-spacer-similarity=0.8\nmax-cycles=600\n')
configured = a.output / 'configured60'
run('configured60', [*common, '--output-folder', configured, '--settings', config])
params=json.loads((configured / 'parameters.json').read_text())['dna_sequence']
assert params['repeat_support']==.60 and params['variant_similarity']==.95 and params['max_spacer_similarity']==.7 and params['max_repeat_spacer_similarity']==.8 and params['max_cycles']==600

override = a.output / 'override85'
run('override85', [*common, '--output-folder', override, '--settings', config, '--repeat-support', '.85','--spacer-variant-similarity','.90','--max-spacer-similarity','.60','--max-repeat-spacer-similarity','.60','--max-cycles','500'])
assert json.loads((override / 'parameters.json').read_text())['dna_sequence']['repeat_support'] == .85
assert (override / 'CRISPR_graph_paths.tsv').read_bytes() == (default / 'CRISPR_graph_paths.tsv').read_bytes()
report = {'passed': len(checks), 'checks': checks, 'assertions': 'default; CLI errors; config; CLI overrides; full pipeline vs replay equality'}
(a.output / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
print(f'{len(checks)} CLI/pipeline cases passed')
