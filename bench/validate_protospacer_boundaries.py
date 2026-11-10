#!/usr/bin/env python3
"""Real-SDBG boundary/junction regression; extension requires no read files."""
import argparse,csv,hashlib,json,random,subprocess
from pathlib import Path

p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('output',type=Path)
a=p.parse_args();binary=a.binary.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
if any(out.iterdir()):raise SystemExit('Use an empty output directory')
rng=random.Random(111031)
def dna(n):return ''.join(rng.choice('ACGT') for _ in range(n))
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def rows(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def run(name,args):
    with (out/(name+'.log')).open('w') as log:
        subprocess.run([str(binary),*map(str,args)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
def hashes(prefix):return {f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in prefix.parent.glob(prefix.name+'.*')}
repeat='CCACGCACCCCGCGAGGGGTGCGAC';tail='TTTCAAT'
cores=['ACGT'[i%4]+dna(33)+'ACGT'[(i+1)%4] for i in range(8)]
spacers=[s+tail for s in cores]
targets=[dna(1000)+s+dna(1000) for s in cores]
# Two indistinguishable target contexts for one spacer. Do not force a pairing.
ambiguous=dna(1000)+cores[0]+dna(1000)
hosts=[dna(250)+tail+repeat+s+repeat+dna(250) for s in spacers]
reads=out/'construction_reads.fasta'
with reads.open('w') as f:
    for genome,sequence in enumerate(targets+[ambiguous]+hosts):
        for pos in range(len(sequence)-149):
            for copy in range(4):f.write(f'>{genome}_{pos}_{copy}\n{sequence[pos:pos+150]}\n')
arrays=out/'arrays';arrays.mkdir()
text='>Array_forward\n'+repeat+'\n'+''.join('        '+'-'*len(repeat)+'\t'+s+'\n' for s in spacers)
text+='>Array_reverse\n'+rc(repeat)+'\n'+''.join('        '+'-'*len(repeat)+'\t'+rc(s)+'\n' for s in spacers)
report=arrays/'CRISPR_Arrays_1.txt';report.write_text(text);array_hash=hashlib.sha256(report.read_bytes()).hexdigest()
graph_run=out/'graph_build'
run('build',['--input-files',reads,'--output-folder',graph_run,'--threads','2','--ram','2G','--threshold-multiplicity','1000000','--autoclean','false'])
graph=graph_run/'graph/graph';before=hashes(graph)
common=['phage-curate','--arrays',arrays,'--graph',graph,'--branch-guard','true','--resolve-graph-ties','true']
run('default',[*common,'--output',out/'default'])
result=out/'default';queries=rows(result/'protospacer_queries.tsv');regions=rows(result/'protospacer_regions.tsv')
links=rows(result/'protospacer_associations.tsv');by_source={q['source_id']:q for q in queries};by_region={r['region_id']:r for r in regions}
params=json.loads((result/'protospacer_parameters.json').read_text())
assert params['schema']==6 and params['branch_guard_active'] and params['boundary_corrected_source_observations']==16
assert params['max_flank_length']==750 and params['boundary_support']==.85
for q in queries:
    assert q['boundary_status']=='conserved_terminal_reassigned'
    assert int(q['trimmed_left_bp'])+int(q['trimmed_right_bp'])==7
    assert len(q['search_spacer'])==35 and len(q['spacer'])==42
    assert len(q['search_repeat'])==len(repeat)+7
    if q['spacer'] in [spacers[0],rc(spacers[0])]:
        assert q['mapping_status']=='unresolved', 'indistinguishable source contexts forced into a target'
    else:
        assert q['mapping_status']=='candidate', 'unambiguous true target lost'
assert len(regions)==7
for row in regions:
    seq=row['sequence'];assert any(seq in t or rc(seq) in t for t in targets[1:])
    assert len(seq)==len(row['vertices'].split(','))+22
    assert len(row['edges'].split(','))+1==len(row['vertices'].split(','))
for link in links:
    q=by_source[link['source_id']];seq=by_region[link['region_id']]['sequence']
    lo,hi=int(link['seed_begin']),int(link['seed_end'])
    assert lo<=750 and len(seq)-hi<=750, 'per-side cap exceeded'
    expected=q['spacer'][int(link['query_begin']):int(link['query_end'])]
    observed=seq[lo:hi];assert (rc(observed) if link['strand']=='-' else observed)==expected
    assert int(link['expanded_states'])<=200000
assert any(len(r['sequence'])==1535 for r in regions), 'default 750+35+750 extension not retained'
checks=rows(result/'protospacer_branch_checks.tsv')
assert any(r['decision']=='ambiguous_branch' for r in checks)
assert any(r['status']=='source_array_context' for r in checks)
assert all(int(r['probe_steps'])<=75 for r in checks)
run('threads24',[*common,'--output',out/'threads24','--threads','24'])
for file in result.iterdir():
    if file.name!='protospacer_parameters.json':assert file.read_bytes()==(out/'threads24'/file.name).read_bytes(),file.name
run('disabled',[*common,'--output',out/'disabled','--correct-spacer-boundaries','false'])
assert all(q['search_spacer']==q['spacer'] and q['boundary_status']=='disabled' for q in rows(out/'disabled/protospacer_queries.tsv'))
# Exercise TSV input and the same conserved-boundary grouping.
tsv=out/'tsv_arrays';tsv.mkdir()
with (tsv/'CRISPR_graph_paths.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['array_id','observation_id','repeat','spacer'])
    for i,s in enumerate(spacers):w.writerow(['Array_forward',i,repeat,s])
run('tsv',[*common[:1],'--arrays',tsv,'--graph',graph,'--output',out/'tsv','--branch-guard','true','--resolve-graph-ties','true'])
assert {r['sequence'] for r in rows(out/'tsv/protospacer_regions.tsv')}=={r['sequence'] for r in regions}
assert hashes(graph)==before and hashlib.sha256(report.read_bytes()).hexdigest()==array_hash
summary=dict(corrected_observations=16,true_targets_retained=7,ambiguous_query_withheld=1,
             default_region_length=1535,threads_1_and_24_identical=True,original_coordinates_verified=True,
             graph_and_array_files_unchanged=True)
(out/'validation.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
