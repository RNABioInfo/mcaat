#!/usr/bin/env python3
"""Production negative/positive controls for the source-array leakage regression.

Uses a reported repeat/spacer and the flanking repeat variants observed in the
failing run. Separate cases cover indels and exact copies. All remaining host
and target bases are deterministic synthetic DNA; no biological accuracy claim.
"""
import argparse,csv,json,random,subprocess,time
from pathlib import Path

parser=argparse.ArgumentParser()
parser.add_argument('binary',type=Path)
parser.add_argument('output',type=Path)
args=parser.parse_args()
binary=args.binary.resolve();out=args.output.resolve()
out.mkdir(parents=True,exist_ok=True)
if any(out.iterdir()):raise SystemExit('Use an empty output directory')
rng=random.Random(92141)
def dna(n):return ''.join(rng.choice('ACGT') for _ in range(n))
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def table(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def run(name,arguments):
    start=time.monotonic()
    with (out/(name+'.log')).open('w') as log:
        subprocess.run([str(binary),*map(str,arguments)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=180)
    return time.monotonic()-start

repeat='CGGTTCATCCCCGCTTGCGCGGGGAACTC'
spacer='ACTACTTCGCTGGAAACGCTAGAGGTGCGTTTT'
left='CGGTTCATCCCCGCTGGCGCGGGGAACTC'
right='CGGGTAATCCCCGCTTGCGCGGGGAACTC'
indel_spacer=dna(35);exact_spacer=dna(34)
cases=[(spacer,left,right),(indel_spacer,repeat[:10]+repeat[11:],repeat[:15]+'A'+repeat[15:]),(exact_spacer,repeat,repeat)]
hosts=[dna(1500)+l+s+r+dna(1500) for s,l,r in cases]
target=dna(1500)+spacer+dna(1500)
arrays=out/'arrays';arrays.mkdir()
report=''
for i,(s,l,r) in enumerate(cases,1):report+=f'>Array_{i}\n{repeat}\n        '+len(repeat)*'-'+'\t'+s+'\n'
report+=f'>Array_RC\n{rc(repeat)}\n        '+len(repeat)*'-'+'\t'+rc(spacer)+'\n'
(arrays/'CRISPR_Arrays_1.txt').write_text(report)
positive_arrays=out/'positive_arrays';positive_arrays.mkdir()
(positive_arrays/'CRISPR_Arrays_1.txt').write_text(
    f'>Array_1\n{repeat}\n        '+len(repeat)*'-'+'\t'+spacer+'\n'+
    f'>Array_RC\n{rc(repeat)}\n        '+len(repeat)*'-'+'\t'+rc(spacer)+'\n')
results={}
for name,genomes in [('negative',hosts),('positive',[hosts[0],target])]:
    reads=out/(name+'.fasta')
    with reads.open('w') as f:
        for genome,sequence in enumerate(genomes):
            for start in range(len(sequence)-149):
                for copy in range(12):f.write(f'>{genome}_{start}_{copy}\n{sequence[start:start+150]}\n')
    graph_run=out/(name+'_graph')
    run(name+'_build',['--input-files',reads,'--output-folder',graph_run,'--threads','2','--ram','2G','--threshold-multiplicity','1000000','--autoclean','false'])
    result=out/name
    active_arrays=arrays if name=='negative' else positive_arrays
    elapsed=run(name,['phage-curate','--arrays',active_arrays,'--graph',graph_run/'graph/graph','--output',result,'--threads','2'])
    queries=table(result/'protospacer_queries.tsv');links=table(result/'protospacer_associations.tsv')
    regions=table(result/'protospacer_regions.tsv');seeds=table(result/'protospacer_candidates.tsv')
    assert all(q['mapping_status']!='not_fully_mapped' for q in queries),'control input missing from graph'
    if name=='negative':
        assert not regions and not seeds and not links,'host-only CRISPR array reported as a target'
        assert all(q['mapping_status']=='no_external_flank' for q in queries)
        assert not (result/'Protospacers.fasta').read_text()
        assert any(r['reason']=='source_array_context' for r in table(result/'protospacer_decisions.tsv'))
        # Raising identity to 1 disables approximate exclusion and must reproduce
        # the original leak for variant copies. Exact copies stay excluded.
        strict=out/'exact_only'
        run('exact_only',['phage-curate','--arrays',arrays,'--graph',graph_run/'graph/graph','--output',strict,'--array-repeat-identity','1'])
        assert any(q['spacer']==spacer and q['mapping_status']=='candidate' for q in table(strict/'protospacer_queries.tsv'))
    else:
        assert any(q['spacer']==spacer and q['mapping_status']=='candidate' for q in queries),'true target lost together with array path'
        assert all(q['mapping_status']=='no_external_flank' for q in queries if q['spacer'] in (indel_spacer,exact_spacer))
        assert regions and all(r['sequence'] in target or rc(r['sequence']) in target for r in regions),'source/chimeric context survived in positive control'
        assert any(int(r['length'])>=1000 for r in regions),'target reduced to tiny fragment'
        assert {r['array_id'] for r in links}=={'Array_1','Array_RC'},'opposite-strand source lost'
        assert all(r['status']!='seed_only' for r in links)
    results[name]=dict(unique_regions=len(regions),seed_intervals=len(seeds),elapsed_seconds=elapsed)
(out/'validation.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
