#!/usr/bin/env python3
"""Compare saved-cycle filtering with annotated local reference spacer sequences.
Exact strand-independent sequence recovery is reported, not locus precision.
"""
import argparse,csv,json,subprocess,time
from collections import defaultdict
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('old',type=Path);p.add_argument('new',type=Path);p.add_argument('data',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
a.output.mkdir(parents=True,exist_ok=True)
def records(file):
    name=None;seq=[]
    for line in file.read_text().splitlines():
        if line.startswith('>'):
            if name is not None:yield name,''.join(seq).upper()
            name=line[1:];seq=[]
        else:seq.append(line.strip())
    if name is not None:yield name,''.join(seq).upper()
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def key(s):return min(s,rc(s))
def circles(s):return min(s[i:]+s[:i] for s in (s,rc(s)) for i in range(len(s)))
def table(file):
    with file.open() as f:return list(csv.DictReader(f,delimiter='\t'))
results=[]
for cycles in sorted(a.data.glob('*/cycles/str_paths.txt')):
    base=cycles.parent.parent;graph=base/'graph/graph'
    if not graph.with_suffix('.sdbg_info').exists():continue
    refs={}
    for file in sorted((base/'reference').glob('*')):
        if not file.is_file():continue
        rec=list(records(file));spacers=frozenset(key(s) for name,s in rec if name.lower().startswith('spacer') and s)
        if len(spacers)>=2:refs.setdefault(spacers,file.name)
    known=set().union(*refs) if refs else set()
    result={'dataset':base.name,'reference_arrays':len(refs),'reference_spacers':len(known),'versions':{}}
    evidence={}
    for version,binary in [('before',a.old),('after',a.new)]:
        out=a.output/base.name/version;out.mkdir(parents=True,exist_ok=True);log=out/'run.log';started=time.monotonic()
        with log.open('w') as f:
            run=subprocess.run([str(binary.resolve()),'postprocess-arrays','--graph',str(graph),'--cycles',str(cycles),'--output-folder',str(out)],stdout=f,stderr=subprocess.STDOUT)
        if run.returncode:raise RuntimeError(f'{base.name} {version}: see {log}')
        paths=table(out/'CRISPR_graph_paths.tsv');sequences={key(r['spacer']) for r in paths};evidence[version]=sequences
        per_array=defaultdict(set)
        for r in paths:per_array[r['array_id']].add(key(r['spacer']))
        shared=defaultdict(set)
        for r in paths:shared[circles(r['repeat']+r['spacer'])].add(r['array_id'])
        duplicates=sum(len(v)>1 for v in shared.values())
        result['versions'][version]={'seconds':time.monotonic()-started,'arrays':len(per_array),'observations':len(paths),'distinct_spacers':len(sequences),'exact_reference_spacers':len(known&sequences),'reference_arrays_with_two_exact_spacers':sum(any(len(s&ref)>=2 for s in per_array.values()) for ref in refs),'shared_circles':duplicates}
        if version=='after':
            assert duplicates==0,(base.name,'duplicate accepted cycle')
            variants=table(out/'CRISPR_spacer_variants.tsv')
            assert {(r['array_id'],r['variant']) for r in variants}=={(r['array_id'],r['spacer']) for r in paths}
            review=table(out/'CRISPR_graph_alternatives.tsv');decisions=table(out/'CRISPR_graph_decisions.tsv')
            assert len(paths)+len(review)==sum(int(r['original_assigned']) for r in decisions)
            observed={(r['array_id'],r['spacer']) for r in paths}
            assert all((r['array_id'],r['representative']) in observed for r in variants)
    result['lost_exact_reference_spacers']=sorted(known&evidence['before']-evidence['after'])
    results.append(result);(a.output/'reference_validation.json').write_text(json.dumps(results,indent=2)+'\n')
    print(base.name, result['versions']['before']['exact_reference_spacers'],'->',result['versions']['after']['exact_reference_spacers'],'known spacers;',result['versions']['before']['arrays'],'->',result['versions']['after']['arrays'],'candidates',flush=True)
