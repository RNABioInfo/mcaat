#!/usr/bin/env python3
"""Production SDBG: preserve supplied detections independently of extension quality."""
import argparse,csv,hashlib,json,random,subprocess
from pathlib import Path
from check_retained_protospacers import check as check_outputs
p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('output',type=Path)
a=p.parse_args();binary=a.binary.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
if any(out.iterdir()):raise SystemExit('Use an empty output directory')
rng=random.Random(20261111)
def dna(n):return ''.join(rng.choice('ACGT') for _ in range(n))
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def rows(path):return list(csv.DictReader(path.open(),delimiter='\t'))
def fasta(path):
 result={};name=None
 for line in path.read_text().splitlines():
  if line.startswith('>'):
   name=line[1:].split()[0];assert name not in result;result[name]=''
  else:result[name]+=line.strip()
 return result
def run(name,args,success=True):
 with (out/(name+'.log')).open('w') as f:
  r=subprocess.run([str(binary),*map(str,args)],stdout=f,stderr=subprocess.STDOUT,timeout=180)
 assert (r.returncode==0)==success,name

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
seeds={k:dna(35) for k in ['straight','branch','adapter','dead_end','absent','unknown_source']}
seeds['below_k']=dna(17)
seeds['variant']=dna(37)
variant_allele=list(seeds['variant'])
for position in (0,2,25,35):variant_allele[position]='ACGT'[('ACGT'.index(variant_allele[position])+1)%4]
variant_allele=''.join(variant_allele)
variant_truth=dna(1000)+variant_allele+dna(1000)
repeat=dna(27);adapter='AGATCGGAAGAGCACACGT'
left,right=dna(1000),dna(1000);straight=left+seeds['straight']+right
branch_left,common=dna(900),dna(20)
arm_a='A'+dna(8);arm_b='C'+dna(8)
true_branch=branch_left+seeds['branch']+common+arm_a+'A'+dna(900)
branch_other=seeds['branch']+common+arm_b+'A'+dna(900)
# Both first continuations reach another fork after nine edges. The old local
# gate stops early; complete-path scoring can see the unstable alternatives.
branch_a_other=seeds['branch']+common+arm_a+'C'+dna(900)
branch_b_other=seeds['branch']+common+arm_b+'C'+dna(900)
adapter_left=dna(800);adapter_target=adapter_left+seeds['adapter']+adapter+dna(800)
sequences=[(straight,4),(true_branch,4),(branch_other,4),(branch_a_other,4),(branch_b_other,4),(adapter_target,4),(seeds['dead_end'],4),(variant_truth,4)]
seeds['masked_variant']='AGGAGGGGGCTCTGCCGCACGGCGC'
masked_allele='AGCAGGGGGCTCTGCTGCACCGCGC'
masked_truth=dna(1000)+masked_allele+dna(1000)
sequences.append((masked_truth,4))
reads=out/'construction.fasta'
with reads.open('w') as f:
 for num,(seq,copies) in enumerate(sequences):
  if len(seq)<150:
   for j in range(copies):f.write(f'>{num}_short_{j}\n{seq}\n')
   continue
  for pos in range(len(seq)-149):
   n=copies
   if num in (2,3,4) and pos>len(seeds['branch'])+len(common):n=1 if (pos//55)%2 else 18
   for j in range(n):f.write(f'>{num}_{pos}_{j}\n{seq[pos:pos+150]}\n')
arrays=out/'arrays';arrays.mkdir()
source_sequences={}
with (arrays/'CRISPR_Arrays_1.txt').open('w') as f:
 for name,s in seeds.items():
  if name=='unknown_source':continue
  full=dna(8)+s+dna(9) if name=='below_k' else s
  source_sequences['CRISPR_Arrays_1.txt:Array_'+name+':1']=full
  f.write('>Array_'+name+'\n'+repeat+'\n        '+'-'*len(repeat)+'\t'+full+'\n')
seeds['reverse']=rc(seeds['straight']);seeds['alias']=seeds['straight']
provided=out/'supplied.fasta';provided.write_text(''.join('>'+n+'\n'+s+'\n' for n,s in seeds.items()))
input_hash=digest(provided);array_hash=digest(arrays/'CRISPR_Arrays_1.txt')
run('build',['--input-files',reads,'--output-folder',out/'graph_build','--threads','2','--ram','2G','--threshold-multiplicity','1000000','--autoclean','false'])
graph=out/'graph_build/graph/graph';graph_hash={f.name:digest(f) for f in graph.parent.glob('graph.*')}
common_args=['phage-curate','--arrays',arrays,'--graph',graph,'--protospacers',provided]
run('default',[*common_args,'--output',out/'default'])
result=out/'default';extensions={r['seed_id']:r for r in rows(result/'protospacer_extensions.tsv')}
check_outputs(provided,result)
assert fasta(result/'Protospacers.fasta')==seeds
assert list(fasta(result/'Protospacers.fasta'))==list(seeds), 'input IDs/order changed'
assert set(extensions)==set(seeds)
assert all(r['retention_status']=='retained_input' for r in extensions.values())
for source in rows(result/'protospacer_sources.tsv'):
 full=source_sequences[source['source_id']]
 observed=full[int(source['query_begin']):int(source['query_end'])]
 assert (rc(observed) if source['strand']=='-' else observed)==seeds[source['seed_id']], 'source coordinates disagree'
assert extensions['below_k']['reason']=='no_supported_variant_extension'
assert extensions['absent']['reason']=='no_supported_variant_extension'
assert extensions['unknown_source']['reason']=='source_context_missing'
assert extensions['dead_end']['extension_status']=='not_extended'
regions={r['region_id']:r for r in rows(result/'protospacer_regions.tsv')}
links=rows(result/'protospacer_associations.tsv')
assert len(links)==len({r['seed_id'] for r in links})
for link in links:
 r=regions[link['region_id']];seq=r['sequence'];lo,hi=int(link['seed_begin']),int(link['seed_end'])
 expected={'variant':variant_allele,'masked_variant':masked_allele}.get(link['seed_id'],seeds[link['seed_id']])
 assert (rc(seq[lo:hi]) if link['strand']=='-' else seq[lo:hi])==expected
 assert lo<=750 and len(seq)-hi<=750 and len(seq)<=1550
 assert len(seq)==len(r['vertices'].split(','))+22
 assert len(r['edges'].split(','))+1==len(r['vertices'].split(','))
 assert int(extensions[link['seed_id']]['expanded_states'])<=100000
for name in ['straight','reverse','alias']:
 e=extensions[name];assert int(e['length'])==1535,(name,e)
 seq=regions[e['region_id']]['sequence'];assert seq in straight or rc(seq) in straight
assert int(extensions['adapter']['length'])>=700
assert all(adapter not in r['sequence'] and rc(adapter) not in r['sequence'] for r in regions.values())
assert int(extensions['branch']['length'])>=1000,extensions['branch']
branch_sequence=regions[extensions['branch']['region_id']]['sequence']
assert branch_sequence in true_branch or rc(branch_sequence) in true_branch, 'longer-context ranking chose unstable branch'
variant=extensions['variant'];assert variant['extension_status'].startswith('extended_variant')
variant_sequence=regions[variant['region_id']]['sequence']
assert len(variant_sequence)==1537 and (variant_sequence in variant_truth or rc(variant_sequence) in variant_truth)
masked=extensions['masked_variant'];assert masked['extension_status'].startswith('extended_variant')
masked_sequence=regions[masked['region_id']]['sequence']
assert len(masked_sequence)==1525 and (masked_sequence in masked_truth or rc(masked_sequence) in masked_truth)
assert any(r['seed_id']=='masked_variant' and r['phase']=='masked_extend' for r in rows(result/'protospacer_placement_checks.tsv'))
run('exact_only',[*common_args,'--output',out/'exact_only','--target-variant-fallback','false'])
exact_only={r['seed_id']:r for r in rows(out/'exact_only/protospacer_extensions.tsv')}
assert exact_only['below_k']['reason']=='shorter_than_graph_k'
assert exact_only['absent']['reason']=='seed_not_mapped'
assert exact_only['variant']['reason']=='seed_not_mapped'
assert exact_only['masked_variant']['reason']=='seed_not_mapped'
check_outputs(provided,out/'exact_only')
# Twenty-four threads must yield exactly the same records, paths and diagnostics.
run('threads24',[*common_args,'--output',out/'threads24','--threads','24'])
for f in result.iterdir():
 if f.name!='protospacer_parameters.json':assert f.read_bytes()==(out/'threads24'/f.name).read_bytes(),f.name
# Even user-requested strict clipping / tiny budgets cannot delete the inputs.
run('strict',[*common_args,'--output',out/'strict','--branch-guard','true','--resolve-graph-ties','true'])
assert fasta(out/'strict/Protospacers.fasta')==seeds
strict={r['seed_id']:r for r in rows(out/'strict/protospacer_extensions.tsv')}
assert int(strict['branch']['length'])<int(extensions['branch']['length']), 'fixture does not reproduce premature stopping'
run('limited',[*common_args,'--output',out/'threads24','--threads','24','--max-query-states','2'])
assert fasta(out/'threads24/Protospacers.fasta')==seeds
assert not (out/'threads24/ContiguousPaths.fasta').read_text(), 'seed-only paths or stale extensions exported'
check_outputs(provided,out/'threads24')
snapshot={f.name:f.read_bytes() for f in result.iterdir()}
duplicate=out/'duplicate.fasta';duplicate.write_text('>same\n'+seeds['straight']+'\n>same\n'+seeds['branch']+'\n')
run('invalid',[*common_args,'--protospacers',duplicate,'--output',result],success=False)
assert snapshot=={f.name:f.read_bytes() for f in result.iterdir()}, 'failed run replaced successful results'
run('input_overwrite',[*common_args,'--output',out],success=False)
# Real saved 81-match baseline. This tests retention against an unrelated graph,
# deliberately not the original metagenome's reconstruction accuracy.
baseline=Path(__file__).resolve().parent/'fixtures/protospacer-baseline-81/Protospacers.fasta'
run('baseline',[*common_args,'--protospacers',baseline,'--output',out/'baseline81','--threads','24'])
assert len(fasta(out/'baseline81/Protospacers.fasta'))==81
assert (out/'baseline81/Protospacers.fasta').read_bytes()==baseline.read_bytes()
check_outputs(baseline,out/'baseline81')
assert digest(provided)==input_hash and digest(arrays/'CRISPR_Arrays_1.txt')==array_hash
assert graph_hash=={f.name:digest(f) for f in graph.parent.glob('graph.*')}
summary=dict(supplied_records=len(seeds),retained_records=len(seeds),baseline_records_preserved=81,
 straight_length=int(extensions['straight']['length']),branch_length=int(extensions['branch']['length']),
 old_strict_branch_length=int(strict['branch']['length']),adapter_length=int(extensions['adapter']['length']),
 true_branch_recovered=True,variant_target_recovered=True,variant_length=len(variant_sequence),
 masked_target_recovered=True,masked_length=len(masked_sequence),
 exact_only_policy_checked=True,threads_1_24_identical=True,failed_extensions_keep_protospacers=True,
 original_metagenome_graph_tested=False)
(out/'validation.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
