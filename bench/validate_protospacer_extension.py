#!/usr/bin/env python3
"""Build a controlled real SDBG and check the production protospacer CLI.

Truth includes a host array, two target contexts (one one-sided), a host-only
spacer, a missing query, opposite-strand source associations, and read/PAM input.
150-base synthetic reads build the SDBG; separate read pairs exercise evidence
scoring. This is software validation, not biological accuracy.
"""
import argparse,csv,gzip,hashlib,json,random,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('output',type=Path)
p.add_argument('--graph-builder',type=Path);a=p.parse_args()
a.binary=a.binary.resolve();a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=True)
if any(a.output.iterdir()):raise SystemExit('Use an empty validation output directory')
rng=random.Random(77031)
def dna(n):return ''.join(rng.choice('ACGT') for _ in range(n))
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def table(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def run(name,args,success=True):
    started=time.monotonic();log=a.output/(name+'.log')
    with log.open('w') as f:r=subprocess.run([str(a.binary),*map(str,args)],stdout=f,stderr=subprocess.STDOUT,timeout=180)
    if (r.returncode==0)!=success:raise AssertionError(f'{name}: exit {r.returncode}; inspect {log}')
    return time.monotonic()-started
repeat=dna(36);s1=dna(34);s2=s1[:23]+dna(12);s3=dna(36);missing=dna(33)
host=dna(100)+repeat+s1+repeat+s2+repeat+s3+repeat+dna(100)
left=dna(140);right='AGG'+dna(157);target=left+s1+right
target2=dna(130)+s2
genomes={'host':(host,30),'target':(target,12),'one_sided_target':(target2,8)}
truth={name:s for name,(s,_) in genomes.items()}
(a.output/'truth.json').write_text(json.dumps(truth,indent=2)+'\n')
reads=a.output/'edge_reads.fasta'
with reads.open('w') as f:
    for name,(seq,depth) in genomes.items():
        for pos in range(len(seq)-149):
            for copy in range(depth):f.write(f'>{name}_{pos}_{copy}\n{seq[pos:pos+150]}\n')
arrays=a.output/'arrays';arrays.mkdir()
report=arrays/'CRISPR_Arrays_1.txt'
report.write_text('>Array_1\n'+repeat+'\n'+''.join('        '+'-'*len(repeat)+'\t'+s+'\n' for s in (s1,s2,s3,missing))+
                  '>Array_2\n'+rc(repeat)+'\n        '+rc(repeat)+'\t'+rc(s1)+'\n')
input_hash=hashlib.sha256(report.read_bytes()).hexdigest()
builder=(a.graph_builder or a.binary).resolve();graph_run=a.output/'graph_build'
with (a.output/'graph_build.log').open('w') as log:
    subprocess.run([str(builder),'--input-files',str(reads),'--output-folder',str(graph_run),'--threads','2','--ram','2G','--threshold-multiplicity','1000000','--autoclean','false'],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
graph=graph_run/'graph/graph'
if not graph.with_suffix('.sdbg_info').exists():raise AssertionError('graph not built')
common=['phage-curate','--arrays',arrays,'--graph',graph,'--max-length','500','--coverage-window','3']
out=a.output/'default';runtime=run('default',[*common,'--output',out])
queries=table(out/'protospacer_queries.tsv');seeds=table(out/'protospacer_candidates.tsv');regions=table(out/'protospacer_regions.tsv');links=table(out/'protospacer_associations.tsv')
assert len(queries)==5 and len({r['source_id'] for r in queries})==5
by_source={r['source_id']:r for r in queries};by_seed={r['seed_id']:r for r in seeds};by_region={r['region_id']:r for r in regions}
assert any(r['spacer']==missing and r['mapping_status']=='not_fully_mapped' for r in queries)
assert any(r['spacer']==s3 and r['mapping_status']=='no_external_flank' for r in queries)
assert links and regions
assert all(r['status'] in ('one_sided','two_sided') for r in links),'unextended seed leaked into detections'
params=json.loads((out/'protospacer_parameters.json').read_text())
assert params['matched_queries']==len({r['query_id'] for r in seeds})
assert len(seeds)==params['matched_queries'] and len(regions)<=len(seeds),'more than one final prediction per spacer'
selection=table(out/'protospacer_selection.tsv')
assert len(selection)==len(seeds) and len({r['query_id'] for r in selection})==len(selection)
assert any(int(r['compared_regions'])>1 for r in selection) or table(out/'protospacer_branch_checks.tsv'), 'no alternatives inspected during traversal or final selection'
assert all(sum(r['source_id']==q['source_id'] for r in links)<=1 for q in queries),'multiple predictions for a source spacer'
assert (out/'MatchedSpacers.fasta').read_text().count('>')==params['matched_queries']
for r in regions:
    assert int(r['length'])==len(r['sequence'])<=500
    assert len(r['sequence'])==len(r['vertices'].split(','))+22
    assert len(r['edges'].split(','))==len(r['vertices'].split(','))-1
assert len({min(r['sequence'],rc(r['sequence'])) for r in regions})==len(regions),'duplicate strand regions'
for link in links:
    original=by_source[link['source_id']]['spacer'];region=by_region[link['region_id']]['sequence']
    matched=original[int(link['query_begin']):int(link['query_end'])]
    observed=region[int(link['seed_begin']):int(link['seed_end'])]
    assert matched==(rc(observed) if link['strand']=='-' else observed),'wrong seed/source coordinates'
    assert link['read_evidence']=='not_provided'
    assert link['pam_status']=='unknown_no_model'
sources_s1={r['source_id'] for r in queries if r['spacer'] in (s1,rc(s1))}
assert len(sources_s1)==2
assert any(sources_s1<={r['source_id'] for r in links if r['region_id']==region['region_id']} for region in regions),'shared target lost source associations'
for source in [s1,s2]:assert any(r['spacer']==source and r['mapping_status']=='candidate' for r in queries),'same first kmer overwrote a spacer'
assert any(r['status']=='one_sided' and by_source[r['source_id']]['spacer']==s2 for r in links)
truth_sequences=[s for s in truth.values()]+[rc(s) for s in truth.values()]
exact_regions=sum(any(r['sequence'] in truth for truth in truth_sequences) for r in regions)
assert exact_regions==len(regions), 'junction ambiguity produced a recombinant in this fixture'
# This deliberately ambiguous graph shares a 23-mer between different spacers.
# Both retained regions must now match a known source in this fixture. This is
# a regression check, not a claim that graph-only scoring proves genomic origin.
assert all(len(r['sequence'])>len(by_seed[next(x['seed_id'] for x in links if x['region_id']==r['region_id'])]['sequence']) for r in regions)
# Deterministic output, no append duplicates, no edits to array inputs.
snapshot={f.name:f.read_bytes() for f in out.iterdir()}
run('repeat',[*common,'--output',out]);assert snapshot=={f.name:f.read_bytes() for f in out.iterdir()}
parallel=a.output/'parallel';run('parallel',[*common,'--output',parallel,'--threads','2'])
for name,content in snapshot.items():
    if name!='protospacer_parameters.json':assert (parallel/name).read_bytes()==content,name
assert hashlib.sha256(report.read_bytes()).hexdigest()==input_hash
# Optional gzipped paired reads: exact inward-facing matches and observed PAM.
r1=a.output/'R1.fastq.gz';r2=a.output/'R2.fastq.gz'
with gzip.open(r1,'wt') as f,gzip.open(r2,'wt') as g:
    for i,pos in enumerate(range(60,151,10)):
        x=target[pos:pos+75];y=rc(target[pos+100:pos+175]);assert len(x)==len(y)==75
        f.write(f'@pair{i}/1\n{x}\n+\n'+75*'I'+'\n');g.write(f'@pair{i}/2\n{y}\n+\n'+75*'I'+'\n')
pam=a.output/'pam.tsv';pam.write_text('array_id\tsubtype\tstrand\tside\tmotif\nArray_1\tsynthetic-NGG\t+\t3prime\tNGG\nArray_2\tsynthetic-NGG\t-\t3prime\tNGG\n')
supported=a.output/'supported'
run('supported',[*common,'--output',supported,'--reads1',r1,'--reads2',r2,'--min-insert','150','--max-insert','200','--pam-rules',pam])
supported_links=table(supported/'protospacer_associations.tsv')
assert len(table(supported/'protospacer_candidates.tsv'))==len({r['query_id'] for r in table(supported/'protospacer_candidates.tsv')})
assert any(int(r['pair_support'])>0 for r in supported_links),'mate evidence absent'
assert any(r['pam_status']=='match' for r in supported_links),'PAM evidence absent'
assert any(int(r['read_support'])>0 for r in supported_links),'junction read evidence absent'
supported_regions={r['region_id']:r['sequence'] for r in table(supported/'protospacer_regions.tsv')}
full_seeds={r['seed_id'] for r in seeds if float(r['query_coverage'])==1 and r['sequence'] in (s1,rc(s1))}
for source in sources_s1:
    top=[r for r in supported_links if r['source_id']==source and r['seed_id'] in full_seeds and r['alternative_rank']=='1']
    assert top and all(int(r['pair_support'])>0 for r in top),'supported target not ranked first'
    assert all(supported_regions[r['region_id']] in target or rc(supported_regions[r['region_id']]) in target for r in top),'top supported full-spacer path is chimeric'
# Both current graph reports and original text reports use the production command.
tsv=a.output/'tsv_arrays';tsv.mkdir()
(tsv/'CRISPR_graph_paths.tsv').write_text('array_id\tobservation_id\trepeat\tspacer\n'+''.join(f"{q['array_id']}\t{q['observation_id']}\t{rc(repeat) if q['array_id']=='Array_2' else repeat}\t{q['spacer']}\n" for q in queries))
run('tsv',['protospacer-curate','--arrays',tsv,'--graph',graph,'--output',a.output/'tsv','--max-length','500','--coverage-window','3'])
assert (a.output/'tsv/ContiguousPaths.fasta').read_bytes()==(out/'ContiguousPaths.fasta').read_bytes(),'input formats give different regions'
limited=a.output/'limited';run('limited',[*common,'--output',limited,'--max-query-states','2'])
assert table(limited/'protospacer_unresolved.tsv'),'cap not disclosed separately'
assert not table(limited/'protospacer_associations.tsv'),'unresolved seed counted as a detection'
assert not (limited/'Protospacers.fasta').read_text() and (limited/'UnresolvedSeeds.fasta').read_text()
assert any(r['mapping_status']=='unresolved' for r in table(limited/'protospacer_queries.tsv'))
short=a.output/'short';run('short',[*common,'--output',short,'--max-length','36'])
short_pending=table(short/'protospacer_unresolved.tsv')
assert all(any(r['sequence'] in (s,rc(s)) for r in short_pending) for s in (s1,s2,s3)), 'complete-spacer seeds lack enough flank context at this length'
short_decisions=table(short/'protospacer_decisions.tsv')
for row in short_pending:
    assert any(d['seed_id']==row['seed_id'] and
               {'length_limit','ambiguous_branch'} & {d['left_stop'],d['right_stop']}
               for d in short_decisions), 'short-context failure lacks length or early-ambiguity stop'
assert json.loads((short/'protospacer_parameters.json').read_text())['full_length_matched_queries']==0
# A shorter seed can still have enough bases to establish a genuinely different
# flank within 36 bp. Do not turn that positive partial match into a negative.
assert all(r['status'] in ('one_sided','two_sided') for r in table(short/'protospacer_associations.tsv'))
for i,flags in enumerate([['--beam-width','0'],['--coverage-model','garbage'],['--max-coverage-penalty','nan'],['--reads2',r2],['--threads','-1'],['--max-length','10'],['--array-repeat-identity','1.1'],
                          ['--max-flank-length','0'],['--branch-lookahead','1001'],['--branch-score-margin','nan'],['--prediction-score-margin','-1'],
                          ['--boundary-support','.5'],['--boundary-min-spacers','1'],['--correct-spacer-boundaries','yes']]):
    run(f'invalid-{i}',[*common,'--output',a.output/f'invalid-{i}',*flags],False)
bad=a.output/'bad_R2.fq';bad.write_text('@wrong/2\nACGT\n+\nIIII\n')
run('bad-mates',[*common,'--output',supported,'--reads1',r1,'--reads2',bad],False)
assert table(supported/'protospacer_associations.tsv')==supported_links,'failed run replaced successful results'
run('retention-limit',[*common,'--output',supported,'--max-retained-bases','1'],False)
assert table(supported/'protospacer_associations.tsv')==supported_links,'memory limit replaced successful results'
run('read-work-limit',[*common,'--output',supported,'--reads1',r1,'--reads2',r2,'--max-read-work','1'],False)
assert table(supported/'protospacer_associations.tsv')==supported_links,'read-work limit replaced successful results'
empty=a.output/'host_only';empty.mkdir();(empty/'CRISPR_Arrays_1.txt').write_text(f'>Array_1\n{repeat}\n        {repeat}\t{s3}\n')
run('empty-rerun',['phage-curate','--arrays',empty,'--graph',graph,'--output',parallel])
assert (parallel/'ContiguousPaths.fasta').read_text()=='' and (parallel/'Protospacers.fasta').read_text()=='','empty rerun left stale candidates'
run('resolved-rerun',[*common,'--output',limited])
assert (limited/'UnresolvedSeeds.fasta').read_bytes()==(out/'UnresolvedSeeds.fasta').read_bytes(),'rerun left stale unresolved intervals'
summary={'checks':'production SDBG mapping, exact reconstruction, partial targets, overlapping sources, reverse complements, all coordinates, deterministic reruns and threads, read pairs, PAM, budgets, input errors',
         'sources':len(queries),'seeds':len(seeds),'regions':len(regions),'regions_wholly_matching_a_truth_sequence':exact_regions,
         'regions_not_wholly_matching_a_truth_sequence':len(regions)-exact_regions,'associations':len(links),'runtime_seconds':runtime,
         'one_prediction_per_distinct_spacer':True,'shortlisted_paths_compared':sum(int(r['compared_paths']) for r in selection),
         'top_read_supported_full_spacer_target_matches_truth':True,
         'interpretation':'Controlled software fixture. The multiplicity-only winner can be chimeric; graph connectivity alone is not viral confirmation.'}
(a.output/'validation.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
