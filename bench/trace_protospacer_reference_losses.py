#!/usr/bin/env python3
"""Locate witness paths in logging-only production traces; reference scoring stays external."""
import argparse,csv,collections,math
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--trace',type=Path,required=True);p.add_argument('--oracle',type=Path,required=True);p.add_argument('--metrics',type=Path,required=True);p.add_argument('--repeat-identity',type=float,default=.85);a=p.parse_args()
if not 0<a.repeat_identity<=1:p.error('--repeat-identity must be in (0, 1]')
trace=a.trace;oracle=a.oracle
def rows(p):return list(csv.DictReader(p.open(),delimiter='\t'))
contexts=collections.defaultdict(list)
for r in rows(trace/'production_contexts.tsv'):contexts[r['seed_id']].append(r)
def boundary_context(witness):
 # Independently evaluate the complete-repeat edit rule on the witness string.
 # A short connected prefix can still end before source/target disambiguation.
 pending=False;outward=witness['sequence_outward']
 for c in contexts[witness['seed_id']]:
  back=witness['side']=='left'
  remainder=c['spacer'][:int(c['begin'])] if back else c['spacer'][int(c['end']):]
  repeat=c['repeat'][::-1] if back else c['repeat']
  if not repeat:continue
  if len(outward)<len(remainder):pending=True;continue
  edits=math.floor((1-a.repeat_identity)*len(repeat)+1e-9);available=len(outward)-len(remainder)
  row=list(range(len(repeat)+1))
  for i,base in enumerate(outward[len(remainder):len(remainder)+len(repeat)+edits],1):
   nxt=[i]
   for j in range(1,len(row)):nxt.append(min(row[j]+1,nxt[j-1]+1,row[j-1]+(base!=repeat[j-1])))
   row=nxt
   if i+edits>=len(repeat) and row[-1]<=edits:return 'array'
  if available<len(repeat)+edits and min(row)<=edits:pending=True
 return 'pending' if pending else 'external'
def hashes(nodes):
 h=14695981039346656037
 for n in nodes:
  h=((h^n)*1099511628211)&((1<<64)-1);yield h
metrics={(r['case_id'],r['side'],int(r['checkpoint_ref_bp'])):r for r in rows(a.metrics)}
paths=rows(oracle/'oracle_paths.tsv');watch=collections.defaultdict(list);longwatch=collections.defaultdict(list)
for i,r in enumerate(paths):
 ns=list(map(int,r['vertices'].split(',')));r['_nodes']=ns;r['_loss']=None
 for depth,h in enumerate(hashes(ns)):
  watch[r['seed_id'],r['side'],depth,h].append((i,ns[depth+1] if depth+1<len(ns) else None))
 longwatch[r['seed_id'],r['side'],len(ns)-1,list(hashes(ns))[-1]].append(i)
for e in csv.DictReader((trace/'production_events.tsv').open(),delimiter='\t'):
 if not e.get('reason'):continue
 if e['event'] in {'terminal_pruned','side_shortlist_pruned'} and int(e['depth'])==0:continue
 key=e['seed_id'],e['side'],int(e['depth']),int(e['path_hash'])
 for i,nxt in watch.get(key,[]):
  r=paths[i]
  if r['_loss'] is not None:continue
  matched=(e['event']=='reject' and nxt is not None and e['target']!='-' and int(e['target'])==nxt)
  matched|=e['event'] in {'beam_pruned','terminal_pruned','side_shortlist_pruned'}
  matched|=e['event']=='stopping_prefix' and nxt is not None and e['reason'] in {'work_limit','length_limit'}
  if matched:r['_loss']=dict(e)
 if e['event'] in {'beam_pruned','terminal_pruned','side_shortlist_pruned'}:
  for x in e.get('prefix_hashes','').split(','):
   if not x:continue
   d,h=map(int,x.split(':'))
   for i in longwatch.get((e['seed_id'],e['side'],d,h),[]):
    if paths[i]['_loss'] is None:paths[i]['_loss']=dict(e)
byseed=collections.defaultdict(list)
for r in paths:byseed[r['seed_id']].append(r)
for c in csv.DictReader((trace/'production_candidates.tsv').open(),delimiter='\t'):
 if not c.get('seed_end') or not c['sequence']:continue
 lo,hi=int(c['seed_begin']),int(c['seed_end']);parts={'left':c['sequence'][:lo][::-1],'right':c['sequence'][hi:]}
 for r in byseed[c['seed_id']]:
  if parts[r['side']].startswith(r['sequence_outward']):
   r.setdefault('_candidate_stages',{})[c['stage']]=c
out=[]
for r in paths:
 m=metrics[r['case_id'],r['side'],int(r['checkpoint_ref_bp'])];e=r['_loss'];stages=r.get('_candidate_stages',{})
 context=boundary_context(r)
 if m['qualifies_90pct']=='1':status='selected_reference_consistent_prefix'
 elif 'eligible' in stages:status='candidate_ranked_below_selected'
 elif e:status=e['event']+':'+e['reason']
 elif 'after_refine' in stages:status='candidate_removed_after_refinement'
 elif 'after_refresh' in stages:status='flank_removed_by_joint_check'
 elif 'before_refresh' in stages:status='flank_removed_by_context_refresh'
 elif 'joined' in stages:status='joined_shortlist_loss'
 elif 'join_rejected' in stages:status='join_rejected'
 elif context=='pending':status='witness_ends_before_source_context_resolved'
 elif context=='array':status='source_context_exclusion'
 else:status='not_located_in_trace'
 if status=='reject:source_array_context':assert context=='array',(r['case_id'],context)
 c=stages.get('eligible',{})
 out.append(dict(seed_id=r['seed_id'],case_id=r['case_id'],side=r['side'],reference_bp=r['checkpoint_ref_bp'],graph_bp=r['query_bp'],status=status,witness_source_context=context,max_visits_in_flank=max(collections.Counter(r['_nodes']).values()),
   loss_at_graph_bp=e['depth'] if e else '',vertex=e['vertex'] if e else '',expected_edge=e['edge'] if e else '',target=e['target'] if e else '',multiplicity=e['count'] if e else '',local_median=e['local_median'] if e else '',gate_penalty=e['penalty'] if e else '',beam_loss=e['beam_loss'] if e else '',eligible_rank=c.get('rank',''),candidate_score=c.get('score',''),selected_flank_bp=m['selected_flank_bp'],selected_prefix_cost=m['weighted_cost']))
with (trace/'verified_path_losses.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
print(collections.Counter(r['status'] for r in out))
for sid in sorted({r['seed_id'] for r in out}):
 rs=[r for r in out if r['seed_id']==sid]
 print(sid,collections.Counter(r['status'] for r in rs))
 for side in ['left','right']:
  choices=[r for r in rs if r['side']==side and r['status']!='selected_reference_consistent_prefix']
  if choices:
   r=max(choices,key=lambda r:int(r['reference_bp']));print(' ',side,r)
