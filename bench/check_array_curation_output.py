#!/usr/bin/env python3
"""Verify final curation preserves all assigned input observations and variants."""
import argparse,csv,json
from collections import Counter,defaultdict
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('before',type=Path);p.add_argument('after',type=Path);a=p.parse_args()
def rows(file):
    with file.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def reverse(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def canonical(s):return min(t[i:]+t[:i] for t in (s,reverse(s)) for i in range(len(t)))
old_systems=rows(a.before/'CRISPR_graph_systems.tsv');old_ids={r['array_id']:str(i+1) for i,r in enumerate(old_systems)}
new_systems=rows(a.after/'CRISPR_graph_systems.tsv');new_ids={r['array_id']:r['candidate_id'] for r in new_systems}
original=rows(a.before/'CRISPR_graph_paths.tsv');accepted=rows(a.after/'CRISPR_graph_paths.tsv');review=rows(a.after/'CRISPR_graph_alternatives.tsv')
fields=[f for f in original[0] if f not in ('array_id','observation_id')] if original else []
def evidence(row,candidate):return (candidate,*(row[f] for f in fields))
expected=Counter(evidence(r,old_ids[r['array_id']]) for r in original)
actual=Counter(evidence(r,new_ids[r['array_id']]) for r in accepted)
actual.update(evidence(r,r['candidate_id']) for r in review)
assert expected==actual,'observations lost or invented during curation'
metrics={r['metric']:r['value'] for r in rows(a.after/'CRISPR_graph_summary.tsv')}
strict=metrics.get('array_filter_mode','strict')=='strict'
owned=defaultdict(set)
for row in accepted:
    key=canonical(row['repeat']+row['spacer'])
    assert row['array_id'] not in owned[key],'circular sequence repeated within candidate'
    if strict:assert not owned[key],'strict output shares a circular sequence'
    owned[key].add(row['array_id'])
    if 'repeat_min_length' in metrics:
        assert int(metrics['repeat_min_length'])<=len(row['repeat'])<=int(metrics['repeat_max_length'])
        assert int(metrics['spacer_min_length'])<=len(row['spacer'])<=int(metrics['spacer_max_length'])
shared={key:owners for key,owners in owned.items() if len(owners)>1}
if (a.after/'CRISPR_graph_shared_cycles.tsv').exists():
    reported={r['canonical_circular_sequence']:set(r['array_ids'].split(',')) for r in rows(a.after/'CRISPR_graph_shared_cycles.tsv')}
    assert reported==shared,'shared evidence not fully disclosed'

variants=rows(a.after/'CRISPR_spacer_variants.tsv')
assert {(r['array_id'],r['variant']) for r in variants}=={(r['array_id'],r['spacer']) for r in accepted}
observed={(r['array_id'],r['spacer']) for r in accepted}
assert all((r['array_id'],r['representative']) in observed for r in variants),'invented representative'
assert len({(r['array_id'],r['variant']) for r in variants})==len(variants),'variant belongs to multiple families'
by_array=defaultdict(list)
for r in accepted:by_array[r['array_id']].append(r)
for r in new_systems:
    paths=by_array[r['array_id']];assert len(paths)==int(r['assigned_cycles'])
    support=sum(int(p['supports_repeat_path']) for p in paths)
    assert support==int(r['support_cycles']) and support+1e-9>=float(r['repeat_support_threshold'])*int(r['candidate_cycles'])
    assert len({p['spacer'] for p in paths})>=2
    if strict:assert len({v['family_id'] for v in variants if v['array_id']==r['array_id']})>=2
# The text report must expose every retained spacer, not only family medoids.
text_observed=set();current=None
for file in a.after.glob('CRISPR_Arrays_*.txt'):
    for line in file.read_text().splitlines():
        if line.startswith('>'):current=line[1:].split()[0]
        elif line.startswith('        ') and '\t' in line:text_observed.add((current,line.split('\t')[-1]))
assert text_observed==observed,'main text report hides observed spacer variants'
redundancy_checks=0
if 'redundant_candidates_removed' in metrics:
    per_candidate=defaultdict(set)
    for row in accepted:per_candidate[new_ids[row['array_id']]].add(canonical(row['repeat']+row['spacer']))
    # Independently verify each removed set has ONE surviving parent; no union
    # of several unrelated candidates may justify removal.
    old_evidence=defaultdict(set)
    for row in original:
        if not int(metrics['repeat_min_length'])<=len(row['repeat'])<=int(metrics['repeat_max_length']):continue
        if not int(metrics['spacer_min_length'])<=len(row['spacer'])<=int(metrics['spacer_max_length']):continue
        if not int(metrics['cycle_min_length'])<=len(row['cycle_vertices'].split(','))<=int(metrics['cycle_max_length']):continue
        old_evidence[old_ids[row['array_id']]].add(canonical(row['repeat']+row['spacer']))
    for decision in rows(a.after/'CRISPR_graph_decisions.tsv'):
        if decision['status']!='redundant_circular_subset':continue
        child=old_evidence[decision['candidate_id']]
        assert child and any(parent in per_candidate and child<=per_candidate[parent]
                             for parent in decision['competing_candidates'].split(',')),'removed set lacks a surviving superset'
        redundancy_checks+=1
    assert redundancy_checks==int(metrics['redundant_candidates_removed'])
    indexed=defaultdict(set)
    for name,sequence_set in per_candidate.items():
        for key in sequence_set:indexed[key].add(name)
    for name,sequence_set in per_candidate.items():
        neighbors=min((indexed[key] for key in sequence_set),key=len)
        assert not any(name!=other and sequence_set<=per_candidate[other] for other in neighbors),'redundant set remains accepted'
print(json.dumps({'original_observations':len(original),'accepted_observations':len(accepted),'review_observations':len(review),'accepted_arrays':len(new_systems),'unique_accepted_circles':len(owned),'missing_or_invented_observations':0,'shared_circular_sequences':len(shared),'duplicate_circles_within_candidate':0,'verified_redundant_removals':redundancy_checks},indent=2))
