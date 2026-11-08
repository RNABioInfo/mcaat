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
owned={}
for row in accepted:
    key=canonical(row['repeat']+row['spacer'])
    assert key not in owned,'accepted circular sequence appears more than once'
    owned[key]=row['array_id']
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
    assert len({v['family_id'] for v in variants if v['array_id']==r['array_id']})>=2
print(json.dumps({'original_observations':len(original),'accepted_observations':len(accepted),'review_observations':len(review),'accepted_arrays':len(new_systems),'unique_accepted_circles':len(owned),'missing_or_invented_observations':0,'duplicate_accepted_circles':0},indent=2))
