#!/usr/bin/env python3
"""Evaluate a completed native evidence run; the reference never selects paths.

Outputs: updated extensions.tsv and summary.json, REPORT.md, and two plots.
The existing target FASTA and viral contig FASTA are never changed.
"""
import argparse
import csv
import json
import os
from pathlib import Path

import seeded_extension_benchmark as bench


def load_hit_cache(paths, reference):
    """Read frozen development evaluations only after prediction has completed."""
    hits, provenance = {}, []
    for path in paths:
        data = json.loads(path.read_text())
        prediction = path.with_name(path.name.replace('.evaluation.json', '.json'))
        assert data['prediction_sha256'] == bench.digest(prediction)
        assert data['reference'] == str(reference)
        assert data['metrics'] == dict(identity=.9, query_coverage=.5, single_hsp=True)
        hits.update(data['hits'])
        provenance.append(dict(path=str(path), sha256=bench.digest(path)))
    return hits, provenance


def figures(output, rows, selected, hits):
    """Plot measured top-K matches and actual source-masked multiplicities."""
    os.environ.setdefault('MPLCONFIGDIR', '/tmp/mcaat-matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), constrained_layout=True)
    matches = [hits[bench.canonical(r['sequence'])]['passes'] for r in selected]
    axes[0].plot(range(1, len(matches)+1), np.cumsum(matches), color='#27618c', lw=2)
    axes[0].plot([0, 50], [0, 25], '--', color='#999999', lw=1, label='25 / 50 target')
    axes[0].set(xlabel='Saved contigs, in prediction rank order', ylabel='Reference matches',
                xlim=(0, max(50, len(matches))))
    axes[0].legend(frameon=False)
    supporters = np.arange(1, 17)
    axes[1].plot(supporters, np.log2(supporters), color='#236d5e', lw=2, marker='o', ms=3)
    axes[1].set(xlabel='Distinct targets or source arrays supporting a word', ylabel='Word support: log₂(n)',
                xticks=[1, 2, 4, 8, 16], xlim=(1, 16))
    axes[1].set_title('Each target / array contributes once per word')
    fig.savefig(output / 'selection.png', dpi=180)
    plt.close(fig)
    groups = [(label, [r for r in selected if hits[bench.canonical(r['sequence'])]['passes'] == passed])
              for label, passed in [('Matched', True), ('No qualifying match', False)]]
    fig, axes = plt.subplots(2, 2, figsize=(10, 6), sharey='row', constrained_layout=True)
    maximum = max((float(r[side+'_counts'].split(',')[int(p)-1]) for r in selected
                   for side in ('left','right') for p in r[side+'_scored_positions'].split(',') if p), default=1.)
    bins = np.linspace(0, np.log2(1+maximum), 31)
    for j, (label, group) in enumerate(groups):
        counts = []
        for r in group:
            for side in ('left', 'right'):
                values = [float(v) for v in r[side+'_counts'].split(',') if v]
                positions = [int(v) for v in r[side+'_scored_positions'].split(',') if v]
                counts.extend(values[p-1] for p in positions)
        if counts:
            axes[0, j].hist(np.log2(1+np.array(counts)), bins=bins, density=True, color='#27618c', alpha=.85)
        axes[0, j].set(title=f'{label} ({len(group)})', xlabel='log₂(1 + external multiplicity)', ylabel='Position density', xlim=(0,bins[-1]))
        if group:
            r = group[0]  # first prediction-ranked example, not best alignment
            for side, sign, color in [('left', -1, '#27618c'), ('right', 1, '#a95834')]:
                values = [float(v) for v in r[side+'_counts'].split(',') if v]
                positions = [int(v) for v in r[side+'_scored_positions'].split(',') if v]
                axes[1, j].plot([sign*p for p in positions], [values[p-1] for p in positions], color=color, lw=1)
            axes[1, j].set_title(r['query_id'] + ': first ranked example')
        axes[1, j].set(xlabel='Bases from protospacer boundary', ylabel='Observed multiplicity', yscale='symlog')
    fig.savefig(output / 'multiplicities.png', dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--baseline', type=Path, required=True, help='Completed lookahead comparison run')
    parser.add_argument('--reuse-evaluation', action='append', type=Path, default=[])
    parser.add_argument('--cached-hits', action='append', type=Path, default=[])
    parser.add_argument('--blastn', default='blastn')
    parser.add_argument('--threads', type=int, default=2)
    args = parser.parse_args()
    out = args.output.resolve()
    summary = json.loads((out/'summary.json').read_text())
    assert summary['status'] == 'complete', 'Do not evaluate an incomplete run'
    rows = bench.rows(out/'extensions.tsv')
    assert len(rows) == summary['input_starts'] == len({r['query_id'] for r in rows})
    targets = bench.fasta(out/'protospacers.fasta')
    seeds = Path(summary['seed_directory'])
    assert targets == bench.fasta(seeds/'CandidateProtospacers.fasta')
    selected = sorted((r for r in rows if r['reported']=='1' and r['relation_to_representative']=='representative'),
                      key=lambda r: int(r['representative_rank']))
    contigs = bench.fasta(out/'viral_contigs_evidence.fasta')
    assert len(selected) == len(contigs) <= len(targets)
    output_limit = summary.get("max_reported_contigs", 50)
    assert not output_limit or len(contigs) <= output_limit
    assert len({r['target_id'] for r in selected}) == len(selected)
    assert set(contigs.values()) == {r['sequence'] for r in selected}
    for r in rows:
        b, e = int(r['seed_begin']), int(r['seed_end'])
        assert bench.canonical(r['sequence'][b:e]) == targets[r['target_id']]
        flank = summary["max_flank_per_side"]
        assert b <= flank and len(r['sequence'])-e <= flank
        if r['accepted']=='1': assert len(r['sequence']) >= summary['min_contig_length']
    frozen = {name:bench.digest(out/name) for name in ['extensions.tsv','candidate_paths.tsv.gz',
              'protospacers.fasta','viral_contigs_evidence.fasta']}
    cache, provenance = bench.evaluation_cache(args.reuse_evaluation, args.reference)
    extra, extra_provenance = load_hit_cache(args.cached_hits, args.reference)
    cache.update(extra)
    sequences = sorted({bench.canonical(r['sequence']) for r in rows if r['accepted']=='1'})
    missing = {f'C{i+1:06d}':s for i,s in enumerate(sequences) if s not in cache}
    print(f'Evaluating {len(selected)} saved contigs; {len(missing)} uncached candidate sequences', flush=True)
    fresh, command = bench.align(missing, args)
    for key,s in missing.items(): cache[s] = fresh[key]
    target_hits, target_command = bench.align(targets, args, target=True)
    matched = {r['query_id'] for r in selected if cache[bench.canonical(r['sequence'])]['passes']}
    old = [r for r in bench.rows(args.baseline/'extensions.tsv') if r['reported']=='1' and r['relation_to_representative']=='representative']
    previous = {r['query_id'] for r in old if r['reference_pass']=='1'}
    old_summary = json.loads((args.baseline/'summary.json').read_text())['evaluation']
    assert old_summary['identity']==.9 and old_summary['extension_coverage']==.5 and old_summary['single_hsp']
    assert Path(old_summary['reference']).resolve() == args.reference.resolve()
    result = dict(processed_inputs=len(rows), protospacers=len(targets), matching_protospacers=sum(h['passes'] for h in target_hits.values()),
                  saved_contigs=len(selected), matching_contigs=len(matched), baseline_saved=len(old), baseline_matches=len(previous),
                  gained_queries=sorted(matched-previous), lost_queries=sorted(previous-matched), retained_queries=sorted(matched&previous),
                  length_min=min(map(len,contigs.values()), default=0), length_max=max(map(len,contigs.values()), default=0),
                  reference_loci=bench.reference_loci([cache[bench.canonical(r['sequence'])] for r in selected]))
    top50_matched = sum(cache[bench.canonical(r['sequence'])]['passes'] for r in selected[:50])
    result['top50'] = dict(saved=min(50, len(selected)), matching=top50_matched)
    for r in rows:
        h = cache.get(bench.canonical(r['sequence'])) if r['accepted']=='1' else None
        r.update(reference_evaluated=int(h is not None), reference_pass=int(bool(h and h['passes'])),
                 reference_record=h['reference'] if h else '', reference_identity=h['identity'] if h else '',
                 reference_coverage=h['coverage'] if h else '', reference_bitscore=h.get('bitscore',0) if h else '',
                 reference_query_begin=h.get('query_begin',0) if h else '', reference_query_end=h.get('query_end',0) if h else '',
                 reference_begin=h.get('reference_begin',0) if h else '', reference_end=h.get('reference_end',0) if h else '',
                 reference_query_reverse_complemented=int(r['sequence']!=bench.canonical(r['sequence'])))
    with (out/'extensions.tsv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
    summary['evaluation']=dict(reference=str(args.reference), identity=.9, extension_coverage=.5, target_coverage=.85,
        single_hsp=True, reference_used_for_selection=False, alignment_query_orientation='canonical(sequence)', frozen_prediction_sha256=frozen, results=result,
        commands=[command,target_command], target_hits=target_hits, cache_provenance=provenance+extra_provenance)
    figures(out, rows, selected, cache)
    summary['artifact_sha256'] = {name: bench.digest(out/name) for name in
        ['extensions.tsv', 'candidate_paths.tsv.gz', 'protospacers.fasta',
         'viral_contigs_evidence.fasta', 'selection.png', 'multiplicities.png']}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'REPORT.md').write_text(f'''# Native protospacer extension: SRR14017553

Processed **{len(rows)} starts**; preserved **{len(targets)} protospacer sequences**.
Saved **{len(selected)} contigs**, of which **{len(matched)} match** the local reference.
The first {min(50, len(selected))} ranked contigs include **{top50_matched} matches**.
The previous run saved {len(old)} contigs with {len(previous)} matches.

A match requires one HSP with at least **90% identity and 50% full-query coverage**.
Selected contigs are **{result['length_min']}–{result['length_max']} bp**; each flank is at most 750 bp.
The qualifying best hits span {len(result['reference_loci'])} overlapping reference intervals. These counts describe contigs and reference intervals; they do not establish the number of distinct phages.

For inputs represented by a saved matching contig: **{len(matched-previous)} gained, {len(previous-matched)} lost, {len(matched&previous)} retained**.
Exact IDs and all identity/coverage values are in `summary.json` and `extensions.tsv`.
The unchanged targets include {result['matching_protospacers']} reference matches at 90% identity / 85% query coverage.

## Method

1. Extend from the actual saved graph nodes using bounded Huber and lookahead searches; preserve their candidate routes.
2. Add seed-containing windows at the minimum contig length and halfway between that minimum and the maximum flank length (500 and 625 bp here). Keep the complete candidates too.
3. Average the percentile ranks of raw branch support, stop-free coding coherence, and external sequence shared across distinct source arrays. Use raw branch support to break ties.
4. Consider five distinct route families per input, then save nonredundant regions, with at most one contig per target and graph start, subject to the optional output limit. Preserve every input protospacer independently of extension success.
5. Freeze predictions, then align them for evaluation. No reference, read aligner or external assembler chooses paths.

Shorter windows change the amount of retained sequence: this result does not imply recovery of both full 750-base flanks. Scores are ranking heuristics, not calibrated probabilities. The reference was used to compare methods during development after each prediction set was frozen; it is never an input to the predictor. This is development performance on one dataset, not independent validation.
''')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    main()
