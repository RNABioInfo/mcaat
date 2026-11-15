# Native workflow v2: SRR14017553 evaluation

The connected-path mapper improves detection over the preceding local prototype.
The extension method is **not yet a reliable replacement**: only 5 of 58 saved
extensions meet the specified reference criterion. Four of those five are shorter
than 750 bp. This is an experimental implementation, not a demonstrated solution
to the extension-quality problem.

## Input and procedure

- One existing k=23 graph: SRR14017553. No graph rebuilding or read downloads.
- The complete array report: 49 arrays, 530 spacer observations, 398 distinct
  corrected queries. All observations remain in `sources.tsv`; no 81-input
  whitelist or supplied-target preservation pass was used.
- Array-report SHA-256:
  `996237f15be63221091946ba066e7e501200857e4b7a679984225b2208a32832`.
- Two CPU threads; beam 64; maximum 750 added bases per side; minimum accepted
  flank 200 bases. Other parameters are the documented v2 defaults.
- Prediction used only the arrays and graph. Predictions were frozen before
  reference evaluation. No reference alignment selected a target or branch.
- Afterwards, local BLAST compared all outputs and all 398 corrected queries with
  the existing curated database: 10,021 sequences, 250,720,859 bases.
- Target criterion: at least 90% identity and 85% query coverage. Extension criterion:
  at least 90% identity and 50% query coverage. Each requires one qualifying HSP;
  disjoint hits are not added together. Exact identity uses `nident / alignment length`.

The comparison below uses the preceding **full-input local v1 prototype**, which had
failed to produce a reference-matching accepted extension. It is not a comparison
with the historical 81-input experiment or an earlier server's intermediate output.

## Results

| Measurement | Earlier local v1 | Final native v2 |
|---|---:|---:|
| Original spacer observations | 530 | 530 |
| Distinct corrected queries | 398 | 398 |
| Unique detected target intervals | 56 | 144 |
| Reference-matching target intervals | 13 | 48 |
| Unique accepted extensions | 1 | 58 |
| Reference-matching extensions | 0 | 5 |
| Queries reaching a search limit | 398 | 268 |
| Total native wall time, two threads | 95.2 s | 325.8 s |

The final 144 target intervals are associated with 148 distinct input queries.
The 58 extensions are associated with 60 queries. These are different counting
units: several inputs can support the same saved sequence.

Comparing the **same query IDs and sequences**:

- Matching target results: 39 query gains, 1 loss, 12 retained; 13 → 51 queries.
  The lost query, `Q000393`, is dominated by a poly-G run and now fails the explicit
  complexity check. It is still recorded as an input and rejection.
- Matching extension results: 7 query gains, 0 losses; these seven queries support
  the five distinct matching extensions. This gain is against v1's zero, not against
  any historical set of previously reported successful extensions.
- Of the 140 input queries that themselves match the database, 50 have a detected
  target and 12 have an accepted extension. An input reference match does not establish
  that its corresponding flanks exist in this sample graph.
- Six queries have a reference-matching detected target but their accepted extension
  fails the extension criterion: `Q000007`, `Q000111`, `Q000206`, `Q000235`, `Q000331`,
  `Q000382`. The present test does not establish which alternative graph paths could
  replace them.

## The five matching extensions

| Extension | Length | Reference | Identity | Query coverage |
|---|---:|---|---:|---:|
| E000001 | 256 bp | OTU_1556 | 100% | 100% |
| E000009 | 442 bp | OTU_2041 | 97.96% | 100% |
| E000012 | 550 bp | OTU_2041 | 98.36% | 99.82% |
| E000028 | 1,536 bp | OTU_2191 | 100% | 100% |
| E000035 | 247 bp | OTU_1556 | 100% | 100% |

All 58 saved extensions range from 226 to 1,536 bp, with median 776 bp.
Thirty of the 60 extended queries reach 750 added bases on at least one side;
three reach 750 on both sides. The remainder stop earlier. None of these counts
means complete phage genomes.

## What this experiment establishes

The new mapper finds more reference-matching target intervals, but takes longer.
The widest-path score and ambiguity filtering do not yet select enough useful
extensions to justify treating this as an accurate phage-extension method.
The 5/58 figure is a reference-match rate (8.6%), not a measured biological precision:
an unmatched sequence is not automatically false, and a matched sequence does not
prove an entire graph route was uniquely resolvable.

The final optimization improved seed searching and removed repeated alignment work
without changing thresholds. The first v2 pass produced 139 targets / 55 extensions,
with 47 / 5 reference matches. The final pass produced 144 / 58 with 48 / 5 matches.
It therefore did not improve the number of matching extensions.

Search limits remain material: 137 of the 148 detected-query results have a mapping
or extension limit flag. A selected limited result must not be called globally unique.
No exhaustive reference-to-graph recoverability audit was performed in this experiment.
No claim is made that unmatched continuations are absent from the graph.

## Verification and measurements

- Native release build and arithmetic/context/ambiguity checks passed.
- 192 exact/mutated seed searches recovered sampled known graph labels.
- Native output construction verified the selected graph edges and sequence spelling.
- Independent Python checks passed for every saved association, target identity,
  query coverage, fragment coordinate, count and reverse-complement deduplication.
- All 28 original array/graph files retained their original SHA-256 checksums.
- Prediction binary and source matched the frozen implementation snapshot.
- Peak sampled native RSS: approximately 2.52 GiB; no memory-pressure stop.
- Native wall times: graph load 70.05 s; index checks 24.51 s; spacer preparation
  0.005 s; detection/extension/filtering 231.07 s; consolidation 0.017 s; writing 0.094 s.
- Batch checks passed for authoritative input selection, empty valid output, duplicate
  ID rejection, successful logs and failed-stage accounting. The other two requested
  datasets were listed in a dry run only; their graphs were not run locally.

## Reproduce

See [workflow options and function reference](protospacer-extension.md) for prediction.
To evaluate a frozen result directory separately:

```bash
python3 tests/evaluate_protospacer_workflow.py \
  --srr SRR14017553 \
  --results /path/to/native_flow_results \
  --output /path/to/new_evaluation_directory \
  --database /path/to/all_curated_db \
  --threads 2
```

The evaluator saves exact commands, runtime, all hit tables, every sequence's best
coverage, validation and checksums. BLAST is optional evaluation only and is not
part of the native prediction workflow.
