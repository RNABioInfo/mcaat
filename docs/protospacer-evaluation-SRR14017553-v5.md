# Source-crossing repair: SRR14017553, 16 November 2026

The complete 398-query workflow now saves **20 reference-matching extensions out of 182**, versus **6 out of 146** before this repair. That is **14 new matching saved sequences and zero lost**. All 146 earlier extension sequences remain exactly unchanged. Nineteen of the 20 matches are at least 750 bp.

The separate test of 81 supplied targets improves from **5/26 to 16/58** matching saved sequences. This repair does **not** restore the historical 44/65 result. Of its 44 matching input IDs, 16 still have a matching prediction, 28 do not, and two other input IDs now match. Twenty of those 28 missing historical matches have a different nonmatching selected extension; eight have no accepted extension. They have not been proven unreachable in the graph.

## What changed

1. Discovery still requires independent context. Traversing the source array cannot create a new detected target.
2. An already detected target, or an explicitly supplied target, first uses the strict extension pass. An accepted extension is kept.
3. When that pass cannot extend the target, allow one flagged source-region passage. Prevent source-node revisits and source-region re-entry; keep external repeats bounded by the existing visit limit.
4. Exclude source-region multiplicities and the overlapping k-mer boundary from downstream scoring. Keep a bounded share of unscored crossing candidates alive until they leave the source region. The existing continuity/flow score is unchanged.
5. Require supplied targets to match the graph exactly over their complete sequence. Shorter-than-k targets use explicit graph padding and verified target coordinates. All supplied records remain in a separate FASTA, even without an extension.
6. Consolidate strict-pass outputs first. A longer recovery candidate cannot replace an already accepted saved sequence. Exact and reverse-complement duplication is still consolidated.

The implementation remains one native C++ class with plain structs. The graph is read-only. There is no reference input, external aligner, or assembler in prediction. BLAST is used only afterward to evaluate frozen outputs.

## Equal-input comparisons

A match requires at least **90% identity and 50% coverage of the complete saved extension in one BLAST HSP**. Identity uses exact identical-base/alignment-length counts; coverage uses query coordinates. All comparisons use the same local curated reference database and BLAST settings.

| Test | Supplied inputs | Saved extensions before → after | Matching saved extensions before → after | Unique matching sequences gained / lost |
|---|---:|---:|---:|---:|
| Full workflow | 530 observations / 398 distinct queries | 146 → 182 | 6 → 20 | +14 / 0 |
| Supplied-target extension | 81 records / 75 canonical queries | 26 → 58 | 5 → 16 | +11 / 0 |

The full workflow retains the same 205 detected queries and 200 unique target intervals; 76 target intervals meet the separate 90% identity/85% coverage target criterion. Its target FASTA is byte-for-byte unchanged. The 36 added extension sequences comprise 14 additional reference matches and 22 that do not meet the agreed criterion. A nonmatch is not proof of a false positive; a matching extension is not a count of distinct complete phage genomes.

Input-level counts differ from unique FASTA counts because several inputs can support one saved region. In the full workflow, matching input queries increase 9 → 26: **17 gained, zero lost**. In the supplied test they increase 5 → 18: **13 gained, zero lost**.

The 81 historical supplied records were previously selected by reference matching. They test extension of a fixed target cohort, not end-to-end discovery accuracy. Their sequences are held fixed in both compared runs. The 398-query workflow starts from the full array report. Neither run uses the reference to choose graph paths. This is a dataset-specific regression test, not held-out validation.

Matching extension lengths are 245–1,536 bp (median 785) in the full run, with 19/20 at least 750 bp. The supplied-target run has lengths 245–1,530 bp (median 785), with 15/16 at least 750 bp. The requested 750 bases per flank remains a maximum; a shorter supported prefix can be saved under the unchanged 200-base minimum.

## Validation and limits

- All 81 supplied records and every emitted full target sequence passed exact retention/coordinate checks, including all five inputs shorter than k=23.
- All 398 query sequences and 530 source observations are identical to the full-workflow baseline.
- All 28 original input files (27 graph files plus the array report) retain their original SHA256 values.
- The tiny native-graph test checks source-only discovery rejection, one permitted passage, source revisit/re-entry stops, exclusion of source counts, exact 19-base target coordinates on both strands, equal minimum/maximum flank lengths, and unique audit path IDs.
- Native self-tests cover preservation during containment consolidation and consistent short-target length qualification. Batch checks and Python compilation passed.
- The full run took 493.34 seconds with two threads; peak observed RSS was 2.332 GiB. The supplied control/recovery runs took 140.00/215.45 seconds. No timeout or guard interruption occurred.
- Search limits were reached for 272/398 full-workflow queries and 74/75 supplied recovery queries. These are bounded predictions, not exhaustive searches or proofs that alternative routes are absent.

The final source differs from the binary used for the three large runs only in the optional audit-ID phase marker. Those runs used `save_search=false`, so the changed line was not executed. The final binary passed the full native fixture with auditing enabled. The frozen prediction source, binary, manifests, command lines, raw BLAST results and comparisons remain in the audit bundle.

## Reproduce the server workflow

```bash
./build/mcaat phage-curate \
  --arrays /vol/d/data/copsac/SRR14017553 \
  --graph /vol/d/data/copsac/SRR14017553/graph/graph \
  --output /vol/d/data/copsac/SRR14017553/protospacer_source_crossing_v5 \
  --source-array-crossing true \
  --threads 24
```

Use a new output directory. This reuses and retains the existing graph. The output is still experimental; the remaining historical losses are recorded rather than presented as solved.

## Audit bundle

`outputs/protospacer-crossing-repair-2026-11-16/` in the working workspace contains:

- `FINAL_SUMMARY.json`: exact counts and provenance.
- `comparison_full398_vs_revision5/`: every full-workflow query and its before/after result.
- `comparison_v4_crossing_change/`: matched control/recovery comparison for all 81 supplied IDs.
- `comparison_v4_true_vs_historical/`: all historical gains and losses.
- `v4_full398_crossing/predictions/`: the complete native outputs.
- `v4_full398_crossing/evaluation/ContiguousPaths_coverage.tsv`: every saved extension's coverage and identity.
- `input_integrity.json`, `final_run_validation.json`, `implementation-v4/`, and `implementation-final/`: integrity checks and frozen implementations.
