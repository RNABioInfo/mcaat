# Recall correction and inclusive cycle lengths

The strict default introduced in `21f8f8d` was too aggressive to use as the
sensitivity default for a metagenome with roughly 750 annotated arrays. Its
small local reference regression did not calibrate the sequence cutoffs or
establish recall on those complete annotations. Reported candidate counts and
annotated arrays must not be equated.

## Current defaults

- Cycle length: **47–100 inclusive**, in distinct nodes on the closed simple
  cycle. The returning start is counted once. The former lower-bound comparison
  excluded the stated minimum; this is fixed. The maximum is checked explicitly.
- Array filter mode: **conservative**. Existing structural checks, 85% repeat
  support, and graph-compatible deduplication remain. The older pairwise
  similarity heuristics and remaining shared-cycle evidence remain flags.
- Native string filtering: **on** in both normal execution and saved-cycle
  processing. The later motif and conserved-spacer rules now reject candidates
  automatically before final exact-set deduplication. Their observed strings and
  reasons remain in alternatives; all scores are recorded. See
  [the native string-filter implementation](string-structure-audit.md) for the
  separate settings and available-reference validation. `--string-filter off`
  reproduces behavior before this addition.
- Both strands of a circle count once within a candidate. Shared evidence
  between candidates is disclosed in `CRISPR_graph_shared_cycles.tsv` and their
  `quality_flags`. Candidates consisting entirely of assigned circles already
  covered by one retained candidate are now removed as redundant reports;
  their assignments remain in review. Partial overlaps retain their unique
  evidence. See [the later exact-set fix](circular-sequence-deduplication.md).
- Every distinct observed spacer appears in `CRISPR_Arrays_N.txt`; `spacers=`
  counts these strings. Similarity families, observed representatives and every
  variant are also listed in `CRISPR_spacer_variants.tsv`. Family grouping is
  not genomic spacer ordering and does not assert sequencing errors.
- `--array-filter-mode strict` opts into the previous sequence rejections and
  exclusive cross-candidate evidence assignment. It is retained for comparison,
  not presented as a calibrated high-confidence classifier.
- The enumeration cap still exits immediately when reached. A merged candidate
  at the cap is flagged in conservative mode and its quadratic similarity
  grouping is skipped; its observed sequences remain available.

`--cycle-min-length`, `--cycle-max-length` and `--array-filter-mode` work in the
normal pipeline and in `postprocess-arrays`. Settings files accept the same
names without `--`. Invalid, fractional and reversed cycle bounds fail before
graph loading. The chosen range and mode are written to output metadata.

## What the supplied reports can establish

The measurements in this section describe commit `2fec173`, before the later
exact-set deduplication fix. They concern the older 906-candidate reports.

The supplied 906 candidates contain assigned cycles of length **46–77**, with
4,014 observations exactly at 77. They cannot reveal cycles of length 78–100
that an earlier search never enumerated. Re-filtering them is not equivalent to
running a fresh search through 100.

Applying the corrected conservative curation and 47–100 bounds to those
previously retained candidates gives:

| Measurement | Result |
|---|---:|
| Retained candidates | 906 |
| Distinct observed spacer strings, summed per candidate | 23,990 |
| Similarity families, summed per candidate | 13,633 |
| Retained cycle observations | 23,991 |
| Duplicate opposite-strand observations removed within candidates | 16 |
| Observations outside the requested length range | 2 |
| Original observations preserved across accepted/review tables | 24,009 / 24,009 |

These 906 are **candidates, not confirmed true positives**. The count was not
chosen to match 750. All uncertain candidates and their overlap flags remain
available for matching against the annotations. Counts of globally unique
spacer strings and circular sequences are reported separately from totals
summed per candidate.

## Validation

The production cycle-enumeration body is exercised at lengths 46, 47, 48, 99,
100 and 101: exactly the requested interval is retained. The cap-exit regression
still passes. The candidate and graph suites test conservative preservation,
strict opt-in rejection, symmetric disclosure of shared evidence, support,
observed variants and within-candidate duplicate removal. Sanitizers pass.

The CLI checks cover defaults, invalid bounds/modes, settings and CLI overrides,
and equality between a normal run and reprocessing its saved cycles.

Replaying 24 local saved-cycle fixtures with the previous strict default versus
the revised conservative mode and length bounds increased exact reference-spacer
matches from **1,031 to 1,105**, with no prior matches lost. Reference sets with at
least two exact spacer matches in one candidate increased from 66 to 67.
Candidate counts changed from 79 to 115. These counts are summed over fixtures,
which include repeated accessions/excerpts; they are not independent genomes,
locus-level precision, or the complete 750-array benchmark. The comparison
changes both curation policy and input length filtering, so it does not isolate
their effects individually.

Fresh searches on the local CP003198.1 and CP025990.1 graphs used the same
conservative policy, minimum 47, and maxima of 77 versus 100. CP003198.1 retained
185 exact reference-spacer matches in both runs. CP025990.1 produced 55 cycles
longer than 77 and increased exact reference-spacer matches from 97 to 98, with
no prior exact matches lost. Candidate counts need not grow monotonically:
additional cycles can change repeat inference and compatible merges. These are
small local comparisons, not a measurement of the complete annotated
metagenome's recall.

## Run from the existing graph

A fresh cycle search is necessary to recover missing longer cycles:

```bash
build/mcaat --graph /path/to/existing/graph/graph \
  --threads 24 --cycle-min-length 47 --cycle-max-length 100 \
  --array-filter-mode conservative --repeat-support 0.85 \
  --autoclean false --output-folder results/cycles47_100_conservative
```

This reuses the graph and reruns cycle search. Use a fresh output directory.
Do not rebuild the graph or replay the old capped cycle file expecting it to
contain newly discovered lengths.

When those new cycles have been saved, strict/conservative comparisons can use
`postprocess-arrays` on that **same** input to isolate filtering losses from
search losses. Compare annotated loci and all observed spacer sequences—not
just the number of candidates or displayed similarity families.
