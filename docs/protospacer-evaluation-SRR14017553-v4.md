# SRR14017553: placement and coverage experiments

The retained change tries another qualifying spacer placement when the better-aligned placement cannot be extended. On the complete dataset it saved **146 unique extensions**, of which **6** met the reference criterion. Compared with the preceding local v3 run, matching input-query results gained **1** and lost **0**.

**No additional distinct matching extension was recovered.** The matched set stayed at the same six sequences while total saved extensions increased from 130 to 146. The extension-quality goal remains unresolved.

This remains an experimental extension method. A code fix and a completed benchmark do not establish reliable phage prediction. The broader coverage-ranking changes were rejected because they shortened previously matching regions; one trial increased the match count while losing most of the previously recovered sequence.

## Inputs and evaluation

- One unchanged k=23 SRR14017553 graph: 49 arrays, 530 spacer observations, 398 distinct corrected queries. No historical 81-input whitelist.
- Native prediction used only arrays and graph, with two threads, beam 64, maximum 750 added bases per side, and minimum accepted flank 200 bases. Every source observation remains recorded.
- Predictions were frozen before BLAST against the existing 10,021-sequence, 250,720,859-base curated database.
- Targets: at least 90% identity and 85% query coverage. Extensions: at least 90% identity and 50% query coverage in one HSP. The thresholds were not lowered.
- References informed development and evaluation, but never selected prediction branches. This is a development benchmark, not an independent validation. SRR14017554 and SRR14017569 were not run.

## Results

| Implementation | Target intervals / matches | Extensions / matches | Matching sequences ≥750 bp | Median matching length | Sum of matching HSP query spans |
|---|---:|---:|---:|---:|---:|
| Published v2 baseline | 144 / 48 | 58 / 5 | 1 | 442 bp | 3030 bp |
| Previous local v3 | 203 / 77 | 130 / 6 | 5 | 1157.5 bp | 6411 bp |
| Average support + distance buckets | 238 / 90 | 106 / 7 | 0 | 258 bp | 1920 bp |
| Average support + previous traversal | 230 / 89 | 159 / 5 | 0 | 267 bp | 1316 bp |
| Placement continuation only (v4) | 200 / 76 | 146 / 6 | 5 | 1157.5 bp | 6411 bp |

Each FASTA sequence is counted once after exact, reverse-complement and containment consolidation. Matching input queries are a separate count because several inputs can support the same saved sequence. HSP spans can overlap in the reference; their sum is not unique reference-genome coverage.

## Query gains, losses and preservation

- Gained matching query IDs: Q000305.
- Lost matching query IDs: none.
- Previous matching queries whose complete saved sequence is still present on either strand: 8/8.
- Previous matched sequences shortened, changed or missing: none.

The broader scoring trial retained eight matching query IDs but preserved the complete old sequence for only one of them. This is why a match-count increase alone was not accepted as improvement.

## Retained implementation

Mapping tiers are ordered by query coverage and identity. The earlier best tier retains its original work allocation. A detected target without an accepted extension no longer ends the search: lower tiers receive the remaining budget. Selection prefers an extended placement, then alignment quality. Equally aligned different targets in the same outcome category remain unresolved.

The preceding unitig flow and continuity ranking, source-array checks, adapter/complexity checks, bounded external repeats, ambiguity handling and one-primary-result policy are retained. The experimental average-log-support score and distance-bucket search are absent from the retained code. No additional prediction tool or dependency was added.

## Targeted failures and controls

The six traced failure queries were tested separately. Placement continuation recovered a 245-base matching fragment for Q000305 (99.59% identity, 100% coverage). That sequence already appeared in the previous full run through another input. It is additional query support, not by itself a new matching phage region. The other five targeted cases did not yield accepted reference-matching extensions. The known alternative graph walks have not all been recovered or correctly selected.

The fixed 398 dinucleotide-preserving shuffled controls yielded 7 targets / 5 extensions before, and 7 targets / 5 extensions after. These are background controls, not certified biological negatives or a calibrated false-discovery rate.

## Every matching saved extension

| ID | Length | Reference | Identity | Query coverage |
|---|---:|---|---:|---:|
| E000021 | 786 bp | OTU_2041 | 98.09% | 100.00% |
| E000025 | 1529 bp | OTU_1556 | 100.00% | 100.00% |
| E000063 | 1536 bp | OTU_2191 | 100.00% | 100.00% |
| E000095 | 245 bp | OTU_815 | 99.59% | 100.00% |
| E000106 | 786 bp | OTU_2041 | 98.22% | 99.87% |
| E000145 | 1530 bp | OTU_1556 | 100.00% | 100.00% |

## Validation and practical limits

- Native wall time: 385.8 seconds with two threads; peak sampled RSS: 2.46 GiB. Queue waits are excluded. The memory guard did not stop the run.
- 268/398 queries reached a mapping, beam, visit or work limit. Completion is not an exhaustive search guarantee.
- Native graph-edge/sequence-spelling checks and the independent coordinate, identity, count and association validator passed.
- All 28 original graph/array files retain their original SHA-256 checksums. The workspace binary and source match the frozen tested implementation.
- Build, native self-checks, batch handling checks, Python compilation and whitespace checks are recorded separately in the local validation log.
- Reference-match rate is not biological precision: an unmatched sequence is not automatically false, and these fragments are not complete phage genomes.

See the [workflow and every-function reference](protospacer-extension.md) and [sequencing evidence](protospacer-sequencing-evidence.md).
