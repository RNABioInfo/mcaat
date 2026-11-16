# SRR14017553 evaluation: v8 extension coverage filter

V7 could save a selected extension whose multiplicity changes substantially across longer stretches. V8 selects the same route first, then checks complete blocks of downstream counts. If block medians differ by more than `--max-coverage-fold 2`, the extension is withheld. Its protospacer, input associations and discovery evidence are retained. No alternative route is promoted after this failure.

## Full discovery run

Inputs: the same existing k=23 graph, 49 arrays, 530 spacer observations and 398 distinct corrected spacers. Comparison is against published v7 (`c69f35b8`). Identical sequences and reverse complements share one saved sequence; input associations remain available.

| Saved output | V7 | V8 |
|---|---:|---:|
| Candidate protospacer regions | 159 | 159 |
| Reference-matching protospacer regions | 53 | 53 |
| Extensions | 142 | 106 |
| Reference-matching extensions | 50 | 49 |
| Nonmatching extensions | 92 | 57 |

The filter removes **35 nonmatching extensions and 1 matching extension**. The reference-matching fraction changes from 35.2% to 46.2%. This is a reference-match fraction on this development dataset, not a biological false-positive rate.

The lost extension belongs to Q000191 and Q000288. It previously aligned to OTU_2363 at 100% identity and 82.9% query coverage. Its complete block medians are 7, 29, 26, 22 and 36 (maximum/minimum = 5.14), so it fails the twofold rule. Matching spacer queries change from 62 to 60, with no gains. Both queries share the one lost extension.

All 106 accepted extension strings were already present before filtering. Accepted routes retain their sequence and length. All 40 withheld query routes were checked against the earlier graph-node paths and counts; the target FASTA and discovery-context table are byte-identical to v7.

## Separate provided-target preservation check

The supplied FASTA contains 81 records and 75 distinct sequences after reverse-complement grouping. All 81 records, IDs and sequences remain unchanged. This reference-selected input set tests preservation, not discovery sensitivity.

| Saved extensions | V7 | V8 |
|---|---:|---:|
| Total | 69 | 54 |
| Reference matches | 43 | 43 |
| Nonmatches | 26 | 11 |

All 43 previously matching extension strings and all 48 input IDs associated with a matching extension remain. No previously matching input is lost.

## Shuffled-spacer control

The same 398 shuffled queries produce 8 candidate targets before and 8 after. Accepted extensions change from 8 to 6; reference-matching extensions change from 1 to 1. This is a background control on the same graph, not an estimate of the real sample's false-discovery rate.

## Method and evaluation

For configured read length 150 and k=23, each block contains 127 scored edge counts. Source-array passage counts are excluded. At least two complete blocks are required; a trailing incomplete block is ignored. The twofold value is an empirical, adjustable setting, not a probability or a published significance threshold. Genuine abundance variation can fail it. `--max-coverage-fold 0` disables this final check.

Prediction receives no reference. Frozen outputs are evaluated afterward against `all_curated_db`: targets require ≥90% identity and ≥85% query coverage; extensions require ≥90% identity and ≥50% coverage of the whole saved sequence, in one HSP. The evaluation does not require that HSP to overlap the target. The parameters and output files are unchanged during evaluation.

An earlier experiment applied a fourfold check while traversing and allowed alternatives to replace rejected routes. It saved 142 extensions with 51 reference matches, compared with v7's 142/50. It was not adopted: the final filter makes the sensitivity cost explicit and substantially reduces nonmatching calls.

## Verification and limits

- The weighted-graph regression verifies rejection of the unstable winner, target preservation, retained node/count evidence and absence of fallback to a weaker route.
- Source-array passage, repeat re-entry, adapter branches, pure-spacer membership, short-target boundaries, reverse complements and exact flank lengths pass their regression tests.
- Batch integration passes with the native twofold default and a valid zero-extension output. Graph construction tests use the native `buildlib`/`read2sdbg` stages used by the batch runner.
- The initial tiny fixture exposed a crash in the legacy combined command's subsequent array-cycle search. That separate issue remains unresolved; these tests cover graph construction and the protospacer workflow, not the legacy array-cycle stage.
- All 28 original graph/input component hashes are unchanged. All 398 full-run queries and all 81 supplied records are accounted for.
- No external validation dataset has been tested. Missing reference matches alone do not establish that a prediction is false.

Local wall time was 536.2 s for v7 and 914.1 s for v8, using two threads. These runs do not establish a speed improvement. Full-run mapping work (28,932,847) and extension work (54,593,096) are identical; the complete filter/assembly stage consumed 1.05 summed worker-seconds in v8. Timing tables remain available per query and stage.

`ContiguousPaths.fasta` contains accepted extensions. `withheld_extensions.tsv` contains the excluded selected routes and their block statistics for inspection; these rows are not accepted extensions. See the [workflow and function reference](protospacer-extension.md) for the full method and output definitions.
