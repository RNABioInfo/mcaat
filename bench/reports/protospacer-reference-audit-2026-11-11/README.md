# SRR14017553: all-81 reference and graph audit

All 81 supplied protospacers were checked against the local reference database
and the existing k=23 graph. The earlier claim that almost all 81 reference
extensions could be recovered was unsupported: k-mer presence and endpoint
reachability did not establish the intervening sequence. This audit checks the
actual path from each original seed and records the exact missing edge or search
decision. No reference sequence enters production path generation or ranking.

The final code preserves the 81 FASTA records byte-for-byte and retains all eight
previous extension matches at N=16, 32 and 64. **The number of matching extensions
remains eight; this change does not increase it.** At N=32, three matching
predictions become longer, one becomes shorter, and four keep their length.
The shorter prediction and remaining ambiguous decisions are reported below.

## What the graph actually supports

The audit covers 1,041 reference locations: all 1,017 qualifying short-seed HSPs,
plus 24 additional seed-spanning HSP locations from previously matching extensions.
A fresh BLAST of all 81 inputs reproduced the 1,017 seed HSPs exactly.

- **8 inputs:** a connected flank matching at least 100 reference bases at >=90%
  aligned identity was explicitly found and its graph path was spelled and checked.
- **68 inputs:** no >=90%-identity path through the first 24 or first 100 reference
  bases exists at any tested location with that much known reference flank.
  These negatives have completed cost-bound certificates, not merely a timeout.
- **5 inputs:** shorter than k=23 (17, 17, 19, 21 and 22 bases). The current exact-seed
  entry cannot place them as a full graph node. They remain in the output unchanged.

Exact-reference traversal is complete for every available flank, independent of
the approximate-search cap. Five inputs have an exact 750-base flank; the other
three positive inputs reach exact maxima of 42, 120 and 175 bases. The remaining
mapped inputs reach at most 14 exact flank bases at the tested locations.
Three inputs have exact 750-base paths on both sides of the same reference
location, but source-array exclusions can still make one side ineligible.

The approximate oracle is a diagnostic Dijkstra search over graph-node/reference-
position pairs. At reference length L, matches cost 0, substitutions/deletions 10
and insertions 9; cost <= L is equivalent to >=90% aligned identity. Every reported
positive includes a connected node path, actual emitted bases and edge counts.
It found 810 location/side/checkpoint witnesses belonging to eight inputs; these
are **not 810 systems**. There are verified >=250-base flanks for seven inputs
and >=500/750-base flanks for five.

The oracle permits indels and unrestricted graph revisits, so a witness is not
necessarily admissible under production's 750-base flank cap, three-visit cap or
source-array rules. Those distinctions are recorded in the outcome tables.
`edit_cost` in the oracle table is the unit edit count of the witness; the queue's
weighted identity bound is described above.

The million-state cap leaves some longer approximate searches unresolved for
73 inputs, including inputs with shorter verified flanks. **The verified counts
are not an upper bound on every possible approximate reconstruction.** The
anchored flank test also differs from the production benchmark's local BLAST
criterion: a local match elsewhere in a predicted sequence can pass even when
its first 100 bases from the seed do not match that reference.

## Where paths were lost and what changed

The logging-only baseline reproduces all 81 saved production outcomes exactly.
All 810 verified witnesses now have an explicit disposition in
[before_path_outcomes.tsv](before_path_outcomes.tsv) and
[after_path_outcomes.tsv](after_path_outcomes.tsv). No verified witness is left
with an unexplained `not_located_in_trace` result. Outcomes include selected
reference-consistent sequence, a specific rejected edge, beam rank, candidate
ranking, source-context exclusion, or a short witness that ends before the
source context can be resolved. Counts across reference locations are not counts
of independent biological systems.

1. **Stale multiplicity baseline.** For `seed_f00881e5a59e1c4e`, the old search
   rejected a verified continuation after 87 outward bases, at vertex 114005295,
   edge 114005295, target 287048003. Its count was 67; the separated-sample median
   was still 23, producing a loss of 27.636 above the limit of 20. The real path
   had already climbed through counts around 37 to 71. The gate now retains both
   the immediate all-edge window and the separated-observation window, accepting
   support from either. A rejected sampled edge updates neither. The initial
   flank-comparison baseline remains separate.
2. **A short, flatter prefix winning over a supported continuation.** For
   `seed_815985dbbb929d9a` and `seed_c35af5b66b335ecd`, complete candidates survived
   but adapter-trimmed prefixes received slightly better stability scores.
   Supplied-target selection now prefers more sequence within the existing
   fixed score margin, provided within-flank instability also stays within that
   margin. It keeps the support penalty, so an absent flank cannot be rewarded.
   Agreement between mean flank depths cannot hide a large increase in variation.
   This is a graph-only ranking heuristic, not a claim that longer is always correct.
3. **Uncertainty remains explicit.** Source-array repeat exclusions, source-cycle
   rejection, adapter checks and bounded external repeats remain active. Beam
   pruning remains bounded. For example, one previously selected
   `seed_0eb07015d23769af` path was pruned at depth 374 and rank 33 in a diagnostic
   trial. A valid graph path is not guaranteed to win a finite multiplicity beam.
   Multiple reference locations and alternative flanks can be incompatible with
   one another; the tables do not pretend a single chosen path recovers all of them.

The reference-consistent path's multiplicity is not uniformly flat: the verified
250-base `hit_0672` flank includes sustained changes between approximately
23, 71, 55, 13 and 67. Counts alone cannot prove which branches belong together
when the pooled graph has merged shared sequence.

## Final comparison

An extension match means >=90% identity and >=50% single-HSP query-coordinate
coverage against the local Shah/EarlyVir database. HSPs are not pooled. The
81 frozen short-input matches use >=90% identity and >=85% coverage; the two
counts must not be conflated.

| Run | Retained inputs | Extended inputs | Unique regions | Matching inputs | Previous matches lost | Work-limit inputs |
|---|---:|---:|---:|---:|---:|---:|
| baseline | 81 | 45 | 40 | 8 | 0 | 0 |
| N32 | 81 | 45 | 40 | 8 | 0 | 0 |
| N16 | 81 | 43 | 40 | 8 | 0 | 0 |
| N64 | 81 | 45 | 42 | 8 | 0 | 15 |

N=32 remains a practical default for this dataset: it extends 45 inputs without
work-limit stops. N=64 has 15 work-limit inputs with the same 50,000-state cap.
No beam width produces additional matching input IDs. Single-run timings include
graph loading and are saved in [beam_comparison.tsv](beam_comparison.tsv).

| Input | Before bp | After bp, N32 | Identity % | Query coverage % | Best reference |
|---|---:|---:|---:|---:|---|
| seed_541b71503abcacf8 | 1530 | 1530 | 100.000 | 100.00 | OTU_1556 |
| seed_815985dbbb929d9a | 1240 | 1526 | 100.000 | 52.03 | OTU_2046 |
| seed_0e3fc52adc6e1e7e | 1528 | 1528 | 100.000 | 100.00 | OTU_1556 |
| seed_f00881e5a59e1c4e | 117 | 780 | 100.000 | 89.23 | OTU_2363 |
| seed_bdf0cbd0c80d2f18 | 321 | 321 | 96.262 | 100.00 | OTU_585 |
| seed_c35af5b66b335ecd | 1242 | 1526 | 99.372 | 52.16 | OTU_2046 |
| seed_0eb07015d23769af | 995 | 421 | 99.035 | 73.87 | OTU_815 |
| seed_5f99a193e6244b29 | 779 | 779 | 96.662 | 100.00 | OTU_2363 |

`seed_0eb07015d23769af` shortens from 995 to 421 bases. Its best match changes from
90.669% identity over 561 query bases to 99.035% over 311 bases. Thus the match is
retained, but the aligned span is shorter. The 1,526-base predictions have only
about 52% single-HSP coverage; the added length is not evidence that every added
base has been independently validated. Full per-input results are in
[all81_audit.tsv](all81_audit.tsv).

## Inspect the evidence

- [all81_audit.tsv](all81_audit.tsv): every input, exact/approximate flank lengths,
  unresolved-search counts, before/after predictions and match status.
  Left/right maxima can come from different reference locations.
- [exact_missing_edge_examples.tsv](exact_missing_edge_examples.tsv): original
  graph vertex, its 23-mer, required next 23-mer and actually available successors.
- [oracle_paths.tsv.gz](oracle_paths.tsv.gz): verified emitted sequences, graph
  node IDs, edge counts and alignment operations for every witness.
- [verified_seed_anchored_paths.fasta.gz](verified_seed_anchored_paths.fasta.gz):
  original seed plus one verified graph flank, in input-forward orientation.
  These are diagnostic witnesses, not production predictions or invented joins.
- [exact_native_routes.tsv.gz](exact_native_routes.tsv.gz) and
  [oracle_routes.tsv.gz](oracle_routes.tsv.gz): all reference-location outcomes,
  including certified negatives and explicitly unresolved searches.
- The reproducible reference windows and their checksums are in
  [the audit fixture](../../fixtures/protospacer-reference-audit/README.md).

## Validation and rerun

Release tests and AddressSanitizer/UndefinedBehaviorSanitizer tests pass, including
independent string-DP checks of the oracle, exhaustive short branching graphs,
work-cap certificates, gradual depth drift, local valleys, repeat limits,
source exclusions and near-tie selection. Native SDBG fixtures pass with one and
24 threads. On the complete supplied graph, four and 24 threads produce identical
non-parameter output files. Both logging-only traces reproduce their production
selections exactly. The graph files and all 81 input FASTA bytes remain unchanged.

Baseline commit: `4ae6c31b3b8cc935e780a8abf3f1e9fd5cc76eb6`.
Final binary SHA256: `91ad9a6297e17cdae5a685253e6c2f955a2c4edde45d7eb02f2a9d57237ad2b5`.
The production algorithm is schema 9, `equal_depth_beam_local_window_v4`.
Validation records are in [validation.json](validation.json) and
[native_fixture_validation.json](native_fixture_validation.json).

On the server, use the binary from the updated checkout:

```sh
cd /vol/d/mcaat-protospacer/mcaat-protospacer-full
git switch codex/protospacer-extension
git pull --ff-only origin codex/protospacer-extension
cmake --build build --parallel 24
./build/mcaat phage-curate \
  --arrays /vol/d/data/copsac/SRR14017553 \
  --graph /vol/d/data/copsac/SRR14017553/graph/graph \
  --protospacers bench/fixtures/protospacer-baseline-81/Protospacers.fasta \
  --output /vol/d/data/copsac/SRR14017553/protospacer_audited_extensions \
  --threads 24 --beam-width 32 --max-node-visits 3 \
  --max-flank-length 750 --max-length 1550 --max-states 50000
```

This reuses the existing graph and preserves all supplied protospacers. It does
not rebuild or delete the graph. References are used only in the separate audit
and post-prediction evaluation utilities.
