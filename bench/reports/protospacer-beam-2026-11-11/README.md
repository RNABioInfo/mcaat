# SRR14017553: equal-depth beam and bounded external repeats

All 81 supplied protospacers remain byte-for-byte unchanged. At >=90% identity and >=50% single-HSP query coverage against the local Shah/EarlyVir database, each beam width produces eight matching extensions, versus six in the saved baseline. All six previous matching inputs are retained. These are reference matches, not a measured false-positive rate.

| Run | Extended input records | Unique extension regions | Matching extensions | Median length (bp) | Reported work-limit records | Wall time (s) |
|---|---:|---:|---:|---:|---:|---:|
| baseline | 45 | 43 | 6 | 432 | 0 | not re-timed |
| N16 | 43 | 41 | 8 | 432 | 0 | 56.9 |
| N32 | 45 | 40 | 8 | 432 | 0 | 61.0 |
| N64 | 45 | 42 | 8 | 601 | 15 | 63.4 |
| N64_100k_states | 45 | 41 | 8 | 774 | 0 | 65.0 |

Use N=32 as the default for this dataset. It extends the same 45 input records as the baseline, adds the two reference matches, and has no reported work-limit stops. N=16 leaves two additional inputs without an extension, although their protospacer records remain intact. N=64 with the original budget hits work limits for 15 inputs. Doubling its per-side budget to 100,000 removes those reported stops and increases median length, but still yields the same eight matching inputs. Longer unvalidated sequences are not evidence of improved accuracy.

The two additional matching inputs are:

- `seed_5f99a193e6244b29`: 377 bp previously; 779 bp at N=32, with 96.662% identity across 100% of the query. Its previously discarded high-count initial branch is now retained long enough to extend. N=64 with the smaller budget stops at 751 bp; the 100,000-state run restores 779 bp.
- `seed_0eb07015d23769af`: 863 bp previously; 995 bp at N=32, with 90.669% identity and 56.382% query coverage. The unaligned remainder is not established by this match.

Among the six existing matches, one extension shortens from 351 to 321 bp while remaining matched over its full query length. The per-input comparison records these changes explicitly. The 81 retained protospacers must not be confused with 81 reconstructed phages: 36 inputs still have no accepted extension at N=32, including five shorter than graph k=23.

## Search changes

- Advance live alternatives at equal base-pair depth before pruning to N. Keep terminal paths separately; rank provisional paths by complete observed multiplicity profiles.
- Permit at most three visits per new vertex. For a vertex already repeated in the immutable seed, permit at most two further visits. Apply the same bound after joining left and right flanks.
- Reject a closed walk when every edge since the previous visit belongs to cycles reported in the source arrays. Include all spacer observations from those arrays, in the correct orientation. Shared nodes alone remain traversable.
- Score all emitted external edge multiplicities using log-count variance, adjacent changes, left/right level agreement and a soft multiplicity-one fraction penalty. The existing sampled local median/MAD gate, complete source-repeat check and adapter check remain active. Scores are heuristics, not probabilities.
- Select at most one final extension per input. Preserve all input IDs, sequences and order independently of extension success.

## Reproducibility

- Baseline commit: `b59de342aebd1381616a0cdd9ed6016985804583`, schema 7, `supported_successors_local_median_v2`. Existing saved outputs were re-evaluated with the same BLAST commands as the new outputs; baseline runtime was not re-measured.
- New binary: schema 8, `equal_depth_beam_bounded_repeats_v3`; source-cycle scope `all_observations_from_source_arrays`.
- Graph: existing SRR14017553, k=23; 49 source arrays / 530 spacer observations. No graph rebuild or mutation.
- Default comparison: 750 bp maximum per flank, 1550 bp maximum total, 50,000 expansions per side, 200,000 per-query cap, three visits per vertex, four local worker threads. Wall times include graph loading and are single runs, not controlled performance benchmarks.
- Local database: 10,021 Shah/EarlyVir sequences. BLAST uses both strands, no DUST/soft masking, E-value search ceiling 1000, word size 7 for queries <=100 bp and 11 for longer queries. Identity is exact `nident/alignment_length`; coverage is one HSP query-coordinate span divided by full query length. HSPs are not pooled.
- No reference sequences, reference alignments or read files were used during path generation or ranking. BLAST runs only after predictions are written.
- Input FASTA SHA256: `9a7886832f257d1d7e4bc69f0cd64aab9354195dbdad47a456e5c6b1376cff12`.
- Graph metadata SHA256: `94ca5ad4422b97de35f7e0d46e165bdfb8efaaabe283166ecd3ea69360e21bce`.
- Matching array report SHA256: `996237f15be63221091946ba066e7e501200857e4b7a679984225b2208a32832`.
- Tested release binary SHA256: `c9fed639bcfc06a8cdae836df4b797e6ba57c8250a2c7e40b7e8dc213c67b048`.

Run `bench/benchmark_retained_beams.py --help` for the reusable benchmark driver. It verifies retained sequences, coordinates, node visit limits and work/length bounds, and writes raw BLAST alignments, commands, per-input results and a summary. The additional budget check uses `--widths 64 --max-states 100000`.

## Validation

- Both C++ test suites pass in Release and under AddressSanitizer/UndefinedBehaviorSanitizer, including two observed real-graph fragments, equal-depth pruning, external repeat exits, forward/backward source-cycle rejection, joined-path bounds, unchanged seeds and work caps.
- Six native SDBG regression scripts pass: extension, retained protospacers, graph filters, source context, adapters and spacer boundaries. Retained-target and boundary fixtures give byte-identical non-parameter outputs with one and 24 threads.
- All three main real-graph runs and the larger-budget run retain the 81 FASTA records unchanged and pass seed-coordinate and per-vertex visit-limit checks. N=32 emits 13 regions containing revisited vertices; all obey the configured cap.
- The 13 distinct qualifying region alignments across baseline/N16/N32/N64 were independently checked against both the predicted query sequences and the local reference FASTA, including coordinates and exact identity counts.
- Input graph-part sizes/mtimes, graph metadata and source/input checksums were checked. The supplied graph is read-only in the extension module.
- One initial synthetic graph-filter fixture build crashed in the pre-existing `buildlib` stage, before extension. A clean rerun completed and all graph-filter assertions passed; the crash cause was not established. The real-graph extension runs did not crash.
