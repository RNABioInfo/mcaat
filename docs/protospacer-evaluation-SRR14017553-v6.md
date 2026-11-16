# MCAAT v6: completed SRR14017553 evaluation

The completed ordinary-discovery run on **398 source queries** produces **172 unique protospacer targets** and **154 unique extensions**. **52 targets** and **49 extensions** match the reference under their respective criteria. Those matching extensions serve **61 source queries**. Relative to commit `094da53`, 48 queries gain a matching extension and 13 lose one: a net gain of 35 matching-query outcomes.

The separate fixed-input test retains all **81 supplied protospacer records exactly**, produces **69 unique extensions**, and finds matches for **43 of those sequences**, serving **48 original input records**. Shuffled controls produce eight targets and eight extensions from 398 control queries. These are distinct experiments and denominators; they must not be combined into a biological precision estimate. All three runs and their reference evaluations completed. Increasing the pending context pool from 256 to 1,024 produced no matching-input gains or losses and took 73.5% longer; the default remains 256.

## Ordinary discovery on all 398 queries

This test starts from the same 530 source-spacer observations and 398 canonical queries as the previous workflow. Input selection did not use the reference, and the saved comparison confirms unchanged observations and corrected boundaries.

| Quantity | Previous workflow (`094da53`) | Current v6 |
|---|---:|---:|
| Canonical source queries | 398 | 398 |
| Queries with detected targets | 205 | 177 |
| Unique target sequences | 200 | 172 |
| Unique targets matching at 90% identity / 85% coverage | 76 | 52 |
| Queries with accepted extensions | 196 | 173 |
| Unique saved extensions | 182 | 154 |
| Unique extensions matching at 90% identity / 50% coverage | 20 | 49 |
| Queries associated with a matching extension | 26 | 61 |

The 49 matching v6 extensions are **779–792 bp** long; their median is 780 bp. Across all 154 outputs, 151 are 750–799 bp and three are 314–499 bp. The increase in matching extensions is therefore not produced by shortening every prediction to a small matching fragment. The 49/154 reference-match fraction is 31.8%; this is not a biological precision estimate.

The native run took **472.09 seconds with two threads**, including graph loading and excluding the subsequent BLAST evaluation. It completed successfully and passed output validation. Of the 398 queries, 201 recorded a search-limit flag; a completed bounded run is not an exhaustive graph search.

### Exactly what was gained and lost

Of the 26 previously matching extension queries, 13 still have a matching extension. Another **48 gain a match**, while **13 lose a match**, yielding 61 current matching queries.

* **Ten losses occur before extension recovery.** Their formerly eligible strict detection contexts encounter an adapter core within one substitution under the corrected artifact check. Those saved contexts no longer supply sufficient independent target context, so the crossing pass is not entered. This explains rejection of the previously used contexts; it does not prove that every biological alternative is absent from the graph.
* **Three queries remain detected but lose the extension match:** `Q000010`, `Q000074` and `Q000165`. The new one-flank selection retains a different flank from the matching alternative. This is a selection failure, not evidence that their targets disappeared.

Target matching also decreases: 78 matching target-query outcomes become 56, with five gains and 27 losses. Twenty-five losses are associated with rejection of the old strict contexts; two retain a different selected target that does not meet the reference criterion. The extension improvement does not mean every detection measure improved.

## Shuffled-query controls

The control set contains one distinct dinucleotide-preserving shuffle for each of the 398 canonical queries, retaining sequence lengths and source associations. All 398 original/control pairs differ. The audit verifies **zero boundary changes** in the controls: no sequence was moved back toward its original input by boundary correction.

| Control outcome | Result |
|---|---:|
| Canonical shuffled queries | 398 |
| Detected targets | 8 |
| Accepted extensions | 8 |
| Unique targets / unique extensions | 8 / 8 |
| Target reference matches | 1 |
| Extension reference matches | 1 |
| Boundary-corrected observations | 0 |

The target and extension matches come from **different controls**: `Q000044` supplies the matching target, while `Q000392` supplies the matching extension. They are not one completely matching target-plus-extension prediction. The background detection fraction is 8/398, or 2.01%. These shuffled sequences are not certified biological negatives, so this fraction is **not a biological false-discovery rate**. The control run took 199.73 seconds before reference evaluation.

## Completed fixed-input test

The 81 unchanged input records form 75 distinct sequences after reverse-complement grouping. There are 126 source-context associations; these are not extra inputs. The fixed input set was originally selected using reference matches. MCAAT's prediction run did not load the reference or use it to choose graph paths.

| Quantity | Result |
|---|---:|
| Supplied input records retained unchanged | 81/81 |
| Complete supplied sequence verified in every output association | Passed |
| Distinct canonical queries | 75 |
| Queries with accepted extensions | 75/75 |
| Original input records with accepted extensions | 81/81 |
| Unique saved extensions after consolidation | 69 |
| Unique extensions matching the reference | 43/69 |
| Original input records associated with a matching extension | 48/81 |
| Unique extension length range | 769–792 bp |
| Median unique extension length | 784 bp |
| Queries encountering a recorded search limit | 75/75 |

All 69 outputs contain the requested 750-base selected flank, with the original target and any explicit short-target k-mer padding accounted for in the saved coordinates. The 43 matching outputs range from 771 to 788 bp. Matching extensions comprise 62.3% of saved unique sequences under this reference criterion; this is a reference-match fraction, not measured biological precision.

The run completed and its consistency validation passed. The search-limit flags remain material: these are results of the bounded search, not exhaustive enumeration or proof that each selected route is globally unique. Prediction files were fingerprinted and remained unchanged during evaluation.

## Fixed-input gains and losses

An input is counted as matching when at least one of its selected saved extensions passes the same reference threshold. Retained means the input still has a matching extension; it does not mean the output DNA is unchanged.

| Comparison baseline | Matching inputs before | Matching inputs now | Retained matching inputs | Newly matching inputs | No longer matching inputs | Net change |
|---|---:|---:|---:|---:|---:|---:|
| Previous workflow v5 | 18 | 48 | 16 | 32 | 2 | +30 |
| Historical continuity-filter result | 44 | 48 | 38 | 10 | 6 | +4 |

The previous workflow v5 baseline is archived as `v4_supplied81_crossing`; it is distinct from the failed development experiment named `revision5_supplied81`. That experimental revision is not the baseline used for the +32/−2 comparison.

At the unique-sequence level, the previous v5 baseline had 16 matches among 58 saved extensions, and the historical baseline had 44 among 65. The current result is 43 among 69. Therefore the current result improves the number of original inputs served by matching extensions compared with both baselines, but has one fewer unique matching sequence than the historical output. These facts should be reported together rather than substituting one denominator for another.

No complete saved extension sequence is identical, after strand normalization, to the previous v5 or historical output. This repair does not preserve all earlier predicted sequences. It preserves the 81 supplied targets, changes their chosen extensions, gains matching outcomes for some targets and loses them for others.

The two input IDs losing a match relative to workflow v5 are:

* `seed_0e3fc52adc6e1e7e`
* `seed_541b71503abcacf8`

The six losing a match relative to the historical result are:

* `seed_00b0f720c7a86847`
* `seed_03150ef5d6c40be9`
* `seed_0e3fc52adc6e1e7e`
* `seed_5f99a193e6244b29`
* `seed_bfe1e95d60e5270a`
* `seed_c8187729f2e6f27e`

The complete per-input table records every gain and loss, including unchanged original sequences and their associations. These IDs are audit identifiers, not biological names.

## Prediction and evaluation settings

All tests used the existing k=23 graph, a maximum selected flank of 750 bases, at least 200 scored external positions, a mature beam of 64, 256 pending context paths, at most two external-node visits, and 2,000,000 work units for each query pass. Ordinary discovery and controls used the 90% identity / 85% query-coverage mapping defaults; supplied targets required exact full-length mapping. Source-repeat crossing was enabled only for identified targets. Local conditional branch surprisal per scored position was the primary score; count variation broke ties. The reference was used only after prediction.

An extension match requires **at least 90% identity and at least 50% whole-query coverage in one BLAST HSP**. Coverage uses the complete saved extension as its denominator. The criterion does not require that the HSP include the original target interval, and does not establish physical linkage across every graph branch. Reference absence does not prove a false positive. Conversely, finding matching downstream sequence does not prove that every base of the reconstructed route occurred together in one molecule.

All 75 unique supplied targets also matched at the separate target criterion of 90% identity and 85% query coverage. Because this fixed target set was reference-selected, that target result is an input check, not a new discovery-performance estimate. Repeated development on this dataset also means the result is not independent validation.

## Input integrity and one-dataset reproduction

`input_integrity.json` verifies unchanged SHA-256 digests for all **28 original files**: the array report and all 27 graph files. The graph was not rebuilt or modified during these prediction runs. Evaluated prediction files also remained unchanged during reference evaluation.

After building the recorded v6 source on the server, run ordinary discovery on SRR14017553 from the repository root with 24 threads:

```bash
./build/mcaat phage-curate \
  --arrays /vol/d/data/copsac/SRR14017553 \
  --graph /vol/d/data/copsac/SRR14017553/graph/graph \
  --output /vol/d/data/copsac/SRR14017553/protospacer_v6 \
  --threads 24
```

Use a new or empty output directory. The omitted parameters take the recorded defaults above: 750 maximum flank, 200 minimum total flank, 200 minimum external scored positions, beam 64, pending pool 256, read length 150 and source crossing enabled. This command uses the existing graph and processes this one dataset. BLAST evaluation is separate; the benchmark runtime above was measured with two threads, not this 24-thread server command. No outcome for SRR14017554 or SRR14017569 is claimed here.

## Search-bound check

The same frozen implementation was also tested on all 81 supplied inputs with 1,024 pending context paths instead of 256. Both runs retain all 81 records, extend all 75 canonical targets and save 69 unique extensions with 43 reference matches, serving the same 48 original inputs. There are **zero matching-input gains and zero losses**. Five output strings change without changing any input's matching status. The two earlier context-pruning examples, Q000009 and Q000030, retain their same nonmatching selected sequences.

Runtime increases from 189.42 to 328.61 seconds, or 73.5%. The default remains **256 pending paths and 64 mature paths**. Increasing the search bound alone does not fix the remaining selection failures.

## Implementation validation

The native graph fixture checks one permitted source passage, source-node revisit/re-entry guards, exclusion of repeat counts, exact short-target coordinates, one selected flank, source-only discovery rejection, and rejection of an adapter core with one substitution. It also adds an already present downstream segment as another reported spacer without changing the graph: the selected route and all 13 internal k-mers of that segment remain unchanged and scored. This directly tests the distinction between pure spacer membership and repeat masking.

Arithmetic tests verify that branch support and count variation use the same adjacent unmasked positions, that trimming recomputes both summaries, and that path ranking and output consolidation remain consistent. Batch checks cover input selection, invalid/empty output handling and failure accounting.

All benchmark predictions use the frozen binary below. The final source differs from that binary only by four console-output lines reporting the actual target, extension and unextended-query counts after successful completion. The final compiled binary passed the same expanded native fixture; prediction and selection code are unchanged.

## Evidence files

The local audit bundle is `outputs/protospacer-method-development-2026-11-16/` in the development workspace. Paths below are relative to that bundle; graph and reference files are not part of the repository:

* `revision6_supplied81/predictions/run_summary.json` — completed run counts, method and exact settings.
* `revision6_supplied81/evaluation/evaluation.json` — validation, saved-file hashes and completed reference-match counts.
* `revision6_supplied81/evaluation/ContiguousPaths_coverage.tsv` — per-extension coverage outcomes.
* `statistics/revision6_supplied81_comparison.json` — exact cohort counts and input-level gains/losses.
* `statistics/revision6_supplied81_input_comparison.tsv` — complete per-input comparison.
* `revision6_full398/evaluation/evaluation.json` and `timing.json` — completed ordinary-discovery validation, matching counts and native runtime.
* `statistics/revision6_full398_comparison.json` and `revision6_full398_query_comparison.tsv` — identical-input comparison with the previous workflow and all query-level gains/losses.
* `statistics/full398_discovery_loss_summary.json` — audit of the rejected old detection contexts.
* `revision6_shuffled398/evaluation/evaluation.json` — completed control evaluation.
* `statistics/revision6_shuffled398_control_audit.json` and `revision6_shuffled398_match_trace.json` — unchanged control boundaries, generation checks and the two distinct matching controls.
* `input_integrity.json` — unchanged SHA-256 digests for the array report and 27 graph files.
* `revision6-implementation/manifest.json` — frozen prediction binary and source hashes.
* `statistics/context256_vs1024.json` — complete larger-context comparison.
* `final_fixture/validation.json` — regression checks on the final compiled binary.

The frozen MCAAT binary SHA-256 is `0b71afaaec20a53facee05cba0ab693a1591aed7b31efb632ff34a916c043e4c`. The final compiled source and binary fingerprints are recorded in `final-implementation/manifest.json`. The reference database fingerprints are recorded in `reference_manifest.json`.
