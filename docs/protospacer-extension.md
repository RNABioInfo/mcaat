# Native protospacer detection and extension (v8)

MCAAT searches an existing de Bruijn graph for a possible target of an array spacer, extends that target through competing graph paths, filters source-only and technical paths, and saves the selected sequence with its evidence. All of this runs inside the single `ProtospacerWorkflow` class. The graph remains unchanged. No reference sequence or external aligner is consulted during prediction.

## Run on an existing graph

From the repository root:

```bash
build/mcaat phage-curate \
  --arrays /vol/d/data/copsac/SRR14017553 \
  --graph /vol/d/data/copsac/SRR14017553/graph/graph \
  --output /vol/d/data/copsac/SRR14017553/native_protospacer_results \
  --threads 24
```

Use a new or empty output folder. This command reuses and retains the graph; it does not download reads or rebuild it. To extend an explicit known target set instead, add `--protospacers /path/to/targets.fasta` and use a separate output folder.

The batch wrapper `scripts/run_copsac_protospacers.sh` processes eligible SRR folders sequentially. It reuses an existing graph or downloads reads and builds one when needed. It defaults to 28 threads, a 192 GiB graph-build budget, and the same extension settings listed below: 200 scored external positions, a 200–750-base flank, 256 pending paths, 64 mature paths and a read length of 150. The batch invokes source crossing explicitly. Use `--existing-only` to prohibit downloading/building and `--dry-run` to list eligible samples without running them. Final result directories are recorded in each sample's `protospacer_workflow/LATEST.txt`; the batch's stage logs and measurements are under `_mcaat_workflow_logs`.

## Inputs and two distinct modes

In ordinary discovery mode, read the array reports and retain every source-spacer observation. Conservative boundary correction can move shared spacer-end motifs into the repeat: the default requires 85% support among at least five distinct spacers, a motif of at least six bases, and agreement across duplicate and reverse-complement observations. Group corrected sequences with their reverse complements into internal queries while retaining the original source associations.

Map each query to a connected graph path. The defaults require at least 90% identity and 85% query coverage. The mapper permits substitutions and terminal clipping, but not internal gaps. It searches seeds in the succinct index and verifies their connected extensions base by base. Mapping is bounded; a search that reaches a work or placement limit is not exhaustive.

In supplied-target mode (`--protospacers FASTA`), the inputs are already identified targets. Require an exact, full-length graph embedding: 100% identity, 100% query coverage, no clipping. Associate each input with source-array context by substring lookup on both strands; do not replace it with the longer source spacer. A target shorter than graph k is embedded as the exact suffix of a full k-mer. The preceding `anchor_prefix_bp` bases belong to the left flank, not to the supplied target. Saved target coordinates still identify the unchanged original sequence.

Every supplied FASTA record is preserved even if it cannot be extended. Unlike new discovery queries, explicitly supplied targets are allowed to map even when their sequence fails the discovery complexity test; an exact mapped target remains identified. Flank complexity and technical-adapter checks still apply. Supplied targets with no source-array association receive an explicit unresolved outcome. Supplied-input evaluation and ordinary discovery evaluate different populations and must be reported separately.

## Establish a target before allowing a source-repeat passage

First explore the neighboring context with source crossing disabled. The source repeat, in the orientation and position expected from the associated array, cannot establish an external target. Require a continuation incompatible with that source context and sufficient subsequent graph sequence. A shared k-mer alone is not a global rejection rule. This identifies a possible external target; it does not establish that the target is viral.

If crossing is enabled, give already detected placements an additional extension pass with a separate work allowance. Supplied targets are eligible because their identity is provided explicitly and their complete sequence must map exactly. The crossing pass cannot turn a previously rejected source-array-only placement into a newly discovered target. Retain the original independent context in `target_context.tsv` when evaluating a crossing extension.

The crossing pass may traverse the associated repeat region once. It cannot revisit an anchor or source-array node. A repeat/junction occurrence adds a k-mer to the repeat-context set; an occurrence wholly inside a spacer does not. The same k-mer can still enter this set if it also occurs at a repeat or junction elsewhere in the associated arrays. After leaving, entering a known repeat/junction node is rejected. This separation matters because a genuine target can share k-mers with another spacer. A source passage must clear by the expected unmatched spacer/repeat interval, its allowed edit/graph overlap, plus one context window; otherwise it stops with `source_passage_not_exited`. External nodes have a configurable visit cap, default two. Counts from the initial source-repeat passage are excluded from the downstream count profile; a repeat's pooled abundance should not set the extension's ordinary coverage baseline.

Unlike the previous policy, a short accepted strict extension does not prevent this comparison. Strict and permitted-crossing alternatives are compared before selecting the result. Consequently, this revision does not promise to preserve every previous strict result unchanged.

## Extension score: conditional graph-walk support

At a fork, compare the selected edge's multiplicity with the sum over valid alternatives in the traversal direction. If `c(e)` is the selected count and `C` is that sum, define the empirical local transition proportion `p(e | node) = c(e)/C`. Forward extension uses outgoing arcs; backward extension uses predecessor arcs.

```
B = sum over scored forks of -log2(p(e | node))
N = number of recorded, unmasked edge-count positions

path_cost = B / N
```

Lower cost is preferred. An unbranched position contributes zero decision cost. This is **conditional empirical graph-walk surprisal per scored position**: a locally weak choice costs more than a well-supported choice. It does not estimate the probability that an extension is genomically correct or belongs to a phage. Algebraically, `B` is the negative logarithm of the product of the scored local proportions, and `B/N` is the negative logarithm of that product's `N`th root.

The numerator accumulates decision costs; it is not averaged by the number of forks. At a fixed scored length, taking extra easy forks cannot improve the score merely by increasing a fork-count denominator. The normalization still has a length effect: sufficiently many easy added positions can reduce the average cost. Completeness, eligibility and bounded search are handled separately.

### Apply the same source mask to branch and count evidence

Store each raw local branch cost at its actual added-base position, including while source clearance remains unresolved. Once eligible downstream positions are known, include a branch cost only when the preceding and selected edge-count positions are both unmasked and consecutive. This excludes a source-exit transition whose preceding edge is masked, while retaining eligible downstream decisions encountered before clearance was confirmed.

Use the same adjacent-position condition for count variation:

```
E = sum over adjacent scored positions of
        [log2((count[i] + 1) / (count[i-1] + 1))]^2
mean_energy = E / number of those adjacent comparisons
```

`mean_energy` is a tie-breaker when primary path costs are equal. It is not added to `path_cost`, cannot compensate for a worse primary branch score, and has no fitted coefficient. After trimming a path, recompute both summaries from the remaining positions and the same mask. No score coefficient is fitted to reference labels.

The denominator at a fork includes valid local graph arcs before path-specific source, revisit and artifact guards are applied. It is not renormalized over only the branches surviving those guards. Masked transitions are omitted, so this is not the likelihood of the entire reconstruction.

These graph-local proportions are **not probabilities that an extension or phage prediction is correct**. Overlapping k-mer counts are correlated and can pool organisms, strains, repeats and sequencing errors. A count-1 route can be real; stable high counts can be misleading. A route with no scored forks has zero branch cost, which means that no branch comparison challenged it—not certainty or viral origin. When `N = 0`, zero is an unscored placeholder; output eligibility still requires the configured external evidence.

No reference sequence or reference-match label is used during prediction. The score ranks paths retained by a bounded search, so the selected route need not minimize it over the complete graph. Sequence-artifact checks, source/revisit rules, ambiguity reporting and exact target/graph validation remain separate checks.

## Keep bounded alternatives at the same extension distance

Expand one graph edge per round. Competing frontier paths therefore have the same physical distance from the anchor rather than advancing by unequal unitig lengths. Cache local adjacency queries to avoid fetching the same neighborhood repeatedly.

Keep separate bounded pools:

* **Pending:** a source passage or an external path with too little downstream count evidence; default maximum 256 paths.
* **Mature:** an external path with at least `max(k + 1, read_length - k)` recorded count positions; default beam width 64. For k = 23 and read length 150, this threshold is 127 positions.

The maturity threshold is a chosen observation scale related to one read's effective k-mer span. It does not mean those positions are 127 independent reads, and `--read-length` is supplied metadata rather than an estimate inferred from the graph. Keeping separate pools prevents an unscored source passage from being compared immediately with a path whose downstream support is already known.

Within a pool, prefer established external context and a full requested flank, then lower conditional branch cost, lower mean adjacent-count energy, longer retained sequence, and deterministic node order. Reserve different current oriented endpoints and context states before filling remaining slots by score. Rejoined variants therefore cannot occupy every slot ahead of a distinct continuation. Every retained state keeps its complete real path and revisit history; states are not spliced or asserted to be mathematically equivalent. These policies improve search diversity only within the configured bounds; no route-retention guarantee follows.

Store at most four mature beam widths of eligible terminal candidates. Trimming this terminal buffer does not terminate the live frontier. Beam-pruned states are diagnostic paths, not completed predictions. Reaching a work, placement, frontier or terminal-storage limit is recorded.

When any terminal candidate reaches the requested maximum flank, compare those full-length candidates first. Otherwise compare the retained shorter candidates. Retain the best actual graph route and annotate close competitors, their runner-up score and first divergence. Do not shorten a completed path solely because another route is plausible. Short bubbles that rejoin and then share at least k+1 oriented nodes are treated as local sequence alternatives on the same route; this is a heuristic route equivalence, not read-backed phasing.

The default `branch_ratio` is 1.5. A retained competitor is comparable when its conditional branch cost is at most 1.5 times the best cost (plus numerical tolerance). This is a heuristic cost tolerance, not a probability ratio or confidence interval. Source-region branching is also marked unphased even if its other arms were too short to enter the completed-candidate pool. A graph-valid route can connect a target to downstream sequence that never occurred with it in one molecule; the toy fixture deliberately records this limitation. `ambiguous=0` still does not establish genomic truth or exhaustive uniqueness. Competitor annotations concern the retained routes for a particular placement and side, not every possible target embedding or every unexplored path.

## Filter, retain one flank, and save

Stop before the exact shared TruSeq prefix `AGATCGGAAGAGC`, or either configured longer TruSeq core within one substitution. The shorter prefix requires all 13 bases with no mismatches; it catches technical paths whose continuation differs from the two longer cores. Check both strands, both traversal directions, and the target/flank boundary. The prefix is shared by the [Illumina adapter sequences](https://support.illumina.com/ko-kr/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html).

This check participates in discovery as well as extension: an adapter-containing path cannot establish independent target context unless enough clean context precedes the adapter. Stop before the motif, retain that preceding prefix, and continue competing branches. In SRR14017553, the old full-core check missed technical branches supporting 21 accepted targets; all 21 targets and their selected extensions failed the fixed reference criteria. These observations motivate the repair, but reference alignments are never read by prediction.

Also check 64-base sequence-complexity windows every eight added bases. A window is uninformative if triplet entropy is below 2 bits or one homopolymer occupies at least half its length. These finite motif/complexity rules do not detect every technical artifact.

A retained flank must have at least 200 added bases by default, at least 200 scored edge positions outside the source region, positive recorded support and informative sequence. `--min-extension-flank` and `--min-external-bases` control these distinct requirements. A 750-base path whose first 623 bases lie in unresolved source context is not a supported 750-base external reconstruction. The default maximum is 750 added bases. For short targets, the explicit anchor padding counts toward the left-flank minimum and maximum. Source revisits, source re-entry, failure to clear the expected repeat passage, excessive external revisits, saturated recorded multiplicities, technical motifs and search limits can shorten or reject a flank. Competing-route uncertainty is annotated, not used to manufacture a short confident-looking prefix.

When both sides qualify, keep one side using the same path ranking. Pooled edge counts do not establish that two independently chosen flanks belong to the same molecule. The result is an extended target region with one selected flank, not a complete phage genome. Short-target k-mer padding may still precede a right-side extension; it is the exact anchor embedding, not a separately selected second flank.

Select at most one primary placement per canonical query. Prefer detection and an accepted extension, then alignment coverage and identity, full requested flank, conditional branch cost, mean adjacent-count energy, length and deterministic anchor order. Equal-quality distinct target sequences remain unresolved. Verify the final sequence against the graph before writing it.

Consolidate exact/reverse-complement duplicate sequences and wholly contained extensions. Do not stitch partial overlaps into a new sequence. Strict representatives selected in this run take precedence in containment consolidation; this is not a promise that the preceding revision's primary selection is unchanged. Preserve every source/input association and its exact coordinates.

## Limits and reproducibility

| Parameter | Default | Meaning |
|---|---:|---|
| `--max-flank-length` | 750 | Maximum added bases on a side, including short-target left padding. |
| `--min-extension-flank` | 200 | Minimum total retained flank, including short-target padding. |
| `--min-external-bases` | 200 | Minimum scored positions outside the initial source context. |
| `--beam-width` | 64 | Mature frontier width; terminal storage is at most four times this width. |
| `--max-context-paths` | 256 | Pending/context frontier width. |
| `--read-length` | 150 | Used to set the downstream-evidence maturity scale. |
| `--max-placement-states` | 400,000 | Bounded index/connected-mapping work per query. |
| `--max-placements` | 128 | Maximum retained target placements; short-seed hit lists also have a 4,096-node cap. |
| `--max-query-states` | 2,000,000 | Main query work allowance; an enabled crossing pass has a separate allowance of this size. |
| `--max-node-visits` | 2 | Maximum visits to an external node; source/anchor rules are stricter. |

These are search bounds, not wall-clock timeouts. Completion means the bounded procedure finished, not that all graph routes were enumerated. Parameter values, work counts and limit flags are saved. SIGINT/SIGTERM request partial output and a nonzero exit status. Use a new or empty output folder and retain the graph and run summary when comparing versions.

There is no reference-guided branch choice, reference-derived path score or downstream BLAST call in MCAAT. Evaluate frozen predictions afterward with fixed identity/coverage criteria and unchanged inputs. A missing reference hit is not proof of a false positive; a target hit alone does not validate its extension. Report both sequence counts and input-query counts because consolidation can change their relationship.

Repeatedly evaluating method changes on one dataset makes it a development dataset. Independent data are needed before claiming general performance; reference-free prediction does not make repeated development-set evaluation independent.

## Output counts

| Field or file | What it counts or stores |
|---|---|
| `supplied_targets` / `SuppliedProtospacers.fasta` | Every supplied input record, unchanged, including unextended records. Only meaningful in supplied mode. |
| `spacer_observations` / `sources.tsv` | Ordinary mode: original source-spacer observations. Supplied mode: source-context associations, which may exceed the number of supplied input records. |
| `distinct_spacers` / `queries.tsv` | Internal canonical queries after grouping duplicates and reverse complements. |
| `detected_queries` | Queries with a selected possible target, whether extended or not. In supplied mode this is an extension outcome, not a new biological discovery count. |
| `protospacer_regions` / `Protospacers.fasta` | Unique selected target sequences after target deduplication. |
| `extended_queries` | Queries associated with an accepted extension. Several queries may share one saved region. |
| `unique_extensions` / `ContiguousPaths.fasta` | Unique retained extended regions after consolidation. This is not a count of phage genomes. |
| `limited_queries` | Queries that encountered any recorded search limit. |
| `queries_without_extension` | Canonical queries without an accepted extension, including detected-only and unresolved outcomes. |
| `supplied_targets.tsv`, `associations.tsv` | Original record → query → target → extension mappings, strands and exact intervals. |
| `candidates.tsv`, `unresolved_and_rejected.tsv` | Considered placements and reasons; alternative rows are not extra primary predictions. |
| `target_context.tsv`, `graph_paths.tsv` | Independent discovery context and selected connected graph routes. |
| `coverage_evidence.tsv` | Query ID, target/extension IDs, actual scored edge multiplicities and running minimum count. Query IDs distinguish different input queries that map to the same deduplicated target. The historical `path_capacity` column name no longer means segment-median capacity. |
| `path_selection.tsv` | Retained flank scores, phase, selection status, competing-route count, first divergence and unphased source branching. The legacy `independent_bases` field counts scored edge positions, not statistically independent observations or necessarily distinct nodes. `path_cost = branch_loss / independent_bases` for a nonzero denominator; `runner_up_path_cost` uses that same primary score. `mean_energy` is only a tie-breaker. `branch_observations` is an audit count, not the selection denominator. |
| `withheld_extensions.tsv` | One row per withheld selected route: target, reason, full sequence, block statistics, outward graph nodes and scored multiplicities. These are diagnostic routes, not accepted extensions. `count_positions` are one-based from the graph anchor; add `anchor_prefix_bp` for a short target’s left flank. |
| `query_metrics.tsv`, `run_summary.json` | Per-query work/timings, all parameter values, count meanings and completion status. |
| `search_paths/*.tsv.gz` | Optional finished/pruned path audit. These rows are search states, not protospacer or phage counts. `branch_costs` contains raw costs, including masked positions; its unfiltered sum is not `branch_loss`. Only costs whose preceding and selected positions both occur consecutively in `count_positions` contribute. |

Interval coordinates are zero-based and half-open. Supplied-mode query intervals cover the original complete supplied sequence. Left/right paths are oriented to the canonical query; association strand fields relate them back to the original source. Empty extension IDs explicitly mean no accepted extension. Summed worker times can exceed wall time because queries run concurrently.

`first_ambiguous_base` is a zero-based outward position from the graph anchor, interpreted only when `ambiguous=1`. For a short target's left flank, add `anchor_prefix_bp` to express it relative to the target boundary. Coverage-evidence positions are one-based outward added-base positions and already include that padding. The optional search audit's `count_positions` are one-based from the graph anchor and exclude padding; a cost at position `p` is `branch_costs[p-1]`. A `runner_up_path_cost` of -1 means no scored competing completed route was recorded; source-branch uncertainty alone can still set `ambiguous=1`.

## Coverage validation after selection

`--max-coverage-fold 2` validates the one selected extension. Set it to `0` to disable this filter. Search, candidate ranking, spacer detection and source-array handling run first. There is no fallback to a weaker placement, opposite flank or smoother alternative after the selected route fails. This avoids accepting an arbitrary route simply because one of many explored paths passes a smoothness filter.

For configured read length `R` and graph k, use blocks of `max(k+1, R-k)` consecutive scored edges: 127 positions for 150-base reads and k=23. Calculate the median multiplicity of every complete block of the selected flank. With at least two blocks, withhold the extension if the largest median exceeds the smallest by the configured factor. An incomplete trailing block is ignored. Paths with fewer than two complete blocks receive no decision from this check and remain subject to the other extension rules.

Counts from the source-array passage are excluded. Neither repeat multiplicity nor anchor multiplicity sets this baseline. Isolated count spikes or dips generally do not determine a block median. The read-length scale reduces sensitivity to adjacent correlated k-mer counts; it does not make the medians statistically independent.

The original target, source associations and strict discovery context stay intact. A withheld route is absent from `ContiguousPaths.fasta`; `withheld_extensions.tsv` retains its complete sequence, outward graph nodes and scored counts for inspection. Its target has reason `coverage_discontinuity` and an empty extension ID. Accepted routes are not shortened by this filter.

This is a tunable continuity heuristic, not a p-value or a guarantee of genomic linkage. Genuine coverage variation can also fail it. Report matching extensions removed as well as nonmatching ones removed; evaluate the setting on additional samples before treating the development result as general performance.

## Function responsibilities

All prediction functions belong to one class. The former unitig-median scoring helpers have been removed.

| Function | Responsibility |
|---|---|
| `ProtospacerWorkflow` | Own options, graph, observations and per-query results. |
| `main` | Dispatch help/self-test, validate startup, install signals and save failure status. |
| `parse_options` | Parse CLI options and validate files, budgets and output safety. |
| `print_help` | Describe the current options and count/search semantics. |
| `interrupt` | Record an interruption for workers to stop and save partial results. |
| `run` | Execute stages, schedule queries, collect errors and write final status. |
| `reverse_complement`, `canonical`, `dna` | Orient, group and validate DNA sequences. |
| `json`, `numbered` | Escape output strings and generate deterministic readable identifiers. |
| `triplet_entropy`, `informative` | Measure sequence complexity and apply entropy/homopolymer checks. |
| `adapter_position` | Find the exact shared 13-base TruSeq prefix or a longer configured core within one substitution, in the traversal orientation. |
| `repeat_context` | Classify the expected adjacent source repeat as pending, matching or incompatible. |
| `label`, `find_node` | Read a real k-mer label and retrieve a verified exact node. |
| `neighbors`, `added_base`, `spell` | Return real graph arcs/counts, emit their bases and verify a connected path's sequence. |
| `search_range`, `find_starts` | Search bounded colexicographic seed intervals and order candidate starts. |
| `verify_index` | Check exact/mutated seed lookup against sampled real graph nodes. |
| `load_sources` | Read array observations and preserve their identities and repeat variants. |
| `correct_boundaries` | Infer conservative shared-end motifs and save original boundaries. |
| `prepare_source_nodes` | Index both broad source-array nodes and separate repeat/junction nodes on both strands; pure spacer matches do not imply repeat context. |
| `load_supplied_targets` | Retain exact supplied records and associate their source contexts. |
| `group_queries` | Group canonical query sequences while retaining every source/input association. |
| `map_spacer` | Construct bounded connected mappings; require exact full embeddings in supplied mode. |
| `source_context`, `source_clearance` | Orient source context, test it and calculate the excluded initial repeat interval. |
| `trim_walk` | Truncate bases, nodes, positions and score evidence together. |
| `update_profile` | Apply the same consecutive unmasked-position rule to cumulative branch costs and squared log-count changes. |
| `walk_score` | Return cumulative branch loss per scored position; zero is an unscored placeholder when no counts exist. |
| `rank_walks` | Rank completeness, conditional branch cost and mean energy; reserve distinct current endpoints within bounded pools. |
| `save_walks` | Stream finished/pruned diagnostic states with phase-specific unique identifiers. |
| `extension_barrier` | Stop before configured adapter or low-complexity sequence while preserving valid coordinates. |
| `filter_selected_extension` | Measure complete read-span block medians after final selection; withhold an inconsistent extension without fallback, preserving its target, discovery context and route audit. |
| `advance_walk` | Store one real edge and its raw branch cost; apply source/revisit/artifact guards, count masking and the flank-length bound. |
| `retain_completed` | Remove ineligible terminal paths and cap terminal storage without stopping live paths. |
| `same_route` | Recognize bounded rejoining bubbles with shared downstream context. |
| `choose_walk` | Retain the best actual terminal route and annotate close competitors and unphased source branching. |
| `extend_flank` | Expand equal-distance states with separate pending/mature frontiers and record branch scores. |
| `assemble_candidate` | Preserve discovery evidence, verify spelling, filter both sides and retain one flank. |
| `select_candidate` | Choose one placement using alignment and extension evidence or mark target ambiguity. |
| `process_query` | Map, establish strict target context, compare eligible crossing alternatives and record work/timings. |
| `consolidate` | Deduplicate/contain selected sequences without inventing overlap joins or losing associations. |
| `write_outputs`, `write_summary` | Publish validated sequence/evidence tables and atomically save status/counts. |
| `self_test` | Exercise arithmetic, context, ranking, ambiguity and output invariants without loading a dataset. |

## Mathematical background and its limits

If outgoing-base counts are modeled as multinomial observations, `c(e)/sum(c)` is the maximum-likelihood proportion. Applying that model to pooled graph counts and using the proportions to rank a complete route are additional modelling choices: graph counts do not retain which adjacent decisions came from the same read. The score has no fitted weights and does not estimate a sequencing-error distribution.

Squared adjacent differences are a standard smoothness measure; the broader mathematical treatment is described by [Shuman et al. (2013)](https://documents.epfl.ch/users/s/sh/shuman/www/Papers/Journal/Shuman_et_al_SPM_2013.pdf). Here that measure is a descriptive tie-breaker, not an additive likelihood or penalty in the primary score. Overlapping-read coverage is correlated, as discussed by [Evans et al. (2010)](https://arxiv.org/abs/1004.5587); the scored positions must not be counted as independent reads. [Haploflow](https://pmc.ncbi.nlm.nih.gov/articles/PMC8287296/) is related work using local graph coverage for reconstruction, not validation of this specific MCAAT selector.

Local coverage consistency is also used in [MetaCortex (2023)](https://academic.oup.com/bioinformatics/article/39/1/btad020/6986127). Its adjacent-node model is related motivation, not the block-median rule implemented here. [García-López et al. (2015)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4585024/) document coverage variation and reconstruction limitations in viral metagenomes. Neither paper validates MCAAT’s twofold default; its measured trade-off must be reported separately.
