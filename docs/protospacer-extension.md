# Native protospacer detection and seed-local extension

MCAAT searches an existing de Bruijn graph for a possible target of an array spacer, extends that target through competing graph paths, filters source-only and technical paths, and saves the selected sequence with its evidence. All of this runs inside the single `ProtospacerWorkflow` class. The graph remains unchanged. No reference sequence or external aligner is consulted during prediction.

In ordinary discovery mode (v9), final protospacers require an accepted extension. `Protospacers.fasta` and `associations.tsv` contain only targets supported by a saved region in `ContiguousPaths.fasta`. Detection without an accepted extension is retained in `CandidateProtospacers.fasta` and `candidate_associations.tsv` for diagnosis. The candidate FASTA includes the final targets; do not add these two counts together. The new seeded comparison mode below preserves every input candidate independently of extension success.

## Extend from saved starts

This mode starts at the actual saved protospacer graph nodes. It constructs **only local paths, up to 750 added bases on each side**, without whole-graph assembly or rebuilding the graph. A full result contains the target plus up to 1,500 flank bases. The `evidence` mode generates competing routes, scores branch support, coding coherence and shared external sequence across arrays, and preserves all input targets:

```bash
build/mcaat phage-curate \
  --arrays /path/to/array_reports \
  --graph /path/to/graph/graph \
  --extension-seeds /path/to/previous_workflow_output \
  --output /path/to/new_seeded_results \
  --coverage-method evidence --max-flank-length 750 \
  --min-contig-length 500 --beam-width 64 --max-context-paths 256 --threads 24
```

`--extension-seeds` requires the existing `queries.tsv`, `sources.tsv`, `candidates.tsv` and `graph_paths.tsv` together. All detected/extended saved starts are used, including those without an earlier accepted extension. Query sequences and corrected source boundaries must match the array input; every anchor must spell its target on the loaded graph. This format currently supports one saved placement per query. This is an extension comparison, not a new test of discovery sensitivity across every possible graph placement.

For the earlier coverage-model comparison, choose `--coverage-method poisson`, `negbin` or `huber`. `all` runs those three sequentially with identical limits; it does not include the evidence selector. This option acts in seeded mode; ordinary discovery still uses the documented v9 workflow below.

### Select supported extension regions

Use `--coverage-method evidence` with the same saved starts. The native workflow preserves every protospacer, generates local candidates from the existing graph, and saves nonredundant extension regions. No reference sequence, external assembler or gene caller is used to choose a path.

The bounded generator keeps the Huber result and its retained flanks, plus lookahead terminal routes, prefixes and the selected compatible pair. Huber and lookahead both use the requested beam width (64 by default). They have separate query-work budgets. The graph is read-only. Source nodes cannot be revisited; other nodes retain the existing bounded-repeat rule. Actual graph spelling, source masks, technical-sequence barriers, sufficient external sequence and the 750-base per-side limit remain enforced.

Keep complete candidates and add seed-containing windows of 500 and 625 bp at the default limits: one facing left, one facing right and one centered on the target, where those windows fit. More generally, use `min_contig_length` and `(min_contig_length + max_flank_length)/2`. The target and its complete graph anchor must remain inside every window. Shortening a candidate never joins disconnected graph segments.

Rank candidates by external sequence supported by **distinct canonical protospacers** and **distinct source arrays**. Index every canonical `(2*k+5)`-base word outside the source mask (51 bases at k=23). Each target and each array votes once per word, regardless of how many competing path variants contain it. For a candidate's distinct external words W, calculate:

```
target_support = mean(log2(number of distinct targets supporting w))  for w in W
array_support  = mean(log2(number of distinct source arrays supporting w))
```

Average the percentile midranks for lower mean branch surprisal, higher stop-free coding coherence, and higher source-array support. Distinct protospacer support stays in the audit. First retain five distinct candidate families per input; then rank those families globally, keeping at most one representative per distinct target and canonical graph anchor. The direct 90%-shared-graph-run / 95%-overlap-identity grouping rule applies across starts. Equal scores prefer lower mean raw branch penalty, then longer sequence and deterministic sequence/query ordering. Feature values and combined ranks use nine-decimal rounding for stable ties. Other multiplicity features remain diagnostic audit columns.

The rationale is that competing branches reached from different spacer targets can corroborate an external sequence, while counting each target/array once prevents a large candidate pool from manufacturing votes. Multiple spacer hits have also been associated with more reliable virus–host links in [Roux et al. (2026)](https://pmc.ncbi.nlm.nih.gov/articles/PMC13739734/). Our shared-path score is an inference inspired by that idea: it uses graph candidates, not their reference alignments, and is neither that paper's method nor a calibrated probability. Related arrays, host sequence and shared repeats can still create correlated support.

`evidence_rank`, `evidence_score`, `target_support` and `array_support` identify saved choices and their support. Support values are −1 when no globally selected choice was assigned. Eligible but unsaved candidates have `candidate_not_saved` and rank −1; they do not count as FASTA outputs. All original targets remain in `protospacers.fasta`. The candidate audit retains per-position multiplicities, source-mask boundaries and branch/coding features. Retained alternative families are marked as competing routes rather than silently described as unique.

An optional five-column `--candidate-pool` replays frozen graph-generated candidates. Every route is checked against the actual graph. Evaluation columns are rejected. Evidence mode accepts at most 512 cached rows per input; mixture mode retains its 256-row bound. A cache omits original search-limit flags, which is recorded by `original_search_limits_known: false`; replay is not proof of exhaustive graph search.

`--expand-candidates true` is a separate experimental comparison. It preserves the base pool and adds bidirectional-support and local-variation searches, with an endpoint/abundance-stratified beam and postponed twofold coverage gate. The broader pool did not improve the measured result; it is not part of the recommended bounded method.

| Function | Responsibility |
|---|---|
| `mixture_candidates` | Generate the protected Huber and lookahead base pool; Huber uses width B and lookahead uses min(256, 2*B). |
| `add_evidence_windows` | Trim actual graph walks to the fixed seed-containing windows while keeping full candidates. |
| `evidence_candidates` | Add the two optional broader searches, preserving bounded repeats and source-node barriers. |
| `select_evidence` | Replay graph routes, calculate source-masked features and write the compact audit. |
| `coding_evidence` | Compare the longest stop-free run on either strand with a composition-adjusted background; this is not gene calling. |
| `select_evidence_regions` | Count distinct target/array support, rank local families, then save nonredundant globally ranked regions with one contig per target. |
| `same_graph_region` | Apply the unchanged grouping criterion using an exact central-core search; exhaustive-offset tests check equivalence. |

Run `tests/check_seeded_extension.py build/mcaat --lookahead` for route preservation, source masks, replay and window bounds. `--evidence-cap` tests 51 eligible inputs, 50 outputs, a complete omitted-candidate row and preservation of all 51 targets. The original three-model fixtures remain applicable.

Evaluate a completed run with `scripts/evaluate_protospacer_evidence.py`. Its `load_hit_cache` checks frozen prediction hashes before reusing evaluations; `figures` plots top-K matches, the logarithmic shared-word support function and observed multiplicities; `main` validates inputs/FASTA counts, evaluates the frozen outputs, computes gains/losses and writes the report. The predictor never reads the reference or these evaluation columns.

**SRR14017553 archived result (2026-11-17):** 164 starts, 159 preserved target sequences, 50 saved contigs and 19 reference matches. The 19 distinct sequences had best qualifying alignments in nine intervals across seven reference records. A match used at least 90% identity and 50% full-query coverage in one HSP. This configuration permits shorter windows and single-sided paths.

For a different sample, run ordinary discovery against that sample's existing graph. Pass its output directory to the seeded command above. Each dataset needs its own array reports, graph and saved discovery output:

```bash
sample=/vol/d/data/copsac/SRR14017553
build/mcaat phage-curate --arrays "$sample" --graph "$sample/graph/graph" \
  --output "$sample/evidence19_seeds" --threads 24
build/mcaat phage-curate --arrays "$sample" --graph "$sample/graph/graph" \
  --extension-seeds "$sample/evidence19_seeds" --output "$sample/evidence19_extensions" \
  --coverage-method evidence --max-flank-length 750 --min-contig-length 500 \
  --beam-width 64 --max-context-paths 256 --threads 24
```

Both output directories must be new or empty. The seeded FASTA is ordered by rank. Use its first 50 entries for the fixed-count benchmark. Add `--max-reported-contigs 50` to save just that subset. The default `0` saves every eligible representative, bounded by the number of distinct protospacers and graph starts.

| Flag | Meaning |
|---|---|
| `--min-contig-length 500` | Minimum total length including the target. |
| `--max-contig-length 0` | Maximum total length. Zero disables this filter. |
| `--max-flank-length 750` | Maximum bases to search on each side. |
| `--min-flank-per-side 0` | Minimum added bases required on each side. Zero allows a single flank. |
| `--min-branch-support 0` | Optional minimum `2^(-mean branch penalty in bits)` across scored external positions. A path with no scored forks can reach 1, so treat this as a graph statistic. |
| `--max-reported-contigs 0` | Optional FASTA output limit applied after ranking and grouping. |

Length and support filters run before ranking and can change which candidate wins. The optional FASTA limit takes a prefix of the same rank order. `extensions.tsv` reports `branch_support` and `same_start_representative_query`; the latter links input queries with the same graph anchor to the representative query. Use the representative's `reported` field to determine whether it was included in the FASTA.

A direct rerun of the restored search reproduced the 11,295 archived SRR14017553 routes and all 50 ranked contigs. Replaying those routes with the new selector saved 117 representatives with 28 matches. Its first 50 contained the archived 50 sequences and 19 matches, with a small order difference. A pre-ranking cutoff of `--min-branch-support 0.998` saved 87 representatives with 24 matches and retained seven of the original 19 matching sequences. The cutoff is opt-in because path selection changed on the test sample.

### Joint local coverage fitting (experimental)

Use `--coverage-method mixture` with `--extension-seeds`. Keep the ordinary/default workflow unchanged. Native generation protects the Huber candidate and its eligible independent flanks, then adds lookahead terminal routes, prefixes and its selected pair. Huber uses beam width `B` (default 64); lookahead uses `min(256,2*B)`. Each generator has its own `--max-query-states` budget. This is **two bounded searches per start**, not an exhaustive graph search. The 750-bp per-side maximum and 500-bp total minimum still apply.

An optional `--candidate-pool FILE.tsv.gz` replays already generated routes. Its header must be exactly `query_id`, `sequence`, `seed_begin`, `seed_end`, `origin`, separated by tabs. Coordinates are zero-based, half-open. Other columns, including reference annotations, are rejected. Every cached route is replayed through graph edges and traversal barriers, checked for bounded revisits and verified by spelling its connected nodes. The origin string is diagnostic, never a fitting feature. Cache input is limited to 256 records per query. Empty candidate pools leave the original target intact and unextended.

For candidate routes and nearby background routes, let `A[b,j]` be the number of visits by route `j` to an edge in block `b`, and `a[j] >= 0` its fitted abundance. The predicted multiplicity is:

```
mu[b] = sum_j A[b,j] * a[j]
objective = sum_b w[b] * ((observed[b] - mu[b]) / sigma[b])^2
            + log(1 + sum_b w[b]) * number_of_active_routes
```

Edges shared by two routes receive their **sum**, so a 10-copy and a 30-copy route can explain a shared 40-copy segment. This avoids forcing every segment into one flat target baseline. It cannot identify which equal-abundance entrance connects to which exit when both decompositions fit equally well.

For fitting, combine all compatible maximal left and right flanks; shorter output prefixes do not become separate biological molecules. If this exceeds 1,024 model routes, mark the start unresolved instead of pruning competitors and claiming confidence. This matters: comparing one complete two-sided route against only one-sided alternatives can falsely favor that complete route. The local graph includes up to 64 background routes using actual off-candidate incoming/outgoing arms and an existing candidate boundary, with slots distributed round-robin across distinct flanks. Hitting this background limit is recorded in the candidate audit and the selected start’s limit flag. Removing the target alone does not create an artificial background boundary. Outside arms follow at most one read-span of nonbranching graph, with no revisit. Connected source `repeat + spacer + repeat` paths on both strands are additional nuisance components. They are never output as targets solely because they explain shared counts. This is a bounded approximation: some real background routes may be missing.

Consecutive edges with identical route-incidence vectors form blocks, each at most `R=max(k+1,read_length-k)` edges. `observed` is the block median; `w=block_edges/R` gives overlapping k-mers read-span weight rather than counting them as independent reads. `sigma` is the larger of the scaled within-block MAD and `sqrt(max(1, observed + observed^2/dispersion))`, with `--nb-dispersion 10` by default. This is a weighted least-squares approximation using an overdispersed count variance, **not a negative-binomial likelihood**. The route-count penalty is complexity-inspired, not a proven BIC calibration for these correlated data.

A native greedy nonnegative solver adds at most 12 routes, refitting abundances with at most 80 coordinate sweeps per addition. It is not an exact minimum-flow decomposition. Repeat fitting 16 times after deterministic Gaussian block perturbations at scale `sigma`, clipped at zero. Identical design columns share support rather than receiving artificial certainty from tie order. A reported candidate needs support in at least 80% of these trials, positive support in the unperturbed fit, and standardized global fit error at most 2. Support for an exact anchored prefix includes longer routes containing that prefix. These are explicit heuristics, **not probabilities that a contig is correct**.

When both a prefix and its longer continuation pass, keep the longer continuation. Among other candidates, rank stability, fit error, fitted support and the existing tie-breakers. Keep one sequence per start, apply direct representative grouping, and save at most 50 representatives. Preserve all supplied target sequences irrespective of extension outcome. The fitter uses no reference, no read aligner and no external optimizer. Background abundance does not turn a graph route into a proven viral sequence.

This follows the general idea of [inexact flow decomposition under noisy weights](https://pmc.ncbi.nlm.nih.gov/articles/PMC9700332/). It is not an implementation of that paper's exact optimization, nor of [Phables](https://pmc.ncbi.nlm.nih.gov/articles/PMC10563150/), which also uses evidence unavailable in this graph-only experiment.

The native output adds no per-query file tree. `candidate_paths.tsv.gz` stores every model-eligible target route, its origin, stability, fitted support, eligibility, model/background route counts, background-limit flag and observed/fitted blocks intersecting that candidate. `base_abundance` is the sum of abundances of complete model routes containing the candidate, not necessarily one independent molecule. `extensions.tsv` retains one row per start, including unresolved starts. `protospacers.fasta` preserves all input targets; `viral_contigs_mixture.fasta` contains the selected representatives. The benchmark adds reference results **after** these files are frozen.

| Function / struct | Responsibility |
|---|---|
| `SavedExtension` | Sequence, anchor coordinates and diagnostic origin of a cached route. |
| `CoverageBlock` | One contiguous incidence-consistent block, observed median, noise scale, weight and visit vector. |
| `CoverageFit` | Nonnegative route abundances, fitted block counts, objective and normalized error. |
| `load_candidate_pool` | Read only the five allowed cache fields; reject unknown queries, malformed rows and evaluation columns. |
| `replay_candidate` | Verify target coordinates, reconstruct the connected graph walk and multiplicities, enforce source/revisit/artifact barriers. |
| `mixture_candidates` | Generate protected Huber routes and lookahead alternatives with separate work budgets and local mode overrides. |
| `fit_coverage` | Approximate the sparse weighted nonnegative fit using greedy additions and coordinate refits; no graph or reference access. |
| `select_mixture` | Build the local mixture, add background routes, form coverage blocks, perturb/refit, select a stable route and write its audit. |
| `evaluation_cache` | After prediction, reuse exact canonical-sequence evaluations with source hashes and matching metrics/reference path; reject newer reference indexes. Earlier caches lack a reference-content hash, which is explicitly recorded as a limitation. |
| `make_mixture_plots` / `write_mixture_report` | Show fixed-ranking performance, observed/fitted counts, actual multiplicities and exact outcome counts/limitations. |

Native arithmetic tests include a 10+30 shared-segment example and an equal-depth ambiguity. `tests/check_seeded_extension.py --lookahead` also checks mixture live generation versus replay, 750-bp linear extension, an equal-abundance shared repeat that must not be resolved across its split, preserved targets, connected sequence spelling and unchanged graph files. The ordinary three-model checks remain separate.

**SRR14017553 result (2026-11-17):** all **164** starts were tested and all **159** targets preserved. Joint fitting accepted **59** start records (58 distinct targets), giving **36** saved representatives of 519–1,539 bp. **7/36** representatives matched the reference at 90% identity / 50% full-query coverage in one HSP. This is worse than lookahead's **15/50**, and below the requested 25-of-50 goal. The method remains opt-in and is **not recommended as a replacement**.

Per-input matches: **13**, versus lookahead's 26 (**7 gained, 20 lost, 6 retained**) and Huber's 60 (**3 gained, 50 lost, 10 retained**). The candidate audit retains all 148 previously accepted Huber sequences. Matching alternatives exist for 72 inputs; the filter leaves only 16 with an eligible matching alternative (7 inputs fail fit error; 49 fail stability/positive support). Selection loses 3 more, leaving 13. Two matching input sequences are associated with a nonmatching representative during grouping. These are distinct stages, not counts to add together.

The complete bounded model uses at most 700 target-route combinations on this dataset. Background exploration reaches its cap for 135 starts; this is recorded, not treated as exhaustive. The real-data comparison replays the frozen Huber/lookahead routes through the actual graph. Live generation and replay agree on the controlled fixtures, including a shared repeat with equal-abundance entrances/exits that must remain unresolved. This does **not** substitute for a full live candidate-generation rerun of the real dataset. Both existing source searches were already bounded.

The final run used two threads, 2.50 GiB peak sampled RSS and approximately 139 seconds including graph loading; local candidate replay/model fitting took 46.7 seconds. All 27 graph components stayed unchanged. A provisional asymmetric model scored 12/44 but failed the equal-abundance repeat fixture and was discarded; the final comparison includes compatible complete flank combinations and uses the same thresholds. Evaluation reused exact-sequence matches from verified prior output artifacts; the source evaluations used the same local reference and metrics. No reference was passed to MCAAT.

### Branch lookahead and independent prefixes (experimental)

Use `--coverage-method lookahead` with the same saved seeds. The ordinary workflow and `all` remain unchanged. This experiment separates the candidate pool from the final reference-free selection. The real-graph test uses `--beam-width 128 --max-context-paths 256`, with up to 750 added bases per side.

At an unresolved fork, inspect unitigs within `2*R` bases (`R=read_length-k`, bounded below by `k+1`). A breadth-first probe visits at most 256 distinct unitig-entry edges, with at most 512 graph operations per new probe. Each flank can spend at most one eighth of the query budget on these optional probes, preserving most of the budget for actual extension. Exhausted probes fall back to immediate edge counts; they do not stop the walk. Only fully unmasked segments of at least `k+1` edges contribute depth. After calibration, each arm's weight is `exp(-min Huber(log((depth+0.5)/(mu+0.5))))` across observed downstream segments; without such a segment, use the immediate edge. Before calibration, use the largest observed downstream depth, falling back to the immediate count. Normalize these weights across retained arms for the branch-surprisal term. These are ranking weights, **not posterior probabilities**. The first baseline remains the median of that path's first read-span of external counts.

Looking beyond short connections for coverage evidence is related to [metaSPAdes](https://genome.cshlp.org/content/27/5/824). This implementation is a graph-only heuristic; it does not implement metaSPAdes' read-pair support or its complete assembly procedure. Optimistically taking the best future depth can also favor a wrong arm that reaches a similarly covered region.

- Protect alternative families defined by the first two unresolved fork choices, allocating beam slots round-robin within the score-sorted pool. No hard weak-link pruning is applied. Short dead-end tips and clearly dominant reconverging bubbles retain their local rules.
- Save supported checkpoints every read span after two spans and before later forks, as well as actual terminal routes. Twofold block-median continuity and existing source/technical barriers still apply. Keep up to 16 terminal routes and 16 prefixes per side, with family quotas.
- Consider compatible left/right pairs and one-sided alternatives. Require at least 500 total bases and at least two read spans of external coverage. Minimize the worse flank score, using total length to break equal scores. A poorly supported second flank is not mandatory. Source nodes cannot be revisited; other nodes use the configured bounded revisit rule.
- Rank and group the selected per-input sequences without reference data. Save the first 50 representatives. The 50 limit is an output budget, **not a demonstrated accuracy threshold**. All candidates remain in the evidence table; `reported` and `representative_rank` distinguish saved output from the pre-cap pool.

`candidate_paths.tsv.gz` is the only extra native file: up to 32 retained one-sided routes per flank, filtered to at least 500 bp, with sequence, target coordinates, score, stop and limit flags. Each is checked against its actual connected graph walk. The selected paired sequence is in `extensions.tsv`. The pool is bounded, not exhaustive; failure to retain a route does not prove it is absent from the graph. Both final alternative retention and the main beam can truncate exploration.

Evaluate only after freezing predictions. The benchmark records top 10/20/30/40/50 matching representatives, per-input gains/losses against frozen Huber, and candidate-pool availability. A retrospective grouping of matching candidates is explicitly diagnostic; reference hits never choose predictions. `selection.png` shows the measured top-K curve and loss function; `multiplicities.png` shows actual counts. Use `python3 tests/check_seeded_extension.py build/mcaat --lookahead` for the low-count bridge and independent-flank checks.

**SRR14017553 result (2026-11-17):** all 164 saved starts processed; all 159 candidate targets preserved. The fixed top 10/20/30/40/50 contain **5/7/9/12/15 reference matches**, respectively, at 90% identity and 50% full-query coverage in one HSP. The goal of 25 within 50 was **not reached**. The 50 final contigs are 501–871 bp. Per-input matching predictions fell from Huber's 60 to 26: **7 gained, 41 lost, 19 retained**. This should not replace Huber based on these results.

The retained pool contains 3,324 one-sided alternatives plus selected sequences, totaling 3,329 unique sequences. Of these, 541 pass the alignment metric, spread across 42 inputs; they are alternative variants, not 541 independent recoveries. For 16 inputs a matching alternative was retained but not selected. Of Huber's 60 matching inputs, 28 have no qualifying match in this retained pool. Nineteen selected representatives match before the output cap; 15 remain in the top 50. The final matches cover seven overlapping-best-hit reference interval groups, not seven proven distinct phages.

Prediction took 9.60 minutes on two threads with 2.54 GiB peak sampled RSS; all 27 graph components were unchanged. An initial diagnostic exhausted its work budget in optional probes and was corrected **before** reference evaluation; it was not used for accuracy claims. The final search still has 29 flank stops at the work limit, and it is not exhaustive. Native self-tests, the three coverage-model fixtures, the flow fixtures and the new low-count bridge/independent-flank fixtures passed. A first flow-fixture attempt failed inside legacy tiny-graph construction before extension; the repeat passed.

### Flow-guided extension and stopping at ambiguity

Use the same command with `--coverage-method flow`. This experimental mode follows one supported prefix from each saved anchor. It does not use a beam to guess across an unresolved split. `all` retains its original meaning (the three coverage-loss models).

Initial support `F` is the minimum observed multiplicity along the complete saved anchor. A one-node anchor uses the smaller total of its incoming and outgoing multiplicities. Probe nonbranching arms only within the remaining 750-bp flank budget. An arm's mean multiplicity supplies its weight `c[i]`; when it contains masked source-repeat positions, use its first edge count instead. At each junction:

```
leak[i] = sum(c[j] for j != i)
next[i] = min(c[i], F - leak[i])
```

Continue only when **exactly one** `next[i]` is positive. Several positive continuations cause `unresolved_flow_branch`; none cause `insufficient_seed_flow`. Both stops occur before adding the next edge, leaving the supported prefix intact. Set `F` to the chosen residual. Merges and high-copy repeats never replenish it. Genuine short tips are removed locally. A reconverging sequence bubble is simplified only when its strongest arm exceeds the runner-up by `--branch-ratio` (default 1.5); tied variants stay unresolved. Other weak arms consume flow rather than being silently discarded.

The excess-flow idea follows [Sena et al., Flowtigs (2024)](https://doi.org/10.1016/j.isci.2024.111208). Unitig depth estimates, the cap on residual flow, local graph simplifications and requiring a unique continuation are our conservative adaptations. This is **not exhaustive flowtig enumeration or a correctness guarantee for noisy sequencing counts**. Real counts need not obey flow conservation, and several genuine continuations can share one anchor. Either condition can cause early termination. No reference is consulted by this code.

Keep the existing source-crossing, source-revisit, bounded-repeat and technical-sequence barriers. When sides conflict, prefer the side retaining a larger fraction of its seed support, then length. As in the comparison modes, output requires at least 500 total bases and sufficient external coverage. All 159 input targets in the development set remain in the target FASTA regardless of this minimum; shorter attempted extensions remain in the evidence table.

Representative selection now ranks the **worst flank's support** before search truncation, ambiguity, relative-count variation and length. For flow, cost is `1 - remaining_flow / seed_flow`; for the other modes it is branch surprisal per external position. Taking the worse flank prevents a strong side from diluting a weak side. A longer prediction no longer wins merely by being longer. This ranking is evidence-based, not a probability of correctness.

`extensions.tsv` adds initial/residual flow, representative support cost and both decision traces. A trace entry is `position:node:incoming:outgoing_weights:choice:remaining`. Position is the number of added bases before the decision; weights are in retained-arm order; choice −2 means multiple supported arms, −1 means none. The retained prefix excludes that unresolved edge. Trimming restores the prefix's flow evidence rather than retaining decisions from discarded sequence.

Use `python3 tests/check_seeded_extension.py build/mcaat --flow` for four controlled graphs: a linear route, an unresolved split, a supported split and high-copy count inflation. Native self-tests also verify the formula and that support outranks representative length. To compare a completed flow run against a frozen Huber run, add `--previous-seeded /path/to/three_model_results` to the evaluation command below. This reuses the same identity/coverage definition and reports exact gained/lost input IDs.

### Shared graph rules

Follow unitigs directly from the target; at branches, probe only within the remaining flank budget. Apply local versions of the vendored MEGAHIT rules:

- Remove a genuine dead-end arm shorter than `2*k`. A probe stopped by a length/work limit is **not** a dead end.
- For alternatives reconverging at the same oriented node, keep the deepest simple bubble (`k+2` edges), or a longer bubble with at least 98% sequence similarity and the `20*k/0.98` length bound.
- Discard a weak arm only after external baseline calibration: its mean depth must be at most 10% of the neighboring-arm total **and** below this target's baseline by more than `--max-coverage-fold`. This protects low-count source-exit bridges and rare paths compatible with their own local coverage. Depth-based decisions require complete probes and no masked source-repeat counts in the compared arms. Setting the fold limit to zero also disables weak-link pruning.

These are local adaptations of [`bubble_remover.cpp`](../libs/megahit/src/assembly/bubble_remover.cpp), [`bubble_remover.h`](../libs/megahit/src/assembly/bubble_remover.h) and [`weak_link_remover.cpp`](../libs/megahit/src/assembly/weak_link_remover.cpp), not MEGAHIT's full global, iterative, multi-k assembly. Read-pair support is not available from this graph alone. The graph is never edited.

Keep competing unresolved branches in a bounded beam: 64 mature and 256 pending states by default. Each state retains its own actual walk. A source region can be crossed once by an already identified target; source/anchor revisits and source-array re-entry are forbidden. Other nodes permit at most two visits. Mask the source passage when measuring external multiplicity. Stop at technical adapters, low complexity, saturated counts, dead ends, revisit bounds or 750 added bases.

### Coverage loss and acceptance

For each flank, estimate `mu` as the median of its first `R=max(k+1, read_length-k)` unmasked multiplicities and freeze it. At k=23 and read length=150, `R=127`. This initial window calibrates the baseline; count residuals are evaluated afterward. At each scored fork retain branch support `p=c/chosen_direction_total_count`. Before baseline calibration, rank by mean branch surprisal, breaking ties by adjacent log-count variation. After calibration, combine branch support with one of these per-position count losses:

| Model | Loss | Interpretation |
|---|---|---|
| Poisson | `2*[c*ln(c/mu) - c + mu]` | Deviance for a count model with variance `mu`. |
| Negative binomial | `2*[c*ln(c/mu) - (c+r)*ln((c+r)/(mu+r))]` | Allows variance `mu + mu²/r`; default `r=10`, adjustable with `--nb-dispersion`. |
| Huber log-ratio | `z=abs(ln((c+0.5)/(mu+0.5))); z²/2 if z<=1, otherwise z-1/2` | Relative-count loss with reduced influence from extreme deviations. |

Use the continuous limit `c*ln(c/mu)=0` at `c=0`. Count deviances follow the standard [GLM deviance definition](https://www.statsmodels.org/stable/generated/statsmodels.genmod.families.family.NegativeBinomial.deviance.html). Negative-binomial modeling has precedent for overdispersed k-mer spectra in [GenomeScope 2.0](https://pmc.ncbi.nlm.nih.gov/articles/PMC7080791/). The robust loss follows [Huber's location-estimation work](https://doi.org/10.1214/aoms/1177703732); applying it to log multiplicity ratios is our adaptation. These sources do not validate this particular extension workflow.

The final flank score is:

```
S = sum(-ln(p)) / N + a * mean(count_loss_after_calibration) / R
a = 1/2 for Poisson/NB deviance; 1 for Huber loss
N = all unmasked external positions, including calibration
```

Branch terms use the same consecutive external-position mask as the established workflow. The denominator includes all valid neighboring arcs before local pruning. Coverage alone is insufficient: a flat, weakly supported count-1 route can otherwise beat a supported route near count 40. Branch surprisal penalizes that weak choice. The initial window is not also penalized as if it were a stationary test interval; source exits can change abundance before the baseline is established.

`1/R` is an explicit conservative weighting for overlapping counts, **not** an exact effective-sample-size calculation. The median baseline, fixed NB dispersion, beam limits and log transformation are design choices. The combined score is a composite ranking rule, not a calibrated likelihood. **It does not supply probabilities of correct sequence or viral origin.** Each branch may establish a different initial baseline; ambiguous wrong branches remain possible.

Filter terminal candidates before choosing a route: complete read-span block medians must remain within a twofold range (`--max-coverage-fold`). Prefer a qualifying full 750-base flank when available, otherwise a shorter supported flank. When joining sides would revisit forbidden source nodes or join incompatible depth baselines, retain the longer supported side, breaking equal-length ties by loss. Require at least 500 total bases and at least `max(min_external_bases, 2*read_span)` scored external positions. Close alternatives and search truncation are reported rather than presented as certainty.

### Compact outputs

| File | Contents |
|---|---|
| `protospacers.fasta` | Every input candidate target, including those without a retained extension. |
| `viral_contigs_<method>.fasta` | One representative sequence per grouped region; only retained contigs of at least 500 bp by default. These are candidate viral regions, not confirmed phages. |
| `extensions.tsv` | One row per input start per method: outcome, actual sequence, source observations, target coordinates, both count traces, scored positions, limits, stops and representative association. |
| `summary.json` | Parameters, counts, timings and completion status. |

Group directly to a representative when at least 90% of the shorter graph-node path is one identical contiguous run and overlapping sequence identity is at least 95%. Both orientations are checked. Grouping is not transitive. A representative need not contain every associated target variant: the **actual per-start sequence and its target coordinates remain in `extensions.tsv`**. Coordinates are zero-based, end-exclusive in that row's sequence. `*_counts` includes each added edge; `*_scored_positions` is a one-based index into that side's count vector, ordered outward from the seed.

### Function guide and checks

All search logic remains in `ProtospacerWorkflow`; `LocalSegment` and `FlowDecision` only store evidence.

| Function | Responsibility |
|---|---|
| `load_extension_seeds` | Validate input identities, boundaries, node IDs, connected anchors and exact target spelling. |
| `probe_segment` | Inspect one bounded unitig arm, recording multiplicity, junctions and genuine dead ends. |
| `sequence_similarity` | Banded edit-distance comparison for longer local bubbles. |
| `count_loss` | Evaluate the selected mathematical coverage penalty. |
| `score_seeded_walk` | Freeze a flank's initial baseline and accumulate post-calibration loss; reset correctly after trimming. |
| `seeded_walk_score` | Combine branch surprisal with the read-span-weighted coverage loss. |
| `extend_seeded` | Apply local graph rules, explore bounded alternatives, filter and select actual flanks; optionally return retained prefixes and terminal routes. |
| `future_depths` | Probe bounded downstream unitigs for unmasked mean-depth evidence at a fork. |
| `reserve_branch_families` | Reserve beam or audit slots across the first two unresolved branch decisions. |
| `compatible_flanks` | Reject joined sides with incompatible depth baselines or forbidden node revisits. |
| `same_graph_region` | Test direct representative grouping using shared graph-node runs and sequence identity. |
| `run_seeded` | Run the methods on identical starts, validate joined walks, group results and save compact outputs. |
| `flow_choice` | Compute residual support for each arm and require a unique positive continuation. |
| `extend_flow` | Follow a single prefix with local simplifications, flow accounting and explicit ambiguity stops. |
| `support_cost` | Measure the least-supported retained flank without rewarding extra length. |
| `supported_representative` | Rank actual predictions by support, uncertainty, variation, then length. |

`tests/check_seeded_extension.py` builds a small temporary graph and checks both 750-base flanks, target preservation, bubble/tip/weak-link handling, source masking, bounded probes, and unchanged graph files. It also checks that a weak flat count-1 branch cannot displace a supported count-21 continuation after a low-count source exit, and that a rare correct arm survives an abundant neighbor both before and after calibration. `phage-curate --self-test` checks the loss formulas and region grouping. Temporary fixtures are removed automatically.

`scripts/seeded_extension_benchmark.py` evaluates an already completed three-method run against a **local** BLAST database, adds per-start reference evidence to the existing table, and creates `coverage_models.png` and `multiplicities.png`. It never changes prediction paths or uses the reference during search. It compares matched starts with an earlier run and records exact gained/lost query IDs. The criterion is one HSP with at least 90% identity and 50% of the entire extension query covered; targets use 85% coverage. Unmatched candidates are not thereby proven false, and overlapping reference loci are not independent phage discoveries.

The benchmark helpers have narrow roles: `fasta`/`rows` load records and tables, `canonical` normalizes strand, and `digest` fingerprints inputs. `same_region` mirrors native grouping so the baseline is counted consistently and each new association can be checked. `align` streams HSPs and retains the best qualifying match, or the best nonqualifying match when none passes. `reference_loci` groups overlapping best-hit intervals. `make_plots` draws the mathematical losses and observed traces; `write_report` writes the count tables and limitations. `main` verifies the unchanged input set, performs evaluation, and saves the augmented evidence and summary. `make_flow_plots` draws excess-flow curves and measured multiplicities; `write_flow_report` reports strict-flow outcomes and gains/losses against frozen Huber results. `check_flow` in the existing test script builds and removes the four controlled fixtures. `check_lookahead` verifies retention of a true 750-base path through a low-count connection, independent flank selection, target spelling and unchanged graph files. `make_lookahead_plots` draws the fixed ranking curve, loss and measured multiplicities; `write_lookahead_report` keeps final representative counts distinct from per-input and candidate-pool outcomes. `diagnose_lookahead_selection` attributes evaluated misses to per-input selection, representative grouping or the output cap, without altering predictions.

Run the controlled check with `python3 tests/check_seeded_extension.py build/mcaat`. For a completed three-method run, evaluate with:

```bash
python3 scripts/seeded_extension_benchmark.py /path/to/new_seeded_results \
  --baseline /path/to/previous_workflow_output \
  --baseline-evaluation /path/to/previous_blast_evaluation \
  --reference /path/to/local_blast_database_prefix --threads 2
```

### Measured development result: SRR14017553

The v3 comparison below predates the support-first representative change. All three models were tested on the same 164 saved query/start records: 160 distinct oriented anchor paths and 159 candidate target sequences, originating from 398 distinct corrected spacers. All 159 targets were preserved. Extension matches require one HSP at 90% identity and 50% coverage of the whole predicted sequence; the reference was used only after prediction.

| Method | Representative contigs | Matching representatives | Matching input records | Gained / lost records |
|---|---:|---:|---:|---:|
| Previous workflow, same grouping | 30 | 5 | 60 | — |
| Poisson | 50 | 3 | 36 | +0 / −24 |
| Negative binomial | 52 | 6 | 34 | +3 / −29 |
| Huber log-ratio | 51 | 8 | 60 | +6 / −6 |

Huber performed best among the three models but produced **no net increase in matching input records**. Multiple input records can describe the same region; representative contigs are partial candidate regions, not counts of distinct phages. The earlier workflow retained one flank, whereas these methods can join both; this changes the coverage denominator. This is a development comparison on one dataset, not independent validation or a reason to replace the ordinary workflow by default.

All three searches together took 10.24 minutes with 2 threads and 2.54 GiB peak RSS. Accepted contigs were 773–1546 bp. The graph remained unchanged. Controlled tests passed for eight cases under each model. Early weak-link pruning was corrected after tracing 30 lost matching routes: source-exit branches must survive until the external baseline is established. Search limits and ambiguity remain recorded in the output table.

The subsequent strict flow test processed all 164 records and preserved all 159 targets, but only four records reached the output requirements. They produced three representative contigs, one matching the reference, representing two matching input records. Against the frozen Huber run this is **0 gained, 58 lost, 2 retained**. Stops across both sides were mostly exhausted support (192) or multiple supported continuations (70). The method is too restrictive for this graph; it remains opt-in. These results do not establish that the discarded routes are absent: the conservative rule refuses to choose among them.

An isolated representative-selection control reran Huber on all 164 inputs. Every per-input sequence and acceptance decision was identical to v3. Support-first grouping yielded **44 representatives / 8 reference matches**, versus 51 / 8 previously; all 60 matching input records remained. This changes saved representatives, not the underlying predicted paths. Existing reference evidence was reused only after exact sequence verification. The new strict-flow run took 1.48 minutes including initialization, with 2 threads and 2.34 GiB peak sampled RSS.

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

From the updated repository root, build and run the whole server batch:

```bash
git submodule update --init --recursive
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 28
bash scripts/run_copsac_protospacers.sh \
  --root /vol/d/data/copsac \
  --mcaat "$PWD/build/mcaat" \
  --threads 28 --ram-gb 192
```

The runner processes SRR folders containing array reports one by one, retains the graphs and earlier results, and has no wall-clock timeout. Its summary distinguishes candidate regions, final protospacers, candidate-only regions and accepted extensions. `NO_ACCEPTED_EXTENSIONS` means candidates were found but none passed extension filtering; `NO_DETECTIONS` means no candidate was detected. Both are completed outcomes, not failed runs. A changed binary or settings starts a new run rather than reusing an older completed result.

With `--existing-only`, samples without graph metadata are recorded as `SKIPPED_NO_GRAPH`, without creating a sample run, downloading, building or treating the skip as a batch failure. An existing graph with invalid metadata or missing/truncated components still fails validation. Old `protospacer_batch/run_*/discovery` and `extensions` outputs, top-level `ContiguousPaths.fasta`, downloaded reads and `.part` files are not inputs to the new prediction workflow; they remain untouched. Only the selected array report and validated graph are reused.

## Inputs and two distinct modes

In ordinary discovery mode, read the array reports and retain every source-spacer observation. Conservative boundary correction can move shared spacer-end motifs into the repeat: the default requires 85% support among at least five distinct spacers, a motif of at least six bases, and agreement across duplicate and reverse-complement observations. Group corrected sequences with their reverse complements into internal queries while retaining the original source associations.

Map each query to a connected graph path. The defaults require at least 90% identity and 85% query coverage. The mapper permits substitutions and terminal clipping, but not internal gaps. It searches seeds in the succinct index and verifies their connected extensions base by base. Mapping is bounded; a search that reaches a work or placement limit is not exhaustive.

In supplied-target mode (`--protospacers FASTA`), the inputs are already identified targets. Require an exact, full-length graph embedding: 100% identity, 100% query coverage, no clipping. Associate each input with source-array context by substring lookup on both strands; do not replace it with the longer source spacer. A target shorter than graph k is embedded as the exact suffix of a full k-mer. The preceding `anchor_prefix_bp` bases belong to the left flank, not to the supplied target. Saved target coordinates still identify the unchanged original sequence.

Every supplied FASTA record is preserved in `SuppliedProtospacers.fasta` even if it cannot be extended. The final `Protospacers.fasta` contains only those with accepted extensions. Unlike new discovery queries, explicitly supplied targets are allowed to map even when their sequence fails the discovery complexity test; an exact mapped target remains identified in the candidate audit. Flank complexity and technical-adapter checks still apply. Supplied targets with no source-array association receive an explicit unresolved outcome. Supplied-input evaluation and ordinary discovery evaluate different populations and must be reported separately.

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
| `candidate_protospacer_regions` / `CandidateProtospacers.fasta` | All unique detected target sequences, including final targets and detected-only candidates. Diagnostic output. |
| `protospacer_regions` / `Protospacers.fasta` | Unique target sequences with at least one accepted, saved extension. Final protospacer count. |
| `candidate_only_regions` | Candidate target sequences with no accepted extension: candidate count minus final count. |
| `extended_queries` | Queries associated with an accepted extension. Several queries may share one saved region. |
| `unique_extensions` / `ContiguousPaths.fasta` | Unique retained extended regions after consolidation. This is not a count of phage genomes. |
| `limited_queries` | Queries that encountered any recorded search limit. |
| `queries_without_extension` | Canonical queries without an accepted extension, including detected-only and unresolved outcomes. |
| `supplied_targets.tsv` | Every supplied input record mapped to its canonical query, including unsuccessful inputs. |
| `associations.tsv` | Final source → query → target → accepted extension mappings, strands and exact intervals. Every extension ID is populated. Several targets may share a consolidated extension. |
| `candidate_associations.tsv` | All detected source → query → target mappings, including final associations. Empty extension IDs mean no accepted extension for that query. |
| `candidates.tsv`, `unresolved_and_rejected.tsv` | Considered placements and reasons; alternative rows are not extra primary predictions. |
| `target_context.tsv`, `graph_paths.tsv` | Independent discovery context and selected connected graph routes. |
| `coverage_evidence.tsv` | Query ID, target/extension IDs, actual scored edge multiplicities and running minimum count. Query IDs distinguish different input queries that map to the same deduplicated target. The historical `path_capacity` column name no longer means segment-median capacity. |
| `path_selection.tsv` | Retained flank scores, phase, selection status, competing-route count, first divergence and unphased source branching. The legacy `independent_bases` field counts scored edge positions, not statistically independent observations or necessarily distinct nodes. `path_cost = branch_loss / independent_bases` for a nonzero denominator; `runner_up_path_cost` uses that same primary score. `mean_energy` is only a tie-breaker. `branch_observations` is an audit count, not the selection denominator. |
| `withheld_extensions.tsv` | One row per withheld selected route: target, reason, full sequence, block statistics, outward graph nodes and scored multiplicities. These are diagnostic routes, not accepted extensions. `count_positions` are one-based from the graph anchor; add `anchor_prefix_bp` for a short target’s left flank. |
| `query_metrics.tsv`, `run_summary.json` | Per-query work/timings, all parameter values, count meanings and completion status. |
| `search_paths/*.tsv.gz` | Optional finished/pruned path audit. These rows are search states, not protospacer or phage counts. `branch_costs` contains raw costs, including masked positions; its unfiltered sum is not `branch_loss`. Only costs whose preceding and selected positions both occur consecutively in `count_positions` contribute. |

Interval coordinates are zero-based and half-open. Supplied-mode query intervals cover the original complete supplied sequence. Left/right paths are oriented to the canonical query; association strand fields relate them back to the original source. Candidate tables and graph diagnostics retain unsuccessful candidates; use the final FASTAs and `associations.tsv` for downstream analysis. IDs are shared across final and candidate files and may have gaps. Summed worker times can exceed wall time because queries run concurrently.

`first_ambiguous_base` is a zero-based outward position from the graph anchor, interpreted only when `ambiguous=1`. For a short target's left flank, add `anchor_prefix_bp` to express it relative to the target boundary. Coverage-evidence positions are one-based outward added-base positions and already include that padding. The optional search audit's `count_positions` are one-based from the graph anchor and exclude padding; a cost at position `p` is `branch_costs[p-1]`. A `runner_up_path_cost` of -1 means no scored competing completed route was recorded; source-branch uncertainty alone can still set `ambiguous=1`.

## Coverage validation after selection

`--max-coverage-fold 2` validates the one selected extension. Set it to `0` to disable this filter. Search, candidate ranking, spacer detection and source-array handling run first. There is no fallback to a weaker placement, opposite flank or smoother alternative after the selected route fails. This avoids accepting an arbitrary route simply because one of many explored paths passes a smoothness filter.

For configured read length `R` and graph k, use blocks of `max(k+1, R-k)` consecutive scored edges: 127 positions for 150-base reads and k=23. Calculate the median multiplicity of every complete block of the selected flank. With at least two blocks, withhold the extension if the largest median exceeds the smallest by the configured factor. An incomplete trailing block is ignored. Paths with fewer than two complete blocks receive no decision from this check and remain subject to the other extension rules.

Counts from the source-array passage are excluded. Neither repeat multiplicity nor anchor multiplicity sets this baseline. Isolated count spikes or dips generally do not determine a block median. The read-length scale reduces sensitivity to adjacent correlated k-mer counts; it does not make the medians statistically independent.

The original target, source associations and strict discovery context stay intact in the diagnostic outputs. A withheld route is absent from `ContiguousPaths.fasta`; `withheld_extensions.tsv` retains its complete sequence, outward graph nodes and scored counts for inspection. Its candidate has reason `coverage_discontinuity` and an empty extension ID. That query is absent from the final associations; its target is absent from the final FASTA unless another query supplies an accepted extension for the same target. Accepted routes are not shortened by this filter.

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
| `write_outputs`, `write_summary` | Publish final targets with accepted extensions, preserve all candidate evidence separately, and atomically save their distinct counts/status. |
| `self_test` | Exercise arithmetic, context, ranking, ambiguity and output invariants without loading a dataset. |

The batch helper's `validate_predictions` checks completion, final and candidate FASTA counts, the target subset, complete final associations, and the exact target coordinates inside every saved extension. `process` calls it for both new and resumed results before marking a sample complete. Empty final outputs are valid; a final target without a saved extension is an error.

`existing_graph_prefix` searches the original, current-workflow and legacy-batch graph locations in that order. `prepare_graph` validates an existing graph before reuse. `process` records missing-graph skips before creating a run directory when `--existing-only` is set.

## Mathematical background and its limits

If outgoing-base counts are modeled as multinomial observations, `c(e)/sum(c)` is the maximum-likelihood proportion. Applying that model to pooled graph counts and using the proportions to rank a complete route are additional modelling choices: graph counts do not retain which adjacent decisions came from the same read. The score has no fitted weights and does not estimate a sequencing-error distribution.

Squared adjacent differences are a standard smoothness measure; the broader mathematical treatment is described by [Shuman et al. (2013)](https://documents.epfl.ch/users/s/sh/shuman/www/Papers/Journal/Shuman_et_al_SPM_2013.pdf). Here that measure is a descriptive tie-breaker, not an additive likelihood or penalty in the primary score. Overlapping-read coverage is correlated, as discussed by [Evans et al. (2010)](https://arxiv.org/abs/1004.5587); the scored positions must not be counted as independent reads. [Haploflow](https://pmc.ncbi.nlm.nih.gov/articles/PMC8287296/) is related work using local graph coverage for reconstruction, not validation of this specific MCAAT selector.

Local coverage consistency is also used in [MetaCortex (2023)](https://academic.oup.com/bioinformatics/article/39/1/btad020/6986127). Its adjacent-node model is related motivation, not the block-median rule implemented here. [García-López et al. (2015)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4585024/) document coverage variation and reconstruction limitations in viral metagenomes. Neither paper validates MCAAT’s twofold default; its measured trade-off must be reported separately.
