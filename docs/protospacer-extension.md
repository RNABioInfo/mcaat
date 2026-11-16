# Native protospacer workflow (v5)

Experimental: the [complete v5 evaluation](protospacer-evaluation-SRR14017553-v5.md)
on 398 distinct SRR14017553 queries found **20 matching extensions out of 182**, up
from 6/146: **14 new matching sequences, zero lost**. All 146 earlier sequences
remain unchanged; 19/20 matching extensions are at least 750 bp. The separate
81-supplied-target control improves 5/26 to 16/58. The historical 44-match result
has not been fully recovered; the report lists the remaining losses.

This implementation is entirely within MCAAT: one C++ class, `ProtospacerWorkflow`,
and plain structs for options, sources, graph arcs, local unitigs, walks and candidates.
It adds no assembler, aligner, optimization solver or library dependency. The existing
succinct graph is opened read-only and retained. Array detection and CAS code are not
part of this replacement.

## Run on an existing graph

```bash
build/mcaat phage-curate \
  --arrays /vol/d/data/copsac/SRR14017553 \
  --graph /vol/d/data/copsac/SRR14017553/graph/graph \
  --output /vol/d/data/copsac/SRR14017553/native_flow_results \
  --source-array-crossing true \
  --threads 24
```

The output folder must be new or empty. No reads, cycle search or graph rebuilding
are needed for this command. There is no wall-clock timeout. Explicit search budgets
bound work and are recorded; a completed invocation is not necessarily an exhaustive
search. SIGINT/SIGTERM request partial-output saving and return a nonzero exit status.

### Extend supplied targets

To extend a known set of target sequences on the same graph:

```bash
build/mcaat phage-curate \
  --arrays /vol/d/data/copsac/SRR14017553 \
  --graph /vol/d/data/copsac/SRR14017553/graph/graph \
  --protospacers /path/to/known_protospacers.fasta \
  --output /vol/d/data/copsac/SRR14017553/supplied_target_extensions \
  --threads 24
```

This is extension of supplied targets, not a new discovery count. Every supplied
FASTA identifier and sequence is retained in `SuppliedProtospacers.fasta`, including
inputs without an accepted extension. `supplied_targets.tsv` maps input records to
internal queries; query/result tables then identify any accepted extension. Each
supplied sequence is associated with its original array context by substring lookup
on both strands. The full source spacer and unmatched remainder locate the adjacent
repeat; the longer source spacer does not replace the supplied target. There is no
reference-database lookup during this process.

Supplied targets require an exact, full-length graph match: identity and query
coverage are both 100%, with no substitutions or clipping. A target shorter than
graph k is embedded as the exact suffix of a full graph k-mer. Its preceding
padding is explicit in `anchor_prefix_bp`; it is part of the left flank, not part
of the supplied target. Both minimum and maximum flank lengths include that padding.
Saved seed coordinates still identify only the original supplied sequence.

`--source-array-crossing true` is the default. It permits the additional extension
pass described below; it never relaxes strict discovery. Set it to `false` to compare
against the noncrossing extension policy on the same inputs. Keep supplied-input and
full-discovery evaluations separate, even when they use the same graph.

## Detection

1. Read graph-path observations if available; otherwise use numbered array reports,
   falling back to the legacy report. Preserve every source observation.
2. Optionally correct shared spacer-end motifs with at least 85% support among at
   least five distinct spacers. Motifs must be at least six bases; duplicate and
   reverse-complement observations must agree. Record the original boundaries.
3. Group reverse-complement copies of corrected spacers. Reject uninformative queries
   with triplet entropy below 2 or a homopolymer covering at least half the sequence.
4. Partition the minimum covered query interval into short seeds. With at most E
   substitutions among P blocks, at least one block contains at most floor(E/P)
   substitutions. Search these short suffixes in the graph's colexicographic index.
5. From each hit, walk backwards/forwards to align a contiguous graph path to the
   query. Check every base and every edge. Require at least 90% identity and 85%
   query coverage by default. Allow substitutions and terminal clipping; internal
   alignment gaps are not supported. The pigeonhole argument applies before search
   limits and does not make a truncated search exhaustive.
6. Check nearby sequence against the associated source repeat. A continuation of
   that source repeat cannot establish an external target. A shared k-mer alone is
   not globally blacklisted. Require k+1 further bases beyond the initial k/context
   check to call a possible external target. This is a sequence-context criterion,
   not proof of viral origin.

Source-repeat crossing is not evidence for this detection decision. The additional
extension pass runs only after strict discovery has identified a primary target, or
when the user has explicitly supplied the target with `--protospacers`.

The default mapping budget is 400,000 index/graph operations per query and the
placement cap is 128. Seed searches cap their hit list at 4,096 nodes. Reaching one seed’s hit cap does
not disable the remaining seeds. Once the mismatch allowance is spent, the remaining
seed prefix is found with an exact index lookup. Duplicate connected alignments are
visited once per query start. Every reached limit is reported. Runtime index checks independently recover sampled real graph
nodes using exact and mutated seeds before detection begins.

## Native local-unitig extension

The graph is compacted lazily into unbranched segments for each flank. Every
retained sequence follows actual graph arcs. The original graph is not changed.
`path_capacity` is the minimum of the segment-median multiplicities along a path.
A high-count shared segment therefore does not, by itself, reduce this bottleneck.

### One source-repeat passage after detection

In both modes, first run extension with source crossing disabled and keep every
already accepted extension. In discovery mode, a primary candidate must pass strict
target detection before becoming eligible for the additional pass. In supplied-target
mode, the input provides its identity and its entire sequence must match the graph.
Only a target without an accepted extension receives the additional pass, which
can follow one source-repeat passage and continue into the flank. That passage is
marked; it does not establish a newly discovered target.

The additional pass has its own `--max-query-states` allowance. It does not consume
the original pass's remaining budget or replace its accepted extensions. Search
limits still apply and must be recorded. Eligibility to traverse a saved route is
not a guarantee that beam search will retain it or that final filtering will select
it. The repair's reference-match counts must be measured on completed runs.

### Local support and coverage variation

The retained extension algorithm measures counts beyond the seed/source-array
context. `update_profile` forms a robust continuity trace from segment medians:

1. A segment with multiple incoming and outgoing connections is a possible shared
   segment. Exempt its elevated count from the continuity trace only when ordinary
   entry and exit segments are both observed and no higher. Low counts are not exempted.
2. Repeat each remaining segment median according to its number of recorded edges,
   and group this trace into windows of k+1 positions (24 for the tested graph).
3. Take each window's median, including a final partial window if at least half full.
   Record the largest adjacent-window change in `log2(median + 1)`.
4. Prefer continuity within `log2(branch_ratio)` of the best observed value. Within
   that group, flows within a factor `branch_ratio` of the greatest flow are comparable.
   Prefer length within the best group, then flow, then lower deviation.

The default ratio 1.5 is a heuristic tolerance. This is a bounded widest-path and
continuity heuristic, not a fitted count model or a probability of correctness.
Overlapping k-mer counts are dependent. The archived calibration sample contained
only one qualifying long unbranched segment, so it did not support a meaningful
negative-binomial fit. A tested alternative using 150-base windows and average log
support shortened previously matching extensions and was not retained as the default.

Local relative coverage has precedent in
[Haploflow](https://link.springer.com/article/10.1186/s13059-021-02426-8).
The score, tolerances and beam policy here are MCAAT choices. See
[sequencing evidence and limitations](protospacer-sequencing-evidence.md).

### Technical sequences and low complexity

Coverage continuity alone can prefer a stable artifact. By default,
`--filter-extension-artifacts true` adds two native string checks during traversal:

- Stop before an exact 18-base prefix of either [Illumina TruSeq trimming
  sequence](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/CDIndexes.htm),
  on either strand. Check across the anchor/flank boundary too. Queries or mapped
  targets containing these cores are recorded as `adapter_like_query` or
  `adapter_like_target`; their source observations remain in the output.
- Inspect a 64-base window every eight added bases. Stop before that window if its
  triplet entropy is below 2 bits or a homopolymer occupies at least half the window,
  using the same complexity predicate already applied to input queries.

The clean preceding prefix can still qualify; other competing branches continue.
These checks prevent extension through known technical motifs or uninformative
sequence. They do not establish viral origin. The motif policy is explicit and
can be disabled with `--filter-extension-artifacts false`.

Comparable routes that diverge restrict output to their common prefix, with
`ambiguous_coverage` recorded. A shorter prefix of the same route is not treated
as a competing biological path. The implementation does not silently choose an
unsupported long route merely to reach the requested maximum length.

### Search and work limits

- Up to 750 added bases per flank, beam 64, and at most two visits to external nodes.
  Source-anchor revisits remain prohibited. A source-array continuation cannot
  establish a target; one passage is available only in the additional extension pass.
- Expand cached unbranched segments in rounds, retaining a single bounded live
  frontier. Segment lengths can differ; the work cap counts their actual arcs.
  A tested distance-bucket schedule with early bubble pruning lost known routes
  earlier and was withdrawn. The retained search uses the earlier continuity
  and flow ranking throughout.
- Reserve a beam slot for each distinct combination of the first two branch choices
  across retained states. This retains lower-abundance alternatives rather than using
  every slot for similar descendants. Remaining slots follow the ranking above.
- Keep at most four beam widths of terminal candidates. Filling this buffer trims
  the retained candidate list; **it does not terminate the live frontier**. This fixes
  the earlier stop at roughly 226 bases observed while a strongly supported path
  was still advancing. Work, flank-length and visit limits still apply.
- Beam-pruned prefixes are diagnostic records, not completed candidate extensions.
  A work-limited selected prefix remains explicitly marked as limited.
- Visit mapping tiers in decreasing query coverage and identity. Share a tier's
  budget across its placements, reserving half the remaining work for lower tiers
  when they exist. Stop after a tier produces an accepted extension, or the query
  budget is exhausted. Merely detecting a target without extending it no longer
  suppresses the next tier. The preceding best tier keeps its original search budget.
  Mapping itself still has explicit state and placement limits.
- A flank of at least 200 bases may qualify. Counts below the initial source/seed
  context do not define external support. There is no sparse 100-base sampling rule.
- The additional crossing pass applies only to a primary detected but unextended
  candidate, or to an explicitly supplied target. It receives a separate query-work
  allowance. Disabling it does not change strict detection or mapping thresholds.
- If the two flanks conflict by more than fourfold path capacity or reuse external nodes,
  retain the greater-capacity flank instead of joining incompatible contexts.

At most one primary target is selected per corrected query. Prefer an accepted
extension over a detected-only placement; then prefer query coverage and alignment
identity. An already extended better alignment stays preferred. Distinct targets
with identical alignment quality in the same outcome category remain unresolved.
No abundance score overrides the spacer alignment. Exact/RC duplicates and contained
extensions are consolidated; partial overlaps are not joined into an invented
sequence. Strict-pass extensions are consolidated first. Their saved representative
sequences are preserved when new crossing-pass extensions are added; a longer new
extension cannot replace an already accepted strict representative. This also avoids
changing a baseline sequence's reference-coverage denominator solely through
containment consolidation. Every original observation and its association are preserved. Search-limit
flags mean the result is conditional on the explored paths, not globally unique.
In supplied-target mode, also preserve the original input records and their mapping;
an accepted extension is optional and must not determine whether an input is retained.

### Optional path audit

`--save-search true` streams every finished or pruned path to
`search_paths/Q....tsv.gz`, including its sequence, nodes, counts, unitig blocks,
stop reason and coverage deviation. Records are written as the search runs; discarded
paths need not remain in RAM. This is optional diagnostic output and can be large.
It is disabled in ordinary and batch runs.
Path identifiers include the strict or crossing phase, so rerunning extension for
one placement does not reuse an identifier within that query's audit file.

## Files and counts

| File | Meaning |
|---|---|
| `Protospacers.fasta` | Unique detected target intervals, with P IDs. |
| `ContiguousPaths.fasta` | Unique accepted extended regions, with E IDs; not a count of complete phage genomes. |
| `SuppliedProtospacers.fasta` | In supplied-target mode, every original input identifier and sequence, whether extended or not. |
| `supplied_targets.tsv` | In supplied-target mode, original input IDs and sequences mapped to internal query IDs. |
| `sources.tsv` | Discovery: original array-spacer observations, corrected sequences and trim offsets. Supplied mode: source-context associations for the supplied targets. |
| `queries.tsv` | Every distinct query, outcome, work count and search-limit flag. |
| `query_metrics.tsv` | Mapping/extension work and per-query mapping, extension and filter seconds. |
| `candidates.tsv` | Examined placements, identity/coverage, acceptance, stop reasons, left/right coverage deviations, crossing flags and prior target identification. Alternatives are not primary predictions. |
| `unresolved_and_rejected.tsv` | Queries/placements without an accepted extension, including detected but unextended targets. |
| `associations.tsv` | Source/query/target/extension relationships and verified coordinates. |
| `target_context.tsv` | Selected contexts before extension filtering and their original stop reasons. |
| `graph_paths.tsv` | Anchor and retained flank node paths. Left nodes are listed outwards. |
| `coverage_evidence.tsv` | Actual per-edge multiplicities and the associated path capacity. |
| `run_summary.json` | Completion/error status, counts, discovery/supplied mode, crossing policy, parameters, method and stage timings. |

In discovery mode, `spacer_observations` counts original array-spacer observations.
In supplied-target mode, it counts source-context associations, which can exceed
the number of input records because one target can occur in several source contexts.
The exact supplied-input count is `supplied_targets` in the summary; the unchanged
records are in `SuppliedProtospacers.fasta` and `supplied_targets.tsv`.
`distinct_spacers` counts internal canonical queries, so duplicate or reverse-complement
inputs can share one query while remaining separate supplied records.

Association coordinates are zero-based, half-open. Query coordinates refer to the
original untrimmed source spacer. Seed/fragment coordinates refer to the saved
extension. `target_strand` relates the original source to that extension, or to the
canonical target FASTA if no extension exists. Left/right fields follow the canonical
query orientation. Containment does not imply support beyond a source's own recorded
fragment. Empty extension fields mean the target was detected without an accepted
extension. Worker times may sum beyond wall time because workers run concurrently.
For short supplied targets, `anchor_prefix_bp` distinguishes graph-anchor padding
from the target. Reported left-flank lengths and coverage-evidence positions are
measured from the target boundary and include that padding. Supplied-mode query
coordinates cover the complete supplied sequence.

## Every C++ function

| Function | Responsibility |
|---|---|
| `ProtospacerWorkflow` | Own options, graph, source observations and query results. |
| `main` | Dispatch help/self-test, initialize output, install signals and save failure status. |
| `parse_options` | Parse the explicit native CLI and validate inputs, budgets and output safety. |
| `print_help` | Describe options and limits. |
| `interrupt` | Record an interruption for bounded workers to finish/checkpoint. |
| `run` | Execute stages, parallelize queries, propagate errors and save completion status. |
| `reverse_complement` | Reverse/complement strict ACGT DNA. |
| `canonical` | Select the lexical sequence/RC minimum for deduplication. |
| `dna` | Check nonempty ACGT input. |
| `json` | Escape output strings safely. |
| `numbered` | Produce readable deterministic Q/P/E identifiers. |
| `median` | Calculate unitig capacities from count vectors. |
| `triplet_entropy` | Measure sequence complexity from overlapping nucleotide triplets. |
| `informative` | Apply entropy and long-homopolymer checks. |
| `adapter_position` | Find exact trimming cores in the required traversal orientation. |
| `repeat_context` | Edit-compare an adjacent repeat; distinguish pending, matching and incompatible context. |
| `label` | Read a real k-mer label and reject dummy labels. |
| `find_node` | Retrieve and verify an exact k-mer index hit. |
| `neighbors` | Return actual incoming/outgoing arcs with their multiplicities. |
| `added_base` | Read the nucleotide contributed by a forward or backward arc. |
| `spell` | Verify adjacency and spell a connected graph-node path. |
| `search_range` | Search a colexicographic suffix interval within a substitution/work budget. |
| `find_starts` | Retrieve and deterministically order short-seed hits. |
| `verify_index` | Check exact/mutated seed searches against sampled known graph labels. |
| `load_sources` | Read reports, expand repeat variants and validate source identities. |
| `correct_boundaries` | Infer shared-end motifs and record conservative boundary adjustments. |
| `prepare_source_nodes` | Index reported source-array k-mers on both strands for extension checks. |
| `load_supplied_targets` | Preserve explicit FASTA inputs and associate them with source contexts by substring lookup on both strands. |
| `group_queries` | Group corrected strands while retaining all observations and supplied-input associations when applicable. |
| `map_spacer` | Construct bounded connected alignments; require the full exact sequence for supplied targets and explicitly pad short targets to graph k. |
| `source_context` | Orient and test the associated source repeats for a candidate's flank. |
| `source_clearance` | Bound the source-overlap interval using the unmatched spacer remainder, repeat length, allowed edits and graph k. |
| `build_unitig` | Follow an actual unbranched chain up to a junction, cycle or explicit bound. |
| `update_flow` | Record a unitig block and diagnostic bottleneck, then advance its coverage profile. |
| `update_profile` | Track the largest adjacent change between windows of segment-median counts, exempting bracketed shared increases. |
| `select_candidate` | Prefer an extended placement, then alignment quality; retain one primary result or record equal-quality ambiguity. |
| `rank_walks` | Rank continuity/support bands and preserve distinct branch families in the beam. |
| `retain_completed` | Remove ineligible terminal paths, then bound their storage without stopping live paths. |
| `save_walks` | Stream diagnostic paths and evidence to per-query gzip files with unique IDs. |
| `advance_walk` | Apply one arc, anchor/revisit/source checks and count evidence under the strict or permitted extension-only crossing policy. |
| `extension_barrier` | Stop before a technical motif or uninformative window, keeping graph coordinates consistent. |
| `extend_flank` | Explore cached unitigs with bounded competing path states and the selected source-crossing policy. |
| `choose_walk` | Prefer length within comparable continuity/support and retain only a common prefix when routes remain ambiguous. |
| `trim_walk` | Truncate sequence, nodes and count evidence consistently. |
| `assemble_candidate` | Validate spelling, preserve target evidence, filter flanks and check joint compatibility. |
| `process_query` | Run noncrossing extension first in both modes, preserve accepted extensions, and give an eligible unextended primary target a separately budgeted crossing pass; record per-function timings. |
| `consolidate` | Consolidate strict-pass representatives first, then recovered sequences, preserving accepted strict sequences and all associations. |
| `write_outputs` | Stage validated FASTA/TSV files and publish them before completion status. |
| `write_summary` | Atomically save status, counts, parameters and wall times. |
| `self_test` | Check strand/context logic, branch ranking, partial coverage updates, artifact boundaries and bounded terminal storage. |

## Evaluation, separate from prediction

`tests/evaluate_protospacer_workflow.py` checks frozen outputs and can optionally
compare them with a local reference using BLAST. BLAST is only an evaluation tool;
it is not called by MCAAT or required for prediction. The helper permits only the
three requested dataset labels: SRR14017553, SRR14017554 and SRR14017569.

| Function | Responsibility |
|---|---|
| `fasta` | Read strict DNA FASTA and reject duplicate/empty records. |
| `table` | Read TSV and reject inconsistent columns. |
| `reverse_complement` | Independently orient sequences for association checks. |
| `digest` | Fingerprint files with SHA-256. |
| `write_json` | Save readable evaluation records. |
| `validate` | Check all inputs/outcomes, identities, coordinates, retained evidence and deduplication. |
| `blast` | Evaluate frozen sequences, saving the exact command and runtime. |
| `summarize_matches` | Calculate exact identity/coverage and save every sequence's best qualifying coverage. |
| `main` | Validate, freeze/fingerprint, evaluate targets/extensions/all inputs and verify unchanged predictions. |

Default reference criteria: 90% identity with 85% target coverage or 50% extension
coverage in one HSP. These are reference-match rates, not biological precision.
Targets and input spacers use 7-base BLAST seeds. Extensions use 9-base seeds when
every sequence is at least 180 bases: at the specified thresholds, a qualifying
alignment then contains an exact run of at least nine bases. BLAST remains a
heuristic alignment search. Use the same settings for before/after comparisons.
Missing database matches do not prove false positives, and input reference matches
do not prove their flanks exist in the sample graph. Report search limits and compare
the same source-query set when measuring gains and losses.

## Diagnostic helpers

`tests/audit_protospacer_search.py` evaluates saved candidates only after prediction.
It uses SQLite for deduplication and records shorter prefixes separately from full
750-base flanks. It seeks one qualifying reference witness per query; successful
queries need no further searches. For unsuccessful queries, it searches every saved
full-length candidate with word-size-16 BLAST. This gives confirmed witnesses, not
a maximum recoverability estimate. An optional `--sensitive` pass exhausts the
remaining pool with word size 9; repetitive paths can make this very expensive.
The report explicitly records whether that additional pass completed. Neither
mode proves absence of alternative routes in the full graph. Batches rotate among
placements and flank sides. All commands and hit files are retained.

| Function | Responsibility |
|---|---|
| `reverse_complement` | Normalize candidate strands. |
| `prepare` | Stream gzip records into a disk database, count partial paths and export eligible distinct sequences. |
| `screen` | Find positive witnesses with longer exact seeds, stopping each query once one is confirmed. |
| `evaluate` | Search unresolved queries in batches, preserve witnesses and record how many candidates were tested. |
| `main` | Parse paths and start/resume the separate evaluation. |

`tests/make_protospacer_controls.py` creates one shuffled control for each corrected
query, preserving length and exact dinucleotide counts. It keeps repeat associations,
uses a recorded random seed and rejects exact copies of real/control queries. Inputs
are already boundary-corrected, so run them with `--correct-spacer-boundaries false`.
These are background controls, not certified biological negatives or a calibrated FDR.

| Function | Responsibility |
|---|---|
| `reverse_complement` | Preserve orientation and reject duplicate controls on either strand. |
| `shuffle_dinucleotides` | Randomize an Euler trail through nucleotide-pair edges and verify exact pair counts. |
| `main` | Generate controls and their original-query mapping, preserving source observations. |

`tests/compare_protospacer_runs.py` compares two completed, evaluated runs. It first
requires identical query IDs and sequences, then maps reference matches back through
`associations.tsv`. It reports both unique saved sequences and input queries, with
explicit gained/lost/retained query lists, extension lengths, reference-covered bases,
and whether each previous matching sequence was actually preserved. It saves a row
for every query. It never chooses
prediction paths.

| Function | Responsibility |
|---|---|
| `table` | Read tab-separated records. |
| `load_run` | Require completed predictions and resolve saved matches to query IDs. |
| `main` | Check identical cohorts and save exact counts and per-query comparisons. |
