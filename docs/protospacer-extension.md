# Native protospacer workflow (v2)

Experimental: the [full SRR14017553 evaluation](protospacer-evaluation-SRR14017553.md)
found 48/144 target matches and 5/58 extension matches. The mapper improved over the
preceding local prototype; extension quality remains inadequate for a reliable replacement.

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
  --threads 24
```

The output folder must be new or empty. No reads, cycle search or graph rebuilding
are needed for this command. There is no wall-clock timeout. Explicit search budgets
bound work and are recorded; a completed invocation is not necessarily an exhaustive
search. SIGINT/SIGTERM request partial-output saving and return a nonzero exit status.

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

The default mapping budget is 400,000 index/graph operations per query and the
placement cap is 128. Seed searches cap their hit list at 4,096 nodes. Reaching one seed’s hit cap does
not disable the remaining seeds. Once the mismatch allowance is spent, the remaining
seed prefix is found with an exact index lookup. Duplicate connected alignments are
visited once per query start. Every reached limit is reported. Runtime index checks independently recover sampled real graph
nodes using exact and mutated seeds before detection begins.

## Native local-unitig extension

Unbranched paths are compacted lazily in a per-flank cache. A unitig retains actual
arc IDs, bases and all multiplicities. The original graph is never simplified or
rewritten. Every candidate still has a verifiable complete graph-node path.

For each traversed external unitig segment U, its capacity is the median of its
edge counts. The capacity of a candidate path is

`F(path) = min(capacity(U1), capacity(U2), ...)`.

Candidates with larger F are preferred: this is a bounded **widest-path** search.
For example, capacities 20 → 25 → 20 retain capacity 20; a shared repeat's increase
is not treated as a coverage discontinuity. Counts from initial seed/source context
are excluded. A capacity is an abundance-support heuristic, not a correctness
probability, likelihood ratio or independent-read count.

This uses the widest-path idea found in [Haploflow](https://pmc.ncbi.nlm.nih.gov/articles/PMC8287296/).
It is our native, local adaptation. It does not run Haploflow, implement its whole
genome assembler, solve a global minimum-cost flow, or claim exact strain phasing.
The beam and work bounds mean the selected path need not be globally widest.

- Up to 750 added bases on each side, with a beam of 64 candidates by default.
- No return to an anchor node; at most two visits to other nodes.
- Record every external edge count. There is no 100-base sampling gap or hidden
  three-observation minimum.
- A supported flank of 200 bases can qualify. Short dead-end alternatives cannot
  erase a longer eligible extension merely by shortening the comparison horizon.
- For eligible competing routes within the default 1.5 capacity ratio, retain their
  common prefix and mark ambiguity. Prefixes of one route are not separate alternatives.
- Beam/work-truncated candidates are explicitly marked. A limited prediction must
  not be described as uniquely resolved across every possible route.
- If both flanks conflict by more than fourfold capacity or reuse the same external
  nodes, retain the stronger flank rather than joining incompatible contexts.

At most one primary target is selected per corrected query. Equal-quality distinct
target intervals remain unresolved. Exact/RC duplicate extensions are collapsed;
contained extensions point to a longer representative. Partial overlaps are never
joined into an invented sequence. Source-to-target associations remain intact.

## Files and counts

| File | Meaning |
|---|---|
| `Protospacers.fasta` | Unique detected target intervals, with P IDs. |
| `ContiguousPaths.fasta` | Unique accepted extended regions, with E IDs; not a count of complete phage genomes. |
| `sources.tsv` | Every original observation, corrected sequences and trim offsets. |
| `queries.tsv` | Every distinct query, outcome, work count and search-limit flag. |
| `query_metrics.tsv` | Mapping/extension work and per-query mapping, extension and filter seconds. |
| `candidates.tsv` | Examined placements, identity/coverage, acceptance and stop reasons. |
| `unresolved_and_rejected.tsv` | Queries/placements without an accepted extension, including detected but unextended targets. |
| `associations.tsv` | Source/query/target/extension relationships and verified coordinates. |
| `target_context.tsv` | Selected contexts before extension filtering and their original stop reasons. |
| `graph_paths.tsv` | Anchor and retained flank node paths. Left nodes are listed outwards. |
| `coverage_evidence.tsv` | Actual per-edge multiplicities and the associated path capacity. |
| `run_summary.json` | Completion/error status, all counts, parameters, method and stage timings. |

Association coordinates are zero-based, half-open. Query coordinates refer to the
original untrimmed source spacer. Seed/fragment coordinates refer to the saved
extension. `target_strand` relates the original source to that extension, or to the
canonical target FASTA if no extension exists. Left/right fields follow the canonical
query orientation. Containment does not imply support beyond a source's own recorded
fragment. Empty extension fields mean the target was detected without an accepted
extension. Worker times may sum beyond wall time because workers run concurrently.

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
| `group_queries` | Group corrected strands while retaining all observations. |
| `map_spacer` | Seed both portions of covered intervals and construct bounded connected alignments. |
| `source_context` | Orient and test the associated source repeats for a candidate's flank. |
| `build_unitig` | Follow an actual unbranched chain up to a junction, cycle or explicit bound. |
| `update_flow` | Compute a unitig segment's median capacity and update the path bottleneck. |
| `advance_walk` | Apply one arc, anchor/revisit/source checks, and retain count evidence. |
| `extend_flank` | Explore cached unitigs with bounded competing path states. |
| `choose_walk` | Prefer eligible widest paths and preserve ambiguity between comparably supported routes. |
| `trim_walk` | Truncate sequence, nodes and count evidence consistently. |
| `assemble_candidate` | Validate spelling, preserve target evidence, filter flanks and check joint compatibility. |
| `process_query` | Map, extend, select the query's primary result and record function timings. |
| `consolidate` | Assign IDs and collapse exact/RC/contained sequences without losing associations. |
| `write_outputs` | Stage validated FASTA/TSV files and publish them before completion status. |
| `write_summary` | Atomically save status, counts, parameters and wall times. |
| `self_test` | Check strand/context logic, shared-repeat capacity, weak/tied paths and short-dead-end handling. |

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
Missing database matches do not prove false positives, and input reference matches
do not prove their flanks exist in the sample graph. Report search limits and compare
the same source-query set when measuring gains and losses.
