# Graph-based repeat and spacer filtering

MCAAT's array post-processing now works on the ordered graph cycles returned by cycle enumeration. It does not use SPOA consensus to determine repeat/spacer boundaries, rotate a linear string with duplicated overlap, or merge spacers merely because they share a k-mer.

## Repeat support

Use `--repeat-support 0.85` (the default), or set `repeat-support=0.85` in a settings file. A value of `0.60` requests 60%. Values must be finite fractions greater than zero and at most one.

The denominator is the number of **distinct valid cycles in a candidate group**, initially the cycles collected from one start edge. Rotated copies of the same cycle count once. The numerator is the number of those cycles that contain the **entire same contiguous repeat path**, in the same direction, and pass repeat/spacer length bounds. This is neither nucleotide identity nor prevalence across organisms, reads, or all detected arrays.

For example, 17 of 20 cycles meet 85% support; 12 of 20 meet 60%. The minimum is `max(2, ceil(support * cycle_count))`. Final accepted candidates require at least two distinct observed spacer sequences. Strict mode additionally requires two spacer families (defined below). Separate 85% majorities at different nodes cannot be combined into a path supported by only 70% of cycles.

The algorithm enumerates observed contiguous paths meeting this support, retains paths that cannot extend further with sufficient support, and applies the configured length bounds. It chooses the longest qualifying path, then the highest support, then a deterministic graph-ID tie-break. The number of additional eligible repeat paths is reported. It retains the existing minimum median spacer/repeat length ratio of 0.5. It does not crop an overlong conserved path to make it pass the maximum repeat length.

Lowering support can extend boundaries or merge candidate groups, so the final number of arrays or spacers is **not necessarily monotonic**.

## Node identity and boundaries

SDBG stores outgoing edge slots. Several slots can represent transitions from the same k-mer. Post-processing uses a normalized **k-mer vertex identity** to recognize shared repeat nodes, and saves the actual traversed **edge handles** separately. This also handles a repeat of exactly k bases, whose shared repeat path has one vertex but potentially different outgoing edges.

For each observed circular path, the output separates:

- k-mers wholly inside the repeat;
- k-mers wholly inside the spacer;
- k-mers crossing from repeat to spacer;
- k-mers crossing from spacer back to repeat.

A k-mer crossing a boundary is not mislabeled as a pure spacer. If a spacer is shorter than k, it has no pure spacer k-mer; boundary categories can overlap. Its nucleotide sequence and full cycle remain available.

Minority repeat variants retain their **observed** sequence when both repeat boundary anchors are present on their cycle and lengths pass the filters. Cycles without those anchors remain explicitly unresolved. No consensus sequence is invented to fill them in. Repeats shorter than k cannot be isolated as repeat-only k-mer paths by this method.

## Deduplication

Candidate comparisons require at least one shared complete cycle, with rotations and reverse complements accounted for. The candidates must share repeat nodes and must not assign one another's repeat nodes to pure spacers.

- Duplicate and contained detections collapse into one candidate.
- Partly overlapping detections can merge, adding their distinct cycles, only if repeat inference at the same support threshold succeeds and preserves all previously assigned cycles without role conflicts.
- A shared repeat alone, or a shared spacer associated with a different repeat, does not trigger a merge.
- Original start-edge IDs are retained as provenance.

The cycle-owner index is built once. Removing a contained detection removes only
that detection's entries; merging or replacing a system reindexes only the
changed system. Comparisons follow the original start-group order, and changed
systems reschedule earlier neighbors when necessary. Repeat/spacer roles and
canonical cycle sets are cached between changes. Candidate pairs are generated
for one system at a time instead of allocating every pair sharing a cycle.

This removes the previous global index rebuild after every duplicate removal.
It does not skip repeat inference, reduce the support threshold, or cap candidate
comparisons. Incompatible overlapping candidates can still require many
comparisons. Filtering is currently serial; `--threads` does not parallelize this
stage. See [the performance checks](repeat-filter-performance.md) for measured
scope and reproducible commands.

These rules consolidate graph evidence. A pooled graph cannot establish that two identical repeat–spacer combinations came from the same physical locus. Nor does this output infer genomic spacer order. Additional read or assembly linkage is needed for those claims.

## Cycle length range

The default range is **47–100 inclusive**, measured as distinct k-mer nodes on
the simple cycle, counting the returning start only once. Both CLI and settings
accept `cycle-min-length` and `cycle-max-length`. Enumeration includes the stated
minimum (the older implementation incorrectly used `>`), and enforces the
maximum even when returning to the start.

`postprocess-arrays` applies these bounds to saved input as well, before repeat
inference, and reports excluded cycles. It cannot discover cycles missing from
an earlier search. If the saved search stopped at 77, rerun ordinary MCAAT with
`--graph`, minimum 47 and maximum 100; do not rebuild the graph or expect a
post-processing replay to supply lengths 78–100.

## Final quality and competing assignments

`--array-filter-mode conservative` is now the default. It retains candidates
passing the structural, length and 85% support checks, and reports sequence
heuristics as `quality_flags`. Similarity-based grouping annotates the variant table;
it does not assert that nearby sequences are errors or delete candidates with
two distinct observations in one display family. Within-candidate strand copies
still count once. Original graph-compatible deduplication/merging still runs.

Shared cycles across remaining candidates are preserved and explicitly listed
in `CRISPR_graph_shared_cycles.tsv`; both candidates get a `shared_cycle_evidence`
flag. Alternative partitions of the same circle are not proven false merely
because another candidate was ranked first. Large merged candidates are flagged
and quadratic grouping is skipped: each observed string is annotated as its own
family in that case. The enumeration cap itself is unchanged.

The numbered rejection/ownership procedure below describes **opt-in strict mode**
(`--array-filter-mode strict`), introduced in commit `21f8f8d`. Its thresholds
were not calibrated against the complete annotated metagenome, so it must not
be the default used for sensitivity assessment.

After repeat inference and compatible merges, candidates undergo an explicit
curation pass. These are configurable heuristics, not probabilities or a trained
CRISPR classifier. The 85% repeat-support definition above is unchanged.

Sequence similarity is `1 - global_Levenshtein_distance / max(lengths)`, with
unit costs for substitutions, insertions and deletions. It is not BLAST identity.
A Myers bit-vector implementation handles normal spacer lengths; longer
configured sequences use dynamic programming.

1. Sort unique observed spacer strings lexically. Group them deterministically
   so **every pair** within a family meets `--spacer-variant-similarity` (default
   0.90). Similar A–B and B–C alone cannot bridge dissimilar A–C. Choose an
   **observed medoid**, with lexical tie breaking, as the display representative.
   A family can contain sequencing errors or real alleles; neither is inferred.
   Set the threshold to 1 to keep distinct strings in separate families.
2. Require at least two families. Reject candidates whose mean pairwise
   similarity between family representatives exceeds `--max-spacer-similarity`
   (0.60), or whose mean repeat–representative similarity exceeds
   `--max-repeat-spacer-similarity` (0.60). Setting either maximum to 1 disables
   that rejection. The median spacer/repeat length ratio must still be at least
   0.5. Candidates with at least `--max-cycles` distinct candidate cycles (500)
   are rejected before allocating a pairwise distance matrix. This limit also
   applies after merging starts.
3. Rank passing candidates by more spacer families, lower mean inter-family
   similarity, longer repeat path, then original candidate order. Assign each
   **complete circular sequence**, normalized over rotations and both strands,
   to at most one accepted candidate in strict mode. This is a deterministic evidence ownership
   convention, not proof that a chosen boundary or physical locus is correct.
4. Remove already-owned cycles from each later candidate. Opposite-strand
   copies inside one candidate count once too; their alternative paths remain
   in review, and an observation supporting the selected repeat is preferred. Recheck 85% support
   on the residual denominator, at least two families, and the quality limits.
   Unresolved cycles remain in the denominator unless their exact cycle is
   already owned. Keep independently supported residual evidence; otherwise
   retain the competing interpretation only in the review files. A rejected
   candidate never claims cycles. Do not re-infer a new boundary to force a
   residual candidate to pass.

High spacer similarity is an established warning sign in CRISPR detection (for
example [CRISPRdb](https://link.springer.com/article/10.1186/1471-2105-8-172) and
[CRISPRidentify](https://pmc.ncbi.nlm.nih.gov/articles/PMC7913763/)). MCAAT's global
edit-distance metric and family grouping are explicitly defined here; these
thresholds do not reproduce those tools' alignment scores or classification.
The retained variants and rejection reasons allow the user to inspect the effect
of changing the thresholds.

During cycle enumeration, reaching the exclusive cycle limit now **returns
immediately**, discarding that start's cycles without marking them visited.
Previously clearing the results and leaving the closing edge on the search
stack could repeat the same work until the ten-million-step watchdog.

## Output

`CRISPR_Arrays_N.txt` keeps the legacy text structure. Its repeat is an observed supported sequence; every distinct observed spacer is displayed with its observed repeat variant when applicable. `spacers=` counts observed strings and `spacer_families=` counts similarity groups. Rows are sorted for deterministic output, not biological order.

`CRISPR_graph_systems.tsv` has one row per retained candidate: repeat sequence, support numerator, distinct candidate-cycle denominator, assigned-cycle count, threshold, repeat vertices, source start edges, number of alternative repeat paths, original candidate ID, spacer-family count and sequence-variant count.

`CRISPR_graph_paths.tsv` preserves every assigned cycle observation, including repeat variants, with its repeat/spacer strings, the four vertex categories, ordered complete cycle vertices, ordered actual SDBG edge handles, and whether it supports the selected repeat path. IDs are comma-separated; `-` means an empty path.

`CRISPR_graph_candidates.tsv` records each retained system's candidate cycles and whether boundaries were assigned or unresolved. `CRISPR_graph_summary.tsv` records input, rejection, duplicate and merge counts. These files describe retained candidates. `spacers` and `sequence_variants` in the summary count distinct observed strings
summed per candidate; `spacer_families` and `assigned_observations` are separate.
Groups rejected during initial inference are still only counted in the summary.

`CRISPR_spacer_variants.tsv` maps every retained observed spacer string to one
family and an observed representative. CAS/phage graph-path readers still receive
all variants, as does the main text report.

`CRISPR_graph_shared_cycles.tsv` lists normalized circular sequences shared by retained candidates and their IDs. Complete overlap membership is stored once per cycle rather than as a quadratic set of candidate pairs; `competing_candidates` in the decision table lists exclusive-assignment competitors only in strict mode. `unique_circular_sequences` and `unique_spacer_sequences_across_arrays` in the summary separate globally deduplicated evidence counts from per-candidate totals. These are sequence counts, not inferred physical loci.

`CRISPR_graph_decisions.tsv` records every candidate entering final curation,
its original counts and sequence metrics, status, competing candidate IDs,
retained observation count and accepted array ID. Metrics not computed after an
early rejection are `NA`. Family counts are zero when grouping was skipped at
the cycle limit. Candidate IDs refer to this table; array IDs are assigned only
to accepted candidates. The `quality_flags` column describes retained evidence when accepted, and the original candidate otherwise. Original metrics describe the candidate before shared
cycles were removed; retained counts describe the accepted residual candidate.

`CRISPR_graph_alternatives.tsv` preserves every rejected or removed observation
from final curation, with its proposed repeat/spacer assignment and actual graph
paths. It is not a second list of accepted arrays. Every observation entering
curation appears exactly once in either this file or the accepted path table.
`CRISPR_graph_review_cycles.tsv` also preserves non-retained candidate paths,
including paths whose repeat/spacer boundaries could not be assigned. Review
and accepted candidate-cycle files can both mention an unresolved retained path;
the two **observation** tables remain disjoint by original candidate/observation.
Initial inference failures and discarded enumeration tangles are not reconstructed
in these review files.

The CAS and phage directory readers prefer `CRISPR_graph_paths.tsv`. They validate edge continuity, source labels, sequences and boundary partitions against the loaded graph before passing actual edge handles downstream. A mismatched or damaged sidecar fails explicitly. Legacy output directories without a sidecar still use the old text reader. The C++ `getGraphSystems()` accessor preserves candidate identity; the older `getSystems()` map is a sequence-only compatibility view that pools entries with identical repeat strings.

## Re-filter without running cycle enumeration again

Retain the graph and cycles with `--autoclean false`, then run:

```bash
mcaat postprocess-arrays \
  --graph /data/graph/graph \
  --cycles /data/run/cycles.txt \
  --output-folder /data/run-support-60 \
  --repeat-support 0.60
```

The command also accepts `--array-filter-mode`, `--cycle-min-length`, `--cycle-max-length`, all four quality options described above, `--min-repeat-len`, `--max-repeat-len`, `--min-spacer-len`, and `--max-spacer-len`. It validates every k-mer and transition against the graph. Saved cycles must retain their original order and start; hand-rotating individual records changes which start groups they belong to. Normal full-pipeline execution uses in-memory cycle IDs, even when no cycles were found, rather than accidentally reusing a stale `cycles.txt`.

Saved-cycle loading looks up the first k-mer of each record, checks each textual
overlap, and follows the corresponding outgoing edge for subsequent k-mers.
It also checks cycle closure. It no longer performs a global graph lookup for
every k-mer. Start lookups use a bounded cache, and validation of repeated or
rotated cycles is cached within one filter run. Malformed records report their
line number; I/O failures abort the run.

Stage messages distinguish graph loading, cycle loading, repeat inference,
candidate indexing, duplicate comparisons, and output writing. Long phases
report progress approximately every five seconds at record/group boundaries.
Loading reports bytes, records, and k-mers processed. No stage silently switches
to a different support threshold or skips records to finish sooner.

To rerun **cycle search** on a saved graph instead of reusing saved cycles, use
the ordinary pipeline with `--graph`:

```bash
mcaat --graph /data/run/graph/graph --threads 24 \
  --repeat-support 0.85 --autoclean false --output-folder /data/run-restarted
```

This starts from the graph, runs cycle enumeration again, and filters the new
in-memory results. The `postprocess-arrays` command instead reads saved cycles.

## Validation and limits

Synthetic regression tests cover threshold boundaries, duplicate cycles, observed repeat variants, unresolved anchors, exact-k repeats, short spacers, incompatible systems, overlapping candidates, reverse complements and shifted repeat boundaries. The real-SDBG integration test validates output paths, downstream edge handles, deterministic reruns, stale output cleanup and rejection of altered provenance.

This changes the definition of post-processed candidates, so old SPOA output is not expected to match. Graph-path validity is testable, but it is not proof of biological accuracy: shared spacer prefixes/suffixes can move inferred repeat boundaries, incomplete cycle enumeration affects support, and multiple eligible repeat paths can remain ambiguous. Benchmark precision and recall against labeled arrays before drawing biological accuracy conclusions.


See [curation validation](array-curation-validation.md) for the local reference
comparison, regression scope and commands. The report replay helper operates
only on previously retained candidates; it cannot recover cycles or candidates
that were not included in those reports.
