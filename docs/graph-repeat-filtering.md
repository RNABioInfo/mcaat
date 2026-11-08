# Graph-based repeat and spacer filtering

MCAAT's array post-processing now works on the ordered graph cycles returned by cycle enumeration. It does not use SPOA consensus to determine repeat/spacer boundaries, rotate a linear string with duplicated overlap, or merge spacers merely because they share a k-mer.

## Repeat support

Use `--repeat-support 0.85` (the default), or set `repeat-support=0.85` in a settings file. A value of `0.60` requests 60%. Values must be finite fractions greater than zero and at most one.

The denominator is the number of **distinct valid cycles in a candidate group**, initially the cycles collected from one start edge. Rotated copies of the same cycle count once. The numerator is the number of those cycles that contain the **entire same contiguous repeat path**, in the same direction, and pass repeat/spacer length bounds. This is neither nucleotide identity nor prevalence across organisms, reads, or all detected arrays.

For example, 17 of 20 cycles meet 85% support; 12 of 20 meet 60%. The minimum is `max(2, ceil(support * cycle_count))`. At least two distinct spacer sequences are required. Separate 85% majorities at different nodes cannot be combined into a path supported by only 70% of cycles.

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

## Output

`CRISPR_Arrays_N.txt` keeps the legacy text structure. Its repeat is an observed supported sequence; each unique spacer is listed once per candidate, with its observed repeat variant when applicable. Rows are sorted for deterministic output, not biological order.

`CRISPR_graph_systems.tsv` has one row per retained candidate: repeat sequence, support numerator, distinct candidate-cycle denominator, assigned-cycle count, threshold, repeat vertices, source start edges, and number of alternative repeat paths.

`CRISPR_graph_paths.tsv` preserves every assigned cycle observation, including repeat variants, with its repeat/spacer strings, the four vertex categories, ordered complete cycle vertices, ordered actual SDBG edge handles, and whether it supports the selected repeat path. IDs are comma-separated; `-` means an empty path.

`CRISPR_graph_candidates.tsv` records each retained system's candidate cycles and whether boundaries were assigned or unresolved. `CRISPR_graph_summary.tsv` records input, rejection, duplicate and merge counts. These files describe retained candidates; rejected groups are counted in the summary.

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

The command also accepts `--min-repeat-len`, `--max-repeat-len`, `--min-spacer-len`, and `--max-spacer-len`. It validates every k-mer and transition against the graph. Saved cycles must retain their original order and start; hand-rotating individual records changes which start groups they belong to. Normal full-pipeline execution uses in-memory cycle IDs, even when no cycles were found, rather than accidentally reusing a stale `cycles.txt`.

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
