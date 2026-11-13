# Protospacer detection and bounded target extension

This implementation is on `codex/protospacer-extension`, based on
`origin/master` at `4f75041`. It replaces the production `phage-curate`
implementation; `protospacer-curate` is an alias. It does not incorporate the
separate array-filter branch. The old experimental classes remain in the
repository for reference, but the release subcommand uses the new module.

## Extend previously accepted protospacers without redetecting them

Use `--protospacers FILE` to freeze an existing FASTA as input to extension:

```bash
build/protospacer/mcaat phage-curate \
  --arrays /path/to/array-results \
  --graph /path/to/previous-run/graph/graph \
  --protospacers /path/to/accepted/Protospacers.fasta \
  --output /path/to/separate-extension-results \
  --threads 24 --max-flank-length 750 --min-extension-flank 750 --max-length 1550
```

Every supplied ID and sequence is retained in `Protospacers.fasta`, including
unmapped seeds, short seeds and failed extensions. Sequence letters are normalized
to uppercase; IDs and input order remain unchanged. Supplied sequences are already
the desired protospacers: boundary correction applies to source-array context,
**not** to the supplied FASTA. This mode does not rediscover, reference-filter,
shorten, deduplicate or replace its inputs. It is graph-only; read/PAM options
are rejected explicitly. Use a new output directory or one from the same mode;
it cannot overwrite the input FASTA or mix with a discovery run.

Source-array observations identify the complete spacer and repeat around each
supplied interval, including reverse orientation and original source coordinates.
The default quality filter requires the entire input sequence exactly, reranks
the resulting paths, and withholds unresolved competing routes. The
[quality rules below](#reference-free-quality-filtering-of-supplied-protospacers)
describe the current acceptance policy.

With `--path-quality-filter false`, the former workflow tries exact placement
first and preserves its selected extension. Otherwise,
`--target-variant-fallback true` searches the graph
for related target sequences. Its aligned interval must cover at least 85% of
the original query at at least 90% substitution identity. Original FASTA records
remain unchanged; a differing observed target is reported separately, with
its query coordinates, identity, coverage, and `extended_variant` status.

Short queries of at least 15 bases can be placed inside actual graph k-mers.
The surrounding bases come from graph labels and edges; they are neither
invented padding nor counted as matched query sequence. Complete source-array
context must be available around the chosen core. The BOSS suffix lookup uses
the existing colexicographic ordering without rebuilding an index or graph.

Short-target placement first probes up to 75 bases on each side of every candidate, then
extends the best N placements using the existing multiplicity score and source,
cycle, adapter and depth guards. It shares the existing per-query work budget.
`--max-target-placement-states 50000` separately bounds placement work;
`--max-seeds` caps retained placements. Work limits and unextended alternatives
are explicit. With the quality filter disabled, if these anchors produce no extension and budget remains, queries
of at least k bases receive a masked BOSS lookup. This searches beyond the
one-substitution anchor restriction: a valid 90%/85% target may have two
substitutions in every complete k-mer. Candidates must still pass the full
query-interval check and spell actual graph paths. The rescue shares both work
budgets with the earlier attempts and reports `masked_probe`/`masked_extend`
checks. It never replaces an existing exact or first-fallback extension.

Only this final rescue uses a counting-noise tolerance in beam and final profile
ranking. From positive edge counts, estimate the zero-truncated Poisson rate
whose conditional mean equals the observed mean. A delta-method variance for
`log2(1 + count)` supplies a tolerance for within-flank variance and adjacent
changes; the expected fraction of singleton counts supplies the corresponding
singleton tolerance. Excesses remain penalized, as do finite support and residual
between-flank differences. This allows ordinary variation between counts one
and two without automatically preferring a short flat dead end. It is a ranking
heuristic for correlated k-mers, not a fitted probability of phage origin.
Existing predictions retain their original scoring. The policy is recorded as
`masked_rescue_stability` in schema 12 parameters.

The production mapper remains bounded and does not allow target indels. It cannot prove
which organism supplied a short sequence shared by several graph locations.
Use `--target-variant-fallback false` to retain the earlier exact-placement policy.
Missing source context, adapter-containing seeds and work limits receive
separate reasons.

At branches the bounded beam keeps alternatives and selects at most one extension
per supplied record after the quality checks. Similar scores on substantially
different routes now produce a separate unresolved candidate. Similar variants
of the same saved route can still mark a selected result `extended_ambiguous`.
Neither outcome deletes the protospacer. The adapter, coverage, cycle,
incompatible-flank and work-budget checks remain active. Source-repeat crossings are flagged by default, as described below.
Optional `--branch-guard true` and `--resolve-graph-ties true` restore the strict
stopping/clipping policies for explicit comparisons. They still cannot remove a
supplied protospacer. A selected path remains a candidate, not proof of its origin.

| Output | Meaning in supplied-protospacer mode |
|---|---|
| `Protospacers.fasta`, `protospacer_candidates.tsv` | All supplied IDs/sequences with `retained_input` provenance |
| `ContiguousPaths.fasta`, `protospacer_regions.tsv` | Only selected actual extensions; no seed-only placeholders |
| `protospacer_extensions.tsv` | One row per input, independent extension status/reason, chosen region, lengths, stops, scores, ambiguity and work |
| `protospacer_associations.tsv` | Target-to-region coordinates and strand; observed allele and query alignment for variant placements |
| `protospacer_target_placements.tsv` | One row per input, observed target, alignment, placement work and limits |
| `protospacer_placement_checks.tsv` | Diagnostic probe/extension alternatives and scores; these are not additional final predictions |
| `protospacer_path_quality.tsv` | Candidate quality measurements, acceptance, selection and rejection reasons |
| `UnresolvedExtensions.fasta`, `protospacer_unresolved.tsv` | Best provisional path when substantially different routes remain near-tied; excluded from primary predictions |
| `protospacer_sources.tsv` | Matching source-array observations and zero-based, half-open coordinates in their original spacers |
| `protospacer_branch_checks.tsv`, `protospacer_adapter_stops.tsv` | Traversal diagnostics; branch probes are present only if enabled |
| `protospacer_parameters.json` | Mode, all effective algorithm settings, input/retained/extended counts |
| `retained_protospacer_mode.txt` | Marker allowing safe reruns in this output directory |

`extension_status` is `extended`, `extended_ambiguous`, `extended_variant`,
`extended_variant_ambiguous`, or `not_extended`.
`retained_input` records provenance, not a claim that this run validated a target.
Duplicate sequences with distinct supplied IDs keep their records; identical
extension regions can share a FASTA record. No target is silently substituted by
a shorter interval. The repeat and adapter checks do not establish viral origin.

The saved 81-match baseline is in `bench/fixtures/protospacer-baseline-81/`.
Its original reference matches are preserved by keeping exactly those sequences,
independently of extension success. Five intervals are
shorter than k=23. The corresponding manifest records thresholds and the FASTA
hash. No reference sequence is used to choose extension paths.

Verify a finished run independently (nonzero exit on lost/changed sequences,
wrong coordinates, duplicate associations, seed-only extensions or cap violations):

```bash
python3 bench/check_retained_protospacers.py \
  /path/to/accepted/Protospacers.fasta /path/to/separate-extension-results
```

This checker verifies retention and output consistency. It does not assess
reference-match quality of the extended flanks.

### One passage through source-array context

The default `--source-array-crossing true` keeps a branch after its source
repeat is recognized. The path is marked `crosses_source_array`, rather than
rejected at that boundary. A source-array vertex may not be revisited outside
the immutable seed, including across the joined left/right flanks. A second
approximate source-repeat copy is also rejected at the existing
`--array-repeat-identity` threshold: slightly different copies must not evade
the one-passage rule simply by using different k-mer vertices. External repeats
still use `--max-node-visits`. If both flanks cross a source repeat, or their
join would revisit a source vertex, they remain separate one-sided candidates.
The one-passage rule applies to the complete emitted sequence, not separately to
each half.

Source-cycle edges do not contribute to the coverage baseline or beam stability.
The baseline restarts after leaving that region; final rescoring uses the same
mask. This prevents repeat abundance from determining the downstream depth. The
equal-depth beam still keeps `--beam-width` competing paths (default 32), with
750 added bases per side and explicit work limits. A separate
`--source-context-width 256` pool keeps paths that have not yet accumulated k+1
unmasked downstream edges; they cannot be ranked by downstream stability while
still inside the array. That pool is also bounded and reports truncation. Among
scored paths, one best representative of each source-array exit is considered
before filling remaining beam slots with alternatives from the same exit. Once
enough external sequence is established, a return into source-cycle edges stops
that branch. A beam is bounded, so preserving a boundary crossing does not
guarantee every downstream route survives pruning or wins the final score.
References are never used for ranking.

Supplied-target reports include `left_crosses_source_array`,
`right_crosses_source_array`, and the number of scored edges on each side.
Discovery-mode associations include the same crossing flags. A flag records the
path's traversal, not a validation of viral origin. All supplied protospacers
remain unchanged regardless of which extension wins. Set
`--source-array-crossing false` for the previous boundary-rejection policy.
Parameters record `equal_depth_beam_source_crossing_v5` and the chosen policy.

For an explicit audit, add `--report-alternatives true` in supplied-target mode.
`ExtensionAlternatives.fasta` and `protospacer_alternatives.tsv` contain at most
twice `--beam-width` candidates per selected target placement (a left/right
conflict can split one joined candidate into two one-sided candidates), with
ranks, scores, coordinates, selection status and crossing flags. These are
diagnostic alternatives, not additional final predictions.
`ContiguousPaths.fasta` still contains at most one selected extension per
supplied record. The option is off by default and its sequence storage counts
against `--max-retained-bases`. A matching alternative does not imply that the
graph-only score selected it; evaluate and report those two outcomes separately.

### Target-placement correction and full input audit (2026-11-12)

The following paragraph describes the earlier one-substitution-anchor run.
The subsequent complete-placement audit and masked-rescue comparison are in
[the updated report](../bench/reports/protospacer-complete-placement-2026-11-12/README.md).

The graph-only retained-input workflow now selects 10 reference-matching
extensions, compared with 8 before the placement fallback. N=16, 32 and 64
all preserve the previous eight matches and every one of the 81 FASTA records.
The additional matches are 1,537 bp and 133 bp. A third demonstrated matching
path belongs to a short, ambiguous query and still loses the final ranking.
The code does not claim to select that eleventh known matching continuation.

The [per-input report and reproduction command](../bench/reports/protospacer-target-placement-2026-11-12/README.md)
cover all 81 inputs and 1,041 reference locations, including absent reference
alleles, completed negative searches, path witnesses and unresolved alternatives.
References are not used to generate or rank production predictions.

### Equal-depth comparison on SRR14017553 (2026-11-11)

N=16,32,64 each increased local reference-matching extensions from six to eight
at >=90% identity / >=50% single-HSP coverage. All 81 inputs are unchanged and
all six baseline matches survive. N=32 extends 45 input records into 40 unique
regions; it is the current default. Doubling the N=64 work budget removes its
reported work-limit stops but adds no reference match.
See the [measured comparison and limitations](../bench/reports/protospacer-beam-2026-11-11/README.md).

### SRR14017553 graph investigation (2026-11-10)

The downloaded k=23 graph and its matching 49-array/530-spacer report were used
to rerun the fixed 81-protospacer set. Two traversal defects were reproduced:

1. A rejected successor produced a terminal copy of the shared path even when
   another successor passed. Those artificial short candidates could outrank the
   continuing path. Now a rejection terminates the shared prefix only when no
   successor passes; rejected alternatives still cannot be traversed.
2. The fixed initial depth window could stop a continuing path after gradual
   depth drift. The walking gate now uses the median and MAD of the most recent
   accepted observations. It tests an edge **before** updating that window.
   The initial window remains separate for comparing the two flanks.

Both changes apply to ordinary discovery and `--protospacers` extension.
Parameters schema 7 identifies `supported_successors_local_median_v2`.
Source-repeat checks still examine complete repeats: a shared prefix alone can
also occur on a valid target and is not sufficient to discard it.

On the complete real graph, all 81 supplied sequences and IDs remained unchanged.
There were 45 extended input records, represented by 43 unique regions, and 36
without an accepted extension (31 context/adapter/extension failures and five
inputs shorter than k). Local BLAST against the saved Shah reference at >=90%
identity and >=50% single-HSP query coverage retained the same six matching
regions. It did **not** increase that match count. Four improved:

| Input seed suffix | Extension length, before → after | Longest qualifying HSP, before → after |
|---|---:|---:|
| `541b71503abcacf8` | 1315 → 1530 | 1315 → 1530 |
| `0e3fc52adc6e1e7e` | 987 → 1528 | 987 → 1528 |
| `815985dbbb929d9a` | 860 → 1240 | 466 → 799 |
| `c35af5b66b335ecd` | 860 → 1227 | 466 → 786 |

The first two match the reference across the full new sequence. The remaining
extensions are still candidates, and a reference match to the short spacer does
not establish the flanking sequence. Reference tracing also found cases where
the compared reference's next edge was absent from this graph. This does not
exclude another strain or another reference; it limits what that exact route can
reconstruct. No reference sequence participates in the production path scoring.

`tests/fixtures/protospacer-srr14017553/` contains observed graph fragments that
reproduce the premature stopping without distributing the complete graph.
`bench/trace_protospacer_graph.cpp` is the read-only diagnostic used to inspect
k-mer availability, neighboring edges and their multiplicities.

## Build and run on an existing graph

With the repository's normal build dependencies installed:

```bash
git fetch origin
git switch --track origin/codex/protospacer-extension
git submodule update --init --recursive
cmake -S . -B build/protospacer -DCMAKE_BUILD_TYPE=Release
cmake --build build/protospacer -j 24

build/protospacer/mcaat phage-curate \
  --arrays /path/to/array-results \
  --graph /path/to/previous-run/graph/graph \
  --output /path/to/protospacer-results \
  --threads 24 \
  --max-flank-length 750 \
  --max-length 1550
```

If the branch already exists locally, use `git switch codex/protospacer-extension`
and `git pull --ff-only` instead of `git switch --track`.
The graph argument is the SDBG file prefix, without `.sdbg_info`, or its containing
`graph` directory. This command loads the existing graph; it does not build a graph
or rerun array/cycle detection. Use a separate output directory for comparisons.

The array directory can contain either:

- `CRISPR_graph_paths.tsv`, with `array_id`, `observation_id`, `repeat`, `spacer`
  columns. These are the actual sequence observations, not a new consensus.
- Original `CRISPR_Arrays_*.txt` reports. An all-dash repeat field means the
  consensus repeat above it, as in MCAAT's writer. Literal variants are retained.

The TSV takes precedence when both formats exist. Sequences must contain ACGT;
mixed alignment gaps and ambiguous bases fail explicitly. Every source
observation is retained, including observations whose spacers cannot be mapped.
Whole search spacers (after boundary correction) must map exactly to the graph before their matching intervals are
considered. This version does not search for mismatched or indel-containing
protospacers.

## Detection and extension

1. Correct conserved terminal repeat sequence misassigned to spacers using the
   source-array observations (details below). Group identical corrected spacers
   and their reverse complements for processing, retaining
   all source array/observation links and their original strands. Distinct spacers
   sharing the same first k-mer remain separate queries.
2. Trace the complete spacer. Enumerate exact intervals with external graph
   connections. Internal fragments need external connections at both ends;
   one-sided intervals must reach the corresponding corrected spacer end. Longer
   intervals are considered first. The default minimum match is the graph's `k`.
3. Test the local source-array context on each growing flank, on both strands.
   Repeat copies are aligned with substitutions and indels allowed by
   `--array-repeat-identity` (default 0.85). For a partial seed, the unmatched
   spacer length places the expected repeat boundary; differences in those
   unaligned spacer bases do not by themselves establish an external target.
   A repeat-like flank is stopped and removed. Its opposite flank can still be
   retained if it establishes an external context. A shared repeat prefix is
   traversable until the path diverges or qualifies as a repeat copy. This is a
   per-path decision, not a global deletion of repeat edges or an entire spacer.
4. Extend at most 750 bp per side by default. Advance all live paths one base
   before pruning to `--beam-width` alternatives: paths of different progress
   never compete in the live beam. Terminal candidates are kept separately.
   The provisional beam score uses all traversed counts, including unresolved
   boundary context, so an untested one-base branch is not compared with a
   completed flank. Final ranking excludes source-context observations.
   Strict local lookahead is available with `--branch-guard true`.
   Allow up to `--max-node-visits 3` visits per vertex outside source-array cycles.
   For a vertex already repeated in the immutable input, permit at most two
   additional visits. A closed segment consisting entirely of edges from the
   source arrays' repeat/spacer cycles is rejected, even after the path previously
   left array context. All spacer observations from the same source arrays are
   included in this guard. A shared vertex alone is never blacklisted. Recheck
   both visit limits and source-cycle closures across the joined left/seed/right
   path. Existing repeats wholly within the supplied seed remain untouched.
   Check sequencing adapters during the walk. Roll back the complete detected
   adapter and stop only that branch; retain other continuations and the clean
   opposite flank. Recompute coverage on the retained sequence before ranking.
5. Retain a flank only when it establishes an external context and contains at
   least `k+1` bases outside the matched seed (one complete external graph edge).
   Shorter tails are removed; a supported opposite side remains eligible.
   Compare left/right multiplicity baselines and split incompatible pairs into
   separate one-sided candidates. Failed joins, work-limited seeds with no usable
   flank, and other unextended intervals go into separate diagnostic outputs.
   They are excluded from the primary detection FASTAs and associations. A
   completed search showing only array context is `no_external_flank`; a query
   with no accepted result and unfinished searches is `unresolved`.
6. Score multiplicity stability on the emitted flanks, including variation
   within each flank, adjacent changes and agreement between the two flanks.
   Optionally incorporate reads/read pairs and supplied PAM models. Retain up
   to `max-paths` candidates per seed **internally** for final comparison.
7. Consolidate nested intervals only when they produced the identical retained
   extension sequence at the same query alignment. Different paths/alignments
   remain candidates. Then compare candidates across **all seed intervals of
   the spacer** and select at most one prediction. Save its detected protospacer
   independently of extension selection. If the optional strict ambiguity rule
   withholds the extension, the already detected interval remains in the
   protospacer FASTA and receives a separate `not_extended` status.
   Write only selected predictions
   to the primary FASTA/TSV outputs. Identical final sequences, including reverse
   complements, share a region without losing their source associations.

Connectivity is evidence for a candidate path, not proof that its complete
sequence occurs in a single genome. Shared sequences can connect unrelated
flanks. Selecting the best available prediction does not establish its true
genomic origin. Comparison counts, score ties and search limits remain visible.

## Spacer boundary correction

By default, each source array must have at least five **distinct** spacer sequences.
For each end, find the longest exact terminal motif present in at least 85% of
those sequences; motifs shorter than six bases are ignored. Trim that motif only
from observations that contain it exactly. Duplicate observations do not inflate
support. For the same full original spacer, including its reverse complement,
all source arrays must agree on the proposed canonical boundaries; conflicting
proposals leave it unchanged. A corrected complete spacer must remain at least
graph `k` bases long. This does not impose a new length rule on the array detector.

If the reported spacer is `A X B`, move the conserved prefix `A` and suffix `B`
into the repeat: `repeat + A X B` becomes `(B + repeat + A) + X`. This preserves
the cyclic sequence and supplies consistent repeat boundaries for the subsequent
source-array context test. It is a conservation heuristic, not independent proof
of the biological boundary. Conserved genuine spacer ends can be mistaken for
repeat sequence; sparse arrays and divergent terminal repeats can remain uncorrected.
The choice of boundaries does not use reference matches.

Use `--boundary-support`, `--boundary-min-spacers`, and `--boundary-min-bases`
to adjust the evidence requirements, or `--correct-spacer-boundaries false`
to disable correction. The original files are never rewritten.
`protospacer_queries.tsv` preserves `spacer` and `original_repeat`, and adds
`search_spacer`, `search_repeat`, `trimmed_left_bp`, `trimmed_right_bp` and
`boundary_status`. Queries are deduplicated **after** correction, so several
different original observations can now share one search query. Association
coordinates continue to refer to the original source spacer.

## Short extensions and junction decisions

`--max-flank-length 750` caps each side independently. A complete 35-bp match can
therefore produce at most `750 + 35 + 750 = 1535` bp. A blocked left flank does
not allow 1500 bp on the right. `--max-length 1550` is an additional total cap,
including the seed; larger seeds can reduce the available flank budget. These
are maximum lengths: dead ends, adapters, ambiguity and other checks can stop
earlier. For up to 900 bp per side, use `--max-flank-length 900 --max-length 1850`.

Every outgoing fork and incoming merge supplies a coverage observation, even
between the ordinary spaced samples (`--junction-check true`). At an outgoing
fork, the optional `--branch-guard true` examines each continuation for at most 75 edges,
stopping earlier at the next fork, a dead end, a repeat context or an adapter.
These probes consume the same edge-expansion budget as traversal; they do not
launch recursive searches. Source-array and adapter branches are rejected.
When at least three prior baseline observations exist, the probe also rejects
an edge exceeding the usual Poisson/NB coverage-loss threshold against that
fixed baseline.

Eligible external counts in a probe are compared using:

```text
branch_score = mean coverage loss + variance(log2(1 + count))
             + mean squared adjacent change of log2(1 + count)
```

Lower is better. Without an established baseline, coverage loss is zero and
only internal stability is compared. Each competing branch needs at least
three eligible counts; an unresolved short/context-pending alternative prevents
claiming a winner. Choose the best only if the runner-up is worse by more than
`--branch-score-margin 1`. If only one branch survives the context/adapter/coverage
checks, follow it. Otherwise stop at the shared junction and retain only the
already supported prefix. Equal stable branches at different depths cannot be
resolved from internal stability alone. The mean score avoids an automatic
reward for longer probes, but correlated k-mers do not become independent
evidence: neither the score nor its margin is a confidence probability.

`protospacer_branch_checks.tsv` records each considered edge, probe length,
observation count, score, status, choice and decision. A score is usable only
with `status=scored`; zero on rejected/incomplete probes is not a good score.
Rows use the canonical search-spacer orientation. At most 4096 rows per query
are stored, with uncapped counts and stored-row counts in the parameters JSON.
`--branch-guard false` is the default and uses beam traversal at junctions. Supplying optional
read files or PAM rules also disables this early branch choice, preserving
alternatives for their existing reranking; `branch_guard_active` records this.

## Equal-depth beam and bounded repeats (schema 9)

`equal_depth_beam_local_window_v4` applies to both discovery and supplied-target
extension. The graph stays read-only. Length limits, work limits, source-repeat
edit checks, adapter checks and flank compatibility checks remain enforced.
A repeated external vertex can represent a real repeat; its copy number is still
ambiguous without additional linkage evidence. The visit cap bounds that
ambiguity rather than declaring a biologically correct copy number.

The beam score is log-count variance plus squared adjacent changes, the
multiplicity-one fraction, and `1/(observations+1)`. Final selection uses the
emitted external flanks and additionally compares their levels, as below.
Overlapping k-mers are correlated: these are ranking scores, not probabilities.
Wider beams need more work and can hit the same `--max-states` cap earlier in
base-pair progress. `search_limited` discloses pruning and budget exhaustion.

`bench/benchmark_retained_beams.py` runs N=16,32,64 sequentially against an existing
graph, checks every supplied ID/sequence and seed coordinate, then performs local
BLAST on predictions and the saved baseline. Reference sequences are never
provided to the graph search. Its comparison separates matching input records
from unique extension regions and reports baseline matches gained/lost, lengths,
work, revisited nodes and wall time including graph loading.

The [all-81 reference audit](../bench/reports/protospacer-reference-audit-2026-11-11/README.md)
records exact missing edges, verified connected witnesses, production path losses,
before/after matches and unresolved bounded searches for the supplied SRR14017553 graph.

## Selecting one final prediction

Use every emitted edge multiplicity outside the seed and known source-array context.
The local stopping gate retains its spaced samples; final profile scoring does not.
Its local median/MAD window advances on every accepted external edge, including
unitig interiors. A second moving window retains the separated observations.
At a sampled edge, accept support from either scale; stop only when the count is
incompatible with both. This handles gradual drift without letting a short spike
or valley erase the broader depth context. A sampled edge is checked before its
count can update either window. The initial baseline remains separately stored
for the left/right compatibility check. These are robust stopping heuristics,
not independent probability tests.
For each flank, transform each count `c` to `x = log2(1 + c)` and record its mean,
population variance and mean squared difference between consecutive observations.
Initial baseline observations participate in this calculation. All statistics are
recomputed on the emitted sequence, after any trimming.

```text
within  = (variance_left + variance_right) / 2
jumps   = (mean_squared_step_left + mean_squared_step_right) / 2
between = (mean_left - mean_right)^2          # when both sides have observations
support = 1/(n_left + 1) + 1/(n_right + 1)
weak = (fraction_of_count_1_edges_left + fraction_of_count_1_edges_right) / 2
stability_score = -(within + jumps + between + support + weak)
match_penalty = 1 - matched_bases / corrected_search_spacer_length
total_score = stability_score - match_penalty + optional_read_pair_PAM_score
```

Higher is better. The weak-edge fraction adds a soft penalty for multiplicity-one
runs; isolated count-one edges remain allowed. Above that floor, a uniformly
high-count path has no automatic advantage over a uniformly lower-count path. The between-flank term penalizes joining a stable low-depth flank
to a stable high-depth flank. The support penalty prevents tiny or unobserved
flanks from receiving a perfect score; one-sided reconstructions remain eligible.
Missing-side variance/step terms are zero, and missing observations contribute
through `support`. These are explicit ranking heuristics, not posterior
probabilities, p-values or calibrated confidence intervals.

With supplied read pairs and positive `pair-weight`, pair-supported candidates
still rank ahead of unsupported candidates. Within that tier, compare total
score; exact ties prefer more matched spacer bases, more observed external edges,
then a longer region and finally sequence order for reproducibility. Choose one
winner for the complete spacer query, including all its partial seed intervals.
`--max-paths 1` alone would not do this, because it limits candidates per seed.

In supplied-protospacer mode, the default additionally uses the existing
`prediction-margin` to choose among near ties. First sort by the ordinary score.
Among candidates within that fixed best-score margin (including the support
term), prefer the longest validated sequence, but only if its
`within + jumps + weak` loss is also within the same margin of the original
winner. Agreement between left/right mean depths cannot compensate for a large
deterioration within a flank. Equal lengths retain the original score ordering.
This prevents a slightly flatter adapter-trimmed prefix automatically beating
its supported continuation. It does not claim the longest path is correct.
Near ties remain explicitly ambiguous. `--resolve-graph-ties true` retains the
conservative common-prefix behavior; a zero margin disables length preference.
For a longer near-tied selection, the signed score gap can be negative because
the numerically highest-scoring shorter alternative remains the runner-up.

`protospacer_selection.tsv` gives the selected seed/region, score components,
sample counts, multiplicity levels, number of compared shortlisted candidates, and the
runner-up score for a different region sequence. `score_gap` is reported only
within the same pair-support tier and before ambiguity trimming. Numerical score
ties describe the original shortlist. After retaining only common sequence,
the old runner score is kept for diagnostics but `score_gap` is `NA`.
Low coverage evidence is
labelled `insufficient`, rather than presented as demonstrated stability.

## Graph-only acceptance checks

External-sequence and incompatible-flank checks are active in ordinary runs.
The strict ambiguity trimming below is opt-in. No read files or reference
database are required.

- **External sequence:** default `--min-external-flank 0` means `k+1` bases
  wholly outside the matched seed. At `k=23`, the minimum is 24 bp. This is a
  structural minimum, not a confidence threshold. If neither flank qualifies,
  keep the interval in unresolved diagnostics instead of the detection FASTA.
  A positive value overrides the minimum; `1` restores previous short-tail
  acceptance. Removing a short reconstruction does not establish that the
  underlying spacer match is false.
- **Incompatible flanks:** with at least three baseline observations per side,
  compare each side's median multiplicity against the other using the existing
  Poisson/NB loss. Use the greater observed overdispersion of the two sides and
  the existing `--max-coverage-penalty` threshold (20 by default). A conflict
  splits the joined candidate into two one-sided alternatives; final selection
  still produces at most one output. Too few observations cannot establish a
  conflict. `--joint-flank-check false` disables this gate.
- **Optional near-best trimming:** without supplied read or PAM evidence,
  `--resolve-graph-ties true` compares the winner with alternatives having the
  stability-plus-match evidence within `--prediction-score-margin 0.1` below
  the winner (or better), and the same number of retained sides.
  The sample-count reward is excluded from this ambiguity comparison: one extra
  observation cannot identify an otherwise indistinguishable branch. Align
  candidates by their corrected search-spacer coordinates and retain only identical
  contiguous sequence outward from the seed. Recheck context, minimum length,
  path coordinates and coverage on that sequence. If no external side remains,
  report the extension as unresolved while keeping the detected protospacer.
  This checks the bounded shortlist, not every path in the
  graph. The margin is a heuristic, not calibrated confidence; a value of zero
  restores exact-evidence ties only. `--resolve-graph-ties false` restores deterministic winner-only
  behavior; existing optional read/PAM reranking also keeps that behavior.

Selection flags `short_flank_trimmed`, `flank_multiplicity_conflict` and
`ambiguity_trimmed` disclose changes to the selected candidate.
`joint_flank_loss_before_split` records the pre-split comparison. Decisions
also include candidates withheld completely, so the selection table alone is
not a count of all rejected alternatives. Graph ties that are shortened carry
association ambiguity `shared_sequence_of_near_best_predictions`.

## Coverage model

Each flank has its own local baseline: the median of the most recent
`coverage-window` accepted eligible counts, ordinarily sampled at least `k+2` edges apart,
with additional observations at outgoing forks and incoming merges. Counts on
known array boundaries, the matched spacer, and the prefix before source-context
divergence are excluded. The robust variance
estimate is `(1.4826 * MAD)^2`; the overdispersion estimate is
`alpha = max(0, (variance - median) / median^2)`.
The first window is also retained unchanged for the separate two-flank comparison.
An observation is checked against the existing local window before being admitted;
a rejected jump cannot shift that window. Window storage is bounded.

- `auto`: negative binomial when estimated dispersion is positive, otherwise
  Poisson.
- `poisson`: Poisson throughout.
- `negative-binomial`: negative binomial, with an alpha floor of 0.05 when the
  local estimate is smaller.

For count `x`, baseline `mu`, and NB size `r = 1/alpha`, the penalty is the
log-likelihood loss relative to a saturated mean:

```text
Poisson: x log(x/mu) + mu - x
NB:      x log(x/mu) - (x+r) log((x+r)/(mu+r))
```

The zero-count logarithm contributes zero. After collecting the baseline, an
individual observation above `max-coverage-penalty` rejects that successor before
the incompatible edge. Other accepted successors continue, without creating a
terminal copy of their shared prefix. The diagnostic `coverage_score` is minus the mean loss per side, summed
over the two emitted flanks. Trimmed parts do not contribute. A side with too few
observations to finish its baseline has score zero and no scored samples; zero
is **not** a confidence statement. These are ranking heuristics, not calibrated
probabilities. Nearby k-mer counts remain correlated, and systematic coverage
changes can still stop a true extension. Final selection uses the stability score
above; the Poisson/NB model supplies the branch-stopping test.

## Optional read and read-pair evidence

Append to the command:

```bash
  --reads1 /path/to/sim_R1.fastq \
  --reads2 /path/to/sim_R2.fastq \
  --min-insert 150 --max-insert 600
```

Choose the insert bounds for the actual library; the defaults are 0–1000 bp.
FASTA/FASTQ and gzip are supported. A single read file is also supported. Reads
are streamed once; the candidate index is bounded by `max-read-index`. Mate names
must agree after removing `/1` and `/2` suffixes.

- Junction evidence counts exact `k+2`-mers spanning graph junctions, on either
  strand, once per word per read. `read_support` sums junction/read observations;
  it is not a count of distinct reads across all junctions.
- Pair evidence requires both full reads to match a retained region exactly,
  facing inward and within the specified outer insert span. Each pair contributes
  at most once per region. `unique_pair_support` means only one **retained region**
  supports that pair, not uniqueness in the complete metagenome.
- With positive `pair-weight`, pair-supported candidates rank ahead of
  unsupported candidates across all seeds of the same spacer. Within that tier, the score is
  stability score - match penalty + `read-weight * fraction_supported_junctions` +
  `pair-weight * log(1 + pair_support)` + any PAM bonus. Supporting several
  overlapping regions does not erase a pair's evidence.

These reads rerank graph-generated alternatives; they cannot recover an already
pruned path. Full-read matching does not tolerate sequencing errors or perform
alignment. A pair supporting part of a region does not validate its unspanned
ends. All candidates compete for the single final prediction per spacer.

## Optional PAM rules

Pass `--pam-rules rules.tsv`, with tab-separated columns:

```text
array_id  subtype         strand  side    motif
Array_1   supplied-model  +       3prime  NGG
Array_2   supplied-model  ?       5prime  AAG
Array_3   supplied-model  ?       none    -
```

The examples illustrate the format, not a built-in subtype catalogue. Supply
rules appropriate to the system and experimental convention:

- `strand=+` evaluates a sequence oriented like the original reported spacer;
  `-` evaluates its reverse complement; `?` means orientation is unknown.
- `5prime` and `3prime` specify the adjacent side in that oriented sequence.
- IUPAC motifs are accepted. `side=none`, `motif=-` explicitly records a model
  for which PAM testing is not applicable.

Full-spacer matches can be `match` or `mismatch`; incomplete targets, missing
flanks, unknown orientations and absent models have separate `unknown_*`
statuses. PAM mismatch is not a hard exclusion. A match adds `pam-weight`.
Several models for one array are allowed; a matching model wins. Ranking is
shared by sources of the same query and uses a match from any source. The
reported association's `total_score` and PAM status use only that source's model,
so per-source total scores can differ from the score used for the shared selection.
No subtype, orientation or PAM model is inferred from an unlabelled graph.

## Source-array repeat similarity

`--array-repeat-identity 0.85` permits an edit distance of at most
`floor((1 - 0.85) * repeat_length)` to a reported repeat. Substitutions,
insertions and deletions each cost one; alignment is anchored at the expected
boundary and cannot freely clip the repeat. Both repeat orientations and all
repeat observations belonging to the spacer are considered. This is a sequence
similarity setting, distinct from the array detector's cycle-support fraction.

Lowering the value recognizes more divergent copies as array context; raising
it can let variant array paths through. A value of 1 tests exact copies only.
This is a configurable exclusion heuristic, not proof of origin. Large spacer
length differences, repeat divergence beyond the tolerance, incomplete array
input, and unrelated graph recombinations remain limitations. The threshold is
applied to each path; another external path for the same spacer is preserved.

## Adapter-mediated extension failure and guard

An audit of 217 exported regions found the Illumina adapter core in 191 regions.
All 191 also contain confirming adapter context: an exact 13-base core followed
by six bases with at most two differences from a TruSeq adapter. In one
2,500-base region, a 387-base exact viral-reference match ends immediately before
the adapter. Multiplicity stability did not prevent this walk: a contaminating
adapter can itself have stable coverage and create connections between inserts.
This establishes a concrete failure mode, not an explanation for every unmatched
region. Reference absence and unrelated graph ambiguities remain possible.

The default guard recognizes both strands of these patterns:

| Sequence | Required evidence |
|---|---|
| `AGATCGGAAGAGCACACGT` | Exact first 13 bases; next six model bases within two edits |
| `AGATCGGAAGAGCGTCGTG` | Exact first 13 bases; next six model bases within two edits |
| `CTGTCTCTTATACACATCT` | Exact 19-base Nextera sequence |

The confirming bases matter: the 13-base Illumina core alone also occurs in
biological genomes. See the [Cutadapt TruSeq documentation](https://cutadapt.readthedocs.io/en/stable/guide.html#illumina-truseq)
for adapter sequences and that limitation. The anchored edit allowance here
handles substitutions and indels in the confirming context: 4–8 observed bases
can align to the six model bases with at most two edits. The exact core cannot
be edited. This recognizes deletion/insertion variants observed in the newer
exports that escaped the previous fixed-length substitution-only matcher.
It is a heuristic, not a posterior probability of contamination. An isolated
core without confirming sequence, changes in the core, and variants beyond the
allowance can escape it. Reads should still be adapter-trimmed before building
new graphs; this guard also works on an existing graph without changing its files.

Checks include sequence crossing the seed/flank boundary. If an adapter overlaps
the matched seed, that branch combination is rejected rather than outputting
adapter bases as a protospacer. A complete input spacer containing a recognized
adapter is reported with `mapping_status=adapter_in_spacer`. Adapter-free branches
for other queries remain searchable. The algorithm still selects one final
prediction per matched spacer; it does not label a retained region as confirmed
viral DNA merely because it passed the guard.

`--adapter-sequences SEQ[,SEQ...]` replaces the defaults with exact user-supplied
motifs, checked on both strands. Motifs must be 12–256 ACGT bases; at most 64 are
allowed. `--adapter-sequences none` explicitly disables the guard for comparisons.
The effective patterns, exact anchor coordinates, `max_edits`, and `allow_indels` are
saved in `protospacer_parameters.json`.

`protospacer_adapter_stops.tsv` records each blocked branch's seed, left/right
side in canonical query orientation, observed adapter motif, clean flank length,
whether it overlaps the seed, last clean endpoint vertex, and vertex completing
detection. These are actual graph vertex IDs from the current run. At most 4,096
diagnostic rows are stored per query; the JSON total `adapter_blocked_branches`
is not capped, and `adapter_diagnostic_rows` discloses any truncation. Repeated
interval searches and speculative junction probes can encounter the same adapter
and are counted separately. Diagnostic counts are encounters, not unique
contamination events or emitted paths.

Production regression (including a pre-fix comparison when supplied):

```bash
python3 bench/validate_protospacer_adapters.py \
  build/protospacer/mcaat build/adapter-validation \
  --before-binary /path/to/pre-fix/mcaat
```

This builds a real SDBG from synthetic reads, including the observed failing
spacer/adapter boundary. It checks both extension directions, a genuine competing
branch, Nextera, a contaminated input spacer, coordinate/path consistency,
deterministic threaded output, and unchanged graph files. It is not a rerun of
the original metagenomic graph, which was not available in the local audit data.

## Bounds and outputs

| Option | Default | Meaning |
|---|---:|---|
| `--max-length` | 1550 | Total emitted bp, including the matched seed and both flanks |
| `--max-flank-length` | 750 | Independent maximum bp on each side; unused capacity cannot exceed this cap |
| `--correct-spacer-boundaries` | true | Reassign conserved spacer ends to the repeat before searching |
| `--boundary-support` | 0.85 | Required fraction of distinct source-array spacers; allowed range (0.5, 1] |
| `--boundary-min-spacers` | 5 | Minimum distinct spacers in an array for boundary correction |
| `--boundary-min-bases` | 6 | Minimum length of a conserved terminal motif |
| `--junction-check` | true | Observe multiplicity at forks/merges between ordinary samples |
| `--branch-guard` | false | Stop or choose at forks using local lookahead; inactive with optional read/PAM evidence |
| `--branch-lookahead` | 75 | Maximum inspected edges per branch probe; range 1–1000 |
| `--branch-score-margin` | 1 | Required improvement over the next scored branch; lower score is better |
| `--prediction-score-margin` | 0.1 | Final graph-only alternatives this close to the winner retain common sequence |
| `--beam-width` | 32 | Live alternatives at equal extension depth; also bounds terminal/joined shortlists |
| `--max-node-visits` | 3 | Visits per new vertex outside source-array cycles; range 1–64; existing seed repeats remain intact |
| `--max-states` | 50000 | Edge expansions per seed side, including lookahead |
| `--max-query-states` | 200000 | Edge expansions across all intervals of one canonical spacer |
| `--max-seeds` | 128 | Retained intervals per canonical spacer |
| `--min-match-bp` | graph k | Minimum exact matched interval; CLI value 0 means graph k |
| `--min-external-flank` | 0 (= graph k+1) | Minimum retained bases wholly outside the matched seed |
| `--joint-flank-check` | true | Split incompatible left/right multiplicity contexts |
| `--resolve-graph-ties` | false | Retain common sequence for indistinguishable graph-only alternatives |
| `--max-paths` | 16 | Internal shortlist per seed; final output remains one prediction per spacer |
| `--coverage-window` | 7 | Baseline count observations per flank; minimum 3 |
| `--array-repeat-identity` | 0.85 | Source-repeat edit similarity threshold, allowed range 0.5–1 |
| `--adapter-sequences` | Confirmed Illumina / exact Nextera patterns above | Replace with comma-separated exact motifs, or `none` to disable |
| `--max-coverage-penalty` | 20 | Per-observation branch stopping threshold |
| `--max-retained-bases` | 50000000 | Global intermediate region sequence budget; error if exceeded |
| `--max-read-index` | 5000000 | Indexed candidate positions; error if exceeded |
| `--max-read-work` | 1000000 | Exact-placement comparisons and mate combinations per read pair; error if exceeded |
| `--read-weight / --pair-weight / --pam-weight` | 1 / 2 / 1 | Evidence score weights; zero disables that contribution |
| `--threads` | 1 | Independent-query OpenMP workers |

The retained-base limit is not a total RAM limit: graph storage, vertex paths
(about eight bytes per vertex), containers, threads and read indexes also consume
memory. Join combinations are bounded to four times the beam width. Search and
output pruning set `search_limited`; they do not claim an exhaustive result.

| Output | Contents |
|---|---|
| `protospacer_queries.tsv` | Every source observation, canonical query ID, boundary provenance and separate detection/extension status |
| `protospacer_candidates.tsv`, `Protospacers.fasta` | One detected exact interval per matched spacer, retained even if subsequent extension selection is unresolved |
| `MatchedSpacers.fasta` | One complete corrected search spacer per matched canonical query; partial-only queries are counted separately in JSON |
| `protospacer_regions.tsv`, `ContiguousPaths.fasta` | Unique accepted extended/partial region sequences, vertices and edges |
| `protospacer_associations.tsv` | Source → seed → region links, coordinates, strand, scores, rank, ambiguity, stop reasons and evidence |
| `protospacer_unresolved.tsv`, `UnresolvedSeeds.fasta` | Unresolved intervals without a retained detection, including unperformed searches and failed joins |
| `protospacer_decisions.tsv` | Rejected source contexts, unresolved path outcomes, and redundant intervals with their retained replacement seed |
| `protospacer_selection.tsv` | One row per selected prediction: stability components, comparison counts, runner-up score and ties |
| `protospacer_adapter_stops.tsv` | Blocked adapter branches and their actual graph endpoint/detection vertices |
| `protospacer_extensions.tsv` | Independent extension outcome for every retained detected interval |
| `protospacer_branch_checks.tsv` | Junction alternatives, local scores, actual vertex/edge IDs, choices and stop decisions |
| `protospacer_parameters.json` | Effective settings, input paths and result counts |

Output schema 6 separates detection retention from extension outcome.
In discovery mode, `extended_queries` and `detected_without_extension` split
`matched_queries`; the protospacer count can exceed the extension count.
`mapping_status` is the legacy combined status; use the explicit `detection_status`
and `extension_status` fields for the independent outcomes. Adapter patterns
use `max_edits`/`allow_indels`.
It declares `selection=one_prediction_per_distinct_spacer` and distinguishes `matched_queries`, `full_length_matched_queries`,
`partial_only_matched_queries`, `unresolved_queries` and
`unresolved_seed_intervals`. The matched-query denominator is `unique_queries`,
with identical spacers and reverse complements processed together. Counting
`Protospacers.fasta` gives selected intervals; counting `MatchedSpacers.fasta` gives
distinct corrected search spacers with a selected prediction. Both counts are now equal;
the number of unique sequences in `ContiguousPaths.fasta` cannot exceed that count.
A full corrected search spacer
in `MatchedSpacers.fasta` is not a claim that its entire length matched: consult
the seed table's coverage. No `seed_only` rows occur in primary associations.

Coordinates are **zero-based, half-open**. Seed-table query coordinates refer to
the corrected canonical search-spacer orientation; association query coordinates
refer to the original source spacer (including its trim offset). Decision and
unresolved-table query coordinates also refer to the corrected canonical spacer. Association `seed_begin/end`, left/right evidence and
stop reasons refer to the stored region sequence. `strand` tells whether that
region's seed matches the source directly or as a reverse complement.
`query_coverage` uses corrected search-spacer length and is a fraction, not a percentage; `identity=1` describes an exact
interval, not full query coverage. Stable IDs are sequence-derived hashes.

Stop reasons distinguish dead ends, cycles, source-array context, `sequencing_adapter`,
`insufficient_external_sequence`, `flank_multiplicity_conflict`,
`ambiguous_graph_extension`, `ambiguous_branch`, `no_supported_branch`, coverage breaks,
length limits, work limits and incompatible flanks. `coverage_samples` counts
scored observations after baseline collection. `expanded_states` is the total
edge work for that seed search, repeated on its associations, not additive across
rows. A single selected path may still have pruned alternatives. The association
`ambiguity` reports `selected_best_of_alternatives` or `selected_score_tie` when
appropriate; `alternative_rank` remains 1 for the selected result. Selecting one
path does not imply that it was the only feasible path.
`left_context` and `right_context` distinguish `external`, `array`, `not_retained`, and
`incomplete` context tests in the stored region orientation. `pending` with zero
flank bases denotes a side where extension did not start. Repeat-compatible
parts trimmed from an accepted one-sided region do not contribute coverage.
The decision report can contain multiple outcomes for one seed; a successful
alternative and a rejected array alternative can coexist. The unresolved FASTA
only contains intervals with no viable detection; a lower-ranked valid candidate
does not become unresolved because another candidate won. Completely skipped intervals
beyond `max-seeds` are disclosed by `eligible_intervals` and `seed_limit_reached`.
The decision report's interval-consolidation references describe the internal
shortlist; a referenced seed need not be the final winner.

Successful reruns replace the discovery mode's fifteen named files, including empty FASTAs;
they do not append old results. Computation/input failures before publication
leave previous results intact. Files are staged then individually renamed; this
is not a filesystem-wide atomic transaction. Array and graph inputs are unchanged.

`viral_status=not_assessed` is deliberate. These outputs are candidate target
regions; identifying a viral origin, estimating completeness or claiming a full
phage requires separate reference/composition/viral-marker validation. No external
viral classifier is silently run or bundled.

## Validation

```bash
cmake -S tests -B build/protospacer-tests -DCMAKE_BUILD_TYPE=Release
cmake --build build/protospacer-tests -j 2
ctest --test-dir build/protospacer-tests --output-on-failure

cmake -S tests -B build/protospacer-sanitizers -DPROTO_SANITIZERS=ON
cmake --build build/protospacer-sanitizers -j 2
ctest --test-dir build/protospacer-sanitizers --output-on-failure

# The output directory must be empty. Builds a real SDBG from 150-bp reads.
python3 bench/validate_protospacer_extension.py \
  build/protospacer/mcaat build/protospacer-validation

# Variant-repeat regression: host-only negatives and a target for the same spacer.
python3 bench/validate_protospacer_context.py \
  build/protospacer/mcaat build/protospacer-context-validation

# Graph-only acceptance: short tails, incompatible depths and tied branches.
python3 bench/validate_protospacer_graph_filters.py \
  build/protospacer/mcaat build/protospacer-graph-filter-validation

# Explicit strict policy: conserved ends, junction ambiguity and 750-bp caps.
python3 bench/validate_protospacer_boundaries.py \
  build/protospacer/mcaat build/protospacer-boundary-validation

# Current defaults: frozen detections, longer branch comparison and failure retention.
python3 bench/validate_retained_protospacers.py \
  build/protospacer/mcaat build/retained-protospacer-validation
```

Core checks cover overlapping seeds, bidirectional/one-sided reconstruction,
coverage compatibility, branch and work limits, intersecting flanks, repeated
seed k-mers, strand-aware reads/PAMs, pair ambiguity and invalid parameters.
The production fixture checks both array input formats, actual SDBG traversal,
source/region coordinates, reverse-complement deduplication, deterministic
threads/reruns, explicit caps and failed/empty rerun behavior.

In the general controlled fixture, five source observations yield two selected
seed intervals, two regions and three source associations. With the explicitly enabled strict junction guard, both
selected regions are contained in known source sequences; the previous guard
version retained one true segment and one recombinant in this fixture.
The fixture deliberately shares a 23-mer between different spacers. Resolving
this example does not mean that stability can resolve every graph ambiguity. With the supplied
pairs, the selected full-spacer target for the two opposite-strand source
observations is a true target segment. This demonstrates selection and physical
evidence handling, **not biological sensitivity or precision**. Metagenome-scale
runtime and accuracy require server testing.

The separate context regression uses a repeat/spacer context observed in a
failing output. Host-only loci with substituted, indel-containing and exact
repeat copies now produce zero detections. Adding an independent target for the
same spacer retains one matched interval and one selected region, contained in
the true target sequence, with both source strands preserved.
Setting repeat identity to 1 reproduces variant-repeat leakage, exercising the
configured threshold. Unit tests also cover a genuine target sharing more than
`k` bases with a repeat prefix, short unresolved contexts, spacer variation in
partial matches, and interval consolidation. Selection tests compare matching
against crossed low/high-multiplicity flanks, penalize within-flank instability,
and verify that absent or short coverage profiles do not get a perfect score.

The graph-only fixture retains the unambiguous 1,000-bp target. It withholds
10-bp tails, removes the ambiguous part beyond a shared branch, and withholds
an arbitrary left/right pairing at a shared seed. The incompatible-depth case
is split into one-sided candidates and then withheld as ambiguous under the explicitly enabled strict policy. A separate `post_selection_only` run disables the new junction checks
and near-score margin to exercise the previous splitting/tie gates independently.
A zero-margin control restores the slightly better one-sided depth candidate
with a nonzero score gap, verifying that near-score handling is configurable.
With all gates disabled, the short tails and arbitrary tied paths return; the
depth case was already one-sided under the previous ranking, so it is not evidence
of a corrected previous final false call. These fixtures verify deterministic
one/two-thread outputs and unchanged graph files. They do not establish sensitivity
or precision in the original metagenome.

The explicit strict-policy boundary/junction fixture builds a real SDBG with eight target spacers,
source-array contexts and an indistinguishable alternate target for one spacer.
Sixteen source observations (both strands) have their conserved seven-base suffix
reassigned. The seven unambiguous targets are retained and the ambiguous query
is withheld. Output includes 1,535-bp reconstructions (750 + 35 + 750), original
source-coordinate checks, unchanged graph/array hashes, both input formats,
explicit boundary opt-out, and byte-identical non-parameter files with one and
24 threads. Unit checks additionally exercise off-grid junction coverage jumps,
near-equal versus clearly different lookahead scores, conflicting boundaries,
minimum spacer lengths, distinct-spacer support and no transfer beyond the
opposite 750-bp cap.

Production boundary normalization was also compared with the earlier independent
sequence audit: all 199 audited original queries have identical corrected bounds
when loaded from the complete 530-observation source report. This checks boundary
normalization only, not traversal of the unavailable original metagenomic graph.

The release build also excludes an unused experimental debug translation unit
with obsolete interfaces and guards Linux-only `malloc.h`/`malloc_trim` calls.
Those portability fixes do not change the array-search/filter algorithms.
The graph builder's memory conversion is also corrected: settings store GiB,
while the underlying MEGAHIT `--host_mem` argument takes bytes. Existing-graph
`phage-curate` runs do not invoke that builder.

After the selection change, local smoke checks used existing graphs CP058556.1 and
CP071793.1 with two threads, a 500-bp output limit, 2,000 expansions per side and
8,000 per query. Of 317 and 277 distinct input spacers, respectively, 22 and 51
had selected predictions. The selector compared 229 and 784 shortlisted paths;
the final FASTAs contain 18 and 46 distinct region sequences because different
spacers can select the same region. All emitted source coordinates, sequence/path lengths,
distinct-spacer counts and edge-work bounds passed checks; no unextended seeds
occurred in detection outputs. These counts do not establish target accuracy:
the array input itself includes candidate detections, and these genome fixtures
are much smaller than a full pooled metagenome.

After enabling the graph-only acceptance checks above, the same bounded smoke
runs retain 21 predictions / 17 regions for CP058556.1 and 45 predictions /
40 regions for CP071793.1. Six and 22 selected predictions, respectively, are
shortened to common sequence because of indistinguishable alternatives. One
selected CP071793.1 prediction has incompatible flanks split before selection.
All source/seed coordinates, graph path-length accounting and one-prediction
invariants pass. Lower counts alone are not evidence of improved accuracy.

Replaying the shared C++ context classifier on the saved failing run recognized
all 694 previously identified repeat-flanked regions as array context on both
sides. This checked saved sequences, not a rerun of the original graph search;
the original graph was not available locally. The current core suite passes
321 checks in both release and AddressSanitizer/UndefinedBehaviorSanitizer builds.

With boundary correction and the new junction/near-score guards, the same bounded
genome smoke runs retain 2 predictions / 2 regions for CP058556.1 and 12 predictions /
11 regions for CP071793.1. No boundary changes qualify in these source reports.
Source coordinates, path accounting, independent flank caps and work bounds pass.
The drop from 21/45 predictions is a material sensitivity tradeoff: unresolved
junctions now stop earlier, and valid extensions can be lost as well as incorrect
ones. These runs do not measure precision. Test the defaults on the original
server graph before interpreting a reduced count as an improvement.


## Retention regression and current policy

### Reference-free quality filtering of supplied protospacers

The supplied-target CLI now defaults to `--path-quality-filter true`. This is
applied to `--protospacers` runs, not array discovery. Every supplied ID and DNA
sequence remains in `Protospacers.fasta`; a rejected extension is reported as
`not_extended`, with no sequence fabricated to meet a length target. The former
selection behavior is available with `--path-quality-filter false`.

In quality mode, the whole supplied target must occur exactly in the graph.
Targets shorter than graph k are embedded in observed graph k-mers without
changing or clipping their bases. The approximate-allele and masked-k-mer rescue
are not used in this mode. The historical variant tests below describe the
quality-filter-disabled workflow.

Search still uses a bounded, equal-depth beam. Low-complexity sequence adds a
soft beam-ranking penalty. A single multiplicity jump no longer terminates the
branch immediately in quality mode; final checks assess sustained changes.
Source-array passage and node-visit limits remain enforced. Search/work limits
are reported and do not prove that other paths are absent.

After searching, inspect each retained complete path and its individual left
and right flanks. Reject a candidate if any existing flank has:

- More than 20% low-complexity windows: 64 bases per window, step 16, triplet
  entropy below 3 bits.
- More than 20% repeated k-mer occurrences within the flank. Bounded external
  repeats below that fraction remain permitted.
- Fewer than 64 downstream edge-count observations after excluding the seed
  overlap and known source-array context.
- More than 50% singleton edge counts.
- A sustained depth step between adjacent, nonoverlapping 64-edge blocks:
  median high is at least four times median low (low floored at one), and
  high minus low is at least `3 * sqrt(high + low)`.

Also require at least `--min-extension-flank` added bases on **one** side
(default 250; use 750 for a full-750 experiment), and the unchanged full input
target. Two shorter flanks are not added together to satisfy this requirement.

Rank passing candidates by the worse existing flank's loss: variance of
`log2(1 + block median)` about the median, plus singleton fraction,
low-complexity fraction, repeated-k-mer fraction, and `1/(observations + 1)`.
An absent flank contributes no loss. Exact score ties prefer longer sequence,
then lexical order. If another passing candidate is within the configured
prediction-score margin (default 0.1) and follows a substantially different
same-side route, withhold the primary prediction. Compare flanks of at least
128 bases using their distinct `(k+1)`-mers: less than 80% overlap relative to
the smaller set indicates a different route. Small bubbles and an opposite
flank alone do not trigger this rule. The best unresolved candidate is retained
in `UnresolvedExtensions.fasta` with reasons in `protospacer_unresolved.tsv`.
This is ambiguity within the saved beam, not proof of uniqueness in the graph.

Output at most one primary extension per input. Near-score variants of the same
route remain marked ambiguous; passing these heuristics is not proof
of viral origin or a calibrated correctness probability. No reference sequence,
BLAST result, or trained reference label is an input to these decisions.

`protospacer_path_quality.tsv` records each assessed candidate, its acceptance,
selection, measured support, sequence, and rejection reasons. Optional
`ExtensionAlternatives.fasta` remains a diagnostic file containing passing and
rejected paths, not extra final predictions. `ContiguousPaths.fasta` contains
only selected passing outputs. Parameters record the quality-policy version
and minimum requested flank. Use `protospacer_extensions.tsv` to account for
every original input, including those for which nothing passed.

The 719632c strict defaults caused a large drop in retained predictions. The
subsequent source-crossing release restored bounded beam traversal and one-best-path selection. Strict
local stopping and common-sequence clipping are available only when requested.
The historical smoke counts above describe the strict version, not the current
supplied-protospacer workflow.

The production retention fixture keeps all ten supplied records, including
unmapped, short, source-context-missing, reverse, duplicate-sequence and variant
target records. Straight and branching targets extend to their known 1,535-bp
sequences; the variant target extends to its known 1,537-bp sequence while the
original input remains unchanged. Exact-only mode rejects that unmapped original
query. Strict options reduce the branching target to 785 bp. Adapter continuation
is excluded while a 785-bp clean fragment remains. One- and 24-thread outputs
match; tiny-budget reruns remove stale extensions while retaining inputs, and
invalid reruns leave prior output intact. The 81-record baseline is also preserved
byte for byte on an unrelated synthetic graph. Separate full-graph measurements
are linked above; this synthetic fixture tests retention and implementation,
not accuracy on the original metagenome.
