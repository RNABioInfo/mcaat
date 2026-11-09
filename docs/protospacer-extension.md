# Protospacer detection and bounded target extension

This implementation is on `codex/protospacer-extension`, based on
`origin/master` at `4f75041`. It replaces the production `phage-curate`
implementation; `protospacer-curate` is an alias. It does not incorporate the
separate array-filter branch. The old experimental classes remain in the
repository for reference, but the release subcommand uses the new module.

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
  --threads 24
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
Whole spacers must map exactly to the graph before their matching intervals are
considered. This version does not search for mismatched or indel-containing
protospacers.

## Detection and extension

1. Group identical spacers and their reverse complements for processing, retaining
   all source array/observation links and their original strands. Distinct spacers
   sharing the same first k-mer remain separate queries.
2. Trace the complete spacer. Enumerate exact intervals with external graph
   connections. Internal fragments need external connections at both ends;
   one-sided intervals must reach the corresponding original spacer end. Longer
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
4. Extend both sides with a bounded, coverage-ranked beam. Walk nonbranching
   stretches in place; keep alternatives at branches. New flanks cannot revisit
   seed or flank vertices. Repeated k-mers already present in the input spacer
   are allowed. Joining left and right flanks rejects intersections between them.
5. Retain one-sided and short reconstructions only when at least one emitted
   flank establishes an external context. A short flank still compatible with a
   source repeat is unresolved. Failed joins, work-limited seeds with no usable
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
   the spacer** and select exactly one prediction. Write only selected predictions
   to the primary FASTA/TSV outputs. Identical final sequences, including reverse
   complements, share a region without losing their source associations.

Connectivity is evidence for a candidate path, not proof that its complete
sequence occurs in a single genome. Shared sequences can connect unrelated
flanks. Selecting the best available prediction does not establish its true
genomic origin. Comparison counts, score ties and search limits remain visible.

## Selecting one final prediction

Use sampled edge multiplicities outside the seed and known source-array context.
For each flank, transform each count `c` to `x = log2(1 + c)` and record its mean,
population variance and mean squared difference between consecutive observations.
Initial baseline observations participate in this calculation. All statistics are
recomputed on the emitted sequence, after any trimming.

```text
within  = (variance_left + variance_right) / 2
jumps   = (mean_squared_step_left + mean_squared_step_right) / 2
between = (mean_left - mean_right)^2          # when both sides have observations
support = 1/(n_left + 1) + 1/(n_right + 1)
stability_score = -(within + jumps + between + support)
match_penalty = 1 - matched_bases / input_spacer_length
total_score = stability_score - match_penalty + optional_read_pair_PAM_score
```

Higher is better. Stability rewards consistent multiplicity, not higher absolute
multiplicity. The between-flank term penalizes joining a stable low-depth flank
to a stable high-depth flank. The support penalty prevents tiny or unobserved
flanks from receiving a perfect score; one-sided reconstructions remain eligible.
Missing-side variance/step terms are zero, and missing observations contribute
through `support`. These are explicit ranking heuristics, not posterior
probabilities, p-values or calibrated confidence intervals.

With supplied read pairs and positive `pair-weight`, pair-supported candidates
still rank ahead of unsupported candidates. Within that tier, compare total
score; exact ties prefer more matched spacer bases, more sampled observations,
then a longer region and finally sequence order for reproducibility. Choose one
winner for the complete spacer query, including all its partial seed intervals.
`--max-paths 1` alone would not do this, because it limits candidates per seed.

`protospacer_selection.tsv` gives the selected seed/region, score components,
sample counts, multiplicity levels, number of compared shortlisted candidates, and the
runner-up score for a different region sequence. `score_gap` is reported only
within the same pair-support tier. Numerical score ties are disclosed even though
the deterministic tie-break chooses one sequence. Low coverage evidence is
labelled `insufficient`, rather than presented as demonstrated stability.

## Coverage model

Each flank has its own frozen local baseline: the median of the first
`coverage-window` eligible counts, sampled at least `k+2` edges apart. Counts on
known array boundaries, the matched spacer, and the prefix before source-context
divergence are excluded. The robust variance
estimate is `(1.4826 * MAD)^2`; the overdispersion estimate is
`alpha = max(0, (variance - median) / median^2)`.

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
individual observation above `max-coverage-penalty` stops that branch before the
incompatible edge. The diagnostic `coverage_score` is minus the mean loss per side, summed
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

## Bounds and outputs

| Option | Default | Meaning |
|---|---:|---|
| `--max-length` | 2500 | Total emitted bp, including the matched seed and both flanks |
| `--beam-width` | 32 | Retained search frontier and joined alternatives before read scoring |
| `--max-states` | 50000 | Edge expansions per seed side |
| `--max-query-states` | 200000 | Edge expansions across all intervals of one canonical spacer |
| `--max-seeds` | 128 | Retained intervals per canonical spacer |
| `--min-match-bp` | graph k | Minimum exact matched interval; CLI value 0 means graph k |
| `--max-paths` | 16 | Internal shortlist per seed; final output remains one prediction per spacer |
| `--coverage-window` | 7 | Baseline count observations per flank; minimum 3 |
| `--array-repeat-identity` | 0.85 | Source-repeat edit similarity threshold, allowed range 0.5–1 |
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
| `protospacer_queries.tsv` | Every source observation, its canonical query ID and mapping status |
| `protospacer_candidates.tsv`, `Protospacers.fasta` | One selected exact seed interval per matched spacer, with an external flank and graph vertices |
| `MatchedSpacers.fasta` | One complete input spacer per matched canonical query; partial-only queries are counted separately in JSON |
| `protospacer_regions.tsv`, `ContiguousPaths.fasta` | Unique accepted extended/partial region sequences, vertices and edges |
| `protospacer_associations.tsv` | Source → seed → region links, coordinates, strand, scores, rank, ambiguity, stop reasons and evidence |
| `protospacer_unresolved.tsv`, `UnresolvedSeeds.fasta` | Unresolved intervals without a retained detection, including unperformed searches and failed joins |
| `protospacer_decisions.tsv` | Rejected source contexts, unresolved path outcomes, and redundant intervals with their retained replacement seed |
| `protospacer_selection.tsv` | One row per selected prediction: stability components, comparison counts, runner-up score and ties |
| `protospacer_parameters.json` | Effective settings, input paths and result counts |

Output schema 3 declares `selection=one_prediction_per_distinct_spacer` and distinguishes `matched_queries`, `full_length_matched_queries`,
`partial_only_matched_queries`, `unresolved_queries` and
`unresolved_seed_intervals`. The matched-query denominator is `unique_queries`,
with identical spacers and reverse complements processed together. Counting
`Protospacers.fasta` gives selected intervals; counting `MatchedSpacers.fasta` gives
distinct input spacers with a selected prediction. Both counts are now equal;
the number of unique sequences in `ContiguousPaths.fasta` cannot exceed that count.
A full input spacer
in `MatchedSpacers.fasta` is not a claim that its entire length matched: consult
the seed table's coverage. No `seed_only` rows occur in primary associations.

Coordinates are **zero-based, half-open**. Seed-table query coordinates refer to
the canonical spacer orientation; association query coordinates refer to the
original source spacer. Association `seed_begin/end`, left/right evidence and
stop reasons refer to the stored region sequence. `strand` tells whether that
region's seed matches the source directly or as a reverse complement.
`query_coverage` is a fraction, not a percentage; `identity=1` describes an exact
interval, not full query coverage. Stable IDs are sequence-derived hashes.

Stop reasons distinguish dead ends, cycles, source-array context, coverage breaks,
length limits, work limits and incompatible flanks. `coverage_samples` counts
scored observations after baseline collection. `expanded_states` is the total
edge work for that seed search, repeated on its associations, not additive across
rows. A single selected path may still have pruned alternatives. The association
`ambiguity` reports `selected_best_of_alternatives` or `selected_score_tie` when
appropriate; `alternative_rank` remains 1 for the selected result. Selecting one
path does not imply that it was the only feasible path.
`left_context` and `right_context` distinguish `external`, `array`, and
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

Successful reruns replace the module's twelve named files, including empty FASTAs;
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
```

Core checks cover overlapping seeds, bidirectional/one-sided reconstruction,
coverage compatibility, branch and work limits, intersecting flanks, repeated
seed k-mers, strand-aware reads/PAMs, pair ambiguity and invalid parameters.
The production fixture checks both array input formats, actual SDBG traversal,
source/region coordinates, reverse-complement deduplication, deterministic
threads/reruns, explicit caps and failed/empty rerun behavior.

In the general controlled fixture, five source observations yield two selected
seed intervals, two regions and three source associations. Without reads, one
selected region is contained in a known source sequence and one is recombinant.
The fixture deliberately shares a 23-mer between different spacers: stable
multiplicity alone cannot resolve every graph ambiguity. With the supplied
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

Replaying the shared C++ context classifier on the saved failing run recognized
all 694 previously identified repeat-flanked regions as array context on both
sides. This checked saved sequences, not a rerun of the original graph search;
the original graph was not available locally. The core suite passes 158 checks
in both the release and AddressSanitizer/UndefinedBehaviorSanitizer builds.
