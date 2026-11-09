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
3. Mark the observed `repeat|spacer|repeat` boundary context on both strands.
   Pure-repeat edges are excluded from extension. Mixed boundary k-mers remain
   traversable: a real target flank can initially share bases with the repeat.
   Known array-context counts never set the target coverage baseline. A side
   containing only that context is trimmed from the reported target region.
4. Extend both sides with a bounded, coverage-ranked beam. Walk nonbranching
   stretches in place; keep alternatives at branches. New flanks cannot revisit
   seed or flank vertices. Repeated k-mers already present in the input spacer
   are allowed. Joining left and right flanks rejects intersections between them.
5. Retain one-sided and short reconstructions. If search limits or an incompatible
   join prevent reconstruction, retain an unresolved seed with an explicit reason.
   A completed search showing only the original array context is reported as
   `no_external_flank`, not as a new target.
6. Optionally score reads/read pairs and PAM models, then retain up to `max-paths`
   alternatives **per seed**. Deduplicate complete region sequences, including
   reverse complements, without discarding their source associations.

Connectivity is evidence for a candidate path, not proof that its complete
sequence occurs in a single genome. Shared sequences can connect unrelated
flanks. Alternative paths, unresolved limits and all source links remain visible.

## Coverage model

Each flank has its own frozen local baseline: the median of the first
`coverage-window` eligible counts, sampled at least `k+2` edges apart. Counts on
known array boundaries and the matched spacer are excluded. The robust variance
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
incompatible edge. The coverage score is minus the mean loss per side, summed
over the two emitted flanks. Trimmed parts do not contribute. A side with too few
observations to finish its baseline has score zero and no scored samples; zero
is **not** a confidence statement. These are ranking heuristics, not calibrated
probabilities. Nearby k-mer counts remain correlated, and systematic coverage
changes can still stop a true extension.

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
- With positive `pair-weight`, pair-supported alternatives rank ahead of
  unsupported alternatives within the same seed. Within that tier, the score is
  coverage score + `read-weight * fraction_supported_junctions` +
  `pair-weight * log(1 + pair_support)` + any PAM bonus. Supporting several
  overlapping regions does not erase a pair's evidence.

These reads rerank graph-generated alternatives; they cannot recover an already
pruned path. Full-read matching does not tolerate sequencing errors or perform
alignment. A pair supporting part of a region does not validate its unspanned
ends. Unsupported candidates remain alternatives unless removed by output limits.

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
so per-source total scores need not follow the shared `alternative_rank`.
No subtype, orientation or PAM model is inferred from an unlabelled graph.

## Bounds and outputs

| Option | Default | Meaning |
|---|---:|---|
| `--max-length` | 2500 | Total emitted bp, including the matched seed and both flanks |
| `--beam-width` | 32 | Retained search frontier and joined alternatives before read scoring |
| `--max-states` | 50000 | Edge expansions per seed side |
| `--max-query-states` | 200000 | Edge expansions across all intervals of one canonical spacer |
| `--max-seeds` | 128 | Retained intervals per canonical spacer |
| `--min-match-bp` | graph k | Minimum exact matched interval; CLI value 0 means graph k |
| `--max-paths` | 16 | Final alternatives per seed |
| `--coverage-window` | 7 | Baseline count observations per flank; minimum 3 |
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
| `protospacer_candidates.tsv`, `Protospacers.fasta` | Exact matched seed intervals and graph vertices |
| `protospacer_regions.tsv`, `ContiguousPaths.fasta` | Unique extended/partial region sequences, vertices and edges |
| `protospacer_associations.tsv` | Source → seed → region links, coordinates, strand, scores, rank, ambiguity, stop reasons and evidence |
| `protospacer_parameters.json` | Effective settings, input paths and result counts |

Coordinates are **zero-based, half-open**. Seed-table query coordinates refer to
the canonical spacer orientation; association query coordinates refer to the
original source spacer. Association `seed_begin/end`, left/right evidence and
stop reasons refer to the stored region sequence. `strand` tells whether that
region's seed matches the source directly or as a reverse complement.
`query_coverage` is a fraction, not a percentage; `identity=1` describes an exact
interval, not full query coverage. Stable IDs are sequence-derived hashes.

Stop reasons distinguish dead ends, cycles, pure-array context, coverage breaks,
length limits, work limits and incompatible flanks. `coverage_samples` counts
scored observations after baseline collection. `expanded_states` is the total
edge work for that seed search, repeated on its associations, not additive across
rows. A single retained path may still have pruned alternatives.

Successful reruns replace the module's seven named files, including empty FASTAs;
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
```

Core checks cover overlapping seeds, bidirectional/one-sided reconstruction,
coverage compatibility, branch and work limits, intersecting flanks, repeated
seed k-mers, strand-aware reads/PAMs, pair ambiguity and invalid parameters.
The production fixture checks both array input formats, actual SDBG traversal,
source/region coordinates, reverse-complement deduplication, deterministic
threads/reruns, explicit caps and failed/empty rerun behavior.

In the controlled fixture, five source observations yield four seed candidates,
16 regions and 24 source associations. Eight regions are entirely contained in
the known source sequences; eight graph-compatible alternatives are recombined.
With the supplied pairs, the first-ranked full-spacer target for the two opposite
strand source observations is a true target segment. This demonstrates retention,
ranking and ambiguity handling, **not biological sensitivity or precision**.
The fixture deliberately keeps alternatives rather than presenting them as
confirmed phages. Metagenome-scale runtime and accuracy require server testing.

The release build also excludes an unused experimental debug translation unit
with obsolete interfaces and guards Linux-only `malloc.h`/`malloc_trim` calls.
Those portability fixes do not change the array-search/filter algorithms.

Additional local smoke checks used existing graphs CP058556.1 and CP071793.1
with two threads, a 500-bp output limit, 2,000 expansions per side and 8,000 per
query. They processed 317 and 277 source spacers in approximately 1.7 and 4.7
seconds respectively, including graph loading. All emitted source coordinates,
sequence/path lengths and edge-work bounds passed checks. These runs do not
establish target accuracy: the array input itself includes candidate detections,
and these genome fixtures are much smaller than a full pooled metagenome.
