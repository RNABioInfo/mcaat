<p align="center">
  <img src="icon.png?v=2" alt="MCAAT" width="330"/>
</p>

<p align="center">
  <a href="https://www.gnu.org/licenses/gpl-3.0"><img src="https://img.shields.io/badge/License-GPLv3-blue.svg" alt="License: GPL v3"/></a>
  <a href="https://hub.docker.com/r/feeka94/mcaat"><img src="https://img.shields.io/docker/v/feeka94/mcaat?label=Docker%20Hub&color=2496ED" alt="Docker Hub"/></a>
  <a href="https://hub.docker.com/r/feeka94/mcaat"><img src="https://img.shields.io/docker/pulls/feeka94/mcaat" alt="Docker Pulls"/></a>
  <a href="https://github.com/RNABioInfo/mcaat/releases"><img src="https://img.shields.io/github/v/release/RNABioInfo/mcaat" alt="GitHub Release"/></a>
</p>

# metagenomic CRISPR analysis tool - MCAAT v1.0.0


****
Finds CRISPR arrays in raw, un-assembled metagenomic reads. Builds a succinct de Bruijn graph and detects multicycles - the structural signature of CRISPR repeat-spacer arrays - without any prior assembly step.

Outperforms assembly-based workflows and other assembly-free CRISPR detectors on synthetic and real metagenomes.

## Requirements

- CMake ≥ 3.12, C++17, zlib, OpenMP, BZip2
- Docker (recommended for production use)

## Build

```bash
git clone --recurse-submodules https://github.com/RNABioInfo/mcaat.git
cd mcaat

chmod +x ./install.sh
./install.sh
```

The `mcaat` binary will be at `build/mcaat`.

Optional flags:

```bash
./install.sh --install   # also installs to system
./install.sh --clean     # clean build artifacts
```

### Docker

> [!NOTE]
> A pre-built image is available on [Docker Hub](https://hub.docker.com/r/feeka94/mcaat) — no manual dependency setup required.
> ```bash
> docker pull feeka94/mcaat:1.0.0
> ```

```bash
docker build -t mcaat .

docker run --rm -v $(pwd):/data mcaat \
  --input-files /data/reads_R1.fastq /data/reads_R2.fastq \
  --output-folder /data/results
```

The image is based on `debian:bookworm-slim` and ships only the `mcaat` binary and runtime libs (`libomp5`, `zlib1g`).

## Usage

Detailed usage of the tool is outlined: [rnabioinfo.github.io/mcaat](https://rnabioinfo.github.io/mcaat/)

Exactly one input source is required — either raw reads or a pre-built graph:

```bash
# From reads (builds the graph internally)
mcaat --input-files <file1> [file2] [options]

# From a pre-built graph (skips graph construction)
mcaat --graph <path> [options]
```

**Required (one of):**

| Flag | Description |
|---|---|
| `--input-files <file1> [file2]` | One or two FASTA/FASTQ files — plain or gzipped. One file = single-end, two = paired-end |
| `--graph <path>` | Pre-built SDBG graph directory (or file prefix) from a previous run (skips graph construction) |

**Optional:**

| Flag | Default | Description |
|---|---|---|
| `--output-folder <path>` | `mcaat_run_YYYY-MM-DD_HH-MM-SS/` | Output directory |
| `--ram <amount>` | 95% of system RAM | Memory cap. Units: `B`, `K`, `M`, `G` (e.g. `--ram 8G`) |
| `--threads <num>` | CPU cores − 2 | Thread count |
| `--cycle-max-length <int>` | `100` | Inclusive maximum number of distinct nodes per cycle |
| `--cycle-min-length <int>` | `47` | Inclusive minimum number of distinct nodes per cycle |
| `--threshold-multiplicity <int>` | `20` | Min edge multiplicity for cycle start nodes |
| `--repeat-support <fraction>` | `0.85` | Minimum fraction of distinct cycles supporting the entire repeat path; range `(0, 1]` |
| `--array-filter-mode <conservative\|strict>` | `conservative` | Flag uncertain sequence/overlap evidence; opt into strict heuristic rejection |
| `--spacer-variant-similarity <fraction>` | `0.90` | Minimum normalized edit similarity within each spacer family; every variant is preserved |
| `--max-spacer-similarity <fraction>` | `0.60` | Flag mean similarity above this value; reject only in strict mode |
| `--max-repeat-spacer-similarity <fraction>` | `0.60` | Flag repeat/spacer similarity above this value; reject only in strict mode |
| `--max-cycles <int>` | `500` | Exclusive enumeration limit; merged candidates at this limit are flagged (rejected in strict mode) |
| `--low-abundance <true\|false>` | `true` | Enable low-abundance mode |
| `--autoclean <true\|false>` | `true` | Remove intermediate graph/cycle files after run. Set to `false` to keep them |
| `--settings <path>` | — | Key=value settings file (CLI flags override it) |
| `--help`, `-h` | — | Show usage and exit |

## Output

```
<output-folder>/
├── CRISPR_Arrays_1.txt  # detected arrays (split into numbered files if large)
├── CRISPR_graph_systems.tsv    # repeat support and original start IDs per candidate
├── CRISPR_graph_paths.tsv      # observed cycles, repeat/spacer/boundary vertices and edges
├── CRISPR_graph_candidates.tsv # assigned and unresolved candidate cycles
├── CRISPR_graph_summary.tsv    # filtering and deduplication counts
├── CRISPR_spacer_variants.tsv   # observed variants and their family representatives
├── CRISPR_graph_decisions.tsv   # acceptance/rejection reasons for final candidates
├── CRISPR_graph_alternatives.tsv # rejected or competing observations, including paths
├── CRISPR_graph_review_cycles.tsv # non-retained assigned/unresolved paths
├── CRISPR_graph_shared_cycles.tsv # shared evidence retained across competing candidates
├── graph/                     # graph files, unless autoclean removes them
└── cycles.txt                 # ordered k-mers, unless autoclean removes them
```

Each `CRISPR_Arrays_N.txt` file has a short header followed by one block per array:

```
# MCAAT graph-based CRISPR array candidates
# Repeat support: 0.85 (distinct cycles per candidate group; not sequence identity)
# Rows: all distinct observed spacers, without inferred genomic order. Similarity families and paths are in TSV files.

>Array_1  spacers=3  spacer_families=3  flags=-
ATCGATCGATCGATCGATCGATCG
        --------------------    AACCCGGTTAATCGATCGTTTCGAGC
        --------------------    TTGGCCAATCGATCGATCAAAACGGG
        ATCGATCGATCGATCTATCG    GGAATTCCAATCGATCGAATACCCAC   ← repeat variant
```

The supported observed repeat sequence is on its own line. Each distinct observed spacer is shown with its observed repeat variant (or dashes when it matches the supported repeat) followed by the spacer sequence. All observed variants remain in the variant and graph-path tables. `spacers=` counts distinct observed spacer strings; `spacer_families=` gives the separate similarity-group count. The text structure stays compatible with existing readers. CAS and phage commands prefer the graph path table and validate its sequences, edges and boundaries against the loaded graph.

Default curation preserves candidates that pass structural and repeat-support checks, while flagging uncertain sequence features and shared evidence. It does not assign each shared cycle exclusively to one candidate. Read [recall and cycle lengths](docs/recall-and-cycle-lengths.md) before comparing candidate counts with annotated arrays.

See [graph-based repeat/spacer filtering](docs/graph-repeat-filtering.md) for the support definition, deduplication rules, output fields and limitations.

To change support without enumerating cycles again, retain the graph and `cycles.txt` (`--autoclean false`):

```bash
mcaat postprocess-arrays --graph /data/graph/graph \
  --cycles /data/run/cycles.txt --output-folder /data/run-support-60 \
  --repeat-support 0.60
```

## Settings file

Pass a `key=value` file with `--settings`. CLI flags override any value from the file.

```
input-files=/data/R1.fastq /data/R2.fastq
ram=128G
threads=26
output-folder=results/run_1
cycle-max-length=100
cycle-min-length=47
threshold-multiplicity=20
low-abundance=true
repeat-support=0.85
array-filter-mode=conservative
autoclean=true
```

`input-files` accepts one or two paths separated by spaces, commas, or semicolons.

## v2.0.0 (planned)
- **CAS detection**: identify and annotate CAS genes flanking detected CRISPR arrays
- **Protospacer detection**: map spacers back to reads/contigs to find protospacer sequences and PAM sites

## Citation

If you use MCAAT please cite: https://academic.oup.com/microlife/article/doi/10.1093/femsml/uqaf016/8205558

Contact: Please write an issue on our GitHub page if any problems occur.

## Optional cached start-node precheck

`--cycle-precheck unitig` enables cached traversal of non-branching paths during
start-node selection. The default remains `bfs`. The settings-file equivalent is
`cycle-precheck=unitig`; the selected mode is saved in the run configuration.

The cache stores each path's endpoint, original step count and minimum
multiplicity. It is bounded to 16,384 entries per worker and lives only during
start-node selection, after graph pruning. A weighted bidirectional search tracks
distances in both directions and checks their sum when the two searches meet.
Cycle enumeration, parallel spacer paths and the original SDBG sequences are
unchanged. Do not reuse a search object after mutating its graph without clearing
its cache.

For compatibility, the CLI mode preserves the existing bidirectional precheck's
effective bound of `2 * (floor(cycle_max_length / 2) + 1)` steps. In particular,
the existing default of 77 permits a 78-step precheck cycle. The underlying
`UnitigCycleSearch::has_cycle` API uses an exact bound; correcting the existing
precheck's length convention is a separate change.

Build and run the independent correctness tests and the actual-SDBG benchmark:

```sh
cmake -S tests -B build/search -DCMAKE_BUILD_TYPE=Release
cmake --build build/search -j 4
ctest --test-dir build/search --output-on-failure
build/search/benchmark_unitig_search /path/to/graph results.csv 0 5 77 quick
```

`0` tests every qualifying candidate; a positive value selects an evenly spaced
sample. `quick` compares the original BFS with the compatible cached search.
Omit it to also measure exact-bound search with jumps disabled, uncached jumps and
cached jumps. Each trial starts with an empty cache, includes its construction
cost, and validates answers. Loading, pruning and scanning times are reported
separately. The benchmark wraps the original precheck and pruning functions from
the current source, rather than reimplementing them.

To compare full single-thread cycle and array outputs, run
`python3 bench/check_pipeline.py build/mcaat /path/to/graph /path/to/new-output`.
The comparison ignores only cycle-line ordering and generated timestamps.

For sanitizer checks without the SDBG benchmark, configure another build with
`-DBUILD_SDBG_BENCHMARK=OFF -DSEARCH_SANITIZERS=ON`.

