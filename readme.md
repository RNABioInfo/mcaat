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
| `--cycle-max-length <int>` | `77` | Maximum cycle length to search |
| `--cycle-min-length <int>` | `27` | Minimum cycle length to search |
| `--threshold-multiplicity <int>` | `20` | Min edge multiplicity for cycle start nodes |
| `--low-abundance <true\|false>` | `true` | Enable low-abundance mode |
| `--autoclean <true\|false>` | `true` | Remove intermediate graph/cycle files after run. Set to `false` to keep them |
| `--settings <path>` | — | Key=value settings file (CLI flags override it) |
| `--help`, `-h` | — | Show usage and exit |

## Output

```
<output-folder>/
├── CRISPR_Arrays_1.txt  # detected arrays (split into numbered files if large)
├── graph/               # succinct de Bruijn graph files
└── cycles/              # raw cycle data
```

Each `CRISPR_Arrays_N.txt` file has a short header followed by one block per array:

```
# MCAAT — CRISPR Array Output
# Generated : 2026-05-12 10:30:21
# Arrays    : 42
# Spacers   : 312

>Array_1  spacers=8
ATCGATCGATCGATCGATCGATCG
        --------------------    AACCCGGTTAATCGATCGTTTCGAGC
        --------------------    TTGGCCAATCGATCGATCAAAACGGG
        ATCGATCGATCGATCTATCG    GGAATTCCAATCGATCGAATACCCAC   ← repeat variant
```

The consensus repeat sequence is on its own line. Each spacer entry shows the repeat variant (or dashes when it matches the consensus exactly) followed by the spacer sequence.

## Settings file

Pass a `key=value` file with `--settings`. CLI flags override any value from the file.

```
input-files=/data/R1.fastq /data/R2.fastq
ram=128G
threads=26
output-folder=results/run_1
cycle-max-length=77
cycle-min-length=27
threshold-multiplicity=20
low-abundance=true
autoclean=true
```

`input-files` accepts one or two paths separated by spaces, commas, or semicolons.

## v2.0.0 (planned)
- **CAS detection**: identify and annotate CAS genes flanking detected CRISPR arrays.
- **Protospacer workflow (experimental)**: one native class performs connected
  spacer alignment, source-context and sequence-artifact filtering, bounded beam
  extension and output. Local branch support ranks actual graph routes; count
  variation breaks ties. Source-repeat crossing is available only for an already
  detected or explicitly supplied target. The workflow saves one selected flank,
  records competing-route uncertainty and search limits, and preserves input
  associations with an optional path audit. It uses no reference database during
  prediction. Discovery and extension reject the exact shared TruSeq prefix as
  well as longer adapter cores, preventing technical branches from establishing
  false outside-array context. After selection, `--max-coverage-fold 2` withholds
  extensions with inconsistent downstream block medians, preserving their
  diagnostic records in `withheld_extensions.tsv`. Final `Protospacers.fasta`
  now contains only targets with an accepted extension, linked through
  `associations.tsv` to `ContiguousPaths.fasta`. All detected candidates remain
  separately in `CandidateProtospacers.fasta` and `candidate_associations.tsv`.
  Supplied inputs remain unchanged in `SuppliedProtospacers.fasta`. Use
  `--max-coverage-fold 0` to disable the coverage check. On SRR14017553 it removed 35 nonmatching
  extensions and one matching extension; the separate 81-input check retained
  all 81 records and all 43 previously matching extensions. The full discovery
  run now saves 117 final protospacers linked to 106 unchanged accepted
  extensions. See the [output-policy validation](docs/protospacer-evaluation-SRR14017553-v9.md),
  [coverage-filter evaluation](docs/protospacer-evaluation-SRR14017553-v8.md)
  and [workflow and function reference](docs/protospacer-extension.md).

## Citation

If you use MCAAT please cite: https://academic.oup.com/microlife/article/doi/10.1093/femsml/uqaf016/8205558

Contact: Please write an issue on our GitHub page if any problems occur.

