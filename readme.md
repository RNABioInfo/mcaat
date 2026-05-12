<p align="left">
  <img src="icon.png?v=2" alt="MCAAT" width="350"/>
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
git clone --recurse-submodules https://github.com/feeka94/mcaat.git
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

```bash
docker build -t mcaat .

docker run --rm -v $(pwd):/data mcaat \
  --input-files /data/reads_R1.fastq /data/reads_R2.fastq \
  --output-folder /data/results
```

The image is based on `debian:bookworm-slim` and ships only the `mcaat` binary and runtime libs (`libomp5`, `zlib1g`).

## Usage

```bash
mcaat --input-files <file1> [file2] [options]
```

| Flag | Description |
|---|---|
| `--input-files <file1> [file2]` | One or two FASTA/FASTQ files - plain or gzipped. One file = single-end, two files = paired-end |
| `--graph <path>` | Path to a pre-built SDBG graph folder (skips graph construction) |
| `--ram <amount>` | Memory cap. Units: `B`, `K`, `M`, `G` (e.g. `--ram 8G`). Default: 95% of system RAM |
| `--threads <num>` | Thread count. Default: CPU cores − 2 |
| `--output-folder <path>` | Output directory. Default: timestamped folder `mcaat_run_YYYY-MM-DD_HH-MM-SS/` |
| `--benchmark <file>` | File with expected CRISPR sequences (one per line) for evaluation |
| `--cycle-max-length <int>` | Maximum cycle length to search. Default: 77 |
| `--cycle-min-length <int>` | Minimum cycle length to search. Default: 27 |
| `--threshold-multiplicity <int>` | Min edge multiplicity for cycle start nodes. Default: 20 |
| `--low-abundance <true\|false>` | Enable low-abundance mode. Default: true |
| `--settings <path>` | Key=value settings file (CLI flags override it) |
| `--help`, `-h` | Show usage and exit |

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
# ──────────────────────────────────────────────────────

>Array_1  spacers=8
ATCGATCGATCGATCGATCGATCG
        --------------------    AACCCGGTTAATCGATCG
        --------------------    TTGGCCAATCGATCGATC
        ATCGATCGATCGATCTATCG    GGAATTCCAATCGATCGA   ← repeat variant
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
```

`input-files` accepts one or two paths separated by spaces, commas, or semicolons.

## Source

```
CMakeLists.txt
install.sh
Dockerfile
libs/
  megahit/                 MEGAHIT sdbg (submodule)
  spoa/                    sequence alignment (submodule)
  kseqpp/                  FASTA/FASTQ I/O (submodule)
include/
  core/
    settings.h             global settings struct + file parser
    pipeline.h             production pipeline steps (header-only)
    sdbg_build.h
    tmp_utils.h
    main_run_and_debug.h   debug-only scaffolding
  crispr_array/
    cycle_finder.h
    filters.h
    post_processing.h      array consensus + output (header-only)
    evaluation.h
    spacer_ordering.h
    graph_generic_func.h
    sequence_generic_func.h
    node_counter.h
  io/
    reads.h
    io_ops.h
    path_writer.h
  protospacer_detection/
    isolate_protospacers.h
    phage_curator.h
  cas/                     (v2.0.0)
src/
  main.cpp                 CLI + argument parsing
  core/
    sdbg_build.cpp
    tmp_utils.cpp
    main_run_and_debug.cpp
  crispr_array/
    cycle_finder.cpp       multicycle detection (parallel DFS)
    filters.cpp
    evaluation.cpp
    spacer_ordering.cpp
    graph_generic_func.cpp
    sequence_generic_func.cpp
  io/
    reads.cpp
    io_ops.cpp
    path_writer.cpp
  protospacer_detection/
    isolate_protospacers.cpp
    phage_curator.cpp
  cas/                     (v2.0.0)
tests/                     unit tests
docs/                      algorithmic notes and optimization report
```

## Roadmap

### v2.0.0 (planned)
- **CAS detection**: identify and annotate CAS genes flanking detected CRISPR arrays
- **Protospacer detection**: map spacers back to reads/contigs to find protospacer sequences and PAM sites

## Citation

If you use MCAAT please cite: https://academic.oup.com/microlife/article/doi/10.1093/femsml/uqaf016/8205558

Contact: fikrat.talibli@ibmg.uni-stuttgart.de

