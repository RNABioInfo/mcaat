<p align="center">
  <img src="icon.png?v=2" alt="MCAAT" width="350"/>
</p>

# MCAAT - v1.0.0

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
├── CRISPR_Arrays.txt    # detected arrays: repeat + spacers per system
├── graph/               # succinct de Bruijn graph files
└── cycles/              # raw cycle data
```

## Settings file

Pass a `key=value` file with `--settings`. CLI flags override any value from the file.

```
input_files=/data/R1.fastq /data/R2.fastq
ram=128G
threads=26
output_folder=results/run_1
cycle_max_length=77
cycle_min_length=27
threshold_multiplicity=20
low_abundance=true
```

`input_files` accepts one or two paths separated by spaces, commas, or semicolons.

## Source

```
CMakeLists.txt
install.sh
Dockerfile
libs/megahit/              MEGAHIT sdbg (submodule)
libs/spoa/                 sequence alignment (submodule)
libs/kseqpp/               FASTA/FASTQ I/O (submodule)
libs/rapidfuzz-cpp/        fuzzy string matching (submodule)
include/                   headers
src/
  main.cpp                 CLI + argument parsing
  sdbg_build.cpp           de Bruijn graph construction
  cycle_finder.cpp         multicycle detection (parallel DFS)
  post_processing.h        CRISPR array extraction and consensus
  pipeline.h               production pipeline steps
include/                   headers
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

