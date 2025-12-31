## MCAAT — Metagenomic CRISPR Array Analysis Tool

Usage
-----
```
./mcaat --input-files <file1> [file2] [--ram <amount>] [--threads <num>] [--output-folder <path>] [--help]
```

Quick description
-----------------
MCAAT detects CRISPR arrays in unassembled metagenomic reads using de Bruijn graph cycle detection. 

Examples
--------
Paired-end:
```
./mcaat --input-files reads_R1.fastq reads_R2.fastq --ram 8G --threads 12 --output-folder results
```
Single-end:
```
./mcaat --input-files reads.fastq
```

### Installation using docker
#### Docker Build

```bash
docker build -t mcaat .
```

---

#### Run the Tool Using Docker

Mount your working directory to access input/output files:

```bash
docker run --rm -v $(pwd):/data mcaat \
  --input_files /data/reads_R1.fastq /data/reads_R2.fastq \
  --output-folder /data/results
```

---

#### Final Image Size

The final image is based on `debian:bookworm-slim` and includes only:

- The `mcaat` binary
- Runtime libraries: `libomp5`, `zlib1g`

This keeps the image small and portable.

---

#### Clean Up

To remove the image:

```bash
docker rmi mcaat
```

### Compiling the project

#### Build the Project
To allow ./install.sh make changes, we execute following command:
```bash
chmod +x ./install.sh
```
You can build the project and the working version will be saved in the build folder.
```bash
./install.sh
```
It is also possible to install the library by simply putting the --install flag.
```bash
./install.sh --install
```
To clean up you can use --clean flag.

---
### Command-Line Arguments

Required
- `--input_files <file1> [file2]` — One or two FASTA/FASTQ files (single or paired-end).

Optional
- `--ram <amount>` — Max RAM (e.g. 4G). Default: ~95% of system RAM.
- `--threads <num>` — Number of threads. Default: CPU cores - 2.
- `--output-folder <path>` — Output directory (default: timestamped folder).
- `--help`, `-h` — Show help and exit.

---

### Output

The tool creates the following directory structure inside the specified output folder:

```
<output-folder>/
├── CRISPR_Arrays.txt         # Raw CRISPR array output
```

---

### Examples

| Scenario                     | Command                                                                 |
|-----------------------------|-------------------------------------------------------------------------|
| Paired-end input with custom output | `./mcaat --input_files reads_R1.fastq reads_R2.fastq --ram 8G --threads 12 --output-folder results/my_run` |
| Single-end input with default output | `./mcaat --input_files reads.fastq` <br>Creates a folder like `mcaat_run_2025-07-07_15-30-00/` |


---

#### Notes

- Input files must exist and be accessible.
- If RAM is set below 1 GB or above system capacity, the program will exit with an error.
- If only one input file is provided, the tool assumes single-end data.

---

### Settings file

Create a simple key=value text file (one setting per line) and pass it with `--settings /path/to/file`.

CLI flags override file values.

Example `settings.txt` (must include `input_files`):
```
input_files=/data/sample/1.fastq /data/sample/2.fastq
ram=128G
threads=26
output_folder=results/run1
```

Notes:
- `input_files` accepts one or two paths (space, comma, or semicolon separated).
- CLI flags override file values.

---

#### Requirements

- C++17 compiler
- [RapidFuzz](https://github.com/maxbachmann/rapidfuzz-cpp) (for fuzzy string matching)
- Filesystem support (`<filesystem>`)

---

## Support

If you encounter issues or have questions, feel free to open an issue or write us an email: fikrat.talibli@ibmg.uni-stuttgart.de. If you are using this software please cite this paper: https://academic.oup.com/microlife/article/doi/10.1093/femsml/uqaf016/8205558.
