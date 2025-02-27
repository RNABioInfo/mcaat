# MCAAT - metagenomic CRISPR array analysis tool

## Compilation
To compile the software, ensure you have CMake (minimum version 3.12) installed. Then, follow these steps:

```bash
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

After compilation, you can run the tool using:

```bash
./mcaat --input_files <file1.fastq> [file2.fastq] [options]
```

Alternatively, you can build and run the tool using Docker.

### Using Docker
#### Build the Docker Image
```bash
docker build -t mcaat_image .
```

#### Run the Docker Container
```bash
docker run --rm -v $(pwd):/mcaat mcaat_image --input_files <file1.fastq> [file2.fastq] [options]
```

## Required Arguments
- `--input_files <file1.fastq> [file2.fastq]`  
  - You must provide **one or two** `.fastq` input files.  
  - If an input file does not exist, the program will terminate with an error.

## Optional Arguments
- `--ram <value>`  
  - Sets the amount of RAM (in MiB) to be used.  
  - **Default:** 90% of total system memory.

- `--threads <num>`  
  - Number of threads to use.  
  - **Default:** Total available CPU cores minus 2. If the device has only 2 both will be used.

- `--k <k-mer>`  
  - Specifies the k-mer size.

- `--output-folder <path>`  
  - Specifies the directory for storing output files.
> **IMPORTANT**
> MCAAT saves the graph under output_folder/graph. If a user stops the software early(by pressing CMD+C), the folder will not be deleted.   

## Example Usage
```bash
./mcaat --input_files sample1.fastq sample2.fastq --ram 4000 --threads 8 --output-folder results/
```

Or using Docker:
```bash
docker run --rm -v $(pwd):/mcaat mcaat_image --input_files sample1.fastq sample2.fastq --ram 4000 --threads 8 --output-folder results/
```

This command:
- Uses `sample1.fastq` and `sample2.fastq` as input files.
- Allocates 4000 MiB of RAM.
- Runs with 8 threads.
- Saves outputs to the `results/` directory.

If an option is not provided, **mcaat** applies default values where applicable.

