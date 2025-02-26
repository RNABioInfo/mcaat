# Command-Line Tool Usage

## Compilation

After compiling the software, you can run the tool using the following command:

```bash
./your_program --input_files <file1> [file2] [options]
```

## Required Arguments

- `--input_files <file1> [file2]`
  - You must provide **one or two** input files - **paired-end or signle-end** reads.
  - If an input file does not exist, the program will terminate with an error.

## Optional Arguments

- `--ram <value>`

  - Sets the amount of RAM (in MiB) to be used.
  - **Default:** 90% of total system memory.

- `--threads <num>`

  - Number of threads to use.
  - **Default:** Total available CPU cores minus 2. If you have only 2 cores software takes both.

- `--graph-folder <path>`

  - Specifies the directory where graph-related files will be stored. In future might be omitted.

- `--cycles-folder <path>`

  - Specifies the directory where cycle-related files will be stored. In future might be omitted.

- `--output-folder <path>`

  - Specifies the directory for storing output files.

## Example Usage

```bash
./mcaat --input_files R1.fastq R2.fastq --ram 4000 --threads 8 --output-folder results/
```

This command:

- Uses `data1.fastq` and `data2.fastq` as input files.
- Allocates 4000 MiB of RAM.
- Runs with 8 threads.
- Saves outputs to the `results/` directory.

If an option is not provided, the tool applies default values where applicable.

