# Sequential native MCAAT batch

Read the [experimental method’s measured results](protospacer-evaluation-SRR14017553-v5.md)
before using it broadly. A working batch does not establish prediction accuracy.

The batch runs **one native v5 workflow invocation per SRR**. There is no separate
supplied-protospacer preservation pass. Prediction uses no third-party assembler,
aligner or solver. Python standard-library code orchestrates the bundled MCAAT graph
builder; curl is used only for ENA downloads if a graph is missing.

For the 28-thread / 256-GB server:

```bash
bash scripts/run_copsac_protospacers.sh \
  --root /vol/d/data/copsac --threads 28 --ram-gb 192
```

For existing graphs only, with explicit dataset selection:

```bash
bash scripts/run_copsac_protospacers.sh \
  --root /vol/d/data/copsac --threads 28 --existing-only \
  SRR14017553 SRR14017554 SRR14017569
```

Pass `--mcaat /absolute/path/to/build/mcaat` if the binary is elsewhere. The script
rejects older implementations. `--dry-run` lists inputs without starting work.
Do not infer cross-dataset accuracy from the single local SRR14017553 test.

## Lifecycle

1. Find SRR folders with array observations. Graph-path TSV takes precedence over
   numbered reports; legacy reports are a fallback.
2. Validate and reuse the existing graph when available. If absent, download ENA
   FASTQs with size/checksum validation, then call MCAAT's bundled `buildlib` and
   `read2sdbg`. `--existing-only` disables fetching/building.
3. Snapshot arrays into a readable timestamped job directory and record provenance.
4. Call native `phage-curate` once for detection, extension, filtering and saving.
5. Verify final status/counts and unchanged graph metadata/shards. Retain graph,
   reads, logs, inputs and results. Failed stages preserve partial files.
6. Continue sequentially to the next SRR. Identical completed jobs are skipped;
   `--rerun` requests a fresh job. Changed binary, helper, arrays, graph or settings
   prevent reuse of an old completed job.

No wall-clock stage timeout is imposed. Network connection/stall checks allow failed
transfers to retry. Native search budgets are finite and reported; they are not timeouts.
A root-level lock prevents two batches using the same folders concurrently. Signals
are forwarded to the active stage and recorded as interruptions.

## Where results and logs go

- `SRR.../protospacer_workflow/YYYY-MM-DD_HH-MM-SS/predictions/`: final native files.
- `SRR.../protospacer_workflow/LATEST.txt`: latest successful prediction folder.
- `job.json`: input/binary/helper identity and complete/failed/interrupted status.
- Per-stage `.command.txt` and `.log`: exact commands and output.
- `_mcaat_workflow_logs/<timestamp>/summary.tsv`: one status/count row per sample.
- `measurements.tsv`: per-process wall/CPU time, peak RSS, exit status and log location.
- Native `query_metrics.tsv` and `run_summary.json`: per-query and stage measurements.

`Protospacers.fasta` counts target intervals; `ContiguousPaths.fasta` counts extended
regions. Neither count is a count of complete phage genomes. Empty FASTAs with a
complete native summary mean no detections; a failed/interrupted job is different.

## Every helper function

| Function | Responsibility |
|---|---|
| `utc` | Format UTC timestamps. |
| `digest` | Incrementally hash inputs or downloads. |
| `write_json` | Atomically replace a JSON status/provenance file. |
| `append_tsv` | Append a row, adding the header for a new file. |
| `array_files` | Select one authoritative array-input representation. |
| `graph_info` | Validate k=23 metadata, bucket extents and complete shard sizes; record identity. |
| `read_fasta` | Validate DNA records and count final sequences. |
| `free_space` | Check the configured disk reserve before major work. |
| `fetch_reads` | Obtain ENA file reports, resume transfers and verify sizes/MD5 checksums. |
| `library_config` | Describe paired and unpaired downloaded reads for bundled graph construction. |
| `interrupt` | Forward a signal to the active stage. |
| `run_stage` | Run and account for one process without a wall timeout. |
| `prepare_graph` | Reuse a complete graph or build a new retained graph. |
| `process` | Snapshot one SRR, resume compatible completed work, run one native workflow and validate outputs. |
| `main` | Validate CLI settings, lock the root, iterate samples and write the batch summary. |

The shell entry point only locates and executes this Python helper. No helper classes
or secondary prediction engines are used.

`tests/check_workflow_batch.py` has one `main` function. It checks authoritative
input selection, valid empty results, invalid duplicate IDs, successful stage logs,
failed-stage propagation and per-stage measurements without downloads or graph builds.
