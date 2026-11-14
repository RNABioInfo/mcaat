# Sequential COPSAC protospacer batch

Run from the repository on the Linux server after rebuilding this branch:

```bash
cd /vol/d/mcaat-protospacer/mcaat-protospacer-full
git pull --ff-only origin codex/protospacer-extension
cmake --build build -j28
nohup bash scripts/run_copsac_protospacers.sh --threads 28 --ram-gb 192 \
  > /vol/d/data/copsac/protospacer_batch.console.log 2>&1 < /dev/null &
```

The Bash entry point uses a Python 3.9+ standard-library helper and `curl`.
No separate MEGAHIT or SRA Toolkit installation is needed. The rebuilt MCAAT
exposes its bundled `buildlib` and `read2sdbg` functions as subcommands; these
never enter array/cycle detection or the normal pipeline's autocleanup.

Use `--dry-run` to list eligible folders without downloading or changing files.
To restrict the batch, append IDs, for example `SRR14017550 SRR14017552`.
`--root` defaults to `/vol/d/data/copsac`; `--mcaat` defaults to `build/mcaat`
relative to this repository, independent of the current working directory.

After the discovery correction, rerun with the rebuilt binary and updated script.
Old jobs are preserved for comparison; their detections must not be reused as
supplied targets. The changed code creates a fresh results job and reuses complete
saved graphs. Downloads and graph construction are skipped for those samples.

## What happens for each SRR

1. Select immediate `SRR[digits]` folders containing a nonempty
   `CRISPR_Arrays_*.txt` or `CRISPR_Arrays.txt`. Existing `FAILED` markers do not
   exclude a sample with arrays. If numbered reports exist, the unnumbered
   report is not duplicated. `CRISPR_graph_paths.tsv`, when present, takes
   precedence inside MCAAT, matching its normal input semantics.
2. Validate and reuse `SRR.../graph/graph` if complete and k=23. Validation checks
   every metadata bucket and shard byte extent; it does not establish which
   reads originally generated a user-provided graph. If that graph is missing
   or incomplete, retain it and use/build `SRR.../protospacer_batch/graph/graph`.
3. For missing graphs, fetch the SRR FASTQs over HTTPS from ENA, including any
   unpaired file as well as the paired files. Resume partial transfers and
   verify published sizes and MD5 checksums. An existing unverified FASTQ in
   the original SRR directory is not silently substituted for the archive data.
   Complete existing graphs require no download.
4. Build the read library and k=23, minimum-edge-frequency=1 graph. These are
   graph-only commands: array detection and cycle search are not rerun.
   Keep graphs, downloaded reads, library files, logs and failed partial work.
5. Discover protospacers using the saved arrays with `--source-array-crossing false`.
   A path through the source repeat does not establish a separate target, even
   if it subsequently reaches a long, stable flank. Then feed its
   `Protospacers.fasta` into a separate `--protospacers` invocation using quality
   filter version 3 and `--source-array-crossing true`. Only this supplied-target
   extension pass permits a single source-array passage. The CLI rejects crossing
   in discovery mode. This second pass is necessary because array-discovery mode
   does not apply the supplied-target quality filter. Retaining supplied IDs is
   not itself evidence that they are separate targets.
6. Verify that every detected ID and sequence survives the second pass, report
   extension counts, and move to the next SRR. Zero detections is a successful
   `NO_DETECTIONS` outcome; the empty FASTA is not submitted for extension.

The batch uses 28 threads within one SRR and processes SRRs sequentially. Defaults are
beam 128, source-context pool 256, up to 750 added bases per side, and at least
750 added bases on one side for a filtered extension. Search is bounded at
200,000 states per side and 800,000 per query. The sequence-retention budget is
500,000,000 bases per MCAAT invocation; this is not a process memory cap.
Prediction does not use reference sequences or BLAST.

Graph construction defaults to 70% of available host RAM. The launch command
above uses a fixed 192 GiB budget for the user's 256 GB server. Set `--ram-gb`
to choose another budget. This setting applies to graph
construction, not extension-process RSS. Use an explicit budget under scheduler
or container memory limits. The batch checks for 20 GiB of free disk before
major stages, and stops an active stage if free disk falls below 1 GiB.
`--min-free-gb` changes the pre-stage reserve. Optional
`--stage-timeout-hours 6` stops a stage after six hours; the default has no
overall stage timeout. Network operations have their own timeouts/retries.

## Outputs and measurements

Each sample has:

```text
SRR.../protospacer_batch/
  reads/                      downloaded gzip files, ENA metadata, MD5 manifest
  graph/graph.sdbg_info        newly built graph, when needed
  graph/graph.sdbg.*           graph shards, kept permanently by this script
  LATEST.txt                  absolute path of the latest completed result
  run_<input-and-settings-hash>/
    arrays/                   snapshot of supplied array reports
    discovery/Protospacers.fasta
    extensions/Protospacers.fasta
    extensions/ContiguousPaths.fasta
    extensions/protospacer_path_quality.tsv
    extensions/protospacer_extensions.tsv
    job.json
    summary.json
    DONE.json
```

`extensions/` is absent when discovery found zero protospacers. Existing files
such as `ContiguousPaths.fasta`, `protospacer_results*` and the original array
reports are not overwritten.

Batch logs are under `/vol/d/data/copsac/_mcaat_batch_logs/<UTC-time>-<PID>/`:

- `batch.log`: human-readable progress and errors.
- `summary.tsv`: per-SRR status, detected/retained targets, extended inputs,
  distinct extension sequences, unextended inputs, search-limit counts and paths.
  The `unresolved_inputs` column counts all inputs without an emitted extension,
  not just ambiguity cases. These are predictions, not BLAST match counts.
- `metrics.tsv`: stage/function wall seconds, user CPU seconds, system CPU
  seconds, peak RSS in MiB, exit code and log location. CPU time can exceed wall
  time with multiple threads. RSS is the OS's per-process high-water accounting,
  not a sum of all simultaneously live processes.
- `<SRR>/download_reads.log`, `build_read_library.log`, `build_graph.log`,
  `discover_protospacers.log`, `filter_extensions.log`: exact executed command
  plus complete stdout/stderr. The graph builder also reports internal sorting
  and construction stage timings.

The measured functions are these five pipeline stages. The script does not
claim to profile every internal C++ call; discovery and filtering logs include
MCAAT's native progress and diagnostic tables.

To monitor:

```bash
tail -f /vol/d/data/copsac/protospacer_batch.console.log
cat "$(cat /vol/d/data/copsac/_mcaat_batch_logs/LATEST.txt)/summary.tsv"
```

Rerun the same command to resume. Completed output hashes are verified before
skipping stages; changed array files, executable, helper or search settings get
a separate result directory. Graphs are reused across runs. Interrupted outputs
are retained under `.previous-*`; failed download chunks stay as `.part` files.
A failed sample is logged and the next SRR runs. Low disk space stops the batch.
A root-level advisory lock prevents two copies of this batch from overlapping.

## Verification

`python3 bench/test_copsac_batch.py` exercises discovery of eligible samples,
complete/truncated graph validation, checksum/resume behavior, paired plus
unpaired library construction, completion hashes and process accounting.
`python3 bench/validate_copsac_batch.py build/mcaat /path/to/new/test-output`
builds a small native graph and exercises both prediction stages, a failed first
sample followed by a successful one, rerun/resume, existing-graph reuse and an
unnumbered legacy array report. Its FASTQs are synthetic and cached locally;
no real metagenomic reads are downloaded by this test. It also requires exactly
one detection from three source spacers: only one has a separate synthetic target.
The other two must be rejected despite having traversable source-array flanks.
