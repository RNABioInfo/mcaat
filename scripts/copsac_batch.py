#!/usr/bin/env python3
"""Sequential, resumable array -> graph -> protospacer/extension batch.

Only Python's standard library, curl and the repository's rebuilt mcaat are
required. ENA API: https://ena-docs.readthedocs.io/en/latest/retrieval/
programmatic-access/file-reports.html . No reference/BLAST is used for prediction.
"""
import argparse
import csv
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time
from urllib.parse import urlencode, urlparse

HERE = Path(__file__).resolve().parent
GIB = 1024 ** 3
METRICS = ["srr", "stage", "started_utc", "wall_seconds", "user_cpu_seconds",
           "system_cpu_seconds", "peak_rss_mib", "exit_code", "log"]
SUMMARY = ["srr", "status", "spacer_observations", "distinct_spacers", "candidate_protospacer_regions",
           "protospacer_regions", "candidate_only_regions",
           "unique_extensions", "limited_queries", "graph", "results", "detail"]


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def digest(path, kind="sha256"):
    h = hashlib.new(kind)
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, data):
    path = Path(path)
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    temp.replace(path)


def append_tsv(path, fields, row):
    new = not path.exists()
    with path.open("a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        if new:
            w.writeheader()
        w.writerow(row)


def array_files(sample):
    numbered = sorted(p for p in sample.glob("CRISPR_Arrays_*.txt") if p.is_file() and p.stat().st_size)
    legacy = sample / "CRISPR_Arrays.txt"
    chosen = numbered or ([legacy] if legacy.is_file() and legacy.stat().st_size else [])
    # Match the native reader: graph paths take precedence over text reports.
    paths = sample / "CRISPR_graph_paths.tsv"
    if paths.is_file():
        return [paths]
    return chosen


def graph_info(prefix):
    """Validate every bucket's extent against actual shard sizes, without loading graph RAM."""
    prefix = Path(prefix)
    meta = Path(str(prefix) + ".sdbg_info")
    with meta.open() as f:
        headers = {}
        for expected in ("k", "words_per_tip_label", "num_buckets", "num_files"):
            name, value = f.readline().split()
            if name != expected:
                raise ValueError("invalid graph header: " + str(meta))
            headers[name] = int(value)
        if headers["k"] != 23 or not 0 < headers["num_files"] <= 4096:
            raise ValueError("expected a complete k=23 graph: " + str(prefix))
        if headers["words_per_tip_label"] != 2 or not 0 < headers["num_buckets"] <= 65536:
            raise ValueError("unexpected graph layout: " + str(prefix))
        offsets, seen, items, rows = {}, set(), 0, 0
        for line in f:
            bucket, fid, start, count, tips, large = map(int, line.split())
            rows += 1
            if bucket == 2 ** 64 - 1 and fid == 2 ** 64 - 1:
                continue
            if not (0 <= bucket < headers["num_buckets"] and 0 <= fid < headers["num_files"]):
                raise ValueError("invalid bucket/file id")
            if bucket in seen or min(start, count, tips, large) < 0 or max(tips, large) > count:
                raise ValueError("invalid graph bucket")
            if start != offsets.get(fid, 0):
                raise ValueError("noncontiguous graph shard")
            seen.add(bucket)
            # SdbgItem/mul_t are uint16_t; tip-label words are uint32_t.
            offsets[fid] = start + 2 * count + 2 * large + 4 * headers["words_per_tip_label"] * tips
            items += count
        if rows != headers["num_buckets"] or not items:
            raise ValueError("truncated or empty graph metadata")
    shards = {}
    for fid, size in offsets.items():
        path = Path(str(prefix) + ".sdbg." + str(fid))
        st = path.stat()
        if st.st_size != size:
            raise ValueError("incomplete graph shard: " + str(path))
        shards[path.name] = [st.st_size, st.st_mtime_ns]
    return {"prefix": str(prefix.resolve()), "k": 23, "items": items,
            "metadata_sha256": digest(meta), "shards_size_mtime": shards}


def read_fasta(path):
    records, name = {}, None
    with Path(path).open() as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in records:
                    raise ValueError("duplicate FASTA ID: " + name)
                records[name] = ""
            elif line.strip():
                if name is None or re.search("[^ACGTacgt]", line.strip()):
                    raise ValueError("invalid DNA FASTA: " + str(path))
                records[name] += line.strip()
    if any(not s for s in records.values()):
        raise ValueError("empty FASTA sequence: " + str(path))
    return records


def free_space(path, minimum):
    free = shutil.disk_usage(path).free
    if free < minimum:
        raise RuntimeError(f"only {free/GIB:.1f} GiB free; need {minimum/GIB:.1f} GiB reserve")


def fetch_reads(srr, folder, reserve):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    report = folder / "ena.json"
    if not report.exists():
        url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urlencode({
            "accession": srr, "result": "read_run", "format": "json",
            "fields": "run_accession,fastq_ftp,fastq_md5,fastq_bytes"})
        result = subprocess.run(["curl", "--fail", "--location", "--silent", "--show-error",
                                 "--retry", "5", "--retry-delay", "5", "--connect-timeout", "30",
                                 url], check=True, stdout=subprocess.PIPE, text=True)
        records = json.loads(result.stdout)
    else:
        records = json.loads(report.read_text())
    if len(records) != 1 or records[0].get("run_accession") != srr:
        raise ValueError("ENA returned no unique matching run: " + srr)
    r = records[0]
    urls, hashes, sizes = (r.get(k, "").split(";") for k in ("fastq_ftp", "fastq_md5", "fastq_bytes"))
    if not urls[0] or not len(urls) == len(hashes) == len(sizes):
        raise ValueError("ENA has no complete FASTQ file report for " + srr)
    entries = []
    for url, checksum, size in zip(urls, hashes, sizes):
        url = "https://" + url.removeprefix("ftp://").removeprefix("https://")
        parsed = urlparse(url)
        name = Path(parsed.path).name
        if parsed.hostname != "ftp.sra.ebi.ac.uk" or not re.fullmatch(r"SRR\d+(?:_[12])?\.fastq\.gz", name):
            raise ValueError("unexpected ENA FASTQ location: " + url)
        if not name.startswith(srr + ".") and not name.startswith(srr + "_"):
            raise ValueError("FASTQ accession differs from requested run")
        if not re.fullmatch(r"[0-9a-fA-F]{32}", checksum) or int(size) <= 0:
            raise ValueError("invalid ENA checksum/size")
        entries.append({"url": url, "md5": checksum.lower(), "bytes": int(size), "path": str(folder / name)})
    if len({e["path"] for e in entries}) != len(entries):
        raise ValueError("duplicate ENA FASTQ files")
    write_json(report, records)
    for e in entries:
        dest = Path(e["path"])
        if dest.exists():
            if dest.stat().st_size == e["bytes"] and digest(dest, "md5") == e["md5"]:
                print("Verified cached FASTQ:", dest, flush=True)
                continue
            raise ValueError("existing FASTQ failed checksum; preserved: " + str(dest))
        partial = dest.with_name(dest.name + ".part")
        for attempt in range(1, 4):
            needed = max(0, e["bytes"] - (partial.stat().st_size if partial.exists() else 0))
            free_space(folder, reserve + needed)
            print(f"Downloading {srr}: {dest.name}, attempt {attempt}", flush=True)
            cmd = ["curl", "--fail", "--location", "--silent", "--show-error", "--retry", "5",
                   "--retry-delay", "5", "--retry-connrefused", "--connect-timeout", "30",
                   "--speed-limit", "1024", "--speed-time", "120",
                   "--continue-at", "-", "--output", str(partial),
                   "--write-out", "HTTP %{http_code}; downloaded %{size_download} bytes; %{speed_download} bytes/s\n", e["url"]]
            result = subprocess.run(cmd)
            valid = partial.exists() and partial.stat().st_size == e["bytes"] and digest(partial, "md5") == e["md5"]
            if valid:
                partial.replace(dest)
                print("MD5 verified:", dest.name, flush=True)
                break
            if result.returncode == 0 or result.returncode in (33, 36) or (partial.exists() and partial.stat().st_size >= e["bytes"]):
                if partial.exists():
                    partial.rename(partial.with_name(partial.name + f".bad-{time.time_ns()}"))
            if attempt == 3:
                raise RuntimeError("download/checksum failed; partial file preserved: " + str(dest))
        e["verified_utc"] = utc()
    write_json(folder / "reads.json", entries)


def library_config(reads_folder, destination):
    entries = json.loads((Path(reads_folder) / "reads.json").read_text())
    paths = [Path(e["path"]).resolve() for e in entries]
    if any(not p.is_file() or re.search(r"\s", str(p)) for p in paths):
        raise ValueError("graph-library input paths must exist and contain no whitespace")
    by_name = {p.name: p for p in paths}
    used, lines = set(), []
    for p in paths:
        if p.name.endswith("_1.fastq.gz"):
            mate = by_name.get(p.name.replace("_1.fastq.gz", "_2.fastq.gz"))
            if mate:
                lines += ["# paired reads", f"pe {p} {mate}"]
                used.update((p, mate))
    for p in paths:
        if p not in used:
            lines += ["# single/unpaired reads", f"se {p}"]
    Path(destination).write_text("\n".join(lines) + "\n")


ACTIVE_PID = None


def interrupt(signum, frame):
    """Forward interruption to the current stage and preserve its partial output."""
    if ACTIVE_PID is not None:
        try:
            os.killpg(ACTIVE_PID, signum)
        except ProcessLookupError:
            pass
    raise KeyboardInterrupt


def run_stage(command, folder, metrics, srr, stage):
    """Run one stage without a wall timeout; capture its own CPU/RSS accounting."""
    global ACTIVE_PID
    command = [str(x) for x in command]
    log_path = folder / (stage + ".log")
    (folder / (stage + ".command.txt")).write_text(shlex.join(command) + "\n")
    started, stamp = time.monotonic(), utc()
    interrupted = False
    with log_path.open("w") as log:
        child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        ACTIVE_PID = child.pid
        try:
            _, status, usage = os.wait4(child.pid, 0)
        except KeyboardInterrupt:
            interrupted = True
            _, status, usage = os.wait4(child.pid, 0)
        finally:
            ACTIVE_PID = None
        child.returncode = os.waitstatus_to_exitcode(status)
    peak = usage.ru_maxrss / (1024**2 if sys.platform == "darwin" else 1024)
    append_tsv(metrics, METRICS, dict(srr=srr, stage=stage, started_utc=stamp,
        wall_seconds=time.monotonic()-started, user_cpu_seconds=usage.ru_utime,
        system_cpu_seconds=usage.ru_stime, peak_rss_mib=peak, exit_code=child.returncode, log=str(log_path)))
    if interrupted:
        raise KeyboardInterrupt
    if child.returncode:
        raise RuntimeError(f"{stage} exited {child.returncode}; see {log_path}")


def existing_graph_prefix(sample):
    """Locate graph metadata in original, current-workflow or legacy-batch storage."""
    for prefix in (sample/"graph/graph", sample/"protospacer_workflow/graph/graph", sample/"protospacer_batch/graph/graph"):
        if Path(str(prefix)+".sdbg_info").exists():
            return prefix
    return None


def prepare_graph(sample, job, options, metrics):
    """Reuse a validated graph, or download reads and build one with bundled MCAAT code."""
    base = sample / "protospacer_workflow"
    prefix = existing_graph_prefix(sample)
    if prefix is not None:
        return graph_info(prefix)
    if options.existing_only:
        raise ValueError("no existing graph for " + sample.name)
    reads = base/"reads"
    run_stage([sys.executable, __file__, "_fetch", sample.name, reads, int(options.min_free_gb*GIB)],
              job, metrics, sample.name, "fetch_reads")
    free_space(base, options.min_free_gb*GIB)
    build = job/"graph_build"
    build.mkdir()
    library_config(reads, build/"reads.lib")
    run_stage([options.mcaat, "buildlib", build/"reads.lib", build/"reads"],
              job, metrics, sample.name, "build_read_library")
    run_stage([options.mcaat, "read2sdbg", "--read_lib_file", build/"reads", "--output_prefix", build/"graph",
               "--host_mem", int(options.ram_gb*GIB), "--num_cpu_threads", options.threads, "-k", 23, "-m", 1],
              job, metrics, sample.name, "build_graph")
    graph_info(build/"graph")
    if (base/"graph").exists():
        raise ValueError("existing incomplete graph directory preserved: " + str(base/"graph"))
    build.rename(base/"graph")
    return graph_info(base/"graph/graph")


def validate_predictions(folder):
    """Require complete v9 output and one accepted extension for every final target."""
    summary = json.loads((folder / "run_summary.json").read_text())
    if summary["status"] != "complete" or summary["implementation"] != "protospacer-workflow-v9" or \
       summary.get("protospacer_output_policy") != "accepted_extensions_only":
        raise ValueError("native workflow did not complete with the final-target policy")
    targets = read_fasta(folder / "Protospacers.fasta")
    candidates = read_fasta(folder / "CandidateProtospacers.fasta")
    extensions = read_fasta(folder / "ContiguousPaths.fasta")
    if len(targets) != summary["protospacer_regions"] or len(extensions) != summary["unique_extensions"] or \
       len(candidates) != summary["candidate_protospacer_regions"] or \
       len(candidates) - len(targets) != summary["candidate_only_regions"]:
        raise ValueError("FASTA count differs from native summary")
    if any(candidates.get(key) != sequence for key, sequence in targets.items()):
        raise ValueError("final target is absent from candidate audit")
    with (folder / "associations.tsv").open() as stream:
        links = list(csv.DictReader(stream, delimiter="\t"))
    if {r["target_id"] for r in links} != set(targets) or \
       {r["extension_id"] for r in links} != set(extensions) or \
       len({r["query_id"] for r in links}) != summary["extended_queries"]:
        raise ValueError("final target/extension associations are incomplete")
    for row in links:
        sequence = extensions[row["extension_id"]]
        begin, end = int(row["seed_begin"]), int(row["seed_end"])
        fragment = sequence[begin:end]
        canonical = min(fragment, fragment.translate(str.maketrans("ACGT", "TGCA"))[::-1])
        if not 0 <= begin < end <= len(sequence) or canonical != targets[row["target_id"]]:
            raise ValueError("final target is not present at its extension coordinates")
    return summary


def process(sample, arrays, options, logdir, binary_hash):
    """Snapshot arrays and execute exactly one native detection/extension/filtering invocation."""
    if options.existing_only and existing_graph_prefix(sample) is None:
        return dict(srr=sample.name, status="SKIPPED_NO_GRAPH", results="",
                    detail="No existing graph; downloading and building are disabled.")
    free_space(options.root, options.min_free_gb*GIB)
    base = sample/"protospacer_workflow"
    base.mkdir(exist_ok=True)
    job = base/dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d_%H-%M-%S")
    job.mkdir()
    manifest = {"status": "running", "sample": sample.name, "implementation": "protospacer-workflow-v9",
                "binary_sha256": binary_hash, "helper_sha256": digest(__file__), "arrays": {p.name: digest(p) for p in arrays}}
    write_json(job/"job.json", manifest)
    try:
        graph = prepare_graph(sample, job, options, logdir/"measurements.tsv")
        manifest.update(graph=graph, threads=options.threads, beam_width=options.beam_width,
                        max_flank_length=options.max_flank_length, min_extension_flank=options.min_extension_flank,
                        min_external_bases=options.min_external_bases, read_length=options.read_length,
                        max_context_paths=options.max_context_paths, max_coverage_fold=options.max_coverage_fold)
        # Resume only a completed invocation with exactly the same inputs/settings.
        identity = {k: v for k, v in manifest.items() if k != "status"}
        if not options.rerun:
            for previous in sorted(base.glob("*/job.json"), reverse=True):
                if previous.parent == job:
                    continue
                old = json.loads(previous.read_text())
                if old.get("status") == "complete" and {k: v for k, v in old.items() if k != "status"} == identity:
                    result = previous.parent/"predictions/run_summary.json"
                    summary = validate_predictions(result.parent)
                    manifest["status"] = "skipped_identical_complete"
                    write_json(job/"job.json", manifest)
                    return dict(summary, srr=sample.name, status="SKIPPED_COMPLETE", results=str(result.parent))
        snapshot = job/"arrays"
        snapshot.mkdir()
        for path in arrays:
            shutil.copy2(path, snapshot/path.name)
            if digest(snapshot/path.name) != manifest["arrays"][path.name]:
                raise ValueError("array input changed while copying")
        write_json(job/"job.json", manifest)
        if digest(options.mcaat) != binary_hash:
            raise ValueError("mcaat binary changed during batch; completed work preserved")
        run_stage([options.mcaat, "phage-curate", "--arrays", snapshot, "--graph", graph["prefix"],
                   "--output", job/"predictions", "--threads", options.threads, "--beam-width", options.beam_width,
                   "--source-array-crossing", "true", "--max-context-paths", options.max_context_paths,
                   "--read-length", options.read_length, "--min-external-bases", options.min_external_bases,
                   "--max-coverage-fold", options.max_coverage_fold,
                   "--max-flank-length", options.max_flank_length, "--min-extension-flank", options.min_extension_flank],
                  job, logdir/"measurements.tsv", sample.name, "protospacer_workflow")
        summary = validate_predictions(job/"predictions")
        if graph_info(graph["prefix"]) != graph:
            raise ValueError("graph changed during prediction")
        manifest["status"] = "complete"
        write_json(job/"job.json", manifest)
        (base/"LATEST.txt").write_text(str(job/"predictions") + "\n")
        status = "DONE" if summary["protospacer_regions"] else (
            "NO_ACCEPTED_EXTENSIONS" if summary["candidate_protospacer_regions"] else "NO_DETECTIONS")
        return dict(summary, srr=sample.name, status=status,
                    results=str(job/"predictions"))
    except BaseException as error:
        manifest["status"] = "interrupted" if isinstance(error, KeyboardInterrupt) else "failed"
        manifest["error"] = str(error)
        write_json(job/"job.json", manifest)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("srr", nargs="*", help="SRRs to process; default all folders with array reports")
    parser.add_argument("--root", type=Path, default=Path("/vol/d/data/copsac"))
    parser.add_argument("--mcaat", type=Path, default=HERE.parent/"build/mcaat")
    parser.add_argument("--threads", type=int, default=28)
    parser.add_argument("--ram-gb", type=float, default=192, help="graph-build budget in GiB")
    parser.add_argument("--min-free-gb", type=float, default=20)
    parser.add_argument("--beam-width", type=int, default=64)
    parser.add_argument("--max-flank-length", type=int, default=750)
    parser.add_argument("--min-extension-flank", type=int, default=200)
    parser.add_argument("--min-external-bases", type=int, default=200)
    parser.add_argument("--max-context-paths", type=int, default=256)
    parser.add_argument("--read-length", type=int, default=150)
    parser.add_argument("--max-coverage-fold", type=float, default=2.,
                        help="maximum ratio of read-span coverage medians; 0 disables")
    parser.add_argument("--existing-only", action="store_true", help="reuse graphs and skip samples without one; do not download or build")
    parser.add_argument("--rerun", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    args.root, args.mcaat = args.root.resolve(), args.mcaat.resolve()
    if not args.root.is_dir() or any(not re.fullmatch(r"SRR\d+", s) for s in args.srr):
        parser.error("root must exist; sample names must be SRR accessions")
    if not 1 <= args.threads <= 1024 or not 1 <= args.beam_width <= 256 or args.ram_gb <= 1 or args.min_free_gb < 0:
        parser.error("invalid resource settings")
    if not 1 <= args.min_extension_flank <= args.max_flank_length <= 10000:
        parser.error("invalid flank lengths")
    if not 1 <= args.min_external_bases <= args.max_flank_length or not 1 <= args.max_context_paths <= 4096 or not 25 <= args.read_length <= 10000:
        parser.error("invalid external evidence or context settings")
    if args.max_coverage_fold != 0 and not 1 < args.max_coverage_fold <= 1000000:
        parser.error("coverage fold must exceed one, or be zero to disable")
    samples = [p for p in sorted(args.root.iterdir()) if p.is_dir() and re.fullmatch(r"SRR\d+", p.name)
               and (not args.srr or p.name in args.srr)]
    missing = set(args.srr)-{p.name for p in samples}
    if missing:
        parser.error("missing folders: " + ", ".join(sorted(missing)))
    eligible = [(p, array_files(p)) for p in samples if array_files(p)]
    if args.dry_run:
        for sample, arrays in eligible:
            print(sample.name + "\t" + ",".join(p.name for p in arrays))
        return 0
    help_result = subprocess.run([str(args.mcaat), "phage-curate", "--help"], text=True, capture_output=True)
    if help_result.returncode or "workflow v9" not in help_result.stdout:
        parser.error("rebuild MCAAT with native protospacer workflow v9")
    logs = args.root/"_mcaat_workflow_logs"
    logs.mkdir(exist_ok=True)
    for signum in (signal.SIGINT, signal.SIGTERM):
        signal.signal(signum, interrupt)
    with (logs/"batch.lock").open("a+") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error("another workflow batch is already using this root")
        logdir = logs/dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d_%H-%M-%S")
        logdir.mkdir()
        binary_hash = digest(args.mcaat)
        write_json(logdir/"configuration.json", {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()})
        failures = 0
        for sample, arrays in eligible:
            print("Processing " + sample.name, flush=True)
            try:
                row = process(sample, arrays, args, logdir, binary_hash)
            except KeyboardInterrupt:
                append_tsv(logdir/"summary.tsv", SUMMARY, dict(srr=sample.name, status="INTERRUPTED", detail="partial outputs preserved"))
                return 130
            except Exception as error:
                failures += 1
                row = dict(srr=sample.name, status="FAILED", detail=str(error))
            append_tsv(logdir/"summary.tsv", SUMMARY, row)
            print(f"{sample.name}: {row['status']}; candidates={row.get('candidate_protospacer_regions', '-')}; "
                  f"final protospacers={row.get('protospacer_regions', '-')}; extensions={row.get('unique_extensions', '-')}", flush=True)
        print("Summary: " + str(logdir/"summary.tsv"), flush=True)
        return int(failures > 0)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "_fetch":
        fetch_reads(sys.argv[2], sys.argv[3], int(sys.argv[4]))
    else:
        sys.exit(main())
