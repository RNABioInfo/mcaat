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
SUMMARY = ["srr", "status", "detected_protospacers", "retained_protospacers",
           "extended_inputs", "unique_extensions", "unresolved_inputs",
           "search_limited_inputs", "graph", "results", "detail"]


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
    if chosen and paths.is_file():
        chosen.append(paths)
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
                                 "--max-time", "180", url], check=True, stdout=subprocess.PIPE, text=True)
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
                   "--max-time", "3600", "--speed-limit", "1024", "--speed-time", "120",
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


class StageError(RuntimeError):
    pass


class Runner:
    def __init__(self, folder, root, timeout):
        self.folder, self.root, self.timeout = folder, root, timeout
        self.active = None
        for sig in (signal.SIGINT, signal.SIGTERM):
            signal.signal(sig, self.interrupt)

    def interrupt(self, signum, frame):
        if self.active is not None:
            try:
                os.killpg(self.active, signum)
            except ProcessLookupError:
                pass
        raise KeyboardInterrupt

    def say(self, message):
        line = f"[{utc()}] {message}"
        print(line, flush=True)
        with (self.folder / "batch.log").open("a") as f:
            f.write(line + "\n")

    def run(self, srr, stage, command):
        logs = self.folder / srr
        logs.mkdir(exist_ok=True)
        logfile = logs / (stage + ".log")
        self.say(f"{srr}: {stage} — log: {logfile}")
        started, clock, stopped, interrupted = utc(), time.monotonic(), None, False
        with logfile.open("w") as log:
            log.write("Command: " + shlex.join(list(map(str, command))) + "\n")
            log.flush()
            p = subprocess.Popen(list(map(str, command)), stdout=log, stderr=subprocess.STDOUT,
                                 stdin=subprocess.DEVNULL, start_new_session=True)
            self.active = p.pid
            try:
                while True:
                    pid, status, usage = os.wait4(p.pid, os.WNOHANG)
                    if pid:
                        break
                    elapsed = time.monotonic() - clock
                    if stopped is None and ((self.timeout and elapsed > self.timeout) or shutil.disk_usage(self.root).free < GIB):
                        stopped = time.monotonic()
                        log.write("Stopping: stage timeout or less than 1 GiB disk space remains.\n"); log.flush()
                        os.killpg(p.pid, signal.SIGTERM)
                    if stopped is not None and time.monotonic() - stopped > 10:
                        try:
                            os.killpg(p.pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                    time.sleep(0.2)
            except KeyboardInterrupt:
                interrupted = True
                try:
                    os.killpg(p.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                deadline = time.monotonic() + 10
                while True:
                    pid, status, usage = os.wait4(p.pid, os.WNOHANG)
                    if pid:
                        break
                    if time.monotonic() > deadline:
                        try:
                            os.killpg(p.pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                    time.sleep(0.2)
            finally:
                self.active = None
            rc = os.waitstatus_to_exitcode(status)
            p.returncode = rc
        elapsed = time.monotonic() - clock
        rss = usage.ru_maxrss / (1024 ** 2 if sys.platform == "darwin" else 1024)
        append_tsv(self.folder / "metrics.tsv", METRICS, dict(srr=srr, stage=stage,
            started_utc=started, wall_seconds=f"{elapsed:.3f}", user_cpu_seconds=f"{usage.ru_utime:.3f}",
            system_cpu_seconds=f"{usage.ru_stime:.3f}", peak_rss_mib=f"{rss:.2f}", exit_code=rc, log=str(logfile)))
        self.say(f"{srr}: {stage} finished: exit={rc}, {elapsed:.1f}s, peak RSS={rss:.1f} MiB")
        if interrupted:
            raise KeyboardInterrupt
        if rc or stopped is not None:
            raise StageError(f"{stage} failed (exit {rc}); see {logfile}")


def checkpoint(directory, filenames):
    write_json(directory / "DONE.json", {"completed_utc": utc(),
        "sha256": {name: digest(directory / name) for name in filenames}})


def completed(directory):
    try:
        hashes = json.loads((directory / "DONE.json").read_text())["sha256"]
        return bool(hashes) and all(digest(directory / name) == h for name, h in hashes.items())
    except (OSError, ValueError, KeyError):
        return False


def archive_incomplete(directory):
    if directory.exists():
        directory.rename(directory.with_name(directory.name + f".previous-{time.time_ns()}"))


def process(sample, arrays, a, runner, binary_hash):
    srr = sample.name
    base = sample / "protospacer_batch"
    base.mkdir(exist_ok=True)
    graph = None
    for candidate in (sample / "graph/graph", base / "graph/graph"):
        if Path(str(candidate) + ".sdbg_info").exists():
            try:
                graph = graph_info(candidate)
                runner.say(f"{srr}: reusing complete graph {candidate}; no read download needed")
                break
            except (OSError, ValueError) as e:
                runner.say(f"{srr}: existing graph is incomplete/incompatible; preserving it: {e}")
    if graph is None:
        free_space(base, a.min_free_gb * GIB)
        reads = base / "reads"
        runner.run(srr, "download_reads", [sys.executable, __file__, "_fetch", srr, reads, int(a.min_free_gb * GIB)])
        free_space(base, a.min_free_gb * GIB)
        stage = base / ("graph_build_" + runner.folder.name)
        stage.mkdir()
        library_config(reads, stage / "reads.lib")
        runner.run(srr, "build_read_library", [a.mcaat, "buildlib", stage / "reads.lib", stage / "reads"])
        free_space(base, a.min_free_gb * GIB)
        runner.run(srr, "build_graph", [a.mcaat, "read2sdbg", "--read_lib_file", stage / "reads",
            "--output_prefix", stage / "graph", "--host_mem", int(a.ram_gb * GIB),
            "--num_cpu_threads", a.threads, "-k", 23, "-m", 1])
        graph_info(stage / "graph")
        archive_incomplete(base / "graph")
        stage.rename(base / "graph")
        graph = graph_info(base / "graph/graph")
        write_json(base / "graph/build_manifest.json", {"graph": graph, "binary_sha256": binary_hash,
            "reads": json.loads((reads / "reads.json").read_text()), "ram_gb": a.ram_gb, "threads": a.threads})
    # File identity, settings and exact array contents identify a resumable job.
    inputs = {p.name: digest(p) for p in arrays}
    config = {"version": 1, "arrays": inputs, "binary_sha256": binary_hash, "helper_sha256": digest(__file__),
              "graph": graph, "threads": a.threads, "beam_width": a.beam_width,
              "max_states": a.max_states, "max_query_states": a.max_query_states,
              "max_retained_bases": a.max_retained_bases, "max_flank_length": 750, "min_extension_flank": 750}
    signature = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()[:16]
    job = base / ("run_" + signature)
    job.mkdir(exist_ok=True)
    if completed(job) and completed(job / "discovery") and (
            not (job / "extensions").exists() or completed(job / "extensions")):
        row = json.loads((job / "summary.json").read_text())
        row["status"] = "SKIPPED_COMPLETE"
        runner.say(f"{srr}: already completed with identical inputs/settings: {job}")
        return row
    snapshot = job / "arrays"
    snapshot.mkdir(exist_ok=True)
    for path in arrays:
        name = "CRISPR_Arrays_legacy.txt" if path.name == "CRISPR_Arrays.txt" else path.name
        shutil.copy2(path, snapshot / name)
        if digest(snapshot / name) != inputs[path.name]:
            raise RuntimeError("array input changed while copying: " + str(path))
    write_json(job / "job.json", config)
    common = [a.mcaat, "phage-curate", "--arrays", snapshot, "--graph", graph["prefix"],
        "--threads", a.threads, "--beam-width", a.beam_width, "--source-context-width", 256,
        "--max-node-visits", 3, "--max-length", 1550, "--max-flank-length", 750,
        "--max-states", a.max_states, "--max-query-states", a.max_query_states,
        "--max-retained-bases", a.max_retained_bases, "--source-array-crossing", "true",
        "--correct-spacer-boundaries", "true", "--resolve-graph-ties", "false"]
    discovery, filtered = job / "discovery", job / "extensions"
    if not completed(discovery):
        free_space(base, a.min_free_gb * GIB)
        archive_incomplete(discovery)
        runner.run(srr, "discover_protospacers", common + ["--output", discovery])
        read_fasta(discovery / "Protospacers.fasta")
        checkpoint(discovery, ["Protospacers.fasta", "protospacer_parameters.json", "protospacer_queries.tsv"])
    else:
        runner.say(f"{srr}: resuming after completed discovery")
    detected = read_fasta(discovery / "Protospacers.fasta")
    row = dict(srr=srr, status="DONE" if detected else "NO_DETECTIONS", detected_protospacers=len(detected),
               retained_protospacers=len(detected), extended_inputs=0, unique_extensions=0,
               unresolved_inputs=0, search_limited_inputs=0, graph=graph["prefix"], results=str(job), detail="")
    if detected:
        if not completed(filtered):
            free_space(base, a.min_free_gb * GIB)
            archive_incomplete(filtered)
            runner.run(srr, "filter_extensions", common + ["--output", filtered,
                "--protospacers", discovery / "Protospacers.fasta", "--path-quality-filter", "true",
                "--min-extension-flank", 750, "--prediction-score-margin", 0.1])
            if read_fasta(filtered / "Protospacers.fasta") != detected:
                raise RuntimeError("retention check failed: protospacer IDs/sequences changed")
            params = json.loads((filtered / "protospacer_parameters.json").read_text())
            if params.get("path_quality_version") != 3 or not params.get("path_quality_filter"):
                raise RuntimeError("the rebuilt binary did not use quality filter version 3")
            checkpoint(filtered, ["Protospacers.fasta", "ContiguousPaths.fasta", "protospacer_parameters.json",
                                  "protospacer_extensions.tsv", "protospacer_path_quality.tsv"])
        params = json.loads((filtered / "protospacer_parameters.json").read_text())
        row.update(retained_protospacers=len(read_fasta(filtered / "Protospacers.fasta")),
                   extended_inputs=params["extended_records"], unique_extensions=len(read_fasta(filtered / "ContiguousPaths.fasta")),
                   unresolved_inputs=params["not_extended_records"], search_limited_inputs=params["search_limited_records"])
    if graph_info(graph["prefix"]) != graph:
        raise RuntimeError("graph files changed during prediction")
    write_json(job / "summary.json", row)
    files = ["summary.json", "job.json", "discovery/DONE.json", "discovery/Protospacers.fasta"]
    if detected:
        files += ["extensions/DONE.json", "extensions/Protospacers.fasta", "extensions/ContiguousPaths.fasta"]
    checkpoint(job, files)
    (base / "LATEST.txt").write_text(str(job) + "\n")
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("srr", nargs="*", help="optional SRR IDs; default: all eligible folders")
    parser.add_argument("--root", type=Path, default=Path("/vol/d/data/copsac"))
    parser.add_argument("--mcaat", type=Path, default=HERE.parent / "build/mcaat")
    parser.add_argument("--threads", type=int, default=24)
    parser.add_argument("--ram-gb", type=float, help="graph builder's memory budget; default 70%% of available RAM")
    parser.add_argument("--min-free-gb", type=float, default=20, help="free-disk reserve before each major stage")
    parser.add_argument("--stage-timeout-hours", type=float, default=0, help="0 disables stage timeout")
    parser.add_argument("--beam-width", type=int, default=128)
    parser.add_argument("--max-states", type=int, default=200000)
    parser.add_argument("--max-query-states", type=int, default=800000)
    parser.add_argument("--max-retained-bases", type=int, default=500000000,
                        help="native retained-sequence budget per invocation (not an RSS limit)")
    parser.add_argument("--dry-run", action="store_true")
    a = parser.parse_args()
    a.root, a.mcaat = a.root.resolve(), a.mcaat.resolve()
    if not a.root.is_dir():
        parser.error("root directory does not exist")
    if not 2 <= a.threads <= 1024 or not 1 <= a.beam_width <= 1024 or min(a.min_free_gb, a.stage_timeout_hours) < 0:
        parser.error("invalid thread, beam, disk or timeout setting")
    if not all(0 < n <= 1000000000 for n in (a.max_states, a.max_query_states, a.max_retained_bases)):
        parser.error("search and retention budgets must be between 1 and 1000000000")
    if any(not re.fullmatch(r"SRR\d+", s) for s in a.srr):
        parser.error("sample arguments must be SRR accessions")
    requested = set(a.srr)
    samples = [p for p in sorted(a.root.iterdir()) if p.is_dir() and re.fullmatch(r"SRR\d+", p.name)
               and (not requested or p.name in requested)]
    missing = requested - {p.name for p in samples}
    if missing:
        parser.error("sample directories missing: " + ", ".join(sorted(missing)))
    eligible = [(p, array_files(p)) for p in samples if array_files(p)]
    print(f"Eligible samples: {len(eligible)}; skipped without nonempty array reports: {len(samples)-len(eligible)}", flush=True)
    if a.dry_run:
        for p, arrays in eligible:
            print(p.name + "\t" + ",".join(f.name for f in arrays))
        return 0
    if not eligible:
        return 0
    if not a.mcaat.is_file() or not os.access(a.mcaat, os.X_OK):
        parser.error("mcaat binary missing; rebuild or pass --mcaat /absolute/path/to/mcaat")
    if not shutil.which("curl"):
        parser.error("curl is required")
    check = subprocess.run([str(a.mcaat), "buildlib", "--help"], capture_output=True, text=True)
    if check.returncode or "Usage: mcaat buildlib" not in check.stdout:
        parser.error("mcaat is outdated: pull this branch and rebuild before starting the batch")
    if a.ram_gb is None:
        try:
            mem = dict(line.split(":", 1) for line in Path("/proc/meminfo").read_text().splitlines())
            a.ram_gb = int(mem["MemAvailable"].split()[0]) * 1024 * .7 / GIB
        except (OSError, KeyError):
            parser.error("pass --ram-gb explicitly on systems without /proc/meminfo")
    if not 1 < a.ram_gb < 100000:
        parser.error("graph builder requires --ram-gb greater than 1")
    logs = a.root / "_mcaat_batch_logs"
    logs.mkdir(exist_ok=True)
    with (logs / "batch.lock").open("a+") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error("another batch is already using this root")
        run_id = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ") + f"-{os.getpid()}"
        logdir = logs / run_id
        logdir.mkdir()
        (logs / "LATEST.txt").write_text(str(logdir) + "\n")
        runner = Runner(logdir, a.root, a.stage_timeout_hours * 3600)
        write_json(logdir / "configuration.json", {k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()})
        binary_hash = digest(a.mcaat)
        write_json(logdir / "provenance.json", {"binary_sha256": binary_hash, "helper_sha256": digest(__file__)})
        runner.say(f"Processing {len(eligible)} SRRs sequentially, {a.threads} threads each; graph budget {a.ram_gb:.1f} GiB")
        failures = 0
        for sample, arrays in eligible:
            try:
                free_space(a.root, a.min_free_gb * GIB)
            except RuntimeError as e:
                runner.say(f"Batch stopped; completed work preserved: {e}")
                return 1
            try:
                if digest(a.mcaat) != binary_hash:
                    raise RuntimeError("binary changed during the batch; restart with a fixed binary")
                row = process(sample, arrays, a, runner, binary_hash)
            except KeyboardInterrupt:
                append_tsv(logdir / "summary.tsv", SUMMARY, dict(srr=sample.name, status="INTERRUPTED", detail="partial files preserved"))
                runner.say("Interrupted; partial downloads and outputs preserved. Rerun the same command to resume.")
                return 130
            except Exception as e:
                failures += 1
                row = dict(srr=sample.name, status="FAILED", detail=str(e))
                runner.say(f"{sample.name}: FAILED: {e}; continuing with next SRR")
            append_tsv(logdir / "summary.tsv", SUMMARY, row)
            runner.say(f"{sample.name}: {row['status']}; detected={row.get('detected_protospacers', '-')}, unique extensions={row.get('unique_extensions', '-')}")
        runner.say(f"Batch finished: {len(eligible)-failures} completed/skipped, {failures} failed. Summary: {logdir/'summary.tsv'}")
        return int(failures > 0)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "_fetch":
        fetch_reads(sys.argv[2], sys.argv[3], int(sys.argv[4]))
    else:
        sys.exit(main())
