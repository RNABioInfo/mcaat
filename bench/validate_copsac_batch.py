#!/usr/bin/env python3
"""Native graph build + two-pass batch, failure continuation, restart and legacy inputs.

Usage: python3 bench/validate_copsac_batch.py build/mcaat /path/to/new/test-output
Uses cached synthetic ENA FASTQs; never downloads real sequencing reads.
"""
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import shutil
import subprocess

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("batch", REPO / "scripts/copsac_batch.py")
batch = importlib.util.module_from_spec(spec); spec.loader.exec_module(batch)
p = argparse.ArgumentParser(); p.add_argument("binary", type=Path); p.add_argument("output", type=Path)
a = p.parse_args(); binary = a.binary.resolve(); out = a.output.resolve()
out.mkdir(parents=True, exist_ok=True)
if any(out.iterdir()): raise SystemExit("use an empty output directory")
root = out / "copsac"; root.mkdir()
rng = random.Random(314159)
def dna(n): return "".join(rng.choice("ACGT") for _ in range(n))
repeat = dna(29); spacers = [dna(35) for _ in range(3)]
report = ">Array_1\n" + repeat + "\n" + "".join("        " + "-"*len(repeat) + "\t" + s + "\n" for s in spacers)
bad, good, skip = (root / s for s in ("SRR000001", "SRR000002", "SRR000003"))
for d in (bad, good, skip): d.mkdir()
(bad / "CRISPR_Arrays.txt").write_text(report)
(bad / "FAILED").write_text("previous run failed\n")
(good / "CRISPR_Arrays_1.txt").write_text(report)
(skip / "CRISPR_Arrays_1.txt").touch()
for sample in (bad, good): (sample / "protospacer_batch/reads").mkdir(parents=True)
(bad / "protospacer_batch/reads/ena.json").write_text(json.dumps([dict(
    run_accession=bad.name, fastq_ftp="", fastq_md5="", fastq_bytes="")]))
seqs = [dna(120)+repeat+repeat.join(spacers)+repeat+dna(120), dna(1000)+spacers[0]+dna(1000)]
reads = good / "protospacer_batch/reads/SRR000002.fastq.gz"
with gzip.open(reads, "wt") as f:
    for number, seq in enumerate(seqs):
        for pos in range(len(seq)-149):
            for copy in range(3):
                f.write(f"@r{number}_{pos}_{copy}\n{seq[pos:pos+150]}\n+\n" + "I"*150 + "\n")
(reads.parent / "ena.json").write_text(json.dumps([dict(run_accession=good.name,
    fastq_ftp="ftp.sra.ebi.ac.uk/vol1/SRR000002.fastq.gz", fastq_md5=batch.digest(reads,"md5"),
    fastq_bytes=str(reads.stat().st_size))]))
cmd = ["bash", str(REPO / "scripts/run_copsac_protospacers.sh"), "--root", str(root),
       "--mcaat", str(binary), "--threads", "2", "--ram-gb", "2", "--min-free-gb", "1",
       "--beam-width", "16", "--max-states", "20000", "--max-query-states", "80000"]
def run(name, expected):
    with (out / (name + ".log")).open("w") as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, timeout=180)
    assert result.returncode == expected, (name, result.returncode)
    logs = Path((root / "_mcaat_batch_logs/LATEST.txt").read_text().strip())
    with (logs / "summary.tsv").open() as f: summary = {r["srr"]:r for r in csv.DictReader(f, delimiter="\t")}
    return logs, summary
logs, rows = run("first", 1)
assert rows[bad.name]["status"] == "FAILED" and rows[good.name]["status"] == "DONE"
assert skip.name not in rows
job = Path(rows[good.name]["results"])
detected = batch.read_fasta(job / "discovery/Protospacers.fasta")
assert len(detected) > 0 and batch.read_fasta(job / "extensions/Protospacers.fasta") == detected
assert len(batch.read_fasta(job / "extensions/ContiguousPaths.fasta")) > 0
assert not list(root.rglob("cycles.txt")), "graph-only build must not run cycle search"
graph_dir = good / "protospacer_batch/graph"
before = {f.name: batch.digest(f) for f in graph_dir.glob("graph.*")}
logs2, rows2 = run("restart", 1)
assert rows2[good.name]["status"] == "SKIPPED_COMPLETE"
assert not (logs2 / good.name).exists(), "completed sample must not rerun stages"
# Repair the first sample by providing an existing graph. The bad cached
# download is irrelevant now, and CRISPR_Arrays.txt must work through snapshotting.
shutil.copytree(graph_dir, bad / "graph")
logs3, rows3 = run("existing_graph_and_legacy_arrays", 0)
assert rows3[bad.name]["status"] == "DONE" and rows3[good.name]["status"] == "SKIPPED_COMPLETE"
assert not (logs3 / bad.name / "download_reads.log").exists()
assert (bad / "FAILED").read_text() == "previous run failed\n"
assert (bad / "CRISPR_Arrays.txt").read_text() == report
assert {f.name: batch.digest(f) for f in graph_dir.glob("graph.*")} == before
empty = root / "SRR000004"; empty.mkdir()
(empty / "CRISPR_Arrays_1.txt").write_text(">Array_absent\n" + dna(29) + "\n        " + "-"*29 + "\t" + dna(45) + "\n")
shutil.copytree(graph_dir, empty / "graph")
logs4, rows4 = run("zero_detections", 0)
assert rows4[empty.name]["status"] == "NO_DETECTIONS"
assert rows4[empty.name]["detected_protospacers"] == "0"
assert not (logs4 / empty.name / "filter_extensions.log").exists()
with (logs / "metrics.tsv").open() as f: metrics = list(csv.DictReader(f, delimiter="\t"))
assert {r["stage"] for r in metrics if r["srr"] == good.name} == {
    "download_reads", "build_read_library", "build_graph", "discover_protospacers", "filter_extensions"}
result = dict(status="passed", detected=len(detected),
    extensions=len(batch.read_fasta(job / "extensions/ContiguousPaths.fasta")),
    checks=["graph-only construction", "two-pass native prediction", "all detections retained",
            "failed sample does not stop next sample", "checksum-verified cached download",
            "restart skips unchanged completed outputs", "existing graph skips download/build",
            "unnumbered legacy array report", "old FAILED marker retried", "graph bytes unchanged",
            "five stage metrics present", "zero detections skip extension successfully"])
(out / "validation.json").write_text(json.dumps(result, indent=2)+"\n")
print(json.dumps(result, indent=2))
