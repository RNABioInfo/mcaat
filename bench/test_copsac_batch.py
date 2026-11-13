#!/usr/bin/env python3
"""Small offline tests for overnight batch failure/resume behavior."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location("batch", Path(__file__).resolve().parents[1] / "scripts/copsac_batch.py")
batch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(batch)


class BatchTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name).resolve()
        self.handlers = {s: signal.getsignal(s) for s in (signal.SIGINT, signal.SIGTERM)}

    def tearDown(self):
        for s, handler in self.handlers.items():
            signal.signal(s, handler)
        self.tmp.cleanup()

    def test_reports_do_not_duplicate_legacy_and_retry_failed(self):
        p = self.root / "SRR123"; p.mkdir()
        (p / "FAILED").touch()
        legacy = p / "CRISPR_Arrays.txt"; legacy.write_text("legacy\n")
        self.assertEqual(batch.array_files(p), [legacy])
        numbered = p / "CRISPR_Arrays_1.txt"; numbered.write_text("numbered\n")
        self.assertEqual(batch.array_files(p), [numbered])
        paths = p / "CRISPR_graph_paths.tsv"; paths.write_text("header\n")
        self.assertEqual(batch.array_files(p), [numbered, paths])
        empty = self.root / "SRR456"; empty.mkdir(); (empty / "CRISPR_Arrays_1.txt").touch()
        self.assertEqual(batch.array_files(empty), [])

    def test_graph_extent_checks_catch_partial_shard(self):
        prefix = self.root / "graph"
        meta = self.root / "graph.sdbg_info"
        meta.write_text("k 23\nwords_per_tip_label 2\nnum_buckets 2\nnum_files 1\n0 0 0 3 1 1\n18446744073709551615 18446744073709551615 0 0 0 0\n")
        shard = self.root / "graph.sdbg.0"; shard.write_bytes(bytes(16))
        self.assertEqual(batch.graph_info(prefix)["items"], 3)
        shard.write_bytes(bytes(15))
        with self.assertRaises(ValueError): batch.graph_info(prefix)
        shard.write_bytes(bytes(16)); meta.write_text(meta.read_text().rsplit("\n", 2)[0] + "\n")
        with self.assertRaises(ValueError): batch.graph_info(prefix)

    def test_download_resumes_and_checks_md5(self):
        payload = b"fixture downloaded bytes"
        folder = self.root / "reads"; folder.mkdir()
        dest = folder / "SRR123.fastq.gz"
        partial = folder / "SRR123.fastq.gz.part"; partial.write_bytes(payload[:7])
        report = [{"run_accession": "SRR123", "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/SRR123.fastq.gz",
                   "fastq_md5": hashlib.md5(payload).hexdigest(), "fastq_bytes": str(len(payload))}]
        (folder / "ena.json").write_text(json.dumps(report))
        def curl(cmd, **kw):
            self.assertIn("--continue-at", cmd)
            self.assertEqual(partial.read_bytes(), payload[:7])
            partial.write_bytes(payload)
            return subprocess.CompletedProcess(cmd, 0)
        with patch.object(batch.subprocess, "run", side_effect=curl) as mock:
            batch.fetch_reads("SRR123", folder, 0)
            self.assertEqual(mock.call_count, 1)
        self.assertEqual(dest.read_bytes(), payload)
        with patch.object(batch.subprocess, "run", side_effect=AssertionError("cached file should not download")):
            batch.fetch_reads("SRR123", folder, 0)
        dest.write_bytes(b"corrupted")
        with self.assertRaises(ValueError): batch.fetch_reads("SRR123", folder, 0)
        self.assertEqual(dest.read_bytes(), b"corrupted", "user-visible corrupted file must be preserved")

    def test_bad_download_is_never_promoted(self):
        folder = self.root / "reads"; folder.mkdir()
        expected = b"correct"
        (folder / "ena.json").write_text(json.dumps([dict(run_accession="SRR123",
            fastq_ftp="ftp.sra.ebi.ac.uk/vol1/SRR123.fastq.gz", fastq_md5=hashlib.md5(expected).hexdigest(),
            fastq_bytes=str(len(expected)))]))
        def curl(cmd, **kw):
            Path(cmd[cmd.index("--output")+1]).write_bytes(b"wrong!!")
            return subprocess.CompletedProcess(cmd, 0)
        with patch.object(batch.subprocess, "run", side_effect=curl):
            with self.assertRaises(RuntimeError): batch.fetch_reads("SRR123", folder, 0)
        self.assertFalse((folder / "SRR123.fastq.gz").exists())
        self.assertEqual(len(list(folder.glob("*.bad-*"))), 3)

    def test_library_retains_unpaired_reads_alongside_pair(self):
        paths = [self.root / name for name in ("SRR123_1.fastq.gz", "SRR123_2.fastq.gz", "SRR123.fastq.gz")]
        for p in paths: p.write_bytes(b"fixture")
        (self.root / "reads.json").write_text(json.dumps([dict(path=str(p)) for p in paths]))
        output = self.root / "reads.lib"; batch.library_config(self.root, output)
        self.assertEqual(output.read_text().splitlines(), ["# paired reads", f"pe {paths[0]} {paths[1]}",
                                                          "# single/unpaired reads", f"se {paths[2]}"])

    def test_completion_requires_unchanged_outputs(self):
        p = self.root / "output.txt"; p.write_text("original")
        batch.checkpoint(self.root, [p.name])
        self.assertTrue(batch.completed(self.root))
        p.write_text("modified")
        self.assertFalse(batch.completed(self.root))

    def test_failed_command_has_metrics_and_log(self):
        runner = batch.Runner(self.root, self.root, 0)
        with self.assertRaises(batch.StageError):
            runner.run("SRR123", "failure", [sys.executable, "-c", "print('test failure'); raise SystemExit(7)"])
        with (self.root / "metrics.tsv").open() as f: rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(rows[0]["exit_code"], "7")
        self.assertGreater(float(rows[0]["wall_seconds"]), 0)
        self.assertIn("test failure", (self.root / "SRR123/failure.log").read_text())


if __name__ == "__main__":
    unittest.main()
