#!/usr/bin/env python3
"""Infrastructure checks only; no biological dataset, downloads, or graph builds."""
import importlib.util
from pathlib import Path
import sys
import tempfile


def main():
    spec = importlib.util.spec_from_file_location("batch", Path(__file__).parents[1]/"scripts/copsac_batch.py")
    batch = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(batch)
    with tempfile.TemporaryDirectory(prefix="mcaat-workflow-check-") as directory:
        root = Path(directory)
        legacy, numbered, paths = (root/name for name in ("CRISPR_Arrays.txt", "CRISPR_Arrays_1.txt", "CRISPR_graph_paths.tsv"))
        legacy.write_text("legacy report\n")
        assert batch.array_files(root) == [legacy]
        numbered.write_text("numbered report\n")
        assert batch.array_files(root) == [numbered]
        paths.write_text("graph paths\n")
        assert batch.array_files(root) == [paths], "input versions must never be combined"
        empty = root/"empty.fasta"
        empty.write_text("")
        assert batch.read_fasta(empty) == {}, "completed zero-detection output is valid"
        bad = root/"duplicate.fasta"
        bad.write_text(">P000001\nACGT\n>P000001\nACGT\n")
        try:
            batch.read_fasta(bad)
        except ValueError:
            pass
        else:
            raise AssertionError("duplicate output identifiers must fail")
        metrics = root/"metrics.tsv"
        batch.run_stage([sys.executable, "-c", "print('stage output')"], root, metrics, "infrastructure", "success")
        assert (root/"success.log").read_text() == "stage output\n"
        try:
            batch.run_stage([sys.executable, "-c", "raise SystemExit(7)"], root, metrics, "infrastructure", "failure")
        except RuntimeError as error:
            assert "exited 7" in str(error)
        else:
            raise AssertionError("stage failure must propagate")
        assert len(metrics.read_text().splitlines()) == 3, "both successful and failed stages need measurements"
    print("Batch input selection, empty/invalid output handling and per-stage failure accounting passed.")


if __name__ == "__main__":
    main()
