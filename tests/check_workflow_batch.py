#!/usr/bin/env python3
"""Infrastructure checks only; no biological dataset, downloads, or graph builds."""
import importlib.util
import json
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
        predictions = root / "predictions"
        predictions.mkdir()
        summary = dict(status="complete", implementation="protospacer-workflow-v9",
                       protospacer_output_policy="accepted_extensions_only", protospacer_regions=2,
                       candidate_protospacer_regions=3, candidate_only_regions=1,
                       unique_extensions=1, extended_queries=2)
        (predictions / "run_summary.json").write_text(json.dumps(summary))
        (predictions / "CandidateProtospacers.fasta").write_text(">P1\nACG\n>P2\nGTA\n>P3\nCCC\n")
        (predictions / "Protospacers.fasta").write_text(">P1\nACG\n>P2\nGTA\n")
        (predictions / "ContiguousPaths.fasta").write_text(">E1\nAACGTAC\n")
        header = "query_id\ttarget_id\textension_id\tseed_begin\tseed_end\n"
        links = predictions / "associations.tsv"
        links.write_text(header + "Q1\tP1\tE1\t1\t4\nQ2\tP2\tE1\t3\t6\n")
        assert batch.validate_predictions(predictions) == summary, "two targets may share one saved extension"
        links.write_text(header + "Q1\tP1\tE1\t1\t4\n")
        try:
            batch.validate_predictions(predictions)
        except ValueError:
            pass
        else:
            raise AssertionError("final target without an accepted extension must fail")
        summary.update(protospacer_regions=0, unique_extensions=0, extended_queries=0, candidate_only_regions=3)
        (predictions / "run_summary.json").write_text(json.dumps(summary))
        (predictions / "Protospacers.fasta").write_text("")
        (predictions / "ContiguousPaths.fasta").write_text("")
        links.write_text(header)
        assert batch.validate_predictions(predictions) == summary, "candidate-only output is valid with zero final calls"
    print("Batch input selection, empty/invalid output handling and per-stage failure accounting passed.")


if __name__ == "__main__":
    main()
