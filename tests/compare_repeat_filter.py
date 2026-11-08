#!/usr/bin/env python3
"""Compare complete filter results from two benchmark binaries, using seeded fixtures."""
import argparse
import json
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("before", type=Path)
    parser.add_argument("after", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--seeds", type=int, default=30)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    cases = []
    for mode in ("identical", "mixed", "variants", "ambiguous"):
        for support in ("0.60", "0.85", "1.0"):
            for seed in range(1, args.seeds + 1):
                stem = f"{mode}-{support}-{seed}"
                outputs = []
                for label, binary in (("before", args.before), ("after", args.after)):
                    path = args.output / f"{stem}-{label}.txt"
                    command = [str(binary.resolve()), "2", "6", str(path), str(seed), mode, support]
                    result = subprocess.run(command, check=True, capture_output=True, text=True, timeout=60)
                    outputs.append((path, result.stdout.strip()))
                if outputs[0][0].read_bytes() != outputs[1][0].read_bytes():
                    raise AssertionError(f"Different results for {stem}; inspect {args.output}")
                cases.append({"case": stem, "before": outputs[0][1], "after": outputs[1][1]})
                # Keep mismatch files, but avoid retaining hundreds of identical dumps.
                for path, _ in outputs:
                    path.unlink()
            print(f"Matched: {mode}, support={support}, {args.seeds} seeds", flush=True)
    (args.output / "results.json").write_text(json.dumps(cases, indent=2) + "\n")
    print(f"All {len(cases)} cases matched exactly (systems, paths, sequences, sources and original statistics).")


if __name__ == "__main__":
    main()
