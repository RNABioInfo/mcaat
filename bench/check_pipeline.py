"""Compare complete single-thread MCAAT array runs on an existing SDBG."""
import argparse
import hashlib
import json
import pathlib
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("executable", type=pathlib.Path)
parser.add_argument("graph", type=pathlib.Path)
parser.add_argument("output", type=pathlib.Path)
parser.add_argument("--compare-only", action="store_true")
args = parser.parse_args()
exe = args.executable.resolve()
args.output.mkdir(parents=True, exist_ok=True)
help_result = subprocess.run([str(exe), "--help"], capture_output=True, text=True, check=True)
assert "--cycle-precheck" in help_result.stdout
summaries = {}
for mode in ("bfs", "unitig"):
    dest = args.output / mode
    if dest.exists() and not args.compare_only:
        raise SystemExit(f"Use a fresh output directory: {dest} already exists")
    if not args.compare_only:
        with (args.output / f"{mode}.log").open("w") as log:
            subprocess.run([
                str(exe), "--graph", str(args.graph), "--output-folder", str(dest),
                "--threads", "1", "--ram", "2G", "--autoclean", "false",
                "--cycle-precheck", mode,
            ], stdout=log, stderr=subprocess.STDOUT, check=True)
    parameters = json.loads((dest / "parameters.json").read_text())
    assert parameters["cycle_finder"]["cycle_precheck"] == mode
    cycles = (dest / "cycles.txt").read_text().splitlines()
    # Exclude only the known wall-clock header; preserve all sequence, ID,
    # ordering and count fields. Different runs necessarily have different times.
    arrays = {p.name: "\n".join(line for line in p.read_text().splitlines()
              if not line.startswith("# Generated :")).encode()
              for p in dest.glob("CRISPR_Arrays_*.txt")}
    summaries[mode] = {
        "cycles": len(cycles),
        "canonical_cycles_sha256": hashlib.sha256("\n".join(sorted(cycles)).encode()).hexdigest(),
        "array_files": {name: hashlib.sha256(data).hexdigest() for name, data in arrays.items()},
    }
    if not cycles or not arrays:
        raise SystemExit(f"{mode}: choose a fixture producing cycles and arrays")
assert summaries["bfs"] == summaries["unitig"], "Cycle or array output changed"
result = {"graph": str(args.graph), "matched": True,
          "normalization": "Sort cycle lines; omit only array '# Generated :' timestamps",
          "outputs": summaries}
(args.output / "comparison.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
