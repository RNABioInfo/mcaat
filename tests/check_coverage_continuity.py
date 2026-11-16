#!/usr/bin/env python3
"""Verify final coverage validation withholds an unstable selected route.

Build a graph from counted 24-mers. Two 750-base continuations share 200 bases;
the dominant branch jumps from multiplicity 10 to 80, while its competitor
stays at 10. The final filter must withhold the dominant route without choosing
the weaker alternative just because it appears smoother. Prediction receives the graph and target only. Check the saved
sequence, coordinates and counted evidence against fixture truth afterward.
Usage: python3 tests/check_coverage_continuity.py BINARY NEW_OUTPUT
"""
import argparse
import json
import random
import subprocess
from pathlib import Path
from check_source_crossing import canonical, digest, fasta, table, write_fasta, graph_commands


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    binary, output = args.binary.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    rng = random.Random(2026111608)
    dna = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    repeat, target, left, shared = dna(32), dna(35), dna(80), dna(200)
    good, bad = dna(550), dna(550)
    while good[0] == bad[0]:
        bad = dna(550)
    molecules = {"stable": left + target + shared + good,
                 "jump": left + target + shared + bad}
    counts = {}
    for name, molecule in molecules.items():
        for start in range(len(molecule) - 23):
            word = canonical(molecule[start:start + 24])
            value = 80 if name == "jump" and start + 24 > len(left + target + shared) else 10
            counts[word] = max(counts.get(word, 0), value)
    with (output / "counts.fasta").open("w") as stream:
        for index, (word, copies) in enumerate(sorted(counts.items())):
            for copy in range(copies):
                stream.write(f">edge_{index}_{copy}\n{word}\n")
    arrays = output / "arrays"
    arrays.mkdir()
    (arrays / "CRISPR_Arrays_1.txt").write_text(f">Array_1\n{repeat}\n{'-' * len(repeat)}\t{target}\n")
    write_fasta(output / "supplied.fasta", {"known_target": target})
    base = [str(binary), "phage-curate", "--arrays", str(arrays), "--graph", str(output / "graph_build/graph/graph"),
            "--protospacers", str(output / "supplied.fasta"), "--threads", "2", "--max-flank-length", "750",
            "--source-array-crossing", "false", "--correct-spacer-boundaries", "false", "--save-search", "true"]
    commands = {**graph_commands(output, binary, output / "counts.fasta"), "disabled": base + ["--max-coverage-fold", "0", "--output", str(output / "disabled")],
                "guarded": base + ["--max-coverage-fold", "4", "--output", str(output / "guarded")]}
    (output / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
    (output / "truth.json").write_text(json.dumps(dict(molecules=molecules, target=target, counts=counts), indent=2) + "\n")
    graph_hashes = {}
    (output / "graph_build/graph").mkdir(parents=True, exist_ok=True)
    for name, command in commands.items():
        with (output / (name + ".log")).open("w") as stream:
            subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, check=True)
        if name == "build_graph":
            graph_hashes = {str(p): digest(p) for p in (output / "graph_build/graph").glob("graph.sdbg*")}
    disabled = output / "disabled"
    saved = fasta(disabled / "ContiguousPaths.fasta")
    expected = canonical(target + shared + bad)
    assert len(saved) == 1 and next(iter(saved.values())) == expected, saved
    for link in table(disabled / "associations.tsv"):
        sequence = saved[link["extension_id"]]
        assert canonical(sequence[int(link["seed_begin"]):int(link["seed_end"])]) == canonical(target)
    guarded = output / "guarded"
    assert fasta(guarded / "SuppliedProtospacers.fasta") == {"known_target": target}
    assert set(fasta(guarded / "Protospacers.fasta").values()) == {canonical(target)}
    assert not fasta(guarded / "ContiguousPaths.fasta"), "unstable route was replaced by a permissive fallback"
    withheld = table(guarded / "withheld_extensions.tsv")
    assert len(withheld) == 1 and canonical(withheld[0]["sequence"]) == expected
    assert float(withheld[0]["block_max"]) / float(withheld[0]["block_min"]) == 8
    assert len(withheld[0]["counts"].split(",")) >= 700
    assert all(not row["extension_id"] for row in table(guarded / "associations.tsv"))
    assert (disabled / "target_context.tsv").read_bytes() == (guarded / "target_context.tsv").read_bytes()
    assert graph_hashes and all(digest(Path(p)) == value for p, value in graph_hashes.items())
    result = dict(passed=True, supplied_target_retained=True, disabled_selects_coverage_jump=True,
                  guarded_withholds_unstable_winner=True, weaker_alternative_not_promoted=True, withheld_flank_length=750,
                  graph_unchanged=True, binary_sha256=digest(binary))
    (output / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
