#!/usr/bin/env python3
"""Regression fixture: an extension must not follow a CRISPR array into the
phage of a neighbouring spacer.

crossing      Array_1 = R1 S1 R1 S2 R1 S3 R1 (reported). S1 has only 120 bp of
              its own phage context per side; S2's phage is well covered. The
              path S1 -> R1 -> S2 -> S2's phage is graph-valid but chimeric.
second_array  S3 also sits, unreported, in a host array of a different repeat
              R2 next to S5, whose phage is well covered. R2 is reported through
              Array_2 = R2 S6 R2 S7 R2. The path S3 -> R2 -> S5 -> S5's phage is
              not source context for S3, so only a repeat barrier stops it.

Both server stages run: discovery, then seeded evidence at branch support 0.998.
A saved extension containing an 18-base word of R1 or R2 is an array chimera.

Usage: python3 tests/check_array_chimera.py BINARY OUTPUT [--expect-chimera]
  default           the binary's defaults must save no array chimera, S2 must still
                    extend inside its own phage, and turning the repeat barrier
                    off must bring the second-array chimera back.
  --expect-chimera  for an older binary: require at least one chimera (shows the
                    fixture discriminates); prints the same report.
Builds its graph with the binary's bundled tools; nothing is downloaded.
"""
import argparse
import json
import random
import subprocess
from pathlib import Path

from check_source_crossing import canonical, fasta, graph_commands, reverse_complement, table

WORD = 18


def words(sequence):
    return {o[i:i + WORD] for o in (sequence, reverse_complement(sequence)) for i in range(len(o) - WORD + 1)}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("binary", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--expect-chimera", action="store_true")
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    rng = random.Random(2026111901)
    dna = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    R1, R2 = dna(32), dna(32)
    S = {name: dna(35) for name in ("S1", "S2", "S3", "S5", "S6", "S7")}
    phage = {"S1": dna(120) + S["S1"] + dna(120),   # weak own context
             "S2": dna(1050) + S["S2"] + dna(1050),
             "S3": dna(120) + S["S3"] + dna(120),   # weak own context
             "S5": dna(1050) + S["S5"] + dna(1050)}
    molecules = {
        "array_1": dna(300) + R1 + S["S1"] + R1 + S["S2"] + R1 + S["S3"] + R1 + dna(300),
        "unreported_host_array": dna(300) + R2 + S["S3"] + R2 + S["S5"] + R2 + dna(300),
        "array_2": dna(300) + R2 + S["S6"] + R2 + S["S7"] + R2 + dna(300),
        **{f"phage_{k}": v for k, v in phage.items()}}
    arrays = root / "arrays"
    arrays.mkdir()
    gap = "-" * 32
    (arrays / "CRISPR_Arrays_1.txt").write_text(
        f">Array_1\n{R1}\n" + "".join(f"{gap}\t{S[s]}\n" for s in ("S1", "S2", "S3")) +
        f">Array_2\n{R2}\n" + "".join(f"{gap}\t{S[s]}\n" for s in ("S6", "S7")))
    reads = root / "reads.fasta"
    reads.write_text("".join(f">{name}_{i}\n{seq}\n" for name, seq in molecules.items() for i in range(20)))
    (root / "graph_build/graph").mkdir(parents=True)
    log = root / "test.log"

    def run(command):
        with log.open("a") as stream:
            stream.write(" ".join(map(str, command)) + "\n")
            stream.flush()
            result = subprocess.run([str(c) for c in command], stdout=stream, stderr=subprocess.STDOUT)
        if result.returncode:
            raise AssertionError(log.read_text()[-6000:])

    for command in graph_commands(root, args.binary, reads).values():
        run(command)
    graph = root / "graph_build/graph/graph"
    common = [args.binary, "phage-curate", "--arrays", arrays, "--graph", graph, "--threads", "2",
              "--correct-spacer-boundaries", "false", "--min-identity", "1", "--min-query-coverage", "1"]
    help_text = subprocess.run([str(args.binary), "phage-curate", "--help"], text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout
    configurations = {"defaults": []}
    if "--crispr-repeat-barrier" in help_text:
        configurations["no_barrier"] = ["--crispr-repeat-barrier", "false"]
    repeat_words = words(R1) | words(R2)
    spacer_of = {canonical(seq): name for name, seq in S.items()}
    report = {}
    for label, extra in configurations.items():
        seeds, support = root / f"{label}_seeds", root / f"{label}_support"
        run(common + extra + ["--output", seeds])
        run(common + extra + ["--extension-seeds", seeds, "--output", support, "--coverage-method", "evidence",
                              "--min-contig-length", "500", "--max-flank-length", "750", "--beam-width", "64",
                              "--max-context-paths", "256", "--min-branch-support", "0.998",
                              "--max-reported-contigs", "0"])
        saved = []
        contigs = fasta(seeds / "ContiguousPaths.fasta")
        targets = fasta(seeds / "CandidateProtospacers.fasta")
        for row in table(seeds / "associations.tsv"):
            if row["extension_id"]:
                saved.append(("discovery", targets[row["target_id"]], contigs[row["extension_id"]]))
        targets = fasta(support / "protospacers.fasta")
        for row in table(support / "extensions.tsv"):
            if row["reported"] == "1" and row["relation_to_representative"] == "representative":
                saved.append(("support", targets[row["target_id"]], row["sequence"]))
        outcome = []
        for stage, target, sequence in saved:
            spacer = spacer_of.get(canonical(target), "other")
            own = phage.get(spacer, "")
            outcome.append(dict(stage=stage, spacer=spacer, length=len(sequence),
                                array_chimera=any(sequence[i:i + WORD] in repeat_words
                                                  for i in range(len(sequence) - WORD + 1)),
                                inside_own_phage=bool(own) and (sequence in own or reverse_complement(sequence) in own)))
        report[label] = outcome
    summary = {label: dict(saved=len(o), array_chimeras=sum(x["array_chimera"] for x in o),
                           chimeric_spacers=sorted({f"{x['stage']}:{x['spacer']}" for x in o if x["array_chimera"]}),
                           extended_inside_own_phage=sorted({f"{x['stage']}:{x['spacer']}" for x in o if x["inside_own_phage"]}))
               for label, o in report.items()}
    (root / "report.json").write_text(json.dumps(dict(summary=summary, records=report), indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    defaults = summary["defaults"]
    if args.expect_chimera:
        assert defaults["array_chimeras"] > 0, "fixture produced no chimera with this binary"
        print("PASS (expected chimera): fixture reproduces the array chimera with this binary.")
        return
    assert defaults["array_chimeras"] == 0, ("array chimera saved with defaults", defaults["chimeric_spacers"])
    # S5 is only in the unreported host array, so it is never a query; S2 is
    # the reported neighbour whose own phage extension must survive.
    for stage in ("discovery", "support"):
        assert f"{stage}:S2" in defaults["extended_inside_own_phage"], (stage, "S2 true extension lost")
    for x in report["defaults"]:
        if x["spacer"] in ("S1", "S3"):
            assert x["inside_own_phage"], ("weak-context target extended outside its own phage", x)
    if "no_barrier" in summary:
        assert summary["no_barrier"]["array_chimeras"] > 0, \
            "fixture no longer shows the second-array chimera without the barrier"
    print("PASS: defaults save no array chimera; S2 still extends inside its own phage; "
          "without the repeat barrier the second-array chimera returns.")


if __name__ == "__main__":
    main()
