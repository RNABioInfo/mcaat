#!/usr/bin/env python3
"""Build a small graph, detect starts, then test all three local extension models.

No reference is supplied to prediction. All temporary reads, graphs and outputs
are deleted on exit. The expected molecules below are used only in assertions.
"""
import argparse
import csv
import hashlib
import gzip
import json
import random
import subprocess
import tempfile
from pathlib import Path

from check_source_crossing import canonical, fasta, graph_commands, reverse_complement, table


def check_flow(binary):
    """Check resolved flow, excluded ambiguous edges, and count inflation."""
    rng = random.Random(2026111704)
    dna = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    with tempfile.TemporaryDirectory(prefix="mcaat-flow-test-") as tmp:
        root = Path(tmp)
        arrays = root / "arrays"
        arrays.mkdir()
        cases, records = {}, []
        for name in ("linear", "ambiguous", "supported", "merge"):
            repeat, spacer, left, right = dna(32), dna(35), dna(1000), dna(1000)
            seq = left + spacer + right
            cases[name] = dict(spacer=spacer, repeat=repeat, sequence=seq)
            records.extend((f"{name}_{i}", seq) for i in range(20))
            p = len(left) + len(spacer) + 250
            if name == "ambiguous":
                alternative = seq[:p] + next(b for b in "ACGT" if b != seq[p]) + dna(1000)
                records.extend((f"alternate_{i}", alternative) for i in range(20))
            if name == "supported":
                records.append(("weak_alternative", seq[p - 23:p] + next(b for b in "ACGT" if b != seq[p]) + dna(1000)))
            if name == "merge":
                records.extend((f"high_copy_{i}", seq[p:p + 100]) for i in range(100))
        (arrays / "CRISPR_Arrays_1.txt").write_text("".join(
            f">Array_{name}\n{c['repeat']}\n{'-' * 32}\t{c['spacer']}\n" for name, c in cases.items()))
        reads = root / "construction.fasta"
        reads.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))
        (root / "graph_build/graph").mkdir(parents=True)
        def run(command):
            result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            assert result.returncode == 0, result.stdout[-10000:]
        for command in graph_commands(root, binary, reads).values():
            run(command)
        graph = root / "graph_build/graph/graph"
        hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in graph.parent.glob("graph.sdbg*")}
        base = [str(binary), "phage-curate", "--arrays", str(arrays), "--graph", str(graph),
                "--threads", "2", "--max-flank-length", "750", "--correct-spacer-boundaries", "false",
                "--min-identity", "1", "--min-query-coverage", "1", "--filter-extension-artifacts", "false"]
        run(base + ["--output", str(root / "seeds")])
        run(base + ["--extension-seeds", str(root / "seeds"), "--coverage-method", "flow",
                    "--output", str(root / "result")])
        result = root / "result"
        targets = fasta(result / "protospacers.fasta")
        rows = table(result / "extensions.tsv")
        assert len(targets) == len(cases) == len(rows), rows
        names = {canonical(c["spacer"]): name for name, c in cases.items()}
        for row in rows:
            seq = row["sequence"]
            begin, end = int(row["seed_begin"]), int(row["seed_end"])
            assert canonical(seq[begin:end]) == targets[row["target_id"]]
            name = names[canonical(seq[begin:end])]
            assert row["accepted"] == "1", (name, row)
            assert seq in cases[name]["sequence"] or reverse_complement(seq) in cases[name]["sequence"], (name, row)
            flanks = sorted([int(row["left_bp"]), int(row["right_bp"])])
            if name == "ambiguous":
                assert flanks == [250, 750], row
                assert "unresolved_flow_branch" in (row["left_stop"], row["right_stop"]), row
            else:
                assert flanks == [750, 750], (name, row)
            assert row["search_limited"] == "0" and row["ambiguous"] == "0", row
            for side in ("left", "right"):
                initial = float(row[f"{side}_seed_flow"])
                assert 0 <= float(row[f"{side}_remaining_flow"]) <= initial
                last = initial
                for entry in row[f"{side}_flow_decisions"].split(";"):
                    if not entry: continue
                    position, node, incoming, outgoing, chosen, remaining = entry.split(":")
                    assert float(incoming) <= last + 1e-5
                    if int(chosen) >= 0:
                        last = float(remaining)
                        assert last > 0
                    else:
                        assert int(position) == int(row[f"{side}_bp"]), row
        assert all(hashlib.sha256(p.read_bytes()).hexdigest() == sha for p, sha in hashes.items())
        print("PASS: four flow graphs; supported branch, stop BEFORE unresolved split, no flow inflation, both 750-bp flanks, preserved targets and unchanged graph.")



def check_lookahead(binary):
    """Keep a low-count bridge, preserve alternative routes and choose flank prefixes."""
    rng = random.Random(2026111705)
    dna = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    with tempfile.TemporaryDirectory(prefix="mcaat-lookahead-test-") as tmp:
        root = Path(tmp)
        arrays = root / "arrays"
        arrays.mkdir()
        cases, records = {}, []
        for name in ("linear", "bridge", "depth_change"):
            repeat, spacer, left, right = dna(32), dna(35), dna(1050), dna(1050)
            seq = left + spacer + right
            cases[name] = dict(spacer=spacer, repeat=repeat, sequence=seq, start=len(left))
            depth = 1 if name == "bridge" else 20
            records.extend((f"{name}_{i}", seq) for i in range(depth))
            if name == "bridge":
                p = len(left) + len(spacer) + 350
                # Coverage 20 before the fork, 1 across a 40-base connection,
                # and 20 downstream. The competing arm has depth 20 throughout.
                records.extend((f"bridge_prefix_{i}", seq[:p]) for i in range(19))
                records.extend((f"bridge_tail_{i}", seq[p + 40:]) for i in range(19))
                alternate = seq[p - 23:p] + next(b for b in "ACGT" if b != seq[p]) + dna(1100)
                records.extend((f"bridge_alternate_{i}", alternate) for i in range(20))
            if name == "depth_change":
                # An incompatible high-depth left flank must not spoil the
                # stable right flank. Right changes sharply after 550 bases.
                records.extend((f"deep_left_{i}", left) for i in range(60))
                p = len(left) + len(spacer) + 550
                records.extend((f"deep_tail_{i}", seq[p:]) for i in range(180))
        (arrays / "CRISPR_Arrays_1.txt").write_text("".join(
            f">Array_{name}\n{c['repeat']}\n{'-' * 32}\t{c['spacer']}\n" for name, c in cases.items()))
        reads = root / "construction.fasta"
        reads.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))
        (root / "graph_build/graph").mkdir(parents=True)
        def run(command):
            result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            assert result.returncode == 0, result.stdout[-10000:]
        for command in graph_commands(root, binary, reads).values():
            run(command)
        graph = root / "graph_build/graph/graph"
        hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in graph.parent.glob("graph.sdbg*")}
        base = [str(binary), "phage-curate", "--arrays", str(arrays), "--graph", str(graph),
                "--threads", "2", "--max-flank-length", "750", "--correct-spacer-boundaries", "false",
                "--min-identity", "1", "--min-query-coverage", "1", "--filter-extension-artifacts", "false"]
        run(base + ["--output", str(root / "seeds")])
        run(base + ["--extension-seeds", str(root / "seeds"), "--coverage-method", "lookahead",
                    "--output", str(root / "result")])
        result = root / "result"
        targets = fasta(result / "protospacers.fasta")
        evidence = table(result / "extensions.tsv")
        with gzip.open(result / "candidate_paths.tsv.gz", "rt") as stream:
            audit = list(csv.DictReader(stream, delimiter="\t"))
        assert len(targets) == len(cases) == len(evidence)
        names = {canonical(c["spacer"]): name for name, c in cases.items()}
        correct_bridge = False
        for row in evidence + audit:
            seq = row["sequence"]
            b, e = int(row["seed_begin"]), int(row["seed_end"])
            assert canonical(seq[b:e]) == targets[row["target_id"]]
            assert b <= 750 and len(seq) - e <= 750 and len(seq) >= 500
            name = names[canonical(seq[b:e])]
            if name == "bridge" and "side" in row:
                if max(b, len(seq) - e) == 750 and (seq in cases[name]["sequence"] or reverse_complement(seq) in cases[name]["sequence"]):
                    # Both true directions must exist, including the fork side.
                    if seq in cases[name]["sequence"][cases[name]["start"]:] or reverse_complement(seq) in cases[name]["sequence"][cases[name]["start"]:]:
                        correct_bridge = True
        assert correct_bridge, "low-count bridge lost its true full-length continuation"
        for row in evidence:
            name = names[targets[row["target_id"]]]
            assert row["accepted"] == row["reported"] == "1"
            if name == "linear":
                assert int(row["left_bp"]) == int(row["right_bp"]) == 750, row
            if name == "depth_change":
                assert min(int(row["left_bp"]), int(row["right_bp"])) == 0, row
        assert len(fasta(result / "viral_contigs_lookahead.fasta")) <= 50
        assert all(hashlib.sha256(p.read_bytes()).hexdigest() == sha for p, sha in hashes.items())
        print("PASS: lookahead retains the correct 750-base low-depth bridge route; independent sides, targets, bounded paths and graph verified.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("--flow", action="store_true", help="Run the flow-specific controlled graphs")
    parser.add_argument("--lookahead", action="store_true")
    args = parser.parse_args()
    binary = args.binary.resolve()
    if args.lookahead:
        check_lookahead(binary)
        return
    if args.flow:
        check_flow(binary)
        return
    rng = random.Random(2026111703)
    dna = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    with tempfile.TemporaryDirectory(prefix="mcaat-seeded-test-") as tmp:
        root = Path(tmp)
        arrays = root / "arrays"
        arrays.mkdir()
        cases = {}
        for name in ("two_flanks", "bubble_and_tip", "bounded_probe", "weak_link", "source_passage", "coverage_trap", "rare_early", "rare_mature"):
            repeat, spacer, left, right = dna(32), dna(35), dna(1050), dna(1050)
            cases[name] = dict(repeat=repeat, spacer=spacer, left=left, right=right)
        (arrays / "CRISPR_Arrays_1.txt").write_text("".join(
            f">Array_{name}\n{c['repeat']}\n{'-' * 32}\t{c['spacer']}\n"
            for name, c in cases.items()))
        records = []
        for name, c in cases.items():
            molecule = c["left"] + c["spacer"] + (c["repeat"] if name == "source_passage" else "") + c["right"]
            c["molecule"] = molecule
            # Long full molecules give constant interior edge multiplicity.
            records.extend((f"{name}_main_{i}", molecule) for i in range(1 if name == "coverage_trap" else 20))
            records.extend((f"{name}_array_{i}", c["repeat"] + c["spacer"] + c["repeat"])
                           for i in range(80))
        c = cases["bubble_and_tip"]
        molecule = c["molecule"]
        start = len(c["left"]) + len(c["spacer"])
        # One low-depth substitution makes a reconverging simple bubble.
        p = start + 300
        alternative = molecule[:p] + next(b for b in "ACGT" if b != molecule[p]) + molecule[p + 1:]
        records.append(("low_depth_bubble", alternative))
        # A genuine ten-edge dead end is a short tip, not a missing probe tail.
        p = start + 400
        records.append(("short_tip", molecule[p - 23:p] + next(b for b in "ACGT" if b != molecule[p]) + dna(9)))
        # Both branches near the 750-base bound continue beyond the window.
        # Neither may be declared a short tip because its probe was truncated.
        c = cases["bounded_probe"]
        p = len(c["left"]) + len(c["spacer"]) + 740
        alternate_tail = next(b for b in "ACGT" if b != c["molecule"][p]) + dna(400)
        records.extend((f"long_alternative_{i}", c["molecule"][:p] + alternate_tail) for i in range(20))
        # A low-depth 70-edge arm is long enough not to be a short tip.
        # The high-depth arm ends at a second fork, so both probes are complete.
        c = cases["weak_link"]
        p = len(c["left"]) + len(c["spacer"]) + 350
        records.append(("weak_arm", c["molecule"][p - 23:p]
                        + next(b for b in "ACGT" if b != c["molecule"][p]) + dna(69)))
        p += 150
        records.append(("weak_case_bubble", c["molecule"][:p]
                        + next(b for b in "ACGT" if b != c["molecule"][p]) + c["molecule"][p + 1:]))
        # A source-exit path at count 1 forks into a well-supported continuation
        # at count 21 and a long, perfectly flat count-1 alternative. Calibrate
        # beyond the source exit and use branch support: flatness is not truth.
        c = cases["coverage_trap"]
        p = len(c["left"]) + len(c["spacer"]) + 60
        records.extend((f"supported_tail_{i}", c["molecule"][p - 23:]) for i in range(20))
        records.extend((f"supported_left_{i}", c["left"][:-60]) for i in range(20))
        records.append(("flat_wrong_tail", c["molecule"][p - 23:p]
                        + next(b for b in "ACGT" if b != c["molecule"][p]) + dna(1400)))
        # Rare true arms must survive a much more abundant neighboring unitig:
        # once before calibration, once after a matching local baseline exists.
        for name, offset in (("rare_early", 60), ("rare_mature", 350)):
            c = cases[name]
            p = len(c["left"]) + len(c["spacer"]) + offset
            high = c["molecule"][p - 23:p] + next(b for b in "ACGT" if b != c["molecule"][p]) + dna(89)
            records.extend((f"{name}_high_{i}", high) for i in range(1000))
            p = len(c["left"]) + len(c["spacer"]) + 500
            records.append((f"{name}_bubble", c["molecule"][:p]
                            + next(b for b in "ACGT" if b != c["molecule"][p]) + c["molecule"][p + 1:]))
        reads = root / "construction.fasta"
        reads.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))
        (root / "graph_build/graph").mkdir(parents=True)
        log = root / "test.log"

        def run(command):
            with log.open("a") as stream:
                result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
            if result.returncode:
                raise AssertionError(log.read_text()[-10000:])

        for command in graph_commands(root, binary, reads).values():
            run(command)
        graph = root / "graph_build/graph/graph"
        hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in graph.parent.glob("graph.sdbg*")}
        base = [str(binary), "phage-curate", "--arrays", str(arrays), "--graph", str(graph),
                "--threads", "2", "--max-flank-length", "750", "--correct-spacer-boundaries", "false",
                "--min-identity", "1", "--min-query-coverage", "1", "--filter-extension-artifacts", "false"]
        run(base + ["--output", str(root / "seeds")])
        seed_summary = json.loads((root / "seeds/run_summary.json").read_text())
        assert seed_summary["detected_queries"] == len(cases), seed_summary
        run(base + ["--extension-seeds", str(root / "seeds"), "--coverage-method", "all",
                    "--output", str(root / "result")])
        result = root / "result"
        summary = json.loads((result / "summary.json").read_text())
        assert summary["status"] == "complete" and summary["input_starts"] == len(cases), summary
        targets = fasta(result / "protospacers.fasta")
        rows = table(result / "extensions.tsv")
        assert len(targets) == len(cases) and len(rows) == 3 * len(cases)
        names = {canonical(c["spacer"]): name for name, c in cases.items()}
        for row in rows:
            assert row["accepted"] == "1", row
            seq = row["sequence"]
            begin, end = int(row["seed_begin"]), int(row["seed_end"])
            assert canonical(seq[begin:end]) == targets[row["target_id"]]
            name = names[canonical(seq[begin:end])]
            case = cases[name]
            assert begin == 750 and len(seq) - end == 750, (name, row)
            if name != "bounded_probe":
                assert seq in case["molecule"] or reverse_complement(seq) in case["molecule"], (name, seq)
            else:
                assert int(row["tips_removed"]) == 0, "length-bounded probe falsely called a tip"
            if name == "bubble_and_tip":
                assert int(row["bubbles_removed"]) >= 1, row
                assert int(row["tips_removed"]) >= 1, row
            if name == "weak_link":
                assert int(row["weak_links_removed"]) >= 1, row
            if name == "source_passage":
                assert max(float(row["left_baseline"]), float(row["right_baseline"])) < 30, row
            for side in ("left", "right"):
                counts = [float(x) for x in row[f"{side}_counts"].split(",") if x]
                positions = [int(x) for x in row[f"{side}_scored_positions"].split(",") if x]
                assert len(counts) == 750 and all(1 <= p <= 750 for p in positions)
                assert positions == sorted(set(positions))
        for method in ("poisson", "negbin", "huber"):
            representatives = fasta(result / f"viral_contigs_{method}.fasta")
            assert len(representatives) == len(cases)
            assert {r["viral_contig_id"] for r in rows if r["method"] == method} == set(representatives)
        assert all(hashlib.sha256(p.read_bytes()).hexdigest() == sha for p, sha in hashes.items()), "graph modified"
        print("PASS: all three models; unchanged targets/graph; 750 bases on each side; bubble, tip, weak-link, source-mask, probe-bound, flat-wrong-branch and rare-arm protection tests.")


if __name__ == "__main__":
    main()
