#!/usr/bin/env python3
"""Exercise source crossing on a reproducible tiny native succinct graph.

Usage: python3 tests/check_source_crossing.py BINARY NEW_OUTPUT [--prepare-only]
Preparation writes synthetic reads, arrays, supplied targets, truth and commands.
It starts no process. A prepared directory may subsequently be run without the
flag. Native construction and predictions use two threads; no reference tool is
invoked. All biological truth checks are performed after prediction.
"""
import argparse
import csv
import gzip
import hashlib
import json
import random
import subprocess
from collections import Counter
from pathlib import Path


K = 23


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def canonical(sequence):
    return min(sequence, reverse_complement(sequence))


def fasta(path):
    result = {}
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            assert name not in result, (path, "duplicate identifier", name)
            result[name] = ""
        else:
            result[name] += line.strip()
    return result


def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_fasta(path, records):
    path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in records.items()))


def make_fixture(output, binary):
    rng = random.Random(2026111601)
    dna = lambda length: "".join(rng.choice("ACGT") for _ in range(length))
    cases = {}
    for name, flank_length in (("single_pass", 500), ("reentry", 300), ("variant_reentry", 300)):
        repeat, spacer, leader, flank = dna(32), dna(35), dna(80), dna(flank_length)
        first_repeat = repeat
        if name == "variant_reentry":
            # The central substitution changes every full repeat k-mer. The
            # first passage is still recognized at 85% identity; later exact
            # repeat nodes are new to this walk and exercise re-entry rather
            # than the simpler source-node revisit guard.
            at = 15
            first_repeat = repeat[:at] + "ACGT"[("ACGT".index(repeat[at]) + 1) % 4] + repeat[at + 1:]
        molecule = leader + repeat + spacer + first_repeat + flank
        if name != "single_pass":
            molecule += repeat
        cases[name] = dict(repeat=repeat, spacer=spacer, molecule=molecule,
                           downstream=flank, source_sequence=repeat + spacer + repeat)
    # A source sequence occurs only in its source array: it is not discovery
    # evidence. Its known target can nevertheless be extended through one repeat.
    # The second component returns to that repeat after 300 external bases.
    arrays = output / "arrays"
    arrays.mkdir(exist_ok=True)
    (arrays / "CRISPR_Arrays_1.txt").write_text("".join(
        f">Array_{name}\n{case['repeat']}\n{'-' * len(case['repeat'])}\t{case['spacer']}\n"
        for name, case in cases.items()))
    short = cases["single_pass"]["spacer"][8:27]
    supplied = {"known_single_pass": cases["single_pass"]["spacer"],
                "known_reentry": cases["reentry"]["spacer"],
                "known_variant_reentry": cases["variant_reentry"]["spacer"],
                "short_19": short, "short_19_reverse": reverse_complement(short)}
    assert len(short) == 19
    write_fasta(output / "supplied.fasta", supplied)
    # Sliding reads guarantee all intended adjacencies. Extra source-only reads
    # create a deliberate source multiplicity spike; those source nodes must not
    # enter the exported downstream coverage trace.
    with (output / "construction.fasta").open("w") as stream:
        for name, case in cases.items():
            sequence = case["molecule"]
            for position in range(len(sequence) - 149):
                for copy in range(2):
                    stream.write(f">{name}_{position}_{copy}\n{sequence[position:position + 150]}\n")
            for copy in range(240):
                stream.write(f">{name}_source_boost_{copy}\n{case['source_sequence']}\n")
    graph = output / "graph_build/graph/graph"
    base = [str(binary), "phage-curate", "--arrays", str(arrays), "--graph", str(graph),
            "--threads", "2", "--max-flank-length", "600", "--min-extension-flank", "200",
            "--min-identity", "1", "--min-query-coverage", "1", "--correct-spacer-boundaries", "false",
            "--beam-width", "64", "--save-search", "true"]
    commands = {
        "build": [str(binary), "--input-files", str(output / "construction.fasta"),
                  "--output-folder", str(output / "graph_build"), "--threads", "2", "--ram", "2G",
                  "--threshold-multiplicity", "1000000", "--autoclean", "false"],
        "discovery": base + ["--output", str(output / "discovery"), "--source-array-crossing", "true"],
        "supplied_strict": base + ["--protospacers", str(output / "supplied.fasta"),
                                   "--output", str(output / "supplied_strict"), "--source-array-crossing", "false"],
        "supplied_crossing": base + ["--protospacers", str(output / "supplied.fasta"),
                                     "--output", str(output / "supplied_crossing"), "--source-array-crossing", "true"],
        "supplied_exact_length": base + ["--protospacers", str(output / "supplied.fasta"),
                                         "--output", str(output / "supplied_exact_length"),
                                         "--source-array-crossing", "true",
                                         "--max-flank-length", "500", "--min-extension-flank", "500"],
    }
    fixture = {"k": K, "random_seed": 2026111601, "cases": cases, "supplied": supplied,
               "binary": str(binary), "binary_sha256": digest(binary),
               "reference_used_for_prediction": False}
    (output / "fixture.json").write_text(json.dumps(fixture, indent=2) + "\n")
    (output / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
    return fixture, commands


def source_words(case):
    sequence = case["source_sequence"]
    return {oriented[i:i + K] for oriented in (sequence, reverse_complement(sequence))
            for i in range(len(oriented) - K + 1)}


def canonical_fragment(link, query, sequences):
    sequence = sequences[link["extension_id"]]
    begin, end = int(link["seed_begin"]), int(link["seed_end"])
    low, high = int(link["supported_fragment_begin"]), int(link["supported_fragment_end"])
    if sequence[begin:end] != query:
        assert reverse_complement(sequence[begin:end]) == query, "target coordinates changed"
        length = len(sequence)
        sequence = reverse_complement(sequence)
        begin, end, low, high = length - end, length - begin, length - high, length - low
    return sequence[low:high], begin - low, end - low


def check_exact_supplied_associations(run, supplied):
    """Every exported target association covers its complete unchanged input."""
    targets = fasta(run / "Protospacers.fasta")
    sequences = fasta(run / "ContiguousPaths.fasta")
    seen = set()
    for link in table(run / "associations.tsv"):
        input_id = link["source_id"].rsplit(":", 1)[0]
        original = supplied[input_id]
        assert int(link["query_begin"]) == 0
        assert int(link["query_end"]) == len(original)
        assert float(link["identity"]) == float(link["query_coverage"]) == 1.0
        assert canonical(targets[link["target_id"]]) == canonical(original)
        if link["extension_id"]:
            sequence = sequences[link["extension_id"]]
            begin, end = int(link["seed_begin"]), int(link["seed_end"])
            observed = sequence[begin:end]
            if link["target_strand"] == "-":
                observed = reverse_complement(observed)
            assert observed == original, (input_id, "reported target coordinates disagree")
        seen.add(input_id)
    assert seen == set(supplied), ("missing supplied target associations", set(supplied) - seen)


def check_outputs(output, fixture):
    discovery = output / "discovery"
    summary = json.loads((discovery / "run_summary.json").read_text())
    assert summary["status"] == "complete"
    assert summary["detected_queries"] == summary["unique_extensions"] == 0, (
        "source-only paths became new targets", summary)
    assert not fasta(discovery / "Protospacers.fasta")
    assert not fasta(discovery / "ContiguousPaths.fasta")

    supplied = fixture["supplied"]
    for name in ("supplied_strict", "supplied_crossing", "supplied_exact_length"):
        run = output / name
        assert fasta(run / "SuppliedProtospacers.fasta") == supplied, "supplied records changed"
        mapping = table(run / "supplied_targets.tsv")
        assert {row["input_id"]: row["sequence"] for row in mapping} == supplied
        queries = {row["query_id"]: row["spacer"] for row in table(run / "queries.tsv")}
        assert all(queries[row["query_id"]] == canonical(row["sequence"]) for row in mapping)
        check_exact_supplied_associations(run, supplied)
    assert not fasta(output / "supplied_strict/ContiguousPaths.fasta"), (
        "strict extension unexpectedly traversed a source repeat")

    run = output / "supplied_crossing"
    mapping = {row["input_id"]: row for row in table(run / "supplied_targets.tsv")}
    queries = {row["query_id"]: row["spacer"] for row in table(run / "queries.tsv")}
    sequences = fasta(run / "ContiguousPaths.fasta")
    links = table(run / "associations.tsv")
    candidates = {row["query_id"]: row for row in table(run / "candidates.tsv") if row["status"] == "extended"}
    evidence = table(run / "coverage_evidence.tsv")
    checked = {}
    for input_id, input_sequence in supplied.items():
        query_id = mapping[input_id]["query_id"]
        own_links = [row for row in links if row["query_id"] == query_id and row["extension_id"]
                     and row["source_id"].startswith(input_id + ":")]
        assert own_links, (input_id, "known source passage was not extended")
        case_name = input_id.removeprefix("known_") if input_id.startswith("known_") else "single_pass"
        case = fixture["cases"][case_name]
        words = source_words(case)
        candidate = candidates[query_id]
        assert candidate["identified_before_extension"] == "1"
        assert "1" in (candidate["left_crosses_source_array"], candidate["right_crosses_source_array"])
        for link in own_links:
            fragment, begin, end = canonical_fragment(link, queries[query_id], sequences)
            assert fragment in case["molecule"] or reverse_complement(fragment) in case["molecule"], (
                input_id, "output is not a contiguous intended graph route")
            assert fragment[begin:end] == canonical(input_sequence), (input_id, "seed was altered")
            assert end - begin == len(input_sequence)
            assert max(begin, len(fragment) - end) >= 200
            assert max(begin, len(fragment) - end) <= 600
            occurrences = Counter(fragment[i:i + K] for i in range(len(fragment) - K + 1))
            assert all(count <= 1 for word, count in occurrences.items() if word in words), (
                input_id, "source node revisited in saved prediction")
            if "reentry" in input_id:
                # The second repeat first becomes a source k-mer after K bases.
                # At most its first K-1 bases can remain before re-entry stops.
                assert max(begin, len(fragment) - end) <= 32 + 300 + K - 1
            recorded = [row for row in evidence if row["target_id"] == link["target_id"]
                        and row["extension_id"] == link["extension_id"]]
            assert recorded, (input_id, "extension has no downstream counts")
            padding = int(candidate.get("anchor_prefix_bp", "0"))
            for row in recorded:
                position = int(row["added_base_position"])
                if row["side"] == "left":
                    start = begin - position
                    word = fragment[start:start + K]
                else:
                    stop = end + position
                    word = fragment[stop - K:stop]
                assert len(word) == K
                assert word not in words, (input_id, "source multiplicity leaked into score", row)
            checked[input_id] = dict(query_id=query_id, extension_id=link["extension_id"],
                                    seed_length=end - begin, length=len(fragment),
                                    coverage_observations=len(recorded), anchor_prefix_bp=padding)
    assert checked["short_19"]["seed_length"] == 19
    assert checked["short_19"]["anchor_prefix_bp"] == K - 19
    assert mapping["short_19"]["query_id"] == mapping["short_19_reverse"]["query_id"]
    exact = output / "supplied_exact_length"
    exact_mapping = {row["input_id"]: row for row in table(exact / "supplied_targets.tsv")}
    exact_sequences = fasta(exact / "ContiguousPaths.fasta")
    exact_short_links = [row for row in table(exact / "associations.tsv")
                         if row["source_id"].startswith("short_19:") and row["extension_id"]]
    assert exact_short_links, "19-base target was lost when minimum and maximum flank both equal 500"
    for link in exact_short_links:
        fragment, begin, end = canonical_fragment(link, canonical(supplied["short_19"]), exact_sequences)
        assert int(link["left_bp"]) == begin == 500, ("padding not counted consistently in flank bounds", link)
        assert end - begin == 19
        assert link["query_id"] == exact_mapping["short_19"]["query_id"]
    # Confirm the intended rejection was actually exercised, not just avoided.
    stops = Counter()
    audit_records = {}
    for name in ("discovery", "supplied_strict", "supplied_crossing", "supplied_exact_length"):
        for path in (output / name / "search_paths").glob("*.tsv.gz"):
            path_ids = set()
            with gzip.open(path, "rt") as stream:
                for row in csv.DictReader(stream, delimiter="\t"):
                    assert row["path_id"] not in path_ids, (
                        "strict and recovery audit records reused an identifier", path, row["path_id"])
                    path_ids.add(row["path_id"])
                    if name == "supplied_crossing":
                        stops[row["stop"]] += 1
            audit_records[str(path.relative_to(output))] = len(path_ids)
    assert stops["source_array_reentry"] > 0, stops
    assert stops["source_node_revisit"] > 0, stops
    assert stops["anchor_revisit"] > 0, stops
    return dict(strict_discovery_rejects_source_only=True, supplied_records_retained=len(supplied),
                single_pass_extension=True, source_reentry_exercised=True,
                source_counts_excluded=True, short_target_coordinates_preserved=True,
                exact_full_supplied_associations=True, short_minimum_equals_maximum_flank=500,
                per_query_audit_path_ids_unique=True, audit_records=audit_records,
                checked_predictions=checked, diagnostic_stops=dict(stops))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    binary, output = args.binary.resolve(), args.output.resolve()
    assert binary.is_file(), binary
    if output.exists() and any(output.iterdir()):
        assert (output / "fixture.json").is_file() and not (output / "graph_build").exists(), (
            "Use a new directory, or one created with --prepare-only")
    output.mkdir(parents=True, exist_ok=True)
    fixture, commands = make_fixture(output, binary)
    if args.prepare_only:
        print(json.dumps({"prepared": str(output), "native_processes_started": 0}, indent=2))
        return
    graph_hashes = None
    for name, command in commands.items():
        with (output / (name + ".log")).open("w") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        if name == "build":
            graph_hashes = {str(path): digest(path) for path in (output / "graph_build/graph").glob("graph.sdbg*")}
    result = check_outputs(output, fixture)
    assert graph_hashes and all(digest(Path(path)) == sha for path, sha in graph_hashes.items()), "graph changed"
    result["graph_sha256"] = graph_hashes
    result["binary_sha256"] = fixture["binary_sha256"]
    (output / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
