#!/usr/bin/env python3
"""Reconcile per-location certificates with unchanged production predictions.

This script never labels a work-limited search negative. The initial repeat
filter's results are deliberately excluded from the main classification.
"""
import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


def table(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path, rows):
    assert rows
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fasta(path):
    result = {}
    name = None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            assert name not in result
            result[name] = ""
        elif line:
            assert name is not None
            result[name] += line.strip()
    return result


def rc(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def weighted_alignment(reference, sequence):
    """Independent full string DP, with an explicit alignment traceback."""
    matrix = [[0] * (len(reference) + 1) for _ in range(len(sequence) + 1)]
    for i in range(len(sequence) + 1):
        matrix[i][0] = 9 * i
    for j in range(len(reference) + 1):
        matrix[0][j] = 10 * j
    for i, query_base in enumerate(sequence, 1):
        for j, ref_base in enumerate(reference, 1):
            matrix[i][j] = min(matrix[i - 1][j] + 9, matrix[i][j - 1] + 10,
                               matrix[i - 1][j - 1] + (0 if query_base == ref_base else 10))
    i, j = len(sequence), len(reference)
    operations = []
    while i or j:
        if i and j and matrix[i][j] == matrix[i - 1][j - 1] + (0 if sequence[i - 1] == reference[j - 1] else 10):
            operations.append("M" if sequence[i - 1] == reference[j - 1] else "X")
            i, j = i - 1, j - 1
        elif i and matrix[i][j] == matrix[i - 1][j] + 9:
            operations.append("I")
            i -= 1
        else:
            operations.append("D")
            j -= 1
    operations.reverse()
    return matrix[-1][-1], "".join(operations)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("audit", type=Path)
    parser.add_argument("previous", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--allow-unresolved", action="store_true")
    args = parser.parse_args()
    root, previous, output = args.audit, args.previous, args.output
    targets = table(root / "adapter-free/targets.tsv")
    seeds = {row["seed_id"]: row["sequence"] for row in targets}
    assert len(seeds) == len(targets) == 81
    cases = {row["case_id"]: row for row in table(root / "adapter-free/reference_cases.tsv")}
    assert len(cases) == 1041
    expected = {(case, side) for case in cases for side in ("left", "right")}
    negative = "no_path_at_90pct_anchored_identity"
    positive = "verified_graph_path"
    chosen = {}
    # Later stages use the same or fewer graph walks. No initial source-repeat
    # sequence filter is used in this classification.
    stages = ["adapter-free", "adapter-deep", "strong-c-rich-no-prefix", "cycle-only", "unitig-cycle", "unitig-cycle-deep", "strong-three", "rolling-three"]
    for stage in stages:
        path = root / stage / "oracle_routes.tsv"
        if not path.exists():
            continue
        rows = table(path)
        assert len({(row["case_id"], row["side"]) for row in rows}) == len(rows)
        for original in rows:
            row = dict(original, evidence_stage=stage)
            key = row["case_id"], row["side"]
            assert key in expected
            if key in chosen and chosen[key]["status"] == negative:
                assert row["status"] != positive, ("contradictory certificate", key, stage)
                continue
            if row["status"] == "search_limit_unresolved" and key in chosen:
                continue
            chosen[key] = row
    assert set(chosen) == expected
    exact = table(root / "all-locations/exact_native_routes.tsv")
    exact_by_seed = defaultdict(list)
    for row in exact:
        exact_by_seed[row["seed_id"]].append(row)

    selected = {row["seed_id"]: row for row in table(previous / "final-trace-v8/production_selected.tsv")}
    metrics = table(previous / "final-trace-v8/selected_reference_metrics.tsv")
    eligible_metrics = defaultdict(list)
    for row in metrics:
        if row["qualifies_90pct"] == "1" and int(row["checkpoint_ref_bp"]) >= 100:
            eligible_metrics[row["seed_id"]].append(row)
    run_rows = table(previous / "final/per_seed.tsv")
    runs = {(row["run"], row["seed_id"]): row for row in run_rows}
    regions = fasta(previous / "final/N32/ContiguousPaths.fasta")
    witnesses = []
    for seed_id, values in eligible_metrics.items():
        metric = min(values, key=lambda row: (int(row["checkpoint_ref_bp"]), row["case_id"], row["side"]))
        prediction = selected[seed_id]
        sequence = prediction["sequence"]
        begin, end = int(prediction["seed_begin"]), int(prediction["seed_end"])
        assert sequence[begin:end] == seeds[seed_id]
        emitted = regions[runs[("N32", seed_id)]["region_id"]]
        assert sequence in (emitted, rc(emitted))
        flank = sequence[:begin][::-1] if metric["side"] == "left" else sequence[end:]
        length, query_length = int(metric["checkpoint_ref_bp"]), int(metric["best_prefix_query_bp"])
        reference = cases[metric["case_id"]][metric["side"] + "_outward"][:length]
        cost, operations = weighted_alignment(reference, flank[:query_length])
        assert cost == int(metric["weighted_cost"]) <= length
        assert operations.count("M") * 100 >= len(operations) * 90
        assert chosen[(metric["case_id"], metric["side"])]["status"] != negative
        witnesses.append(dict(seed_id=seed_id, case_id=metric["case_id"], reference=cases[metric["case_id"]]["reference"],
                              side=metric["side"], reference_bp=length, query_bp=query_length, weighted_cost=cost,
                              identity=100 * operations.count("M") / len(operations), operations=operations,
                              observed_flank=flank[:query_length], reference_flank=reference))
    production_positive = {row["seed_id"] for row in witnesses}

    for row in chosen.values():
        case = cases[row["case_id"]]
        reference = case[row["side"] + "_outward"]
        known = len(reference)
        for i, base in enumerate(reference):
            if base not in "ACGT":
                known = i
                break
        row["known_reference_bp"] = known
        row["minimum_tested_bp"] = 100
        if known < 100:
            row["status"] = "insufficient_known_reference_flank"
        if row["status"] == negative:
            assert row["exhausted"] == "1" and row["limited"] == "0"
    by_seed = defaultdict(list)
    for row in chosen.values():
        by_seed[row["seed_id"]].append(row)

    records, unresolved = [], []
    for seed_id, sequence in seeds.items():
        rows = by_seed[seed_id]
        usable = [row for row in rows if int(row["known_reference_bp"]) >= 100]
        assert usable
        unknown = [row for row in usable if row["status"] not in (negative, positive)]
        oracle_positive = any(row["status"] == positive for row in usable)
        if seed_id in production_positive:
            status = "recovered_by_current_production"
        elif oracle_positive:
            status = "additional_graph_path_needs_production_fix"
            unresolved.append(seed_id)
        elif unknown:
            status = "unresolved"
            unresolved.append(seed_id)
        else:
            status = "no_qualifying_adapter_free_source_cycle_free_path"
        longest = max(exact_by_seed[seed_id], key=lambda row: int(row["exact_bp"]))
        records.append(dict(seed_id=seed_id, seed_length=len(sequence), classification=status,
                            reference_locations=len(rows) // 2, known_flanks_tested=len(usable),
                            reference_flanks_shorter_than_100=len(rows)-len(usable), unresolved_alternative_flanks=len(unknown),
                            longest_exact_flank_bp=int(longest["exact_bp"]), exact_reference=longest["reference"],
                            exact_side=longest["side"], exact_stop_vertex=longest["vertex"], expected_next_base=longest["expected_base"],
                            available_edges=longest["available_edges"],
                            N16_length=runs[("N16", seed_id)]["length"], N32_length=runs[("N32", seed_id)]["length"],
                            N64_length=runs[("N64", seed_id)]["length"],
                            N16_local_reference_match=runs[("N16", seed_id)]["match_90_50"],
                            N32_local_reference_match=runs[("N32", seed_id)]["match_90_50"],
                            N64_local_reference_match=runs[("N64", seed_id)]["match_90_50"]))
    preserved = {}
    expected_hash = "9a7886832f257d1d7e4bc69f0cd64aab9354195dbdad47a456e5c6b1376cff12"
    for run in ("N16", "N32", "N64"):
        path = previous / "final" / run / "Protospacers.fasta"
        assert fasta(path) == seeds
        preserved[run] = hashlib.sha256(path.read_bytes()).hexdigest()
        assert preserved[run] == expected_hash
    summary = dict(inputs=81, reference_locations=len(cases), classification_counts=dict(Counter(row["classification"] for row in records)),
                   unresolved_input_ids=unresolved, classifications_complete=not unresolved,
                   criterion="At least one seed-anchored flank with 100..750 known reference bases and >=90% aligned identity; unchanged input seed; no default adapters or source-array-only closed walk.",
                   production_positive_inputs=len(production_positive), production_positive_ids=sorted(production_positive),
                   initial_repeat_filter_used_in_main_classification=False,
                   reference_use="Offline diagnosis and evaluation only; production scoring has no reference input.",
                   all_inputs_preserved_sha256=preserved,
                   reference_gap_note="Each location is tested only up to the first unknown reference base or reference end. Unknown DNA beyond a reference gap is not asserted absent.")
    print(json.dumps(summary, indent=2))
    if unresolved and not args.allow_unresolved:
        raise SystemExit("Refusing to publish a complete classification while input IDs remain unresolved")
    output.mkdir(parents=True, exist_ok=True)
    write_table(output / "all81_recovery.tsv", records)
    write_table(output / "current_prediction_witnesses.tsv", witnesses)
    write_table(output / "reference_location_certificates.tsv", [chosen[key] for key in sorted(chosen)])
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
