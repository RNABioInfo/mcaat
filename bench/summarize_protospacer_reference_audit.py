#!/usr/bin/env python3
"""Summarize all supplied seeds; never convert a bounded search into an absence claim."""
import argparse
import collections
import csv
import json
from pathlib import Path


def table(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_table(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--audit", type=Path, required=True)
    p.add_argument("--oracle", type=Path, required=True)
    p.add_argument("--before", type=Path, required=True)
    p.add_argument("--after", type=Path, required=True)
    p.add_argument("--benchmark", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    targets = table(a.audit / "targets.tsv")
    cases = table(a.oracle / "reference_cases.tsv")
    exact = table(a.oracle / "exact_native_routes.tsv")
    routes = table(a.oracle / "oracle_routes.tsv")
    witnesses = table(a.oracle / "oracle_paths.tsv")
    for r in routes:
        assert r.get("lowest_unexplored_cost"), "oracle output is incomplete"
    assert len(exact) == 2 * len(cases), "not every reference side was audited"
    by_seed = lambda rows: {t["seed_id"]: [r for r in rows if r["seed_id"] == t["seed_id"]] for t in targets}
    exacts, outcomes = by_seed(exact), by_seed(routes)
    before = {r["seed_id"]: r for r in table(a.before / "production_selected.tsv")}
    after = {r["seed_id"]: r for r in table(a.after / "production_selected.tsv")}
    metrics = {}
    for name, directory in [("before", a.before), ("after", a.after)]:
        metrics[name] = by_seed(table(directory / "selected_reference_metrics.tsv"))
    comparisons = table(a.benchmark / "per_seed.tsv")
    matches = {(r["run"], r["seed_id"]): r for r in comparisons}
    losses = by_seed(table(a.before / "verified_path_losses.tsv"))
    rows = []
    missing_examples = []
    for t in targets:
        sid, seq = t["seed_id"], t["sequence"]
        es, rs = exacts[sid], outcomes[sid]
        found = [r for r in rs if r["status"] == "verified_graph_path"]
        mapped = any(r["seed_mapped"] == "1" for r in es)
        best = max(es, key=lambda r: int(r["exact_bp"]))
        if len(seq) < 23:
            finding = "input_shorter_than_graph_k"
        elif not mapped:
            finding = "original_seed_not_mapped"
        elif any(int(r["checkpoint_ref_bp"]) >= 24 for r in found):
            finding = "connected_reference_consistent_extension_verified"
        else:
            short = [r for r in rs if r["checkpoint_ref_bp"] == "24"]
            finding = "no_90pct_24bp_flank_at_tested_locations" if short and all(r["status"] == "no_path_at_90pct_anchored_identity" for r in short) else "no_verified_24bp_flank_search_unresolved"
        row = dict(seed_id=sid, seed_length=len(seq), reference_locations=len(es)//2,
                   graph_seed_mapped=mapped, finding=finding,
                   longest_exact_left=max(int(r["exact_bp"]) for r in es if r["side"] == "left"),
                   longest_exact_right=max(int(r["exact_bp"]) for r in es if r["side"] == "right"))
        for side in ["left", "right"]:
            row["longest_verified_90pct_" + side] = max((int(r["checkpoint_ref_bp"]) for r in found if r["side"] == side), default=0)
        row.update(unresolved_reference_checkpoints=sum(r["status"] == "search_limit_unresolved" for r in rs),
                   longest_exact_case=best["case_id"], longest_exact_reference=best["reference"],
                   longest_exact_side=best["side"], exact_first_stop=best["first_stop"],
                   exact_stop_after_bp=best["exact_bp"], exact_stop_vertex=best["vertex"],
                   exact_required_base=best["expected_base"], exact_available_edges=best["available_edges"])
        for label, selected, run in [("before", before, "baseline"), ("after", after, "N32")]:
            r = selected[sid]
            row[label + "_length"] = len(r["sequence"])
            row[label + "_reason"] = r["reason"]
            row[label + "_reference_match_90id_50cov"] = matches[run, sid]["match_90_50"]
            ms = [m for m in metrics[label][sid] if m["qualifies_90pct"] == "1"]
            for side in ["left", "right"]:
                row[label + "_verified_ref_" + side] = max((int(m["checkpoint_ref_bp"]) for m in ms if m["side"] == side), default=0)
        row["before_observed_path_outcomes"] = ";".join(sorted({r["status"] for r in losses[sid]}))
        rows.append(row)
        if best["first_stop"] == "exact_edge_absent":
            case = next(c for c in cases if c["case_id"] == best["case_id"])
            n = int(best["exact_bp"])
            if best["side"] == "right":
                sequence = seq + case["right_outward"][:n]
                vertex = sequence[-23:]
                required = vertex[1:] + best["expected_base"]
            else:
                sequence = case["left_outward"][:n][::-1] + seq
                vertex = sequence[:23]
                required = best["expected_base"] + vertex[:-1]
            missing_examples.append(dict(seed_id=sid, case_id=best["case_id"], reference=best["reference"], side=best["side"],
                                         matched_flank_bp=n, vertex=best["vertex"], vertex_kmer=vertex,
                                         required_next_kmer=required, available_edges=best["available_edges"]))
    write_table(a.output / "all81_audit.tsv", rows)
    write_table(a.output / "exact_missing_edge_examples.tsv", missing_examples)
    summary = dict(inputs=len(targets), reference_locations=len(cases), sides=len(exact),
                   findings=dict(collections.Counter(r["finding"] for r in rows)),
                   confirmed_inputs_by_reference_flank_length={n: len({r["seed_id"] for r in routes if r["status"] == "verified_graph_path" and int(r["checkpoint_ref_bp"]) >= n}) for n in [24,100,250,500,750]},
                   exact_inputs_by_flank_length={n: len({r["seed_id"] for r in exact if int(r["exact_bp"]) >= n}) for n in [24,100,250,500,750]},
                   route_statuses=dict(collections.Counter(r["status"] for r in routes)),
                   unresolved_inputs=len({r["seed_id"] for r in routes if r["status"] == "search_limit_unresolved"}),
                   reference_used_for_production_ranking=False)
    (a.output / "audit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    # Emit actual graph-spelled sequence in the original input orientation.
    # Each record is ONE flank, not an unverified join of two reference paths.
    seeds = {t["seed_id"]: t["sequence"] for t in targets}
    with (a.output / "verified_graph_flanks.fasta").open("w") as f, (a.output / "verified_seed_anchored_paths.fasta").open("w") as anchored:
        for r in witnesses:
            flank = r["sequence_outward"][::-1] if r["side"] == "left" else r["sequence_outward"]
            name = r["case_id"] + "|" + r["seed_id"] + "|" + r["side"] + "|ref_bp=" + r["checkpoint_ref_bp"] + "|input_forward"
            f.write(">" + name + "\n" + flank + "\n")
            seed = seeds[r["seed_id"]]
            seq = flank + seed if r["side"] == "left" else seed + flank
            anchored.write(">" + name + "\n" + seq + "\n")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
