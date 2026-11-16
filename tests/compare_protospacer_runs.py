#!/usr/bin/env python3
"""Compare complete runs on identical input queries; never select prediction paths."""
import argparse
import csv
import json
import statistics
from pathlib import Path

from evaluate_protospacer_workflow import fasta, reverse_complement


def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def load_run(predictions, evaluation):
    """Resolve unique sequence matches back to all associated input query IDs."""
    report = json.loads((evaluation / "evaluation.json").read_text())
    summary = json.loads((predictions / "run_summary.json").read_text())
    if summary["status"] != "complete":
        raise ValueError("comparison requires complete native predictions")
    queries = {r["query_id"]: r for r in table(predictions / "queries.tsv")}
    matches, coverage = {}, {}
    for name, column in (("Protospacers", "target_id"), ("ContiguousPaths", "extension_id")):
        rows = {r["sequence_id"]: r for r in table(evaluation / (name + "_coverage.tsv"))}
        passing = {key for key, value in rows.items() if value["passes_threshold"] == "1"}
        matches[name] = set()
        for link in table(predictions / "associations.tsv"):
            key = link[column]
            if key in passing:
                matches[name].add(link["query_id"])
            if name == "ContiguousPaths" and key:
                coverage[link["query_id"]] = rows[key]
    counts = {key: summary[key] for key in ("spacer_observations", "distinct_spacers", "detected_queries",
                                          "protospacer_regions", "extended_queries", "unique_extensions")}
    counts.update(matching_targets=report["Protospacers"]["matching_sequences"],
                  matching_extensions=report["ContiguousPaths"]["matching_sequences"],
                  matching_extension_queries=len(matches["ContiguousPaths"]))
    # Count lengths once per saved sequence, not once per associated input query.
    passing_extensions = [r for r in rows.values() if r["passes_threshold"] == "1"]
    lengths = [int(r["query_length"]) for r in passing_extensions]
    counts.update(matching_extensions_at_least_750_bp=sum(n >= 750 for n in lengths),
                  median_matching_extension_length=statistics.median(lengths) if lengths else None,
                  matching_extension_aligned_bases=sum(abs(int(r["query_end"])-int(r["query_begin"]))+1
                                                       for r in passing_extensions))
    return queries, matches, coverage, counts, fasta(predictions / "ContiguousPaths.fasta")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("baseline", "baseline-evaluation", "revised", "revised-evaluation", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    old = load_run(args.baseline, args.baseline_evaluation)
    new = load_run(args.revised, args.revised_evaluation)
    if {q: r["spacer"] for q, r in old[0].items()} != {q: r["spacer"] for q, r in new[0].items()}:
        raise ValueError("input query identities differ; cohorts cannot be compared")
    if args.output.exists():
        raise ValueError("comparison output already exists")
    args.output.mkdir(parents=True)
    comparison = dict(baseline=old[3], revised=new[3], identical_input_queries=True,
                      counting_unit="Query IDs for gains/losses; unique saved sequences for FASTA totals.")
    for kind in old[1]:
        a, b = old[1][kind], new[1][kind]
        comparison[kind] = dict(gained=sorted(b-a), lost=sorted(a-b), retained=sorted(a & b),
                                gained_count=len(b-a), lost_count=len(a-b), net_gain=len(b)-len(a))
    rows = []
    preserved = []
    for query in old[0]:
        row = dict(query_id=query, spacer=old[0][query]["spacer"])
        for prefix, run in (("baseline", old), ("revised", new)):
            hit = run[2].get(query, {})
            row.update({prefix + "_status": run[0][query]["status"],
                        prefix + "_extension_matches": int(query in run[1]["ContiguousPaths"]),
                        prefix + "_extension_id": hit.get("sequence_id", ""),
                        prefix + "_reference": hit.get("reference", ""),
                        prefix + "_extension_length": hit.get("query_length", ""),
                        prefix + "_identity": hit.get("identity", ""),
                        prefix + "_coverage": hit.get("coverage", "")})
        before = old[4].get(old[2].get(query, {}).get("sequence_id", ""), "")
        after = new[4].get(new[2].get(query, {}).get("sequence_id", ""), "")
        contained = bool(before and (before in after or reverse_complement(before) in after))
        row["baseline_extension_sequence_preserved"] = int(contained) if before else ""
        if contained and query in old[1]["ContiguousPaths"]:
            preserved.append(query)
        rows.append(row)
    comparison["matching_extension_sequence_preservation"] = dict(
        counting_unit="Baseline reference-matching input queries",
        definition="The complete previous saved extension is contained on either strand in the new saved extension.",
        preserved=sorted(preserved),
        shortened_changed_or_missing=sorted(old[1]["ContiguousPaths"]-set(preserved)))
    (args.output / "comparison.json").write_text(json.dumps(comparison, indent=2) + "\n")
    with (args.output / "query_comparison.tsv").open("w") as stream:
        writer = csv.DictWriter(stream, list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(comparison, indent=2))


if __name__ == "__main__":
    main()
