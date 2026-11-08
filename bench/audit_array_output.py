#!/usr/bin/env python3
"""Audit reported arrays without changing them. Requires Biopython.

Sequence similarity is 1 - global Levenshtein distance / max(sequence lengths).
Diagnostic flags are not a trained classifier or a replacement for ground truth.
"""
import argparse
from collections import Counter, defaultdict
import csv
import itertools
import json
from pathlib import Path
import statistics

from Bio.Align import PairwiseAligner


def reverse(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def rotation(sequence):
    return min(sequence[i:] + sequence[:i] for i in range(len(sequence)))


def cycle_key(sequence):
    return min(rotation(sequence), rotation(reverse(sequence)))


def kmers(sequence, k):
    return {sequence[i:i+k] for i in range(len(sequence)-k+1)}


def ids(text):
    return tuple(map(int, text.split(","))) if text != "-" else ()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.input / "CRISPR_graph_systems.tsv").open() as handle:
        systems = {r["array_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
    with (args.input / "CRISPR_graph_summary.tsv").open() as handle:
        run_summary = {r["metric"]: r["value"] for r in csv.DictReader(handle, delimiter="\t")}
    k = int(run_summary["k"])
    observations = defaultdict(list)
    cycle_owners = defaultdict(dict)
    labels, edges = {}, {}
    inconsistencies = Counter()
    with (args.input / "CRISPR_graph_paths.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            name = row["array_id"]
            path = ids(row["cycle_vertices"])
            repeat, spacer = row["repeat"], row["spacer"]
            sequence = repeat + spacer
            padded = sequence + sequence[:k-1]
            inconsistencies["cycle_length"] += len(sequence) != len(path)
            inconsistencies["non_simple_cycle"] += len(set(path)) != len(path)
            for field, expected in (
                ("repeat_vertices", path[:len(repeat)-k+1]),
                ("spacer_vertices", tuple(v for i, v in enumerate(path) if i >= len(repeat) and i+k <= len(path))),
                ("repeat_to_spacer_vertices", tuple(v for i, v in enumerate(path) if i < len(repeat) and i+k > len(repeat))),
                ("spacer_to_repeat_vertices", tuple(v for i, v in enumerate(path) if i+k > len(path))),
            ):
                inconsistencies[field] += ids(row[field]) != expected
            for pos, node in enumerate(path):
                label = padded[pos:pos+k]
                inconsistencies["node_label"] += node in labels and labels[node] != label
                labels[node] = label
            for pos, edge in enumerate(ids(row["cycle_edges"])):
                transition = (path[pos], path[(pos+1) % len(path)])
                inconsistencies["edge_transition"] += edge in edges and edges[edge] != transition
                edges[edge] = transition
            pure_r = {min(w, reverse(w)) for w in kmers(repeat, k)}
            pure_s = {min(w, reverse(w)) for w in kmers(spacer, k)}
            key = cycle_key(sequence)
            cycle_owners[key][name] = (pure_r, pure_s)
            observations[name].append((repeat, spacer, int(row["supports_repeat_path"])))
    for name, rows in observations.items():
        inconsistencies["assigned_count"] += len(rows) != int(systems[name]["assigned_cycles"])
        inconsistencies["support_count"] += sum(r[2] for r in rows) != int(systems[name]["support_cycles"])
    if set(systems) != set(observations):
        raise ValueError("Systems table and path table contain different candidate IDs")
    inconsistencies["different_vertex_ids_with_same_label"] = len(labels) - len(set(labels.values()))

    parent = {name: name for name in systems}

    def root(name):
        while parent[name] != name:
            parent[name] = parent[parent[name]]
            name = parent[name]
        return name

    shared_counts = Counter()
    conflicting_counts = Counter()
    shared_rows = []
    for sequence, owners in cycle_owners.items():
        if len(owners) < 2:
            continue
        names = sorted(owners, key=lambda s: int(s.split("_")[1]))
        conflicts = []
        for left, right in itertools.combinations(names, 2):
            parent[root(right)] = root(left)
            lr, ls = owners[left]
            rr, rs = owners[right]
            if lr & rs or rr & ls:
                conflicts.append(f"{left}:{right}")
                conflicting_counts[left] += 1
                conflicting_counts[right] += 1
        for name in names:
            shared_counts[name] += 1
        shared_rows.append({"circular_sequence": sequence, "arrays": ",".join(names),
                            "repeat_spacer_conflicts": ",".join(conflicts)})
    components = defaultdict(list)
    for name in systems:
        components[root(name)].append(name)

    aligner = PairwiseAligner(mode="global", match_score=0, mismatch_score=-1,
                              open_gap_score=-1, extend_gap_score=-1)

    def similarity(a, b):
        return 1 + aligner.score(a, b) / max(len(a), len(b))

    rows_out = []
    for name, system in systems.items():
        spacers = sorted({r[1] for r in observations[name]})
        repeat = system["repeat"]
        nearest = [0.] * len(spacers)
        total, pairs = 0., 0
        for i, a in enumerate(spacers):
            for j in range(i):
                score = similarity(a, spacers[j])
                total += score; pairs += 1
                nearest[i] = max(nearest[i], score)
                nearest[j] = max(nearest[j], score)
        repeat_similarities = [similarity(repeat, s) for s in spacers]
        repeat_words = kmers(repeat, 7)
        motif_fractions = [len(repeat_words & kmers(s, 7)) / len(repeat_words) for s in spacers]
        near = sum(x >= .90 - 1e-12 for x in nearest)
        repeat_like = sum(x >= .80 - 1e-12 for x in motif_fractions)
        flags = []
        if pairs and total/pairs >= .6:
            flags.append("high_mean_spacer_similarity")
        if near >= len(spacers)/2:
            flags.append("majority_have_near_identical_variant")
        if repeat_like >= len(spacers)/2:
            flags.append("majority_share_repeat_motif")
        if shared_counts[name]:
            flags.append("cycle_evidence_shared_with_other_candidates")
        rows_out.append({"array_id": name, "repeat": repeat, "repeat_length": len(repeat),
                         "spacer_sequences": len(spacers), "near_identical_neighbor_90": near,
                         "mean_pairwise_spacer_similarity": round(total/pairs, 6) if pairs else 0,
                         "mean_repeat_spacer_similarity": round(statistics.mean(repeat_similarities), 6),
                         "spacers_sharing_80pct_repeat_7mers": repeat_like,
                         "shared_cycles": shared_counts[name],
                         "conflicting_repeat_spacer_assignments": conflicting_counts[name],
                         "overlap_component_size": len(components[root(name)]),
                         "diagnostic_flags": ";".join(flags)})
    summary = {
        "arrays": len(systems), "observations": sum(map(len, observations.values())),
        "spacer_sequences_summed_per_array": sum(r["spacer_sequences"] for r in rows_out),
        "unique_circular_sequences_including_reverse_complements": len(cycle_owners),
        "circular_sequences_shared_by_multiple_candidates": len(shared_rows),
        "arrays_with_shared_cycles": len(shared_counts),
        "arrays_with_repeat_spacer_role_conflicts": len(conflicting_counts),
        "cycle_overlap_components": len(components),
        "largest_cycle_overlap_component": max(map(len, components.values()), default=0),
        "spacers_with_near_identical_neighbor_90": sum(r["near_identical_neighbor_90"] for r in rows_out),
        "arrays_with_high_mean_spacer_similarity": sum("high_mean_spacer_similarity" in r["diagnostic_flags"] for r in rows_out),
        "arrays_with_majority_repeat_like_spacers": sum(r["spacers_sharing_80pct_repeat_7mers"] >= r["spacer_sequences"]/2 for r in rows_out),
        "internal_inconsistencies": dict(inconsistencies),
        "interpretation": "Diagnostic properties only; overlap components are not physical loci and flags are not ground-truth false positives."
    }
    for filename, rows, empty_fields in [
        ("array_audit.tsv", rows_out, ["array_id", "repeat", "spacer_sequences", "diagnostic_flags"]),
        ("shared_cycles.tsv", shared_rows, ["circular_sequence", "arrays", "repeat_spacer_conflicts"])
    ]:
        with (args.output / filename).open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else empty_fields, delimiter="\t")
            writer.writeheader(); writer.writerows(rows)
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (args.output / "overlap_components.json").write_text(json.dumps(sorted(components.values(), key=len, reverse=True), indent=2) + "\n")
    print(json.dumps(summary, indent=2))



if __name__ == "__main__":
    main()
