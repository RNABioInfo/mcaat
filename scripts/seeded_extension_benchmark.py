#!/usr/bin/env python3
"""Evaluate completed seed-local runs and draw two compact scientific figures.

BLAST is evaluation only. This script never selects or modifies prediction paths.
Temporary query FASTAs are deleted; best-hit evidence is added to extensions.tsv.
"""
import argparse
import csv
import hashlib
import json
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

METHODS = ("poisson", "negbin", "huber")
LABELS = {"poisson": "Poisson", "negbin": "Negative binomial", "huber": "Huber log-ratio"}
COLORS = {"poisson": "#245b91", "negbin": "#b66925", "huber": "#6d5295"}


def fasta(path):
    records, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            assert name not in records, (path, name)
            records[name] = ""
        elif line.strip():
            assert name is not None, path
            records[name] += line.strip()
    return records


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def canonical(sequence):
    return min(sequence, sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1])


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024**2), b""):
            result.update(block)
    return result.hexdigest()


def same_region(a, b, k=23):
    """Mirror the native direct-representative rule for a fair baseline count."""
    if min(len(a), len(b)) < k:
        return a == b
    index = defaultdict(list)
    for i in range(len(a) - k + 1):
        index[a[i:i + k]].append(i)
    for sequence in (b, b.translate(str.maketrans("ACGT", "TGCA"))[::-1]):
        offsets = {i - j for j in range(len(sequence) - k + 1) for i in index.get(sequence[j:j + k], [])}
        for offset in offsets:
            ab, bb = max(0, offset), max(0, -offset)
            overlap = min(len(a) - ab, len(sequence) - bb)
            if 10 * overlap < 9 * min(len(a), len(sequence)):
                continue
            identical = longest = run = 0
            for x, y in zip(a[ab:ab + overlap], sequence[bb:bb + overlap]):
                identical += x == y
                run = run + 1 if x == y else 0
                longest = max(longest, run)
            if 100 * identical >= 95 * overlap and longest >= k and 10 * (longest - k + 1) >= 9 * (min(len(a), len(sequence)) - k + 1):
                return True
    return False


def align(records, args, target=False):
    """Stream HSPs; prefer a passing HSP, then coverage, identity and bitscore."""
    if not records:
        return {}, []
    minimum = 85 if target else 50
    with tempfile.TemporaryDirectory(prefix="mcaat-evaluation-") as temporary:
        query = Path(temporary) / "queries.fasta"
        query.write_text("".join(f">{key}\n{sequence}\n" for key, sequence in records.items()))
        command = [args.blastn, "-query", str(query), "-db", str(args.reference),
                   "-task", "blastn-short", "-word_size", "7" if target else "9",
                   "-dust", "no", "-soft_masking", "false", "-evalue", "1000",
                   "-perc_identity", "90", "-max_target_seqs", "1000000",
                   "-num_threads", str(args.threads), "-outfmt",
                   "6 qseqid sseqid qlen qstart qend sstart send length nident pident evalue bitscore"]
        best = {}
        with (Path(temporary) / "stderr").open("w+") as error:
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=error, text=True)
            for line in process.stdout:
                q, s, qlen, qb, qe, sb, se, length, same, identity, evalue, bits = line.split()
                qlen, qb, qe, sb, se, length, same = map(int, (qlen, qb, qe, sb, se, length, same))
                assert qlen == len(records[q])
                span = abs(qe - qb) + 1
                passed = 10 * same >= 9 * length and 100 * span >= minimum * qlen
                item = dict(passes=passed, reference=s, identity=same / length,
                            coverage=span / qlen, query_begin=qb, query_end=qe,
                            reference_begin=sb, reference_end=se, aligned_bases=span,
                            alignment_length=length, identical_bases=same, bitscore=float(bits))
                rank = (passed, span, same / length, float(bits))
                if q not in best or rank > best[q][0]:
                    best[q] = (rank, item)
            status = process.wait()
            if status:
                error.seek(0)
                raise RuntimeError(error.read())
        command[command.index("-query") + 1] = "<temporary deduplicated sequences>"
        return {key: best[key][1] if key in best else dict(passes=False, reference="", identity=0.,
                coverage=0., query_begin=0, query_end=0, reference_begin=0, reference_end=0,
                aligned_bases=0, alignment_length=0, identical_bases=0, bitscore=0.) for key in records}, command


def reference_loci(hits):
    """Group overlapping best qualifying HSP intervals; these are not phage counts."""
    intervals = defaultdict(list)
    for hit in hits:
        if hit["passes"]:
            intervals[hit["reference"]].append(sorted((hit["reference_begin"], hit["reference_end"])))
    loci = []
    for reference, locations in sorted(intervals.items()):
        for begin, end in sorted(locations):
            if loci and loci[-1]["reference"] == reference and begin <= loci[-1]["end"]:
                loci[-1]["end"] = max(end, loci[-1]["end"])
                loci[-1]["predictions"] += 1
            else:
                loci.append(dict(reference=reference, begin=begin, end=end, predictions=1))
    return loci


def make_plots(output, evidence, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titleweight": "bold", "figure.facecolor": "white"})
    x = np.linspace(0, 4, 401)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    r = summary["nb_dispersion"]
    for ax, method in zip(axes, METHODS):
        for mean, style in ((10, "-"), (50, "--")):
            count = x * mean
            term = np.zeros_like(count)
            nonzero = count > 0
            term[nonzero] = count[nonzero] * np.log(count[nonzero] / mean)
            if method == "poisson":
                loss = 2 * (term - count + mean)
            elif method == "negbin":
                loss = 2 * (term - (count + r) * np.log((count + r) / (mean + r)))
            else:
                z = np.abs(np.log((count + .5) / (mean + .5)))
                loss = np.where(z <= 1, .5 * z**2, z - .5)
            ax.plot(x, loss, style, color=COLORS[method], lw=2.3, label=fr"Baseline $\mu={mean}$")
        ax.axvline(1, color="#a4acb4", lw=1)
        ax.set(title=LABELS[method], xlabel=r"Observed multiplicity / baseline  ($c/\mu$)",
               ylabel="Penalty per scored position", xlim=(0, 4), ylim=(0, None))
        ax.grid(alpha=.16)
        ax.legend(frameon=False, fontsize=10)
    fig.suptitle("Three coverage losses: lower penalty means closer to the local baseline", fontsize=15)
    fig.savefig(output / "coverage_models.png", dpi=180)
    plt.close(fig)

    representatives = {m: [r for r in evidence if r["method"] == m and r["relation_to_representative"] == "representative"] for m in METHODS}
    by_query = defaultdict(dict)
    for row in evidence:
        by_query[row["query_id"]][row["method"]] = row
    # Choose a disagreement deterministically for illustration, after prediction.
    options = sorted(q for q, group in by_query.items() if len(group) == 3
                     and all(row["accepted"] == "1" for row in group.values()))
    differing = [q for q in options if len({row["sequence"] for row in by_query[q].values()}) > 1]
    different_match = [q for q in differing if len({r["reference_pass"] for r in by_query[q].values()}) > 1]
    any_match = [q for q in differing if any(r["reference_pass"] == "1" for r in by_query[q].values())]
    chosen = (different_match or any_match or differing or options or sorted(by_query))[0]
    summary["plot_example_query"] = chosen
    summary["plot_example_policy"] = "First sorted method disagreement with differing reference outcomes, otherwise any match, otherwise first disagreement; selected after prediction."
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.3), constrained_layout=True)
    for column, method in enumerate(METHODS):
        ax = axes[0, column]
        populations = {True: [], False: []}
        for row in representatives[method]:
            for side in ("left", "right"):
                counts = [float(x) for x in row[f"{side}_counts"].split(",") if x]
                positions = [int(x) for x in row[f"{side}_scored_positions"].split(",") if x]
                populations[row["reference_pass"] == "1"].extend(counts[p - 1] for p in positions)
        maximum = max([v for values in populations.values() for v in values] or [2])
        bins = np.geomspace(1, max(2, maximum + 1), 35)
        for passed, color in ((False, "#949da6"), (True, "#16806b")):
            values = populations[passed]
            if values:
                weights = np.ones(len(values)) / len(values)
                ax.hist(values, bins=bins, weights=weights, histtype="step", linewidth=2,
                        color=color, label="Reference match" if passed else "No qualifying match")
        ax.set_xscale("log")
        ax.set(title=LABELS[method], xlabel="External edge multiplicity (log scale)", ylabel="Fraction of positions / bin")
        ax.legend(frameon=False, fontsize=9)
        ax.grid(alpha=.15)
        ax = axes[1, column]
        row = by_query[chosen][method]
        seed_length = int(row["seed_end"]) - int(row["seed_begin"])
        for side, color in (("left", "#245b91"), ("right", "#b66925")):
            counts = np.array([float(x) for x in row[f"{side}_counts"].split(",") if x])
            positions = np.array([int(x) for x in row[f"{side}_scored_positions"].split(",") if x], dtype=int)
            padding = int(row["left_bp"]) - len(counts) if side == "left" else 0
            xs = -padding - np.arange(1, len(counts) + 1) if side == "left" else seed_length + np.arange(len(counts))
            ax.plot(xs, counts, color="#c8cdd3", linewidth=.8)
            if len(positions):
                ax.plot(xs[positions - 1], counts[positions - 1], color=color, linewidth=1.4, label=f"{side.capitalize()} scored")
                baseline = float(row[f"{side}_baseline"])
                ax.hlines(baseline, min(xs), max(xs), color=color, linestyle="--", alpha=.7)
        ax.axvspan(0, seed_length, color="#64845b", alpha=.18)
        ax.set(xlim=(-750, seed_length + 750), xlabel="Position relative to target start (bp)", ylabel="Edge multiplicity",
               title=f"{chosen} · best HSP covers {100 * float(row['reference_coverage'] or 0):.1f}%")
        ax.set_yscale("symlog", linthresh=1)
        ax.set_ylim(bottom=0)
        ax.grid(alpha=.15)
        ax.legend(frameon=False, fontsize=9)
    fig.suptitle("Measured multiplicity distributions and extension paths", fontsize=15)
    fig.supxlabel("Top: representative paths only; reference labels applied afterward. Bottom: grey = masked context; dashed = frozen baseline.", fontsize=10)
    fig.savefig(output / "multiplicities.png", dpi=180)
    plt.close(fig)


def write_report(output, summary):
    evaluation = summary["evaluation"]
    baseline = evaluation["baseline"]
    lines = ["# Seed-local protospacer extension: SRR14017553", "",
             f"Tested all **{summary['input_starts']} saved query/start records**, representing "
             f"**{summary['input_target_sequences']} candidate target sequences** from "
             f"{summary['original_corrected_spacers']} corrected distinct spacers. No reference-based input selection.", "",
             f"These records use {summary['distinct_oriented_anchor_paths']} distinct oriented anchor paths. "
             "Counts of matched starts and gains/losses below refer to input records; multiple records can share a target or path.", "",
             "Each method starts at the saved target nodes, applies local MEGAHIT-derived tip/bubble/weak-link rules, "
             "and searches up to **750 added bases on each side**. Prediction uses no reference. "
             "There is no whole-metagenome assembly and the existing graph is retained.", "",
             "## Results", "",
             "| Method | Accepted starts | Matched starts | Representative contigs | Matched representatives | Gained / lost starts |",
             "|---|---:|---:|---:|---:|---:|"]
    for method in METHODS:
        result = evaluation["methods"][method]
        lines.append(f"| {LABELS[method]} | {result['accepted_starts']} | {result['matched_starts']} | "
                     f"{result['representative_contigs']} | {result['matched_representatives']} | "
                     f"+{len(result['gained_starts'])} / −{len(result['lost_starts'])} |")
    best_count = max(r["matched_starts"] for r in evaluation["methods"].values())
    if best_count <= baseline["matched_starts"]:
        lines.extend(["", "**Decision:** none of these models improves the previous matched-start count under "
                      "the requested full-contig metric. Keep this as an experimental comparison; it should not replace "
                      "the current default on the strength of this test."])
    if "initial_trial" in summary:
        trial = summary["initial_trial"]
        counts = ", ".join(f"{LABELS[m]} {trial['evaluated_methods'][m]['matched_starts']}" for m in METHODS)
        lines.extend(["", f"The first, coverage-only trial failed ({counts} matching starts). "
                      "It favored flat weakly supported paths and penalized coverage establishment near source exits. "
                      "The corrected comparison above retains branch support and scores coverage after calibration; "
                      "all starts were rerun. Initial counts and losses remain in `summary.json`."])
    if "premature_pruning_audit" in summary:
        audit = summary["premature_pruning_audit"]
        lines.extend(["", f"A direct graph trace then found premature weak-link pruning on "
                      f"{audit['matching_starts_pruned']} previously matching starts, before external coverage was calibrated. "
                      "The final run above protects those uncalibrated bridges and rare arms consistent with the local baseline. "
                      "The second comparison and all audited branch positions/counts are retained in `summary.json`."])
    lines.extend(["", "**Match criterion:** one BLAST HSP with at least 90% identity and at least 50% coverage "
                  "of the entire extension sequence. HSPs are not added together. "
                  "The alignment need not cover the target, matching the previously agreed evaluation.", "",
                  f"Previous run: {baseline['extensions']} sequence records; {baseline['extensions_at_least_500']} "
                  f"were at least 500 bp, including all {baseline['matched_extensions']} matching records. "
                  f"Those matches represented **{baseline['matched_starts']} input starts**. "
                  "Gains/losses above compare the same query IDs before grouping. "
                  "The previous run retained one flank; this run can join both, so a longer query can lower the coverage fraction.", "",
                  f"Applying the same region-grouping rule to the baseline leaves "
                  f"**{baseline['grouped_representatives']} representative contigs, {baseline['matched_grouped_representatives']} matching**. "
                  "This is the appropriate record-count comparison with the new grouped outputs; the old 49 matching records included repeated regions.", "",
                  f"Target-only evaluation: **{evaluation['target_matches']}/{evaluation['target_count']}** "
                  "match at 90% identity and 85% coverage. Every input target remains in `protospacers.fasta`, "
                  "including targets without accepted extensions.", "",
                  "Representative grouping removes redundant region records, not source associations. "
                  "The original per-start sequence, target coordinates and group assignment remain in `extensions.tsv`. "
                  "A representative can differ from an associated target variant. Counts of these partial regions are not counts of distinct phages.", "",
                  "## Lengths, limits and measured work", "",
                  "| Method | Contig length range | Starts with both full flanks | Accepted starts with a search limit | Best-hit reference loci | Search seconds |",
                  "|---|---:|---:|---:|---:|---:|"])
    for measured in summary["methods"]:
        method = measured["method"]
        result = evaluation["methods"][method]
        lines.append(f"| {LABELS[method]} | {result['representative_length_min']}–{result['representative_length_max']} bp | "
                     f"{result['full_two_flank_starts']} | {result['accepted_limited_starts']} | "
                     f"{len(result['best_hit_reference_loci'])} | {measured['extension_seconds']:.1f} |")
    if "resource_measurements" in summary:
        resources = summary["resource_measurements"]
        lines.extend(["", f"All three searches together: **{resources['elapsed_seconds'] / 60:.2f} minutes**, "
                      f"including graph loading/checks; **{resources['peak_rss_gib']:.2f} GiB peak RSS**, "
                      f"{summary['threads']} threads. No memory-guard stop. BLAST is measured separately from prediction."])
    lines.extend(["", "Best-hit loci group overlapping qualifying alignments on the same reference record. "
                  "This is a reporting aid, not an exhaustive catalog of all reference mappings. "
                  "Beam/work limits are recorded; this is not an exhaustive search. "
                  "Unmatched sequences are not automatically proven false, and flat coverage alone cannot establish viral origin. "
                  "This dataset was used during method development; these results are not independent validation.", "",
                  "## What the methods measure", "",
                  "- **Poisson:** penalizes departures from a frozen local count baseline using Poisson deviance.",
                  "- **Negative binomial:** uses an overdispersed count model, with dispersion fixed at 10 before evaluation.",
                  "- **Huber log-ratio:** measures relative count changes while reducing the effect of large deviations.", "",
                  "Weak-link removal additionally requires an established local baseline and a drop below that baseline. "
                  "The baseline is the median of the first 127 unmasked positions per flank, then frozen; "
                  "coverage residuals are scored after that calibration window. Branch support is retained separately. "
                  "Source-repeat counts are excluded. Full read-span medians must stay within a twofold range. "
                  "The beam retains 64 mature and 256 pending states. All three use the same source/revisit/artifact rules, "
                  "a minimum 500-bp contig, and at least 254 scored external positions. "
                  "The final score is mean branch surprisal plus mean coverage loss weighted by 1/127 "
                  "(deviances are divided by two; Huber is used directly). This conservative read-span weighting "
                  "avoids counting overlapping k-mer observations as independent evidence. It remains a design choice, "
                  "not a calibrated likelihood or probability of correct assembly.", "",
                  "See the repository's `docs/protospacer-extension.md` for the exact formulas, local rule adaptations and function guide.", "",
                  "## Files and figures", "",
                  "- `protospacers.fasta`: all input candidate target sequences.",
                  "- `viral_contigs_poisson.fasta`, `viral_contigs_negbin.fasta`, `viral_contigs_huber.fasta`: accepted representatives.",
                  "- `extensions.tsv`: all starts and outcomes, actual selected sequences, multiplicity traces, source associations and reference evidence.",
                  "- `summary.json`: parameters, timings, counts, exact gained/lost IDs, reference loci and provenance.",
                  "- `run.log`: native stage progress.",
                  "- `coverage_models.png`: mathematical penalty versus multiplicity ratio; each panel has its own penalty scale.",
                  "- `multiplicities.png`: measured distributions and the same example start under all three methods.", "",
                  "![Coverage losses](coverage_models.png)", "", "![Observed multiplicities](multiplicities.png)", "",
                  "BLAST query coordinates in `extensions.tsv` refer to the canonical evaluation sequence; "
                  "`reference_query_orientation` identifies its orientation relative to that row's actual graph sequence. "
                  "Native target coordinates are zero-based/end-exclusive. Raw counts are ordered outward on each flank; "
                  "scored positions index those counts from one.", "",
                  "## Validation", "",
                  "Release build and arithmetic/context self-tests passed. The controlled graph test covers all three losses, "
                  "both full flanks, target preservation, simple bubbles, short tips, weak links, source-repeat masking, "
                  "bounded probes, a flat weak-branch trap, and unchanged graph files. Every real output route is re-spelled from connected graph nodes "
                  "and checked against its saved target before writing.", ""])
    (output / "REPORT.md").write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--baseline", type=Path, required=True, help="Previous prediction directory")
    parser.add_argument("--baseline-evaluation", type=Path, required=True, help="Previous BLAST evaluation directory")
    parser.add_argument("--reference", type=Path, required=True, help="Local BLAST database prefix")
    parser.add_argument("--blastn", default="blastn")
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    output = args.output.resolve()
    summary = json.loads((output / "summary.json").read_text())
    assert summary["status"] == "complete", "Evaluate only a completed run"
    evidence = rows(output / "extensions.tsv")
    assert len(evidence) == 3 * summary["input_starts"]
    assert set(row["method"] for row in evidence) == set(METHODS)
    targets = fasta(output / "protospacers.fasta")
    seed_dir = Path(summary["seed_directory"])
    saved_anchors = {r["query_id"]: r["anchor_nodes"] for r in rows(seed_dir / "graph_paths.tsv")}
    summary["distinct_oriented_anchor_paths"] = len(set(saved_anchors.values()))
    source_rows = rows(seed_dir / "sources.tsv")
    summary["source_observations"] = len(source_rows)
    summary["source_arrays"] = len({r["array_id"] for r in source_rows})
    assert targets == fasta(seed_dir / "CandidateProtospacers.fasta"), "Input targets changed"
    for name in ("queries.tsv", "sources.tsv"):
        assert digest(seed_dir / name) == digest(args.baseline / name), "Baseline has different inputs"
    unique, sequence_ids = {}, {}
    for row in evidence:
        seq = row["sequence"]
        begin, end = int(row["seed_begin"]), int(row["seed_end"])
        assert canonical(seq[begin:end]) == targets[row["target_id"]]
        assert begin == int(row["left_bp"]) <= 750 and len(seq) - end == int(row["right_bp"]) <= 750
        if row["accepted"] == "1":
            assert len(seq) >= summary["min_contig_length"]
            sequence = canonical(seq)
            if sequence not in sequence_ids:
                key = f"S{len(unique) + 1:06d}"
                unique[key] = sequence
                sequence_ids[sequence] = key
    print(f"Evaluating {len(unique)} distinct accepted extension sequences; paths already fixed.", flush=True)
    hits, extension_command = align(unique, args)
    print(f"Evaluating {len(targets)} unchanged input targets.", flush=True)
    target_hits, target_command = align(targets, args, target=True)
    for row in evidence:
        hit = hits.get(sequence_ids.get(canonical(row["sequence"])), {}) if row["accepted"] == "1" else {}
        row["reference_evaluated"] = row["accepted"]
        row["reference_pass"] = str(int(hit.get("passes", False)))
        for field in ("reference", "identity", "coverage", "query_begin", "query_end", "reference_begin", "reference_end", "bitscore"):
            row["reference_" + field if field != "reference" else "reference_record"] = str(hit.get(field, ""))
    # Query coordinates in BLAST refer to the canonical evaluation sequence.
    for row in evidence:
        row["reference_query_orientation"] = "+" if row["sequence"] == canonical(row["sequence"]) else "-"
    baseline_sequences = fasta(args.baseline / "ContiguousPaths.fasta")
    baseline_manifest = json.loads((args.baseline_evaluation / "evaluation.json").read_text())
    baseline_metric = baseline_manifest["ContiguousPaths"]
    assert baseline_metric["identity_threshold"] == .90 and baseline_metric["single_hsp_query_coverage_threshold"] == .50, "Baseline thresholds differ"
    assert baseline_manifest["input_sha256"]["ContiguousPaths.fasta"] == digest(args.baseline / "ContiguousPaths.fasta"), "Baseline evaluation used different sequences"
    baseline_hits = {r["sequence_id"]: r for r in rows(args.baseline_evaluation / "ContiguousPaths_coverage.tsv")}
    baseline_eligible = {name for name, seq in baseline_sequences.items() if len(seq) >= summary["min_contig_length"]}
    baseline_passing = {key for key in baseline_eligible if baseline_hits[key]["passes_threshold"] == "1"}
    baseline_queries = {r["query_id"] for r in rows(args.baseline / "associations.tsv") if r["extension_id"] in baseline_passing}
    tested_queries = {r["query_id"] for r in evidence}
    assert baseline_queries <= tested_queries
    baseline_groups = {}
    for key in sorted(baseline_eligible, key=lambda key: (-len(baseline_sequences[key]), baseline_sequences[key], key)):
        representative = next((rep for rep in baseline_groups if same_region(baseline_sequences[rep], baseline_sequences[key], summary["graph_k"])), key)
        baseline_groups.setdefault(representative, []).append(key)
    evaluation = dict(reference=str(args.reference.resolve()), identity=0.90, extension_coverage=0.50,
                      target_coverage=0.85, single_hsp=True, reference_used_for_selection=False,
                      baseline=dict(extensions=len(baseline_sequences), extensions_at_least_500=len(baseline_eligible),
                                    matched_extensions=len(baseline_passing), matched_starts=len(baseline_queries),
                                    grouped_representatives=len(baseline_groups),
                                    matched_grouped_representatives=sum(key in baseline_passing for key in baseline_groups),
                                    groups=baseline_groups),
                      target_matches=sum(h["passes"] for h in target_hits.values()),
                      target_count=len(targets), target_hits=target_hits,
                      commands=[extension_command, target_command], methods={})
    for method in METHODS:
        own = [r for r in evidence if r["method"] == method]
        assert len({r["query_id"] for r in own}) == summary["input_starts"]
        accepted = [r for r in own if r["accepted"] == "1"]
        matched_queries = {r["query_id"] for r in accepted if r["reference_pass"] == "1"}
        reps = fasta(output / f"viral_contigs_{method}.fasta")
        rep_rows = [r for r in own if r["relation_to_representative"] == "representative"]
        assert reps == {r["viral_contig_id"]: r["sequence"] for r in rep_rows}
        assert len(set(map(canonical, reps.values()))) == len(reps), "Duplicate representative sequences"
        assert all(same_region(reps[r["viral_contig_id"]], r["sequence"], summary["graph_k"]) for r in accepted), "Invalid native representative association"
        rep_hits = [hits[sequence_ids[canonical(seq)]] for seq in reps.values()]
        passing_reps = [hit for hit in rep_hits if hit["passes"]]
        lengths = [len(seq) for seq in reps.values()]
        evaluation["methods"][method] = dict(
            accepted_starts=len(accepted), accepted_target_sequences=len({r["target_id"] for r in accepted}),
            representative_contigs=len(reps), matched_representatives=len(passing_reps),
            reference_match_fraction=len(passing_reps) / len(reps) if reps else 0.,
            matched_starts=len(matched_queries), gained_starts=sorted(matched_queries - baseline_queries),
            matched_oriented_anchor_paths=len({saved_anchors[q] for q in matched_queries}),
            lost_starts=sorted(baseline_queries - matched_queries), retained_starts=sorted(baseline_queries & matched_queries),
            representative_length_min=min(lengths, default=0), representative_length_max=max(lengths, default=0),
            full_two_flank_starts=sum(int(r["left_bp"]) == 750 and int(r["right_bp"]) == 750 for r in accepted),
            accepted_limited_starts=sum(r["search_limited"] == "1" for r in accepted),
            accepted_ambiguous_starts=sum(r["ambiguous"] == "1" for r in accepted),
            best_hit_reference_loci=reference_loci(rep_hits),
            stop_reasons=dict(Counter(r[f"{side}_stop"] for r in own for side in ("left", "right"))),
            local_rule_events={field: sum(int(r[field]) for r in own) for field in ("tips_removed", "bubbles_removed", "weak_links_removed")})
    summary["evaluation"] = evaluation
    summary["input_and_prediction_hashes"] = {str(p.resolve()): digest(p) for p in
        [seed_dir / "queries.tsv", seed_dir / "sources.tsv", seed_dir / "graph_paths.tsv", seed_dir / "candidates.tsv",
         output / "protospacers.fasta", *(output / f"viral_contigs_{m}.fasta" for m in METHODS)]}
    make_plots(output, evidence, summary)
    with (output / "extensions.tsv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(evidence[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(evidence)
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_report(output, summary)
    print(json.dumps({"baseline": evaluation["baseline"], "target_matches": evaluation["target_matches"],
                      "methods": evaluation["methods"]}, indent=2))


if __name__ == "__main__":
    main()
