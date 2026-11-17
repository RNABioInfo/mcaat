#!/usr/bin/env python3
"""Evaluate completed seed-local runs and draw two compact scientific figures.

BLAST is evaluation only. This script never selects or modifies prediction paths.
Temporary query FASTAs are deleted; best-hit evidence is added to extensions.tsv.
"""
import argparse
import csv
import hashlib
import gzip
import json
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

METHODS = ("poisson", "negbin", "huber")
LABELS = {"poisson": "Poisson", "negbin": "Negative binomial", "huber": "Huber log-ratio", "flow": "Flow-guided", "lookahead": "Branch lookahead"}
COLORS = {"poisson": "#245b91", "negbin": "#b66925", "huber": "#6d5295", "flow": "#16806b", "lookahead": "#245b91"}


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


def make_flow_plots(output, evidence, summary):
    """Show the flow rule and observed results after paths have been fixed."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    leakage = np.linspace(0, 60, 301)
    for initial, color in ((10, "#b66925"), (30, "#245b91"), (50, "#16806b")):
        axes[0].plot(leakage, initial - leakage, color=color, lw=2, label=f"Initial flow {initial}")
    axes[0].axhline(0, color="#555555", ls="--", lw=1)
    axes[0].axhspan(-60, 0, color="#eeeeee")
    axes[0].set(xlabel="Cumulative competing-exit flow", ylabel="Excess flow = initial − leakage",
                title="Continue only with one positive continuation", ylim=(-20, 52), xlim=(0, 60))
    axes[0].legend(frameon=False)
    representatives = [r for r in evidence if r["relation_to_representative"] == "representative"]
    for passed, color in ((False, "#949da6"), (True, "#16806b")):
        subset = [r for r in representatives if (r["reference_pass"] == "1") == passed]
        fractions = [1 - float(r["support_cost"]) for r in subset]
        if fractions:
            axes[1].hist(fractions, bins=np.linspace(0, 1, 21), histtype="step", lw=2, color=color,
                         label="Reference match" if passed else "No qualifying match")
    axes[1].set(xlabel="Worst retained flank: remaining / initial flow", ylabel="Representative contigs",
                xlim=(0, 1), title="Measured support; references evaluated afterward")
    if representatives: axes[1].legend(frameon=False)
    for ax in axes: ax.grid(alpha=.15)
    fig.suptitle("Flow-guided extension: competing exits consume seed support", fontsize=14)
    fig.savefig(output / "flow_support.png", dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    populations = {True: [], False: []}
    for row in representatives:
        for side in ("left", "right"):
            counts = [float(x) for x in row[f"{side}_counts"].split(",") if x]
            positions = [int(x) for x in row[f"{side}_scored_positions"].split(",") if x]
            populations[row["reference_pass"] == "1"].extend(counts[p - 1] for p in positions)
    maximum = max([v for values in populations.values() for v in values] or [2])
    for passed, color in ((False, "#949da6"), (True, "#16806b")):
        values = populations[passed]
        if values:
            axes[0].hist(values, bins=np.geomspace(1, max(2, maximum + 1), 35),
                         weights=np.ones(len(values)) / len(values), histtype="step", lw=2,
                         color=color, label="Reference match" if passed else "No qualifying match")
    axes[0].set_xscale("log")
    axes[0].set(xlabel="External edge multiplicity", ylabel="Fraction of positions / bin",
                title="Observed multiplicities of saved representatives")
    if representatives: axes[0].legend(frameon=False)
    example = min(representatives or evidence, key=lambda r: r["query_id"])
    summary["plot_example_query"] = example["query_id"]
    summary["plot_example_policy"] = "First saved representative by query ID; no reference-dependent example selection."
    seed_length = int(example["seed_end"]) - int(example["seed_begin"])
    for side, color in (("left", "#245b91"), ("right", "#b66925")):
        counts = np.array([float(x) for x in example[f"{side}_counts"].split(",") if x])
        padding = int(example["left_bp"]) - len(counts) if side == "left" else 0
        xs = -padding - np.arange(1, len(counts) + 1) if side == "left" else seed_length + np.arange(len(counts))
        if len(counts): axes[1].plot(xs, counts, color=color, lw=1.5, label=f"{side.capitalize()} flank")
    axes[1].axvspan(0, seed_length, color="#64845b", alpha=.18)
    axes[1].set(xlabel="Position relative to target start (bp)", ylabel="Edge multiplicity",
                title=f"{example['query_id']}: actual retained path", xlim=(-750, seed_length + 750))
    axes[1].set_yscale("symlog", linthresh=1)
    axes[1].set_ylim(bottom=0)
    if example["left_counts"] or example["right_counts"]: axes[1].legend(frameon=False)
    for ax in axes: ax.grid(alpha=.15)
    fig.suptitle("Observed graph counts; the 750-bp limit is a maximum", fontsize=14)
    fig.savefig(output / "multiplicities.png", dpi=180)
    plt.close(fig)


def write_flow_report(output, summary):
    e = summary["evaluation"]
    r = e["methods"]["flow"]
    previous = e.get("previous_huber", {})
    delta = r.get("versus_previous_huber", {})
    lines = ["# Flow-guided extension: SRR14017553", "",
        f"Processed **{summary['input_starts']} input records**, preserving all **{summary['input_target_sequences']} candidate protospacers**. "
        "Prediction starts at the saved graph anchors. Reference sequences are used only for evaluation after prediction.", "",
        "| Method | Saved representative contigs | Matching representatives | Matching input records |",
        "|---|---:|---:|---:|",
        f"| Previous Huber | {previous.get('representative_contigs', '—')} | {previous.get('matched_representatives', '—')} | {previous.get('matched_starts', '—')} |",
        f"| Flow-guided | {r['representative_contigs']} | {r['matched_representatives']} | {r['matched_starts']} |", "",
        f"Compared with Huber: **{len(delta.get('gained_starts', []))} gained, {len(delta.get('lost_starts', []))} lost**, "
        f"{len(delta.get('retained_starts', []))} retained matching input records. Exact IDs are in `summary.json`.", "",
        f"{r['accepted_starts']}/{summary['input_starts']} input records have an accepted extension, covering "
        f"{r['accepted_target_sequences']} distinct candidate targets. Saved contigs: {r['representative_length_min']}–{r['representative_length_max']} bp. "
        f"{r['full_two_flank_starts']} records reach 750 added bases on both sides. Shorter prefixes remain in `extensions.tsv` even when they fail the 500-bp output minimum.", "",
        "Matches require one HSP with at least 90% identity and 50% coverage of the entire query. "
        "The alignment need not cover the target. Matching input records can share a region; representative contigs are not distinct phage counts. "
        "An unmatched sequence is not automatically proven false.", "",
        f"Unchanged target sequences: {e['target_matches']}/{e['target_count']} match at 90% identity and 85% coverage.", "",
        "## Implemented rule", "",
        "1. Start with the minimum observed flow along the saved anchor. A one-node anchor uses the smaller incident directional total.",
        "2. Probe nonbranching arms within the remaining flank budget. Use their mean multiplicity; masked source-repeat arms use the first edge count instead.",
        "3. For each continuation, subtract the other retained arms from remaining seed flow and cap at that arm's depth. Continue only if exactly one result is positive. Otherwise stop BEFORE the split.",
        "4. Never replenish flow at a merge or source repeat. Keep the existing source-revisit, repeat-bound and artifact barriers. Short tips and clearly dominant reconverging bubbles are local simplifications; tied bubbles remain unresolved.",
        "5. Save qualifying prefixes, up to 750 bp per side. Select representatives by worst-flank support before uncertainty, coverage variation and length. No reference score enters this choice.", "",
        "This is a seeded, conservative adaptation of excess-flow reasoning, not exhaustive flowtig enumeration or a correctness guarantee on noisy counts. "
        "No beam searches through unresolved forks. Low seed support, competing genuine routes and nonconserved observed counts can make it stop early. "
        "Support lost at a stopping edge is not attributed to the retained prefix.", "",
        "## Stopping reasons (left and right counted separately)", ""]
    lines.extend(f"- {reason}: {count}" for reason, count in sorted(r['stop_reasons'].items()))
    diagnosis = summary.get("lost_match_diagnosis")
    if diagnosis:
        lines.extend(["", f"All lost matching inputs stop below the output minimum: retained prefixes are "
            f"{diagnosis['retained_prefix_length_min']}–{diagnosis['retained_prefix_length_max']} bp "
            f"(median {diagnosis['retained_prefix_length_median']:g} bp). "
            "For Q000008, the right side stops after one base: incoming flow 52 and outgoing weights 3 and 49 "
            "leave two positive continuations. This is a refusal to choose, not proof that longer graph routes do not exist.", "",
            "**Decision:** the strict unique-flow rule loses too much sensitivity on this dataset and should not replace Huber."])
    control = summary.get("representative_selection_control")
    if control:
        lines.extend(["", "## Representative selection isolated from flow stopping", "",
            f"Reran Huber on all {control['input_records']} records using only the new representative ranking. "
            "Every per-input sequence and acceptance decision was unchanged; prior BLAST evidence was reused only after exact sequence verification. "
            f"The saved output changed from {control['old_representatives']} representatives / {control['old_matching_representatives']} matches "
            f"to {control['new_representatives']} representatives / {control['new_matching_representatives']} matches. "
            f"Matching input records remain {control['matching_input_records']}: no gains or losses of predicted paths. "
            "The optional `viral_contigs_huber_control.fasta` preserves these control representatives separately from flow predictions."])
    if "resource_measurements" in summary:
        m = summary["resource_measurements"]
        lines.extend(["", f"Prediction: {m['elapsed_seconds'] / 60:.2f} min, {m['peak_rss_gib']:.2f} GiB peak sampled RSS, "
                      f"{m['threads']} threads. Memory guard stop: {m['memory_guard_stop']}."])
    lines.extend(["", "## Files", "",
        "`protospacers.fasta` preserves all inputs. `viral_contigs_flow.fasta` contains accepted representatives. "
        "`extensions.tsv` keeps all outcomes, actual sequences, target positions, multiplicities and flow decisions; "
        "`summary.json` holds parameters, exact gains/losses, measurements and provenance. `run.log` contains stage progress.", "",
        "Flow decisions are `position:node:incoming:outgoing_weights:choice:remaining`, separated by semicolons. "
        "Position is the number of added bases before the decision; choices index the retained arms (−2 means multiple supported arms, −1 means none). "
        "Native target coordinates are zero-based/end-exclusive. BLAST coordinates refer to the canonical evaluation sequence; orientation is recorded.", "",
        "![Flow support](flow_support.png)", "", "![Multiplicities](multiplicities.png)", "",
        "Validation: native arithmetic and representative-selection tests plus controlled flow graphs. "
        "Each output sequence is checked against its connected graph walk and saved target. "
        "The original graph is read-only. This remains a single-dataset development test.", ""])
    (output / "REPORT.md").write_text("\n".join(lines))



def make_lookahead_plots(output, evidence, summary):
    """Plot frozen ranking performance, mathematical loss and actual multiplicities."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
    reps = sorted((r for r in evidence if r["relation_to_representative"] == "representative" and r["reported"] == "1"),
                  key=lambda r: int(r["representative_rank"]))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    ranks = np.arange(1, len(reps) + 1)
    axes[0].step(ranks, np.cumsum([r["reference_pass"] == "1" for r in reps]), where="post", lw=2.5, color="#245b91")
    axes[0].plot([0, 50], [0, 25], "--", color="#8b9198", label="50% benchmark")
    axes[0].scatter([50], [25], marker="*", s=140, color="#16806b", label="Goal: 25 of 50")
    axes[0].set(xlabel="Saved representatives in prediction rank order", ylabel="Cumulative reference matches", xlim=(0, 51), ylim=(0, 30), title="Reference checked only after ranking")
    axes[0].legend(frameon=False)
    ratio = np.linspace(.02, 4, 400)
    z = np.abs(np.log(ratio))
    loss = np.where(z <= 1, z*z/2, z-.5)
    axes[1].plot(ratio, loss, lw=2.3, label="Huber loss", color="#245b91")
    axes[1].plot(ratio, np.exp(-loss), lw=2.3, label="Lookahead weight: exp(−loss)", color="#b66925")
    axes[1].axvline(1, color="#a4acb4", lw=1)
    axes[1].set(xlabel="Smoothed multiplicity ratio: (depth + 0.5) / (baseline + 0.5)", ylabel="Loss or relative branch weight", title="Ranking quantities, not correctness probabilities")
    axes[1].legend(frameon=False)
    for ax in axes: ax.grid(alpha=.15)
    fig.savefig(output / "selection.png", dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    populations = {True: [], False: []}
    for row in reps:
        for side in ("left", "right"):
            counts = [float(x) for x in row[f"{side}_counts"].split(",") if x]
            positions = [int(x) for x in row[f"{side}_scored_positions"].split(",") if x]
            populations[row["reference_pass"] == "1"].extend(counts[p-1] for p in positions)
    high = max([v for values in populations.values() for v in values] or [2])
    bins = np.geomspace(1, max(2, high + 1), 35)
    for passed, color in ((False, "#949da6"), (True, "#16806b")):
        values = populations[passed]
        if values: axes[0].hist(values, bins=bins, weights=np.ones(len(values))/len(values), histtype="step", lw=2,
            color=color, label="Reference match" if passed else "No qualifying match")
    axes[0].set_xscale("log")
    axes[0].set(xlabel="External edge multiplicity (log scale)", ylabel="Fraction of positions per bin", title="Selected representatives")
    axes[0].legend(frameon=False)
    example = (reps or evidence)[0]
    summary["plot_example_query"] = example["query_id"]
    summary["plot_example_policy"] = "First prediction-ranked representative; no reference-dependent choice."
    seed_length = int(example["seed_end"]) - int(example["seed_begin"])
    for side, color in (("left", "#245b91"), ("right", "#b66925")):
        counts = np.array([float(x) for x in example[f"{side}_counts"].split(",") if x])
        positions = np.array([int(x) for x in example[f"{side}_scored_positions"].split(",") if x], dtype=int)
        padding = int(example["left_bp"]) - len(counts) if side == "left" else 0
        xs = -padding - np.arange(1, len(counts)+1) if side == "left" else seed_length+np.arange(len(counts))
        if len(counts): axes[1].plot(xs, counts, color="#c8cdd3", lw=.8)
        if len(positions): axes[1].plot(xs[positions-1], counts[positions-1], color=color, lw=1.5, label=side.capitalize())
    axes[1].axvspan(0, seed_length, color="#64845b", alpha=.2)
    axes[1].set(xlabel="Position relative to target start (bp)", ylabel="Edge multiplicity", title=f"{example['query_id']}: highest-ranked output", xlim=(-750, seed_length+750))
    axes[1].set_yscale("symlog", linthresh=1)
    axes[1].set_ylim(bottom=0)
    axes[1].legend(frameon=False)
    for ax in axes: ax.grid(alpha=.15)
    fig.savefig(output / "multiplicities.png", dpi=180)
    plt.close(fig)


def write_lookahead_report(output, summary):
    """Keep candidate availability, selected output and input-record changes distinct."""
    e = summary["evaluation"]
    r = e["methods"]["lookahead"]
    pool = e["candidate_pool"]
    delta = r.get("versus_previous_huber", {})
    lines = ["# Branch lookahead: SRR14017553", "",
        f"**Final output: {r['representative_contigs']} representative contigs; {r['matched_representatives']} match the reference.** "
        + ("The target of at least 25 matches within 50 outputs was reached." if r['goal_25_of_50_met'] else "The target of at least 25 matches within 50 outputs was not reached."), "",
        f"All **{summary['input_starts']} input records** were processed and all **{summary['input_target_sequences']} candidate protospacers** preserved. "
        f"These inputs contain {summary['distinct_oriented_anchor_paths']} distinct oriented anchors; multiple inputs can describe the same region. "
        "The existing graph was used directly, with at most 750 added bases per side and a 500-bp minimum saved contig.", "",
        "## Fixed ranking evaluated afterward", "", "| Top requested | Saved | Reference matches |", "|---:|---:|---:|"]
    lines.extend(f"| {v['requested']} | {v['saved']} | {v['matching']} |" for v in r["top_k"])
    lines.extend(["", f"Before the 50-output cap: {r['pre_cap_representatives']} representatives from {r['accepted_starts']} accepted input records. "
        f"Final lengths: {r['representative_length_min']}–{r['representative_length_max']} bp.", "",
        "## What the search retained", "",
        f"The compressed audit contains {pool['audit_records']} one-sided alternatives. Together with accepted per-input predictions, "
        f"this is {pool['unique_sequences']} distinct sequences, of which {pool['matching_unique_sequences']} match. "
        f"At least one matching retained sequence exists for {pool['matching_inputs']} input records. "
        f"Grouping only the matching sequences retrospectively gives {pool['retrospective_matching_sequence_groups']} sequence groups. Several groups can come from the same input or reference region; this is not that many independent recoveries. "
        "That last number is a diagnostic of this saved pool, not reference-free output and not a proof that all recoverable graph routes were found.", "",
        "## Change against Huber", "",
        f"The selected per-input sequences in the audit have {r['matched_starts']} matching input records: "
        f"{len(delta.get('gained_starts', []))} gained, {len(delta.get('lost_starts', []))} lost, "
        f"{len(delta.get('retained_starts', []))} retained versus the frozen Huber run. "
        "These are input-record counts before grouping and the output cap; they must not be reported as distinct matching contigs. Exact IDs are in `summary.json`.", "",
        f"The unchanged short targets have {e['target_matches']}/{e['target_count']} reference matches at 90% identity and 85% coverage.", "",
        "## Method", "",
        "1. Extend each saved anchor directly in the graph. Preserve competing families defined by the first two unresolved forks; use a bounded beam, local short-tip handling and clearly dominant reconverging bubbles.",
        "2. At a fork, inspect downstream unitig depths within two read spans, at most 256 unitig-entry edges and 512 graph operations per probe. Probe work is capped at one eighth of the query budget per flank; exhausted probes fall back to immediate counts. After calibration, compare those depths with the flank baseline using Huber log-ratio loss; before calibration, use the largest available downstream depth. Short low-count connections are not hard-pruned.",
        "3. Keep terminal routes and supported prefixes separately. Retain at most 16 of each per side, with family quotas. Select compatible left/right prefixes, allowing either side alone, by the worse flank score; length breaks score ties.",
        "4. Rank direct region representatives by graph-derived support and save at most 50. Save all input targets, per-input outcomes and a single compressed alternative-path audit.",
        "5. Freeze predictions, then BLAST against the local curated database. No reference score or hit status selects a path, flank length or representative.", "",
        "This graph-only coverage heuristic is inspired by looking beyond short connections in metagenomic assembly. It lacks read-pair linkage and is not equivalent to a complete assembler or a calibrated probability model. "
        "A stable multiplicity does not establish viral origin. Similar-coverage paths can remain indistinguishable.", "",
        f"Search limits affected {r['accepted_limited_starts']} accepted inputs; {r['accepted_ambiguous_starts']} have competing routes flagged. "
        "The retained pool is bounded, not exhaustive. Source-array nodes cannot be revisited; non-source nodes retain the configured visit bound. "
        "All saved sequences were checked against their connected graph walks.", "",
        f"The final matching representatives cover {len(r['best_hit_reference_loci'])} overlapping-best-hit reference interval groups; these are not distinct phage counts.", "",
        "Matches use one HSP with at least 90% identity and 50% of the full query covered, including its target. "
        "The HSP need not cover the target. Unmatched sequences are not automatically false; representative regions are not counts of distinct phages. This is a single-dataset development experiment.", "",
        "## Files", "",
        "`protospacers.fasta`: all 159 candidate targets. `viral_contigs_lookahead.fasta`: final selected representatives. "
        "`extensions.tsv`: all 164 input outcomes with sequences, counts, scores, ranks and reference evidence. "
        "`candidate_paths.tsv.gz`: retained one-sided alternatives with reference evaluation added afterward. "
        "`summary.json`: thresholds, input hashes, measurements and exact gains/losses. `run.log`: native progress.", "",
        "![Selection curve and coverage loss](selection.png)", "", "![Observed multiplicities](multiplicities.png)", ""])
    diagnosis = e.get("selection_diagnosis")
    if diagnosis:
        lines.extend(["## Where retained matches were lost", "",
            f"{len(diagnosis['previous_huber_matches_missing_from_retained_pool'])} previously matching Huber inputs have no qualifying match among this retained pool. "
            f"For {len(diagnosis['matching_pool_inputs_not_chosen'])} inputs, a retained alternative matched but the chosen per-input sequence did not. "
            f"Before the output cap, {diagnosis['matching_representatives_before_cap']} representative sequences matched; "
            f"{len(diagnosis['matching_representatives_below_cap'])} were ranked below the saved top 50. "
            f"{len(diagnosis['matching_input_paths_grouped_under_nonmatching_representative'])} matching input sequences were grouped under saved representatives that did not pass the same whole-query coverage criterion. "
            "These are distinct stages, and their counts should not be added together as distinct regions. Exact IDs are in the summary.", ""])
    if "resource_measurements" in summary:
        m = summary["resource_measurements"]
        lines.extend([f"Native prediction: {m['elapsed_seconds']/60:.2f} minutes, {m['peak_rss_gib']:.2f} GiB peak sampled RSS, {m['threads']} threads. "
                      f"Memory guard stop: {m['memory_guard_stop']}.", ""])
    (output / "REPORT.md").write_text("\n".join(lines))



def diagnose_lookahead_selection(output, evidence, summary):
    """Attribute saved-pool misses after evaluation, without changing any prediction."""
    e = summary["evaluation"]
    reps = {r["viral_contig_id"]: r for r in evidence if r["relation_to_representative"] == "representative"}
    matched = {r["query_id"] for r in evidence if r["accepted"] == "1" and r["reference_pass"] == "1"}
    available = set(e["candidate_pool"]["matching_input_ids"])
    old = set(e.get("previous_huber", {}).get("matching_input_ids", []))
    e["selection_diagnosis"] = dict(
        previous_huber_matches_missing_from_retained_pool=sorted(old - available),
        previous_huber_matches_present_in_retained_pool=sorted(old & available),
        new_matching_inputs_available_in_pool=sorted(available - old),
        matching_pool_inputs_not_chosen=sorted(available - matched),
        matching_representatives_before_cap=sum(r["reference_pass"] == "1" for r in reps.values()),
        matching_representatives_below_cap=[key for key, r in reps.items() if r["reference_pass"] == "1" and r["reported"] == "0"],
        matching_input_paths_excluded_by_cap=[r["query_id"] for r in evidence if r["query_id"] in matched and r["reported"] == "0"],
        matching_input_paths_grouped_under_nonmatching_representative=[r["query_id"] for r in evidence
            if r["query_id"] in matched and r["reported"] == "1" and reps[r["viral_contig_id"]]["reference_pass"] == "0"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--baseline", type=Path, required=True, help="Previous prediction directory")
    parser.add_argument("--baseline-evaluation", type=Path, required=True, help="Previous BLAST evaluation directory")
    parser.add_argument("--previous-seeded", type=Path, help="Frozen seeded run to compare against (Huber)")
    parser.add_argument("--reference", type=Path, required=True, help="Local BLAST database prefix")
    parser.add_argument("--blastn", default="blastn")
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    output = args.output.resolve()
    summary = json.loads((output / "summary.json").read_text())
    assert summary["status"] == "complete", "Evaluate only a completed run"
    evidence = rows(output / "extensions.tsv")
    methods = tuple(m["method"] for m in summary["methods"])
    assert len(evidence) == len(methods) * summary["input_starts"]
    assert set(row["method"] for row in evidence) == set(methods)
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
    audit = []
    if methods == ("lookahead",):
        summary.setdefault("frozen_prediction_sha256", {name: digest(output / name) for name in
            ("extensions.tsv", "candidate_paths.tsv.gz", "viral_contigs_lookahead.fasta", "protospacers.fasta")})
        with gzip.open(output / "candidate_paths.tsv.gz", "rt") as stream:
            audit = list(csv.DictReader(stream, delimiter="\t"))
        audit.sort(key=lambda row: (row["query_id"], row["side"], int(row["index"])))
        for row in audit:
            seq = row["sequence"]
            b, e = int(row["seed_begin"]), int(row["seed_end"])
            assert canonical(seq[b:e]) == targets[row["target_id"]]
            assert 0 <= b <= 750 and 0 <= len(seq) - e <= 750 and len(seq) >= summary["min_contig_length"]
            sequence = canonical(seq)
            if sequence not in sequence_ids:
                key = f"S{len(unique) + 1:06d}"
                unique[key] = sequence
                sequence_ids[sequence] = key
    print(f"Evaluating {len(unique)} distinct retained extension sequences; paths already fixed.", flush=True)
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
    for method in methods:
        own = [r for r in evidence if r["method"] == method]
        assert len({r["query_id"] for r in own}) == summary["input_starts"]
        accepted = [r for r in own if r["accepted"] == "1"]
        matched_queries = {r["query_id"] for r in accepted if r["reference_pass"] == "1"}
        reps = fasta(output / f"viral_contigs_{method}.fasta")
        rep_rows = [r for r in own if r["relation_to_representative"] == "representative"]
        all_reps = {r["viral_contig_id"]: r["sequence"] for r in rep_rows}
        reported_rows = [r for r in rep_rows if r.get("reported", "1") == "1"]
        assert reps == {r["viral_contig_id"]: r["sequence"] for r in reported_rows}
        assert len(set(map(canonical, reps.values()))) == len(reps), "Duplicate representative sequences"
        assert all(same_region(all_reps[r["viral_contig_id"]], r["sequence"], summary["graph_k"]) for r in accepted), "Invalid native representative association"
        rep_hits = [hits[sequence_ids[canonical(seq)]] for seq in reps.values()]
        passing_reps = [hit for hit in rep_hits if hit["passes"]]
        lengths = [len(seq) for seq in reps.values()]
        evaluation["methods"][method] = dict(
            pre_cap_representatives=len(all_reps), reported_starts=sum(r.get("reported", "1") == "1" for r in accepted),
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
    if args.previous_seeded:
        previous = json.loads((args.previous_seeded / "summary.json").read_text())
        previous_rows = [r for r in rows(args.previous_seeded / "extensions.tsv") if r["method"] == "huber"]
        assert {r["query_id"] for r in previous_rows} == tested_queries, "Different starting population"
        assert fasta(args.previous_seeded / "protospacers.fasta") == targets, "Different target population"
        assert previous["evaluation"]["identity"] == .90 and previous["evaluation"]["extension_coverage"] == .50
        for name, expected in previous.get("artifact_sha256", {}).items():
            assert digest(args.previous_seeded / name) == expected, "Previous evidence changed"
        old = {r["query_id"] for r in previous_rows if r["accepted"] == "1" and r["reference_pass"] == "1"}
        evaluation["previous_huber"] = dict(directory=str(args.previous_seeded.resolve()),
            matched_starts=len(old), matching_input_ids=sorted(old), representative_contigs=previous["evaluation"]["methods"]["huber"]["representative_contigs"],
            matched_representatives=previous["evaluation"]["methods"]["huber"]["matched_representatives"])
        for method in methods:
            new = {r["query_id"] for r in evidence if r["method"] == method and r["accepted"] == "1" and r["reference_pass"] == "1"}
            evaluation["methods"][method]["versus_previous_huber"] = dict(
                gained_starts=sorted(new - old), lost_starts=sorted(old - new), retained_starts=sorted(new & old))
    if methods == ("lookahead",):
        own = evaluation["methods"]["lookahead"]
        rep_rows = sorted((r for r in evidence if r["relation_to_representative"] == "representative"),
                          key=lambda r: int(r["representative_rank"]))
        own["top_k"] = [dict(requested=n, saved=min(n, len(rep_rows)),
            matching=sum(r["reference_pass"] == "1" for r in rep_rows[:n])) for n in (10, 20, 30, 40, 50)]
        own["goal_25_of_50_met"] = own["representative_contigs"] <= 50 and own["matched_representatives"] >= 25
        accepted = [r for r in evidence if r["accepted"] == "1"]
        pool_queries = {r["query_id"] for r in audit if hits[sequence_ids[canonical(r["sequence"]) ]]["passes"]}
        pool_queries.update(r["query_id"] for r in accepted if r["reference_pass"] == "1")
        matching_unique = [seq for key, seq in unique.items() if hits[key]["passes"]]
        # Availability diagnostic only: matching paths are grouped retrospectively.
        # This is not a predictor and not proof of a recovery ceiling in the graph.
        groups = []
        for seq in sorted(matching_unique, key=lambda seq: (-len(seq), seq)):
            if not any(same_region(rep, seq, summary["graph_k"]) for rep in groups): groups.append(seq)
        evaluation["candidate_pool"] = dict(audit_records=len(audit), unique_sequences=len(unique),
            matching_unique_sequences=len(matching_unique), matching_inputs=len(pool_queries),
            matching_input_ids=sorted(pool_queries), retrospective_matching_sequence_groups=len(groups),
            retrospective_grouping_policy="Matching sequences only, longest-first direct similarity grouping; variants from one input or reference region can form multiple groups. Diagnostic, not independent recoveries or selection.",
            exhaustive=False)
        summary["frozen_pool_blast_input_sha256"] = hashlib.sha256("".join(f">{key}\n{seq}\n" for key, seq in unique.items()).encode()).hexdigest()
        if audit:
            for row in audit:
                hit = hits[sequence_ids[canonical(row["sequence"])]]
                row.update(reference_pass=int(hit["passes"]), reference_record=hit["reference"],
                    reference_identity=hit["identity"], reference_coverage=hit["coverage"],
                    reference_query_begin=hit["query_begin"], reference_query_end=hit["query_end"],
                    reference_begin=hit["reference_begin"], reference_end=hit["reference_end"],
                    reference_query_orientation="+" if row["sequence"] == canonical(row["sequence"]) else "-")
            with gzip.open(output / "candidate_paths.tsv.gz", "wt") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(audit[0]), delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(audit)
    summary["evaluation"] = evaluation
    if methods == ("lookahead",):
        diagnose_lookahead_selection(output, evidence, summary)
    summary["input_and_prediction_hashes"] = {str(p.resolve()): digest(p) for p in
        [seed_dir / "queries.tsv", seed_dir / "sources.tsv", seed_dir / "graph_paths.tsv", seed_dir / "candidates.tsv",
         output / "protospacers.fasta", *(output / f"viral_contigs_{m}.fasta" for m in methods)]}
    if methods == ("lookahead",):
        make_lookahead_plots(output, evidence, summary)
    elif methods == ("flow",):
        make_flow_plots(output, evidence, summary)
    else:
        make_plots(output, evidence, summary)
    with (output / "extensions.tsv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(evidence[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(evidence)
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    if methods == ("lookahead",):
        write_lookahead_report(output, summary)
    elif methods == ("flow",):
        write_flow_report(output, summary)
    else:
        write_report(output, summary)
    print(json.dumps({"baseline": evaluation["baseline"], "target_matches": evaluation["target_matches"],
                      "methods": evaluation["methods"]}, indent=2))


if __name__ == "__main__":
    main()
