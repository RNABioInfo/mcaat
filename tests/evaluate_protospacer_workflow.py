#!/usr/bin/env python3
"""Audit frozen workflow outputs; reference searches never feed back into prediction."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import shutil
import statistics
import subprocess
import time


def fasta(path):
    records, name = {}, None
    with Path(path).open() as stream:
        for line in stream:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in records:
                    raise ValueError(f"duplicate FASTA ID: {name}")
                records[name] = ""
            elif line.strip():
                if name is None or set(line.strip()) - set("ACGT"):
                    raise ValueError(f"invalid DNA in {path}")
                records[name] += line.strip()
    if any(not seq for seq in records.values()):
        raise ValueError(f"empty sequence in {path}")
    return records


def table(path):
    with Path(path).open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if any(None in row or None in row.values() for row in rows):
        raise ValueError(f"inconsistent TSV columns in {path}")
    return rows


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def validate(folder):
    """Check all source/query outcomes, sequence identities, and association coordinates."""
    summary = json.loads((folder / "run_summary.json").read_text())
    if summary["implementation"] not in ("protospacer-workflow-v1", "protospacer-workflow-v2") or summary["status"] != "complete":
        raise ValueError("evaluation requires a complete workflow run")
    sources = {r["source_id"]: r for r in table(folder / "sources.tsv")}
    queries = {r["query_id"]: r for r in table(folder / "queries.tsv")}
    targets, extensions = fasta(folder / "Protospacers.fasta"), fasta(folder / "ContiguousPaths.fasta")
    links = table(folder / "associations.tsv")
    candidates = table(folder / "candidates.tsv")
    accepted = {r["query_id"]: r for r in candidates if r["status"] in ("detected", "extended")}
    contexts = table(folder / "target_context.tsv")
    paths = table(folder / "graph_paths.tsv")
    evidence = table(folder / "coverage_evidence.tsv")
    assert len(sources) == summary["spacer_observations"]
    assert len(queries) == summary["distinct_spacers"]
    assert sum(int(q["source_count"]) for q in queries.values()) == len(sources)
    assert all(q["status"] not in ("not_processed", "interrupted") for q in queries.values())
    assert len(targets) == summary["protospacer_regions"]
    assert len(extensions) == summary["unique_extensions"]
    assert len(accepted) == summary["detected_queries"] == len(contexts) == len(paths)
    assert len({r["query_id"] for r in links if r["extension_id"]}) == summary["extended_queries"]
    for sequences in (targets, extensions):
        assert len({min(s, reverse_complement(s)) for s in sequences.values()}) == len(sequences)
    for r in links:
        source, query = sources[r["source_id"]], queries[r["query_id"]]
        candidate = accepted[r["query_id"]]
        aligned = source["original_spacer"][int(r["query_begin"]):int(r["query_end"])]
        if r["extension_id"]:
            seq = extensions[r["extension_id"]]
            begin, end = int(r["seed_begin"]), int(r["seed_end"])
            fragment_begin, fragment_end = int(r["supported_fragment_begin"]), int(r["supported_fragment_end"])
            assert 0 <= fragment_begin <= begin < end <= fragment_end <= len(seq)
            target = seq[begin:end]
            assert min(target, reverse_complement(target)) == targets[r["target_id"]]
            assert fragment_end-fragment_begin == len(target)+int(r["left_bp"])+int(r["right_bp"])
        else:
            target = targets[r["target_id"]]
            assert not any(r[k] for k in ("seed_begin", "seed_end", "supported_fragment_begin", "supported_fragment_end"))
        oriented = reverse_complement(target) if r["target_strand"] == "-" else target
        assert len(aligned) == len(oriented)
        identity = sum(a == b for a, b in zip(aligned, oriented)) / len(aligned)
        coverage = len(aligned) / len(source["spacer"])
        assert abs(identity-float(r["identity"])) < 1e-5
        assert abs(coverage-float(r["query_coverage"])) < 1e-5
        assert identity+1e-9 >= summary["min_identity"] and coverage+1e-9 >= summary["min_query_coverage"]
        assert candidate["target_id"] == r["target_id"]
    for r in contexts:
        # These are the selected flanks before the extension length filter;
        # a detected-but-unextended target must still retain this evidence.
        assert any(len(r[side+"_outward"]) >= max(summary["graph_k"], int(r[side+"_clear_at"]))+summary["graph_k"]+1
                   for side in ("left", "right"))
    for r in paths:
        c = accepted[r["query_id"]]
        assert len(r["anchor_nodes"].split(",")) == len(c["sequence"])-summary["graph_k"]+1
        for side in ("left", "right"):
            nodes = r[side+"_nodes_outward"].split(",") if r[side+"_nodes_outward"] else []
            assert len(nodes) in (int(c[side+"_bp"])+1, 0)
            assert int(c[side+"_bp"]) <= summary["max_flank_length"]
    lengths = sorted(map(len, extensions.values()))
    return {"validation": "passed", "source_observations": len(sources), "distinct_corrected_spacers": len(queries),
            "detected_queries": len(accepted), "protospacer_regions": len(targets), "unique_extensions": len(extensions),
            "extension_lengths": {"minimum": min(lengths, default=0), "median": statistics.median(lengths) if lengths else 0,
                                  "maximum": max(lengths, default=0)},
            "query_outcomes": dict(Counter(q["status"] for q in queries.values())),
            "accepted_search_limited_queries": sum(c["search_limited"] == "1" for c in accepted.values()),
            "all_search_limited_queries": sum(q.get("search_limited", "0") == "1" for q in queries.values()),
            "left_stops": dict(Counter(c["left_stop"] for c in accepted.values())),
            "right_stops": dict(Counter(c["right_stop"] for c in accepted.values())),
            "coverage_observations": len(evidence)}


def blast(fasta_path, database, destination, threads):
    """Search frozen predictions with the same BLAST settings for every comparison."""
    destination.mkdir(parents=True, exist_ok=True)
    output = destination / "hits.tsv"
    command = [shutil.which("blastn") or "/usr/local/bin/blastn", "-query", str(fasta_path), "-db", str(database),
               "-task", "blastn-short", "-word_size", "7", "-dust", "no", "-soft_masking", "false",
               "-evalue", "1000", "-perc_identity", "90", "-max_target_seqs", "1000000", "-num_threads", str(threads),
               "-outfmt", "6 qseqid sseqid qlen qstart qend sstart send length nident pident evalue bitscore", "-out", str(output)]
    write_json(destination / "command.json", command)
    started = time.monotonic()
    if fasta(fasta_path):
        with (destination / "blast.log").open("w") as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    else:
        output.write_text("")
    write_json(destination / "timing.json", {"seconds": time.monotonic()-started})
    return output


def summarize_matches(sequences, hits, minimum_coverage, destination):
    """Use exact nident/alignment-length and query coordinates, not rounded BLAST percentages."""
    best, passing, rows = {}, {}, 0
    with hits.open() as stream:
        for line in stream:
            row = line.rstrip().split("\t")
            qid, reference = row[:2]
            qlen, qb, qe, sb, se, length, nident = map(int, row[2:9])
            assert qid in sequences and qlen == len(sequences[qid])
            identity = nident/length
            coverage = (abs(qe-qb)+1)/qlen
            rows += 1
            if identity+1e-12 < .9:
                continue
            value = dict(sequence_id=qid, reference=reference, query_length=qlen, identity=identity,
                         coverage=coverage, query_begin=qb, query_end=qe, reference_begin=sb, reference_end=se,
                         bitscore=float(row[11]))
            if qid not in best or (coverage, identity, value["bitscore"]) > (best[qid]["coverage"], best[qid]["identity"], best[qid]["bitscore"]):
                best[qid] = value
            if coverage+1e-12 >= minimum_coverage:
                passing[qid] = True
    fields = ["sequence_id", "query_length", "passes_threshold", "reference", "identity", "coverage", "query_begin", "query_end", "reference_begin", "reference_end", "bitscore"]
    with destination.open("w") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t")
        writer.writeheader()
        for qid, sequence in sequences.items():
            row = best.get(qid, {"sequence_id": qid, "query_length": len(sequence)})
            row["passes_threshold"] = int(qid in passing)
            writer.writerow(row)
    return {"sequences": len(sequences), "matching_sequences": len(passing),
            "reference_match_fraction": len(passing)/len(sequences) if sequences else None,
            "identity_threshold": .9, "single_hsp_query_coverage_threshold": minimum_coverage,
            "alignment_rows": rows, "matching_sequence_ids": sorted(passing)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--srr", required=True, choices=["SRR14017553", "SRR14017554", "SRR14017569"])
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--database", type=Path)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--baseline-hits", type=Path)
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        parser.error("evaluation output must be new or empty")
    args.output.mkdir(parents=True, exist_ok=True)
    report = {"dataset": args.srr, "reference_used_for_selection": False,
              "validation": validate(args.results), "run_summary": json.loads((args.results/"run_summary.json").read_text()),
              "input_sha256": {p.name: digest(p) for p in sorted(args.results.iterdir()) if p.is_file()}}
    write_json(args.output / "evaluation.json", report)
    if args.database:
        dbcmd = shutil.which("blastdbcmd") or "/usr/local/bin/blastdbcmd"
        (args.output/"reference_info.txt").write_text(subprocess.check_output([dbcmd, "-db", str(args.database), "-info"], text=True))
        for name, coverage in (("Protospacers", .85), ("ContiguousPaths", .50)):
            path = args.results / (name+".fasta")
            hitfile = blast(path, args.database, args.output/name, args.threads)
            report[name] = summarize_matches(fasta(path), hitfile, coverage, args.output/(name+"_coverage.tsv"))
            write_json(args.output / "evaluation.json", report)
        queries = table(args.results/"queries.tsv")
        all_inputs = args.output/"all_input_spacers.fasta"
        with all_inputs.open("w") as stream:
            for row in queries:
                stream.write(">"+row["query_id"]+"\n"+row["spacer"]+"\n")
        hits = blast(all_inputs, args.database, args.output/"all_input_spacers", args.threads)
        inputs = summarize_matches(fasta(all_inputs), hits, .85, args.output/"all_input_spacers_coverage.tsv")
        known = set(inputs["matching_sequence_ids"])
        candidates = table(args.results/"candidates.tsv")
        accepted = {r["query_id"] for r in candidates if r["status"] in ("detected", "extended")}
        inputs["reference_matching_queries_detected"] = len(known & accepted)
        inputs["reference_matching_query_outcomes"] = dict(Counter(q["status"] for q in queries if q["query_id"] in known))
        report["input_spacers"] = inputs
    if args.baseline:
        if not args.baseline_hits and not args.database:
            parser.error("baseline needs existing hits or a database")
        hits = args.baseline_hits or blast(args.baseline, args.database, args.output/"baseline", args.threads)
        report["baseline_intermediate_discovery"] = summarize_matches(fasta(args.baseline), hits, .50, args.output/"baseline_coverage.tsv")
        report["baseline_intermediate_discovery"]["warning"] = "Intermediate discovery, not old final filtered output; not a like-for-like algorithm comparison."
        report["baseline_sha256"] = digest(args.baseline)
    for name, before in report["input_sha256"].items():
        assert digest(args.results/name) == before, "prediction changed during evaluation"
    write_json(args.output / "evaluation.json", report)
    print(json.dumps({k: v for k, v in report.items() if k not in ("input_sha256", "run_summary")}, indent=2))


if __name__ == "__main__":
    main()
