#!/usr/bin/env python3
"""Score saved extensions by whether the reference match is around the target.

The usual criterion counts a saved contig as a match when one BLAST HSP has
>= 80% identity and covers >= 50% of the contig. That HSP can lie entirely in a
flank that is unrelated to the target, for example sequence reached through a
CRISPR repeat and a neighbouring spacer. This script keeps that criterion and
adds a same-region criterion:

  spans_target                 a qualifying contig HSP covers >= 90% of the
                               target interval inside the contig, or
  target_at_predicted_position the target itself aligns (>= 80% identity,
                               >= 85% of its length) to the same reference
                               record, same strand, within 50 bp of where the
                               contig HSP places it.

Other matching contigs are classified by where (if anywhere) their target
lands on the reference records the contigs match. The reference is used only
after prediction, on frozen outputs.

Layout: ROOT/SRR*/<stage folder>/ with either
  discovery: run_summary.json, ContiguousPaths.fasta, associations.tsv,
             CandidateProtospacers.fasta
  seeded:    summary.json, viral_contigs_evidence.fasta, extensions.tsv,
             protospacers.fasta

Example:
  python3 evaluate_target_region.py --root /vol/d/data/copsac \
      --db /path/all_curated_db --out region_eval \
      --stage old_seed=evidence19_seeds --stage old_support=evidence19_support_0.998 \
      --stage new_seed=nocross_seeds --stage new_support=nocross_support_0.998
"""
import argparse
import collections
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

COMPLEMENT = str.maketrans('ACGT', 'TGCA')
CLASSES = ['spans_target', 'target_at_predicted_position', 'same_record_nearby',
           'same_record_far_or_opposite_strand', 'target_on_other_record', 'target_not_found']
SAME_REGION = {'spans_target', 'target_at_predicted_position'}


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def canonical(sequence):
    return min(sequence, reverse_complement(sequence))


def read_fasta(path):
    records, name = {}, None
    with open(path) as stream:
        for line in stream:
            line = line.strip()
            if line.startswith('>'):
                name = line[1:].split()[0]
                records[name] = []
            elif name is not None:
                records[name].append(line)
    return {k: ''.join(v).upper() for k, v in records.items()}


def read_tsv(path):
    with open(path, newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def load_stage(folder):
    """Return {contig_id: (sequence, [(target_id, begin, end), ...])} for one completed stage."""
    if (folder / 'run_summary.json').is_file():
        summary = json.loads((folder / 'run_summary.json').read_text())
        if summary.get('status') != 'complete':
            return None
        contigs = read_fasta(folder / 'ContiguousPaths.fasta')
        targets = read_fasta(folder / 'CandidateProtospacers.fasta')
        intervals = collections.defaultdict(set)
        for row in read_tsv(folder / 'associations.tsv'):
            if row['extension_id'] and row['seed_begin'] != '':
                intervals[row['extension_id']].add((row['target_id'], int(row['seed_begin']), int(row['seed_end'])))
    elif (folder / 'summary.json').is_file():
        summary = json.loads((folder / 'summary.json').read_text())
        if summary.get('status') != 'complete':
            return None
        contigs = read_fasta(folder / 'viral_contigs_evidence.fasta')
        targets = read_fasta(folder / 'protospacers.fasta')
        intervals = collections.defaultdict(set)
        for row in read_tsv(folder / 'extensions.tsv'):
            if row.get('reported') == '1' and row.get('relation_to_representative') == 'representative':
                intervals[row['viral_contig_id']].add((row['target_id'], int(row['seed_begin']), int(row['seed_end'])))
    else:
        return None
    out = {}
    for contig, sequence in contigs.items():
        checked = []
        for target, begin, end in sorted(intervals.get(contig, ())):
            if not 0 <= begin < end <= len(sequence):
                raise ValueError(f'{folder}: {contig} target {target} outside contig')
            if target in targets and canonical(sequence[begin:end]) != canonical(targets[target]):
                raise ValueError(f'{folder}: {contig} target {target} sequence differs at its coordinates')
            checked.append((target, begin, end))
        out[contig] = (sequence, checked)
    return out


def tool(name, directory):
    found = shutil.which(name, path=directory) if directory else shutil.which(name)
    if not found:
        sys.exit(f'cannot find {name}; pass --blast-bin')
    return found


def run(command):
    subprocess.run([str(c) for c in command], check=True)


def blast_table(path):
    hits = collections.defaultdict(list)
    if not path.is_file():
        return hits
    with open(path) as stream:
        for line in stream:
            f = line.rstrip('\n').split('\t')
            # qseqid sseqid length nident qstart qend sstart send qlen
            hits[f[0]].append((f[1], int(f[2]), int(f[3]), int(f[4]), int(f[5]), int(f[6]), int(f[7]), int(f[8])))
    return hits


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--root', type=Path, required=True, help='folder holding SRR*/ sample folders')
    parser.add_argument('--stage', action='append', required=True, metavar='LABEL=FOLDER',
                        help='stage label and folder name inside each sample; repeat for each stage')
    parser.add_argument('--db', required=True, help='BLAST nucleotide database prefix (reference catalogue)')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--blast-bin', default=None, help='folder with blastn, blastdbcmd, makeblastdb')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--samples', nargs='*', default=None, help='restrict to these sample names')
    parser.add_argument('--include-incomplete', action='store_true',
                        help='also score samples where some stages are missing (default: all stages complete)')
    parser.add_argument('--identity', type=float, default=0.80)
    parser.add_argument('--contig-coverage', type=float, default=0.50)
    parser.add_argument('--target-coverage', type=float, default=0.85)
    parser.add_argument('--target-span', type=float, default=0.90)
    parser.add_argument('--position-tolerance', type=int, default=50)
    parser.add_argument('--nearby', type=int, default=5000)
    parser.add_argument('--paper-blast', action='store_true',
                        help='contig search with blastn-short, word 9, E 1000 (manuscript settings; slow)')
    args = parser.parse_args()

    stages = []
    for item in args.stage:
        label, _, folder = item.partition('=')
        if not label or not folder:
            sys.exit(f'--stage needs LABEL=FOLDER, got {item}')
        stages.append((label, folder))
    blastn, blastdbcmd, makeblastdb = (tool(n, args.blast_bin) for n in ('blastn', 'blastdbcmd', 'makeblastdb'))
    args.out.mkdir(parents=True, exist_ok=True)

    # Load every stage of every sample; keep samples where all stages completed.
    data, skipped = {}, {}
    for sample in sorted(p for p in args.root.iterdir() if p.is_dir() and p.name.startswith('SRR')):
        if args.samples and sample.name not in args.samples:
            continue
        loaded = {label: load_stage(sample / folder) for label, folder in stages}
        missing = [label for label, value in loaded.items() if value is None]
        if missing and not args.include_incomplete:
            skipped[sample.name] = missing
            continue
        data[sample.name] = {k: v for k, v in loaded.items() if v is not None}
    if not data:
        sys.exit('no sample has all requested stages complete')
    print(f'{len(data)} samples scored; {len(skipped)} skipped for incomplete stages', flush=True)

    # One BLAST per distinct contig sequence.
    sequences = {}
    for stages_of_sample in data.values():
        for contigs in stages_of_sample.values():
            for sequence, _ in contigs.values():
                sequences.setdefault(sequence, f'C{len(sequences) + 1:07d}')
    contig_fasta = args.out / 'contigs.fa'
    with open(contig_fasta, 'w') as stream:
        for sequence, key in sequences.items():
            stream.write(f'>{key}\n{sequence}\n')
    fields = '6 qseqid sseqid length nident qstart qend sstart send qlen'
    contig_hits_path = args.out / 'contig_hits.tsv'
    search = ['-task', 'blastn-short', '-word_size', '9', '-evalue', '1000'] if args.paper_blast else \
             ['-task', 'blastn', '-evalue', '1e-10']
    print(f'BLAST {len(sequences)} distinct contig sequences', flush=True)
    run([blastn, '-query', contig_fasta, '-db', args.db, *search, '-dust', 'no', '-soft_masking', 'false',
         '-perc_identity', str(int(args.identity * 100)), '-qcov_hsp_perc', str(int(args.contig_coverage * 100)),
         '-max_target_seqs', '500', '-num_threads', args.threads, '-outfmt', fields, '-out', contig_hits_path])
    contig_hits = {}
    for key, rows in blast_table(contig_hits_path).items():
        passing = [r for r in rows if r[2] / r[1] >= args.identity and (r[4] - r[3] + 1) / r[7] >= args.contig_coverage]
        if passing:
            contig_hits[key] = passing

    # Target intervals of matching contigs, searched against the records those contigs match.
    targets = {}
    for stages_of_sample in data.values():
        for contigs in stages_of_sample.values():
            for sequence, intervals in contigs.values():
                if sequences[sequence] in contig_hits:
                    for _, begin, end in intervals:
                        targets.setdefault(sequence[begin:end], f'T{len(targets) + 1:07d}')
    subjects = sorted({r[0] for rows in contig_hits.values() for r in rows})
    target_hits = {}
    if subjects and targets:
        (args.out / 'subjects.txt').write_text('\n'.join(subjects) + '\n')
        with open(args.out / 'subjects.fa', 'w') as stream:
            subprocess.run([blastdbcmd, '-db', args.db, '-entry_batch', str(args.out / 'subjects.txt')],
                           check=True, stdout=stream)
        run([makeblastdb, '-in', args.out / 'subjects.fa', '-dbtype', 'nucl', '-parse_seqids',
             '-out', args.out / 'subjects'])
        with open(args.out / 'targets.fa', 'w') as stream:
            for sequence, key in targets.items():
                stream.write(f'>{key}\n{sequence}\n')
        print(f'BLAST {len(targets)} distinct target intervals against {len(subjects)} matched records', flush=True)
        run([blastn, '-query', args.out / 'targets.fa', '-db', args.out / 'subjects', '-task', 'blastn-short',
             '-word_size', '9', '-dust', 'no', '-soft_masking', 'false', '-evalue', '1000',
             '-perc_identity', str(int(args.identity * 100)), '-qcov_hsp_perc', str(int(args.target_coverage * 100)),
             '-max_target_seqs', '100000', '-max_hsps', '1000', '-num_threads', args.threads,
             '-outfmt', fields, '-out', args.out / 'target_hits.tsv'])
        for key, rows in blast_table(args.out / 'target_hits.tsv').items():
            target_hits[key] = [r for r in rows if r[2] / r[1] >= args.identity and
                                (r[4] - r[3] + 1) / r[7] >= args.target_coverage]

    def classify(sequence, intervals):
        hsps = contig_hits.get(sequences[sequence])
        if not hsps:
            return 'no_match', {}
        if not intervals:
            return 'no_target_interval', {}
        best, detail = len(CLASSES), {}
        for target, begin, end in intervals:
            for subject, length, nident, qs, qe, ss, se, qlen in hsps:
                if max(0, min(qe, end) - max(qs - 1, begin)) / (end - begin) >= args.target_span and best > 0:
                    best, detail = 0, dict(target=target, reference=subject)
            found = target_hits.get(targets.get(sequence[begin:end]), [])
            if found and best > 4:
                best, detail = 4, dict(target=target, reference=found[0][0])
            for subject, length, nident, qs, qe, ss, se, qlen in hsps:
                strand = 1 if se >= ss else -1
                predicted = ss + (begin + 1 - qs) * strand
                for t_subject, _, _, tqs, _, tss, tse, _ in found:
                    if t_subject != subject:
                        continue
                    t_strand = 1 if tse >= tss else -1
                    offset = abs(tss - (tqs - 1) * t_strand - predicted)
                    level = 3 if t_strand != strand or offset > args.nearby else \
                        1 if offset <= args.position_tolerance else 2
                    if level < best:
                        best, detail = level, dict(target=target, reference=subject, offset=offset)
        return (CLASSES[best] if best < len(CLASSES) else CLASSES[-1]), detail

    records = []
    for sample, stages_of_sample in data.items():
        for label, contigs in stages_of_sample.items():
            for contig, (sequence, intervals) in contigs.items():
                category, detail = classify(sequence, intervals)
                records.append(dict(sample=sample, stage=label, contig=contig, length=len(sequence),
                                    targets=';'.join(f'{t}:{b}-{e}' for t, b, e in intervals),
                                    contig_match=int(category != 'no_match'),
                                    same_region=int(category in SAME_REGION), category=category,
                                    reference=detail.get('reference', ''), offset=detail.get('offset', ''),
                                    sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest()))
    with open(args.out / 'records.tsv', 'w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)

    labels = [label for label, _ in stages]
    per_sample = collections.defaultdict(lambda: collections.Counter())
    pooled = {label: collections.Counter() for label in labels}
    for r in records:
        for counter in (per_sample[(r['sample'], r['stage'])], pooled[r['stage']]):
            counter['contigs'] += 1
            counter['contig_match'] += r['contig_match']
            counter['same_region'] += r['same_region']
            counter[r['category']] += 1
    with open(args.out / 'samples.tsv', 'w', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerow(['sample'] + [f'{label}_{m}' for label in labels for m in ('contigs', 'contig_match', 'same_region')])
        for sample in data:
            writer.writerow([sample] + [per_sample[(sample, label)][m] for label in labels
                                        for m in ('contigs', 'contig_match', 'same_region')])
    lines = [f'# Same-region evaluation: {len(data)} samples', '',
             f'Contig match: one HSP >= {args.identity:.0%} identity covering >= {args.contig_coverage:.0%} of the contig. '
             f'Same region: that HSP covers >= {args.target_span:.0%} of the target, or the target aligns '
             f'(>= {args.identity:.0%} identity, >= {args.target_coverage:.0%} of its length) to the same record, same strand, '
             f'within {args.position_tolerance} bp of the position the contig HSP predicts. Counts are FASTA records.', '',
             '| Stage | Contigs | Contig matches | Same-region matches | Same-region / contigs |',
             '|---|---:|---:|---:|---:|']
    for label in labels:
        c = pooled[label]
        lines.append(f"| {label} | {c['contigs']:,} | {c['contig_match']:,} | {c['same_region']:,} | "
                     f"{100 * c['same_region'] / max(1, c['contigs']):.2f}% |")
    lines += ['', '| Stage | ' + ' | '.join(CLASSES) + ' |', '|---|' + '---:|' * len(CLASSES)]
    for label in labels:
        lines.append(f'| {label} | ' + ' | '.join(f'{pooled[label][c]:,}' for c in CLASSES) + ' |')
    if skipped:
        lines += ['', 'Skipped (incomplete stages): ' + ', '.join(f'{s} ({",".join(m)})' for s, m in skipped.items())]
    (args.out / 'REPORT.md').write_text('\n'.join(lines) + '\n')
    (args.out / 'summary.json').write_text(json.dumps(dict(
        samples=sorted(data), skipped=skipped, stages=dict(stages), pooled={k: dict(v) for k, v in pooled.items()},
        settings=dict(identity=args.identity, contig_coverage=args.contig_coverage, target_coverage=args.target_coverage,
                      target_span=args.target_span, position_tolerance=args.position_tolerance, nearby=args.nearby,
                      contig_search=search, database=args.db)), indent=2) + '\n')
    print('\n'.join(lines), flush=True)


if __name__ == '__main__':
    main()
