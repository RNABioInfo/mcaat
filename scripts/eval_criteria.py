#!/usr/bin/env python3
"""Score the boundary_walk criteria against a local reference.

For each <criterion>.fasta written by boundary_walk, reports:

  matching      one alignment with at least --identity and --coverage of the contig
  loci          those alignments, separate unless two start within --tolerance
  >=Nbp         one alignment with at least --identity and --absolute aligned bases,
                independent of contig length
  over spacer   the alignment covers at least --spacer-overlap of the spacer, whose
                position comes from the FASTA header rather than being inferred

  python3 eval_criteria.py --dir BW_OUT --db all_curated_db --out DIR
"""
import argparse
import collections
import glob
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


def loci_set(positions, tolerance=100):
    """Group alignment starts on each record into loci.

    Two alignments on one record are the same locus only if their starts are
    within tolerance.  Applied to a run of starts each within tolerance of the
    next, that makes the whole run one locus, which is the only reading that
    does not depend on which alignment you start from.  Returns the set of
    (record, first start of the group) pairs, so two runs can be compared.
    """
    out = set()
    for record, starts in positions.items():
        starts = sorted(starts)
        if not starts:
            continue
        anchor = starts[0]
        out.add((record, anchor))
        for a, b in zip(starts, starts[1:]):
            if b - a > tolerance:
                anchor = b
                out.add((record, anchor))
    return out


def count_loci(positions, tolerance=100):
    """How many loci the alignment starts fall into; see loci_set."""
    return len(loci_set(positions, tolerance))


def read_fasta(path):
    records, name = {}, None
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('>'):
            fields = line[1:].split()
            name = fields[0]
            meta = dict(f.split('=', 1) for f in fields[1:] if '=' in f)
            records[name] = [meta, []]
        elif name is not None:
            records[name][1].append(line.strip())
    return {k: (m, ''.join(s).upper()) for k, (m, s) in records.items()}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--dir', required=True)
    p.add_argument('--db', required=True)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--blast-bin', default=None)
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--identity', type=float, default=.80)
    p.add_argument('--coverage', type=float, default=.50)
    p.add_argument('--absolute', type=int, default=400)
    p.add_argument('--tolerance', type=int, default=100)
    p.add_argument('--spacer-overlap', type=float, default=.90)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    runs = {}
    for path in sorted(glob.glob(os.path.join(args.dir, '*.fasta'))):
        runs[Path(path).stem] = read_fasta(path)
    if not runs:
        sys.exit('no FASTA files in ' + args.dir)

    seqs = {}
    for records in runs.values():
        for _, seq in records.values():
            seqs.setdefault(seq, 'C%07d' % (len(seqs) + 1))
    hits = args.out / 'hits.tsv'
    if not hits.is_file():
        fa = args.out / 'contigs.fa'
        with open(fa, 'w') as f:
            for seq, key in seqs.items():
                f.write('>%s\n%s\n' % (key, seq))
        exe = shutil.which('blastn', path=args.blast_bin) if args.blast_bin else shutil.which('blastn')
        if not exe:
            sys.exit('blastn not found; pass --blast-bin')
        print('aligning %d distinct contigs' % len(seqs), flush=True)
        subprocess.run([exe, '-task', 'blastn', '-query', str(fa), '-db', args.db, '-dust', 'no',
                        '-soft_masking', 'false', '-evalue', '1e-10',
                        '-perc_identity', str(int(args.identity * 100)),
                        '-max_target_seqs', '50', '-num_threads', str(args.threads), '-outfmt',
                        '6 qseqid sseqid length nident qstart qend sstart send qlen',
                        '-out', str(hits)], check=True)

    best, long_best = {}, {}
    spans = collections.defaultdict(list)
    for line in open(hits):
        f = line.rstrip('\n').split('\t')
        length, nident, qstart, qend, sstart, send, qlen = (int(x) for x in f[2:9])
        if nident / length < args.identity:
            continue
        if length >= args.absolute and (f[0] not in long_best or nident > long_best[f[0]][0]):
            long_best[f[0]] = (nident, f[1], min(sstart, send))
        if (qend - qstart + 1) / qlen < args.coverage:
            continue
        spans[f[0]].append((min(qstart, qend), max(qstart, qend)))
        if f[0] not in best or nident > best[f[0]][0]:
            best[f[0]] = (nident, f[1], min(sstart, send))

    header = ('criterion', 'contigs', 'matching', 'loci', 'records',
              '>=%dbp' % args.absolute, 'loci', 'over spacer', 'median bp')
    print('%-22s %8s %9s %6s %8s %8s %6s %12s %10s' % header)
    rows = []
    for name, records in runs.items():
        placed = collections.defaultdict(list)
        long_placed = collections.defaultdict(list)
        matching = long_matching = over = 0
        lengths = []
        for meta, seq in records.values():
            lengths.append(len(seq))
            key = seqs[seq]
            lb = long_best.get(key)
            if lb:
                long_matching += 1
                long_placed[lb[1]].append(lb[2])
            b = best.get(key)
            if not b:
                continue
            matching += 1
            placed[b[1]].append(b[2])
            left, spacer = int(meta.get('left', 0)), int(meta.get('spacer', 0))
            lo, hi = left + 1, left + spacer
            if spacer and any(min(e, hi) - max(s, lo) + 1 >= args.spacer_overlap * spacer
                              for s, e in spans[key]):
                over += 1
        median = sorted(lengths)[len(lengths) // 2] if lengths else 0
        loci = count_loci(placed, args.tolerance)
        long_loci = count_loci(long_placed, args.tolerance)
        print('%-22s %8d %9d %6d %8d %8d %6d %12d %10d' %
              (name, len(records), matching, loci, len(placed), long_matching, long_loci, over, median))
        rows.append(dict(criterion=name, contigs=len(records), matching=matching, loci=loci,
                         records=len(placed), matching_absolute=long_matching, loci_absolute=long_loci,
                         over_spacer=over, median_length=median))
    (args.out / 'criteria.json').write_text(json.dumps(rows, indent=2) + '\n')


if __name__ == '__main__':
    main()
