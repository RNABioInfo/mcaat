#!/usr/bin/env python3
"""Of the >=400bp alignments, how many actually cover the spacer?

A protospacer detection has to include the spacer itself; an alignment that
only covers a flank is the flank matching something else.  The spacer's
position inside the contig comes from the FASTA header, not from the hit.
"""
import argparse
import collections
import csv
import glob
import os
from pathlib import Path

from eval_criteria import count_loci, read_fasta


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--dir', required=True)
    p.add_argument('--eval', required=True)
    p.add_argument('--criterion', default='harmonic_mean')
    p.add_argument('--top', type=int, default=0)
    p.add_argument('--identity', type=float, default=.80)
    p.add_argument('--absolute', type=int, default=400)
    p.add_argument('--overlap', type=float, default=.90)
    args = p.parse_args()

    runs = {}
    for path in sorted(glob.glob(os.path.join(args.dir, '*.fasta'))):
        runs[Path(path).stem] = read_fasta(path)
    seqs = {}
    for records in runs.values():
        for _, seq in records.values():
            seqs.setdefault(seq, 'C%07d' % (len(seqs) + 1))

    spans = collections.defaultdict(list)
    where = {}
    for line in open(os.path.join(args.eval, 'hits.tsv')):
        f = line.rstrip('\n').split('\t')
        length, nident, qstart, qend = (int(x) for x in f[2:6])
        if nident / length < args.identity or length < args.absolute:
            continue
        spans[f[0]].append((min(qstart, qend), max(qstart, qend)))
        sstart, send = int(f[6]), int(f[7])
        if f[0] not in where or nident > where[f[0]][2]:
            where[f[0]] = (f[1], min(sstart, send), nident)

    scores = {}
    with open(os.path.join(args.dir, 'summary.tsv')) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['criterion'] == args.criterion:
                scores[row['source_id']] = float(row['score'])

    items = [(scores.get(h[len(args.criterion) + 1:], float('-inf')), m, s)
             for h, (m, s) in runs[args.criterion].items()]
    items.sort(key=lambda r: -r[0])
    if args.top:
        items = items[:args.top]

    aligned = over = partial = 0
    placed = collections.defaultdict(list)
    for _, meta, seq in items:
        key = seqs[seq]
        if key not in spans:
            continue
        aligned += 1
        left, spacer = int(meta.get('left', 0)), int(meta.get('spacer', 0))
        lo, hi = left + 1, left + spacer
        if not spacer:
            continue
        best = max(min(e, hi) - max(s, lo) + 1 for s, e in spans[key])
        if best >= args.overlap * spacer:
            over += 1
            r = where[key]
            placed[r[0]].append(r[1])
        elif best > 0:
            partial += 1
    print('%-17s %-5s  %4d  %5d  %4d  %4d  %5d'
          % (args.criterion, args.top or len(items), aligned, over,
             count_loci(placed, 100), len({r for r in placed}), aligned - over - partial))


if __name__ == '__main__':
    main()
