#!/usr/bin/env python3
"""Apply the harmonic mean rule to a boundary_walk run and score what it keeps.

The rule reads the harmonic_mean column boundary_walk records per contig, so
nothing here depends on a rank or a count.

  python3 rule_table.py --dir BW_OUT --eval EVAL_OUT [--thresholds 1 1.5 2]
"""
import argparse
import collections
import csv
import glob
import os
from pathlib import Path

from eval_criteria import count_loci, read_fasta


def load(directory, eval_dir, identity=.80, absolute=400):
    runs = {}
    for p in sorted(glob.glob(os.path.join(directory, '*.fasta'))):
        runs[Path(p).stem] = read_fasta(p)
    seqs = {}
    for recs in runs.values():
        for _, s in recs.values():
            seqs.setdefault(s, 'C%07d' % (len(seqs) + 1))
    spans, where = collections.defaultdict(list), {}
    for line in open(os.path.join(eval_dir, 'hits.tsv')):
        f = line.rstrip('\n').split('\t')
        length, nident, qs, qe, ss, se, _ = (int(x) for x in f[2:9])
        if nident / length < identity or length < absolute:
            continue
        spans[f[0]].append((min(qs, qe), max(qs, qe)))
        if f[0] not in where or nident > where[f[0]][2]:
            where[f[0]] = (f[1], min(ss, se), nident)
    harmonic = collections.defaultdict(dict)
    with open(os.path.join(directory, 'summary.tsv')) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            # A run from another pipeline has no such column; it then has no
            # rule to apply and every contig passes.
            harmonic[row['criterion']][row['source_id']] = float(row.get('harmonic_mean') or 0.)
    items = collections.defaultdict(list)
    for crit, recs in runs.items():
        for h, (meta, seq) in recs.items():
            key = seqs[seq]
            sp, left = int(meta.get('spacer', 0)), int(meta.get('left', 0))
            lo, hi = left + 1, left + sp
            covers = bool(sp and key in spans and
                          max(min(e, hi) - max(s, lo) + 1 for s, e in spans[key]) >= .9 * sp)
            items[crit].append(dict(h=harmonic[crit].get(h[len(crit) + 1:], 0.),
                                    covers=covers, place=where.get(key), n=len(seq)))
    return items


def row(label, sel):
    placed = collections.defaultdict(list)
    for i in sel:
        if i['covers'] and i['place']:
            placed[i['place'][0]].append(i['place'][1])
    m = sum(1 for i in sel if i['covers'])
    lengths = sorted(i['n'] for i in sel)
    print('%-22s %6d %8d %7s %6d %6d %9d'
          % (label, len(sel), m, '%.1f%%' % (100 * m / len(sel)) if sel else '-',
             count_loci(placed), len(placed), lengths[len(lengths) // 2] if lengths else 0))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--dir', required=True)
    p.add_argument('--eval', required=True)
    p.add_argument('--thresholds', nargs='*', type=float, default=[0])
    p.add_argument('--criterion', default=None)
    args = p.parse_args()
    items = load(args.dir, args.eval)
    print('%-22s %6s %8s %7s %6s %6s %9s'
          % ('run', 'kept', 'matches', 'rate', 'loci', 'recs', 'median bp'))
    for t in args.thresholds:
        if len(args.thresholds) > 1:
            print('-- harmonic mean >= %g' % t)
        for crit in sorted(items):
            if args.criterion and crit != args.criterion:
                continue
            row(crit if len(args.thresholds) > 1 else '%s H>=%g' % (crit, t),
                [i for i in items[crit] if i['h'] >= t])


if __name__ == '__main__':
    main()
