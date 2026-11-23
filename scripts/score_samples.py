#!/usr/bin/env python3
"""Score boundary_walk runs over many samples and total them per criterion.

One sample cannot settle which criterion to use: on SRR14017565 harmonic_mean
found 93 strict matches against locally_smooth's 29, and on SRR14017570 that
order reversed. This aligns every sample, reports each one, and totals them, so
a criterion is judged on all of them rather than on whichever was looked at
first.

A strict match is one alignment with at least --identity over at least
--absolute aligned bases whose aligned part covers at least --overlap of the
spacer, whose position is read from the FASTA header rather than guessed. Two
alignments on one reference record are the same locus only if their starts are
within --tolerance. "both" counts the matches where that same alignment also
reaches --margin bases past each end of the spacer, which is the direct test of
the two-sided definition.

  python3 score_samples.py --root DIR --label bw --db all_curated_db --out DIR
  python3 score_samples.py --dirs a/bw b/bw --db ... --out DIR

With --baseline the production ContiguousPaths.fasta of each sample is scored
the same way, for the comparison the whole thing is about. It needs the
pipeline's own associations.tsv, which records how far each flank got.
"""
import argparse
import collections
import csv
import glob
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

from eval_criteria import loci_set, read_fasta


def blast(fasta_paths, out_dir, db, blast_bin, threads, identity):
    """One BLAST over the distinct sequences of every file, keyed by sequence."""
    out_dir.mkdir(parents=True, exist_ok=True)
    hits = out_dir / 'hits.tsv'
    seqs = {}
    for path in fasta_paths:
        for _, seq in read_fasta(path).values():
            seqs.setdefault(seq, 'C%07d' % (len(seqs) + 1))
    if not hits.is_file():
        query = out_dir / 'contigs.fa'
        with open(query, 'w') as f:
            for seq, key in seqs.items():
                f.write('>%s\n%s\n' % (key, seq))
        exe = shutil.which('blastn', path=blast_bin) if blast_bin else shutil.which('blastn')
        if not exe:
            sys.exit('blastn not found; pass --blast-bin')
        print('aligning %d distinct contigs' % len(seqs), flush=True)
        subprocess.run([exe, '-task', 'blastn', '-query', str(query), '-db', db, '-dust', 'no',
                        '-soft_masking', 'false', '-evalue', '1e-10',
                        '-perc_identity', str(int(identity * 100)), '-max_target_seqs', '50',
                        '-num_threads', str(threads), '-outfmt',
                        '6 qseqid sseqid length nident qstart qend sstart send qlen',
                        '-out', str(hits)], check=True)
    return seqs, hits


def read_hits(hits, identity, absolute):
    spans = collections.defaultdict(list)
    where = {}
    for line in open(hits):
        f = line.rstrip('\n').split('\t')
        length, nident, qs, qe, ss, se, _ = (int(x) for x in f[2:9])
        if nident / length < identity or length < absolute:
            continue
        spans[f[0]].append((min(qs, qe), max(qs, qe)))
        if f[0] not in where or nident > where[f[0]][2]:
            where[f[0]] = (f[1], min(ss, se), nident)
    return spans, where


def score(records, seqs, spans, where, overlap, margin):
    """Strict matches, how many reach both flanks, and the loci they fall in."""
    kept = len(records)
    matches = both = 0
    placed, placed_both = collections.defaultdict(list), collections.defaultdict(list)
    for meta, seq in records.values():
        key = seqs[seq]
        spacer, left = int(meta.get('spacer', 0)), int(meta.get('left', 0))
        if not spacer or key not in spans:
            continue
        lo, hi = left + 1, left + spacer
        picked = None
        for s, e in spans[key]:
            if min(e, hi) - max(s, lo) + 1 >= overlap * spacer:
                picked = (s, e)
                break
        if not picked:
            continue
        matches += 1
        record, start, _ = where[key]
        placed[record].append(start)
        if lo - picked[0] >= margin and picked[1] - hi >= margin:
            both += 1
            placed_both[record].append(start)
    return dict(kept=kept, matches=matches, both=both, placed=placed, placed_both=placed_both)


COMPLEMENT = str.maketrans('ACGTN', 'TGCAN')


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def baseline_records(sample_dir, seeds):
    """The production contigs, with the spacer located inside each one.

    associations.tsv records left_bp and right_bp, but in the target's
    orientation rather than the contig's: on SRR14017565, 70 of 286 contigs
    come out with the two swapped, which moves the spacer to the wrong end and
    makes the spacer-coverage test meaningless. Searching the contig for the
    spacer is checkable instead — the piece at the position this returns is a
    real spacer by construction.

    A contig whose spacer cannot be found is left out rather than guessed at,
    and the count is reported.
    """
    contigs = Path(sample_dir) / seeds / 'ContiguousPaths.fasta'
    table = Path(sample_dir) / seeds / 'sources.tsv'
    if not contigs.is_file() or not table.is_file():
        return None, 0
    spacers = set()
    with open(table) as f:
        head = f.readline().rstrip('\n').split('\t')
        cols = [head.index(c) for c in ('spacer', 'original_spacer') if c in head]
        for line in f:
            row = line.rstrip('\n').split('\t')
            for c in cols:
                if len(row) > c and row[c]:
                    spacers.add(row[c].upper())
    # Longest first, so a spacer that contains a shorter one is not mistaken for it.
    ordered = sorted(spacers, key=len, reverse=True)

    out, missing = {}, 0
    for name, (_, seq) in read_fasta(str(contigs)).items():
        seq = seq.upper()
        found = None
        for sp in ordered:
            for oriented in (sp, revcomp(sp)):
                at = seq.find(oriented)
                if at >= 0:
                    found = (at, len(oriented))
                    break
            if found:
                break
        if not found:
            missing += 1
            continue
        left, length = found
        out[name] = ({'left': str(left), 'right': str(len(seq) - left - length),
                      'spacer': str(length)}, seq)
    return out, missing


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--root', help='directory holding SRR*/ sample folders')
    p.add_argument('--label', default='bw', help='the run folder inside each sample')
    p.add_argument('--seeds', default='evidence19_seeds')
    p.add_argument('--dirs', nargs='*', default=[], help='run folders, instead of --root')
    p.add_argument('--db', required=True)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--blast-bin', default=None)
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--baseline', action='store_true')
    p.add_argument('--identity', type=float, default=.80)
    p.add_argument('--absolute', type=int, default=400)
    p.add_argument('--overlap', type=float, default=.90)
    p.add_argument('--margin', type=int, default=50)
    p.add_argument('--tolerance', type=int, default=100)
    args = p.parse_args()

    runs = list(args.dirs)
    if args.root:
        runs += sorted(d for d in glob.glob(os.path.join(args.root, '*', args.label))
                       if os.path.isdir(d))
    if not runs:
        sys.exit('no run folders found')

    # Sample name, run folder, and the production contigs beside it.
    samples = []
    for run in runs:
        run = run.rstrip(os.sep).rstrip('/')
        sample = os.path.basename(os.path.dirname(run)) or os.path.basename(run)
        fastas = sorted(glob.glob(os.path.join(run, '*.fasta')))
        if not fastas:
            print('%s: no fasta in %s, skipped' % (sample, run), file=sys.stderr)
            continue
        base = None
        if args.baseline:
            base, missing = baseline_records(os.path.dirname(run), args.seeds)
            if base is None:
                print('%s: no production contigs to compare with' % sample, file=sys.stderr)
            elif missing:
                print('%s: %d production contigs left out, no spacer of the sample found in them'
                      % (sample, missing), file=sys.stderr)
        samples.append((sample, run, fastas, base))

    every = [f for _, _, fs, _ in samples for f in fs]
    seqs, hits = blast(every, args.out, args.db, args.blast_bin, args.threads, args.identity)
    for _, _, _, base in samples:
        if base:
            for _, seq in base.values():
                seqs.setdefault(seq, 'B%07d' % (len(seqs) + 1))
    if any(b for _, _, _, b in samples):
        base_files = []
        for sample, run, _, base in samples:
            if not base:
                continue
            path = args.out / ('%s_baseline.fasta' % sample)
            with open(path, 'w') as f:
                for name, (meta, seq) in base.items():
                    f.write('>%s left=%s right=%s spacer=%s\n%s\n'
                            % (name, meta['left'], meta['right'], meta['spacer'], seq))
            base_files.append(str(path))
        bseqs, bhits = blast(base_files, args.out / 'baseline', args.db, args.blast_bin,
                             args.threads, args.identity)
        seqs.update(bseqs)
        bspans, bwhere = read_hits(bhits, args.identity, args.absolute)
    spans, where = read_hits(hits, args.identity, args.absolute)

    rows = []
    totals = collections.defaultdict(lambda: dict(kept=0, matches=0, both=0, loci=0, loci2=0))
    header = '%-14s %-22s %6s %8s %7s %6s %6s %6s' % (
        'sample', 'run', 'kept', 'matches', 'rate', 'both', 'loci', 'loci2')
    print(header)
    print('-' * len(header))
    for sample, run, fastas, base in samples:
        entries = [(Path(f).stem, read_fasta(f), spans, where) for f in fastas]
        if base:
            entries.append(('PRODUCTION', base, bspans, bwhere))
        for name, records, sp, wh in entries:
            r = score(records, seqs, sp, wh, args.overlap, args.margin)
            loci = len(loci_set(r['placed'], args.tolerance))
            loci2 = len(loci_set(r['placed_both'], args.tolerance))
            print('%-14s %-22s %6d %8d %7s %6d %6d %6d'
                  % (sample, name, r['kept'], r['matches'],
                     '%.1f%%' % (100. * r['matches'] / r['kept']) if r['kept'] else '-',
                     r['both'], loci, loci2))
            t = totals[name]
            t['kept'] += r['kept']
            t['matches'] += r['matches']
            t['both'] += r['both']
            t['loci'] += loci
            t['loci2'] += loci2
            rows.append(dict(sample=sample, run=name, kept=r['kept'], matches=r['matches'],
                             both=r['both'], loci=loci, loci_both=loci2))
    print()
    print('%-14s %-22s %6s %8s %7s %6s %6s %6s'
          % ('TOTAL', 'run', 'kept', 'matches', 'rate', 'both', 'loci', 'loci2'))
    print('-' * len(header))
    for name in sorted(totals, key=lambda k: -totals[k]['matches']):
        t = totals[name]
        print('%-14s %-22s %6d %8d %7s %6d %6d %6d'
              % ('', name, t['kept'], t['matches'],
                 '%.1f%%' % (100. * t['matches'] / t['kept']) if t['kept'] else '-',
                 t['both'], t['loci'], t['loci2']))
    (args.out / 'samples.json').write_text(json.dumps(rows, indent=2) + '\n')


if __name__ == '__main__':
    main()
