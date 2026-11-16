#!/usr/bin/env python3
"""Evaluate saved competing paths after prediction; SQLite bounds memory use."""
import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import sqlite3
import subprocess
import time


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def prepare(results, output, flank_length):
    summary = json.loads((results / 'run_summary.json').read_text())
    if summary['status'] != 'complete' or not summary.get('save_search'):
        raise ValueError('A complete run with --save-search true is required')
    db = sqlite3.connect(output / 'paths.sqlite')
    db.executescript('''
        CREATE TABLE sequences (id INTEGER PRIMARY KEY, sequence TEXT UNIQUE NOT NULL);
        CREATE TABLE paths (path_id TEXT PRIMARY KEY, query_id TEXT, placement INTEGER,
                            side TEXT, flank_bp INTEGER, stop TEXT, limited INTEGER, sequence_id INTEGER);
        CREATE INDEX paths_sequence ON paths(sequence_id);
        CREATE INDEX paths_query ON paths(query_id);
    ''')
    total, eligible, shorter, stops = 0, 0, 0, Counter()
    for path in sorted((results / 'search_paths').glob('Q*.tsv.gz')):
        with gzip.open(path, 'rt') as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                total += 1
                stops[row['stop']] += 1
                flank = int(row['flank_bp'])
                if row['external'] != '1' or flank < max(summary['min_extension_flank'],
                        max(summary['graph_k'], int(row['clear_at'])) + summary['graph_k'] + 1):
                    continue
                if flank < flank_length:
                    shorter += 1
                    continue
                seq = row['sequence']
                assert set(seq) <= set('ACGT') and len(seq) == flank + int(row['query_end']) - int(row['query_begin'])
                seq = min(seq, reverse_complement(seq))
                db.execute('INSERT OR IGNORE INTO sequences(sequence) VALUES (?)', (seq,))
                sequence_id = db.execute('SELECT id FROM sequences WHERE sequence=?', (seq,)).fetchone()[0]
                db.execute('INSERT INTO paths VALUES (?,?,?,?,?,?,?,?)',
                    (row['path_id'], row['query_id'], int(row['placement']), row['side'], flank,
                     row['stop'], int(row['limited']), sequence_id))
                eligible += 1
        db.commit()
    with (output / 'candidates.fasta').open('w') as stream:
        for number, seq in db.execute('SELECT id,sequence FROM sequences ORDER BY id'):
            stream.write(f'>C{number:08d}\n{seq}\n')
    count = db.execute('SELECT COUNT(*) FROM sequences').fetchone()[0]
    db.close()
    report = dict(total_saved_paths=total, eligible_one_sided_paths=eligible,
                  distinct_eligible_sequences=count, stop_counts=dict(stops),
                  minimum_audited_flank=flank_length, shorter_external_paths_not_in_this_comparison=shorter,
                  interpretation='Saved external paths reaching the stated flank length. Shorter prefixes are counted separately. Not exhaustive graph reachability.',
                  reference_used_for_prediction=False)
    (output / 'preparation.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


def screen(output, database, threads):
    """Find positive witnesses; a long-seed nonmatch is not proof of absence."""
    folder = output / 'fast_witness_screen'
    if (folder/'summary.json').exists():
        return json.loads((folder/'summary.json').read_text())
    folder.mkdir(exist_ok=True)
    start = time.monotonic()
    db = sqlite3.connect(output/'paths.sqlite', timeout=60)
    db.execute('CREATE TABLE IF NOT EXISTS matches (sequence_id INTEGER PRIMARY KEY, reference TEXT, identity REAL, coverage REAL)')
    db.execute('CREATE TABLE IF NOT EXISTS screen_tested (sequence_id INTEGER PRIMARY KEY)')
    queries = [r[0] for r in db.execute('SELECT DISTINCT query_id FROM paths ORDER BY query_id')]
    round_number = len(list(folder.glob('round_*')))
    passed = set()
    while True:
        found = {r[0] for r in db.execute('SELECT DISTINCT query_id FROM paths JOIN matches USING(sequence_id)')}
        selected = set()
        for query in queries:
            if query in found:
                continue
            selected.update(r[0] for r in db.execute('''SELECT sequence_id FROM (
                SELECT sequence_id,placement,side,ROW_NUMBER() OVER (
                    PARTITION BY placement,side ORDER BY sequence_id) AS priority
                FROM paths WHERE query_id=? AND sequence_id NOT IN (SELECT sequence_id FROM screen_tested)
                GROUP BY placement,side,sequence_id)
                ORDER BY priority,placement,side LIMIT ?''', (query, min(128,8*2**min(round_number,4)))))
        if not selected:
            break
        batch = folder/f'round_{round_number+1:04d}'
        batch.mkdir()
        with (batch/'candidates.fasta').open('w') as stream:
            for number in sorted(selected):
                sequence = db.execute('SELECT sequence FROM sequences WHERE id=?',(number,)).fetchone()[0]
                stream.write(f'>C{number:08d}\n{sequence}\n')
        command = ['/usr/local/bin/blastn' if Path('/usr/local/bin/blastn').exists() else 'blastn',
            '-query', str(batch/'candidates.fasta'), '-db', str(database), '-task', 'blastn-short',
            '-word_size', '16', '-dust', 'no', '-soft_masking', 'false', '-evalue', '1000',
            '-perc_identity', '90', '-qcov_hsp_perc', '49', '-max_target_seqs', '1000000',
            '-num_threads', str(threads), '-outfmt',
            '6 qseqid sseqid qlen qstart qend sstart send length nident pident evalue bitscore',
            '-out', str(batch/'hits.tsv')]
        (batch/'command.json').write_text(json.dumps(command, indent=2)+'\n')
        print(f'Fast witness round {round_number+1}: {len(found)} queries resolved; {len(selected)} candidates',flush=True)
        with (batch/'blast.log').open('w') as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
        with (batch/'hits.tsv').open() as stream:
            for line in stream:
                r = line.rstrip().split('\t')
                identity = int(r[8])/int(r[7])
                coverage = (abs(int(r[4])-int(r[3]))+1)/int(r[2])
                if identity+1e-12 >= .9 and coverage+1e-12 >= .5:
                    number = int(r[0][1:])
                    passed.add(number)
                    db.execute('INSERT OR IGNORE INTO matches VALUES (?,?,?,?)', (number,r[1],identity,coverage))
        db.executemany('INSERT OR IGNORE INTO screen_tested VALUES (?)',[(i,) for i in selected])
        db.commit()
        round_number += 1
    db.close()
    report = dict(seconds=time.monotonic()-start, verified_matching_sequences=len(passed),
                  interpretation='Positive witnesses only. A negative fast screen never marks a candidate sensitively tested; optional --sensitive uses word size 9.')
    (folder/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


def evaluate(output, database, threads, sensitive=False):
    start = time.monotonic()
    db = sqlite3.connect(output / 'paths.sqlite')
    db.executescript('''CREATE TABLE IF NOT EXISTS matches
        (sequence_id INTEGER PRIMARY KEY, reference TEXT, identity REAL, coverage REAL);
        CREATE TABLE IF NOT EXISTS tested (sequence_id INTEGER PRIMARY KEY);''')
    screen(output, database, threads)
    all_queries = [r[0] for r in db.execute('SELECT DISTINCT query_id FROM paths ORDER BY query_id')]
    round_number = len(list(output.glob('round_*')))
    # This is an existence audit: stop testing a query once a qualifying witness
    # is found. For queries without a witness, visit every saved full-length path.
    while sensitive:
        found = {r[0] for r in db.execute('SELECT DISTINCT query_id FROM paths JOIN matches USING(sequence_id)')}
        per_query = min(128, 4 * 2**min(round_number, 5))
        selected = set()
        for query in all_queries:
            if query in found:
                continue
            selected.update(r[0] for r in db.execute('''SELECT sequence_id FROM (
                SELECT sequence_id,placement,side,ROW_NUMBER() OVER (
                    PARTITION BY placement,side ORDER BY sequence_id) AS priority
                FROM paths WHERE query_id=? AND sequence_id NOT IN (SELECT sequence_id FROM tested)
                GROUP BY placement,side,sequence_id)
                ORDER BY priority,placement,side LIMIT ?''', (query, per_query)))
        if not selected:
            break
        folder = output / f'round_{round_number+1:04d}'
        folder.mkdir()
        with (folder/'candidates.fasta').open('w') as stream:
            for number in sorted(selected):
                seq = db.execute('SELECT sequence FROM sequences WHERE id=?', (number,)).fetchone()[0]
                stream.write(f'>C{number:08d}\n{seq}\n')
        # >=90% identity across >=90 alignment columns implies an exact run of
        # at least nine bases. BLAST itself remains a heuristic alignment search.
        command = ['/usr/local/bin/blastn' if Path('/usr/local/bin/blastn').exists() else 'blastn',
            '-query', str(folder/'candidates.fasta'), '-db', str(database), '-task', 'blastn-short',
            '-word_size', '9', '-dust', 'no', '-soft_masking', 'false', '-evalue', '1000',
            '-perc_identity', '90', '-qcov_hsp_perc', '49', '-max_target_seqs', '1000000',
            '-num_threads', str(threads), '-outfmt',
            '6 qseqid sseqid qlen qstart qend sstart send length nident pident evalue bitscore',
            '-out', str(folder/'hits.tsv')]
        (folder/'command.json').write_text(json.dumps(command, indent=2)+'\n')
        print(f'Audit round {round_number+1}: {len(found)} queries have witnesses; testing {len(selected)} paths', flush=True)
        with (folder/'blast.log').open('w') as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
        with (folder/'hits.tsv').open() as stream:
            for line in stream:
                r = line.rstrip().split('\t')
                coverage = (abs(int(r[4])-int(r[3]))+1)/int(r[2])
                identity = int(r[8])/int(r[7])
                if identity+1e-12 >= .9 and coverage+1e-12 >= .5:
                    db.execute('INSERT OR IGNORE INTO matches VALUES (?,?,?,?)',
                               (int(r[0][1:]), r[1], identity, coverage))
        db.executemany('INSERT OR IGNORE INTO tested VALUES (?)', [(i,) for i in selected])
        db.commit()
        round_number += 1
    queries = [r[0] for r in db.execute('SELECT DISTINCT query_id FROM paths JOIN matches USING(sequence_id) ORDER BY query_id')]
    full = [r[0] for r in db.execute('SELECT DISTINCT query_id FROM paths JOIN matches USING(sequence_id) WHERE flank_bp=750 ORDER BY query_id')]
    with (output / 'matching_paths.tsv').open('w') as stream:
        cursor = db.execute('SELECT paths.*,reference,identity,coverage FROM paths JOIN matches USING(sequence_id) ORDER BY query_id,path_id')
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow([c[0] for c in cursor.description])
        writer.writerows(cursor)
    report = dict(status='complete', evaluation_scope='sensitive_pool_audit' if sensitive else 'long_seed_witness_screen',
                  sensitive_audit_complete=sensitive, queries_with_full_candidates=len(all_queries),
                  matching_queries=queries, matching_queries_count=len(queries),
                  queries_with_matching_750_base_flank=full, matching_750_base_queries_count=len(full),
                  verified_matching_sequences=db.execute('SELECT COUNT(*) FROM matches').fetchone()[0],
                  tested_sequences=db.execute('SELECT COUNT(*) FROM tested').fetchone()[0],
                  fast_screened_sequences=db.execute('SELECT COUNT(*) FROM screen_tested').fetchone()[0],
                  queries_without_a_verified_witness=sorted(set(all_queries)-set(queries)),
                  blast_seconds=time.monotonic()-start, identity_threshold=.9, coverage_threshold=.5,
                  interpretation='One qualifying witness suffices per query. Other candidates for a successful query need not be tested. The long-seed screen uses word size 16; optional sensitive fallback uses 9. These are confirmed witnesses, not a maximum recoverability estimate or proof that other graph routes are absent.',
                  reference_used_for_prediction=False)
    db.close()
    (output / 'evaluation.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--database', type=Path)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--flank-length', type=int, default=750)
    parser.add_argument('--prepared', action='store_true')
    parser.add_argument('--sensitive', action='store_true', help='Also exhaust the remaining candidate pool with word-size-9 BLAST; can be very expensive on repetitive paths.')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if not args.prepared:
        if any(args.output.iterdir()):
            parser.error('new audit output must be empty')
        print(json.dumps(prepare(args.results, args.output, args.flank_length), indent=2), flush=True)
    if args.database:
        print(json.dumps(evaluate(args.output, args.database, args.threads, args.sensitive), indent=2), flush=True)


if __name__ == '__main__':
    main()
