#!/usr/bin/env python3
"""Create composition-matched background queries without reference sequences."""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import random


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def shuffle_dinucleotides(sequence, rng):
    """Randomize an Euler trail through nucleotide-pair edges; preserve pair counts."""
    edges = defaultdict(list)
    for a, b in zip(sequence, sequence[1:]):
        edges[a].append(b)
    for outgoing in edges.values():
        rng.shuffle(outgoing)
    stack, trail = [sequence[0]], []
    while stack:
        if edges[stack[-1]]:
            stack.append(edges[stack[-1]].pop())
        else:
            trail.append(stack.pop())
    result = ''.join(reversed(trail))
    assert len(result) == len(sequence)
    assert Counter(zip(result, result[1:])) == Counter(zip(sequence, sequence[1:]))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seed', type=int, default=20261115)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        parser.error('control output must be new or empty')
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.results/'queries.tsv').open() as stream:
        queries = list(csv.DictReader(stream, delimiter='\t'))
    with (args.results/'sources.tsv').open() as stream:
        sources = list(csv.DictReader(stream, delimiter='\t'))
    rng = random.Random(args.seed)
    used = {q['spacer'] for q in queries}
    controls, records = {}, []
    for query in queries:
        original = query['spacer']
        control = ''
        for attempt in range(1000):
            shuffled = shuffle_dinucleotides(original, rng)
            canonical = min(shuffled, reverse_complement(shuffled))
            if canonical not in used:
                control = shuffled
                used.add(canonical)
                controls[original] = shuffled
                break
        records.append(dict(original_query=query['query_id'], original=original,
                            control=control, status='generated' if control else 'no_distinct_shuffle_found'))
    fields = ['array_id', 'observation_id', 'repeat', 'spacer']
    written = 0
    with (args.output/'CRISPR_graph_paths.tsv').open('w') as stream:
        writer = csv.DictWriter(stream, fields, delimiter='\t'); writer.writeheader()
        for source in sources:
            original = source['spacer']
            canonical = min(original, reverse_complement(original))
            if canonical not in controls:
                continue
            shuffled = controls[canonical]
            if original != canonical:
                shuffled = reverse_complement(shuffled)
            writer.writerow(dict(array_id='control:'+source['array_id'], observation_id=source['source_id'],
                                 repeat=source['repeat'], spacer=shuffled))
            written += 1
    with (args.output/'control_pairs.tsv').open('w') as stream:
        writer = csv.DictWriter(stream, list(records[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(records)
    report = dict(seed=args.seed, original_queries=len(queries), distinct_controls=len(controls),
                  source_observations=written, unavailable_controls=len(queries)-len(controls),
                  preserved='sequence length and exact dinucleotide counts',
                  boundary_correction=False, reference_used=False,
                  interpretation='Background controls, not certified biological negatives or a calibrated false-discovery rate.')
    (args.output/'controls.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
