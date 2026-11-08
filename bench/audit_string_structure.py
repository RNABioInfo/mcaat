#!/usr/bin/env python3
"""Non-destructive, alignment-based audit of MCAAT repeat/spacer strings.

MAFFT aligns the repeat and spacer sections separately: graph-derived boundaries
remain anchored and every observed string is recovered exactly after removing
gaps. Existing variant families receive equal total weight, so a family with
many observed variants does not dominate the conservation profile. The rules
are exploratory review flags, NOT a calibrated CRISPR classifier.
"""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time


def reverse(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def canonical(s):
    return min(s, reverse(s))


def fasta_records(text):
    name = None
    seq = []
    for line in text.splitlines():
        if line.startswith('>'):
            if name is not None:
                yield name, ''.join(seq).upper()
            name, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())
    if name is not None:
        yield name, ''.join(seq).upper()


def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def align_sequences(sequences, executable):
    if not sequences or any(not s or set(s) - set('ACGT') for s in sequences):
        raise ValueError('alignment requires nonempty A/C/G/T strings')
    unique = list(dict.fromkeys(sequences))
    if len(unique) == 1:
        return list(sequences)
    with tempfile.TemporaryDirectory(prefix='mcaat-string-align-') as tmp:
        source = Path(tmp) / 'sequences.fa'
        source.write_text(''.join(f'>s{i}\n{s}\n' for i, s in enumerate(unique)))
        command = [executable, '--nuc', '--quiet', '--thread', '1', '--retree', '1',
                   '--maxiterate', '0', '--inputorder', str(source)]
        result = subprocess.run(command, text=True, capture_output=True, timeout=120)
        if result.returncode:
            raise RuntimeError('MAFFT failed: ' + result.stderr[-2000:])
        aligned = dict(fasta_records(result.stdout))
    if set(aligned) != {f's{i}' for i in range(len(unique))}:
        raise ValueError('MAFFT changed sequence identifiers')
    if len({len(s) for s in aligned.values()}) != 1:
        raise ValueError('alignment rows have inconsistent widths')
    lookup = {}
    for i, sequence in enumerate(unique):
        row = aligned[f's{i}']
        if row.replace('-', '') != sequence:
            raise ValueError('alignment changed an observed sequence')
        lookup[sequence] = row
    return [lookup[s] for s in sequences]


def family_weights(families):
    if not families:
        raise ValueError('at least one family is required')
    counts = Counter(families)
    return [1 / (len(counts) * counts[f]) for f in families]


def column_profile(aligned, weights):
    if not aligned or len(aligned) != len(weights):
        raise ValueError('invalid alignment/weights')
    if len({len(s) for s in aligned}) != 1:
        raise ValueError('unaligned strings supplied to column_profile')
    total = sum(weights)
    if any(not math.isfinite(w) or w <= 0 for w in weights):
        raise ValueError('weights must be finite and positive')
    if any(set(s) - set('ACGT-') for s in aligned):
        raise ValueError('alignment contains unsupported symbols')
    result = []
    for column in zip(*aligned):
        mass = defaultdict(float)
        for base, weight in zip(column, weights):
            mass[base] += weight / total
        bases = {base: mass[base] for base in 'ACGT'}
        best = max(bases, key=lambda base: (bases[base], -ord(base)))
        occupancy = sum(bases.values())
        entropy = -sum(value * math.log2(value) for value in mass.values() if value > 0)
        result.append(dict(base=best, conservation=bases[best], occupancy=occupancy,
                           entropy_with_gaps=entropy))
    return result


def motif_fit(sequence):
    """Best substitution-only short motif with at least three nominal copies.

    This is descriptive identity, not a p-value. The motif is a diagnostic
    template and is never emitted as a reconstructed biological sequence.
    """
    best = (0.0, 0, '')
    for period in range(1, len(sequence) // 3 + 1):
        buckets = [Counter(sequence[i::period]) for i in range(period)]
        motif = ''.join(min(c, key=lambda b: (-c[b], b)) for c in buckets)
        identity = sum(max(c.values()) for c in buckets) / len(sequence)
        if identity > best[0] + 1e-12:
            best = (identity, period, motif)
    return best


def reference_index(root):
    if root is None:
        return [], {}
    if not root.is_dir():
        raise ValueError('reference directory does not exist: ' + str(root))
    references = {}
    for file in sorted(root.glob('*/reference/*')):
        if not file.is_file():
            continue
        records = list(fasta_records(file.read_text()))
        spacers = frozenset(canonical(s) for name, s in records
                            if name.lower().startswith('spacer') and s and set(s) <= set('ACGT'))
        if len(spacers) < 2:
            continue
        # Sequence-set deduplication supplies positive controls, not locus counts.
        references.setdefault(spacers, []).append(str(file))
    items = [dict(spacers=sorted(spacers), sources=sources)
             for spacers, sources in references.items()]
    index = defaultdict(set)
    for i, item in enumerate(items):
        for spacer in item['spacers']:
            index[spacer].add(i)
    return items, index


def score_candidate(name, system, rows, family_map, reference_lookup, args):
    families = [family_map[row['spacer']] for row in rows]
    weights = family_weights(families)
    repeat_alignment = align_sequences([r['repeat'] for r in rows], args.mafft)
    spacer_alignment = align_sequences([r['spacer'] for r in rows], args.mafft)
    rp, sp = column_profile(repeat_alignment, weights), column_profile(spacer_alignment, weights)
    rmean = sum(c['conservation'] for c in rp) / len(rp)
    smean = sum(c['conservation'] for c in sp) / len(sp)
    conserved = [c['conservation'] + 1e-12 >= args.column_support for c in sp]
    left = 0
    while left < len(sp) and conserved[left]:
        left += 1
    right = 0
    while right < len(sp) - left and conserved[len(sp) - 1 - right]:
        right += 1
    periodic = [motif_fit(r['repeat'] + r['spacer']) for r in rows]
    periodic_fraction = sum(w for w, (identity, _, _) in zip(weights, periodic)
                            if identity + 1e-12 >= args.motif_identity)
    conserved_fraction = sum(conserved) / len(sp)
    family_count = len(set(families))
    motif_flag = family_count >= 2 and periodic_fraction + 1e-12 >= args.motif_support
    conservation_flag = (family_count >= 3 and conserved_fraction + 1e-12 >= args.conserved_spacer_fraction
                         and rmean - smean <= args.max_contrast + 1e-12)
    inner_end = len(sp) - right if right else len(sp)
    interior = sp[left:inner_end]
    flank_flag = (max(left, right) >= 3 and len(interior) >= 5
                  and sum(c['conservation'] for c in interior) / len(interior) < args.column_support)
    counts = Counter()
    matched = set()
    for spacer in {canonical(r['spacer']) for r in rows}:
        if spacer in reference_lookup:
            matched.add(spacer)
            for reference in reference_lookup[spacer]:
                counts[reference] += 1
    best_reference, best_matches = min(counts.items(), key=lambda x: (-x[1], x[0])) if counts else (-1, 0)
    flags = []
    if motif_flag:
        flags.append('short_motif_across_unit')
    if conservation_flag:
        flags.append('conservation_across_spacer')
    if flank_flag:
        flags.append('conserved_spacer_flank')
    if family_count == 1:
        flags.append('single_similarity_family')
    score = dict(array_id=name, repeat=system['repeat'], observations=len(rows), spacer_families=family_count,
                 repeat_alignment_columns=len(rp), spacer_alignment_columns=len(sp),
                 repeat_conservation=rmean, spacer_conservation=smean, conservation_contrast=rmean-smean,
                 conserved_spacer_fraction=conserved_fraction, variable_spacer_columns=len(sp)-sum(conserved),
                 conserved_spacer_left_columns=left, conserved_spacer_right_columns=right,
                 short_motif_fraction=periodic_fraction,
                 mean_short_motif_identity=sum(w * p[0] for w, p in zip(weights, periodic)),
                 short_motif_rule=int(motif_flag), conserved_spacer_rule=int(conservation_flag),
                 proposed_review=int(motif_flag or conservation_flag),
                 exact_reference_spacers=len(matched), max_matches_one_reference=best_matches,
                 best_reference_id=best_reference, input_flags=system.get('quality_flags', '-'),
                 string_flags=','.join(flags) or '-')
    alignment = ['# Observed strings; no genomic spacer order inferred.',
                 '# Repeat and spacer aligned separately; | preserves the reported boundary.',
                 '# Weights give each existing similarity family equal total mass.',
                 '# obs\tfamily\tweight\taligned repeat | aligned spacer']
    alignment.extend(f"{r['observation_id']}\t{family}\t{w:.8g}\t{rep} | {spacer}"
                     for r, family, w, rep, spacer in zip(rows, families, weights, repeat_alignment, spacer_alignment))
    alignment.append('# Position-wise A/C/G/T majority is a diagnostic profile, not an inferred sequence.')
    return score, rp, sp, '\n'.join(alignment) + '\n', matched


def write_table(path, rows, fields=None):
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]), delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: f'{v:.8g}' if isinstance(v, float) else v for k, v in row.items()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--references', type=Path)
    parser.add_argument('--mafft', default=shutil.which('mafft'))
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--column-support', type=float, default=.85)
    parser.add_argument('--motif-identity', type=float, default=.90)
    parser.add_argument('--motif-support', type=float, default=.85)
    parser.add_argument('--conserved-spacer-fraction', type=float, default=.80)
    parser.add_argument('--max-contrast', type=float, default=.15)
    args = parser.parse_args()
    if args.input.resolve() == args.output.resolve():
        parser.error('output must differ from input')
    if not args.mafft:
        parser.error('MAFFT is required; use --mafft /path/to/mafft')
    if args.workers < 1:
        parser.error('workers must be positive')
    for value in [args.column_support, args.motif_identity, args.motif_support, args.conserved_spacer_fraction, args.max_contrast]:
        if not math.isfinite(value) or not 0 <= value <= 1:
            parser.error('score parameters must be finite fractions in [0, 1]')
    args.output.mkdir(parents=True, exist_ok=True)
    alignment_dir = args.output / 'alignments'
    alignment_dir.mkdir(exist_ok=True)
    files = ['CRISPR_graph_summary.tsv', 'CRISPR_graph_systems.tsv', 'CRISPR_graph_paths.tsv', 'CRISPR_spacer_variants.tsv']
    fingerprints = {}
    for name in files:
        path = args.input / name
        fingerprints[name] = dict(bytes=path.stat().st_size, sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    system_rows = table(args.input / 'CRISPR_graph_systems.tsv')
    systems = {r['array_id']: r for r in system_rows}
    if len(systems) != len(system_rows):
        parser.error('duplicate candidate IDs in systems table')
    if not systems:
        parser.error('no candidates to audit')
    if any(not re.fullmatch(r'Array_[0-9]+', name) for name in systems):
        parser.error('candidate IDs must have the form Array_N')
    observations = defaultdict(list)
    for row in table(args.input / 'CRISPR_graph_paths.tsv'):
        if set(row['repeat'] + row['spacer']) - set('ACGT'):
            raise ValueError('non-ACGT observed string')
        observations[row['array_id']].append(row)
    families = defaultdict(dict)
    for row in table(args.input / 'CRISPR_spacer_variants.tsv'):
        name, spacer, family = row['array_id'], row['variant'], row['family_id']
        if spacer in families[name] and families[name][spacer] != family:
            raise ValueError('variant belongs to multiple families')
        families[name][spacer] = family
    if set(systems) != set(observations) or set(systems) != set(families):
        raise ValueError('inconsistent candidate IDs across input tables')
    for name, rows in observations.items():
        if len(rows) != int(systems[name]['assigned_cycles']):
            raise ValueError('assigned observation count mismatch: ' + name)
        if set(families[name]) != {r['spacer'] for r in rows}:
            raise ValueError('variant/path table mismatch: ' + name)
    references, ref_index = reference_index(args.references)
    (args.output / 'reference_sets.json').write_text(json.dumps(references, indent=2) + '\n')
    started = time.monotonic()
    scores, profiles, matched_by_array = [], [], {}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs = {pool.submit(score_candidate, name, system, observations[name], families[name], ref_index, args): name
                for name, system in systems.items()}
        for i, future in enumerate(as_completed(jobs), 1):
            score, repeat_profile, spacer_profile, alignment, matched = future.result()
            name = jobs[future]
            scores.append(score)
            matched_by_array[name] = matched
            (alignment_dir / (name + '.txt')).write_text(alignment)
            for section, profile in [('repeat', repeat_profile), ('spacer', spacer_profile)]:
                for j, column in enumerate(profile):
                    profiles.append(dict(array_id=name, section=section, column=j+1, **column))
            if i % 50 == 0 or i == len(jobs):
                print(f'Scored {i}/{len(jobs)} candidates in {time.monotonic()-started:.1f} s', flush=True)
    def number(name):
        return int(name.split('_')[-1])
    scores.sort(key=lambda row: number(row['array_id']))
    profiles.sort(key=lambda row: (number(row['array_id']), row['section'], row['column']))
    write_table(args.output / 'candidate_scores.tsv', scores)
    write_table(args.output / 'column_profiles.tsv', profiles)
    all_matches = set().union(*matched_by_array.values()) if matched_by_array else set()
    impacts = []
    for rule in ['short_motif_rule', 'conserved_spacer_rule', 'proposed_review']:
        flagged = [row for row in scores if row[rule]]
        kept_matches = set().union(*(matched_by_array[row['array_id']] for row in scores if not row[rule]))
        lost = sorted(all_matches - kept_matches)
        impacts.append(dict(rule=rule, flagged_candidates=len(flagged), unflagged_candidates=len(scores)-len(flagged),
                            flagged_candidates_with_two_same_reference_spacers=sum(row['max_matches_one_reference'] >= 2 for row in flagged),
                            exact_reference_spacers_before=len(all_matches), exact_reference_spacers_after=len(kept_matches),
                            exact_reference_spacers_lost=len(lost)))
        (args.output / (rule + '_lost_reference_spacers.txt')).write_text('\n'.join(lost) + ('\n' if lost else ''))
    write_table(args.output / 'rule_impact.tsv', impacts)
    for name, fingerprint in fingerprints.items():
        if hashlib.sha256((args.input / name).read_bytes()).hexdigest() != fingerprint['sha256']:
            raise RuntimeError('input changed during audit; rerun into a fresh directory: ' + name)
    version = subprocess.run([args.mafft, '--version'], capture_output=True, text=True, timeout=15)
    summary = dict(input=str(args.input.resolve()), candidates=len(scores), observations=sum(map(len, observations.values())),
                   local_reference_sequence_sets=len(references), exact_reference_spacers_detected=len(all_matches),
                   candidates_with_two_spacers_from_one_local_reference=sum(r['max_matches_one_reference'] >= 2 for r in scores),
                   rules=impacts, seconds=time.monotonic()-started, parameters={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                   input_fingerprints=fingerprints,
                   mafft_version=(version.stdout + version.stderr).strip(),
                   interpretation='Exploratory string review flags. Unknown annotation matches are not negatives. Local reference coverage is incomplete; this is not sensitivity against all ~750 annotated arrays. No production outputs are changed.')
    (args.output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k:summary[k] for k in ['candidates','local_reference_sequence_sets','exact_reference_spacers_detected','candidates_with_two_spacers_from_one_local_reference','rules','seconds']}, indent=2))


if __name__ == '__main__':
    main()
