#!/usr/bin/env python3
"""String-audit regressions; no graph build or biological labels are required."""
import importlib.util
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True

spec = importlib.util.spec_from_file_location(
    'audit', Path(__file__).resolve().parents[1] / 'bench/audit_string_structure.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class StringAuditTests(unittest.TestCase):
    def test_each_family_has_equal_weight(self):
        weights = audit.family_weights(['a'] * 9 + ['b'])
        self.assertAlmostEqual(sum(weights[:9]), .5)
        self.assertAlmostEqual(weights[-1], .5)
        profile = audit.column_profile(['A'] * 9 + ['C'], weights)
        self.assertAlmostEqual(profile[0]['conservation'], .5)

    def test_gaps_are_not_conserved_bases(self):
        profile = audit.column_profile(['A-', '--', '--', '--'], [.25] * 4)
        self.assertEqual(profile[0]['conservation'], .25)
        self.assertEqual(profile[1]['conservation'], 0)
        self.assertEqual(profile[1]['occupancy'], 0)

    def test_invalid_profiles_rejected(self):
        for rows, weights in [(['AC', 'A'], [1, 1]), (['N'], [1]),
                              (['A'], [0]), (['A'], [float('nan')])]:
            with self.assertRaises(ValueError):
                audit.column_profile(rows, weights)

    def test_short_motif_and_substitution(self):
        self.assertEqual(audit.motif_fit('ACGT' * 8), (1, 4, 'ACGT'))
        altered = 'TCGT' + 'ACGT' * 7
        self.assertEqual(audit.motif_fit(altered), (31/32, 4, 'ACGT'))
        self.assertLess(audit.motif_fit('AGTCGATGCATCAGGATCCGATGACCTGTAGC')[0], .9)

    def test_reference_sets_deduplicated_including_reverse_strands(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            folder = root / 'genome/reference'
            folder.mkdir(parents=True)
            a, b = 'AACCGTAC', 'GGCATAGA'
            (folder / 'one.fa').write_text(f'>array\nAAAA\n>spacer 1\n{a}\n>spacer 2\n{b}\n')
            (folder / 'two.fa').write_text(
                f'>array\nAAAA\n>spacer 1\n{audit.reverse(b)}\n>spacer 2\n{audit.reverse(a)}\n')
            references, index = audit.reference_index(root)
            self.assertEqual(len(references), 1)
            self.assertEqual(len(references[0]['sources']), 2)
            self.assertEqual(index[audit.canonical(a)], {0})

    @unittest.skipUnless(shutil.which('mafft'), 'MAFFT not installed')
    def test_real_alignment_preserves_every_input_and_duplicates(self):
        sequences = ['ATGCGTACCGTA', 'ATGCGTACGTA', 'ATGCGTACCGTA', 'ATGCGTTACCGTA']
        rows = audit.align_sequences(sequences, shutil.which('mafft'))
        self.assertEqual([r.replace('-', '') for r in rows], sequences)
        self.assertEqual(len(set(map(len, rows))), 1)
        self.assertEqual(rows[0], rows[2])

    def test_aligner_must_not_change_observed_strings(self):
        response = SimpleNamespace(returncode=0, stderr='', stdout='>s0\nAAAA\n>s1\nGGGG\n')
        with patch.object(audit.subprocess, 'run', return_value=response):
            with self.assertRaisesRegex(ValueError, 'changed an observed'):
                audit.align_sequences(['AAAA', 'CCCC'], 'mafft')

    def test_empty_or_non_dna_input_rejected(self):
        for sequences in [[], [''], ['ACTN']]:
            with self.assertRaises(ValueError):
                audit.align_sequences(sequences, 'mafft')

    def score(self, repeat, spacers):
        args = SimpleNamespace(mafft='unused', column_support=.85, motif_identity=.9,
                               motif_support=.85, conserved_spacer_fraction=.8, max_contrast=.15)
        rows = [dict(observation_id=str(i), repeat=repeat, spacer=s) for i, s in enumerate(spacers)]
        with patch.object(audit, 'align_sequences', side_effect=lambda strings, _: list(strings)):
            return audit.score_candidate('Array_1', {'repeat': repeat}, rows,
                                         {s: str(i) for i, s in enumerate(spacers)}, {}, args)[0]

    def test_variable_spacer_and_conserved_repeat_not_flagged(self):
        score = self.score('AGTCGATGCATCAGGATCCGATGA',
                           ['ATGCAACGTTGC', 'CGATGGCATACG', 'GTCATTAGCGTA', 'TACGCTTACGAT'])
        self.assertEqual(score['repeat_conservation'], 1)
        self.assertEqual(score['proposed_review'], 0)

    def test_boundary_extension_is_reviewed_separately(self):
        score = self.score('AGTCGATGCATCAGGATCCGATGA',
                           ['GCTATGCAACGTTGC', 'GCTCGATGGCATACG', 'GCTGTCATTAGCGTA', 'GCTTACGCTTACGAT'])
        self.assertEqual(score['conserved_spacer_left_columns'], 3)
        self.assertIn('conserved_spacer_flank', score['string_flags'])
        self.assertEqual(score['proposed_review'], 0)

    def test_conservation_across_spacer_is_flagged(self):
        score = self.score('AGTCGATGCATCAGGATCCGATGA',
                           ['ACTGCTAGCTAGCTAGCTAA', 'ACTGCTAGCTAGCTAGCTAC', 'ACTGCTAGCTAGCTAGCTAG'])
        self.assertEqual(score['conserved_spacer_rule'], 1)

    @unittest.skipUnless(shutil.which('mafft'), 'MAFFT not installed')
    def test_cli_flags_motif_and_preserves_positive_control(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'input'
            source.mkdir()
            ref = root / 'reference/G/reference'
            ref.mkdir(parents=True)
            repeats = ['ACGT' * 6, 'AGTCGATGCATCAGGATCCGATGA']
            spacers = [['ACGT' * 6, 'ACGT' * 5 + 'ACGA'],
                       ['ATGCAACGTTGC', 'CGATGGCATACG', 'GTCATTAGCGTA']]

            def table(name, rows):
                with (source / name).open('w') as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
                    writer.writeheader()
                    writer.writerows(rows)

            table('CRISPR_graph_summary.tsv', [{'arrays': 2}])
            table('CRISPR_graph_systems.tsv', [dict(array_id=f'Array_{i+1}', repeat=r,
                  assigned_cycles=len(ss)) for i, (r, ss) in enumerate(zip(repeats, spacers))])
            table('CRISPR_graph_paths.tsv', [dict(array_id=f'Array_{i+1}', observation_id=j,
                  repeat=r, spacer=s) for i, (r, ss) in enumerate(zip(repeats, spacers))
                  for j, s in enumerate(ss)])
            table('CRISPR_spacer_variants.tsv', [dict(array_id=f'Array_{i+1}', family_id=j,
                  variant=s) for i, ss in enumerate(spacers) for j, s in enumerate(ss)])
            (ref / 'control.fa').write_text(''.join(
                f'>spacer {i}\n{s}\n' for i, s in enumerate(spacers[1])))
            subprocess.run([sys.executable, spec.origin, str(source), str(root / 'output'),
                            '--references', str(root / 'reference')],
                           check=True, capture_output=True, text=True, timeout=30)
            data = json.loads((root / 'output/summary.json').read_text())
            self.assertEqual(data['candidates'], 2)
            self.assertEqual(data['observations'], 5)
            self.assertEqual(data['exact_reference_spacers_detected'], 3)
            impact = {r['rule']: r for r in data['rules']}
            self.assertEqual(impact['short_motif_rule']['flagged_candidates'], 1)
            self.assertEqual(impact['proposed_review']['exact_reference_spacers_lost'], 0)
            self.assertTrue(data['mafft_version'])
            self.assertEqual(len(data['input_fingerprints']), 4)


if __name__ == '__main__':
    unittest.main()
