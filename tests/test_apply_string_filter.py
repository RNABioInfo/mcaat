#!/usr/bin/env python3
"""Test applying cached decisions, input binding, and production serialization."""
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import random
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True
parser = argparse.ArgumentParser()
parser.add_argument('--replay-binary', type=Path, required=True)
arguments, remaining = parser.parse_known_args()
BINARY = arguments.replay_binary.resolve()
spec = importlib.util.spec_from_file_location(
    'apply_filter', Path(__file__).resolve().parents[1] / 'bench/apply_string_filter.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def table(path, records):
    with path.open('w', newline='') as stream:
        # CRLF intentionally exercises report and rejection-manifest parsing.
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(records)


def fixture(root):
    source, audit = root / 'input', root / 'audit'
    source.mkdir()
    audit.mkdir()
    rng = random.Random(1829)
    def dna(n):
        return ''.join(rng.choice('ACGT') for _ in range(n))
    k = 9
    vertices, edges = {}, {}
    systems, paths, variants, scores = [], [], [], []
    for array in range(1, 3):
        repeat = dna(28)
        repeat_path = None
        for obs in range(1, 4):
            while True:
                spacer = dna(30)
                sequence = repeat + spacer
                padded = sequence + sequence[:k-1]
                labels = [padded[i:i+k] for i in range(len(sequence))]
                if len(set(labels)) == len(labels):
                    break
            ids = [vertices.setdefault(s, len(vertices)+1) for s in labels]
            links = [edges.setdefault((ids[i], ids[(i+1) % len(ids)]), len(edges)+1)
                     for i in range(len(ids))]
            repeat_path = ','.join(map(str, ids[:len(repeat)-k+1]))
            paths.append(dict(array_id=f'Array_{array}', observation_id=obs,
                              repeat=repeat, spacer=spacer, repeat_vertices=repeat_path,
                              cycle_vertices=','.join(map(str, ids)), cycle_edges=','.join(map(str, links)),
                              supports_repeat_path='1'))
            variants.append(dict(array_id=f'Array_{array}', family_id=obs, variant=spacer))
        systems.append(dict(array_id=f'Array_{array}', repeat=repeat, repeat_vertices=repeat_path,
                            source_start_edges=str(array), support_cycles=3, candidate_cycles=3,
                            assigned_cycles=3, alternative_repeat_paths=0))
        # These synthetic flags test applying explicit decisions, not the scorer.
        scores.append(dict(array_id=f'Array_{array}', short_motif_rule=int(array == 1),
                           conserved_spacer_rule=0, proposed_review=int(array == 1)))
    table(source / 'CRISPR_graph_systems.tsv', systems)
    table(source / 'CRISPR_graph_paths.tsv', paths)
    table(source / 'CRISPR_spacer_variants.tsv', variants)
    metrics = dict(k=k, repeat_support=.85, cycle_min_length=47, cycle_max_length=100,
                   repeat_min_length=23, repeat_max_length=55, spacer_min_length=23, spacer_max_length=50)
    table(source / 'CRISPR_graph_summary.tsv', [dict(metric=k, value=v) for k, v in metrics.items()])
    table(audit / 'candidate_scores.tsv', scores)
    data = dict(observations=6, parameters={'synthetic_test': True}, input_fingerprints={
        name: dict(bytes=(source / name).stat().st_size, sha256=module.digest(source / name))
        for name in module.INPUT_FILES})
    (audit / 'summary.json').write_text(json.dumps(data))
    (audit / 'reference_sets.json').write_text(json.dumps([
        {'spacers': [r['spacer'] for r in paths if r['array_id'] == 'Array_2']}]))
    return source, audit


class ApplyTests(unittest.TestCase):
    def test_motif_conservation_union_and_boundary_only(self):
        scores = [dict(array_id='Array_1', short_motif_rule='1', conserved_spacer_rule='1', proposed_review='1'),
                  dict(array_id='Array_2', short_motif_rule='0', conserved_spacer_rule='0', proposed_review='0',
                       string_flags='conserved_spacer_flank')]
        self.assertEqual(module.rejection_rows(scores), [dict(
            array_id='Array_1', reason='short_motif_across_unit+conservation_across_spacer')])

    def test_inconsistent_rule_union_fails(self):
        with self.assertRaises(ValueError):
            module.rejection_rows([dict(array_id='Array_1', short_motif_rule='1',
                                        conserved_spacer_rule='0', proposed_review='0')])

    def test_changed_input_fails_before_creating_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, audit = fixture(root)
            with (source / 'CRISPR_graph_paths.tsv').open('a') as stream:
                stream.write('\n')
            with self.assertRaisesRegex(ValueError, 'differs from the audited input'):
                module.apply_filter(source, audit, root / 'output', BINARY)
            self.assertFalse((root / 'output').exists())

    def test_candidate_mismatch_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, audit = fixture(root)
            records = module.rows(audit / 'candidate_scores.tsv')
            records[0]['array_id'] = 'Array_999'
            table(audit / 'candidate_scores.tsv', records)
            with self.assertRaisesRegex(ValueError, 'candidate IDs'):
                module.apply_filter(source, audit, root / 'output', BINARY)

    def test_production_filter_archives_rejections_and_keeps_controls(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, audit = fixture(root)
            output = root / 'output'
            result = module.apply_filter(source, audit, output, BINARY)
            self.assertEqual(result['string_rejected_candidates'], 1)
            self.assertEqual(result['retained_candidates'], 1)
            self.assertEqual(result['retained_spacer_entries'], 3)
            self.assertEqual(result['alternative_observations'], 3)
            self.assertEqual(result['local_reference_matches_before'], 3)
            self.assertEqual(result['local_reference_matches_after'], 3)
            self.assertEqual(result['lost_local_reference_spacers'], [])
            removed = module.rows(output / 'CRISPR_graph_alternatives.tsv')
            original = module.rows(source / 'CRISPR_graph_paths.tsv')
            expected = {(r['repeat'], r['spacer'], r['cycle_vertices']) for r in original if r['array_id'] == 'Array_1'}
            self.assertEqual({(r['repeat'], r['spacer'], r['cycle_vertices']) for r in removed}, expected)
            self.assertTrue(all(r['reason'] == 'string_filter:short_motif_across_unit' for r in removed))
            module.verify_input(source, json.loads((audit / 'summary.json').read_text()))
            with self.assertRaisesRegex(ValueError, 'new or empty'):
                module.apply_filter(source, audit, output, BINARY)


if __name__ == '__main__':
    unittest.main(argv=[sys.argv[0], *remaining])
