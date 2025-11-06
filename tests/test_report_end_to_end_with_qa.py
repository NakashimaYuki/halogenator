# -*- coding: ascii -*-
"""End-to-end tests tying QA loader into finalize and checks."""

import os
import json
import unittest
import tempfile

from src.halogenator.report import (
    _init_incremental_stats,
    _update_incremental_stats,
    _finalize_incremental_stats,
    _load_qa_stats_with_fallback,
)


class TestReportE2EWithQA(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = self.tmp.name
        # Prepare a products table path and parent pairs
        self.products_dir = os.path.join(self.root, 'out')
        os.makedirs(self.products_dir, exist_ok=True)
        self.products_table = os.path.join(self.products_dir, 'products_k1.parquet')
        with open(self.products_table, 'w') as f:
            f.write('')
        self.parent_pairs = [('CCO', 'p1'), ('CCC', 'p2')]

    def _write_json(self, path, payload):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as f:
            json.dump(payload, f)

    def _build_minimal_stats(self):
        stats = _init_incremental_stats()
        # Two products: one F and one Cl
        records = [
            {'rule': 'R1', 'halogen': 'F', 'k': 1, 'sanitize_ok': True, 'parent_smiles': 'CCO'},
            {'rule': 'R1', 'halogen': 'Cl', 'k': 1, 'sanitize_ok': True, 'parent_smiles': 'CCC'},
        ]
        for r in records:
            _update_incremental_stats(stats, r)
        return stats

    def test_with_qa_runs_mutex_and_surfaces_source(self):
        # Place QA in primary location with matching products and valid mutex
        qa_path = os.path.join(self.products_dir, 'qa_summary.json')
        qa_payload = {'total': {'attempts': 3, 'products': 2, 'no_product_matches': 1, 'template_unsupported': 0}}
        self._write_json(qa_path, qa_payload)

        stats = self._build_minimal_stats()
        qa_stats, qa_source = _load_qa_stats_with_fallback(self.products_table, {})

        cfg = {'subset': 'test'}
        # Store qa_source before finalize, as generate_summary_report does
        stats['qa_source'] = qa_source
        _finalize_incremental_stats(stats, self.parent_pairs, cfg, qa_stats)
        result = stats['qa_consistency']

        # three_way_mutex is checked and consistent
        self.assertTrue(result['three_way_mutex_checked'])
        self.assertIn('three_way_mutex', result['checks'])
        self.assertTrue(result['checks']['three_way_mutex']['consistent'])

        # Global conservation consistent
        self.assertIn('global_product_conservation', result['checks'])
        self.assertTrue(result['checks']['global_product_conservation']['consistent'])

        # qa_source surfaced
        self.assertIn('qa_source', result)
        self.assertEqual(result['qa_source'].get('path_type'), 'primary')

        # product_count_source is fixed
        self.assertEqual(result.get('product_count_source'), 'report:product_count')

    def test_without_qa_skips_mutex(self):
        stats = self._build_minimal_stats()
        cfg = {'subset': 'test'}
        _finalize_incremental_stats(stats, self.parent_pairs, cfg, None)
        result = stats['qa_consistency']
        # Flag must be false and no three_way_mutex check present
        self.assertFalse(result['three_way_mutex_checked'])
        self.assertNotIn('three_way_mutex', result['checks'])
        # Global conservation still present
        self.assertIn('global_product_conservation', result['checks'])


if __name__ == '__main__':
    unittest.main()
