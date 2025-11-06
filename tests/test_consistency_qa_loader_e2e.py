# -*- coding: ascii -*-
"""End-to-end tests for QA loader fallback, schema parsing, and metadata."""

import os
import json
import unittest
import tempfile

from src.halogenator.report import _load_qa_stats_with_fallback


class TestQALoaderE2E(unittest.TestCase):
    """Exercise file discovery and schema parsing for QA loader."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = self.tmp.name
        # Simulate products table path
        self.products_dir = os.path.join(self.root, 'out', 'p0')
        os.makedirs(self.products_dir, exist_ok=True)
        self.products_table = os.path.join(self.products_dir, 'products_k1.parquet')
        # Touch the products file path (not actually read by loader)
        with open(self.products_table, 'w') as f:
            f.write('')

    def _write_json(self, path, payload):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as f:
            json.dump(payload, f)

    def test_accepts_attempts_zero(self):
        """Loader accepts QA when attempts == 0 and canonicalizes fields."""
        qa_payload = {
            'total': {
                'attempts': 0,
                'products': 0,
                'no_product_matches': 0,
                'template_unsupported': 0
            }
        }
        primary_path = os.path.join(self.products_dir, 'qa_summary.json')
        self._write_json(primary_path, qa_payload)

        stats, source = _load_qa_stats_with_fallback(self.products_table, {})
        self.assertEqual(stats['attempts'], 0)
        self.assertEqual(stats['products'], 0)
        self.assertEqual(stats['no_product_matches'], 0)
        self.assertEqual(stats['template_unsupported'], 0)
        self.assertEqual(source['path_type'], 'primary')
        self.assertEqual(source['schema'], 'total')

    def test_config_path_priority_over_primary(self):
        """Config-provided QA path takes precedence over primary path."""
        # Write different QA files to config and primary
        qa_primary = {'total': {'attempts': 1, 'products': 1}}
        qa_config = {'total': {'attempts': 9, 'products': 9}}
        primary_path = os.path.join(self.products_dir, 'qa_summary.json')
        self._write_json(primary_path, qa_primary)

        cfg_path = os.path.join(self.root, 'explicit', 'qa_summary.json')
        self._write_json(cfg_path, qa_config)

        stats, source = _load_qa_stats_with_fallback(self.products_table, {'qa_summary_path': cfg_path})
        self.assertEqual(source['path_type'], 'config')
        self.assertEqual(stats['attempts'], 9)
        self.assertEqual(stats['products'], 9)

    def test_schema_variants_are_detected(self):
        """All supported schemas are parsed and reported with proper schema tag."""
        schemas = {
            'total': {'total': {'attempts': 2, 'products': 1, 'no_product_matches': 1}},
            'enumeration.totals': {'enumeration': {'totals': {'attempts': 3, 'products': 2, 'no_product_matches': 1}}},
            'totals': {'totals': {'attempts': 4, 'products': 3}},
            'flat': {'attempts': 5, 'products': 4}
        }

        for expected_schema, payload in schemas.items():
            # Write to primary location per iteration
            qa_path = os.path.join(self.products_dir, 'qa_summary.json')
            self._write_json(qa_path, payload)
            stats, source = _load_qa_stats_with_fallback(self.products_table, {})
            self.assertEqual(source['path_type'], 'primary')
            self.assertEqual(source['schema'], expected_schema)
            # Canonicalization: ensure all keys exist as ints
            for key in ('attempts', 'products', 'no_product_matches', 'template_unsupported'):
                self.assertIsInstance(stats.get(key, 0), int)

    def test_benchmarks_resolved_relative_to_output_dir(self):
        """benchmarks/qa_summary.json is resolved under configured output_dir when present."""
        outdir = os.path.join(self.root, 'report_out')
        qa_path = os.path.join(outdir, 'benchmarks', 'qa_summary.json')
        self._write_json(qa_path, {'total': {'attempts': 7, 'products': 6}})

        stats, source = _load_qa_stats_with_fallback(self.products_table, {'output_dir': outdir})
        self.assertEqual(source['path_type'], 'benchmarks')
        self.assertEqual(stats['attempts'], 7)

    def test_degrade_mode_accepts_products_only(self):
        """When degrade is enabled, accept products-only schema with attempts=0."""
        qa_payload = {'total': {'products': 4}}
        primary_path = os.path.join(self.products_dir, 'qa_summary.json')
        self._write_json(primary_path, qa_payload)

        stats, source = _load_qa_stats_with_fallback(self.products_table, {'qa_loader_degrade': True})
        self.assertEqual(stats['attempts'], 0)
        self.assertEqual(stats['products'], 4)


if __name__ == '__main__':
    unittest.main()
