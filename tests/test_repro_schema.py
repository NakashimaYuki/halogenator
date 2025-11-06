# -*- coding: ascii -*-
"""Tests for verify_reproducibility schema handling."""

import sys
import pathlib
import unittest

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts import verify_reproducibility as repro  # noqa: E402


def _make_sample(schema_version, include_new_fields=True):
    """Construct a minimal acceptance artifact structure."""
    overall = {
        'p0_h2_total_samples': 11,
        'p0_h2_total_passes': 11,
        'product_total': 100
    }
    if include_new_fields:
        overall['accepted_via_score_count'] = 2
        overall['degraded_due_to_no_score_count'] = 0

    sample_entry = {
        'analysis': {
            'off_products': 10,
            'heuristic_products': 8,
            'off_attempts': 20,
            'heuristic_attempts': 18,
            'reduction_percentage': 20.0,
            'overall_pass': True,
            'sugar_events': {'total': 2}
        },
        'heuristic_mode': {
            'mask_size': 4,
            'degraded': False
        },
        'p0_h2_assertion': {
            'assertion_pass': True
        }
    }

    artifact = {
        'test_metadata': {
            'schema_version': schema_version,
            'random_seed': 12345,
            'samples_file_hash': 'deadbeef',
            'halogens': ['F', 'Cl'],
            'rules': ['R1'],
            'k_max': 2,
            'fast_mode': False,
            'sugar_masking_config': {'sugar_ring_score_threshold': 0.5},
            'rdkit_seed_report': {}
        },
        'overall_assessment': overall,
        'samples': {
            'sample-1': sample_entry
        }
    }
    return artifact


class TestReproSchemaCompatibility(unittest.TestCase):
    """Verify schema-aware comparison logic."""

    def test_ignores_missing_keys_when_not_strict(self):
        run1 = repro.extract_key_metrics(_make_sample('p1.1', include_new_fields=True))
        run2 = repro.extract_key_metrics(_make_sample('p1.0', include_new_fields=False))

        differences, ignored = repro.compare_results(run1, run2, strict=False)
        self.assertFalse(differences)
        self.assertIn('overall_assessment.accepted_via_score_count', ignored)

    def test_strict_mode_reports_missing_keys(self):
        run1 = repro.extract_key_metrics(_make_sample('p1.1', include_new_fields=True))
        run2 = repro.extract_key_metrics(_make_sample('p1.0', include_new_fields=False))

        differences, ignored = repro.compare_results(run1, run2, strict=True)
        self.assertTrue(any('missing in run2' in diff for diff in differences))
        self.assertEqual(ignored, [])

    def test_value_difference_detected(self):
        run1 = repro.extract_key_metrics(_make_sample('p1.1', include_new_fields=True))
        altered = _make_sample('p1.1', include_new_fields=True)
        altered['overall_assessment']['product_total'] = 101
        run2 = repro.extract_key_metrics(altered)

        differences, ignored = repro.compare_results(run1, run2, strict=False)
        self.assertTrue(any('product_total' in diff for diff in differences))
        self.assertFalse(ignored)

    def test_identical_metrics_pass(self):
        run1 = repro.extract_key_metrics(_make_sample('p1.1', include_new_fields=True))
        run2 = repro.extract_key_metrics(_make_sample('p1.1', include_new_fields=True))

        differences, ignored = repro.compare_results(run1, run2, strict=True)
        self.assertFalse(differences)
        self.assertFalse(ignored)

if __name__ == '__main__':
    unittest.main()
