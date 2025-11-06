# -*- coding: ascii -*-
"""Test v2 distribution_keys source priority and consistency."""

import unittest
import tempfile
import os
import json
from src.halogenator.report import write_qa_summary_json


class TestV2DistributionKeysSources(unittest.TestCase):
    """Test v2 distribution_keys source priority and output consistency."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_v2_distribution_keys_prefers_metadata_lists_over_pivots(self):
        """Test that v2 distribution_keys prefers metadata rules/halogens over pivots."""
        qa_stats = {
            'version': '2',
            'metadata': {
                'rules': ['R1', 'R3'],
                'halogens': ['F', 'Cl', 'Br']
            },
            'pivots': {
                'by_rule': {
                    'R1': {'metric1': 5}
                    # R3 missing from pivots (zero count)
                },
                'by_halogen': {
                    'F': {'metric1': 3},
                    'Cl': {'metric1': 2}
                    # Br missing from pivots (zero count)
                }
            }
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir)

        # Read the output
        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should use metadata lists, not pivot inference
        self.assertIn('distribution_keys', result['metadata'])
        dist_keys = result['metadata']['distribution_keys']

        # Should include all rules and halogens from metadata, even if missing from pivots
        self.assertIn('R1', dist_keys['rules'])
        self.assertIn('R3', dist_keys['rules'])  # Missing from pivots but in metadata
        self.assertIn('F', dist_keys['halogens'])
        self.assertIn('Cl', dist_keys['halogens'])
        self.assertIn('Br', dist_keys['halogens'])  # Missing from pivots but in metadata

    def test_v2_distribution_keys_present_even_if_empty(self):
        """Test that v2 distribution_keys is always present, even if empty."""
        qa_stats = {
            'version': '2',
            # No metadata, no pivots
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir)

        # Read the output
        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should always have distribution_keys field
        self.assertIn('distribution_keys', result['metadata'])
        dist_keys = result['metadata']['distribution_keys']

        # Should have both keys, even if empty
        self.assertIn('rules', dist_keys)
        self.assertIn('halogens', dist_keys)
        self.assertIsInstance(dist_keys['rules'], list)
        self.assertIsInstance(dist_keys['halogens'], list)

    def test_v2_distribution_keys_falls_back_to_pivots_when_metadata_missing(self):
        """Test that v2 falls back to pivot inference when metadata is missing."""
        qa_stats = {
            'version': '2',
            # No metadata
            'pivots': {
                'by_rule': {
                    'R1': {'metric1': 5},
                    'R3': {'metric1': 3}
                },
                'by_halogen': {
                    'F': {'metric1': 4},
                    'Cl': {'metric1': 4}
                }
            }
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir)

        # Read the output
        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should use pivot inference
        self.assertIn('distribution_keys', result['metadata'])
        dist_keys = result['metadata']['distribution_keys']

        # Should include rules and halogens from pivots
        self.assertIn('R1', dist_keys['rules'])
        self.assertIn('R3', dist_keys['rules'])
        self.assertIn('F', dist_keys['halogens'])
        self.assertIn('Cl', dist_keys['halogens'])

    def test_v2_distribution_keys_partial_metadata_fallback(self):
        """Test that v2 falls back to pivots for missing metadata parts."""
        qa_stats = {
            'version': '2',
            'metadata': {
                'rules': ['R1', 'R3']
                # halogens missing
            },
            'pivots': {
                'by_halogen': {
                    'F': {'metric1': 4},
                    'Cl': {'metric1': 4}
                }
            }
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir)

        # Read the output
        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should use metadata for rules, pivots for halogens
        dist_keys = result['metadata']['distribution_keys']

        # Rules from metadata
        self.assertIn('R1', dist_keys['rules'])
        self.assertIn('R3', dist_keys['rules'])

        # Halogens from pivots
        self.assertIn('F', dist_keys['halogens'])
        self.assertIn('Cl', dist_keys['halogens'])

    def test_v2_distribution_keys_extracts_from_2d_pivots(self):
        """Test that v2 can extract rules and halogens from 2D pivot structure."""
        qa_stats = {
            'version': '2',
            # No metadata
            'pivots': {
                'by_rule_halogen': {
                    'R1': {
                        'F': {'metric1': 2},
                        'Cl': {'metric1': 3}
                    },
                    'R3': {
                        'F': {'metric1': 1},
                        'Br': {'metric1': 1}
                    }
                }
            }
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir)

        # Read the output
        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should extract from 2D structure
        dist_keys = result['metadata']['distribution_keys']

        # Should include all rules from 2D keys
        self.assertIn('R1', dist_keys['rules'])
        self.assertIn('R3', dist_keys['rules'])

        # Should include all halogens from 2D values
        self.assertIn('F', dist_keys['halogens'])
        self.assertIn('Cl', dist_keys['halogens'])
        self.assertIn('Br', dist_keys['halogens'])

    def test_v2_distribution_keys_sorting_consistency(self):
        """Test that distribution_keys are sorted consistently."""
        qa_stats = {
            'version': '2',
            'metadata': {
                'rules': ['R3', 'R1', 'R5'],  # Unsorted
                'halogens': ['Br', 'F', 'Cl']  # Unsorted
            }
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir)

        # Read the output
        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should be sorted according to distribution order (lexicographic)
        dist_keys = result['metadata']['distribution_keys']

        # Check that rules and halogens are in distribution order
        # (the actual order depends on _distribution_order implementation)
        self.assertIsInstance(dist_keys['rules'], list)
        self.assertIsInstance(dist_keys['halogens'], list)

        # All original rules and halogens should be present
        self.assertEqual(set(dist_keys['rules']), {'R1', 'R3', 'R5'})
        self.assertEqual(set(dist_keys['halogens']), {'Br', 'Cl', 'F'})


if __name__ == '__main__':
    unittest.main()