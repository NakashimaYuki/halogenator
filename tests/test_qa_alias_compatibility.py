# -*- coding: ascii -*-
"""
Test QA key alias compatibility.

This test ensures that legacy alias keys are properly merged into standard keys
for backward compatibility with existing scripts and data.
"""

import unittest
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from halogenator.report import _normalize_qa_paths_with_aliases
from halogenator.schema import ensure_qa_paths_compatibility
from halogenator.enumerate_k import EnumConfig, enumerate_with_stats


class TestQAAliasCompatibility(unittest.TestCase):
    """Test QA key alias compatibility and merging."""

    def test_normalize_qa_paths_with_legacy_keys_only(self):
        """Test normalization when only legacy keys are present."""
        qa_paths = {
            'sugar_filtered_reaction_matches': 5,
            'post_guard_blocked': 3,
            'isotope_miss': 2
        }

        normalized = _normalize_qa_paths_with_aliases(qa_paths)

        expected = {
            'sugar_mask_filtered': 5,
            'sugar_post_guard_blocked': 3,
            'isotope_miss': 2
        }

        self.assertEqual(normalized, expected)

    def test_normalize_qa_paths_with_standard_keys_only(self):
        """Test normalization when only standard keys are present."""
        qa_paths = {
            'sugar_mask_filtered': 7,
            'sugar_post_guard_blocked': 4,
            'isotope_miss': 1
        }

        normalized = _normalize_qa_paths_with_aliases(qa_paths)

        # Should remain unchanged
        self.assertEqual(normalized, qa_paths)

    def test_normalize_qa_paths_with_mixed_keys(self):
        """Test normalization when both legacy and standard keys are present."""
        qa_paths = {
            'sugar_filtered_reaction_matches': 3,  # Legacy
            'sugar_mask_filtered': 4,              # Standard
            'post_guard_blocked': 2,               # Legacy
            'sugar_post_guard_blocked': 1,         # Standard
            'isotope_miss': 5                      # Standard (no alias)
        }

        normalized = _normalize_qa_paths_with_aliases(qa_paths)

        expected = {
            'sugar_mask_filtered': 7,              # 3 + 4 (merged)
            'sugar_post_guard_blocked': 3,         # 2 + 1 (merged)
            'isotope_miss': 5
        }

        self.assertEqual(normalized, expected)

    def test_normalize_qa_paths_empty_dict(self):
        """Test normalization with empty dictionary."""
        qa_paths = {}
        normalized = _normalize_qa_paths_with_aliases(qa_paths)
        self.assertEqual(normalized, {})

    def test_normalize_qa_paths_no_aliases(self):
        """Test normalization when no alias keys are present."""
        qa_paths = {
            'isotope_unavailable': 2,
            'atommap_used': 1,
            'heuristic_used': 3
        }

        normalized = _normalize_qa_paths_with_aliases(qa_paths)

        # Should remain unchanged
        self.assertEqual(normalized, qa_paths)

    def test_normalized_qa_paths_preserves_original_dict(self):
        """Test that normalization doesn't modify the original dictionary."""
        original_qa_paths = {
            'sugar_filtered_reaction_matches': 5,
            'isotope_miss': 2
        }

        # Keep a copy for comparison
        original_copy = dict(original_qa_paths)

        normalized = _normalize_qa_paths_with_aliases(original_qa_paths)

        # Original should be unchanged
        self.assertEqual(original_qa_paths, original_copy)

        # Normalized should have the expected changes
        expected = {
            'sugar_mask_filtered': 5,
            'isotope_miss': 2
        }
        self.assertEqual(normalized, expected)

    def test_ensure_qa_paths_compatibility_legacy_off(self):
        """Test ensure_qa_paths_compatibility with emit_legacy_keys=False."""
        qa_paths = {
            'dedup_hits_inchi': 5,
            'isotope_unavailable': 2,
            'atommap_used': 1
        }

        result = ensure_qa_paths_compatibility(qa_paths, emit_legacy_keys=False)

        # Should contain the new key
        self.assertIn('dedup_hits_inchi', result)
        self.assertEqual(result['dedup_hits_inchi'], 5)

        # Legacy key may be absent or equal to new key if present
        if 'pruned_inchikey_dupe' in result:
            self.assertEqual(result['pruned_inchikey_dupe'], result['dedup_hits_inchi'])

    def test_ensure_qa_paths_compatibility_legacy_on(self):
        """Test ensure_qa_paths_compatibility with emit_legacy_keys=True."""
        qa_paths = {
            'dedup_hits_inchi': 3,
            'isotope_unavailable': 1
        }

        result = ensure_qa_paths_compatibility(qa_paths, emit_legacy_keys=True)

        # Should contain both keys
        self.assertIn('dedup_hits_inchi', result)
        self.assertIn('pruned_inchikey_dupe', result)

        # Legacy key should equal new key
        self.assertEqual(result['pruned_inchikey_dupe'], result['dedup_hits_inchi'])

    def test_ensure_qa_paths_compatibility_zero_values(self):
        """Test ensure_qa_paths_compatibility with zero dedup_hits_inchi."""
        qa_paths = {
            'dedup_hits_inchi': 0,
            'isotope_unavailable': 2
        }

        # Test with legacy keys on
        result_on = ensure_qa_paths_compatibility(qa_paths, emit_legacy_keys=True)
        self.assertEqual(result_on['dedup_hits_inchi'], 0)
        # Legacy key behavior with zero should be consistent
        if 'pruned_inchikey_dupe' in result_on:
            self.assertEqual(result_on['pruned_inchikey_dupe'], 0)

        # Test with legacy keys off
        result_off = ensure_qa_paths_compatibility(qa_paths, emit_legacy_keys=False)
        self.assertEqual(result_off['dedup_hits_inchi'], 0)

    def test_enumerate_with_stats_legacy_alias_control(self):
        """Test that enumerate_with_stats respects engine_cfg.emit_legacy_keys."""
        # Test with emit_legacy_keys=False
        cfg_off = EnumConfig(k_max=1, halogens=('F',), rules=('R6',),
                             engine_cfg={'budget_mode': 'ops', 'emit_legacy_keys': False})
        _, qa_off = enumerate_with_stats('Cc1ccc(C)cc1', cfg_off)
        qp_off = qa_off.get('qa_paths', {})

        self.assertIn('dedup_hits_inchi', qp_off)
        # Legacy key may be absent or equal to new key
        if 'pruned_inchikey_dupe' in qp_off:
            self.assertEqual(qp_off['pruned_inchikey_dupe'], qp_off['dedup_hits_inchi'])

        # Test with emit_legacy_keys=True
        cfg_on = EnumConfig(k_max=1, halogens=('F',), rules=('R6',),
                            engine_cfg={'budget_mode': 'ops', 'emit_legacy_keys': True})
        _, qa_on = enumerate_with_stats('Cc1ccc(C)cc1', cfg_on)
        qp_on = qa_on.get('qa_paths', {})

        self.assertIn('dedup_hits_inchi', qp_on)
        self.assertIn('pruned_inchikey_dupe', qp_on)
        self.assertEqual(qp_on['pruned_inchikey_dupe'], qp_on['dedup_hits_inchi'])


if __name__ == '__main__':
    unittest.main(verbosity=2)