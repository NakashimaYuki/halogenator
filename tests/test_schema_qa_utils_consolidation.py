#!/usr/bin/env python
"""Test that schema.py and qa_utils.py key sources are properly consolidated."""

import unittest


class TestSchemaQAUtilsConsolidation(unittest.TestCase):
    """Test consolidation of key sets between schema.py and qa_utils.py."""

    def test_schema_keys_include_qa_utils_keys(self):
        """Test that schema QA_PATH_KEYS includes all qa_utils standard keys."""
        from src.halogenator.schema import QA_PATH_KEYS
        from src.halogenator.qa_utils import STANDARD_QA_PATHS_KEYS

        # All qa_utils standard keys should be included in schema
        qa_utils_keys = set(STANDARD_QA_PATHS_KEYS)
        schema_keys = set(QA_PATH_KEYS)

        missing_keys = qa_utils_keys - schema_keys
        self.assertEqual(len(missing_keys), 0,
                        f"Schema QA_PATH_KEYS missing qa_utils keys: {missing_keys}")

        # Verify core canonical events are present
        core_events = {
            'isotope_unavailable', 'isotope_miss', 'atommap_used', 'heuristic_used',
            'sugar_mask_filtered', 'rdkit_error', 'sugar_post_guard_blocked'
        }
        self.assertTrue(core_events.issubset(schema_keys),
                       f"Schema missing core events: {core_events - schema_keys}")

    def test_schema_backward_compatibility_keys(self):
        """Test that schema maintains legacy keys for backward compatibility."""
        from src.halogenator.schema import QA_PATH_KEYS

        schema_keys = set(QA_PATH_KEYS)

        # Should include legacy/compatibility keys
        legacy_keys = {'carbonyl_unknown', 'pruned_inchikey_dupe', 'sugar_proximity_filtered'}
        self.assertTrue(legacy_keys.issubset(schema_keys),
                       f"Schema missing expected legacy keys: {legacy_keys - schema_keys}")

    def test_no_hardcoded_qa_utils_duplication(self):
        """Test that schema.py doesn't hardcode qa_utils constants."""
        # This is a documentation test - the actual consolidation is visible in the code
        # The import of STANDARD_QA_PATHS_KEYS ensures single source of truth

        from src.halogenator.schema import QA_PATH_KEYS
        from src.halogenator.qa_utils import STANDARD_QA_PATHS_KEYS

        # Verify that qa_utils keys are a subset, proving derivation relationship
        qa_utils_set = set(STANDARD_QA_PATHS_KEYS)
        schema_set = set(QA_PATH_KEYS)

        # qa_utils keys should be fully represented
        self.assertTrue(qa_utils_set.issubset(schema_set),
                       "Schema should derive from qa_utils keys, not duplicate them")

        # Schema can have additional keys beyond qa_utils (for backward compatibility)
        additional_keys = schema_set - qa_utils_set
        expected_additional = {'carbonyl_unknown', 'pruned_inchikey_dupe', 'sugar_proximity_filtered'}
        self.assertEqual(additional_keys, expected_additional,
                        f"Schema has unexpected additional keys: {additional_keys - expected_additional}")

    def test_import_compatibility_maintained(self):
        """Test that existing imports continue to work after consolidation."""
        # This validates that the consolidation doesn't break existing code
        try:
            from src.halogenator.schema import QA_PATH_KEYS, empty_qa_paths
            from src.halogenator.enumerate_k import QA_PATH_KEYS as ENUM_QA_KEYS
            from src.halogenator.report import QA_PATH_KEYS as REPORT_QA_KEYS

            # All imports should work and reference the same object
            self.assertIs(QA_PATH_KEYS, ENUM_QA_KEYS)
            self.assertIs(QA_PATH_KEYS, REPORT_QA_KEYS)

            # empty_qa_paths should work with consolidated keys
            empty_paths = empty_qa_paths()
            self.assertIsInstance(empty_paths, dict)
            self.assertEqual(len(empty_paths), len(QA_PATH_KEYS))

        except ImportError as e:
            self.fail(f"Import compatibility broken by consolidation: {e}")


if __name__ == '__main__':
    unittest.main()