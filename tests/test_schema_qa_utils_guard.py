# -*- coding: ascii -*-
"""
Schema-qa_utils consolidation guard tests.

Ensures that schema.py properly derives from qa_utils.py as the single source of truth
for event key definitions, preventing dual-source-of-truth regressions.
"""

import unittest
from src.halogenator.schema import QA_PATH_KEYS
from src.halogenator.qa_utils import STANDARD_QA_PATHS_KEYS


class TestSchemaQAUtilsGuard(unittest.TestCase):
    """Guard against schema-qa_utils consolidation regressions."""

    def test_schema_keys_superset_of_qa_utils(self):
        """Test that schema QA_PATH_KEYS includes all qa_utils standard keys."""
        qa_utils_set = set(STANDARD_QA_PATHS_KEYS)
        schema_set = set(QA_PATH_KEYS)

        missing_keys = qa_utils_set - schema_set
        self.assertEqual(len(missing_keys), 0,
                        f"Schema QA_PATH_KEYS missing qa_utils keys: {sorted(missing_keys)}")

    def test_schema_derives_from_qa_utils(self):
        """Test that schema properly imports and derives from qa_utils constants."""
        # Schema should include all standard keys
        qa_utils_set = set(STANDARD_QA_PATHS_KEYS)
        schema_set = set(QA_PATH_KEYS)

        self.assertTrue(qa_utils_set.issubset(schema_set),
                       "Schema should derive from qa_utils keys, not duplicate them")

    def test_schema_only_adds_specific_legacy_keys(self):
        """Test that schema only adds expected additional keys beyond qa_utils."""
        qa_utils_set = set(STANDARD_QA_PATHS_KEYS)
        schema_set = set(QA_PATH_KEYS)

        additional_keys = schema_set - qa_utils_set
        expected_additional = {
            'carbonyl_unknown',           # Schema-specific legacy key (used in tests)
            'pruned_inchikey_dupe',       # Legacy alias maintained for backward compatibility
            'sugar_proximity_filtered',   # Non-pivot path event (used in tests)
        }

        self.assertEqual(additional_keys, expected_additional,
                        f"Schema has unexpected additional keys: {additional_keys - expected_additional}")

    def test_no_hardcoded_event_constants_in_schema(self):
        """Test that schema.py doesn't hardcode event constants that should come from qa_utils."""
        import inspect
        import src.halogenator.schema as schema_module

        # Get the source code of schema.py
        schema_source = inspect.getsource(schema_module)

        # Check that schema doesn't define its own event constants
        # (these should come from qa_utils imports)
        forbidden_patterns = [
            "EV_ISOTOPE_UNAVAILABLE",
            "EV_SUGAR_MASK_FILTERED",
            "EV_RDKIT_ERROR",
            "EV_DEDUP_INCHI",
            "EV_SUGAR_POST_GUARD"
        ]

        violations = []
        for pattern in forbidden_patterns:
            if pattern in schema_source and f"= '{pattern.lower().replace('ev_', '').replace('_', '_')}'" in schema_source:
                violations.append(pattern)

        if violations:
            self.fail(f"Schema hardcodes event constants that should come from qa_utils: {violations}")

    def test_qa_utils_import_present(self):
        """Test that schema.py properly imports from qa_utils."""
        import src.halogenator.schema as schema_module
        import inspect

        schema_source = inspect.getsource(schema_module)

        # Should have the consolidation import
        self.assertIn("from .qa_utils import STANDARD_QA_PATHS_KEYS", schema_source,
                     "Schema should import STANDARD_QA_PATHS_KEYS from qa_utils")


if __name__ == '__main__':
    unittest.main()