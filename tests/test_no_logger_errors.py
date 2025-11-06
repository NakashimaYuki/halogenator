# -*- coding: ascii -*-
"""
Tests to enforce no UnboundLocalError and proper R3 rule functionality.

These tests ensure that logger scoping issues are resolved and that
R3 rules work correctly without being short-circuited by logger errors.
"""

import unittest
import io
import sys
import contextlib
from typing import List, Dict, Any

# Import the modules we want to test
try:
    from rdkit import Chem
    from src.halogenator.enumerate_k import enumerate_products, EnumConfig
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


@contextlib.contextmanager
def capture_stderr():
    """Capture stderr to check for logger errors."""
    old_stderr = sys.stderr
    sys.stderr = captured_stderr = io.StringIO()
    try:
        yield captured_stderr
    finally:
        sys.stderr = old_stderr


@unittest.skipUnless(RDKIT_AVAILABLE, "RDKit not available")
class TestNoLoggerErrors(unittest.TestCase):
    """Test that logger errors are resolved and R3 rules work correctly."""

    def test_r3_rule_no_logger_error(self):
        """Test that R3 rule works without UnboundLocalError."""
        # Simple phenolic substrate where R3 should have viable matches
        phenol_smiles = "c1ccc(O)cc1"

        # Configure to only test R3 with single halogen
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R3',),
            sugar_cfg={'mode': 'off'}
        )

        # Capture stderr to check for logger errors
        with capture_stderr() as captured:
            results = list(enumerate_products(phenol_smiles, cfg, return_qa_stats=True))
            stderr_content = captured.getvalue()

        # Assert: no logger assignment errors
        self.assertNotIn("local variable 'logger' referenced before assignment", stderr_content,
                        "Logger scoping error detected in R3 rule execution")
        self.assertNotIn("local variable 'logging' referenced before assignment", stderr_content,
                        "Logging import scoping error detected in R3 rule execution")

        # Assert: enumeration completed successfully (results may be empty)
        self.assertIsInstance(results, list, "R3 enumeration should return a list")

        # Check QA statistics are populated correctly
        if results:
            qa_stats = results[-1][1] if len(results[-1]) > 1 else {}
            attempts = qa_stats.get('attempts', 0)
            # R3 may or may not have attempts depending on the substrate, but should not crash
            self.assertIsInstance(attempts, int, "Attempts should be an integer")

            # Check that sugar statistics keys are present and initialized
            qa_paths = qa_stats.get('qa_paths', {})
            self.assertIn('sugar_filtered_reaction_matches', qa_paths,
                         "sugar_filtered_reaction_matches should be in qa_paths")
            self.assertIn('sugar_post_guard_blocked', qa_paths,
                         "sugar_post_guard_blocked should be in qa_paths")

            # For non-sugar molecule with sugar masking off, these should be 0
            self.assertEqual(qa_paths['sugar_filtered_reaction_matches'], 0,
                           "Non-sugar molecule should have no sugar filtering")
            self.assertEqual(qa_paths['sugar_post_guard_blocked'], 0,
                           "Non-sugar molecule should have no sugar blocking")

    def test_sugar_masking_reduces_matches(self):
        """Test that sugar masking reduces reaction matches for glycosides."""
        # Simple glycoside-like molecule
        glucose_smiles = "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"

        # Test with sugar masking off
        cfg_no_mask = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R3',),
            sugar_cfg={'mode': 'off'}
        )

        # Test with sugar masking on
        cfg_with_mask = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R3',),
            sugar_cfg={'mode': 'heuristic'}
        )

        # Capture stderr for both runs
        with capture_stderr() as captured_no_mask:
            results_no_mask = list(enumerate_products(glucose_smiles, cfg_no_mask, return_qa_stats=True))
            stderr_no_mask = captured_no_mask.getvalue()

        with capture_stderr() as captured_with_mask:
            results_with_mask = list(enumerate_products(glucose_smiles, cfg_with_mask, return_qa_stats=True))
            stderr_with_mask = captured_with_mask.getvalue()

        # Assert: no logger errors in either case
        self.assertNotIn("local variable 'logger' referenced before assignment", stderr_no_mask,
                        "Logger error in no-mask enumeration")
        self.assertNotIn("local variable 'logger' referenced before assignment", stderr_with_mask,
                        "Logger error in with-mask enumeration")

        # Extract QA statistics
        qa_no_mask = results_no_mask[-1][1] if results_no_mask and len(results_no_mask[-1]) > 1 else {}
        qa_with_mask = results_with_mask[-1][1] if results_with_mask and len(results_with_mask[-1]) > 1 else {}

        attempts_no_mask = qa_no_mask.get('attempts', 0)
        attempts_with_mask = qa_with_mask.get('attempts', 0)

        print(f"Glucose enumeration - no mask: {attempts_no_mask} attempts, with mask: {attempts_with_mask} attempts")

        # For pure sugar molecule, masking should reduce or eliminate attempts
        self.assertGreaterEqual(attempts_no_mask, attempts_with_mask,
                               "Sugar masking should reduce or maintain attempt count for sugar molecules")

        # Check that sugar filtering statistics are recorded when masking is enabled
        if results_with_mask:
            qa_paths = qa_with_mask.get('qa_paths', {})
            # Sugar filtering count should be present (may be 0 for this specific molecule/rules)
            self.assertIsInstance(qa_paths.get('sugar_filtered_reaction_matches', 0), int,
                                "sugar_filtered_reaction_matches should be integer")

    def test_qa_stats_schema_completeness_no_errors(self):
        """Test QA statistics schema completeness without logger errors."""
        # Simple molecule for testing
        test_smiles = "c1ccc(O)cc1"

        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R3', 'R4'),
            sugar_cfg={'mode': 'heuristic'}
        )

        with capture_stderr() as captured:
            results = list(enumerate_products(test_smiles, cfg, return_qa_stats=True))
            stderr_content = captured.getvalue()

        # Assert: no logger errors
        self.assertNotIn("local variable 'logger' referenced before assignment", stderr_content,
                        "Logger scoping errors should be eliminated")

        # Check QA schema completeness
        if results:
            qa_stats = results[-1][1] if len(results[-1]) > 1 else {}
            qa_paths = qa_stats.get('qa_paths', {})

            # All expected sugar statistics should be present
            expected_sugar_keys = ['sugar_filtered_reaction_matches', 'sugar_post_guard_blocked']
            for key in expected_sugar_keys:
                self.assertIn(key, qa_paths, f"Missing QA statistic: {key}")
                self.assertIsInstance(qa_paths[key], int, f"QA statistic {key} should be integer")

            # Standard QA path keys should also be present
            standard_keys = ['isotope_unavailable', 'isotope_miss', 'atommap_used', 'heuristic_used']
            for key in standard_keys:
                self.assertIn(key, qa_paths, f"Missing standard QA statistic: {key}")

    def test_all_rules_no_logger_errors(self):
        """Test that all rules work without logger errors."""
        test_smiles = "c1ccc(O)cc1"

        # Test each rule individually to isolate any remaining issues
        for rule in ['R1', 'R2', 'R3', 'R4']:
            with self.subTest(rule=rule):
                cfg = EnumConfig(
                    k_max=1,
                    halogens=('F',),
                    rules=(rule,),
                    sugar_cfg={'mode': 'off'}
                )

                with capture_stderr() as captured:
                    results = list(enumerate_products(test_smiles, cfg, return_qa_stats=True))
                    stderr_content = captured.getvalue()

                # Assert: no logger errors for any rule
                self.assertNotIn("local variable 'logger' referenced before assignment", stderr_content,
                                f"Logger error detected in rule {rule}")
                self.assertNotIn("local variable 'logging' referenced before assignment", stderr_content,
                                f"Logging import error detected in rule {rule}")

                # Each rule should be able to run (may or may not produce results depending on substrate)
                self.assertIsInstance(results, list, f"Rule {rule} should return a list")


if __name__ == '__main__':
    unittest.main()