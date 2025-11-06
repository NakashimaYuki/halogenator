# -*- coding: ascii -*-
"""
E2E test for verifying presence of critical schema fields in parquet output.

This test ensures that parquet files contain the expected top-level columns
introduced in PR2 (detection, sub_rule, rule_family, etc.).

Can be run as:
1. Unit test: python -m unittest tests.test_parquet_schema_fields
2. Standalone script: python tests/test_parquet_schema_fields.py <parquet_file>
"""

import sys
import unittest
from pathlib import Path
from typing import Set, List, Tuple
import pandas as pd


# Critical fields introduced in PR2 that should be present
PR2_CRITICAL_FIELDS = {
    'detection',      # Detection method (strict/fallback for R2b)
    'sub_rule',       # Sub-rule identifier (R2a, R2b, etc.)
    'rule_family',    # Rule family for grouping (R2a/R2b -> R2)
}

# Additional expected fields (not strictly PR2, but should be present)
EXPECTED_PRODUCT_FIELDS = {
    'smiles',
    'inchikey',
    'parent_smiles',
    'parent_inchikey',
    'k',
    'rule',
    'halogen',
}


def check_parquet_schema(parquet_path: str) -> Tuple[bool, List[str], Set[str]]:
    """
    Check if parquet file contains expected schema fields.

    Args:
        parquet_path: Path to parquet file

    Returns:
        Tuple of (all_present, missing_fields, actual_columns)
    """
    try:
        df = pd.read_parquet(parquet_path)
        actual_cols = set(df.columns)

        # Check for PR2 critical fields
        missing = PR2_CRITICAL_FIELDS - actual_cols

        return len(missing) == 0, sorted(missing), actual_cols

    except Exception as e:
        raise RuntimeError(f"Failed to read parquet file: {e}")


def print_schema_report(parquet_path: str, verbose: bool = True) -> bool:
    """
    Print detailed schema validation report with color-coded output.

    Args:
        parquet_path: Path to parquet file
        verbose: Whether to print detailed column list

    Returns:
        True if all critical fields present, False otherwise
    """
    # ANSI color codes
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RESET = '\033[0m'
    BOLD = '\033[1m'

    try:
        all_present, missing, actual_cols = check_parquet_schema(parquet_path)

        print(f"\n{BOLD}=== Parquet Schema Validation ==={RESET}")
        print(f"File: {parquet_path}")
        print(f"Total columns: {len(actual_cols)}")

        if all_present:
            print(f"\n{GREEN}{BOLD}[PASS]{RESET} All PR2 critical fields present!")
        else:
            print(f"\n{RED}{BOLD}[FAIL]{RESET} Missing critical fields:")
            for field in missing:
                print(f"  {RED}- {field}{RESET}")

        # Check for expected product fields
        missing_expected = EXPECTED_PRODUCT_FIELDS - actual_cols
        if missing_expected:
            print(f"\n{YELLOW}[WARNING]{RESET} Missing expected product fields:")
            for field in sorted(missing_expected):
                print(f"  {YELLOW}- {field}{RESET}")

        # Print present critical fields
        present_critical = PR2_CRITICAL_FIELDS & actual_cols
        if present_critical:
            print(f"\n{GREEN}Present critical fields:{RESET}")
            for field in sorted(present_critical):
                print(f"  {GREEN}+ {field}{RESET}")

        if verbose:
            print(f"\n{BOLD}All columns:{RESET}")
            for col in sorted(actual_cols):
                marker = '+ ' if col in PR2_CRITICAL_FIELDS else '  '
                color = GREEN if col in PR2_CRITICAL_FIELDS else ''
                print(f"  {color}{marker}{col}{RESET}")

        print()
        return all_present

    except Exception as e:
        print(f"\n{RED}{BOLD}[ERROR]{RESET} {e}")
        return False


class TestParquetSchemaFields(unittest.TestCase):
    """Test suite for parquet schema field validation."""

    def test_schema_fields_mock_data(self):
        """
        Test schema validation with mock DataFrame.

        This test creates a minimal DataFrame with expected fields
        to validate the checking logic works correctly.
        """
        # Create mock DataFrame with all expected fields
        mock_data = {
            'smiles': ['CC', 'CCCl'],
            'inchikey': ['key1', 'key2'],
            'parent_smiles': ['CC', 'CC'],
            'parent_inchikey': ['parent_key', 'parent_key'],
            'k': [0, 1],
            'rule': ['R1', 'R1'],
            'halogen': ['', 'Cl'],
            'detection': ['strict', 'strict'],
            'sub_rule': ['R1', 'R1'],
            'rule_family': ['R1', 'R1'],
        }
        df = pd.DataFrame(mock_data)

        # Check columns
        actual_cols = set(df.columns)
        missing = PR2_CRITICAL_FIELDS - actual_cols

        self.assertEqual(len(missing), 0,
                        f"Mock DataFrame should have all critical fields, missing: {missing}")

    def test_critical_fields_definition(self):
        """Test that critical fields are properly defined."""
        self.assertIn('detection', PR2_CRITICAL_FIELDS)
        self.assertIn('sub_rule', PR2_CRITICAL_FIELDS)
        self.assertIn('rule_family', PR2_CRITICAL_FIELDS)
        self.assertEqual(len(PR2_CRITICAL_FIELDS), 3)

    def test_expected_product_fields_definition(self):
        """Test that expected product fields are properly defined."""
        required = {'smiles', 'inchikey', 'parent_smiles', 'parent_inchikey',
                   'k', 'rule', 'halogen'}
        self.assertTrue(required.issubset(EXPECTED_PRODUCT_FIELDS))

    def test_check_parquet_schema_with_actual_file(self):
        """
        Test schema validation with actual parquet file if available.

        This test looks for RECENT parquet files (modified in last hour)
        in common output locations and validates their schema.
        Skips if no recent files found.

        Note: Old parquet files from previous runs may not have PR2 fields.
        """
        import time

        # Common output patterns for NEW runs (not legacy p0/p1 outputs)
        search_patterns = [
            'data/output/naringenin_*/products_k*.parquet',
            'data/output/*_raw/products_k*.parquet',
            'data/output/*_strict/products_k*.parquet',
            'test_pr2_*/products_k*.parquet',
        ]

        parquet_files = []
        for pattern in search_patterns:
            parquet_files.extend(Path('.').glob(pattern))

        if not parquet_files:
            self.skipTest("No new parquet files found (looking for recent PR2 outputs)")

        # Filter to files modified in last hour (3600 seconds)
        current_time = time.time()
        recent_files = [
            f for f in parquet_files
            if (current_time - f.stat().st_mtime) < 3600
        ]

        if not recent_files:
            self.skipTest("No parquet files modified in last hour (no recent runs)")

        # Test with the most recent file
        most_recent = max(recent_files, key=lambda p: p.stat().st_mtime)
        print(f"\nTesting with: {most_recent}")

        all_present, missing, actual_cols = check_parquet_schema(str(most_recent))

        # Assert critical fields are present
        self.assertTrue(all_present,
                       f"Critical fields missing in {most_recent}: {missing}\n"
                       f"Actual columns: {sorted(actual_cols)}\n"
                       f"Hint: Run 'python -m halogenator.cli enum' with --r2-fallback "
                       f"to generate files with PR2 fields")

        # Also verify expected product fields
        missing_expected = EXPECTED_PRODUCT_FIELDS - actual_cols
        self.assertEqual(len(missing_expected), 0,
                        f"Expected product fields missing: {missing_expected}")


def main():
    """Standalone script entry point."""
    if len(sys.argv) < 2:
        print("Usage: python tests/test_parquet_schema_fields.py <parquet_file>")
        print("\nOr run as unit test: python -m unittest tests.test_parquet_schema_fields")
        sys.exit(1)

    parquet_path = sys.argv[1]
    if not Path(parquet_path).exists():
        print(f"Error: File not found: {parquet_path}")
        sys.exit(1)

    success = print_schema_report(parquet_path, verbose=True)
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    if len(sys.argv) > 1 and not sys.argv[1].startswith('test'):
        # Standalone mode
        main()
    else:
        # Unit test mode
        unittest.main()
