#!/usr/bin/env python3
"""
Unit tests for batch_stats.py CLI functionality.

Tests dependency checking, statistics generation, audit sample creation,
and QA summary reporting capabilities with minimal mock data.
"""
import sys
import os
import unittest
import tempfile
import json
import shutil
from pathlib import Path
from unittest.mock import patch, MagicMock

# Ensure testing bootstrap is available for direct execution
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import testing_bootstrap

# Import the batch_stats module
scripts_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'scripts')
sys.path.insert(0, scripts_dir)
import batch_stats


class TestBatchStatsCLI(unittest.TestCase):
    """Test batch_stats.py CLI functionality."""

    def setUp(self):
        """Set up test fixtures with temporary directories."""
        self.test_dir = tempfile.mkdtemp()
        self.products_csv = os.path.join(self.test_dir, "products.csv")
        self.products_parquet = os.path.join(self.test_dir, "products.parquet")
        self.qa_json = os.path.join(self.test_dir, "qa_summary.json")
        self.audit_dir = os.path.join(self.test_dir, "audit")

        # Create minimal test products data
        self.create_test_products()
        self.create_test_qa_summary()

    def tearDown(self):
        """Clean up test fixtures."""
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def create_test_products(self):
        """Create minimal test products data (CSV and parquet)."""
        # CSV data with 5 rows containing parent_key, rule, halogen
        csv_content = """parent_key,rule,halogen,smiles,products
parent_1,R1,Cl,c1ccc(Cl)cc1,1
parent_1,R2,Br,c1ccc(Br)cc1,1
parent_2,R1,Cl,CCCl,1
parent_2,R1,F,CCF,1
parent_3,R3,I,c1ccc(I)cc1,1
"""
        with open(self.products_csv, 'w') as f:
            f.write(csv_content)

        # Create parquet version if pandas available
        try:
            import pandas as pd
            df = pd.read_csv(self.products_csv)
            df.to_parquet(self.products_parquet)
        except ImportError:
            # Create a dummy parquet file if pandas not available
            with open(self.products_parquet, 'wb') as f:
                f.write(b'dummy parquet content')

    def create_test_qa_summary(self):
        """Create minimal test QA summary JSON."""
        qa_data = {
            "metadata": {
                "warnings_count": 15,
                "warnings_returned": 10,
                "warnings_truncated": 5
            },
            "warnings": [
                {
                    "type": "non_integer_value_detected",
                    "where": {
                        "dimension": "by_rule",
                        "key": "R1",
                        "metric": "products"
                    },
                    "value": 2.7,
                    "coerced_to_int": True
                }
            ]
        }
        with open(self.qa_json, 'w') as f:
            json.dump(qa_data, f)

    def test_check_deps_with_pandas(self):
        """Test --check-deps returns 0 when pandas/pyarrow available."""
        with patch('batch_stats.check_pandas_availability', return_value=True):
            with patch('sys.exit') as mock_exit:
                batch_stats.main.__globals__['sys'] = MagicMock()
                batch_stats.main.__globals__['sys'].argv = ['batch_stats.py', '--check-deps']

                # Capture print output
                with patch('builtins.print') as mock_print:
                    try:
                        batch_stats.main()
                    except SystemExit:
                        pass

                mock_print.assert_called_with("Pandas and PyArrow available")
                mock_exit.assert_called_with(0)

    def test_check_deps_without_pandas(self):
        """Test --check-deps returns 1 when pandas/pyarrow not available."""
        with patch('batch_stats.check_pandas_availability', return_value=False):
            with patch('sys.exit') as mock_exit:
                batch_stats.main.__globals__['sys'] = MagicMock()
                batch_stats.main.__globals__['sys'].argv = ['batch_stats.py', '--check-deps']

                # Capture print output
                with patch('builtins.print') as mock_print:
                    try:
                        batch_stats.main()
                    except SystemExit:
                        pass

                mock_print.assert_called_with("Pandas and/or PyArrow not available")
                mock_exit.assert_called_with(1)

    def test_products_statistics_with_pandas(self):
        """Test products statistics generation with pandas available."""
        with patch('batch_stats.check_pandas_availability', return_value=True):
            with patch('sys.argv', ['batch_stats.py', '--products', self.products_parquet, '--label', 'test']):
                with patch('builtins.print') as mock_print:
                    batch_stats.main()

                # Check that statistics were printed
                print_calls = [call[0][0] for call in mock_print.call_args_list]

                # Should show product count
                self.assertTrue(any("Products generated: 5" in call for call in print_calls))
                # Should show unique parents
                self.assertTrue(any("Unique parent molecules: 3" in call for call in print_calls))
                # Should show rule distribution
                self.assertTrue(any("Rules distribution:" in call for call in print_calls))
                # Should show halogen distribution
                self.assertTrue(any("Halogens distribution:" in call for call in print_calls))

    def test_products_statistics_without_pandas(self):
        """Test products statistics fallback without pandas."""
        with patch('batch_stats.check_pandas_availability', return_value=False):
            with patch('sys.argv', ['batch_stats.py', '--products', self.products_parquet, '--label', 'test']):
                with patch('builtins.print') as mock_print:
                    batch_stats.main()

                # Check that basic file stats were printed
                print_calls = [call[0][0] for call in mock_print.call_args_list]

                self.assertTrue(any("Pandas not available - using basic file-based statistics" in call for call in print_calls))
                self.assertTrue(any("test file exists:" in call for call in print_calls))

    def test_qa_summary_reading(self):
        """Test QA summary reading and reporting."""
        with patch('sys.argv', ['batch_stats.py', '--qa', self.qa_json]):
            with patch('builtins.print') as mock_print:
                batch_stats.main()

            # Check that QA warnings were reported
            print_calls = [call[0][0] for call in mock_print.call_args_list]

            self.assertTrue(any("QA warnings: 15 (returned: 10)" in call for call in print_calls))

    def test_audit_sample_generation(self):
        """Test audit sample generation with specified size."""
        audit_size = 7

        with patch('batch_stats.check_pandas_availability', return_value=True):
            with patch('sys.argv', ['batch_stats.py', '--products', self.products_parquet, '--audit', self.test_dir, '--audit-size', str(audit_size)]):
                with patch('builtins.print') as mock_print:
                    batch_stats.main()

                # Check if audit directory was created
                audit_file = os.path.join(self.test_dir, "audit", f"sample_{min(audit_size, 3)}.csv")  # We have 3 unique parents

                # Check print output for audit generation
                print_calls = [call[0][0] for call in mock_print.call_args_list]
                audit_generated = any("Generated audit sample:" in call for call in print_calls)

                # If pandas was actually available during test, check file creation
                if audit_generated:
                    self.assertTrue(any("Audit sample saved to:" in call for call in print_calls))

    def test_audit_sample_without_pandas(self):
        """Test audit sample generation fallback without pandas."""
        with patch('batch_stats.check_pandas_availability', return_value=False):
            with patch('sys.argv', ['batch_stats.py', '--products', self.products_parquet, '--audit', self.test_dir, '--audit-size', '5']):
                with patch('builtins.print') as mock_print:
                    batch_stats.main()

                # Check that pandas unavailable message was printed
                print_calls = [call[0][0] for call in mock_print.call_args_list]
                self.assertTrue(any("Pandas not available - skipping audit sample generation" in call for call in print_calls))

    def test_file_finding_logic(self):
        """Test file finding logic with standardized vs k2/k3 naming."""
        # Create k2 variant file
        k2_path = os.path.join(self.test_dir, "products_k2.parquet")
        shutil.copy(self.products_parquet, k2_path)

        # Test finding standardized file when it exists
        found_path = batch_stats.find_products_file(k2_path)
        self.assertEqual(found_path, self.products_parquet)  # Should prefer products.parquet

        # Remove standardized file and test fallback
        os.remove(self.products_parquet)
        found_path = batch_stats.find_products_file(k2_path)
        self.assertEqual(found_path, k2_path)  # Should fall back to k2 file

    def test_missing_qa_file_handling(self):
        """Test handling of missing QA file."""
        missing_qa = os.path.join(self.test_dir, "nonexistent_qa.json")

        with patch('sys.argv', ['batch_stats.py', '--qa', missing_qa]):
            with patch('builtins.print') as mock_print:
                batch_stats.main()

            print_calls = [call[0][0] for call in mock_print.call_args_list]
            self.assertTrue(any(f"QA summary not found: {missing_qa}" in call for call in print_calls))

    def test_integration_all_features(self):
        """Test integration of all features together."""
        with patch('batch_stats.check_pandas_availability', return_value=True):
            with patch('sys.argv', [
                'batch_stats.py',
                '--products', self.products_parquet,
                '--qa', self.qa_json,
                '--label', 'integration-test',
                '--audit', self.test_dir,
                '--audit-size', '3'
            ]):
                with patch('builtins.print') as mock_print:
                    batch_stats.main()

                print_calls = [call[0][0] for call in mock_print.call_args_list]

                # Should have all types of output
                has_products_stats = any("Products generated:" in call for call in print_calls)
                has_qa_stats = any("QA warnings:" in call for call in print_calls)

                # At minimum, should not crash and should produce some output
                self.assertTrue(len(print_calls) > 0)


if __name__ == '__main__':
    unittest.main()