# -*- coding: ascii -*-
"""
Import style guard tests.

Prevents regression to incorrect import patterns that cause module loading issues.
"""

import pathlib
import re
import unittest


class TestImportStyleGuard(unittest.TestCase):
    """Guard against incorrect import patterns in test files."""

    def test_no_plain_package_imports_in_critical_tests(self):
        """Test that critical test files use 'from src.halogenator...' instead of 'from halogenator...'."""
        # Focus on the critical test files we've been working with for T19-T26
        critical_test_files = [
            'test_streaming_vs_snapshot_pivot_consistency.py',
            'test_schema_qa_utils_consolidation.py',
            'test_dual_track_counting_semantics.py',
            'test_alias_audit.py',
            'test_cli_qa_totals_merge.py',
            'test_import_style_guard.py',  # This file itself
            'test_schema_qa_utils_guard.py'
        ]

        tests_dir = pathlib.Path('tests')
        violations = []

        for file_name in critical_test_files:
            py_file = tests_dir / file_name
            if not py_file.exists():
                continue

            try:
                text = py_file.read_text(encoding='utf-8', errors='ignore')
                # Look for problematic import patterns (actual import statements only)
                if re.search(r'^\s*from\s+halogenator\b', text, re.MULTILINE):
                    violations.append(str(py_file))
            except Exception:
                continue

        if violations:
            self.fail(f"Found critical test files using 'from halogenator...' instead of 'from src.halogenator...': {violations}")

    def test_no_relative_imports_crossing_src_boundary(self):
        """Test that test files don't use relative imports to access src modules."""
        tests_dir = pathlib.Path('tests')
        violations = []

        for py_file in tests_dir.rglob('*.py'):
            if py_file.name == '__init__.py':
                continue

            try:
                text = py_file.read_text(encoding='utf-8', errors='ignore')
                # Look for relative imports that might cross boundaries
                if re.search(r'from\s+\.\..*halogenator', text):
                    violations.append(str(py_file))
            except Exception:
                continue

        if violations:
            self.fail(f"Found test files using relative imports to access src: {violations}")


if __name__ == '__main__':
    unittest.main()