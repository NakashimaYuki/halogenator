# -*- coding: ascii -*-
"""Tests for logging debouncing in seed management."""

import io
import logging
import unittest

from src.halogenator.rdkit_seed_utils import (
    set_rdkit_random_seed,
    _reset_seed_manager_state_for_tests,
)


class TestSeedLoggingDebounce(unittest.TestCase):
    """Ensure RDKit seed logging respects debouncing rules."""

    def setUp(self):
        _reset_seed_manager_state_for_tests()
        self.logger = logging.getLogger('src.halogenator.rdkit_seed_utils')
        self.logger.setLevel(logging.DEBUG)
        self.logger.propagate = False
        self.stream = io.StringIO()
        self.handler = logging.StreamHandler(self.stream)
        self.handler.setLevel(logging.DEBUG)
        self.logger.addHandler(self.handler)

    def tearDown(self):
        self.logger.removeHandler(self.handler)
        self.handler.close()
        self.logger.propagate = True
        _reset_seed_manager_state_for_tests()

    def test_info_summary_emit_once(self):
        """Repeated seed calls should log a single INFO summary and deduplicated DEBUG."""
        set_rdkit_random_seed(123)
        set_rdkit_random_seed(123)
        set_rdkit_random_seed(456)

        log_lines = [line for line in self.stream.getvalue().splitlines() if line.strip()]
        info_lines = [line for line in log_lines if 'RDKit seeding summary' in line]
        detail_lines = [line for line in log_lines if 'failed to seed' in line]
        unchanged_lines = [line for line in log_lines if 'details unchanged' in line]

        self.assertEqual(len(info_lines), 1, msg=f'Unexpected INFO lines: {info_lines}')
        self.assertEqual(len(detail_lines), 1, msg=f'Unexpected DEBUG detail lines: {detail_lines}')
        self.assertGreaterEqual(len(unchanged_lines), 1)


if __name__ == '__main__':
    unittest.main()
