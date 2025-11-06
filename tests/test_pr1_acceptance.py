# -*- coding: ascii -*-
"""
PR1 Sugar Masking Acceptance Tests

This module contains automated tests that validate PR1 sugar masking functionality
by running enumeration on test samples and asserting expected behavior.

Test Cases:
1. Glycoside samples show significant product reduction with sugar masking
2. Non-sugar controls show minimal difference between modes
3. QA metrics properly track sugar filtering events
4. No regressions in core enumeration functionality
"""

import json
import logging
import subprocess
import sys
import unittest
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

LOG = logging.getLogger(__name__)


class TestPR1SugarMaskingAcceptance(unittest.TestCase):
    """Acceptance tests for PR1 sugar masking functionality."""

    @classmethod
    def setUpClass(cls):
        """Run the acceptance test script and load results."""
        cls.project_root = Path(__file__).parent.parent
        cls.script_path = cls.project_root / "scripts/run_pr1_acceptance.py"
        cls.results_path = cls.project_root / "benchmarks/pr1_sugar_acceptance.json"

        # Run the acceptance test script
        LOG.info("Running PR1 acceptance test script...")
        try:
            result = subprocess.run([
                sys.executable, str(cls.script_path)
            ], capture_output=True, text=True, timeout=300)  # 5 minute timeout

            if result.returncode != 0:
                raise Exception(f"Acceptance script failed with return code {result.returncode}\\n"
                              f"STDOUT: {result.stdout}\\n"
                              f"STDERR: {result.stderr}")

            LOG.info("Acceptance test script completed successfully")

        except subprocess.TimeoutExpired:
            raise Exception("Acceptance test script timed out after 5 minutes")
        except Exception as e:
            raise Exception(f"Failed to run acceptance test script: {e}")

        # Load results
        if not cls.results_path.exists():
            raise Exception(f"Results file not found: {cls.results_path}")

        with open(cls.results_path, 'r', encoding='utf-8') as f:
            cls.results = json.load(f)

        LOG.info(f"Loaded acceptance test results with {len(cls.results['samples'])} samples")

    def test_glycoside_samples_show_significant_reduction(self):
        """Test that glycoside samples show significant product reduction with sugar masking."""
        glycoside_results = []

        for sample_id, data in self.results['samples'].items():
            sample_info = data['sample_info']
            if sample_info['type'] == 'glycoside':
                analysis = data['analysis']
                glycoside_results.append((sample_id, analysis))

                # Assert significant reduction (?30%) for each glycoside
                self.assertGreaterEqual(
                    analysis['reduction_percentage'], 30.0,
                    f"Glycoside {sample_id} should show ?30% product reduction, "
                    f"got {analysis['reduction_percentage']:.1f}% "
                    f"({analysis['off_products']} ? {analysis['heuristic_products']} products)"
                )

                # Assert some sugar filtering events occurred
                sugar_events = analysis['sugar_events']
                total_sugar_events = (sugar_events['filtered_matches'] +
                                    sugar_events['post_guard_blocked'])
                self.assertGreater(
                    total_sugar_events, 0,
                    f"Glycoside {sample_id} should have sugar filtering events, got 0"
                )

        # Assert we tested at least 3 glycosides
        self.assertGreaterEqual(
            len(glycoside_results), 3,
            f"Should test at least 3 glycoside samples, got {len(glycoside_results)}"
        )

        LOG.info(f"PASS: All {len(glycoside_results)} glycoside samples showed significant reduction")

    def test_control_samples_show_minimal_difference(self):
        """Test that non-sugar control samples show minimal difference between modes."""
        control_results = []

        for sample_id, data in self.results['samples'].items():
            sample_info = data['sample_info']
            if sample_info['type'] == 'aglycone':
                analysis = data['analysis']
                control_results.append((sample_id, analysis))

                # Assert minimal reduction (<30%) for controls
                self.assertLess(
                    analysis['reduction_percentage'], 30.0,
                    f"Control {sample_id} should show <30% product reduction, "
                    f"got {analysis['reduction_percentage']:.1f}% "
                    f"({analysis['off_products']} ? {analysis['heuristic_products']} products)"
                )

                # Assert minimal sugar filtering events
                sugar_events = analysis['sugar_events']
                total_sugar_events = (sugar_events['filtered_matches'] +
                                    sugar_events['post_guard_blocked'])
                self.assertLessEqual(
                    total_sugar_events, 2,  # Allow minimal false positives
                    f"Control {sample_id} should have minimal sugar filtering events, "
                    f"got {total_sugar_events}"
                )

        # Assert we tested at least 1 control
        self.assertGreaterEqual(
            len(control_results), 1,
            f"Should test at least 1 control sample, got {len(control_results)}"
        )

        LOG.info(f"PASS: All {len(control_results)} control samples showed minimal difference")

    def test_qa_metrics_tracking(self):
        """Test that QA metrics properly track sugar filtering events."""
        for sample_id, data in self.results['samples'].items():
            heuristic_stats = data['heuristic_mode']

            # Skip if enumeration failed
            if 'error' in heuristic_stats:
                continue

            qa_stats = heuristic_stats['qa_stats']
            qa_paths = qa_stats.get('qa_paths', {})

            # Verify QA paths structure exists
            self.assertIsInstance(qa_paths, dict, f"qa_paths should be dict for {sample_id}")

            # Verify sugar-specific QA keys are tracked (may be 0 for controls)
            expected_keys = ['sugar_filtered_reaction_matches', 'sugar_post_guard_blocked', 'sugar_mask_degraded']
            for key in expected_keys:
                self.assertIn(key, qa_paths, f"QA key '{key}' should be tracked for {sample_id}")
                self.assertIsInstance(qa_paths[key], int, f"QA key '{key}' should be integer for {sample_id}")

        LOG.info("? QA metrics structure is correct for all samples")

    def test_no_enumeration_regressions(self):
        """Test that core enumeration functionality works without regressions."""
        for sample_id, data in self.results['samples'].items():
            for mode in ['off_mode', 'heuristic_mode']:
                stats = data[mode]

                # Skip if enumeration failed
                if 'error' in stats:
                    self.fail(f"Enumeration failed for {sample_id} in {mode}: {stats['error']}")

                # Verify basic structure
                self.assertIn('total_products', stats, f"Missing total_products for {sample_id} {mode}")
                self.assertIn('qa_stats', stats, f"Missing qa_stats for {sample_id} {mode}")

                # Verify no sanitization failures (this would indicate serious issues)
                qa_stats = stats['qa_stats']
                template_unsupported = qa_stats.get('template_unsupported', 0)
                attempts = qa_stats.get('attempts', 0)

                if attempts > 0:
                    failure_rate = template_unsupported / attempts
                    self.assertLess(
                        failure_rate, 0.5,  # Allow up to 50% template failures
                        f"High template failure rate ({failure_rate:.1%}) for {sample_id} {mode}"
                    )

        LOG.info("? No enumeration regressions detected")

    def test_overall_assessment_criteria(self):
        """Test overall assessment meets acceptance criteria."""
        assessment = self.results['overall_assessment']

        # Verify glycoside performance
        glycoside_success_rate = (assessment['glycoside_significant_reductions'] /
                                assessment['glycoside_samples'])
        self.assertGreaterEqual(
            glycoside_success_rate, 0.8,  # 80% of glycosides should show reduction
            f"Glycoside success rate should be ?80%, got {glycoside_success_rate:.1%} "
            f"({assessment['glycoside_significant_reductions']}/{assessment['glycoside_samples']})"
        )

        # Verify control performance
        control_success_rate = (1 - assessment['control_significant_reductions'] /
                              assessment['control_samples'])
        self.assertGreaterEqual(
            control_success_rate, 0.8,  # 80% of controls should NOT show reduction
            f"Control success rate should be ?80%, got {control_success_rate:.1%} "
            f"({assessment['control_samples'] - assessment['control_significant_reductions']}/"
            f"{assessment['control_samples']})"
        )

        LOG.info(f"PASS: Overall assessment: {glycoside_success_rate:.1%} glycoside success, "
                f"{control_success_rate:.1%} control success")


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    unittest.main(verbosity=2)