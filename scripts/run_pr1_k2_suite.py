#!/usr/bin/env python
# -*- coding: ascii -*-
"""
PR1 k=2 Verification Suite Runner

This script runs the 8-molecule k=2 test suite to verify that:
1. Core enumeration functionality works without regressions
2. No new sanitization errors are introduced
3. Sugar masking works as expected on diverse molecular types
4. Performance remains reasonable for k=2 enumeration
"""

import json
import logging
import subprocess
import sys
import time
import yaml
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from halogenator.enumerate_k import enumerate_with_stats, EnumConfig

LOG = logging.getLogger(__name__)


def get_git_info() -> Dict[str, str]:
    """Get git repository information."""
    try:
        # Get git commit SHA
        sha_result = subprocess.run(['git', 'rev-parse', 'HEAD'],
                                  capture_output=True, text=True, cwd=Path(__file__).parent.parent)
        git_sha = sha_result.stdout.strip() if sha_result.returncode == 0 else 'unknown'

        # Get short SHA
        short_sha = git_sha[:8] if git_sha != 'unknown' else 'unknown'

        # Get branch name
        branch_result = subprocess.run(['git', 'rev-parse', '--abbrev-ref', 'HEAD'],
                                     capture_output=True, text=True, cwd=Path(__file__).parent.parent)
        branch = branch_result.stdout.strip() if branch_result.returncode == 0 else 'unknown'

        # Check if there are uncommitted changes
        status_result = subprocess.run(['git', 'status', '--porcelain'],
                                     capture_output=True, text=True, cwd=Path(__file__).parent.parent)
        has_changes = bool(status_result.stdout.strip()) if status_result.returncode == 0 else True

        return {
            'git_sha': git_sha,
            'git_short_sha': short_sha,
            'git_branch': branch,
            'has_uncommitted_changes': has_changes,
            'git_available': True
        }
    except Exception as e:
        LOG.warning(f"Failed to get git info: {e}")
        return {
            'git_sha': 'unknown',
            'git_short_sha': 'unknown',
            'git_branch': 'unknown',
            'has_uncommitted_changes': True,
            'git_available': False
        }


def load_test_suite(suite_file: str) -> Dict[str, Any]:
    """Load test suite configuration from YAML file."""
    with open(suite_file, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def run_molecule_test(molecule: Dict[str, Any], sugar_mode: str, config_base: Dict[str, Any]) -> Dict[str, Any]:
    """
    Run enumeration test for a single molecule.

    Returns:
        Test results dictionary with statistics and timing
    """
    mol_id = molecule['id']
    smiles = molecule['smiles']

    LOG.info(f"Testing {mol_id} with sugar_mode={sugar_mode}")

    start_time = time.time()

    try:
        # Create enumeration config
        config = EnumConfig(
            k_max=config_base['k_max'],
            halogens=tuple(config_base['halogens']),
            rules=tuple(config_base['rules']),
            sugar_cfg={'mode': sugar_mode, 'audit': True}
        )

        # Run enumeration
        products, qa_stats = enumerate_with_stats(smiles, config)

        end_time = time.time()
        runtime = end_time - start_time

        # Extract key metrics
        result = {
            'molecule_id': mol_id,
            'sugar_mode': sugar_mode,
            'success': True,
            'runtime_seconds': runtime,
            'total_products': len(products),
            'qa_stats': qa_stats,
            'key_metrics': {
                'attempts': qa_stats.get('attempts', 0),
                'products': qa_stats.get('products', 0),
                'template_unsupported': qa_stats.get('template_unsupported', 0),
                'sanitize_errors': qa_stats.get('qa_paths', {}).get('rdkit_error', 0),
                'sugar_filtered': qa_stats.get('qa_paths', {}).get('sugar_filtered_reaction_matches', 0),
                'sugar_blocked': qa_stats.get('qa_paths', {}).get('sugar_post_guard_blocked', 0)
            }
        }

        LOG.info(f"  {mol_id} ({sugar_mode}): {result['total_products']} products, "
                f"{result['runtime_seconds']:.2f}s, "
                f"{result['key_metrics']['sanitize_errors']} sanitize errors")

        return result

    except Exception as e:
        end_time = time.time()
        LOG.error(f"Test failed for {mol_id} with {sugar_mode}: {e}")

        return {
            'molecule_id': mol_id,
            'sugar_mode': sugar_mode,
            'success': False,
            'runtime_seconds': end_time - start_time,
            'error': str(e),
            'total_products': 0,
            'qa_stats': {},
            'key_metrics': {}
        }


def analyze_results(results: List[Dict[str, Any]], suite_config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Analyze test results and generate summary.

    Returns:
        Analysis summary with pass/fail status and recommendations
    """
    total_tests = len(results)
    successful_tests = sum(1 for r in results if r['success'])

    # Group results by molecule and mode
    by_molecule = {}
    by_mode = {'off': [], 'heuristic': []}

    total_sanitize_errors = 0
    total_products = 0
    total_runtime = 0

    for result in results:
        mol_id = result['molecule_id']
        mode = result['sugar_mode']

        if mol_id not in by_molecule:
            by_molecule[mol_id] = {}
        by_molecule[mol_id][mode] = result

        if result['success']:
            by_mode[mode].append(result)
            total_sanitize_errors += result['key_metrics'].get('sanitize_errors', 0)
            total_products += result['total_products']
            total_runtime += result['runtime_seconds']

    # Calculate error rates
    total_attempts = sum(r['key_metrics'].get('attempts', 0) for r in results if r['success'])
    sanitize_error_rate = total_sanitize_errors / total_attempts if total_attempts > 0 else 0

    # Analyze sugar masking effectiveness
    sugar_masking_analysis = {}
    for mol_id, modes in by_molecule.items():
        if 'off' in modes and 'heuristic' in modes and both_successful(modes):
            off_products = modes['off']['total_products']
            heuristic_products = modes['heuristic']['total_products']

            if off_products > 0:
                reduction_pct = (off_products - heuristic_products) / off_products * 100
            else:
                reduction_pct = 0

            sugar_masking_analysis[mol_id] = {
                'off_products': off_products,
                'heuristic_products': heuristic_products,
                'reduction_percentage': reduction_pct,
                'sugar_events': modes['heuristic']['key_metrics'].get('sugar_filtered', 0) +
                               modes['heuristic']['key_metrics'].get('sugar_blocked', 0)
            }

    # Determine overall status
    expected = suite_config.get('test_config', {}).get('expected_behavior', {})
    max_error_rate = expected.get('max_sanitize_error_rate', 0.05)

    issues = []
    if successful_tests < total_tests:
        issues.append(f"{total_tests - successful_tests} test(s) failed")
    if sanitize_error_rate > max_error_rate:
        issues.append(f"Sanitize error rate {sanitize_error_rate:.3f} exceeds threshold {max_error_rate}")

    overall_status = "PASS" if len(issues) == 0 else "FAIL"

    return {
        'overall_status': overall_status,
        'issues': issues,
        'summary': {
            'total_tests': total_tests,
            'successful_tests': successful_tests,
            'total_products': total_products,
            'total_runtime_seconds': total_runtime,
            'total_sanitize_errors': total_sanitize_errors,
            'sanitize_error_rate': sanitize_error_rate,
            'average_products_per_molecule': total_products / len(by_molecule) if by_molecule else 0
        },
        'by_molecule': by_molecule,
        'sugar_masking_analysis': sugar_masking_analysis
    }


def both_successful(modes: Dict[str, Dict[str, Any]]) -> bool:
    """Check if both off and heuristic modes were successful."""
    return (modes.get('off', {}).get('success', False) and
            modes.get('heuristic', {}).get('success', False))


def main():
    """Main test runner."""
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    # Paths
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    suite_file = project_root / "benchmarks/pr1_k2_suite.yaml"
    output_file = project_root / "benchmarks/pr1_k2_summary.json"

    LOG.info("Starting PR1 k=2 verification suite")

    # Load test suite
    suite_config = load_test_suite(suite_file)
    molecules = suite_config['molecules']
    test_config = suite_config['test_config']

    LOG.info(f"Testing {len(molecules)} molecules with k_max={test_config['k_max']}")

    # Run tests for each molecule in both sugar modes
    results = []

    for molecule in molecules:
        for sugar_mode in test_config['sugar_modes']:
            result = run_molecule_test(molecule, sugar_mode, test_config)
            results.append(result)

    # Analyze results
    analysis = analyze_results(results, suite_config)

    # Get environment and git information
    git_info = get_git_info()
    current_time = datetime.now()

    # Generate final report with comprehensive metadata
    report = {
        'test_metadata': {
            'suite_version': suite_config.get('version', '1.0'),
            'timestamp': current_time.isoformat(),
            'timestamp_unix': time.time(),
            'molecules_tested': len(molecules),
            'sugar_modes': test_config['sugar_modes'],
            'config_snapshot': test_config,
            'git_info': git_info,
            'environment': {
                'python_version': sys.version,
                'script_path': str(Path(__file__).resolve()),
                'working_directory': str(Path.cwd()),
                'suite_file': str(suite_file)
            }
        },
        'results': results,
        'analysis': analysis
    }

    # Write results
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    LOG.info(f"Test results written to {output_file}")

    # Print summary
    print(f"\nPR1 k=2 Verification Suite Results")
    print("=" * 40)
    print(f"Overall Status: {analysis['overall_status']}")
    print(f"Tests: {analysis['summary']['successful_tests']}/{analysis['summary']['total_tests']} passed")
    print(f"Total Products: {analysis['summary']['total_products']}")
    print(f"Runtime: {analysis['summary']['total_runtime_seconds']:.1f}s")
    print(f"Sanitize Errors: {analysis['summary']['total_sanitize_errors']} "
          f"(rate: {analysis['summary']['sanitize_error_rate']:.3f})")

    if analysis['issues']:
        print(f"\nIssues Found:")
        for issue in analysis['issues']:
            print(f"  - {issue}")
    else:
        print(f"\nAll tests passed! No regressions detected.")

    # Print sugar masking analysis
    print(f"\nSugar Masking Analysis:")
    for mol_id, data in analysis['sugar_masking_analysis'].items():
        print(f"  {mol_id}: {data['off_products']} -> {data['heuristic_products']} products "
              f"({data['reduction_percentage']:.1f}% reduction, {data['sugar_events']} sugar events)")

    # Exit with appropriate code
    sys.exit(0 if analysis['overall_status'] == 'PASS' else 1)


if __name__ == '__main__':
    main()