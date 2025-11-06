#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
Reproducibility Verification Script for P0-H Acceptance Tests

This script:
1. Runs acceptance tests twice in full mode (default)
2. Or compares existing full mode artifacts (--reuse-latest N or --from-files A B)
3. Compares results for identical numerical values
4. Generates a reproducibility verification report

Used to verify that random seed fixing is effective and results are deterministic.
"""

import argparse
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path


def run_acceptance_test():
    """Run the acceptance test in full mode and return the results file path."""
    script_dir = Path(__file__).parent
    acceptance_script = script_dir / "run_pr1_acceptance.py"

    # Run the acceptance test
    result = subprocess.run([
        sys.executable, str(acceptance_script), "--verbose"
    ], capture_output=True, text=True, cwd=script_dir.parent)

    if result.returncode != 0:
        print(f"Acceptance test failed with return code {result.returncode}")
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)
        return None

    # Find the most recent artifacts file
    artifacts_dir = script_dir.parent / "artifacts/acceptance"
    if not artifacts_dir.exists():
        print("No artifacts directory found")
        return None

    # Find the most recent *-full.json file
    json_files = list(artifacts_dir.glob("*-full.json"))
    if not json_files:
        print("No full mode results found in artifacts")
        return None

    # Return the most recent file
    latest_file = max(json_files, key=lambda f: f.stat().st_mtime)
    return latest_file


def extract_key_metrics(json_data):
    """Extract key numerical metrics for comparison."""
    test_meta = json_data['test_metadata']
    key_metrics = {
        'test_metadata': {
            'schema_version': test_meta.get('schema_version', 'unknown'),
            'random_seed': test_meta['random_seed'],
            'samples_file_hash': test_meta['samples_file_hash'],
            'halogens': test_meta['halogens'],
            'rules': test_meta['rules'],
            'k_max': test_meta['k_max'],
            'fast_mode': test_meta['fast_mode'],
            'sugar_ring_score_threshold': test_meta.get('sugar_masking_config', {}).get('sugar_ring_score_threshold'),
            'rdkit_seed_report': test_meta.get('rdkit_seed_report', {})
        },
        'overall_assessment': json_data['overall_assessment'],
        'sample_results': {}
    }

    # Extract per-sample numerical results
    for sample_id, sample_data in json_data['samples'].items():
        analysis = sample_data['analysis']
        heuristic_mode = sample_data.get('heuristic_mode', {})

        key_metrics['sample_results'][sample_id] = {
            'products_off': analysis['off_products'],
            'products_heur': analysis['heuristic_products'],
            'attempts_off': analysis['off_attempts'],
            'attempts_heur': analysis['heuristic_attempts'],
            'reduction_pct': analysis['reduction_percentage'],
            'overall_pass': analysis['overall_pass'],
            'mask_size': heuristic_mode.get('mask_size', 0),
            'degraded': heuristic_mode.get('degraded', False),
            'sugar_events': analysis.get('sugar_events', {}),
            'p0_h2_assertion_pass': sample_data['p0_h2_assertion']['assertion_pass']
        }

    return key_metrics


def _collect_numeric(path, value, result):
    if isinstance(value, dict):
        for key, subval in value.items():
            next_path = f"{path}.{key}" if path else key
            _collect_numeric(next_path, subval, result)
    elif isinstance(value, list):
        for idx, item in enumerate(value):
            next_path = f"{path}[{idx}]"
            _collect_numeric(next_path, item, result)
    elif isinstance(value, (int, float, bool)):
        result[path] = value


def build_numeric_index(metrics):
    result = {}
    overall = metrics.get('overall_assessment', {})
    _collect_numeric('overall_assessment', overall, result)

    samples = metrics.get('sample_results', {})
    for sample_id, sample_metrics in samples.items():
        base_path = f"sample_results.{sample_id}"
        _collect_numeric(base_path, sample_metrics, result)

    return result


def compare_results(run1_data, run2_data, *, strict=False):
    """Compare numeric intersections of two result sets."""
    index1 = build_numeric_index(run1_data)
    index2 = build_numeric_index(run2_data)

    keys1 = set(index1.keys())
    keys2 = set(index2.keys())
    shared_keys = sorted(keys1 & keys2)

    differences = []
    for key in shared_keys:
        if index1[key] != index2[key]:
            differences.append(f"{key}: {index1[key]} != {index2[key]}")

    ignored_keys = []
    if strict:
        for key in sorted(keys1 - keys2):
            differences.append(f"{key}: missing in run2")
        for key in sorted(keys2 - keys1):
            differences.append(f"{key}: missing in run1")
    else:
        ignored_keys = sorted(keys1 ^ keys2)

    return differences, ignored_keys


def get_latest_artifacts(n: int = 2):
    """Get the latest N full mode artifacts from the artifacts directory."""
    script_dir = Path(__file__).parent
    artifacts_dir = script_dir.parent / "artifacts/acceptance"

    if not artifacts_dir.exists():
        print("ERROR: No artifacts directory found")
        return None

    # Find all *-full.json files
    json_files = list(artifacts_dir.glob("*-full.json"))
    if len(json_files) < n:
        print(f"ERROR: Need {n} artifacts but only found {len(json_files)}")
        return None

    # Sort by modification time (newest first)
    json_files.sort(key=lambda f: f.stat().st_mtime, reverse=True)
    return json_files[:n]


def compare_and_report(run1_file: Path, run2_file: Path, *, strict=False):
    """Compare two artifact files and generate report."""
    print("\nLoading and comparing results...")

    with open(run1_file, 'r', encoding='utf-8') as f:
        run1_data = json.load(f)

    with open(run2_file, 'r', encoding='utf-8') as f:
        run2_data = json.load(f)

    # Extract key metrics
    run1_metrics = extract_key_metrics(run1_data)
    run2_metrics = extract_key_metrics(run2_data)

    # Compare results
    differences, ignored_keys = compare_results(run1_metrics, run2_metrics, strict=strict)

    # Generate report
    utc_timestamp = datetime.utcnow().strftime('%Y%m%d_%H%M%S')
    report_file = run1_file.parent / f"{utc_timestamp}-repro-diff.json"

    report = {
        'timestamp': utc_timestamp,
        'run1_file': run1_file.name,
        'run2_file': run2_file.name,
        'run1_schema_version': run1_metrics['test_metadata'].get('schema_version', 'unknown'),
        'run2_schema_version': run2_metrics['test_metadata'].get('schema_version', 'unknown'),
        'differences_count': len(differences),
        'differences': differences,
        'ignored_keys': ignored_keys,
        'reproducibility_status': 'PASS' if len(differences) == 0 else 'FAIL',
        'run1_summary': {
            'total_samples': run1_data['overall_assessment']['p0_h2_total_samples'],
            'p0_h2_passes': run1_data['overall_assessment']['p0_h2_total_passes'],
            'random_seed': run1_data['test_metadata']['random_seed']
        },
        'run2_summary': {
            'total_samples': run2_data['overall_assessment']['p0_h2_total_samples'],
            'p0_h2_passes': run2_data['overall_assessment']['p0_h2_total_passes'],
            'random_seed': run2_data['test_metadata']['random_seed']
        }
    }

    # Save report
    with open(report_file, 'w', encoding='utf-8') as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    # Print results
    print(f"\nReproducibility Verification Results:")
    print(f"  Run 1: {run1_file.name}")
    print(f"  Run 2: {run2_file.name}")
    print(f"  Differences found: {len(differences)}")
    print(f"  Status: {'PASS' if len(differences) == 0 else 'FAIL'}")

    if ignored_keys and not strict:
        print(f"\n[WARN] Ignored {len(ignored_keys)} key(s) due to schema differences:")
        for key in ignored_keys[:10]:
            print(f"    {key}")
        if len(ignored_keys) > 10:
            remaining = len(ignored_keys) - 10
            print(f"    ... and {remaining} more keys")

    if differences:
        print(f"\nDifferences detected:")
        for diff in differences[:10]:  # Show first 10 differences
            print(f"    {diff}")
        if len(differences) > 10:
            print(f"    ... and {len(differences) - 10} more differences")
        print(f"Reproducibility: FAIL")
    else:
        print(f"\n[PASS] Numeric intersections match between runs")
        print(f"  Random seed control is effective")
        print(f"Reproducibility: PASS")

    print(f"\nDetailed report saved to: {report_file}")

    # Exit with appropriate code
    sys.exit(0 if len(differences) == 0 else 1)


def main():
    """Main reproducibility verification function."""
    parser = argparse.ArgumentParser(
        description="Verify reproducibility of P0-H acceptance tests",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python verify_reproducibility.py                    # Run two new tests
  python verify_reproducibility.py --reuse-latest 2  # Compare latest 2 artifacts
  python verify_reproducibility.py --from-files A B  # Compare specific files
        """
    )

    parser.add_argument('--strict', action='store_true', help='require identical key sets instead of numeric intersections')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--reuse-latest', type=int, metavar='N',
                       help='Compare the latest N artifact files (default: 2)')
    group.add_argument('--from-files', nargs=2, metavar=('FILE1', 'FILE2'),
                       help='Compare two specific artifact files')

    args = parser.parse_args()

    print("Reproducibility Verification for P0-H Acceptance Tests")
    print("=" * 60)

    if args.from_files:
        # Compare specific files
        run1_file = Path(args.from_files[0])
        run2_file = Path(args.from_files[1])

        if not run1_file.exists():
            print(f"ERROR: File not found: {run1_file}")
            sys.exit(1)
        if not run2_file.exists():
            print(f"ERROR: File not found: {run2_file}")
            sys.exit(1)

        print(f"Comparing files:")
        print(f"  File 1: {run1_file}")
        print(f"  File 2: {run2_file}")

        compare_and_report(run1_file, run2_file, strict=args.strict)

    elif args.reuse_latest is not None:
        # Compare latest N artifacts
        n = args.reuse_latest or 2
        print(f"Comparing latest {n} artifacts...")

        artifacts = get_latest_artifacts(n)
        if not artifacts:
            sys.exit(1)

        run1_file, run2_file = artifacts[0], artifacts[1]
        print(f"  Latest:  {run1_file.name}")
        print(f"  Second:  {run2_file.name}")

        compare_and_report(run1_file, run2_file, strict=args.strict)

    else:
        # Default: run two new tests
        print("\nRunning first acceptance test...")
        run1_file = run_acceptance_test()
        if not run1_file:
            print("ERROR: First run failed")
            sys.exit(1)

        print(f"First run completed: {run1_file.name}")

        print("\nRunning second acceptance test...")
        run2_file = run_acceptance_test()
        if not run2_file:
            print("ERROR: Second run failed")
            sys.exit(1)

        print(f"Second run completed: {run2_file.name}")

        compare_and_report(run1_file, run2_file, strict=args.strict)


if __name__ == '__main__':
    main()