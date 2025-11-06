# -*- coding: utf-8 -*-
"""
CI Validation Gate - Unified Entry Point for All Halogenator Data Quality Checks

This script runs all validation checks on Halogenator output and serves as the
CI gate for data quality. All checks must pass for the gate to succeed.

Validation checks:
1. Field Consistency: k==k_ops, len(substitutions_json)==k_ops, k_atoms==Σatom_cost
2. Structural Consistency: structural_halogen_count matches expected count from history
3. Parent Chain Integrity: parent_inchikey references form valid k-1→k-2→...→k=0 chains
4. SDF Purity: k=1/k=2 directories contain only products with correct k values

Usage:
    # Standard mode (checks parquet files + hierarchical SDF if present)
    python scripts/ci_validation_gate.py <output_directory>

    # Paranoid mode (includes expensive structural validation)
    python scripts/ci_validation_gate.py <output_directory> --paranoid

    # Skip SDF checks (only validate parquet)
    python scripts/ci_validation_gate.py <output_directory> --skip-sdf

    # Verbose output
    python scripts/ci_validation_gate.py <output_directory> --verbose

Examples:
    # Check M1 scenario output
    python scripts/ci_validation_gate.py test_m1_output/

    # Check with full paranoid validation
    python scripts/ci_validation_gate.py test_m1_output/ --paranoid --verbose

Exit codes:
    0: All validations passed
    1: One or more validations failed
    2: Script error (missing directory, no parquet files, etc.)
"""

import sys
import subprocess
from pathlib import Path
import argparse


def run_validation(script_name: str, args: list, description: str, verbose: bool = False):
    """
    Run a validation script and return results.

    Args:
        script_name: Name of the validation script (in scripts/)
        args: Arguments to pass to the script
        description: Human-readable description of the check
        verbose: If True, show detailed output

    Returns:
        tuple: (passed: bool, output: str)
    """
    script_path = Path(__file__).parent / script_name
    if not script_path.exists():
        return False, f"[ERROR] Validation script not found: {script_path}"

    cmd = [sys.executable, str(script_path)] + args

    if verbose:
        print(f"\n{'='*80}")
        print(f"Running: {description}")
        print(f"Command: {' '.join(cmd)}")
        print('='*80)

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout per check
        )

        output = result.stdout + result.stderr

        if verbose:
            print(output)

        passed = (result.returncode == 0)

        if not verbose and not passed:
            # Show output only if failed and not in verbose mode
            print(f"\n[FAIL] {description}")
            print(output)

        return passed, output

    except subprocess.TimeoutExpired:
        return False, f"[ERROR] {description} timed out after 10 minutes"
    except Exception as e:
        return False, f"[ERROR] {description} failed with exception: {e}"


def find_parquet_files(directory: Path) -> list:
    """Find all parquet files in directory."""
    return list(directory.rglob('*.parquet'))


def find_hierarchical_sdf_roots(directory: Path) -> list:
    """
    Find hierarchical SDF output roots (directories containing k1/ and/or k2/ subdirs).

    Returns:
        list: Paths to parent directories containing k-level subdirectories
    """
    roots = set()

    for path in directory.rglob('k1'):
        if path.is_dir():
            roots.add(path.parent)

    for path in directory.rglob('k2'):
        if path.is_dir():
            roots.add(path.parent)

    return sorted(roots)


def main():
    parser = argparse.ArgumentParser(
        description="CI Validation Gate - Run all Halogenator data quality checks",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument("directory", help="Output directory to validate")
    parser.add_argument("--paranoid", action="store_true",
                       help="Enable paranoid mode (includes expensive structural validation)")
    parser.add_argument("--skip-sdf", action="store_true",
                       help="Skip SDF purity checks (only validate parquet)")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Show detailed output from all validation scripts")
    parser.add_argument("--skip-field", action="store_true",
                       help="Skip field consistency check (for debugging)")
    parser.add_argument("--skip-structural", action="store_true",
                       help="Skip structural consistency check (for debugging)")
    parser.add_argument("--skip-parent-chain", action="store_true",
                       help="Skip parent chain check (for debugging)")

    args = parser.parse_args()

    directory = Path(args.directory)

    if not directory.exists():
        print(f"[ERROR] Directory not found: {directory}")
        sys.exit(2)

    if not directory.is_dir():
        print(f"[ERROR] Not a directory: {directory}")
        sys.exit(2)

    print(f"{'='*80}")
    print(f"CI Validation Gate")
    print(f"{'='*80}")
    print(f"Directory: {directory}")
    print(f"Paranoid mode: {'ON' if args.paranoid else 'OFF'}")
    print()

    # Find parquet files
    parquet_files = find_parquet_files(directory)
    if not parquet_files:
        print(f"[ERROR] No parquet files found in {directory}")
        sys.exit(2)

    print(f"Found {len(parquet_files)} parquet file(s)")

    # Find hierarchical SDF roots (if any)
    sdf_roots = [] if args.skip_sdf else find_hierarchical_sdf_roots(directory)
    if sdf_roots:
        print(f"Found {len(sdf_roots)} hierarchical SDF output root(s)")
    elif not args.skip_sdf:
        print("[INFO] No hierarchical SDF outputs found (k1/k2 directories)")

    # Run validation checks
    checks = []
    all_passed = True

    # 1. Field Consistency Check (parquet)
    if not args.skip_field:
        for parquet_file in parquet_files:
            check_args = [str(parquet_file)]
            if args.verbose:
                check_args.append('--verbose')

            passed, output = run_validation(
                'validate_field_consistency.py',
                check_args,
                f"Field Consistency: {parquet_file.name}",
                verbose=args.verbose
            )
            checks.append(('Field Consistency', parquet_file.name, passed))
            if not passed:
                all_passed = False

    # 2. Structural Consistency Check (parquet) - optional paranoid mode
    if args.paranoid and not args.skip_structural:
        for parquet_file in parquet_files:
            check_args = [str(parquet_file)]
            if args.verbose:
                check_args.append('--verbose')

            passed, output = run_validation(
                'validate_structural_consistency.py',
                check_args,
                f"Structural Consistency: {parquet_file.name}",
                verbose=args.verbose
            )
            checks.append(('Structural Consistency', parquet_file.name, passed))
            if not passed:
                all_passed = False

    # 3. Parent Chain Integrity Check (parquet)
    if not args.skip_parent_chain:
        for parquet_file in parquet_files:
            check_args = [str(parquet_file)]
            if args.verbose:
                check_args.append('--verbose')

            passed, output = run_validation(
                'validate_parent_chain.py',
                check_args,
                f"Parent Chain Integrity: {parquet_file.name}",
                verbose=args.verbose
            )
            checks.append(('Parent Chain', parquet_file.name, passed))
            if not passed:
                all_passed = False

    # 4. SDF Purity Check (hierarchical SDF outputs)
    if sdf_roots and not args.skip_sdf:
        for sdf_root in sdf_roots:
            check_args = [str(sdf_root), '--check-k2']
            if args.paranoid:
                check_args.append('--check-structural')
            if args.verbose:
                check_args.append('--verbose')

            passed, output = run_validation(
                'check_k1_sdf_purity.py',
                check_args,
                f"SDF Purity: {sdf_root.name}",
                verbose=args.verbose
            )
            checks.append(('SDF Purity', sdf_root.name, passed))
            if not passed:
                all_passed = False

    # Summary Report
    print(f"\n{'='*80}")
    print("Validation Summary")
    print('='*80)

    # Group by check type
    from collections import defaultdict
    by_type = defaultdict(list)
    for check_type, target, passed in checks:
        by_type[check_type].append((target, passed))

    for check_type in ['Field Consistency', 'Structural Consistency', 'Parent Chain', 'SDF Purity']:
        if check_type not in by_type:
            continue

        results = by_type[check_type]
        passed_count = sum(1 for _, p in results if p)
        total_count = len(results)

        status = '[OK]' if passed_count == total_count else '[FAIL]'
        print(f"\n{status} {check_type}: {passed_count}/{total_count} passed")

        if args.verbose or passed_count < total_count:
            for target, passed in results:
                status_icon = '[OK]' if passed else '[FAIL]'
                print(f"  {status_icon} {target}")

    # Final result
    print(f"\n{'='*80}")
    if all_passed:
        print("[SUCCESS] All validation checks passed!")
        print('='*80)
        sys.exit(0)
    else:
        failed_count = sum(1 for _, _, p in checks if not p)
        print(f"[FAILED] {failed_count}/{len(checks)} checks failed")
        print('='*80)
        sys.exit(1)


if __name__ == '__main__':
    main()
