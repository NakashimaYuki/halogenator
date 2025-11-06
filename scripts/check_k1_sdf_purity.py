#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
K-level Purity Validation Script

This script scans hierarchical SDF output directories to detect k-level mixing issues:
- k=1 directories should ONLY contain molecules with k=1
- k=2 directories should ONLY contain molecules with k=2

Any violation indicates a critical bug in hierarchical output generation.

Usage:
    python scripts/check_k1_sdf_purity.py data/output/m1_strict
    python scripts/check_k1_sdf_purity.py data/output/g1_strict --check-k2

Exit codes:
    0: All k-level purity checks passed
    1: One or more k-level mixing violations detected
    2: Script error (file not found, invalid SDF, etc.)
"""

import os
import sys
import argparse
from pathlib import Path
from typing import List, Tuple, Dict


def scan_k_level_purity(root: str, check_k1: bool = True, check_k2: bool = False, check_structural: bool = False) -> List[Tuple[str, int, int, List[str]]]:
    """
    Scan directory tree for k-level purity violations.

    Args:
        root: Root directory to scan (e.g., data/output/m1_strict)
        check_k1: Whether to check k=1 directories (default: True)
        check_k2: Whether to check k=2 directories (default: False)
        check_structural: Whether to also check structural halogen count vs k (default: False)

    Returns:
        List of offenders: (filepath, total_mols, bad_count, k_values_found)
        Empty list if all checks pass.
    """
    # Import RDKit here to avoid import errors if not installed
    try:
        from rdkit import Chem
    except ImportError:
        print("[ERROR] RDKit not installed. Cannot validate SDF files.")
        sys.exit(2)

    offenders = []

    for dirpath, _, files in os.walk(root):
        # Normalize path separators for cross-platform compatibility
        parts = dirpath.replace("\\", "/").split("/")
        path_str = "/" + "/".join(parts) + "/"

        # Determine expected k level for this directory
        expected_k = None
        if check_k1 and "/k1/" in path_str:
            expected_k = 1
        elif check_k2 and "/k2/" in path_str:
            expected_k = 2
        else:
            continue  # Not a k-level directory we're checking

        # Process all SDF files in this directory
        for filename in files:
            if not filename.endswith(".sdf"):
                continue

            filepath = os.path.join(dirpath, filename)

            try:
                suppl = Chem.SDMolSupplier(filepath, sanitize=False, removeHs=False)
                bad_k_values = []
                total_mols = 0

                for mol in suppl:
                    if mol is None:
                        continue

                    total_mols += 1

                    # Skip parent molecule (first entry in SDF with parent=TRUE property)
                    # Parents don't have k values since they're the starting point
                    is_parent = mol.HasProp('parent') and mol.GetProp('parent').upper() == 'TRUE'
                    if is_parent:
                        continue

                    # Get k property from molecule
                    k_str = mol.GetProp('k') if mol.HasProp('k') else None

                    # Convert to int for comparison
                    try:
                        k_int = int(k_str) if k_str is not None else None
                    except (ValueError, TypeError):
                        k_int = None
                        bad_k_values.append(f"invalid:{k_str}")

                    # Check if k matches expected
                    if k_int is None:
                        bad_k_values.append("missing")
                    elif k_int != expected_k:
                        bad_k_values.append(str(k_int))

                    # Optional: Check structural halogen count vs k
                    if check_structural and k_int is not None and k_int > 0:
                        # Count halogen atoms in the molecule
                        halogen_count = sum(1 for atom in mol.GetAtoms()
                                          if atom.GetSymbol() in ('F', 'Cl', 'Br', 'I'))

                        # For non-macro scenarios, halogen_count should equal k
                        # Note: This is a simplified check; macro substitutions (CF3/CCl3)
                        # would have halogen_count != k_int
                        if halogen_count != k_int:
                            # Check if this might be a macro substitution
                            is_likely_macro = (halogen_count % 3 == 0 and halogen_count > k_int)
                            if not is_likely_macro:
                                bad_k_values.append(f"struct_mismatch:k={k_int},hal={halogen_count}")

                # If any violations found, record this file
                if bad_k_values:
                    unique_bad = sorted(set(bad_k_values))
                    offenders.append((filepath, total_mols, len(bad_k_values), unique_bad))

            except Exception as e:
                # Record files that can't be read
                offenders.append((filepath, 0, -1, [f"error: {e}"]))

    return offenders


def format_report(offenders: List[Tuple[str, int, int, List[str]]], check_k1: bool, check_k2: bool) -> str:
    """
    Format offenders into human-readable report.

    Args:
        offenders: List of (filepath, total, bad_count, k_values)
        check_k1: Whether k=1 was checked
        check_k2: Whether k=2 was checked

    Returns:
        Formatted report string
    """
    if not offenders:
        levels_checked = []
        if check_k1:
            levels_checked.append("k=1")
        if check_k2:
            levels_checked.append("k=2")
        return f"[OK] All {' and '.join(levels_checked)} purity checks PASSED."

    report_lines = [
        f"[FAIL] Found {len(offenders)} file(s) with k-level purity violations:",
        ""
    ]

    for filepath, total, bad_count, k_values in offenders:
        if bad_count == -1:
            # Error reading file
            report_lines.append(f"  - {filepath}")
            report_lines.append(f"    ERROR: {k_values[0]}")
        else:
            report_lines.append(f"  - {filepath}")
            report_lines.append(f"    Total molecules: {total}")
            report_lines.append(f"    Violations: {bad_count}")
            report_lines.append(f"    K values found: {k_values}")

        report_lines.append("")

    return "\n".join(report_lines)


def main():
    parser = argparse.ArgumentParser(
        description="Check k-level purity in hierarchical SDF outputs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Check k=1 purity only (default)
  python scripts/check_k1_sdf_purity.py data/output/m1_strict

  # Check both k=1 and k=2 purity
  python scripts/check_k1_sdf_purity.py data/output/m1_strict --check-k2

  # Check k=2 purity only
  python scripts/check_k1_sdf_purity.py data/output/g1_strict --no-check-k1 --check-k2

  # Batch check all scenarios
  for dir in data/output/*/; do
    python scripts/check_k1_sdf_purity.py "$dir" || echo "FAILED: $dir"
  done
        """
    )

    parser.add_argument("root", help="Root directory to scan (e.g., data/output/m1_strict)")
    parser.add_argument("--check-k2", action="store_true",
                       help="Also check k=2 directories for purity (default: False)")
    parser.add_argument("--no-check-k1", action="store_true",
                       help="Skip checking k=1 directories (default: check k=1)")
    parser.add_argument("--check-structural", action="store_true",
                       help="Also validate structural halogen count vs k value (default: False)")
    parser.add_argument("--verbose", action="store_true",
                       help="Print detailed progress information")

    args = parser.parse_args()

    # Validate root directory exists
    if not os.path.exists(args.root):
        print(f"[ERROR] Directory does not exist: {args.root}", file=sys.stderr)
        sys.exit(2)

    if not os.path.isdir(args.root):
        print(f"[ERROR] Not a directory: {args.root}", file=sys.stderr)
        sys.exit(2)

    # Determine what to check
    check_k1 = not args.no_check_k1
    check_k2 = args.check_k2

    if not check_k1 and not check_k2:
        print("[ERROR] Must check at least one k-level (k=1 or k=2)", file=sys.stderr)
        sys.exit(2)

    # Run scan
    if args.verbose:
        levels = []
        if check_k1:
            levels.append("k=1")
        if check_k2:
            levels.append("k=2")
        print(f"[INFO] Scanning {args.root} for {' and '.join(levels)} purity...", file=sys.stderr)
        if args.check_structural:
            print(f"[INFO] Also checking structural halogen count consistency...", file=sys.stderr)

    offenders = scan_k_level_purity(args.root, check_k1=check_k1, check_k2=check_k2, check_structural=args.check_structural)

    # Print report
    report = format_report(offenders, check_k1, check_k2)
    print(report)

    # Exit with appropriate code
    if offenders:
        sys.exit(1)
    else:
        sys.exit(0)


if __name__ == "__main__":
    main()
