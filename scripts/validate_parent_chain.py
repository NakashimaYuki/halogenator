# -*- coding: utf-8 -*-
"""
Parent Chain Validator for Halogenator Parquet Output

Validates that parent_inchikey references form a valid chain:
1. k=2 products: parent_inchikey must exist in k=1 product set
2. k=1 products: parent_inchikey must be the root parent (k=0)
3. No orphan nodes (products with non-existent parents)

Usage:
    python scripts/validate_parent_chain.py <parquet_file_or_directory>

    # Single file
    python scripts/validate_parent_chain.py output/scenario/products.parquet

    # Directory (finds all .parquet files)
    python scripts/validate_parent_chain.py output/scenario/
"""

import sys
from pathlib import Path
from collections import defaultdict
import pandas as pd


def validate_parent_chain(parquet_path: Path, verbose: bool = False):
    """
    Validate parent chain integrity for a single parquet file.

    Args:
        parquet_path: Path to parquet file
        verbose: If True, show detailed output

    Returns:
        dict: Validation results with counts of violations
    """
    print(f"\n{'='*80}")
    print(f"Validating: {parquet_path}")
    print('='*80)

    try:
        df = pd.read_parquet(parquet_path)
    except Exception as e:
        print(f"[ERROR] Failed to read parquet file: {e}")
        return {'error': str(e)}

    print(f"Total records: {len(df)}")

    # Build sets of InChIKeys by k-level
    inchikeys_by_k = defaultdict(set)
    root_parent_inchikeys = set()

    for _, row in df.iterrows():
        k_ops = row.get('k_ops', 0)
        inchikey = row.get('inchikey')
        root_parent_inchikey = row.get('root_parent_inchikey')

        if inchikey:
            inchikeys_by_k[k_ops].add(inchikey)

        if root_parent_inchikey:
            root_parent_inchikeys.add(root_parent_inchikey)

    # Print k-level distribution
    print(f"\nK-level distribution:")
    for k in sorted(inchikeys_by_k.keys()):
        print(f"  k={k}: {len(inchikeys_by_k[k])} unique products")
    print(f"  Root parents: {len(root_parent_inchikeys)} unique")

    # Validation results
    results = {
        'total_records': len(df),
        'k_distribution': {k: len(v) for k, v in inchikeys_by_k.items()},
        'orphan_k1_products': [],
        'orphan_k2_products': [],
        'orphan_k_high_products': [],
        'missing_root_parent': []
    }

    # Validate each record
    for idx, row in df.iterrows():
        record_id = f"idx={idx}, inchikey={row.get('inchikey', 'UNKNOWN')[:14]}"

        k_ops = row.get('k_ops', 0)
        inchikey = row.get('inchikey')
        parent_inchikey = row.get('parent_inchikey')
        root_parent_inchikey = row.get('root_parent_inchikey')
        rule = row.get('rule', 'UNKNOWN')

        # Skip k=0 (root parents have no parent to validate)
        if k_ops == 0:
            continue

        # Validation 1: k=1 products should have parent_inchikey == root_parent_inchikey
        # (or at least root_parent_inchikey should be in root set)
        if k_ops == 1:
            # Check if parent_inchikey is a known root
            if parent_inchikey and parent_inchikey not in root_parent_inchikeys:
                # Could also be in k=0 set if k=0 products are included
                if parent_inchikey not in inchikeys_by_k[0]:
                    results['orphan_k1_products'].append({
                        'record': record_id,
                        'inchikey': inchikey,
                        'parent_inchikey': parent_inchikey,
                        'root_parent_inchikey': root_parent_inchikey,
                        'rule': rule
                    })

            # Check if root_parent_inchikey is provided
            if not root_parent_inchikey:
                results['missing_root_parent'].append({
                    'record': record_id,
                    'k_ops': k_ops,
                    'inchikey': inchikey
                })

        # Validation 2: k=2 products should have parent_inchikey in k=1 set
        elif k_ops == 2:
            if parent_inchikey and parent_inchikey not in inchikeys_by_k[1]:
                results['orphan_k2_products'].append({
                    'record': record_id,
                    'inchikey': inchikey,
                    'parent_inchikey': parent_inchikey,
                    'rule': rule
                })

            # Check if root_parent_inchikey is provided
            if not root_parent_inchikey:
                results['missing_root_parent'].append({
                    'record': record_id,
                    'k_ops': k_ops,
                    'inchikey': inchikey
                })

        # Validation 3: k>2 products (if any) should have parent_inchikey in k-1 set
        elif k_ops > 2:
            if parent_inchikey and parent_inchikey not in inchikeys_by_k[k_ops - 1]:
                results['orphan_k_high_products'].append({
                    'record': record_id,
                    'k_ops': k_ops,
                    'inchikey': inchikey,
                    'parent_inchikey': parent_inchikey,
                    'rule': rule
                })

    # Print results
    print(f"\nValidation Results:")

    violations_found = False

    if results['orphan_k1_products']:
        violations_found = True
        print(f"\n  [FAIL] k=1 orphan products (parent not in root set): {len(results['orphan_k1_products'])} violations")
        if verbose:
            for v in results['orphan_k1_products'][:5]:
                print(f"    - {v['record']}")
                print(f"      parent_inchikey: {v['parent_inchikey'][:14]}...")
                print(f"      root_parent_inchikey: {v['root_parent_inchikey'][:14] if v['root_parent_inchikey'] else 'None'}...")
                print(f"      rule: {v['rule']}")
            if len(results['orphan_k1_products']) > 5:
                print(f"    ... and {len(results['orphan_k1_products']) - 5} more")
    else:
        k1_count = len(inchikeys_by_k[1])
        if k1_count > 0:
            print(f"  [OK] k=1 parent chain valid: All {k1_count} products have valid root parent")

    if results['orphan_k2_products']:
        violations_found = True
        print(f"\n  [FAIL] k=2 orphan products (parent not in k=1 set): {len(results['orphan_k2_products'])} violations")
        if verbose:
            for v in results['orphan_k2_products'][:5]:
                print(f"    - {v['record']}")
                print(f"      parent_inchikey: {v['parent_inchikey'][:14]}... (not found in k=1 set)")
                print(f"      rule: {v['rule']}")
            if len(results['orphan_k2_products']) > 5:
                print(f"    ... and {len(results['orphan_k2_products']) - 5} more")
    else:
        k2_count = len(inchikeys_by_k[2])
        if k2_count > 0:
            print(f"  [OK] k=2 parent chain valid: All {k2_count} products have valid k=1 parent")

    if results['orphan_k_high_products']:
        violations_found = True
        print(f"\n  [FAIL] k>2 orphan products (parent not in k-1 set): {len(results['orphan_k_high_products'])} violations")
        if verbose:
            for v in results['orphan_k_high_products'][:5]:
                print(f"    - {v['record']}")
                print(f"      k_ops={v['k_ops']}, parent_inchikey: {v['parent_inchikey'][:14]}...")
    else:
        k_high_count = sum(len(v) for k, v in inchikeys_by_k.items() if k > 2)
        if k_high_count > 0:
            print(f"  [OK] k>2 parent chain valid: All {k_high_count} products have valid k-1 parent")

    if results['missing_root_parent']:
        violations_found = True
        print(f"\n  [FAIL] Missing root_parent_inchikey: {len(results['missing_root_parent'])} violations")
        if verbose:
            for v in results['missing_root_parent'][:5]:
                print(f"    - {v['record']}: k_ops={v['k_ops']}")
    else:
        checked_count = sum(len(v) for k, v in inchikeys_by_k.items() if k > 0)
        if checked_count > 0:
            print(f"  [OK] root_parent_inchikey present: All {checked_count} products have root parent")

    if not violations_found:
        print(f"\n[SUCCESS] All parent chain validations passed!")
    else:
        total_violations = (len(results['orphan_k1_products']) +
                          len(results['orphan_k2_products']) +
                          len(results['orphan_k_high_products']) +
                          len(results['missing_root_parent']))
        print(f"\n[FAILED] Found {total_violations} violations")

    return results


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    path = Path(sys.argv[1])
    verbose = '--verbose' in sys.argv or '-v' in sys.argv

    if not path.exists():
        print(f"[ERROR] Path not found: {path}")
        sys.exit(1)

    # Find all parquet files
    if path.is_file():
        parquet_files = [path]
    else:
        parquet_files = list(path.rglob('*.parquet'))
        if not parquet_files:
            print(f"[ERROR] No parquet files found in directory: {path}")
            sys.exit(1)

    print(f"Found {len(parquet_files)} parquet file(s) to validate")

    # Validate each file
    all_results = []
    for parquet_file in parquet_files:
        results = validate_parent_chain(parquet_file, verbose=verbose)
        all_results.append((parquet_file, results))

    # Summary
    print(f"\n{'='*80}")
    print("Summary")
    print('='*80)

    total_violations = 0
    for parquet_file, results in all_results:
        if 'error' in results:
            print(f"  {parquet_file.name}: ERROR - {results['error']}")
            continue

        file_violations = (len(results['orphan_k1_products']) +
                         len(results['orphan_k2_products']) +
                         len(results['orphan_k_high_products']) +
                         len(results['missing_root_parent']))
        total_violations += file_violations

        if file_violations == 0:
            k_dist = ", ".join(f"k={k}:{count}" for k, count in sorted(results['k_distribution'].items()))
            print(f"  {parquet_file.name}: [OK] All validations passed ({k_dist})")
        else:
            print(f"  {parquet_file.name}: [FAIL] {file_violations} violations")

    if total_violations == 0:
        print(f"\n[SUCCESS] All files passed parent chain validation!")
        sys.exit(0)
    else:
        print(f"\n[FAILED] Found {total_violations} total violations across all files")
        sys.exit(1)


if __name__ == '__main__':
    main()
