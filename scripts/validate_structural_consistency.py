# -*- coding: utf-8 -*-
"""
Structural Consistency Validator for Halogenator Parquet Output

Validates that the actual number of halogen atoms in the product structure
matches the expected count based on k_ops and substitution history.

For non-macro scenarios:
    structural_halogen_count == k_ops

For macro scenarios (CF3/CCl3):
    structural_halogen_count == Σ(halogen_count_per_step)

Usage:
    python scripts/validate_structural_consistency.py <parquet_file_or_directory>

    # Single file
    python scripts/validate_structural_consistency.py output/scenario/products.parquet

    # Directory (finds all .parquet files)
    python scripts/validate_structural_consistency.py output/scenario/

    # Strict mode (fails on any mismatch)
    python scripts/validate_structural_consistency.py output/scenario/ --strict
"""

import sys
import json
from pathlib import Path
import pandas as pd


def count_halogens_in_smiles(smiles: str) -> dict:
    """
    Count halogen atoms in a SMILES string.

    Returns:
        dict: {'F': count, 'Cl': count, 'Br': count, 'I': count, 'total': count}
    """
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'F': 0, 'Cl': 0, 'Br': 0, 'I': 0, 'total': 0, 'error': 'Invalid SMILES'}

        counts = {'F': 0, 'Cl': 0, 'Br': 0, 'I': 0}
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol in counts:
                counts[symbol] += 1

        counts['total'] = sum(counts.values())
        return counts

    except Exception as e:
        return {'F': 0, 'Cl': 0, 'Br': 0, 'I': 0, 'total': 0, 'error': str(e)}


def calculate_expected_halogen_count(substitutions: list) -> int:
    """
    Calculate expected halogen count from substitution history.

    For each step:
    - If atom_cost is present, use it
    - If type='macro' and label in ('CF3', 'CCl3'), use 3
    - Otherwise, use 1 (single atom substitution)

    Returns:
        int: Expected total halogen count
    """
    total = 0
    for step in substitutions:
        if 'atom_cost' in step:
            total += int(step['atom_cost'])
        elif step.get('type') == 'macro' and step.get('label') in ('CF3', 'CCl3'):
            total += 3
        else:
            # Default: single atom substitution
            total += 1
    return total


def validate_structural_consistency(parquet_path: Path, strict: bool = False, verbose: bool = False):
    """
    Validate structural consistency for a single parquet file.

    Args:
        parquet_path: Path to parquet file
        strict: If True, fail on any mismatch (even for k=0 parent molecules)
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

    # Validation results
    results = {
        'total_records': len(df),
        'structural_mismatch': [],
        'smiles_parse_errors': [],
        'skipped_k0': 0  # Parent molecules (k=0) are typically not halogenated
    }

    # Check each record
    for idx, row in df.iterrows():
        record_id = f"idx={idx}, inchikey={row.get('inchikey', 'UNKNOWN')[:14]}"

        # Extract fields
        k_ops = row.get('k_ops', 0)
        smiles = row.get('smiles', '')
        substitutions_json_str = row.get('substitutions_json', '[]')
        rule = row.get('rule', 'UNKNOWN')

        # Skip k=0 parent molecules (they have no substitutions)
        if k_ops == 0:
            results['skipped_k0'] += 1
            continue

        # Parse substitutions_json
        try:
            substitutions = json.loads(substitutions_json_str) if substitutions_json_str else []
        except json.JSONDecodeError:
            # Already handled by field consistency validator
            continue

        # Count halogens in structure
        halogen_counts = count_halogens_in_smiles(smiles)

        if 'error' in halogen_counts:
            results['smiles_parse_errors'].append({
                'record': record_id,
                'error': halogen_counts['error'],
                'smiles': smiles[:100]
            })
            continue

        structural_halogen_count = halogen_counts['total']

        # Calculate expected halogen count from history
        expected_halogen_count = calculate_expected_halogen_count(substitutions)

        # Check consistency
        if structural_halogen_count != expected_halogen_count:
            results['structural_mismatch'].append({
                'record': record_id,
                'k_ops': k_ops,
                'rule': rule,
                'structural_count': structural_halogen_count,
                'expected_count': expected_halogen_count,
                'halogen_details': {k: v for k, v in halogen_counts.items() if k != 'total'},
                'substitutions': substitutions,
                'smiles': smiles
            })

    # Print results
    print(f"\nValidation Results:")
    print(f"  Total records: {results['total_records']}")
    print(f"  Skipped k=0 parent molecules: {results['skipped_k0']}")

    violations_found = False

    if results['structural_mismatch']:
        violations_found = True
        print(f"\n  [FAIL] Structural halogen count mismatch: {len(results['structural_mismatch'])} violations")
        if verbose:
            for v in results['structural_mismatch'][:5]:
                print(f"    - {v['record']}")
                print(f"      k_ops={v['k_ops']}, rule={v['rule']}")
                print(f"      Structural count: {v['structural_count']} {v['halogen_details']}")
                print(f"      Expected count: {v['expected_count']}")
                print(f"      SMILES: {v['smiles'][:80]}...")
            if len(results['structural_mismatch']) > 5:
                print(f"    ... and {len(results['structural_mismatch']) - 5} more")
    else:
        checked_count = results['total_records'] - results['skipped_k0']
        print(f"  [OK] Structural halogen count matches expected: All {checked_count} records consistent")

    if results['smiles_parse_errors']:
        violations_found = True
        print(f"\n  [FAIL] SMILES parse errors: {len(results['smiles_parse_errors'])} errors")
        if verbose:
            for v in results['smiles_parse_errors'][:3]:
                print(f"    - {v['record']}: {v['error']}")

    if not violations_found:
        print(f"\n[SUCCESS] All structural consistency validations passed!")
    else:
        print(f"\n[FAILED] Found {len(results['structural_mismatch']) + len(results['smiles_parse_errors'])} violations")

    return results


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    path = Path(sys.argv[1])
    strict = '--strict' in sys.argv
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
    if strict:
        print("[INFO] Running in STRICT mode (will fail on any mismatch)")

    # Validate each file
    all_results = []
    for parquet_file in parquet_files:
        results = validate_structural_consistency(parquet_file, strict=strict, verbose=verbose)
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

        file_violations = len(results['structural_mismatch']) + len(results['smiles_parse_errors'])
        total_violations += file_violations

        if file_violations == 0:
            checked_count = results['total_records'] - results['skipped_k0']
            print(f"  {parquet_file.name}: [OK] All validations passed ({checked_count} checked, {results['skipped_k0']} skipped k=0)")
        else:
            print(f"  {parquet_file.name}: [FAIL] {file_violations} violations ({results['total_records']} total records)")

    if total_violations == 0:
        print(f"\n[SUCCESS] All files passed structural consistency validation!")
        sys.exit(0)
    else:
        print(f"\n[FAILED] Found {total_violations} total violations across all files")
        sys.exit(1)


if __name__ == '__main__':
    main()
