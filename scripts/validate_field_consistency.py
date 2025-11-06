# -*- coding: utf-8 -*-
"""
Field Consistency Validator for Halogenator Parquet Output

Validates the following invariants:
1. k == k_ops (operation count consistency)
2. len(substitutions_json) == k_ops (history length consistency)
3. k_atoms == Σ atom_cost from substitutions_json (atom cost consistency)
4. Macro steps (CF3/CCl3) have atom_cost == 3

Usage:
    python scripts/validate_field_consistency.py <parquet_file_or_directory>

    # Single file
    python scripts/validate_field_consistency.py output/scenario/products.parquet

    # Directory (finds all .parquet files)
    python scripts/validate_field_consistency.py output/scenario/
"""

import sys
import json
from pathlib import Path
import pandas as pd


def validate_field_consistency(parquet_path: Path, verbose: bool = False):
    """
    Validate field consistency for a single parquet file.

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
        'k_mismatch': [],
        'history_length_mismatch': [],
        'k_atoms_mismatch': [],
        'macro_atom_cost_mismatch': [],
        'empty_history': [],
        'parse_errors': []
    }

    # Check each record
    for idx, row in df.iterrows():
        record_id = f"idx={idx}, inchikey={row.get('inchikey', 'UNKNOWN')[:14]}"

        # Extract fields
        k = row.get('k')
        k_ops = row.get('k_ops')
        k_atoms = row.get('k_atoms')
        substitutions_json_str = row.get('substitutions_json', '[]')

        # Parse substitutions_json
        try:
            substitutions = json.loads(substitutions_json_str) if substitutions_json_str else []
        except json.JSONDecodeError as e:
            results['parse_errors'].append({
                'record': record_id,
                'error': str(e),
                'raw': substitutions_json_str[:100]
            })
            continue

        # Validation 1: k == k_ops
        if k is not None and k_ops is not None:
            if k != k_ops:
                results['k_mismatch'].append({
                    'record': record_id,
                    'k': k,
                    'k_ops': k_ops,
                    'rule': row.get('rule')
                })

        # Validation 2: len(substitutions_json) == k_ops
        if k_ops is not None:
            if len(substitutions) != k_ops:
                results['history_length_mismatch'].append({
                    'record': record_id,
                    'k_ops': k_ops,
                    'history_length': len(substitutions),
                    'rule': row.get('rule')
                })

        # Validation 3: k_atoms == Σ atom_cost
        if k_atoms is not None and len(substitutions) > 0:
            # Calculate expected k_atoms from history
            expected_k_atoms = 0
            for step in substitutions:
                if 'atom_cost' in step:
                    expected_k_atoms += int(step['atom_cost'])
                elif step.get('type') == 'macro' and step.get('label') in ('CF3', 'CCl3'):
                    # Legacy: infer atom_cost=3 for macro steps
                    expected_k_atoms += 3
                else:
                    # Default: single atom substitution
                    expected_k_atoms += 1

            if k_atoms != expected_k_atoms:
                results['k_atoms_mismatch'].append({
                    'record': record_id,
                    'k_atoms': k_atoms,
                    'expected_k_atoms': expected_k_atoms,
                    'history': substitutions,
                    'rule': row.get('rule')
                })

        # Validation 4: Macro steps have atom_cost == 3
        for step_idx, step in enumerate(substitutions):
            if step.get('type') == 'macro' and step.get('label') in ('CF3', 'CCl3'):
                atom_cost = step.get('atom_cost')
                if atom_cost is not None and atom_cost != 3:
                    results['macro_atom_cost_mismatch'].append({
                        'record': record_id,
                        'step_index': step_idx,
                        'macro_label': step.get('label'),
                        'atom_cost': atom_cost,
                        'expected': 3
                    })

        # Validation 5: Check for empty history when k_ops > 0
        if k_ops is not None and k_ops > 0 and len(substitutions) == 0:
            results['empty_history'].append({
                'record': record_id,
                'k_ops': k_ops,
                'rule': row.get('rule')
            })

    # Print results
    print(f"\nValidation Results:")
    print(f"  Total records: {results['total_records']}")

    violations_found = False

    if results['k_mismatch']:
        violations_found = True
        print(f"\n  [FAIL] k != k_ops: {len(results['k_mismatch'])} violations")
        if verbose:
            for v in results['k_mismatch'][:5]:  # Show first 5
                print(f"    - {v['record']}: k={v['k']}, k_ops={v['k_ops']}, rule={v['rule']}")
            if len(results['k_mismatch']) > 5:
                print(f"    ... and {len(results['k_mismatch']) - 5} more")
    else:
        print(f"  [OK] k == k_ops: All {results['total_records']} records consistent")

    if results['history_length_mismatch']:
        violations_found = True
        print(f"\n  [FAIL] len(substitutions_json) != k_ops: {len(results['history_length_mismatch'])} violations")
        if verbose:
            for v in results['history_length_mismatch'][:5]:
                print(f"    - {v['record']}: k_ops={v['k_ops']}, history_len={v['history_length']}, rule={v['rule']}")
            if len(results['history_length_mismatch']) > 5:
                print(f"    ... and {len(results['history_length_mismatch']) - 5} more")
    else:
        print(f"  [OK] len(substitutions_json) == k_ops: All records consistent")

    if results['k_atoms_mismatch']:
        violations_found = True
        print(f"\n  [FAIL] k_atoms != Σ atom_cost: {len(results['k_atoms_mismatch'])} violations")
        if verbose:
            for v in results['k_atoms_mismatch'][:3]:
                print(f"    - {v['record']}: k_atoms={v['k_atoms']}, expected={v['expected_k_atoms']}, rule={v['rule']}")
            if len(results['k_atoms_mismatch']) > 3:
                print(f"    ... and {len(results['k_atoms_mismatch']) - 3} more")
    else:
        print(f"  [OK] k_atoms == Σ atom_cost: All records consistent")

    if results['macro_atom_cost_mismatch']:
        violations_found = True
        print(f"\n  [FAIL] Macro atom_cost != 3: {len(results['macro_atom_cost_mismatch'])} violations")
        if verbose:
            for v in results['macro_atom_cost_mismatch'][:5]:
                print(f"    - {v['record']}: macro={v['macro_label']}, atom_cost={v['atom_cost']} (expected 3)")
    else:
        print(f"  [OK] Macro atom_cost == 3: All macro steps consistent")

    if results['empty_history']:
        violations_found = True
        print(f"\n  [FAIL] Empty history with k_ops > 0: {len(results['empty_history'])} violations")
        if verbose:
            for v in results['empty_history'][:5]:
                print(f"    - {v['record']}: k_ops={v['k_ops']}, rule={v['rule']}")
    else:
        print(f"  [OK] No empty history for k_ops > 0")

    if results['parse_errors']:
        violations_found = True
        print(f"\n  [FAIL] JSON parse errors: {len(results['parse_errors'])} violations")
        if verbose:
            for v in results['parse_errors'][:3]:
                print(f"    - {v['record']}: {v['error']}")

    if not violations_found:
        print(f"\n[SUCCESS] All field consistency validations passed!")
    else:
        print(f"\n[FAILED] Found {sum(len(v) if isinstance(v, list) else 0 for v in results.values())} total violations")

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
        results = validate_field_consistency(parquet_file, verbose=verbose)
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

        file_violations = sum(len(v) if isinstance(v, list) else 0 for k, v in results.items() if k != 'total_records')
        total_violations += file_violations

        if file_violations == 0:
            print(f"  {parquet_file.name}: [OK] All validations passed ({results['total_records']} records)")
        else:
            print(f"  {parquet_file.name}: [FAIL] {file_violations} violations ({results['total_records']} records)")

    if total_violations == 0:
        print(f"\n[SUCCESS] All files passed field consistency validation!")
        sys.exit(0)
    else:
        print(f"\n[FAILED] Found {total_violations} total violations across all files")
        sys.exit(1)


if __name__ == '__main__':
    main()
