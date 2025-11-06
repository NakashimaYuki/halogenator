#!/usr/bin/env python
"""
Advanced K-Level Mixing Diagnostic Tool

This script performs deep analysis of hierarchical SDF outputs to detect:
1. Structural k-level (actual number of substitutions) vs declared k property
2. Parent molecule contamination in product lists
3. Wrong parent molecules used as SDF headers
4. Detailed substitution analysis per molecule

Usage:
    python scripts/diagnose_k_level_mixing.py data/output/m1_raw_macro/mol_1/k1/F/mol_1_F.sdf
"""

import sys
import argparse
from typing import List, Dict, Tuple
from collections import Counter


def analyze_sdf_detailed(sdf_path: str, expected_k: int = 1) -> Dict:
    """
    Perform detailed analysis of an SDF file.

    Returns:
        Dict with analysis results including:
        - mol_count: total molecules
        - parent_count: molecules marked as parent
        - k_values: counter of k property values
        - halogen_counts: counter of total halogen count per molecule
        - substitutions_json_analysis: analysis of substitution history
        - violations: list of molecules with k != expected_k
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except ImportError:
        print("[ERROR] RDKit not installed")
        sys.exit(2)

    results = {
        'mol_count': 0,
        'parent_count': 0,
        'product_count': 0,
        'k_values': Counter(),
        'halogen_counts': Counter(),
        'substitutions_lengths': Counter(),
        'violations': [],
        'parent_molecules': [],
        'product_molecules': []
    }

    suppl = Chem.SDMolSupplier(sdf_path, sanitize=False, removeHs=False)

    for idx, mol in enumerate(suppl):
        if mol is None:
            print(f"  [WARN] Molecule {idx} failed to parse")
            continue

        results['mol_count'] += 1

        # Check if parent
        is_parent = mol.HasProp('parent') and mol.GetProp('parent').upper() == 'TRUE'

        # Get k property
        k_str = mol.GetProp('k') if mol.HasProp('k') else None
        k_int = int(k_str) if k_str else None

        # Get substitutions_json
        subs_json = mol.GetProp('substitutions_json') if mol.HasProp('substitutions_json') else None
        subs_length = 0
        if subs_json:
            try:
                import json
                subs_list = json.loads(subs_json)
                subs_length = len(subs_list)
            except:
                pass

        # Count halogens in structure
        halogen_count = 0
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol in ['F', 'Cl', 'Br', 'I']:
                halogen_count += 1

        # Get inchikey
        inchikey = mol.GetProp('inchikey') if mol.HasProp('inchikey') else 'UNKNOWN'

        # Get halogen property
        halogen_prop = mol.GetProp('halogen') if mol.HasProp('halogen') else None

        # Get parent_inchikey
        parent_inchikey = mol.GetProp('parent_inchikey') if mol.HasProp('parent_inchikey') else None

        # Get rule
        rule = mol.GetProp('rule') if mol.HasProp('rule') else None

        mol_info = {
            'index': idx,
            'is_parent': is_parent,
            'k_property': k_int,
            'halogen_count': halogen_count,
            'subs_json_length': subs_length,
            'inchikey': inchikey,
            'halogen_prop': halogen_prop,
            'parent_inchikey': parent_inchikey,
            'rule': rule,
            'smiles': Chem.MolToSmiles(mol)
        }

        if is_parent:
            results['parent_count'] += 1
            results['parent_molecules'].append(mol_info)
        else:
            results['product_count'] += 1
            results['product_molecules'].append(mol_info)

            # Update counters
            if k_int is not None:
                results['k_values'][k_int] += 1
            else:
                results['k_values']['missing'] += 1

            results['halogen_counts'][halogen_count] += 1
            results['substitutions_lengths'][subs_length] += 1

            # Check for violations (k != expected_k)
            if k_int != expected_k:
                results['violations'].append(mol_info)

    return results


def print_analysis_report(sdf_path: str, results: Dict, expected_k: int):
    """Print formatted analysis report."""

    print(f"\n{'='*80}")
    print(f"SDF ANALYSIS REPORT: {sdf_path}")
    print(f"Expected k-level: {expected_k}")
    print(f"{'='*80}\n")

    # Basic stats
    print(f"Total molecules: {results['mol_count']}")
    print(f"  - Parent molecules: {results['parent_count']}")
    print(f"  - Product molecules: {results['product_count']}")

    # K property analysis
    print(f"\nK Property Distribution (products only):")
    for k_val, count in sorted(results['k_values'].items()):
        marker = "[OK]" if k_val == expected_k else "[BAD]"
        print(f"  {marker} k={k_val}: {count} molecules")

    # Halogen count analysis
    print(f"\nStructural Halogen Count Distribution (products only):")
    for hal_count, count in sorted(results['halogen_counts'].items()):
        # For k=1, expect 1 halogen; for k=2, expect 2 halogens (usually)
        marker = "[OK]" if hal_count == expected_k else "[WARN]"
        print(f"  {marker} {hal_count} halogens: {count} molecules")

    # Substitutions JSON length analysis
    print(f"\nSubstitutions JSON Length Distribution (products only):")
    for subs_len, count in sorted(results['substitutions_lengths'].items()):
        marker = "[OK]" if subs_len == expected_k else "[WARN]"
        print(f"  {marker} {subs_len} substitutions: {count} molecules")

    # Violations
    if results['violations']:
        print(f"\n{'!'*80}")
        print(f"VIOLATIONS DETECTED: {len(results['violations'])} molecules with k != {expected_k}")
        print(f"{'!'*80}\n")

        # Group violations by k value
        violations_by_k = Counter([v['k_property'] for v in results['violations']])
        for k_val, count in sorted(violations_by_k.items()):
            print(f"  - k={k_val}: {count} molecules (should be k={expected_k})")

        # Show first 5 violations in detail
        print(f"\nFirst 5 violations (detailed):")
        for i, v in enumerate(results['violations'][:5]):
            print(f"\n  [{i+1}] Molecule #{v['index']}:")
            print(f"      InChIKey: {v['inchikey']}")
            print(f"      k property: {v['k_property']} (expected: {expected_k})")
            print(f"      Structural halogens: {v['halogen_count']}")
            print(f"      Substitutions JSON length: {v['subs_json_length']}")
            print(f"      Halogen: {v['halogen_prop']}")
            print(f"      Rule: {v['rule']}")
            print(f"      Parent InChIKey: {v['parent_inchikey']}")
            print(f"      SMILES: {v['smiles'][:80]}...")

    else:
        print(f"\n[OK] NO VIOLATIONS: All product molecules have k={expected_k}")

    # Parent molecule analysis
    if results['parent_molecules']:
        print(f"\nParent Molecule(s):")
        for i, p in enumerate(results['parent_molecules']):
            print(f"\n  [{i+1}] Molecule #{p['index']}:")
            print(f"      InChIKey: {p['inchikey']}")
            print(f"      Structural halogens: {p['halogen_count']}")
            print(f"      SMILES: {p['smiles'][:80]}...")

            # Check if parent has halogens (suspicious)
            if p['halogen_count'] > 0:
                print(f"      [WARN] WARNING: Parent has {p['halogen_count']} halogens!")
                print(f"             This may indicate wrong parent molecule used")

    print(f"\n{'='*80}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Deep analysis of hierarchical SDF for k-level mixing detection",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze k=1 F SDF
  python scripts/diagnose_k_level_mixing.py data/output/m1_raw_macro/mol_1/k1/F/mol_1_F.sdf

  # Analyze k=1 F SDF with expected k=1
  python scripts/diagnose_k_level_mixing.py data/output/m1_raw_macro/mol_1/k1/F/mol_1_F.sdf --expected-k 1

  # Analyze k=2 SDF
  python scripts/diagnose_k_level_mixing.py data/output/m1_raw_macro/mol_1/k2/F/AJPQ.../AJPQ..._F.sdf --expected-k 2
        """
    )

    parser.add_argument("sdf_path", help="Path to SDF file to analyze")
    parser.add_argument("--expected-k", type=int, default=1,
                       help="Expected k-level for this SDF (default: 1)")
    parser.add_argument("--export-violations",
                       help="Export violation details to JSON file")

    args = parser.parse_args()

    # Validate SDF file exists
    import os
    if not os.path.exists(args.sdf_path):
        print(f"[ERROR] File not found: {args.sdf_path}", file=sys.stderr)
        sys.exit(2)

    # Run analysis
    results = analyze_sdf_detailed(args.sdf_path, expected_k=args.expected_k)

    # Print report
    print_analysis_report(args.sdf_path, results, args.expected_k)

    # Export violations if requested
    if args.export_violations and results['violations']:
        import json
        with open(args.export_violations, 'w') as f:
            json.dump(results['violations'], f, indent=2)
        print(f"Violations exported to: {args.export_violations}")

    # Exit with appropriate code
    if results['violations']:
        sys.exit(1)
    else:
        sys.exit(0)


if __name__ == "__main__":
    main()
