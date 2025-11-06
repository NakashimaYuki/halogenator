# -*- coding: utf-8 -*-
"""
Quick validation script for R6_methyl k-level fix.
Runs M1 parent molecule and checks k-level integrity.
"""

import json
from halogenator.enumerate_k import enumerate_with_stats, EnumConfig

# M1 parent: methoxy-flavone
parent_smiles = "COc1cc(-c2cc(=O)c3c(O)cc(O)cc3o2)c(O)cc1O"

print("="*60)
print("R6_methyl K-Level Fix Validation")
print("="*60)
print(f"Parent: {parent_smiles[:50]}...")
print()

# Configuration: M1 raw macro equivalent
config = EnumConfig(
    k_max=2,
    halogens=['F', 'Cl'],  # Only F and Cl for speed
    rules=['R1', 'R6_methyl'],  # Focus on R1 and R6
    constraints={'enable': True, 'per_ring_quota': 2, 'min_graph_distance': 2, 'max_per_halogen': None},
    std_cfg={'do_tautomer': False},
    qc_cfg={'sanitize_strict': True, 'pains': True},
    pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False},
    sugar_cfg={'mode': 'off'},  # Disable sugar mask for simplicity
    rules_cfg={
        'R6_methyl': {
            'enable': True,
            'allowed': ['F', 'Cl'],
            'allow_on_methoxy': True,
            'allow_allylic_methyl': False,
            'macro': {'enable': False}  # Disable macro for speed
        }
    },
    engine_cfg={
        'budget_mode': 'ops',
        'enable_inchi_dedup': True,
        'dedup_stage': 'pre',
        'dedup_policy': 'auto'
    },
    symmetry_cfg={
        'compute_on_masked_subgraph': True
    }
)

# Run enumeration
print("Running enumeration...")
products, qa_stats = enumerate_with_stats(parent_smiles, config)

# Analyze results
print(f"\nTotal products: {len(products)}")

# Group by k_ops and halogen
k1_f = [p for p in products if p['k_ops'] == 1 and p['halogen'] == 'F']
k1_cl = [p for p in products if p['k_ops'] == 1 and p['halogen'] == 'Cl']
k2_f = [p for p in products if p['k_ops'] == 2 and p['halogen'] == 'F']
k2_cl = [p for p in products if p['k_ops'] == 2 and p['halogen'] == 'Cl']

print(f"\nk=1 F products: {len(k1_f)}")
print(f"k=1 Cl products: {len(k1_cl)}")
print(f"k=2 F products: {len(k2_f)}")
print(f"k=2 Cl products: {len(k2_cl)}")

# Check k=1 F products for structural consistency
print(f"\n{'='*60}")
print("K=1 F Products - Structural Validation")
print('='*60)

k1_f_issues = []
for prod in k1_f:
    smiles = prod['smiles']
    f_count = smiles.count('F')
    cl_count = smiles.count('Cl')
    total_hal = f_count + cl_count

    subs_json = json.loads(prod['substitutions_json'])

    if total_hal != 1:
        k1_f_issues.append({
            'inchikey': prod['inchikey'][:14],
            'rule': prod['rule'],
            'k_ops': prod['k_ops'],
            'k_atoms': prod['k_atoms'],
            'structural_halogens': total_hal,
            'subs_len': len(subs_json),
            'smiles': smiles
        })

if k1_f_issues:
    print(f"\n[ERROR] Found {len(k1_f_issues)} k=1 F products with structural inconsistency:")
    for issue in k1_f_issues:
        print(f"  InChIKey: {issue['inchikey']}...")
        print(f"    k_ops={issue['k_ops']}, struct_hals={issue['structural_halogens']}, subs_len={issue['subs_len']}")
        print(f"    rule={issue['rule']}, SMILES={issue['smiles']}")
else:
    print(f"[OK] All {len(k1_f)} k=1 F products have exactly 1 halogen atom")

# Check k=2 F products from R6_methyl
print(f"\n{'='*60}")
print("K=2 F Products from R6_methyl - Validation")
print('='*60)

r6_k2_f = [p for p in k2_f if p['rule'] == 'R6_methyl']
print(f"Found {len(r6_k2_f)} k=2 F products from R6_methyl")

r6_k2_issues = []
for prod in r6_k2_f:
    smiles = prod['smiles']
    f_count = smiles.count('F')
    cl_count = smiles.count('Cl')
    total_hal = f_count + cl_count

    subs_json = json.loads(prod['substitutions_json'])

    issues = []
    if prod['k_ops'] != 2:
        issues.append(f"k_ops={prod['k_ops']} (expected 2)")
    if total_hal != 2:
        issues.append(f"struct_hals={total_hal} (expected 2)")
    if len(subs_json) != 2:
        issues.append(f"subs_len={len(subs_json)} (expected 2)")

    if issues:
        r6_k2_issues.append({
            'inchikey': prod['inchikey'][:14],
            'issues': issues,
            'smiles': smiles
        })

if r6_k2_issues:
    print(f"\n[ERROR] Found {len(r6_k2_issues)} R6_methyl k=2 F products with issues:")
    for issue in r6_k2_issues:
        print(f"  InChIKey: {issue['inchikey']}...")
        for err in issue['issues']:
            print(f"    - {err}")
else:
    print(f"[OK] All {len(r6_k2_f)} R6_methyl k=2 F products correctly labeled")

# Summary
print(f"\n{'='*60}")
print("Validation Summary")
print('='*60)

if not k1_f_issues and not r6_k2_issues:
    print("[SUCCESS] R6_methyl k-level fix validated!")
    print(f"  - {len(k1_f)} k=1 F products: all structurally consistent")
    print(f"  - {len(r6_k2_f)} R6 k=2 F products: all correctly labeled")
    print("\n[OK] Fix confirmed: No k=2 products mislabeled as k=1")
else:
    print("[FAILED] Validation found issues:")
    if k1_f_issues:
        print(f"  - {len(k1_f_issues)} k=1 F products with wrong halogen count")
    if r6_k2_issues:
        print(f"  - {len(r6_k2_issues)} R6 k=2 F products with labeling errors")
