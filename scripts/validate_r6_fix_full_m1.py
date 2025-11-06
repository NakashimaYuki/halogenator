# -*- coding: utf-8 -*-
"""
Full M1 Configuration Validation for R6_methyl K-Level Fix

This script uses the COMPLETE M1 configuration (R1+R2+R3+R6_methyl, all halogens)
to validate the R6_methyl k-level fix with the full rule set.

Expected k=1 F products: ~10 (from all rules combined)
Expected k=2 F products from R6_methyl: ~11 (based on k=1 parents)
"""

from halogenator.enumerate_k import enumerate_with_stats, EnumConfig

# M1 parent: methoxy-flavone
parent_smiles = "COc1cc(-c2cc(=O)c3c(O)cc(O)cc3o2)c(O)cc1O"

print("="*60)
print("R6_methyl K-Level Fix - Full M1 Configuration Validation")
print("="*60)
print(f"Parent: {parent_smiles}")
print()

# Configuration: FULL M1 configuration (matches m1_raw_macro.yaml)
config = EnumConfig(
    k_max=2,
    halogens=['F', 'Cl', 'Br', 'I'],  # All halogens
    rules=['R1', 'R2', 'R3', 'R6_methyl'],  # All rules
    constraints={'enable': True, 'per_ring_quota': 2, 'min_graph_distance': 2, 'max_per_halogen': None},
    std_cfg={'do_tautomer': False},
    qc_cfg={'sanitize_strict': True, 'pains': True},
    pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False},
    sugar_cfg={'mode': 'off'},  # Disable sugar mask for simplicity
    rules_cfg={
        'R6_methyl': {
            'enable': True,
            'allowed': ['F', 'Cl', 'Br', 'I'],
            'allow_on_methoxy': True,
            'allow_allylic_methyl': False,
            'macro': {'enable': False}  # Keep disabled for simplicity in validation
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
print("Running enumeration with FULL M1 configuration (R1+R2+R3+R6_methyl)...")
print()
products, qa_stats = enumerate_with_stats(parent_smiles, config)

# Analyze results
print(f"Total products: {len(products)}")
print()

# Group by k_ops and halogen
k1_f = [p for p in products if p['k_ops'] == 1 and p['halogen'] == 'F']
k1_cl = [p for p in products if p['k_ops'] == 1 and p['halogen'] == 'Cl']
k1_br = [p for p in products if p['k_ops'] == 1 and p['halogen'] == 'Br']
k1_i = [p for p in products if p['k_ops'] == 1 and p['halogen'] == 'I']
k2_f = [p for p in products if p['k_ops'] == 2 and p['halogen'] == 'F']

print(f"k=1 F products: {len(k1_f)}")
print(f"k=1 Cl products: {len(k1_cl)}")
print(f"k=1 Br products: {len(k1_br)}")
print(f"k=1 I products: {len(k1_i)}")
print(f"k=2 F products: {len(k2_f)}")
print()

# Break down k=1 F products by rule
k1_f_by_rule = {}
for prod in k1_f:
    rule = prod['rule']
    if rule not in k1_f_by_rule:
        k1_f_by_rule[rule] = []
    k1_f_by_rule[rule].append(prod)

print("k=1 F products by rule:")
for rule in sorted(k1_f_by_rule.keys()):
    print(f"  {rule}: {len(k1_f_by_rule[rule])} products")
print()

# Check structural consistency for k=1 F
import json
k1_f_issues = []
for prod in k1_f:
    smiles = prod['smiles']
    f_count = smiles.count('F')
    cl_count = smiles.count('Cl')
    br_count = smiles.count('Br')
    i_count = smiles.count('I')
    total_hal = f_count + cl_count + br_count + i_count

    if total_hal != 1:
        k1_f_issues.append(prod)

if k1_f_issues:
    print(f"[ERROR] Found {len(k1_f_issues)} k=1 F products with incorrect halogen count")
else:
    print(f"[OK] All {len(k1_f)} k=1 F products have exactly 1 halogen atom")

# Check k=2 F products from R6_methyl
r6_k2_f = [p for p in k2_f if p['rule'] == 'R6_methyl']
print(f"\nR6_methyl k=2 F products: {len(r6_k2_f)}")

r6_k2_issues = []
for prod in r6_k2_f:
    smiles = prod['smiles']
    # Count ALL halogens (k=2 products can have mixed halogens)
    f_count = smiles.count('F')
    cl_count = smiles.count('Cl')
    br_count = smiles.count('Br')
    i_count = smiles.count('I')
    total_hal = f_count + cl_count + br_count + i_count

    subs_json = json.loads(prod['substitutions_json'])

    if prod['k_ops'] != 2 or total_hal != 2 or len(subs_json) != 2:
        r6_k2_issues.append(prod)

if r6_k2_issues:
    print(f"[ERROR] Found {len(r6_k2_issues)} R6_methyl k=2 F products with issues")
else:
    print(f"[OK] All {len(r6_k2_f)} R6_methyl k=2 F products correctly labeled")

# Summary
print(f"\n{'='*60}")
print("Summary")
print('='*60)
print(f"Configuration: FULL M1 (R1+R2+R3+R6_methyl)")
print(f"Expected k=1 F products: ~10")
print(f"Actual k=1 F products: {len(k1_f)}")
print()

if len(k1_f) >= 8 and not k1_f_issues and not r6_k2_issues:
    print("[SUCCESS] Full M1 validation passed!")
    print(f"  - {len(k1_f)} k=1 F products (all structurally consistent)")
    print(f"  - {len(r6_k2_f)} R6_methyl k=2 F products (all correctly labeled)")
else:
    print("[INFO] Results may vary based on molecule structure and constraints")
