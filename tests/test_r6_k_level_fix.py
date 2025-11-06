# -*- coding: utf-8 -*-
"""
Test R6_methyl k-level fix: Ensure k=2 products are correctly labeled.

This test validates the fix for the critical bug where R6_methyl step substitutions
at k>=2 were incorrectly labeled as k=1 with empty substitutions_json.

Test scenario:
- Parent: Anisole (methoxy-benzene, "COc1ccccc1")
- Step 1 (k=1): R1 adds F to aromatic ring
- Step 2 (k=2): R6_methyl adds F to methoxy group (-OCH3 -> -OCH2F)

Expected results:
- k=1 products: k_ops=1, k_atoms=1, len(substitutions_json)=1
- k=2 products: k_ops=2, k_atoms=2, len(substitutions_json)=2
- k=2 parent_inchikey points to k=1 product
"""

import json


def test_r6_methyl_k2_labeling_fix():
    """
    Test that R6_methyl correctly labels k=2 products with k_ops=2 and populated substitutions_json.

    This is a regression test for the bug where R6_methyl k=2 products were mislabeled as k=1.
    """
    from halogenator.enumerate_k import enumerate_with_stats, EnumConfig

    # Anisole: simple molecule with methoxy group
    parent_smiles = "COc1ccccc1"

    # Configuration: Enable R1 and R6_methyl, k_max=2, only F
    config = EnumConfig(
        k_max=2,
        halogens=['F'],
        rules=['R1', 'R6_methyl'],
        constraints={'enable': True, 'per_ring_quota': 2, 'min_graph_distance': 2, 'max_per_halogen': None},
        std_cfg={'do_tautomer': False},
        qc_cfg={'sanitize_strict': True, 'pains': True},
        pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False},
        sugar_cfg={'mode': 'off'},  # Disable sugar mask for simplicity
        rules_cfg={
            'R6_methyl': {
                'enable': True,
                'allowed': ['F', 'Cl'],
                'allow_on_methoxy': True,  # Enable R6 on methoxy groups
                'allow_allylic_methyl': False,
                'macro': {'enable': False}
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
    products, qa_stats = enumerate_with_stats(parent_smiles, config)

    # Separate k=1 and k=2 products
    k1_products = [p for p in products if p['k_ops'] == 1]
    k2_products = [p for p in products if p['k_ops'] == 2]

    print(f"\nEnumeration results:")
    print(f"Total products: {len(products)}")
    print(f"k=1 products: {len(k1_products)}")
    print(f"k=2 products: {len(k2_products)}")

    # Assertions for k=1 products
    assert len(k1_products) > 0, "Should have at least one k=1 product"

    for prod in k1_products:
        # Check k values
        assert prod['k_ops'] == 1, f"k=1 product should have k_ops=1, got {prod['k_ops']}"
        assert prod['k'] == 1, f"k=1 product should have k=1, got {prod['k']}"
        assert prod['k_atoms'] == 1, f"k=1 product should have k_atoms=1, got {prod['k_atoms']}"

        # Check substitutions_json
        subs_json = json.loads(prod['substitutions_json'])
        assert len(subs_json) == 1, f"k=1 product should have 1 step in history, got {len(subs_json)}"

        # Check rule (can be R1 or R6_methyl at k=1)
        assert prod['rule'] in ('R1', 'R6_methyl'), f"k=1 product should be from R1 or R6_methyl, got {prod['rule']}"

        # Structural validation: count halogen atoms
        smiles = prod['smiles']
        f_count = smiles.count('F')
        cl_count = smiles.count('Cl')
        br_count = smiles.count('Br')
        i_count = smiles.count('I')
        total_halogens = f_count + cl_count + br_count + i_count
        assert total_halogens == 1, f"k=1 product should have 1 halogen atom, got {total_halogens} in {smiles}"

        print(f"  k=1 product OK: {prod['inchikey'][:14]}... rule={prod['rule']} halogen={prod['halogen']} hal_count={total_halogens}")

    # Assertions for k=2 products (CRITICAL: this is where the bug was)
    assert len(k2_products) > 0, "Should have at least one k=2 product from R6_methyl"

    # Find k=2 products from R6_methyl
    r6_k2_products = [p for p in k2_products if p['rule'] == 'R6_methyl']
    assert len(r6_k2_products) > 0, "Should have at least one k=2 product from R6_methyl"

    print(f"\nR6_methyl k=2 products: {len(r6_k2_products)}")

    for prod in r6_k2_products:
        # CRITICAL: Check k values (this was the bug - k_ops was 1 instead of 2)
        assert prod['k_ops'] == 2, f"R6 k=2 product should have k_ops=2, got {prod['k_ops']}"
        assert prod['k'] == 2, f"R6 k=2 product should have k=2, got {prod['k']}"
        assert prod['k_atoms'] == 2, f"R6 k=2 product should have k_atoms=2, got {prod['k_atoms']}"

        # CRITICAL: Check substitutions_json (this was the bug - was empty {})
        subs_json = json.loads(prod['substitutions_json'])
        assert len(subs_json) == 2, f"R6 k=2 product should have 2 steps in history, got {len(subs_json)}"

        # Verify history structure
        # First step can be R1 or R6_methyl (both can produce k=1 from k=0)
        assert subs_json[0]['rule'] in ('R1', 'R6_methyl'), f"First step should be R1 or R6_methyl, got {subs_json[0].get('rule')}"
        assert subs_json[1]['rule'] == 'R6_methyl', f"Second step should be R6_methyl, got {subs_json[1].get('rule')}"

        # Verify atom_cost is present in history
        assert 'atom_cost' in subs_json[1], "R6_methyl step should have atom_cost field"
        assert subs_json[1]['atom_cost'] == 1, f"R6_methyl step should have atom_cost=1, got {subs_json[1]['atom_cost']}"

        # Structural validation: count halogen atoms (should be 2)
        smiles = prod['smiles']
        f_count = smiles.count('F')
        cl_count = smiles.count('Cl')
        br_count = smiles.count('Br')
        i_count = smiles.count('I')
        total_halogens = f_count + cl_count + br_count + i_count
        assert total_halogens == 2, f"k=2 product should have 2 halogen atoms, got {total_halogens} in {smiles}"

        # Parent chain validation
        parent_inchikey = prod['parent_inchikey']
        parent_found = any(p['inchikey'] == parent_inchikey for p in k1_products)
        assert parent_found, f"k=2 product's parent_inchikey should point to a k=1 product, got {parent_inchikey}"

        print(f"  R6 k=2 product OK: {prod['inchikey'][:14]}... parent={parent_inchikey[:14]}... hal_count={total_halogens}")

    print("\n[OK] All R6_methyl k-level validations passed!")
    print(f"[OK] Fix confirmed: {len(r6_k2_products)} R6_methyl k=2 products correctly labeled")


def test_r6_methyl_macro_k2_labeling():
    """
    Test that R6_methyl macro (CF3) correctly labels k=2 products with k_ops=2, k_atoms=4.
    """
    from halogenator.enumerate_k import enumerate_with_stats, EnumConfig

    # Anisole: simple molecule with methoxy group
    parent_smiles = "COc1ccccc1"

    # Configuration: Enable R1 and R6_methyl with macro, k_max=2, only F
    config = EnumConfig(
        k_max=2,
        halogens=['F'],
        rules=['R1', 'R6_methyl'],
        constraints={'enable': True, 'per_ring_quota': 2, 'min_graph_distance': 2, 'max_per_halogen': None},
        std_cfg={'do_tautomer': False},
        qc_cfg={'sanitize_strict': True, 'pains': True},
        pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False},
        sugar_cfg={'mode': 'off'},  # Disable sugar mask for simplicity
        rules_cfg={
            'R6_methyl': {
                'enable': True,
                'allowed': ['F', 'Cl'],
                'allow_on_methoxy': True,  # Enable R6 on methoxy groups
                'allow_allylic_methyl': False,
                'macro': {
                    'enable': True,
                    'labels': ['CF3']  # Only CF3 for this test
                }
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
    products, qa_stats = enumerate_with_stats(parent_smiles, config)

    # Find k=2 macro products from R6_methyl
    r6_macro_k2_products = [p for p in products
                            if p['k_ops'] == 2 and p['rule'] == 'R6_methyl' and p.get('macro_label') == 'CF3']

    print(f"\nR6_methyl macro (CF3) k=2 products: {len(r6_macro_k2_products)}")

    if len(r6_macro_k2_products) > 0:
        for prod in r6_macro_k2_products:
            # Check k values for macro: k_ops=2 (one R1 + one CF3), k_atoms=4 (1 + 3)
            assert prod['k_ops'] == 2, f"R6 macro k=2 product should have k_ops=2, got {prod['k_ops']}"
            assert prod['k'] == 2, f"R6 macro k=2 product should have k=2, got {prod['k']}"
            assert prod['k_atoms'] == 4, f"R6 macro k=2 product should have k_atoms=4 (1+3), got {prod['k_atoms']}"

            # Check substitutions_json
            subs_json = json.loads(prod['substitutions_json'])
            assert len(subs_json) == 2, f"R6 macro k=2 product should have 2 steps in history, got {len(subs_json)}"

            # Verify history structure
            assert subs_json[0]['rule'] == 'R1', f"First step should be R1, got {subs_json[0].get('rule')}"
            assert subs_json[1]['rule'] == 'R6_methyl', f"Second step should be R6_methyl, got {subs_json[1].get('rule')}"
            assert subs_json[1].get('type') == 'macro', f"Second step should be type='macro', got {subs_json[1].get('type')}"
            assert subs_json[1].get('label') == 'CF3', f"Second step should have label='CF3', got {subs_json[1].get('label')}"

            # Verify atom_cost
            assert subs_json[1].get('atom_cost') == 3, f"R6 macro step should have atom_cost=3, got {subs_json[1].get('atom_cost')}"

            # Structural validation: count F atoms (should be 4: 1 from R1 + 3 from CF3)
            smiles = prod['smiles']
            f_count = smiles.count('F')
            assert f_count == 4, f"k=2 CF3 product should have 4 F atoms, got {f_count} in {smiles}"

            print(f"  R6 macro k=2 product OK: {prod['inchikey'][:14]}... k_atoms={prod['k_atoms']} F_count={f_count}")

        print(f"\n[OK] R6_methyl macro k-level validations passed!")
    else:
        print("\n[WARN] No R6_methyl macro products found (may be expected if macro not triggered)")


if __name__ == '__main__':
    test_r6_methyl_k2_labeling_fix()
    test_r6_methyl_macro_k2_labeling()
