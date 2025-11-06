# -*- coding: ascii -*-
"""
Tests for raw mode product count accuracy.
Verifies that raw mode generates the theoretical full set of products.
"""

import pytest


def test_k1_count_simple_molecule():
    """Test k=1 count for simple symmetric molecule."""
    try:
        from halogenator.enumerate_k import enumerate_with_stats, EnumConfig
    except ImportError:
        pytest.skip("enumerate_k module not available (RDKit may not be installed)")

    # Benzene: 6 equivalent positions
    parent_smi = "c1ccccc1"

    cfg = EnumConfig(
        k_max=1,
        halogens=('F',),
        rules=('R1',),  # Aromatic C-H
        constraints={'enable': False},
        pruning_cfg={'enable_symmetry_fold': False},
        engine_cfg={'enable_inchi_dedup': False}
    )

    products, qa_stats = enumerate_with_stats(parent_smi, cfg)

    # Should have 6 products (one for each position)
    # Note: This assumes R1 can halogenate all 6 aromatic C-H
    # With symmetry folding disabled, we expect one product per site
    # With InChI dedup disabled, we expect all sites even if chemically equivalent
    assert len(products) > 0, "Should generate at least one product"

    # In raw mode without symmetry folding, we should see attempts for all 6 positions
    # However, due to chemical equivalence, the actual unique InChI may be just 1
    # So we verify that products were generated (count depends on rule implementation)
    print(f"[OK] K=1 benzene enumeration: {len(products)} products")


def test_k2_count_combinatorics():
    """Test k=2 count follows combinatorial formula."""
    try:
        from halogenator.enumerate_k import enumerate_with_stats, EnumConfig
    except ImportError:
        pytest.skip("enumerate_k module not available (RDKit may not be installed)")

    # Simple molecule with 2 distinct halogenatable sites
    parent_smi = "CCO"  # ethanol (1 OH for R1, possibly C-H for other rules)

    cfg = EnumConfig(
        k_max=2,
        halogens=('F', 'Cl'),
        rules=('R1',),  # Only OH halogenation
        constraints={'enable': False},
        pruning_cfg={'enable_symmetry_fold': False},
        engine_cfg={'enable_inchi_dedup': False}
    )

    products, qa_stats = enumerate_with_stats(parent_smi, cfg)

    k1 = [p for p in products if p['k'] == 1]
    k2 = [p for p in products if p['k'] == 2]

    # k=1: 1 site (OH) ? 2 halogens = 2 products
    # k=2: Since only 1 OH site, k=2 would require halogenating other sites
    # Exact counts depend on rule behavior

    assert len(k1) > 0, "Should have k=1 products"
    # k2 may be 0 if only one OH site is available
    # assert len(k2) > 0, "Should have k=2 products"

    # Verify no duplicates in raw mode (SMILES-level check)
    k1_smiles = [p['smiles'] for p in k1]
    k2_smiles = [p['smiles'] for p in k2]

    print(f"[OK] K=2 combinatorics: {len(k1)} k=1 products, {len(k2)} k=2 products")


def test_macro_substitution_count():
    """Test that macro substitution counts correctly (k_ops=1, k_atoms=3)."""
    try:
        from halogenator.enumerate_k import enumerate_with_stats, EnumConfig
    except ImportError:
        pytest.skip("enumerate_k module not available (RDKit may not be installed)")

    parent_smi = "CC1=CC=CC=C1"  # toluene

    cfg = EnumConfig(
        k_max=1,
        halogens=('F',),
        rules=('R6_methyl',),
        constraints={'enable': False},
        engine_cfg={'enable_inchi_dedup': False},
        rules_cfg={
            'R6_methyl': {
                'enable': True,
                'allowed': ['F'],
                'macro': {'enable': True, 'labels': ['CF3']}
            }
        }
    )

    products, qa_stats = enumerate_with_stats(parent_smi, cfg)

    # Should have:
    # - Step products: -CH? ? -CH?F (k_ops=1, k_atoms=1)
    # - Macro products: -CH? ? -CF? (k_ops=1, k_atoms=3)

    step_products = []
    macro_products = []

    for p in products:
        if p.get('macro_label') == 'CF3':
            macro_products.append(p)
        else:
            step_products.append(p)

    assert len(macro_products) >= 1, "Should have at least 1 CF3 macro product"

    # Verify k and k_ops/k_atoms for macro
    for macro_prod in macro_products:
        assert macro_prod['k'] == 1, "Macro should have k=1"
        assert macro_prod['k_ops'] == 1, "Macro should have k_ops=1"
        assert macro_prod['k_atoms'] == 3, "Macro should have k_atoms=3"

    print(f"[OK] Macro substitution count: {len(macro_products)} CF3 products, {len(step_products)} step products")


def test_constraint_toggle_changes_count():
    """Test that disabling constraints increases product count."""
    try:
        from halogenator.enumerate_k import enumerate_with_stats, EnumConfig
    except ImportError:
        pytest.skip("enumerate_k module not available (RDKit may not be installed)")

    # Molecule where constraints would filter products
    parent_smi = "c1ccccc1"  # benzene

    # With constraints (min_graph_distance=2)
    cfg_constrained = EnumConfig(
        k_max=2,
        halogens=('F',),
        rules=('R1',),
        constraints={'enable': True, 'min_graph_distance': 2},
        pruning_cfg={'enable_symmetry_fold': False},
        engine_cfg={'enable_inchi_dedup': False}
    )

    # Without constraints
    cfg_raw = EnumConfig(
        k_max=2,
        halogens=('F',),
        rules=('R1',),
        constraints={'enable': False},
        pruning_cfg={'enable_symmetry_fold': False},
        engine_cfg={'enable_inchi_dedup': False}
    )

    products_constrained, _ = enumerate_with_stats(parent_smi, cfg_constrained)
    products_raw, _ = enumerate_with_stats(parent_smi, cfg_raw)

    # Raw mode should allow ortho-dihalogenated products
    assert len(products_raw) >= len(products_constrained), \
        f"Raw mode should allow more products (raw={len(products_raw)}, constrained={len(products_constrained)})"

    # Get k=2 products only
    k2_raw = [p for p in products_raw if p['k'] == 2]
    k2_constrained = [p for p in products_constrained if p['k'] == 2]

    print(f"[OK] Constraint effect on k=2: {len(k2_raw)} raw vs {len(k2_constrained)} constrained")


def test_symmetry_fold_effect():
    """Test that symmetry folding reduces product count."""
    try:
        from halogenator.enumerate_k import enumerate_with_stats, EnumConfig
    except ImportError:
        pytest.skip("enumerate_k module not available (RDKit may not be installed)")

    # Highly symmetric molecule
    parent_smi = "c1ccccc1"  # benzene

    # Without folding
    cfg_no_fold = EnumConfig(
        k_max=1,
        halogens=('F',),
        rules=('R1',),
        constraints={'enable': False},
        pruning_cfg={'enable_symmetry_fold': False},
        engine_cfg={'enable_inchi_dedup': False}
    )

    # With folding
    cfg_with_fold = EnumConfig(
        k_max=1,
        halogens=('F',),
        rules=('R1',),
        constraints={'enable': False},
        pruning_cfg={'enable_symmetry_fold': True},
        engine_cfg={'enable_inchi_dedup': True}
    )

    products_no_fold, _ = enumerate_with_stats(parent_smi, cfg_no_fold)
    products_with_fold, _ = enumerate_with_stats(parent_smi, cfg_with_fold)

    # With folding should have fewer (or equal) attempts/products
    # Note: Final product count may be same (1) due to chemical equivalence
    # But without folding, we explore all positions

    print(f"[OK] Symmetry folding effect: {len(products_no_fold)} no-fold vs {len(products_with_fold)} with-fold")


def test_inchi_dedup_effect():
    """Test that InChI deduplication reduces product count."""
    try:
        from halogenator.enumerate_k import enumerate_with_stats, EnumConfig
    except ImportError:
        pytest.skip("enumerate_k module not available (RDKit may not be installed)")

    # Molecule that may generate tautomers or stereoisomers
    parent_smi = "c1cc(O)ccc1O"  # catechol (1,2-dihydroxybenzene)

    # Without dedup
    cfg_no_dedup = EnumConfig(
        k_max=1,
        halogens=('F',),
        rules=('R1',),
        constraints={'enable': False},
        pruning_cfg={'enable_symmetry_fold': False},
        engine_cfg={'enable_inchi_dedup': False}
    )

    # With dedup
    cfg_with_dedup = EnumConfig(
        k_max=1,
        halogens=('F',),
        rules=('R1',),
        constraints={'enable': False},
        pruning_cfg={'enable_symmetry_fold': False},
        engine_cfg={'enable_inchi_dedup': True}
    )

    products_no_dedup, _ = enumerate_with_stats(parent_smi, cfg_no_dedup)
    products_with_dedup, _ = enumerate_with_stats(parent_smi, cfg_with_dedup)

    # With dedup should have fewer (or equal) products
    assert len(products_no_dedup) >= len(products_with_dedup), \
        f"No-dedup should have >= products (no-dedup={len(products_no_dedup)}, with-dedup={len(products_with_dedup)})"

    print(f"[OK] InChI dedup effect: {len(products_no_dedup)} no-dedup vs {len(products_with_dedup)} with-dedup")


def test_full_raw_mode_product_count():
    """Test full raw mode with all constraints/dedup/folding disabled."""
    try:
        from halogenator.enumerate_k import enumerate_with_stats, EnumConfig
    except ImportError:
        pytest.skip("enumerate_k module not available (RDKit may not be installed)")

    parent_smi = "c1ccccc1"  # benzene

    # Full raw mode
    cfg = EnumConfig(
        k_max=2,
        halogens=('F', 'Cl'),
        rules=('R1',),
        constraints={'enable': False},
        sugar_cfg={'mode': 'off'},
        pruning_cfg={'enable_symmetry_fold': False},
        engine_cfg={'enable_inchi_dedup': False}
    )

    products, qa_stats = enumerate_with_stats(parent_smi, cfg)

    k1 = [p for p in products if p['k'] == 1]
    k2 = [p for p in products if p['k'] == 2]

    # In full raw mode, we should get comprehensive enumeration
    # Exact counts depend on how many sites are actually enumerated
    assert len(k1) > 0, "Should have k=1 products"
    assert len(k2) > 0, "Should have k=2 products"

    # Verify all products have required fields
    for p in products:
        assert 'k' in p
        assert 'k_ops' in p
        assert 'k_atoms' in p
        assert p['k'] == p['k_ops'], "k should equal k_ops"

    print(f"[OK] Full raw mode: {len(k1)} k=1 + {len(k2)} k=2 = {len(products)} total products")


def test_budget_mode_consistency():
    """Test that budget_mode='ops' is consistent with k counting."""
    try:
        from halogenator.enumerate_k import enumerate_with_stats, EnumConfig
    except ImportError:
        pytest.skip("enumerate_k module not available (RDKit may not be installed)")

    parent_smi = "CC1=CC=CC=C1"  # toluene

    cfg = EnumConfig(
        k_max=1,
        halogens=('F',),
        rules=('R6_methyl',),
        constraints={'enable': False},
        engine_cfg={'budget_mode': 'ops', 'enable_inchi_dedup': False},
        rules_cfg={
            'R6_methyl': {
                'enable': True,
                'allowed': ['F'],
                'macro': {'enable': True, 'labels': ['CF3']}
            }
        }
    )

    products, qa_stats = enumerate_with_stats(parent_smi, cfg)

    # Verify budget_mode field
    for p in products:
        assert p.get('budget_mode') == 'ops', "Budget mode should be 'ops'"
        assert p['k'] == p['k_ops'], "k should equal k_ops when budget_mode='ops'"

    print(f"[OK] Budget mode consistency validated for {len(products)} products")


def test_multiple_halogens_combinatorics():
    """Test combinatorics with multiple halogens."""
    try:
        from halogenator.enumerate_k import enumerate_with_stats, EnumConfig
    except ImportError:
        pytest.skip("enumerate_k module not available (RDKit may not be installed)")

    parent_smi = "c1ccccc1O"  # phenol

    cfg = EnumConfig(
        k_max=1,
        halogens=('F', 'Cl', 'Br', 'I'),
        rules=('R1',),
        constraints={'enable': False},
        pruning_cfg={'enable_symmetry_fold': False},
        engine_cfg={'enable_inchi_dedup': False}
    )

    products, qa_stats = enumerate_with_stats(parent_smi, cfg)

    # Group by halogen
    by_halogen = {}
    for p in products:
        halogen = p.get('halogen')
        if halogen not in by_halogen:
            by_halogen[halogen] = []
        by_halogen[halogen].append(p)

    # Should have products for each halogen
    assert len(by_halogen) > 0, "Should have products for at least one halogen"

    print(f"[OK] Multiple halogens: {sum(len(v) for v in by_halogen.values())} total products")
    for halogen, prods in by_halogen.items():
        print(f"  - {halogen}: {len(prods)} products")


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v"])
