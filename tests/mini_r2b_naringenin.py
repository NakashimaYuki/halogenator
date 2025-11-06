# -*- coding: ascii -*-
"""
Minimal regression test for R2b site detection on naringenin (flavanone).

This test verifies that the R2b fallback enumeration correctly identifies
the C3 position (sp3 CH2 alpha to carbonyl) in naringenin.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

def test_r2b_naringenin_detection():
    """Test that R2b finds at least 1 site on naringenin (with fallback enabled)."""
    from halogenator.chem_compat import Chem
    from halogenator.sites import c_ring_sp3_CH2_flavanone_sites

    naringenin_smiles = "O=C1CC(c2ccc(O)cc2)Oc2cc(O)cc(O)c12"
    mol = Chem.MolFromSmiles(naringenin_smiles)

    assert mol is not None, "Failed to parse naringenin SMILES"

    # Test with sugar_cfg mode='off' and fallback enabled (default behavior)
    sugar_cfg = {'mode': 'off'}
    rules_cfg = {'R2': {'allow_alpha_as_beta': False, 'fallback': {'enable': True}}}

    r2b_sites = c_ring_sp3_CH2_flavanone_sites(mol, set(), sugar_cfg, rules_cfg)

    print(f"[R2b naringenin test] Found {len(r2b_sites)} sites: {r2b_sites}")
    assert len(r2b_sites) >= 1, f"Expected at least 1 R2b site on naringenin, got {len(r2b_sites)}"
    assert 2 in r2b_sites, f"Expected atom 2 (C3 position) in R2b sites, got {r2b_sites}"

    print("[OK] R2b naringenin detection test passed")


def test_r2b_naringenin_with_fallback():
    """Test that R2b fallback is triggered when strict detection fails."""
    from halogenator.chem_compat import Chem
    from halogenator.sites import (
        c_ring_sp3_CH2_flavanone_sites,
        c_ring_membership_atoms,
        _has_ring_oxygen_neighbor_on_same_c_ring,
        _is_beta_on_same_c_ring
    )

    naringenin_smiles = "O=C1CC(c2ccc(O)cc2)Oc2cc(O)cc(O)c12"
    mol = Chem.MolFromSmiles(naringenin_smiles)

    assert mol is not None, "Failed to parse naringenin SMILES"

    sugar_cfg = {'mode': 'off'}
    rules_cfg = {'R2': {'allow_alpha_as_beta': False, 'fallback': {'enable': True}}}

    # The C3 position (atom 2) should NOT pass strict checks
    # (it's alpha to carbonyl, not beta, and doesn't have ring oxygen neighbor)
    has_o_neighbor = _has_ring_oxygen_neighbor_on_same_c_ring(mol, 2)
    is_beta = _is_beta_on_same_c_ring(mol, 2, allow_alpha_as_beta=False)
    print(f"[Atom 2 strict checks] ring_O_neighbor={has_o_neighbor}, is_beta={is_beta}")
    assert not has_o_neighbor and not is_beta, "Atom 2 should fail strict checks"

    # But fallback should find it
    r2b_sites = c_ring_sp3_CH2_flavanone_sites(mol, set(), sugar_cfg, rules_cfg)
    print(f"[R2b with fallback] Found {len(r2b_sites)} sites: {r2b_sites}")
    assert len(r2b_sites) >= 1, "Fallback should find at least 1 site"
    assert 2 in r2b_sites, "Fallback should find atom 2"

    print("[OK] R2b fallback test passed")


def test_r2b_fallback_can_be_disabled():
    """Test that R2b fallback can be disabled via configuration."""
    from halogenator.chem_compat import Chem
    from halogenator.sites import c_ring_sp3_CH2_flavanone_sites

    naringenin_smiles = "O=C1CC(c2ccc(O)cc2)Oc2cc(O)cc(O)c12"
    mol = Chem.MolFromSmiles(naringenin_smiles)

    assert mol is not None, "Failed to parse naringenin SMILES"

    sugar_cfg = {'mode': 'off'}

    # Test 1: Fallback enabled (default) - should find sites
    rules_cfg_enabled = {
        'R2': {
            'allow_alpha_as_beta': False,
            'fallback': {'enable': True}
        }
    }
    sites_with_fallback = c_ring_sp3_CH2_flavanone_sites(mol, set(), sugar_cfg, rules_cfg_enabled)
    print(f"[Fallback enabled] Found {len(sites_with_fallback)} sites: {sites_with_fallback}")
    assert len(sites_with_fallback) >= 1, "With fallback enabled, should find at least 1 site"

    # Test 2: Fallback disabled - should find 0 sites (strict detection fails on naringenin)
    rules_cfg_disabled = {
        'R2': {
            'allow_alpha_as_beta': False,
            'fallback': {'enable': False}
        }
    }
    sites_without_fallback = c_ring_sp3_CH2_flavanone_sites(mol, set(), sugar_cfg, rules_cfg_disabled)
    print(f"[Fallback disabled] Found {len(sites_without_fallback)} sites: {sites_without_fallback}")
    assert len(sites_without_fallback) == 0, "With fallback disabled and strict detection failing, should find 0 sites"

    print("[OK] R2b fallback configuration test passed")


def test_r2b_detection_field_tracking():
    """Test that detection field correctly tracks 'strict' vs 'fallback' enumeration."""
    from halogenator.chem_compat import Chem
    from halogenator.sites import c_ring_sp3_CH2_flavanone_sites

    naringenin_smiles = "O=C1CC(c2ccc(O)cc2)Oc2cc(O)cc(O)c12"
    mol = Chem.MolFromSmiles(naringenin_smiles)

    assert mol is not None, "Failed to parse naringenin SMILES"

    sugar_cfg = {'mode': 'off'}
    rules_cfg = {
        'R2': {
            'allow_alpha_as_beta': False,
            'fallback': {'enable': True}
        }
    }

    # Request detection information
    sites, used_fallback = c_ring_sp3_CH2_flavanone_sites(
        mol, set(), sugar_cfg, rules_cfg, return_detection=True
    )

    print(f"[Detection tracking] Found {len(sites)} sites, used_fallback={used_fallback}")

    # For naringenin with allow_alpha_as_beta=False, strict should fail and fallback should be used
    assert len(sites) >= 1, "Should find at least 1 site"
    assert used_fallback == True, "Should use fallback detection for naringenin (strict fails)"

    print("[OK] R2b detection field tracking test passed")


def test_by_rule_family_shows_r2():
    """Test that by_rule aggregation by family shows R2 (not just R2b)."""
    import json
    import pandas as pd

    # Read latest products parquet
    products_path = "data/output/p1/products_k2.parquet"
    if not os.path.exists(products_path):
        print(f"[SKIP] Products file not found: {products_path}")
        return

    df = pd.read_parquet(products_path)

    # Test 1: by rule (default) - should show R2b
    if 'rule' in df.columns:
        rule_counts = df['rule'].value_counts()
        print(f"[By rule] Counts: {dict(rule_counts)}")
        assert 'R2b' in rule_counts.index, "Expected R2b in rule counts"

    # Test 2: by rule_family - should show R2
    if 'rule_family' in df.columns:
        family_counts = df['rule_family'].value_counts()
        print(f"[By family] Counts: {dict(family_counts)}")
        assert 'R2' in family_counts.index, "Expected R2 in rule_family counts"
        assert family_counts['R2'] > 0, "Expected R2 family to have >0 products"

    print("[OK] By-rule family aggregation test passed")


def test_k2_metadata_fields_present():
    """
    Test that k>=2 products contain sub_rule and detection fields.

    This test verifies the fix for the issue where k>=2 products were missing
    the sub_rule and detection metadata fields that k=1 products had.
    """
    import pandas as pd

    # Test with any recent naringenin raw output
    products_k2_path = "data/output/naringenin_k2_raw/products_k2.parquet"

    if not os.path.exists(products_k2_path):
        print(f"[SKIP] k=2 parquet not found: {products_k2_path}")
        print("       Run naringenin validation to generate test data")
        return

    df = pd.read_parquet(products_k2_path)

    print(f"[k=2 metadata] Loaded {len(df)} products from {products_k2_path}")
    print(f"[k=2 metadata] Columns: {sorted(df.columns)}")

    # Test 1: sub_rule field must be present
    assert 'sub_rule' in df.columns, "Expected 'sub_rule' field in k=2 parquet"

    # Test 2: detection field must be present
    assert 'detection' in df.columns, "Expected 'detection' field in k=2 parquet"

    # Test 3: Filter to R2 family and check sub_rule distribution
    if 'rule_family' in df.columns:
        r2_df = df[df['rule_family'] == 'R2']
        print(f"[k=2 metadata] R2 family: n={len(r2_df)}")

        if len(r2_df) > 0:
            # Check sub_rule values
            sub_rule_counts = r2_df['sub_rule'].value_counts(dropna=False)
            print(f"[k=2 metadata] R2 sub_rule distribution: {dict(sub_rule_counts)}")
            assert (r2_df['sub_rule'].isin(['R2a', 'R2b'])).all(), \
                "All R2 products should have sub_rule in ['R2a', 'R2b']"

            # Check detection values
            detection_counts = r2_df['detection'].value_counts(dropna=False)
            print(f"[k=2 metadata] R2 detection distribution: {dict(detection_counts)}")
            assert set(r2_df['detection'].unique()) <= {'strict', 'fallback'}, \
                "All R2 products should have detection in ['strict', 'fallback']"

            # For raw mode with naringenin, expect mainly fallback for R2b
            if (r2_df['sub_rule'] == 'R2b').any():
                r2b_df = r2_df[r2_df['sub_rule'] == 'R2b']
                print(f"[k=2 metadata] R2b detection: {dict(r2b_df['detection'].value_counts())}")

    print("[OK] k>=2 metadata fields test passed")


if __name__ == '__main__':
    print("=" * 80)
    print("Minimal R2b Regression Test for Naringenin")
    print("=" * 80)

    try:
        test_r2b_naringenin_detection()
        test_r2b_naringenin_with_fallback()
        test_r2b_fallback_can_be_disabled()
        test_r2b_detection_field_tracking()
        test_by_rule_family_shows_r2()
        test_k2_metadata_fields_present()  # New test for k>=2 metadata

        print("\n" + "=" * 80)
        print("ALL TESTS PASSED")
        print("=" * 80)
    except AssertionError as e:
        print(f"\n[FAIL] {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n[ERROR] {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
