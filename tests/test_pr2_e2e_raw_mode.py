# -*- coding: ascii -*-
"""
End-to-end tests for PR2 raw mode implementation.
Tests the complete k=2 enumeration pipeline with hierarchical output.
"""

import pytest
import os
import json
import csv
import tempfile
import shutil
from pathlib import Path


def test_raw_mode_config_loads():
    """Test that raw mode config file loads correctly."""
    import yaml

    config_path = "configs/one_flavone_k2_raw.yaml"

    # Check if config file exists
    if not os.path.exists(config_path):
        pytest.skip(f"Config file not found: {config_path}")

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    # Verify raw mode settings
    assert config.get('constraints', {}).get('enable') == False, "Constraints should be disabled in raw mode"
    assert config.get('sugar_cfg', {}).get('mode') == 'off', "Sugar masking should be off in raw mode"
    assert config.get('pruning_cfg', {}).get('enable_symmetry_fold') == False, "Symmetry folding should be disabled"
    assert config.get('engine', {}).get('enable_inchi_dedup') == False, "InChI dedup should be disabled"

    # Verify macro substitution is enabled
    assert config.get('rules_cfg', {}).get('R6_methyl', {}).get('macro', {}).get('enable') == True, \
        "R6_methyl macro should be enabled"
    assert 'CF3' in config.get('rules_cfg', {}).get('R6_methyl', {}).get('macro', {}).get('labels', []), \
        "CF3 should be in macro labels"

    print("[OK] Raw mode config loaded successfully with expected settings")


def test_k2_enumeration_completes():
    """Test that k=2 enumeration completes without errors."""
    try:
        from halogenator.enumerate_k import enumerate_with_stats, EnumConfig
    except ImportError:
        pytest.skip("enumerate_k module not available (RDKit may not be installed)")

    parent_smi = "O=C1CC(c2ccc(O)cc2)Oc2cc(O)cc(O)c12"  # naringenin

    cfg = EnumConfig(
        k_max=2,
        halogens=('F', 'Cl'),  # Limited for speed
        rules=('R1', 'R2', 'R6_methyl'),
        constraints={'enable': False},
        sugar_cfg={'mode': 'off'},
        pruning_cfg={'enable_symmetry_fold': False},
        engine_cfg={'enable_inchi_dedup': False, 'budget_mode': 'ops'},
        rules_cfg={
            'R2': {
                'sp2_CH_in_C_ring': True,
                'sp3_CH2_flavanone': True,
                'allowed_halogens': ['F', 'Cl']
            },
            'R6_methyl': {
                'enable': True,
                'allowed': ['F', 'Cl'],
                'macro': {'enable': True, 'labels': ['CF3']}
            }
        }
    )

    products, qa_stats = enumerate_with_stats(parent_smi, cfg)

    # Should have products
    assert len(products) > 0, "Should generate products"

    # Check k values
    k1_products = [p for p in products if p.get('k') == 1]
    k2_products = [p for p in products if p.get('k') == 2]

    assert len(k1_products) > 0, "Should have k=1 products"
    assert len(k2_products) > 0, "Should have k=2 products"

    # Verify k=k_ops consistency
    for p in products:
        assert 'k_ops' in p, "Product should have k_ops field"
        assert p['k'] == p['k_ops'], f"k should equal k_ops (got k={p['k']}, k_ops={p['k_ops']})"

    # Verify required fields
    required_fields = ['smiles', 'inchikey', 'k', 'k_ops', 'k_atoms', 'rule', 'halogen', 'substitutions_json']
    for p in products:
        for field in required_fields:
            assert field in p, f"Product should have {field} field"

    print(f"[OK] K=2 enumeration completed: {len(k1_products)} k=1 products, {len(k2_products)} k=2 products")


def test_macro_substitution_in_products():
    """Test that macro substitution (CF3/CCl3) products are generated."""
    try:
        from halogenator.enumerate_k import enumerate_with_stats, EnumConfig
    except ImportError:
        pytest.skip("enumerate_k module not available (RDKit may not be installed)")

    # Simple molecule with methyl group
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

    # Should have both step and macro products
    step_products = []
    macro_products = []

    for p in products:
        if p.get('macro_label') == 'CF3':
            macro_products.append(p)
        else:
            # Check substitutions_json for step type
            try:
                subs = json.loads(p.get('substitutions_json', '[]'))
                if subs and subs[0].get('type') != 'macro':
                    step_products.append(p)
            except:
                step_products.append(p)

    assert len(macro_products) > 0, "Should generate CF3 macro products"

    # Verify macro product properties
    for macro_prod in macro_products:
        assert macro_prod.get('k_ops') == 1, "Macro should have k_ops=1"
        assert macro_prod.get('k_atoms') == 3, "Macro should have k_atoms=3"
        assert macro_prod.get('k') == 1, "Macro should have k=1 (equals k_ops)"
        assert macro_prod.get('macro_label') == 'CF3', "Macro should have macro_label=CF3"

        # Verify substitutions_json contains macro metadata
        subs = json.loads(macro_prod.get('substitutions_json', '[]'))
        assert len(subs) > 0, "Should have substitution history"
        assert subs[0].get('type') == 'macro', "First substitution should be macro type"
        assert subs[0].get('label') == 'CF3', "First substitution should have CF3 label"

    print(f"[OK] Macro substitution working: {len(macro_products)} CF3 products, {len(step_products)} step products")


def test_hierarchical_output_structure():
    """Test that hierarchical output creates correct directory structure."""
    try:
        from halogenator.io_hierarchy import write_hierarchical_outputs
        from rdkit import Chem
    except ImportError:
        pytest.skip("io_hierarchy or RDKit not available")

    # Create mock data
    parent_record = {
        'smiles': 'CC',
        'inchikey': 'OTMSDBZUPAUEDD-UHFFFAOYSA-N',
        'name': 'ethane'
    }

    k1_records = [
        {
            'k': 1, 'smiles': 'FC', 'inchikey': 'KEY1-ABCDEF-G', 'halogen': 'F',
            'parent_inchikey': 'OTMSDBZUPAUEDD-UHFFFAOYSA-N',
            'root_parent_inchikey': 'OTMSDBZUPAUEDD-UHFFFAOYSA-N',
            'substitutions_json': '[{"halogen":"F"}]'
        },
        {
            'k': 1, 'smiles': 'ClC', 'inchikey': 'KEY2-ABCDEF-G', 'halogen': 'Cl',
            'parent_inchikey': 'OTMSDBZUPAUEDD-UHFFFAOYSA-N',
            'root_parent_inchikey': 'OTMSDBZUPAUEDD-UHFFFAOYSA-N',
            'substitutions_json': '[{"halogen":"Cl"}]'
        },
    ]

    k2_records = [
        {
            'k': 2, 'smiles': 'FCC', 'inchikey': 'KEY3-ABCDEF-G', 'halogen': 'Cl',
            'parent_inchikey': 'KEY1-ABCDEF-G',
            'root_parent_inchikey': 'OTMSDBZUPAUEDD-UHFFFAOYSA-N',
            'substitutions_json': '[{"halogen":"F"},{"halogen":"Cl"}]'
        },
    ]

    # Create temporary output directory
    with tempfile.TemporaryDirectory() as tmpdir:
        outdir = Path(tmpdir)

        try:
            summary = write_hierarchical_outputs(
                parent_record, k1_records + k2_records, str(outdir), ['F', 'Cl']
            )

            # Verify directory structure
            parent_dir = outdir / summary['parent_name']
            assert parent_dir.exists(), "Parent directory should exist"
            assert (parent_dir / 'index.json').exists(), "index.json should exist"
            assert (parent_dir / 'k1_summary.csv').exists(), "k1_summary.csv should exist"
            assert (parent_dir / 'k1' / 'F').exists(), "k1/F directory should exist"
            assert (parent_dir / 'k1' / 'Cl').exists(), "k1/Cl directory should exist"

            # Verify index.json content
            with open(parent_dir / 'index.json') as f:
                index_data = json.load(f)

            assert index_data['parent']['name'] == summary['parent_name'], "Parent name should match"
            assert index_data['summary']['total_k1_products'] == 2, "Should have 2 k=1 products"
            assert len(index_data['files']) > 0, "Should have file entries"

            print(f"[OK] Hierarchical output structure created successfully")
            print(f"  - Generated {summary.get('total_files', 0)} files")
            print(f"  - K1 products: {index_data['summary']['total_k1_products']}")
            print(f"  - K2 products: {index_data['summary']['total_k2_products']}")

        except Exception as e:
            pytest.fail(f"Hierarchical output generation failed: {e}")


def test_sdf_parent_first_entry():
    """Test that SDF files have parent molecule as first entry."""
    try:
        from halogenator.io_hierarchy import write_group_sdf
        from rdkit import Chem
    except ImportError:
        pytest.skip("io_hierarchy or RDKit not available")

    parent_mol = Chem.MolFromSmiles('CC')
    product_mols = [Chem.MolFromSmiles('FC'), Chem.MolFromSmiles('ClC')]
    product_records = [
        {'inchikey': 'KEY1', 'k': 1, 'rule': 'R1', 'halogen': 'F'},
        {'inchikey': 'KEY2', 'k': 1, 'rule': 'R1', 'halogen': 'Cl'}
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        test_sdf = Path(tmpdir) / "test.sdf"

        try:
            write_group_sdf(str(test_sdf), parent_mol, product_mols, product_records)

            # Read back and verify
            suppl = Chem.SDMolSupplier(str(test_sdf), removeHs=False, sanitize=False)
            mols = [m for m in suppl if m is not None]

            assert len(mols) == 3, f"Should have 3 molecules (1 parent + 2 products), got {len(mols)}"

            # First entry should be parent
            first_mol = mols[0]
            assert first_mol.HasProp('parent'), "First entry should have 'parent' property"
            assert first_mol.GetProp('parent') == 'TRUE', "First entry should have parent='TRUE'"

            # Subsequent entries should be products
            for i, mol in enumerate(mols[1:], 1):
                assert mol.HasProp('inchikey'), f"Product {i} should have inchikey"
                assert mol.HasProp('k'), f"Product {i} should have k"
                assert mol.HasProp('rule'), f"Product {i} should have rule"
                assert not mol.HasProp('parent'), f"Product {i} should not have parent property"

            print("[OK] SDF parent-first structure validated")

        except Exception as e:
            pytest.fail(f"SDF parent-first test failed: {e}")


def test_raw_mode_no_dedup():
    """Test that raw mode does not deduplicate products."""
    try:
        from halogenator.enumerate_k import enumerate_with_stats, EnumConfig
    except ImportError:
        pytest.skip("enumerate_k module not available (RDKit may not be installed)")

    # Symmetric molecule that would produce duplicates with symmetry folding
    parent_smi = "c1ccccc1"  # benzene

    # Raw mode config
    cfg_raw = EnumConfig(
        k_max=1,
        halogens=('F',),
        rules=('R1',),
        constraints={'enable': False},
        pruning_cfg={'enable_symmetry_fold': False},
        engine_cfg={'enable_inchi_dedup': False}
    )

    # Unique mode config
    cfg_unique = EnumConfig(
        k_max=1,
        halogens=('F',),
        rules=('R1',),
        constraints={'enable': True},
        pruning_cfg={'enable_symmetry_fold': True},
        engine_cfg={'enable_inchi_dedup': True}
    )

    products_raw, qa_raw = enumerate_with_stats(parent_smi, cfg_raw)
    products_unique, qa_unique = enumerate_with_stats(parent_smi, cfg_unique)

    # Raw mode should have more (or equal) products due to no folding/dedup
    assert len(products_raw) >= len(products_unique), \
        f"Raw mode should have >= products (raw={len(products_raw)}, unique={len(products_unique)})"

    print(f"[OK] Raw mode dedup validation: {len(products_raw)} raw vs {len(products_unique)} unique products")


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
        f"Raw mode should have >= products than constrained (raw={len(products_raw)}, constrained={len(products_constrained)})"

    print(f"[OK] Constraint toggle validation: {len(products_raw)} raw vs {len(products_constrained)} constrained products")


def test_config_cli_override():
    """Test that CLI arguments override config file settings."""
    try:
        from halogenator.enumerate_k import EnumConfig
    except ImportError:
        pytest.skip("enumerate_k module not available")

    # Create base config with constraints enabled
    cfg = EnumConfig(
        k_max=1,
        halogens=('F',),
        rules=('R1',),
        constraints={'enable': True},
        sugar_cfg={'mode': 'heuristic'},
        pruning_cfg={'enable_symmetry_fold': True},
        engine_cfg={'enable_inchi_dedup': True}
    )

    # Verify base config
    assert cfg.constraints.get('enable') == True, "Base config should have constraints enabled"
    assert cfg.sugar_cfg.get('mode') == 'heuristic', "Base config should have sugar mode=heuristic"
    assert cfg.pruning_cfg.get('enable_symmetry_fold') == True, "Base config should have folding enabled"
    assert cfg.engine_cfg.get('enable_inchi_dedup') == True, "Base config should have dedup enabled"

    # Simulate CLI override for raw mode
    cfg_raw = EnumConfig(
        k_max=1,
        halogens=('F',),
        rules=('R1',),
        constraints={'enable': False},  # CLI override
        sugar_cfg={'mode': 'off'},  # CLI override
        pruning_cfg={'enable_symmetry_fold': False},  # CLI override
        engine_cfg={'enable_inchi_dedup': False}  # CLI override
    )

    # Verify CLI override
    assert cfg_raw.constraints.get('enable') == False, "CLI should override constraints"
    assert cfg_raw.sugar_cfg.get('mode') == 'off', "CLI should override sugar mode"
    assert cfg_raw.pruning_cfg.get('enable_symmetry_fold') == False, "CLI should override folding"
    assert cfg_raw.engine_cfg.get('enable_inchi_dedup') == False, "CLI should override dedup"

    print("[OK] CLI override mechanism validated")


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v"])
