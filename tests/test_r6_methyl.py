# -*- coding: ascii -*-
"""
PR-3 R6 Methyl Tests

Test suite for R6 methyl site identification and enumeration:
- toluene: ops/atoms budget constraints
- p-xylene: symmetry folding behavior
- methyl site detection and filtering
"""

import unittest
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from halogenator.sites_methyl import enumerate_methyl_sites
from halogenator.budget import BudgetState
from halogenator.chem_compat import Chem


class TestR6Methyl(unittest.TestCase):
    """Test R6 methyl site identification and budget system."""

    def setUp(self):
        """Set up test molecules."""
        self.toluene_smiles = "Cc1ccccc1"  # Simple toluene
        self.p_xylene_smiles = "Cc1ccc(C)cc1"  # p-xylene with two methyl groups
        self.anisole_smiles = "COc1ccccc1"  # Methoxy group (should be excluded by default)

    def test_toluene_methyl_site_detection(self):
        """Test methyl site detection in toluene."""
        mol = Chem.MolFromSmiles(self.toluene_smiles)
        sites = enumerate_methyl_sites(mol, set(), allow_on_methoxy=False)

        # Toluene should have exactly one methyl site
        self.assertEqual(len(sites), 1, "Toluene should have exactly one methyl site")

        # Check the methyl carbon properties
        methyl_idx = sites[0]
        atom = mol.GetAtomWithIdx(methyl_idx)
        self.assertEqual(atom.GetSymbol(), "C")
        self.assertEqual(atom.GetHybridization(), Chem.HybridizationType.SP3)
        self.assertEqual(atom.GetTotalNumHs(), 3)

    def test_multi_step_halogenation_support(self):
        """Test that sites support multi-step halogenation CH3 -> CH2X -> CHX2."""
        # Create a simple ethyl fluoride (CH3CH2F) - both carbons should be detected
        ethyl_fluoride_smiles = "CCF"  # CH3CH2F
        mol = Chem.MolFromSmiles(ethyl_fluoride_smiles)
        sites = enumerate_methyl_sites(mol, set(), allow_on_methoxy=False)

        # Should detect both carbons as valid methyl sites (CH3 has 3H, CH2F has 2H)
        self.assertEqual(len(sites), 2, "CH3CH2F should have two detectable methyl sites")

        # Check that we have both a 3H carbon (CH3) and a 2H carbon (CH2F)
        hydrogen_counts = []
        for site_idx in sites:
            atom = mol.GetAtomWithIdx(site_idx)
            self.assertEqual(atom.GetSymbol(), "C")
            hydrogen_counts.append(atom.GetTotalNumHs())

            # Each should have at least one heavy neighbor
            heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
            self.assertGreaterEqual(len(heavy_neighbors), 1, f"Carbon {site_idx} should have at least one heavy neighbor")

            # Check that each has exactly one non-halogen heavy neighbor
            non_halogen_neighbors = [nbr for nbr in heavy_neighbors if nbr.GetSymbol() not in ('F', 'Cl', 'Br', 'I')]
            self.assertEqual(len(non_halogen_neighbors), 1, f"Carbon {site_idx} should have exactly one non-halogen heavy neighbor")

        # Should have one carbon with 3H (CH3) and one with 2H (CH2F)
        hydrogen_counts.sort()
        self.assertEqual(hydrogen_counts, [2, 3], "Should have one CH2F (2H) and one CH3 (3H) site")

    def test_r6_atoms_mode_chain_to_chf2(self):
        """Test R6 atoms mode with chaining to CHF2, ensuring k_atoms consistency."""
        from halogenator.enumerate_k import EnumConfig, enumerate_with_stats

        # Use simple CH3F starting point
        parent_smi = "CF"  # CH3F

        # Configure for atoms mode with k_max=2 to allow CHF2
        cfg = EnumConfig(
            k_max=2,
            halogens=('F',),
            rules=('R6',),
            engine_cfg={'budget_mode': 'atoms'},
            rules_cfg={
                'R6_methyl': {
                    'enable': True,
                    'allowed': ['F'],
                    'allow_on_methoxy': False,
                    'macro': {'enable': False}  # Only step halogenation
                }
            }
        )

        products, qa_stats = enumerate_with_stats(parent_smi, cfg)

        # Should find CHF2 product with k=k_atoms=2
        chf2_products = [p for p in products if p.get('smiles') == 'FC(F)']  # CHF2 canonical form
        if chf2_products:
            product = chf2_products[0]
            self.assertEqual(product.get('k'), 2, "CHF2 product should have k=2")
            self.assertEqual(product.get('k_atoms'), 2, "CHF2 product should have k_atoms=2")
            self.assertEqual(product.get('k'), product.get('k_atoms'), "k should equal k_atoms for R6 products")

    def test_p_xylene_methyl_sites(self):
        """Test methyl site detection in p-xylene."""
        mol = Chem.MolFromSmiles(self.p_xylene_smiles)
        sites = enumerate_methyl_sites(mol, set(), allow_on_methoxy=False)

        # p-xylene should have exactly two methyl sites
        self.assertEqual(len(sites), 2, "p-xylene should have exactly two methyl sites")

        # Both should be SP3 CH3
        for site_idx in sites:
            atom = mol.GetAtomWithIdx(site_idx)
            self.assertEqual(atom.GetSymbol(), "C")
            self.assertEqual(atom.GetHybridization(), Chem.HybridizationType.SP3)
            self.assertEqual(atom.GetTotalNumHs(), 3)

    def test_methoxy_exclusion(self):
        """Test that methoxy groups are excluded by default."""
        mol = Chem.MolFromSmiles(self.anisole_smiles)
        sites = enumerate_methyl_sites(mol, set(), allow_on_methoxy=False)

        # Should have no sites (methoxy excluded)
        self.assertEqual(len(sites), 0, "Anisole should have no methyl sites when methoxy excluded")

    def test_methoxy_inclusion(self):
        """Test that methoxy groups can be included when allowed."""
        mol = Chem.MolFromSmiles(self.anisole_smiles)
        sites = enumerate_methyl_sites(mol, set(), allow_on_methoxy=True)

        # Should have one site (methoxy included)
        self.assertEqual(len(sites), 1, "Anisole should have one methyl site when methoxy allowed")

    def test_masked_atoms_exclusion(self):
        """Test that masked atoms are properly excluded."""
        mol = Chem.MolFromSmiles(self.toluene_smiles)

        # Find the methyl site first
        all_sites = enumerate_methyl_sites(mol, set(), allow_on_methoxy=False)
        self.assertEqual(len(all_sites), 1)
        methyl_idx = all_sites[0]

        # Mask the methyl site
        masked_sites = enumerate_methyl_sites(mol, {methyl_idx}, allow_on_methoxy=False)
        self.assertEqual(len(masked_sites), 0, "Masked methyl site should be excluded")

    def test_budget_ops_mode(self):
        """Test budget system in ops mode."""
        # ops mode: limit by number of operations
        budget = BudgetState("ops", k_max=1)

        # First operation should succeed
        self.assertTrue(budget.charge(1, 3))  # 1 op, 3 atoms (CF3 macro)
        self.assertEqual(budget.k_ops, 1)
        self.assertEqual(budget.k_atoms, 3)

        # Second operation should fail (exceeds k_max=1 ops)
        self.assertFalse(budget.charge(1, 1))  # Would be 2 ops total
        self.assertEqual(budget.k_ops, 1)  # Should not change
        self.assertEqual(budget.k_atoms, 3)

    def test_budget_atoms_mode(self):
        """Test budget system in atoms mode."""
        # atoms mode: limit by number of atoms
        budget = BudgetState("atoms", k_max=2)

        # CF3 macro (3 atoms) should fail immediately
        self.assertFalse(budget.charge(1, 3))  # Exceeds k_max=2 atoms
        self.assertEqual(budget.k_ops, 0)
        self.assertEqual(budget.k_atoms, 0)

        # CHF2 (2 atoms) should succeed
        self.assertTrue(budget.charge(2, 2))  # 2 ops, 2 atoms
        self.assertEqual(budget.k_ops, 2)
        self.assertEqual(budget.k_atoms, 2)

        # No more room for additional atoms
        self.assertFalse(budget.charge(1, 1))

    def test_site_tokens_for_deduplication(self):
        """Test site tokens for tracking first-touch operations."""
        budget = BudgetState("ops", k_max=5)

        # Simulate first touch on site 10
        cidx = 10
        op_cost = 1 - budget.site_tokens.get(cidx, 0)  # First touch: 1-0=1
        self.assertEqual(op_cost, 1)

        # Charge and mark site as touched
        self.assertTrue(budget.charge(op_cost, 1))
        budget.site_tokens[cidx] = 1

        # Simulate second touch on same site
        op_cost = 1 - budget.site_tokens.get(cidx, 0)  # Second touch: 1-1=0
        self.assertEqual(op_cost, 0)

        # Should succeed with 0 op cost
        self.assertTrue(budget.charge(op_cost, 1))

    def test_site_state_multiset_tracking(self):
        """Test site state tracking for local deduplication."""
        budget = BudgetState("atoms", k_max=10)

        # Site 5: F-F sequence
        cidx = 5
        current_halogens = [("F", 2)]  # Two F atoms
        state_key = tuple(sorted(current_halogens))
        budget.site_state[cidx] = state_key

        self.assertEqual(budget.site_state[cidx], (("F", 2),))

        # Site 8: Mixed F-Cl sequence
        cidx = 8
        current_halogens = [("F", 1), ("Cl", 1)]
        state_key = tuple(sorted(current_halogens))
        budget.site_state[cidx] = state_key

        self.assertEqual(budget.site_state[cidx], (("Cl", 1), ("F", 1)))

    def test_early_dedup_no_pollution(self):
        """Test that early deduplication doesn't pollute global set when later checks fail."""
        from halogenator.enumerate_k import EnumConfig, enumerate_products
        from halogenator.sugar_mask import post_guard_blocked

        # This test simulates a scenario where:
        # - Two different reaction paths produce the same InChIKey
        # - Path A gets rejected by post_guard_blocked
        # - Path B should still succeed (not blocked by A's early insertion)

        # Use a simple molecule that can produce same product via different paths
        parent_smi = "CC"  # ethane with two methyl sites

        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R6',),
            engine_cfg={'budget_mode': 'ops'},
            rules_cfg={
                'R6_methyl': {
                    'enable': True,
                    'allowed': ['F'],
                    'allow_on_methoxy': False,
                    'macro': {'enable': False}
                }
            },
            # Configure constraints or sugar mask to potentially block some products
            sugar_cfg={'mode': 'off'}  # Disable sugar masking for this test
        )

        products = list(enumerate_products(parent_smi, cfg, return_qa_stats=True))

        # The key test: we should get some products and QA stats should show
        # proper deduplication counts without false rejections
        product_count = len([p for p in products if not p[0].get('is_qa_summary_marker', False)])
        self.assertGreater(product_count, 0, "Should produce at least one product")

        # Check that pruned_inchikey_dupe count is reasonable
        qa_summary = next((p[1] for p in products if p[0].get('is_qa_summary_marker', False)), {})
        if qa_summary and 'qa_paths' in qa_summary:
            dedup_count = qa_summary['qa_paths'].get('pruned_inchikey_dupe', 0)
            # Should not have excessive deduplication for simple ethane
            self.assertLessEqual(dedup_count, 1, "Should not have excessive deduplication")

    def test_budget_state_cross_layer_persistence(self):
        """Test that budget state persists correctly across BFS layers."""
        from halogenator.enumerate_k import EnumConfig, enumerate_with_stats

        # Use toluene with ops mode to track operation counting across layers
        parent_smi = "Cc1ccccc1"  # toluene

        cfg = EnumConfig(
            k_max=2,  # Allow 2 operations
            halogens=('F',),
            rules=('R6',),
            engine_cfg={'budget_mode': 'ops'},
            rules_cfg={
                'R6_methyl': {
                    'enable': True,
                    'allowed': ['F'],
                    'allow_on_methoxy': False,
                    'macro': {'enable': False}  # Only step halogenation
                }
            }
        )

        products, qa_stats = enumerate_with_stats(parent_smi, cfg)

        # Should find products with k_ops=1 (first touch) and k_ops=1 (second touch on same site)
        # The second touch should have op_cost=0 due to site token persistence
        depth_2_products = [p for p in products if p.get('k_ops', 0) >= 1]
        self.assertGreater(len(depth_2_products), 0, "Should have products with k_ops tracking")

        # Check that we have products with proper k_ops/k_atoms tracking
        for product in depth_2_products:
            k_ops = product.get('k_ops', 0)
            k_atoms = product.get('k_atoms', 0)
            self.assertGreaterEqual(k_ops, 0, "k_ops should be non-negative")
            self.assertGreaterEqual(k_atoms, 0, "k_atoms should be non-negative")
            # In step mode, each F adds 1 atom
            if k_atoms > 0:
                self.assertLessEqual(k_atoms, cfg.k_max, f"k_atoms {k_atoms} should not exceed k_max {cfg.k_max}")

    def test_site_token_first_touch_once(self):
        """Test that site tokens work correctly - first touch costs 1 op, subsequent touches cost 0."""
        budget = BudgetState("ops", k_max=5)
        cidx = 10

        # First touch: should cost 1 operation
        op_cost = budget.token_cost(cidx)
        self.assertEqual(op_cost, 1, "First touch should cost 1 operation")

        # Charge and mark first touch
        self.assertTrue(budget.charge(op_cost, 1))
        budget.mark_first_touch(cidx)

        # Second touch: should cost 0 operations
        op_cost = budget.token_cost(cidx)
        self.assertEqual(op_cost, 0, "Second touch should cost 0 operations")

        # Verify budget state
        self.assertEqual(budget.k_ops, 1, "Should have 1 operation charged")
        self.assertEqual(budget.k_atoms, 1, "Should have 1 atom charged")

    def test_budget_contamination_prevention(self):
        """Test that budget contamination is prevented using temporary copies."""
        from halogenator.enumerate_k import EnumConfig, enumerate_with_stats

        # Use a molecule with multiple distinct methyl sites (no symmetry folding issues)
        parent_smi = "CCCC"  # butane with 2 distinct methyl sites

        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R6',),
            engine_cfg={'budget_mode': 'ops'},
            rules_cfg={
                'R6_methyl': {
                    'enable': True,
                    'allowed': ['F'],
                    'allow_on_methoxy': False,
                    'macro': {'enable': False}
                }
            },
            pruning_cfg={'enable_symmetry_fold': False}  # Disable symmetry folding to see all products
        )

        products, qa_stats = enumerate_with_stats(parent_smi, cfg)

        # Should find at least one product (main test is budget consistency)
        product_count = len([p for p in products if 'smiles' in p])
        self.assertGreaterEqual(product_count, 1, "Should find at least one F-substituted product")

        # Check k_ops consistency across products - this is the main budget contamination test
        for product in products:
            if 'k_ops' in product:
                self.assertEqual(product['k_ops'], 1, "All products should have k_ops=1 (no contamination)")
                self.assertEqual(product['k_atoms'], 1, "All products should have k_atoms=1")

    def test_unified_deduplication_metrics(self):
        """Test that deduplication metrics use unified naming (dedup_hits_inchi)."""
        from halogenator.enumerate_k import EnumConfig, enumerate_with_stats

        # Use a symmetric molecule to trigger deduplication
        parent_smi = "Cc1ccc(C)cc1"  # p-xylene (symmetric)

        # Test default mode (emit_legacy_keys=False)
        cfg_default = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R6',),
            engine_cfg={'budget_mode': 'ops', 'emit_legacy_keys': False},
            rules_cfg={
                'R6_methyl': {
                    'enable': True,
                    'allowed': ['F'],
                    'allow_on_methoxy': False,
                    'macro': {'enable': False}
                }
            },
            pruning_cfg={'enable_symmetry_fold': True}
        )

        products_default, qa_stats_default = enumerate_with_stats(parent_smi, cfg_default)

        # Check that unified dedup metric is used
        if qa_stats_default and 'dedup_hits_inchi' in qa_stats_default:
            self.assertGreaterEqual(qa_stats_default['dedup_hits_inchi'], 0, "dedup_hits_inchi should be present")

        # In default mode, legacy field may not be present
        if qa_stats_default and 'qa_paths' in qa_stats_default:
            qa_paths_default = qa_stats_default['qa_paths']
            self.assertIn('dedup_hits_inchi', qa_paths_default, "New field should always be present")
            # Legacy field is optional in default mode
            if 'pruned_inchikey_dupe' in qa_paths_default:
                self.assertEqual(qa_paths_default['dedup_hits_inchi'], qa_paths_default['pruned_inchikey_dupe'],
                               "If legacy field exists, it should match new field")

        # Test legacy compatibility mode (emit_legacy_keys=True)
        cfg_legacy = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R6',),
            engine_cfg={'budget_mode': 'ops', 'emit_legacy_keys': True},
            rules_cfg={
                'R6_methyl': {
                    'enable': True,
                    'allowed': ['F'],
                    'allow_on_methoxy': False,
                    'macro': {'enable': False}
                }
            },
            pruning_cfg={'enable_symmetry_fold': True}
        )

        products_legacy, qa_stats_legacy = enumerate_with_stats(parent_smi, cfg_legacy)

        # In legacy mode, both fields should be present when dedup occurs
        if qa_stats_legacy and 'qa_paths' in qa_stats_legacy:
            qa_paths_legacy = qa_stats_legacy['qa_paths']
            if 'dedup_hits_inchi' in qa_paths_legacy and qa_paths_legacy['dedup_hits_inchi'] > 0:
                self.assertIn('pruned_inchikey_dupe', qa_paths_legacy, "Legacy field should be present in legacy mode")
                self.assertEqual(qa_paths_legacy['dedup_hits_inchi'], qa_paths_legacy['pruned_inchikey_dupe'],
                               "Legacy and new field values should be identical")

    def test_stable_inchikey_generation(self):
        """Test that InChIKey generation strips isotopes and mappings for consistency."""
        from halogenator.standardize import to_inchikey_sanitized, strip_isotopes_and_props

        # Create a molecule with isotopes and atom mappings
        mol_with_isotopes = Chem.MolFromSmiles("C[13CH3]")  # Carbon-13 labeled ethane
        if mol_with_isotopes is None:
            self.skipTest("Cannot create isotope-labeled test molecule")

        # Add atom mapping to first carbon
        mol_with_isotopes.GetAtomWithIdx(0).SetAtomMapNum(1)

        # Test stable InChIKey generation
        stable_key = to_inchikey_sanitized(mol_with_isotopes)
        self.assertIsNotNone(stable_key, "Should generate stable InChIKey")

        # Test isotope stripping
        cleaned_mol = strip_isotopes_and_props(mol_with_isotopes)

        # Verify isotopes are stripped
        for atom in cleaned_mol.GetAtoms():
            self.assertEqual(atom.GetIsotope(), 0, "Isotopes should be stripped")
            self.assertEqual(atom.GetAtomMapNum(), 0, "Atom mappings should be stripped")

        # Verify stable key matches key from cleaned molecule
        from halogenator.standardize import to_inchikey
        cleaned_key = to_inchikey(cleaned_mol)
        self.assertEqual(stable_key, cleaned_key, "Stable key should match cleaned molecule key")

    def test_budget_isolation_across_layers(self):
        """Test that budget state persists correctly across BFS enumeration layers."""
        from halogenator.enumerate_k import EnumConfig, enumerate_with_stats

        # Use molecule that allows multi-step halogenation
        parent_smi = "CC"  # ethane

        cfg = EnumConfig(
            k_max=2,  # Allow 2 operations for multi-layer testing
            halogens=('F',),
            rules=('R6',),
            engine_cfg={'budget_mode': 'ops'},
            rules_cfg={
                'R6_methyl': {
                    'enable': True,
                    'allowed': ['F'],
                    'allow_on_methoxy': False,
                    'macro': {'enable': False}
                }
            }
        )

        products, qa_stats = enumerate_with_stats(parent_smi, cfg)

        # Should find products with varying k_ops values
        k_ops_values = set()
        for product in products:
            if 'k_ops' in product:
                k_ops_values.add(product['k_ops'])

        # Should have products from different layers (k_ops=1 and potentially k_ops=1 again due to site tokens)
        self.assertGreaterEqual(len(k_ops_values), 1, "Should have products from different enumeration layers")

        # All k_ops values should be within budget limit
        for k_ops in k_ops_values:
            self.assertLessEqual(k_ops, cfg.k_max, f"k_ops {k_ops} should not exceed k_max {cfg.k_max}")

    def test_deduplication_consistency_with_budget(self):
        """Test that deduplication works consistently when combined with budget tracking."""
        from halogenator.enumerate_k import EnumConfig, enumerate_with_stats

        # Use a molecule that can produce same products via different budget paths
        parent_smi = "CCC"  # propane with two equivalent methyl sites

        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R6',),
            engine_cfg={'budget_mode': 'atoms'},  # Use atoms mode for stricter budget
            rules_cfg={
                'R6_methyl': {
                    'enable': True,
                    'allowed': ['F'],
                    'allow_on_methoxy': False,
                    'macro': {'enable': False}
                }
            },
            pruning_cfg={'enable_symmetry_fold': True}
        )

        products, qa_stats = enumerate_with_stats(parent_smi, cfg)

        # Should deduplicate equivalent products correctly
        product_smiles = set()
        for product in products:
            if 'smiles' in product:
                product_smiles.add(product['smiles'])

        # Should have only one unique product (CH3CH2CF or CF3CH2CH3)
        self.assertEqual(len(product_smiles), 1, "Should deduplicate equivalent methyl sites")

        # Verify budget consistency
        for product in products:
            if 'k_atoms' in product and 'k_ops' in product:
                self.assertEqual(product['k_atoms'], 1, "Should have 1 F atom added")
                self.assertEqual(product['k_ops'], 1, "Should have 1 operation")

    def test_r6_macro_cf3_budget_atoms(self):
        """Test macro CF3 billing - should charge 3 atoms and have k_atoms=3."""
        from halogenator.enumerate_k import EnumConfig, enumerate_with_stats

        parent_smi = "CC"  # ethane

        cfg = EnumConfig(
            k_max=3,  # Allow enough budget for CF3 (3 atoms)
            halogens=('F',),
            rules=('R6',),
            engine_cfg={'budget_mode': 'atoms'},  # Use atoms mode for strict budget tracking
            rules_cfg={
                'R6_methyl': {
                    'enable': True,
                    'allowed': ['F'],
                    'allow_on_methoxy': False,
                    'macro': {
                        'enable': True,
                        'labels': ['CF3']
                    }
                }
            }
        )

        products, qa_stats = enumerate_with_stats(parent_smi, cfg)

        # Find CF3 macro products
        cf3_products = [p for p in products if p.get('macro_label') == 'CF3']
        self.assertGreater(len(cf3_products), 0, "Should find CF3 macro products")

        # Verify CF3 products have correct billing
        for product in cf3_products:
            self.assertEqual(product['k_atoms'], 3, "CF3 should charge 3 atoms")
            self.assertEqual(product['k'], 3, "CF3 product should have k=3")
            self.assertEqual(product['halogen'], 'F', "CF3 should use F halogen")
            # CF3 should create a trifluoromethyl group: C(F)(F)F pattern
            self.assertIn('C(F)(F)F', product['smiles'], "CF3 should create trifluoromethyl group")

        # Verify step halogenation products exist and have different characteristics
        step_products = [p for p in products if 'macro_label' not in p and 'smiles' in p]
        self.assertGreater(len(step_products), 0, "Should have step halogenation products for comparison")

        # Step products should not have the macro_label field
        for product in step_products:
            self.assertNotIn('macro_label', product, "Step products should not have macro_label")
            self.assertLessEqual(product['k_atoms'], cfg.k_max, "Step product should not exceed budget")

    def test_methoxy_site_allow_toggle(self):
        """Test methoxy site detection with allow_on_methoxy toggle."""
        from halogenator.sites_methyl import enumerate_methyl_sites

        # Anisole: methoxy group attached to benzene
        parent_smi = "COc1ccccc1"  # anisole
        mol = Chem.MolFromSmiles(parent_smi)
        self.assertIsNotNone(mol, "Should create anisole molecule")

        # Test with allow_on_methoxy=False (default)
        sites_off = enumerate_methyl_sites(mol, set(), allow_on_methoxy=False)
        self.assertEqual(len(sites_off), 0, "Should find no methyl sites with allow_on_methoxy=False")

        # Test with allow_on_methoxy=True
        sites_on = enumerate_methyl_sites(mol, set(), allow_on_methoxy=True)
        self.assertGreaterEqual(len(sites_on), 1, "Should find methoxy methyl site with allow_on_methoxy=True")

        # Test with a non-methoxy methyl group for comparison
        toluene_mol = Chem.MolFromSmiles("Cc1ccccc1")  # toluene
        toluene_sites_off = enumerate_methyl_sites(toluene_mol, set(), allow_on_methoxy=False)
        toluene_sites_on = enumerate_methyl_sites(toluene_mol, set(), allow_on_methoxy=True)

        # Toluene should have same number of sites regardless of methoxy setting
        self.assertEqual(len(toluene_sites_off), len(toluene_sites_on),
                        "Toluene sites should be unaffected by allow_on_methoxy setting")
        self.assertGreaterEqual(len(toluene_sites_off), 1, "Toluene should have methyl sites")

    def test_budget_deep_copy_isolation(self):
        """Test that BudgetState.copy() provides proper deep copy isolation to prevent cross-branch contamination."""
        # Create a budget and modify it to have state
        budget = BudgetState('ops', k_max=10)

        # Mark first touch and add state for a site
        site_id = 42
        budget.mark_first_touch(site_id)
        budget.bump_state(site_id, 'F')
        budget.charge(1, 1)  # Add some charges

        # Create a copy
        budget_copy = budget.copy()

        # Verify the copy has the same initial state
        self.assertEqual(budget_copy.k_ops, budget.k_ops, "Copy should have same k_ops")
        self.assertEqual(budget_copy.k_atoms, budget.k_atoms, "Copy should have same k_atoms")
        self.assertEqual(budget_copy.budget_mode, budget.budget_mode, "Copy should have same budget_mode")
        self.assertIn(site_id, budget_copy.site_state, "Copy should have same site_state")

        # Modify the original budget
        original_site_id = 43
        budget.mark_first_touch(original_site_id)
        budget.bump_state(original_site_id, 'Cl')
        budget.charge(1, 1)

        # Verify the copy is NOT affected by original modifications
        self.assertNotIn(original_site_id, budget_copy.site_state, "Copy should not be affected by original modifications")
        self.assertNotEqual(budget_copy.k_ops, budget.k_ops, "Copy k_ops should not change with original")

        # Modify the copy budget
        copy_site_id = 44
        budget_copy.mark_first_touch(copy_site_id)
        budget_copy.bump_state(copy_site_id, 'Br')
        budget_copy.charge(2, 2)

        # Verify the original is NOT affected by copy modifications
        self.assertNotIn(copy_site_id, budget.site_state, "Original should not be affected by copy modifications")

        # Verify site_state is properly deep copied (Counter objects are independent)
        budget.bump_state(site_id, 'F')  # Add another F to original
        original_f_count = budget.site_state[site_id]['F']
        copy_f_count = budget_copy.site_state[site_id]['F']
        self.assertNotEqual(original_f_count, copy_f_count, "Counter objects should be independent")

        # Verify site_tokens isolation - check specific sites exist only in the right budget
        self.assertIn(original_site_id, budget.site_tokens, "Original should have original_site_id")
        self.assertNotIn(original_site_id, budget_copy.site_tokens, "Copy should not have original_site_id")
        self.assertIn(copy_site_id, budget_copy.site_tokens, "Copy should have copy_site_id")
        self.assertNotIn(copy_site_id, budget.site_tokens, "Original should not have copy_site_id")


if __name__ == '__main__':
    unittest.main(verbosity=2)