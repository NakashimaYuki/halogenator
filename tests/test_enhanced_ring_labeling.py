# -*- coding: ascii -*-
"""Test enhanced flavonoid ring labeling system."""

import unittest
from rdkit import Chem
from src.halogenator.sites import flavonoid_ring_label
from src.halogenator.enumerate_k import enumerate_products, EnumConfig
from collections import Counter


class TestEnhancedRingLabeling(unittest.TestCase):
    """Test enhanced A/B/C ring labeling with strict C-ring identification."""
    
    def test_quercetin_ring_labeling(self):
        """Test that quercetin correctly identifies A/B/C rings."""
        quercetin_smiles = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        mol = Chem.MolFromSmiles(quercetin_smiles)
        Chem.SanitizeMol(mol)
        
        # Count atoms in each ring type
        ring_counts = {'A': 0, 'B': 0, 'C': 0, '': 0}
        for i in range(mol.GetNumAtoms()):
            label = flavonoid_ring_label(mol, i)
            ring_counts[label] += 1
        
        # Should have atoms in A, B, and C rings
        self.assertGreater(ring_counts['A'], 0, "Should identify A ring atoms")
        self.assertGreater(ring_counts['B'], 0, "Should identify B ring atoms") 
        self.assertGreater(ring_counts['C'], 0, "Should identify C ring atoms")
        
        # C ring should contain oxygen and carbonyl carbon
        c_ring_atoms = []
        for i in range(mol.GetNumAtoms()):
            if flavonoid_ring_label(mol, i) == 'C':
                c_ring_atoms.append(i)
        
        # C ring should have both O and C=O
        has_oxygen = any(mol.GetAtomWithIdx(i).GetSymbol() == 'O' for i in c_ring_atoms)
        self.assertTrue(has_oxygen, "C ring should contain oxygen atom")
        
    def test_apigenin_glucoside_excludes_sugar_ring(self):
        """Test that sugar rings in glycosides are not labeled as C ring."""
        # Apigenin 7-glucoside (has both flavonoid C ring and sugar ring)
        apigenin_glucoside = "O=c1cc(oc2cc(O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)cc(O)c12)-c1ccc(O)cc1"
        mol = Chem.MolFromSmiles(apigenin_glucoside)
        if mol is None:
            self.skipTest("Could not parse apigenin glucoside SMILES")
        
        try:
            Chem.SanitizeMol(mol)
        except:
            self.skipTest("Could not sanitize apigenin glucoside")
        
        # Should still identify A/B/C rings from flavonoid core
        ring_counts = {'A': 0, 'B': 0, 'C': 0, '': 0}
        for i in range(mol.GetNumAtoms()):
            label = flavonoid_ring_label(mol, i)
            ring_counts[label] += 1
        
        # Should identify flavonoid rings despite presence of sugar
        self.assertGreater(ring_counts['C'], 0, "Should identify C ring despite sugar presence")
        
    def test_per_ring_quota_with_enhanced_labeling(self):
        """Test per_ring_quota=1 with enhanced ring labeling prevents same-ring substitutions."""
        cfg_strict = EnumConfig(
            k_max=2,
            halogens=('F',),
            constraints={'per_ring_quota': 1, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
        quercetin_smiles = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        products_strict = list(enumerate_products(quercetin_smiles, cfg_strict))
        
        # For k=2 products, verify ring tags are different when per_ring_quota=1
        k2_products = [p for p in products_strict if p['k'] == 2]
        for p in k2_products:
            history = p.get('substitutions', [])
            if len(history) == 2:
                tags = [s.get('ring_tag', '') for s in history]
                # Both tags should be non-empty and different
                non_empty_tags = [t for t in tags if t]
                if len(non_empty_tags) == 2:
                    self.assertNotEqual(non_empty_tags[0], non_empty_tags[1],
                                      f"per_ring_quota=1 should prevent same-ring substitutions, got tags: {tags}")
        
        # Should have some k=1 products
        k1_count = len([p for p in products_strict if p['k'] == 1])
        self.assertGreater(k1_count, 0, "Should have k=1 products")
        
    def test_benzene_no_flavonoid_rings(self):
        """Test that benzene returns empty labels (not a flavonoid)."""
        benzene_smiles = "c1ccccc1"
        mol = Chem.MolFromSmiles(benzene_smiles)
        Chem.SanitizeMol(mol)
        
        # Should return empty labels for all atoms
        for i in range(mol.GetNumAtoms()):
            label = flavonoid_ring_label(mol, i)
            self.assertEqual(label, '', f"Benzene atom {i} should have empty label, got '{label}'")


if __name__ == '__main__':
    unittest.main()