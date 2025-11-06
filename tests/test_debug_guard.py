# -*- coding: ascii -*-
"""Debug guard functionality."""

import unittest
import os
from unittest.mock import patch

from src.halogenator.enumerate_k1 import enumerate_k1_with_stats
from src.halogenator.enumerate_k import enumerate_with_stats, EnumConfig


class TestDebugGuard(unittest.TestCase):
    """Debug guard functionality to understand current behavior."""
    
    def setUp(self):
        """Set up guard environment."""
        os.environ['HALO_RDKIT_GUARD'] = '1'
    
    def tearDown(self):
        """Clean up guard environment."""
        if 'HALO_RDKIT_GUARD' in os.environ:
            del os.environ['HALO_RDKIT_GUARD']
    
    def test_guard_enabled_check(self):
        """Test that guard is enabled."""
        from src.halogenator.guard import rdkit_guard_enabled
        self.assertTrue(rdkit_guard_enabled(), "Guard should be enabled")
    
    def test_k1_legacy_with_chem_error(self):
        """Test k=1 legacy with Chem error."""
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1',),
        )
        
        # Mock Chem.MolFromSmiles to raise an exception
        with patch('src.halogenator.chem_compat.Chem.MolFromSmiles') as mock_mol_from_smiles:
            mock_mol_from_smiles.side_effect = Exception("Mock RDKit error")
            
            products, qa_stats = enumerate_k1_with_stats("c1ccccc1", cfg, stream_shape='legacy')
            
            print(f"K=1 Legacy with error - Products: {len(products)}")
            print(f"K=1 Legacy with error - QA stats: {qa_stats}")
            print(f"K=1 Legacy with error - rdkit_error: {qa_stats.get('qa_paths', {}).get('rdkit_error', 0)}")
            print(f"K=1 Legacy with error - template_unsupported: {qa_stats.get('template_unsupported', 0)}")
            
            # Check that guard events are present
            qa_paths = qa_stats.get('qa_paths', {})
            self.assertGreater(qa_paths.get('rdkit_error', 0), 0, "Guard rdkit_error should be > 0")
            self.assertGreater(qa_stats.get('template_unsupported', 0), 0, "template_unsupported should be > 0")


if __name__ == '__main__':
    unittest.main()