# -*- coding: ascii -*-
"""Test that RDKit guard events surface in final QA statistics."""

import unittest
import os
from unittest.mock import patch, MagicMock

from src.halogenator.enumerate_k1 import enumerate_k1_with_stats
from src.halogenator.enumerate_k import enumerate_with_stats, EnumConfig


class TestGuardEventsSurface(unittest.TestCase):
    """Test that guard events appear in final QA stats."""
    
    def setUp(self):
        """Set up test environment."""
        # Enable guard for testing
        os.environ['HALO_RDKIT_GUARD'] = '1'
    
    def tearDown(self):
        """Clean up test environment."""
        if 'HALO_RDKIT_GUARD' in os.environ:
            del os.environ['HALO_RDKIT_GUARD']
    
    def test_guard_events_surface_in_k1_legacy(self):
        """Test guard events surface in k=1 legacy format."""
        from src.halogenator.enumerate_k import EnumConfig
        
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1',),
            constraints={'per_ring_quota': 2},
            std_cfg={'do_tautomer': False},
            qc_cfg={'sanitize_strict': False}
        )
        
        # Mock RDKit to raise exceptions that guard should catch
        with patch('src.halogenator.chem_compat.Chem.MolFromSmiles') as mock_mol_from_smiles:
            # Make MolFromSmiles raise an exception to trigger guard
            mock_mol_from_smiles.side_effect = Exception("Mock RDKit error")
            
            products, qa_stats = enumerate_k1_with_stats("CCO", cfg, stream_shape='legacy')
            
            # Should have version:'1' for legacy
            self.assertEqual(qa_stats.get('version'), '1')
            
            # Guard events should surface in final qa_paths
            qa_paths = qa_stats.get('qa_paths', {})
            self.assertGreater(qa_paths.get('rdkit_error', 0), 0, "Guard rdkit_error should be > 0")
            self.assertGreater(qa_stats.get('template_unsupported', 0), 0, "template_unsupported should be > 0")
    
    def test_guard_events_surface_in_k1_v2(self):
        """Test guard events surface in k=1 v2 format."""
        from src.halogenator.enumerate_k import EnumConfig
        
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1',),
            constraints={'per_ring_quota': 2},
            std_cfg={'do_tautomer': False},
            qc_cfg={'sanitize_strict': False}
        )
        
        # Mock RDKit to raise exceptions that guard should catch
        with patch('src.halogenator.chem_compat.Chem.MolFromSmiles') as mock_mol_from_smiles:
            # Make MolFromSmiles raise an exception to trigger guard
            mock_mol_from_smiles.side_effect = Exception("Mock RDKit error")
            
            products, qa_stats = enumerate_k1_with_stats("CCO", cfg, stream_shape='v2')
            
            # Should have version:'2' for v2
            self.assertEqual(qa_stats.get('version'), '2')
            
            # Guard events should surface in final qa_paths
            qa_paths = qa_stats.get('qa_paths', {})
            self.assertGreater(qa_paths.get('rdkit_error', 0), 0, "Guard rdkit_error should be > 0")
            self.assertGreater(qa_stats.get('template_unsupported', 0), 0, "template_unsupported should be > 0")
    
    def test_guard_events_surface_in_k_gt_1(self):
        """Test guard events surface in k>1 enumeration."""
        cfg = EnumConfig(
            k_max=2,
            halogens=('F',),
            rules=('R1',),
            constraints={'per_ring_quota': 2},
            std_cfg={'do_tautomer': False},
            qc_cfg={'sanitize_strict': False}
        )
        
        # Mock RDKit to raise exceptions that guard should catch
        with patch('src.halogenator.enumerate_k.std_from_smiles') as mock_std:
            # Make standardization raise an exception to trigger guard
            mock_std.side_effect = Exception("Mock RDKit error")
            
            products, qa_stats = enumerate_with_stats("CCO", cfg)
            
            # Should have version:'2' (k>1 always uses v2)
            self.assertEqual(qa_stats.get('version'), '2')
            
            # Guard events should surface in final qa_paths
            qa_paths = qa_stats.get('qa_paths', {})
            self.assertGreater(qa_paths.get('rdkit_error', 0), 0, "Guard rdkit_error should be > 0")


if __name__ == '__main__':
    unittest.main()