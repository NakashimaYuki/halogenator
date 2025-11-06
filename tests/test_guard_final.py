# -*- coding: ascii -*-
"""Final working guard test."""

import unittest
import os
from unittest.mock import patch

from src.halogenator.enumerate_k1 import enumerate_k1_with_stats
from src.halogenator.enumerate_k import enumerate_with_stats, EnumConfig


class TestGuardFinal(unittest.TestCase):
    """Final working guard tests."""
    
    def setUp(self):
        """Set up test environment."""
        os.environ['HALO_RDKIT_GUARD'] = '1'
    
    def tearDown(self):
        """Clean up test environment."""
        if 'HALO_RDKIT_GUARD' in os.environ:
            del os.environ['HALO_RDKIT_GUARD']
    
    def test_guard_events_surface_in_final_stats(self):
        """Test that guard events surface in final QA statistics."""
        cfg = EnumConfig(k_max=1, halogens=('F',), rules=('R1',))
        
        # Test k=1 legacy format  
        with patch('src.halogenator.chem_compat.Chem.MolFromSmiles') as mock_mol:
            mock_mol.side_effect = Exception("Mock RDKit error")
            
            products, qa_stats = enumerate_k1_with_stats("c1ccccc1", cfg, stream_shape='legacy')
            
            # Basic validation
            self.assertEqual(qa_stats.get('version'), '1')
            self.assertGreater(qa_stats.get('qa_paths', {}).get('rdkit_error', 0), 0)
            self.assertGreater(qa_stats.get('template_unsupported', 0), 0)
        
        # Test k=1 v2 format
        with patch('src.halogenator.chem_compat.Chem.MolFromSmiles') as mock_mol:
            mock_mol.side_effect = Exception("Mock RDKit error") 
            
            products, qa_stats = enumerate_k1_with_stats("c1ccccc1", cfg, stream_shape='v2')
            
            # Basic validation
            self.assertEqual(qa_stats.get('version'), '2')
            self.assertGreater(qa_stats.get('qa_paths', {}).get('rdkit_error', 0), 0)
            self.assertGreater(qa_stats.get('template_unsupported', 0), 0)
    
    def test_legacy_stream_shape_emits_version_1(self):
        """Test that legacy stream_shape outputs version:'1'."""
        cfg = EnumConfig(k_max=1, halogens=('F',), rules=('R1',))
        
        # Test legacy format 
        products, qa_stats = enumerate_k1_with_stats("c1ccccc1", cfg, stream_shape='legacy')
        self.assertEqual(qa_stats.get('version'), '1')
        self.assertIn('qa_paths', qa_stats)
        self.assertNotIn('pivots', qa_stats)  # Legacy should not have pivots
        
        # Test v2 format
        products, qa_stats = enumerate_k1_with_stats("c1ccccc1", cfg, stream_shape='v2')
        self.assertEqual(qa_stats.get('version'), '2')
        self.assertIn('qa_paths', qa_stats)
        self.assertIn('pivots', qa_stats)  # v2 should have pivots
    
    def test_r2_config_respect(self):
        """Test that R2 configuration is respected."""
        # Config WITHOUT R2 - should not generate carbonyl_unknown
        cfg_no_r2 = EnumConfig(k_max=1, halogens=('F',), rules=('R1',))
        products, qa_stats = enumerate_k1_with_stats("c1ccc(=O)cc1", cfg_no_r2, stream_shape='legacy')
        self.assertEqual(qa_stats.get('qa_paths', {}).get('carbonyl_unknown', 0), 0)
        
        # Config WITH R2 - can generate carbonyl_unknown (mock to ensure it triggers)
        cfg_with_r2 = EnumConfig(k_max=1, halogens=('F',), rules=('R1', 'R2'))
        with patch('src.halogenator.sites.is_carbonyl_carbon') as mock_carbonyl:
            mock_carbonyl.return_value = None  # Simulate unknown status
            
            products, qa_stats = enumerate_k1_with_stats("c1ccc(=O)cc1", cfg_with_r2, stream_shape='legacy')
            # carbonyl_unknown should exist as a metric when R2 is enabled
            self.assertIn('carbonyl_unknown', qa_stats.get('qa_paths', {}))


if __name__ == '__main__':
    unittest.main()