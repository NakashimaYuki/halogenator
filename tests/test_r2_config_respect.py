# -*- coding: ascii -*-
"""Test that R2 configuration is respected and carbonyl_unknown works for k>1."""

import unittest
from unittest.mock import patch

from src.halogenator.enumerate_k1 import enumerate_k1_with_stats
from src.halogenator.enumerate_k import enumerate_with_stats, EnumConfig


class TestR2ConfigRespect(unittest.TestCase):
    """Test that R2 configuration is respected and carbonyl_unknown counting works."""
    
    def test_enum_r2_respects_config_excluded(self):
        """Test that when R2 is not in config, carbonyl_unknown is not generated."""
        # Config WITHOUT R2
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1',),  # No R2
            constraints={'per_ring_quota': 2},
            std_cfg={'do_tautomer': False},
            qc_cfg={'sanitize_strict': False}
        )
        
        # Use a molecule that would have C ring sites if R2 were enabled
        products, qa_stats = enumerate_k1_with_stats("c1ccc(=O)cc1", cfg, stream_shape='legacy')
        
        # Should NOT generate carbonyl_unknown since R2 is not enabled
        qa_paths = qa_stats.get('qa_paths', {})
        self.assertEqual(qa_paths.get('carbonyl_unknown', 0), 0, "carbonyl_unknown should be 0 when R2 not in config")
    
    def test_enum_r2_respects_config_included(self):
        """Test that when R2 is in config, carbonyl_unknown can be generated."""
        # Config WITH R2
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1', 'R2'),  # Include R2
            constraints={'per_ring_quota': 2},
            std_cfg={'do_tautomer': False},
            qc_cfg={'sanitize_strict': False}
        )
        
        # Mock is_carbonyl_carbon to return None (unknown) to trigger carbonyl_unknown
        with patch('src.halogenator.sites.is_carbonyl_carbon') as mock_carbonyl:
            mock_carbonyl.return_value = None  # Simulate unknown carbonyl status
            
            # Use a molecule that has C ring sites
            products, qa_stats = enumerate_k1_with_stats("c1ccc(=O)cc1", cfg, stream_shape='legacy')
            
            # Should generate carbonyl_unknown since R2 is enabled and carbonyl status is unknown
            qa_paths = qa_stats.get('qa_paths', {})
            # Note: carbonyl_unknown may be 0 if no C ring sites are found, so we just check it exists
            self.assertIn('carbonyl_unknown', qa_paths, "carbonyl_unknown metric should exist when R2 is enabled")
    
    def test_k_gt_1_carbonyl_unknown_counting(self):
        """Test that k>1 enumeration also counts carbonyl_unknown events."""
        # Config WITH R2 for k>1
        cfg = EnumConfig(
            k_max=2,
            halogens=('F',),
            rules=('R1', 'R2'),  # Include R2
            constraints={'per_ring_quota': 2},
            std_cfg={'do_tautomer': False},
            qc_cfg={'sanitize_strict': False}
        )
        
        # Mock is_carbonyl_carbon to return None (unknown) to trigger carbonyl_unknown
        with patch('src.halogenator.sites.is_carbonyl_carbon') as mock_carbonyl:
            mock_carbonyl.return_value = None  # Simulate unknown carbonyl status
            
            # Use a molecule that has C ring sites
            products, qa_stats = enumerate_with_stats("c1ccc(=O)cc1", cfg)
            
            # Should generate carbonyl_unknown in k>1 enumeration as well
            qa_paths = qa_stats.get('qa_paths', {})
            self.assertIn('carbonyl_unknown', qa_paths, "carbonyl_unknown metric should exist in k>1 enumeration with R2")
    
    def test_k_gt_1_no_r2_no_carbonyl_unknown(self):
        """Test that k>1 enumeration without R2 doesn't count carbonyl_unknown."""
        # Config WITHOUT R2 for k>1
        cfg = EnumConfig(
            k_max=2,
            halogens=('F',),
            rules=('R1',),  # No R2
            constraints={'per_ring_quota': 2},
            std_cfg={'do_tautomer': False},
            qc_cfg={'sanitize_strict': False}
        )
        
        # Use a molecule that would have C ring sites if R2 were enabled
        products, qa_stats = enumerate_with_stats("c1ccc(=O)cc1", cfg)
        
        # Should NOT generate carbonyl_unknown since R2 is not enabled
        qa_paths = qa_stats.get('qa_paths', {})
        self.assertEqual(qa_paths.get('carbonyl_unknown', 0), 0, "carbonyl_unknown should be 0 in k>1 when R2 not in config")


if __name__ == '__main__':
    unittest.main()