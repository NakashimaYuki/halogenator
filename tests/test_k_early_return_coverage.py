# -*- coding: ascii -*-
"""Test coverage for k>1 early return branches."""

import os
import unittest
from unittest.mock import patch, MagicMock
from src.halogenator.enumerate_k import enumerate_with_stats, EnumConfig
from tests import CleanLogsTestCase


class TestKEarlyReturnCoverage(CleanLogsTestCase):
    """Test k>1 early return branches."""
    
    def test_early_return_standardization_failure(self):
        """Test early return when std_from_smiles fails."""
        cfg = EnumConfig(k_max=2, rules=('R1', 'R2'))
        
        # Enable RDKit guard to catch exceptions
        with patch.dict(os.environ, {'HALO_RDKIT_GUARD': '1'}):
            # Mock std_from_smiles to raise an exception (simulate RDKit error)
            with patch('src.halogenator.enumerate_k.std_from_smiles') as mock_std:
                mock_std.side_effect = Exception("RDKit standardization error")
                
                # Call enumerate_with_stats to trigger early return handling
                products, qa_stats = enumerate_with_stats("CCO", cfg)
            
            # Verify early return behavior
            self.assertEqual(len(products), 0, "Should produce no products on early return")
            
            # Verify exactly 1 attempt recorded
            self.assertEqual(qa_stats.get('attempts', 0), 1, "Should record exactly 1 attempt")
            
            # Verify template_unsupported count (guard failures recorded here)
            self.assertEqual(qa_stats.get('template_unsupported', 0), 1, "Should record 1 template_unsupported")
            
            # Verify no_product_matches is 0 (guard failures don't count as no_product_matches)
            self.assertEqual(qa_stats.get('no_product_matches', 0), 0, "Should have 0 no_product_matches")
            
            # Verify rdkit_error was recorded in qa_paths
            qa_paths = qa_stats.get('qa_paths', {})
            self.assertGreater(qa_paths.get('rdkit_error', 0), 0, "Should record rdkit_error > 0")
            
            # Verify version is '2' (k>1 only supports v2)
            self.assertEqual(qa_stats.get('version'), '2', "Early return should produce version '2'")
            
            # Verify pivots structure exists (empty but present)
            self.assertIn('pivots', qa_stats, "Should contain pivots structure")
            pivots = qa_stats['pivots']
            self.assertIsInstance(pivots, dict, "Pivots should be a dictionary")
    
    def test_early_return_parent_key_failure(self):
        """Test early return when parent key extraction fails."""
        cfg = EnumConfig(k_max=2, rules=('R1', 'R2'))
        
        # Enable RDKit guard to catch exceptions
        with patch.dict(os.environ, {'HALO_RDKIT_GUARD': '1'}):
            # Mock to_inchikey to raise an exception (key extraction failure)
            with patch('src.halogenator.enumerate_k.to_inchikey') as mock_inchikey:
                mock_inchikey.side_effect = Exception("InChIKey extraction error")
                
                # Call enumerate_with_stats to trigger early return handling
                products, qa_stats = enumerate_with_stats("CCO", cfg)
            
            # Verify early return behavior
            self.assertEqual(len(products), 0, "Should produce no products on early return")
            
            # Verify attempt count 
            self.assertEqual(qa_stats.get('attempts', 0), 1, "Should record exactly 1 attempt")
            
            # Verify template_unsupported count
            self.assertEqual(qa_stats.get('template_unsupported', 0), 1, "Should record 1 template_unsupported")
            
            # Verify version is '2'
            self.assertEqual(qa_stats.get('version'), '2', "Early return should produce version '2'")
            
            # Verify pivots structure exists
            self.assertIn('pivots', qa_stats, "Should contain pivots structure")


if __name__ == '__main__':
    unittest.main()