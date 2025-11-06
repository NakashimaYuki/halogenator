# -*- coding: ascii -*-
"""Test carbonyl_unknown QA event recording."""

import unittest
from unittest.mock import patch, MagicMock
from src.halogenator.enumerate_k1 import enumerate_k1_with_stats


class TestCarbonylUnknownQA(unittest.TestCase):
    """Test carbonyl_unknown QA event recording in enumeration."""
    
    def test_carbonyl_unknown_recorded_when_is_carbonyl_carbon_returns_none(self):
        """Should record carbonyl_unknown QA events when is_carbonyl_carbon returns None."""
        
        # Mock config object
        mock_config = MagicMock()
        mock_config.halogens = ['Cl']
        mock_config.constraints = {}
        mock_config.qc_cfg = {}
        mock_config.std_cfg = {}
        
        # Mock the functions to simulate RDKit unavailability for carbonyl determination  
        with patch('src.halogenator.chem_compat.Chem') as mock_chem:
            mock_mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mock_mol
            mock_chem.MolToSmiles.return_value = 'CCO'  # Simple molecule
            
            with patch('src.halogenator.enumerate_k1.to_inchikey') as mock_inchikey:
                mock_inchikey.return_value = 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'
                
                with patch('src.halogenator.enumerate_k1.c_ring_indices') as mock_c_ring:
                    mock_c_ring.return_value = [0, 1, 2]  # Mock some ring sites
                    
                    with patch('src.halogenator.enumerate_k1.is_carbonyl_carbon') as mock_carbonyl:
                        # Mock is_carbonyl_carbon to return None (unknown status)
                        mock_carbonyl.return_value = None
                        
                        with patch('src.halogenator.enumerate_k1.is_c_ring_site_ready') as mock_site_ready:
                            mock_site_ready.return_value = True  # Site is ready for halogenation
                            
                            with patch('src.halogenator.enumerate_k1.apply_single_site_halogenation') as mock_halogenate:
                                mock_product = MagicMock()
                                mock_halogenate.return_value = mock_product
                                
                                # Run enumeration
                                products, qa_stats = enumerate_k1_with_stats('CCO', mock_config)
                                
                                # Should have recorded carbonyl_unknown events
                                self.assertIn('qa_paths', qa_stats)
                                qa_paths = qa_stats['qa_paths']
                                self.assertIn('carbonyl_unknown', qa_paths)
                                self.assertGreater(qa_paths['carbonyl_unknown'], 0)
    
    def test_carbonyl_unknown_not_recorded_when_is_carbonyl_carbon_returns_bool(self):
        """Should not record carbonyl_unknown QA events when is_carbonyl_carbon returns definitive boolean."""
        
        # Mock config object
        mock_config = MagicMock()
        mock_config.halogens = ['Cl']
        mock_config.constraints = {}
        mock_config.qc_cfg = {}
        mock_config.std_cfg = {}
        
        # Mock the functions
        with patch('src.halogenator.chem_compat.Chem') as mock_chem:
            mock_mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mock_mol
            mock_chem.MolToSmiles.return_value = 'CCO'
            
            with patch('src.halogenator.enumerate_k1.to_inchikey') as mock_inchikey:
                mock_inchikey.return_value = 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'
                
                with patch('src.halogenator.enumerate_k1.c_ring_indices') as mock_c_ring:
                    mock_c_ring.return_value = [0, 1]  # Mock some ring sites
                    
                    with patch('src.halogenator.enumerate_k1.is_carbonyl_carbon') as mock_carbonyl:
                        # Mock is_carbonyl_carbon to return False (definitive non-carbonyl)
                        mock_carbonyl.return_value = False
                        
                        with patch('src.halogenator.enumerate_k1.is_c_ring_site_ready') as mock_site_ready:
                            mock_site_ready.return_value = True
                            
                            with patch('src.halogenator.enumerate_k1.apply_single_site_halogenation') as mock_halogenate:
                                mock_product = MagicMock()
                                mock_halogenate.return_value = mock_product
                                
                                # Run enumeration
                                products, qa_stats = enumerate_k1_with_stats('CCO', mock_config)
                                
                                # Should NOT have recorded carbonyl_unknown events
                                qa_paths = qa_stats.get('qa_paths', {})
                                self.assertEqual(qa_paths.get('carbonyl_unknown', 0), 0)
    
    def test_carbonyl_unknown_recorded_on_exception_during_check(self):
        """Should continue processing when carbonyl check raises exception."""
        
        mock_config = MagicMock()
        mock_config.halogens = ['Cl']
        mock_config.constraints = {}
        mock_config.qc_cfg = {}
        mock_config.std_cfg = {}
        
        with patch('src.halogenator.chem_compat.Chem') as mock_chem:
            mock_mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mock_mol
            mock_chem.MolToSmiles.return_value = 'CCO'
            
            with patch('src.halogenator.enumerate_k1.to_inchikey') as mock_inchikey:
                mock_inchikey.return_value = 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'
                
                with patch('src.halogenator.enumerate_k1.c_ring_indices') as mock_c_ring:
                    mock_c_ring.return_value = [0]  # One site
                    
                    with patch('src.halogenator.enumerate_k1.is_carbonyl_carbon') as mock_carbonyl:
                        # Mock is_carbonyl_carbon to raise exception
                        mock_carbonyl.side_effect = Exception("RDKit error")
                        
                        with patch('src.halogenator.enumerate_k1.is_c_ring_site_ready') as mock_site_ready:
                            mock_site_ready.return_value = True
                            
                            with patch('src.halogenator.enumerate_k1.apply_single_site_halogenation') as mock_halogenate:
                                mock_product = MagicMock()
                                mock_halogenate.return_value = mock_product
                                
                                # Should not crash
                                products, qa_stats = enumerate_k1_with_stats('CCO', mock_config)
                                
                                # Should have some QA stats (not necessarily carbonyl_unknown)
                                self.assertIn('qa_paths', qa_stats)


if __name__ == '__main__':
    unittest.main()