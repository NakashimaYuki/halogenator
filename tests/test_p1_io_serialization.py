# -*- coding: ascii -*-
"""Test P1 I/O serialization and deserialization."""

import unittest
import tempfile
import os
from src.halogenator.enumerate_k import enumerate_products, EnumConfig
from src.halogenator.io_utils import write_table, read_table


class TestP1IOSerialization(unittest.TestCase):
    """Test I/O serialization preserves data types correctly."""
    
    def setUp(self):
        """Set up test configuration."""
        self.cfg = EnumConfig(
            k_max=2,
            halogens=('F', 'Cl'),
            constraints={'per_ring_quota': 2, 'min_graph_distance': 2, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
    def test_parquet_roundtrip_preserves_types(self):
        """Test that Parquet write/read preserves substitutions and constraints_violations types."""
        benzene_smiles = "c1ccccc1"
        products = list(enumerate_products(benzene_smiles, self.cfg))
        
        self.assertGreater(len(products), 0, "Should generate some products")
        
        with tempfile.NamedTemporaryFile(suffix='.parquet', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            # Write to Parquet
            write_table(products, tmp_path)
            
            # Read back from Parquet
            read_products = read_table(tmp_path)
            
            self.assertEqual(len(read_products), len(products), 
                           "Should read same number of records")
            
            # Check that types are preserved
            for i, (original, read_back) in enumerate(zip(products, read_products)):
                with self.subTest(record=i):
                    # substitutions should be list
                    orig_subs = original.get('substitutions', [])
                    read_subs = read_back.get('substitutions', [])
                    
                    self.assertIsInstance(orig_subs, list, "Original substitutions should be list")
                    self.assertIsInstance(read_subs, list, "Read substitutions should be list")
                    self.assertEqual(len(orig_subs), len(read_subs), "Substitutions length should match")
                    
                    # constraints_violations should be dict
                    orig_viol = original.get('constraints_violations', {})
                    read_viol = read_back.get('constraints_violations', {})
                    
                    self.assertIsInstance(orig_viol, dict, "Original violations should be dict")
                    self.assertIsInstance(read_viol, dict, "Read violations should be dict")
                    
        finally:
            os.unlink(tmp_path)
                
    def test_csv_roundtrip_preserves_types(self):
        """Test that CSV write/read preserves substitutions and constraints_violations types."""
        benzene_smiles = "c1ccccc1"
        products = list(enumerate_products(benzene_smiles, self.cfg))
        
        with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            # Write to CSV
            write_table(products, tmp_path)
            
            # Read back from CSV
            read_products = read_table(tmp_path)
            
            self.assertEqual(len(read_products), len(products),
                           "Should read same number of records")
            
            # Check that types are preserved
            for original, read_back in zip(products, read_products):
                # substitutions should be list
                read_subs = read_back.get('substitutions', [])
                self.assertIsInstance(read_subs, list, "Read substitutions should be list")
                
                # constraints_violations should be dict  
                read_viol = read_back.get('constraints_violations', {})
                self.assertIsInstance(read_viol, dict, "Read violations should be dict")
                
        finally:
            os.unlink(tmp_path)
                
    def test_substitutions_structure_preserved(self):
        """Test that substitutions list structure is fully preserved."""
        benzene_smiles = "c1ccccc1"
        products = list(enumerate_products(benzene_smiles, self.cfg))
        
        # Find a k=2 product to test multi-step history
        k2_products = [p for p in products if p.get('k') == 2]
        if not k2_products:
            self.skipTest("No k=2 products generated")
            
        product = k2_products[0]
        original_subs = product['substitutions']
        
        with tempfile.NamedTemporaryFile(suffix='.parquet', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            write_table([product], tmp_path)
            read_products = read_table(tmp_path)
            read_subs = read_products[0]['substitutions']
            
            # Check structure preservation
            self.assertEqual(len(read_subs), len(original_subs), 
                           "Should preserve substitutions count")
            
            for orig_step, read_step in zip(original_subs, read_subs):
                self.assertIsInstance(read_step, dict, "Each step should be a dict")
                
                # Check key fields are preserved
                for key in ['rule', 'halogen', 'depth']:
                    if key in orig_step:
                        self.assertIn(key, read_step, f"Should preserve {key} field")
                        self.assertEqual(orig_step[key], read_step[key], 
                                       f"Should preserve {key} value")
                        
        finally:
            os.unlink(tmp_path)
                
    def test_json_fields_present_in_storage(self):
        """Test that *_json fields are created during storage but not in memory objects."""
        benzene_smiles = "c1ccccc1"
        products = list(enumerate_products(benzene_smiles, self.cfg))
        
        if not products:
            self.skipTest("No products generated")
            
        product = products[0]
        
        # Memory objects should not have _json fields
        self.assertNotIn('substitutions_json', product, "Memory object should not have _json field")
        self.assertNotIn('constraints_violations_json', product, "Memory object should not have _json field")
        
        with tempfile.NamedTemporaryFile(suffix='.parquet', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            from src.halogenator.io_utils import _prepare_records_for_table
            import pandas as pd
            
            # Check that prepared records have _json fields
            prepared = _prepare_records_for_table([product])
            self.assertIn('substitutions_json', prepared[0], "Prepared record should have _json field")
            self.assertIn('constraints_violations_json', prepared[0], "Prepared record should have _json field")
            
            # Check that _json fields are strings
            self.assertIsInstance(prepared[0]['substitutions_json'], str, "_json field should be string")
            self.assertIsInstance(prepared[0]['constraints_violations_json'], str, "_json field should be string")
            
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)


if __name__ == '__main__':
    unittest.main()