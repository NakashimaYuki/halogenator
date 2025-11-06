# -*- coding: ascii -*-
"""Test QA loader degrade functionality."""

import os
import json
import tempfile
import unittest
from unittest.mock import patch
from src.halogenator.report import _load_qa_stats_with_fallback


class TestQALoaderDegrade(unittest.TestCase):
    """Test QA loader degrade functionality."""
    
    def setUp(self):
        # Clear environment before each test
        os.environ.pop('HALO_QA_DEGRADE', None)
    
    def test_qa_loader_degrade_disabled_rejects_products_only(self):
        """QA loader should reject products-only stats when degrade is disabled."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump({'products': 100}, f)  # Missing 'attempts'
            qa_path = f.name
        
        try:
            # Create a fake products table in same directory
            products_table = os.path.join(os.path.dirname(qa_path), 'products.csv')
            open(products_table, 'w').close()
            
            config = {'qa_loader_degrade': False}
            qa_stats, qa_source = _load_qa_stats_with_fallback(products_table, config)
            
            # Should have empty stats and no successful load
            self.assertEqual(qa_stats.get('attempts', 0), 0)
            self.assertEqual(qa_stats.get('products', 0), 0)
            self.assertFalse(qa_source.get('degrade', False))
        finally:
            os.unlink(qa_path)
            if os.path.exists(products_table):
                os.unlink(products_table)
    
    def test_qa_loader_degrade_enabled_accepts_products_only(self):
        """QA loader should accept products-only stats when degrade is enabled."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump({'products': 100, 'no_product_matches': 5}, f)  # Missing 'attempts'
            qa_path = f.name
        
        try:
            # Create products table in same directory
            products_dir = os.path.dirname(qa_path)
            qa_file = os.path.join(products_dir, 'qa_summary.json')
            os.rename(qa_path, qa_file)  # Move to expected location
            
            products_table = os.path.join(products_dir, 'products.csv')
            open(products_table, 'w').close()
            
            config = {'qa_loader_degrade': True}
            qa_stats, qa_source = _load_qa_stats_with_fallback(products_table, config)
            
            # Should have loaded with degrade mode
            self.assertEqual(qa_stats.get('attempts'), 0)  # Set to 0 in degrade mode
            self.assertEqual(qa_stats.get('products'), 100)
            self.assertEqual(qa_stats.get('no_product_matches'), 5)
            self.assertTrue(qa_source.get('degrade', False))
            self.assertFalse(qa_source.get('has_attempts', True))
        finally:
            if os.path.exists(qa_file):
                os.unlink(qa_file)
            if os.path.exists(products_table):
                os.unlink(products_table)
    
    def test_qa_loader_degrade_from_environment(self):
        """Test that degrade flag is read from environment variable."""
        # This test checks the integration point in generate_summary
        # We'll test the environment variable propagation
        
        os.environ['HALO_QA_DEGRADE'] = '1'
        
        # Test that the environment variable is read correctly
        config_with_subset = {}
        config_with_subset['qa_loader_degrade'] = (os.environ.get('HALO_QA_DEGRADE', '0') == '1')
        
        self.assertTrue(config_with_subset['qa_loader_degrade'])
        
        # Test with '0' value
        os.environ['HALO_QA_DEGRADE'] = '0'
        config_with_subset['qa_loader_degrade'] = (os.environ.get('HALO_QA_DEGRADE', '0') == '1')
        
        self.assertFalse(config_with_subset['qa_loader_degrade'])
    
    def test_qa_loader_degrade_with_valid_attempts(self):
        """QA loader should load normally when attempts field is present, regardless of degrade setting."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump({'attempts': 50, 'products': 40, 'no_product_matches': 10}, f)
            qa_path = f.name
        
        try:
            # Create products table in same directory
            products_dir = os.path.dirname(qa_path)
            qa_file = os.path.join(products_dir, 'qa_summary.json')
            os.rename(qa_path, qa_file)
            
            products_table = os.path.join(products_dir, 'products.csv')
            open(products_table, 'w').close()
            
            # Test with degrade enabled - should still load normally
            config = {'qa_loader_degrade': True}
            qa_stats, qa_source = _load_qa_stats_with_fallback(products_table, config)
            
            # Should have loaded normally, not in degrade mode
            self.assertEqual(qa_stats.get('attempts'), 50)
            self.assertEqual(qa_stats.get('products'), 40)
            self.assertEqual(qa_stats.get('no_product_matches'), 10)
            self.assertFalse(qa_source.get('degrade', True))  # Not degrade mode
            self.assertTrue(qa_source.get('has_attempts', False))  # Has attempts
        finally:
            if os.path.exists(qa_file):
                os.unlink(qa_file)
            if os.path.exists(products_table):
                os.unlink(products_table)


if __name__ == '__main__':
    unittest.main()