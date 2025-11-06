# -*- coding: ascii -*-
"""Tests for deep serialization safety of stats structures."""

import json
import tempfile
import unittest
import pandas as pd
from collections import defaultdict, Counter
from src.halogenator.report import _init_incremental_stats, _finalize_incremental_stats


class TestSerializationSafety(unittest.TestCase):
    """Test that stats structures are fully serializable."""
    
    def test_stats_json_serializable(self):
        """Test that finalized stats can be serialized to JSON."""
        # Create mock parent pairs
        parent_pairs = [
            ('c1ccccc1O', 'phenol'),
            ('c1ccc(O)cc1O', 'resorcinol'),
            ('Cc1ccccc1N', 'toluidine')
        ]
        
        # Initialize stats with some problematic structures
        stats = _init_incremental_stats()
        
        # Add some data that includes sets and other non-serializable types
        stats['parents_smiles_with_products'].add('c1ccccc1O')
        stats['parents_smiles_with_products'].add('Cc1ccccc1N')
        stats['parents_inchikey_with_products'].add('ISWSIDIOOBJBQZ-UHFFFAOYSA-N')
        
        # Mock config
        mock_config = {'k_max': 1, 'subset': 'test'}
        
        # Finalize stats (this should convert sets to lists)
        _finalize_incremental_stats(stats, parent_pairs, mock_config)
        
        # Test JSON serialization
        try:
            json_str = json.dumps(stats)
            self.assertIsInstance(json_str, str)
            
            # Test round-trip
            stats_restored = json.loads(json_str)
            self.assertIsInstance(stats_restored, dict)
            
        except TypeError as e:
            self.fail(f"Stats structure not JSON serializable: {e}")
    
    def test_stats_to_dataframe(self):
        """Test that key stats fields can be used in DataFrame construction."""
        parent_pairs = [('c1ccccc1O', 'phenol')]
        
        stats = _init_incremental_stats()
        stats['parents_smiles_with_products'].add('c1ccccc1O')
        stats['product_count'] = 5
        stats['parent_count'] = 1
        
        mock_config = {'k_max': 1, 'subset': 'test'}
        _finalize_incremental_stats(stats, parent_pairs, mock_config)
        
        # Test DataFrame construction with key fields
        try:
            key_stats = {
                'product_count': stats['product_count'],
                'parent_count': stats['parent_count'],
                'unique_parents_with_products': stats['unique_parents_with_products'],
                'parents_smiles_with_products': stats['parents_smiles_with_products']
            }
            
            # Should not raise any serialization errors
            df = pd.DataFrame([key_stats])
            self.assertEqual(len(df), 1)
            
        except Exception as e:
            self.fail(f"Stats not compatible with DataFrame: {e}")
    
    def test_all_nested_structures_are_basic_types(self):
        """Test that all nested structures contain only basic serializable types."""
        parent_pairs = [('c1ccccc1O', 'phenol'), ('c1ccc(N)cc1', 'aniline')]
        
        stats = _init_incremental_stats()
        # Simulate some complex nested data
        stats['parents_smiles_with_products'].add('c1ccccc1O')
        stats['parents_inchikey_with_products'].add('ISWSIDIOOBJBQZ-UHFFFAOYSA-N')
        stats['all_parent_keys'].add('c1ccccc1O')
        stats['parents_with_products'].add('c1ccccc1O')
        
        mock_config = {'k_max': 1, 'subset': 'test'}
        _finalize_incremental_stats(stats, parent_pairs, mock_config)
        
        # Recursively check that all values are basic types
        def check_basic_types(obj, path=""):
            if isinstance(obj, (str, int, float, bool, type(None))):
                return True
            elif isinstance(obj, list):
                for i, item in enumerate(obj):
                    if not check_basic_types(item, f"{path}[{i}]"):
                        return False
                return True
            elif isinstance(obj, dict):
                for key, value in obj.items():
                    if not isinstance(key, (str, int)):
                        self.fail(f"Non-basic key type at {path}.{key}: {type(key)}")
                    if not check_basic_types(value, f"{path}.{key}"):
                        return False
                return True
            else:
                self.fail(f"Non-basic type at {path}: {type(obj)} = {obj}")
        
        # Check the entire stats structure
        self.assertTrue(check_basic_types(stats), "Stats contains non-basic types")
    
    def test_sets_converted_to_sorted_lists(self):
        """Test that sets are converted to sorted lists."""
        parent_pairs = [('c1ccccc1O', 'phenol'), ('c1ccc(N)cc1', 'aniline')]
        
        stats = _init_incremental_stats()
        # Add items in random order to test sorting
        items = ['zebra', 'alpha', 'beta', 'charlie']
        for item in items:
            stats['parents_smiles_with_products'].add(item)
        
        mock_config = {'k_max': 1, 'subset': 'test'}
        _finalize_incremental_stats(stats, parent_pairs, mock_config)
        
        # Check conversion and sorting
        result = stats['parents_smiles_with_products']
        self.assertIsInstance(result, list, "Should be converted to list")
        self.assertEqual(result, sorted(items), "Should be sorted")


if __name__ == '__main__':
    unittest.main()