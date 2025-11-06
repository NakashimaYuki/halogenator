# -*- coding: ascii -*-
"""Schema stability tests for QA statistics contract."""

import unittest
from src.halogenator.enumerate_k import enumerate_with_stats, EnumConfig, QAAggregator


class TestSchemaStability(unittest.TestCase):
    """Test that QA statistics maintain stable schema regardless of data content."""
    
    def setUp(self):
        """Set up test environment."""
        self.cfg = EnumConfig(
            halogens=['F'],
            k_max=1,
            constraints={},
            std_cfg={},
            qc_cfg={},
            pruning_cfg={}
        )
    
    def test_final_qa_stats_v2_schema_stable(self):
        """Test that final_qa_stats_v2 always contains required fields."""
        # Test with valid molecule that should produce some activity
        products, qa_stats = enumerate_with_stats('c1ccccc1', self.cfg)  # benzene
        
        # Version 2 schema requirements
        self.assertEqual(qa_stats.get('version'), '2')
        
        # Pivots dimension - always present (can be empty dict)
        self.assertIn('pivots', qa_stats)
        pivots = qa_stats['pivots']
        self.assertIsInstance(pivots, dict)
        
        # All pivot dimensions must be present
        required_dimensions = ['by_rule', 'by_halogen', 'by_k', 'by_rule_halogen', 'by_rule_halogen_k']
        for dim in required_dimensions:
            self.assertIn(dim, pivots, f"Pivot dimension '{dim}' must always be present")
            self.assertIsInstance(pivots[dim], dict, f"Pivot dimension '{dim}' must be dict")
        
        # Core totals - always present (can be zero)
        required_totals = ['attempts', 'products', 'no_product_matches', 'template_unsupported']
        for field in required_totals:
            self.assertIn(field, qa_stats, f"Total field '{field}' must always be present")
            self.assertIsInstance(qa_stats[field], int, f"Total field '{field}' must be int")
            self.assertGreaterEqual(qa_stats[field], 0, f"Total field '{field}' must be >= 0")
        
        # QA paths - always present with all subkeys (can be zero)
        self.assertIn('qa_paths', qa_stats)
        qa_paths = qa_stats['qa_paths']
        self.assertIsInstance(qa_paths, dict)
        
        required_qa_paths = ['isotope_unavailable', 'isotope_miss', 'atommap_used', 'heuristic_used']
        for qa_key in required_qa_paths:
            self.assertIn(qa_key, qa_paths, f"QA path '{qa_key}' must always be present")
            self.assertIsInstance(qa_paths[qa_key], int, f"QA path '{qa_key}' must be int")
            self.assertGreaterEqual(qa_paths[qa_key], 0, f"QA path '{qa_key}' must be >= 0")
        
        # Legacy dedup stats - always present (can be zero)
        legacy_fields = ['dedup_hits_statesig', 'dedup_hits_inchi']
        for field in legacy_fields:
            self.assertIn(field, qa_stats, f"Legacy field '{field}' must always be present")
            self.assertIsInstance(qa_stats[field], int, f"Legacy field '{field}' must be int")
            self.assertGreaterEqual(qa_stats[field], 0, f"Legacy field '{field}' must be >= 0")
    
    def test_empty_enumeration_stable_schema(self):
        """Test schema stability even when no enumeration occurs."""
        # Test with molecule that produces no results 
        products, qa_stats = enumerate_with_stats('invalid_smiles', self.cfg)
        
        # Should still have stable schema even with no activity
        self.assertEqual(qa_stats.get('version'), '2')
        self.assertIn('pivots', qa_stats)
        
        # All required fields must exist even if zero
        required_fields = [
            'attempts', 'products', 'no_product_matches', 'template_unsupported',
            'qa_paths', 'dedup_hits_statesig', 'dedup_hits_inchi'
        ]
        for field in required_fields:
            self.assertIn(field, qa_stats, f"Field '{field}' must exist even with no activity")
    
    def test_no_sentinel_key_pollution(self):
        """Test that no sentinel keys (like '*_already_counted') appear in output."""
        products, qa_stats = enumerate_with_stats('c1ccccc1', self.cfg)
        
        # Check pivots for any sentinel key pollution
        pivots = qa_stats['pivots']
        
        def check_no_sentinels(obj, path=""):
            """Recursively check for sentinel keys."""
            if isinstance(obj, dict):
                for key, value in obj.items():
                    # No sentinel keys should exist (only check string keys)
                    if isinstance(key, str):
                        self.assertNotIn('_already_counted', key, 
                                       f"Sentinel key '{key}' found at {path}")
                        self.assertNotIn('_seen', key,
                                       f"Sentinel key '{key}' found at {path}")
                    check_no_sentinels(value, f"{path}.{key}")
        
        check_no_sentinels(pivots, "pivots")
        
        # Also check top-level qa_stats
        for key in qa_stats:
            self.assertNotIn('_already_counted', key, f"Sentinel key '{key}' found in qa_stats")
            self.assertNotIn('_seen', key, f"Sentinel key '{key}' found in qa_stats")
    
    def test_aggregator_stable_schema(self):
        """Test QAAggregator produces stable schema."""
        aggregator = QAAggregator()
        
        # Even empty aggregator should produce stable pivots
        pivots = aggregator.to_pivots_dict()
        
        required_dimensions = ['by_rule', 'by_halogen', 'by_k', 'by_rule_halogen', 'by_rule_halogen_k']
        for dim in required_dimensions:
            self.assertIn(dim, pivots, f"Aggregator dimension '{dim}' must always be present")
            self.assertIsInstance(pivots[dim], dict, f"Aggregator dimension '{dim}' must be dict")


if __name__ == '__main__':
    unittest.main()