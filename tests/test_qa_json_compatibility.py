# -*- coding: ascii -*-
"""Tests to ensure QA JSON reading compatibility between legacy and new formats."""

import unittest
import json
import os
import tempfile
from src.halogenator.report import generate_summary_report


class TestQAJSONCompatibility(unittest.TestCase):
    """Test that QA JSON reading supports both legacy and new formats consistently."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Sample QA statistics data
        self.qa_stats = {
            'attempts': 100,
            'products': 75,
            'no_product_matches': 15,
            'template_unsupported': 10,
            'dedup_hits_statesig': 5,
            'dedup_hits_inchi': 3,
            'qa_paths': {
                'isotope_unavailable': 2,
                'isotope_miss': 1,
                'atommap_used': 45,
                'heuristic_used': 27
            }
        }
        
        # Create temporary directory for test files
        self.temp_dir = tempfile.mkdtemp()
        
    def tearDown(self):
        """Clean up test fixtures."""
        # Clean up temporary directory
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_legacy_qa_json_format_reading(self):
        """Test that legacy QA JSON format (root-level totals) is read correctly."""
        # Legacy format: QA statistics at root level
        legacy_qa_data = self.qa_stats.copy()
        
        # Write legacy format QA JSON
        qa_json_path = os.path.join(self.temp_dir, 'qa_summary.json')
        with open(qa_json_path, 'w') as f:
            json.dump(legacy_qa_data, f, indent=2)
        
        # Read and verify using the unified pattern
        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)
        
        # Apply unified reading pattern
        qa_stats = qa_data.get('total', qa_data)
        
        # Verify all fields are accessible
        self.assertEqual(qa_stats['attempts'], 100)
        self.assertEqual(qa_stats['products'], 75)
        self.assertEqual(qa_stats['no_product_matches'], 15)
        self.assertEqual(qa_stats['template_unsupported'], 10)
        self.assertEqual(qa_stats['qa_paths']['atommap_used'], 45)
        
    def test_new_qa_json_format_reading(self):
        """Test that new QA JSON format (nested structure) is read correctly."""
        # New format: QA statistics nested under 'total' key
        new_qa_data = {
            'version': '2',
            'total': self.qa_stats.copy(),
            'pivots': {
                'by_rule': {'R1': 45, 'R2': 30},
                'by_halogen': {'F': 50, 'Cl': 25}
            },
            'metadata': {
                'generated_at': '2025-01-10T12:00:00Z'
            }
        }
        
        # Write new format QA JSON
        qa_json_path = os.path.join(self.temp_dir, 'qa_summary.json')
        with open(qa_json_path, 'w') as f:
            json.dump(new_qa_data, f, indent=2)
        
        # Read and verify using the unified pattern
        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)
        
        # Apply unified reading pattern
        qa_stats = qa_data.get('total', qa_data)
        
        # Verify all fields are accessible (should be same as legacy)
        self.assertEqual(qa_stats['attempts'], 100)
        self.assertEqual(qa_stats['products'], 75)
        self.assertEqual(qa_stats['no_product_matches'], 15)
        self.assertEqual(qa_stats['template_unsupported'], 10)
        self.assertEqual(qa_stats['qa_paths']['atommap_used'], 45)
        
    def test_both_formats_produce_consistent_results(self):
        """Test that both QA JSON formats produce identical results when processed."""
        # Legacy format
        legacy_qa_data = self.qa_stats.copy()
        
        # New format
        new_qa_data = {
            'version': '2',
            'total': self.qa_stats.copy(),
            'pivots': {'by_rule': {'R1': 45}},
            'metadata': {'generated_at': '2025-01-10T12:00:00Z'}
        }
        
        # Process both formats using unified pattern
        legacy_processed = legacy_qa_data.get('total', legacy_qa_data)
        new_processed = new_qa_data.get('total', new_qa_data)
        
        # Results should be identical
        self.assertEqual(legacy_processed['attempts'], new_processed['attempts'])
        self.assertEqual(legacy_processed['products'], new_processed['products'])
        self.assertEqual(legacy_processed['no_product_matches'], new_processed['no_product_matches'])
        self.assertEqual(legacy_processed['template_unsupported'], new_processed['template_unsupported'])
        self.assertEqual(legacy_processed['qa_paths'], new_processed['qa_paths'])
        
    def test_missing_total_key_falls_back_correctly(self):
        """Test that missing 'total' key correctly falls back to root data."""
        # Data without 'total' key (legacy format)
        legacy_data = self.qa_stats.copy()
        
        # Apply unified reading pattern
        result = legacy_data.get('total', legacy_data)
        
        # Should return the entire legacy_data since no 'total' key exists
        self.assertEqual(result, legacy_data)
        self.assertEqual(result['attempts'], 100)
        self.assertEqual(result['products'], 75)
        
    def test_empty_total_key_handling(self):
        """Test that empty or None 'total' key is handled gracefully."""
        # Data with empty 'total' key
        problematic_data = {
            'total': None,
            'attempts': 100,
            'products': 75
        }
        
        # Apply unified reading pattern
        result = problematic_data.get('total', problematic_data)
        
        # Should get None since 'total' exists but is None
        self.assertIsNone(result)
        
        # But the fallback should work when 'total' key doesn't exist
        del problematic_data['total']
        result = problematic_data.get('total', problematic_data)
        self.assertEqual(result['attempts'], 100)
        self.assertEqual(result['products'], 75)


if __name__ == '__main__':
    unittest.main()