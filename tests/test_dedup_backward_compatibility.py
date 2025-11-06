# -*- coding: ascii -*-
"""Test backward compatibility for dedup field names."""

import unittest
import tempfile
import os
from src.halogenator.cli import cmd_enum
from src.halogenator.enumerate_k import QAStats


class TestDedupBackwardCompatibility(unittest.TestCase):
    """Test backward compatibility for old and new dedup field names."""
    
    def setUp(self):
        """Set up test environment with temporary directories."""
        self.temp_dir = tempfile.mkdtemp()
        self.input_file = os.path.join(self.temp_dir, 'test_parents.smi')
        self.output_dir = os.path.join(self.temp_dir, 'output')
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Create test input file with one parent
        with open(self.input_file, 'w') as f:
            f.write('c1ccccc1\tparent1\n')
    
    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_qa_stats_old_field_compatibility(self):
        """Test that QAStats.merge() handles old field names."""
        qa_stats = QAStats()
        
        # Test merging stats with old field names only
        old_stats = {
            'no_product_matches': 2,
            'template_unsupported': 1,
            'statesig_hits': 3,  # Old field name
            'inchi_hits': 1,     # Old field name
            'qa_paths': {
                'isotope_unavailable': 1,
                'isotope_miss': 0,
                'atommap_used': 1,
                'heuristic_used': 0
            }
        }
        
        qa_stats.merge(old_stats)
        
        # Should have merged old fields into new ones
        self.assertEqual(qa_stats.dedup_hits_statesig, 3)
        self.assertEqual(qa_stats.dedup_hits_inchi, 1)
        self.assertEqual(qa_stats.no_product_matches, 2)
        self.assertEqual(qa_stats.template_unsupported, 1)
    
    def test_qa_stats_mixed_fields(self):
        """Test that QAStats.merge() handles both old and new field names."""
        qa_stats = QAStats()
        
        # Test merging stats with both old and new field names
        mixed_stats = {
            'no_product_matches': 1,
            'template_unsupported': 0,
            'dedup_hits_statesig': 2,  # New field name
            'statesig_hits': 1,        # Old field name - should be added to new
            'dedup_hits_inchi': 1,     # New field name  
            'inchi_hits': 2,           # Old field name - should be added to new
            'qa_paths': {
                'isotope_unavailable': 0,
                'isotope_miss': 1,
                'atommap_used': 0,
                'heuristic_used': 1
            }
        }
        
        qa_stats.merge(mixed_stats)
        
        # Should have summed both old and new fields
        self.assertEqual(qa_stats.dedup_hits_statesig, 2 + 1)  # new + old
        self.assertEqual(qa_stats.dedup_hits_inchi, 1 + 2)     # new + old
        self.assertEqual(qa_stats.no_product_matches, 1)
        self.assertEqual(qa_stats.template_unsupported, 0)
    
    def test_qa_stats_new_fields_only(self):
        """Test that QAStats.merge() still works with new field names only."""
        qa_stats = QAStats()
        
        # Test merging stats with new field names only
        new_stats = {
            'no_product_matches': 3,
            'template_unsupported': 2,
            'dedup_hits_statesig': 4,  # New field name only
            'dedup_hits_inchi': 3,     # New field name only
            'qa_paths': {
                'isotope_unavailable': 2,
                'isotope_miss': 1,
                'atommap_used': 1,
                'heuristic_used': 0
            }
        }
        
        qa_stats.merge(new_stats)
        
        # Should work normally with new field names
        self.assertEqual(qa_stats.dedup_hits_statesig, 4)
        self.assertEqual(qa_stats.dedup_hits_inchi, 3)
        self.assertEqual(qa_stats.no_product_matches, 3)
        self.assertEqual(qa_stats.template_unsupported, 2)
    
    def test_cli_backward_compatibility(self):
        """Test that CLI merge handles old dedup field names."""
        call_count = 0
        
        def mock_enumerate_k1_with_stats(smiles, cfg):
            nonlocal call_count
            call_count += 1
            
            # Return stats with old field names to test backward compatibility
            return [{'smiles': 'Fc1ccccc1', 'parent_smiles': 'c1ccccc1'}], {
                'version': '2',
                'attempts': 3,
                'products': 1,
                'no_product_matches': 2,
                'template_unsupported': 0,
                'statesig_hits': 2,  # Old field name
                'inchi_hits': 1,     # Old field name
                'qa_paths': {
                    'isotope_unavailable': 1,
                    'isotope_miss': 0,
                    'atommap_used': 1,
                    'heuristic_used': 0
                },
                'pivots': {
                    'by_rule_halogen_k': {
                        'R1_F_1': {'attempts': 3, 'products': 1, 'no_product_matches': 2, 'isotope_unavailable': 1, 'atommap_used': 1}
                    },
                    'by_rule': {'R1': {'attempts': 3, 'products': 1, 'no_product_matches': 2}},
                    'by_halogen': {'F': {'attempts': 3, 'products': 1, 'no_product_matches': 2}},
                    'by_k': {1: {'attempts': 3, 'products': 1, 'no_product_matches': 2}},
                    'by_rule_halogen': {'R1_F': {'attempts': 3, 'products': 1, 'no_product_matches': 2}}
                }
            }
        
        # Create config for CLI
        config = {
            'io': {
                'smiles_file': self.input_file,
                'products_table': os.path.join(self.output_dir, 'products.parquet')
            },
            'k_max': 1,
            'halogens': ['F']
        }
        
        # Create mock args to specify output directory
        class MockArgs:
            def __init__(self, outdir):
                self.outdir = outdir
                self.subset = 'flavonoids'
                self.k = 1
        
        mock_args = MockArgs(self.output_dir)
        
        # Run CLI with dependency injection
        cmd_enum(config, args=mock_args, enumerate_k1_fn=mock_enumerate_k1_with_stats)
        
        # Load and validate QA stats
        import json
        qa_file = os.path.join(self.output_dir, 'qa_summary.json')
        self.assertTrue(os.path.exists(qa_file))
        
        with open(qa_file, 'r') as f:
            qa_data = json.load(f)
        
        # Should have merged old field values into new field names
        self.assertEqual(qa_data.get('dedup_hits_statesig'), 2)  # From statesig_hits
        self.assertEqual(qa_data.get('dedup_hits_inchi'), 1)     # From inchi_hits
        self.assertEqual(qa_data.get('version'), '2')
        
        # Old field names should not appear in output
        self.assertNotIn('statesig_hits', qa_data)
        self.assertNotIn('inchi_hits', qa_data)


if __name__ == '__main__':
    unittest.main()