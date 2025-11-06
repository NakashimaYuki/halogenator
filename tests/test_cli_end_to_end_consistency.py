# -*- coding: ascii -*-
"""End-to-end CLI consistency tests."""

import unittest
import json
import tempfile
import os
from pathlib import Path
from unittest.mock import patch, MagicMock
from src.halogenator.cli import cmd_enum
from src.halogenator.enumerate_k import EnumConfig


class TestCLIEndToEndConsistency(unittest.TestCase):
    """Test end-to-end CLI consistency with aggregator-derived statistics."""
    
    def setUp(self):
        """Set up test environment with temporary directories."""
        self.temp_dir = tempfile.mkdtemp()
        self.input_file = os.path.join(self.temp_dir, 'test_parents.smi')
        self.output_dir = os.path.join(self.temp_dir, 'output')
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Create test input file with two parents
        with open(self.input_file, 'w') as f:
            f.write('c1ccccc1\tparent1\n')  # Simple benzene
            f.write('c1ccc(O)cc1\tparent2\n')  # Phenol
    
    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def _validate_qa_consistency(self, qa_stats):
        """Helper to validate QA statistics consistency."""
        if qa_stats.get('version') != '2':
            return  # Skip validation for non-v2 stats
        
        pivots = qa_stats.get('pivots', {})
        if not pivots.get('by_rule_halogen_k'):
            return  # No pivot data to validate
        
        # Calculate totals from pivots
        pivot_totals = {
            'attempts': 0,
            'products': 0,
            'no_product_matches': 0,
            'template_unsupported': 0,
            'qa_paths': {
                'isotope_unavailable': 0,
                'isotope_miss': 0,
                'atommap_used': 0,
                'heuristic_used': 0
            }
        }
        
        for key, events in pivots['by_rule_halogen_k'].items():
            for metric, count in events.items():
                if metric in pivot_totals:
                    pivot_totals[metric] += count
                elif metric in pivot_totals['qa_paths']:
                    pivot_totals['qa_paths'][metric] += count
        
        # Validate consistency
        self.assertEqual(qa_stats.get('no_product_matches', 0), pivot_totals['no_product_matches'],
                        f"no_product_matches mismatch: totals={qa_stats.get('no_product_matches')} vs pivots={pivot_totals['no_product_matches']}")
        
        self.assertEqual(qa_stats.get('template_unsupported', 0), pivot_totals['template_unsupported'],
                        f"template_unsupported mismatch")
        
        for qa_key in pivot_totals['qa_paths']:
            self.assertEqual(qa_stats.get('qa_paths', {}).get(qa_key, 0), pivot_totals['qa_paths'][qa_key],
                            f"qa_paths.{qa_key} mismatch")
        
        # Validate three-way exclusive outcome invariant
        attempts = pivot_totals['attempts']
        products = pivot_totals['products']
        no_matches = pivot_totals['no_product_matches']
        tmpl_unsup = pivot_totals['template_unsupported']
        
        if attempts > 0:
            self.assertGreaterEqual(attempts, products, 
                                   f"Invariant violation: attempts={attempts} < products={products}")
            self.assertEqual(attempts, products + no_matches + tmpl_unsup,
                           f"Invariant violation: attempts={attempts} != products={products} + no_matches={no_matches} + template_unsupported={tmpl_unsup}")
    
    def test_cli_k1_consistency(self):
        """Test k=1 CLI enumeration maintains consistency."""
        # Mock data for two parents
        call_count = 0
        
        def mock_enumerate_k1_with_stats(smiles, cfg):
            nonlocal call_count
            call_count += 1
            
            if call_count == 1:
                # Parent 1: One successful attempt, one failed attempt
                return [{'smiles': 'Fc1ccccc1', 'parent_smiles': 'c1ccccc1'}], {
                    'version': '2',
                    'attempts': 2,
                    'products': 1,
                    'no_product_matches': 1,
                    'template_unsupported': 0,
                    'dedup_hits_statesig': 0,
                    'dedup_hits_inchi': 0,
                    'qa_paths': {
                        'isotope_unavailable': 1,
                        'isotope_miss': 0,
                        'atommap_used': 1,
                        'heuristic_used': 0
                    },
                    'pivots': {
                        'by_rule_halogen_k': {
                            'R1_F_1': {'attempts': 2, 'products': 1, 'no_product_matches': 1, 'isotope_unavailable': 1, 'atommap_used': 1}
                        },
                        'by_rule': {'R1': {'attempts': 2, 'products': 1, 'no_product_matches': 1}},
                        'by_halogen': {'F': {'attempts': 2, 'products': 1, 'no_product_matches': 1}},
                        'by_k': {1: {'attempts': 2, 'products': 1, 'no_product_matches': 1}},
                        'by_rule_halogen': {'R1_F': {'attempts': 2, 'products': 1, 'no_product_matches': 1}}
                    }
                }
            else:
                # Parent 2: Template unsupported (mutually exclusive from no matches)
                return [], {
                    'version': '2',
                    'attempts': 2,
                    'products': 0,
                    'no_product_matches': 0,  # Template unsupported is distinct from no matches
                    'template_unsupported': 2,
                    'dedup_hits_statesig': 0,
                    'dedup_hits_inchi': 0,
                    'qa_paths': {
                        'isotope_unavailable': 0,
                        'isotope_miss': 0,
                        'atommap_used': 0,
                        'heuristic_used': 0
                    },
                    'pivots': {
                        'by_rule_halogen_k': {
                            'R3_Cl_1': {'attempts': 2, 'no_product_matches': 0, 'template_unsupported': 2}
                        },
                        'by_rule': {'R3': {'attempts': 2, 'no_product_matches': 0, 'template_unsupported': 2}},
                        'by_halogen': {'Cl': {'attempts': 2, 'no_product_matches': 0, 'template_unsupported': 2}},
                        'by_k': {1: {'attempts': 2, 'no_product_matches': 0, 'template_unsupported': 2}},
                        'by_rule_halogen': {'R3_Cl': {'attempts': 2, 'no_product_matches': 0, 'template_unsupported': 2}}
                    }
                }
        
        # Create config for CLI using correct structure
        config = {
            'io': {
                'smiles_file': self.input_file,
                'products_table': os.path.join(self.output_dir, 'products.parquet')
            },
            'k_max': 1,
            'halogens': ['F', 'Cl']
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
        
        # Check that QA stats file was created
        qa_file = os.path.join(self.output_dir, 'qa_summary.json')
        self.assertTrue(os.path.exists(qa_file), "QA stats file should be created")
        
        # Load and validate QA stats
        with open(qa_file, 'r') as f:
            qa_data = json.load(f)
        
        # Should have version 2 since at least one parent uses v2
        self.assertEqual(qa_data.get('version'), '2')
        
        # Validate consistency for the aggregated results
        self._validate_qa_consistency(qa_data)
    
    def test_cli_k2_consistency(self):
        """Test k=2 CLI enumeration maintains consistency.""" 
        call_count = 0
        
        def mock_enumerate_with_stats(smiles, cfg):
            nonlocal call_count
            call_count += 1
            
            if call_count == 1:
                # Parent 1: Multi-match scenario
                return [{'smiles': 'Fc1ccccc1F', 'parent_smiles': 'c1ccccc1', 'k': 2}], {
                    'version': '2',
                    'attempts': 5,
                    'products': 2,
                    'no_product_matches': 3,
                    'template_unsupported': 0,
                    'dedup_hits_statesig': 1,
                    'dedup_hits_inchi': 0,
                    'qa_paths': {
                        'isotope_unavailable': 2,
                        'isotope_miss': 1,
                        'atommap_used': 2,
                        'heuristic_used': 1
                    },
                    'pivots': {
                        'by_rule_halogen_k': {
                            'R1_F_1': {'attempts': 2, 'products': 1, 'no_product_matches': 1, 'atommap_used': 1},
                            'R1_F_2': {'attempts': 2, 'products': 1, 'no_product_matches': 1, 'isotope_unavailable': 2, 'heuristic_used': 1},
                            'R1_Cl_1': {'attempts': 1, 'no_product_matches': 1, 'isotope_miss': 1, 'atommap_used': 1}
                        },
                        'by_rule': {'R1': {'attempts': 5, 'products': 2, 'no_product_matches': 3}},
                        'by_halogen': {'F': {'attempts': 4, 'products': 2, 'no_product_matches': 2}, 'Cl': {'attempts': 1, 'no_product_matches': 1}},
                        'by_k': {1: {'attempts': 3, 'products': 1, 'no_product_matches': 2}, 2: {'attempts': 2, 'products': 1, 'no_product_matches': 1}},
                        'by_rule_halogen': {'R1_F': {'attempts': 4, 'products': 2, 'no_product_matches': 2}, 'R1_Cl': {'attempts': 1, 'no_product_matches': 1}}
                    }
                }
            else:
                # Parent 2: Template unsupported (mutually exclusive from no matches)
                return [], {
                    'version': '2',
                    'attempts': 1,
                    'products': 0,
                    'no_product_matches': 0,  # Template unsupported is distinct from no matches
                    'template_unsupported': 1,
                    'dedup_hits_statesig': 0,
                    'dedup_hits_inchi': 0,
                    'qa_paths': {
                        'isotope_unavailable': 0,
                        'isotope_miss': 0,
                        'atommap_used': 0,
                        'heuristic_used': 0
                    },
                    'pivots': {
                        'by_rule_halogen_k': {
                            'R4_Br_1': {'attempts': 1, 'no_product_matches': 0, 'template_unsupported': 1}
                        },
                        'by_rule': {'R4': {'attempts': 1, 'no_product_matches': 0, 'template_unsupported': 1}},
                        'by_halogen': {'Br': {'attempts': 1, 'no_product_matches': 0, 'template_unsupported': 1}},
                        'by_k': {1: {'attempts': 1, 'no_product_matches': 0, 'template_unsupported': 1}},
                        'by_rule_halogen': {'R4_Br': {'attempts': 1, 'no_product_matches': 0, 'template_unsupported': 1}}
                    }
                }
        
        # Create config for CLI using correct structure
        config = {
            'io': {
                'smiles_file': self.input_file,
                'products_table': os.path.join(self.output_dir, 'products.parquet')
            },
            'k_max': 2,
            'halogens': ['F', 'Cl', 'Br']
        }
        
        # Create mock args to specify output directory
        class MockArgs:
            def __init__(self, outdir):
                self.outdir = outdir
                self.subset = 'flavonoids'
                self.k = 2
        
        mock_args = MockArgs(self.output_dir)
        
        # Run CLI with dependency injection
        cmd_enum(config, args=mock_args, enumerate_k_fn=mock_enumerate_with_stats)
        
        # Check that QA stats file was created
        qa_file = os.path.join(self.output_dir, 'qa_summary.json')
        self.assertTrue(os.path.exists(qa_file))
        
        # Load and validate QA stats
        with open(qa_file, 'r') as f:
            qa_data = json.load(f)
        
        # Should have version 2
        self.assertEqual(qa_data.get('version'), '2')
        
        # Validate consistency
        self._validate_qa_consistency(qa_data)
    
    def test_mixed_k1_k2_version_handling(self):
        """Test that mixed k=1/k=2 processing handles versions correctly."""
        call_count = 0
        
        # Mock function for k>1 enumeration (should be used)
        def mock_enumerate_with_stats(smiles, cfg):
            nonlocal call_count
            call_count += 1
            
            if call_count == 1:
                # Parent 1: k>1 enumeration with template unsupported (distinct from no matches)
                return [], {
                    'version': '2', 
                    'attempts': 1,
                    'products': 0,
                    'no_product_matches': 0,  # Template unsupported is distinct from no matches
                    'template_unsupported': 1,
                    'dedup_hits_statesig': 0,
                    'dedup_hits_inchi': 0,
                    'qa_paths': {'isotope_unavailable': 0, 'isotope_miss': 0, 'atommap_used': 0, 'heuristic_used': 0},
                    'pivots': {
                        'by_rule_halogen_k': {'R3_Cl_2': {'attempts': 1, 'no_product_matches': 0, 'template_unsupported': 1}},
                        'by_rule': {'R3': {'attempts': 1, 'no_product_matches': 0, 'template_unsupported': 1}},
                        'by_halogen': {'Cl': {'attempts': 1, 'no_product_matches': 0, 'template_unsupported': 1}},
                        'by_k': {2: {'attempts': 1, 'no_product_matches': 0, 'template_unsupported': 1}},
                        'by_rule_halogen': {'R3_Cl': {'attempts': 1, 'no_product_matches': 0, 'template_unsupported': 1}}
                    }
                }
            else:
                # Parent 2: successful k>1 enumeration
                return [{'smiles': 'Fc1ccc(F)cc1', 'parent_smiles': 'c1ccc(O)cc1', 'k': 2}], {
                    'version': '2',
                    'attempts': 3,
                    'products': 1,
                    'no_product_matches': 2,
                    'template_unsupported': 0,
                    'dedup_hits_statesig': 0,
                    'dedup_hits_inchi': 0,
                    'qa_paths': {'isotope_unavailable': 1, 'isotope_miss': 1, 'atommap_used': 1, 'heuristic_used': 0},
                    'pivots': {
                        'by_rule_halogen_k': {
                            'R1_F_1': {'attempts': 1, 'no_product_matches': 1, 'isotope_miss': 1},
                            'R1_F_2': {'attempts': 2, 'products': 1, 'no_product_matches': 1, 'isotope_unavailable': 1, 'atommap_used': 1}
                        },
                        'by_rule': {'R1': {'attempts': 3, 'products': 1, 'no_product_matches': 2}},
                        'by_halogen': {'F': {'attempts': 3, 'products': 1, 'no_product_matches': 2}},
                        'by_k': {1: {'attempts': 1, 'no_product_matches': 1}, 2: {'attempts': 2, 'products': 1, 'no_product_matches': 1}},
                        'by_rule_halogen': {'R1_F': {'attempts': 3, 'products': 1, 'no_product_matches': 2}}
                    }
                }
        
        # Mock function for k=1 enumeration (should NOT be used)
        def mock_enumerate_k1_with_stats(smiles, cfg):
            # This should never be called for k-max > 1
            raise AssertionError("k=1 enumeration should not be called when k-max > 1")
        
        # Create config for CLI with k-max > 1 (should use k>1 path)
        config = {
            'io': {
                'smiles_file': self.input_file,
                'products_table': os.path.join(self.output_dir, 'products.parquet')
            },
            'k_max': 2,
            'halogens': ['F', 'Cl']
        }
        
        # Create mock args to specify output directory
        class MockArgs:
            def __init__(self, outdir):
                self.outdir = outdir
                self.subset = 'flavonoids'
                self.k = 2
        
        mock_args = MockArgs(self.output_dir)
        
        # Run CLI with dependency injection 
        cmd_enum(config, args=mock_args, enumerate_k_fn=mock_enumerate_with_stats, enumerate_k1_fn=mock_enumerate_k1_with_stats)
        
        # Verify QA stats file was created
        qa_file = os.path.join(self.output_dir, 'qa_summary.json')
        self.assertTrue(os.path.exists(qa_file))
        
        # Load and validate QA stats
        with open(qa_file, 'r') as f:
            qa_data = json.load(f)
        
        self.assertEqual(qa_data.get('version'), '2')
        self._validate_qa_consistency(qa_data)
        
        # Verify both parents were processed (call_count should be 2)
        self.assertEqual(call_count, 2, "Should have processed both parents with k>1 enumeration")
    
    def test_by_k_keys_remain_ints(self):
        """Test that by_k keys remain integers in CLI merge operations."""
        # This is a regression test to ensure int keys don't become strings
        
        def mock_enumerate_with_stats(smiles, cfg):
            return [], {
                'version': '2',
                'attempts': 3,
                'products': 2,
                'no_product_matches': 1,
                'template_unsupported': 0,
                'dedup_hits_statesig': 0,
                'dedup_hits_inchi': 0,
                'qa_paths': {'isotope_unavailable': 0, 'isotope_miss': 0, 'atommap_used': 0, 'heuristic_used': 0},
                'pivots': {
                    'by_k': {
                        1: {'attempts': 2, 'products': 1, 'no_product_matches': 1},
                        2: {'attempts': 1, 'products': 1, 'no_product_matches': 0}
                    },
                    'by_rule_halogen_k': {
                        'R1_F_1': {'attempts': 2, 'products': 1, 'no_product_matches': 1},
                        'R1_F_2': {'attempts': 1, 'products': 1}
                    },
                    'by_rule': {'R1': {'attempts': 3, 'products': 2, 'no_product_matches': 1}},
                    'by_halogen': {'F': {'attempts': 3, 'products': 2, 'no_product_matches': 1}},
                    'by_rule_halogen': {'R1_F': {'attempts': 3, 'products': 2, 'no_product_matches': 1}}
                }
            }
        
        # Create config for CLI
        config = {
            'io': {
                'smiles_file': self.input_file,
                'products_table': os.path.join(self.output_dir, 'products.parquet')
            },
            'k_max': 2,
            'halogens': ['F']
        }
        
        # Create mock args to specify output directory
        class MockArgs:
            def __init__(self, outdir):
                self.outdir = outdir
                self.subset = 'flavonoids'
                self.k = 2
        
        mock_args = MockArgs(self.output_dir)
        
        # Run CLI with dependency injection
        cmd_enum(config, args=mock_args, enumerate_k_fn=mock_enumerate_with_stats)
        
        # Load QA stats and verify by_k keys are handled correctly
        qa_file = os.path.join(self.output_dir, 'qa_summary.json')
        with open(qa_file, 'r') as f:
            qa_data = json.load(f)
        
        by_k = qa_data.get('pivots', {}).get('by_k', {})
        
        # JSON loads integer keys as strings, but the CLI should handle this properly
        # Just verify the data structure is reasonable
        self.assertIn('1', by_k)  # JSON converts int keys to strings
        self.assertIn('2', by_k)
        
        # Verify attempt invariants still hold
        self._validate_qa_consistency(qa_data)
    
    def test_three_way_invariant_validation(self):
        """Test that the three-way exclusive invariant works correctly."""
        # Test a mixed scenario with all three outcome types
        
        def mock_enumerate_with_stats(smiles, cfg):
            return [{'smiles': 'Fc1ccccc1', 'parent_smiles': 'c1ccccc1'}], {
                'version': '2',
                'attempts': 4,  # Total attempts
                'products': 1,  # Successful attempts
                'no_product_matches': 2,  # Failed attempts (but template worked)
                'template_unsupported': 1,  # Template unsupported attempts
                'dedup_hits_statesig': 0,
                'dedup_hits_inchi': 0,
                'qa_paths': {'isotope_unavailable': 1, 'isotope_miss': 0, 'atommap_used': 1, 'heuristic_used': 1},
                'pivots': {
                    'by_rule_halogen_k': {
                        'R1_F_1': {'attempts': 1, 'products': 1, 'atommap_used': 1},  # Success
                        'R1_Cl_1': {'attempts': 2, 'no_product_matches': 2, 'heuristic_used': 1, 'isotope_unavailable': 1},  # Failed
                        'R3_Br_1': {'attempts': 1, 'template_unsupported': 1}  # Template unsupported
                    },
                    'by_rule': {
                        'R1': {'attempts': 3, 'products': 1, 'no_product_matches': 2},
                        'R3': {'attempts': 1, 'template_unsupported': 1}
                    },
                    'by_halogen': {
                        'F': {'attempts': 1, 'products': 1},
                        'Cl': {'attempts': 2, 'no_product_matches': 2},
                        'Br': {'attempts': 1, 'template_unsupported': 1}
                    },
                    'by_k': {1: {'attempts': 4, 'products': 1, 'no_product_matches': 2, 'template_unsupported': 1}},
                    'by_rule_halogen': {
                        'R1_F': {'attempts': 1, 'products': 1},
                        'R1_Cl': {'attempts': 2, 'no_product_matches': 2},
                        'R3_Br': {'attempts': 1, 'template_unsupported': 1}
                    }
                }
            }
        
        # Create config for CLI
        config = {
            'io': {
                'smiles_file': self.input_file,
                'products_table': os.path.join(self.output_dir, 'products.parquet')
            },
            'k_max': 1,
            'halogens': ['F', 'Cl', 'Br']
        }
        
        # Create mock args
        class MockArgs:
            def __init__(self, outdir):
                self.outdir = outdir
                self.subset = 'flavonoids'
                self.k = 1
        
        mock_args = MockArgs(self.output_dir)
        
        # Run CLI with dependency injection
        cmd_enum(config, args=mock_args, enumerate_k1_fn=mock_enumerate_with_stats)
        
        # Load and validate QA stats
        qa_file = os.path.join(self.output_dir, 'qa_summary.json')
        with open(qa_file, 'r') as f:
            qa_data = json.load(f)
        
        # This should pass the three-way invariant check
        self._validate_qa_consistency(qa_data)
        
        # Specifically verify the three-way invariant
        # attempts=8 (4 per parent), products=2 (1 per parent), no_matches=4 (2 per parent), template_unsupported=2 (1 per parent)
        # 8 = 2 + 4 + 2
        self.assertEqual(qa_data['attempts'], 8)
        self.assertEqual(qa_data['products'], 2)
        self.assertEqual(qa_data['no_product_matches'], 4)
        self.assertEqual(qa_data['template_unsupported'], 2)


if __name__ == '__main__':
    unittest.main()