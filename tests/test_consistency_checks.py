# -*- coding: ascii -*-
"""Test consistency check validation logic for rule_halogen_k counting."""

import unittest
from collections import defaultdict, Counter
from src.halogenator.report import _validate_rule_halogen_k_consistency, _init_incremental_stats, _update_incremental_stats


class TestConsistencyChecks(unittest.TestCase):
    """Test rule_halogen_k consistency validation with independent cross-validation."""
    
    def setUp(self):
        """Set up test fixtures with sample data."""
        # Sample consistent statistics
        self.consistent_stats = {
            'rule_halogen_k_counts': {
                'R1': {
                    'F': {1: 10, 2: 5},
                    'Cl': {1: 8, 2: 3}
                },
                'R2': {
                    'F': {1: 15, 2: 7},
                    'Br': {1: 6, 2: 2}
                }
            },
            'rule_halogen_k_by_parent': {
                'parent1': {
                    'R1': {'F': {1: 5, 2: 2}, 'Cl': {1: 4, 2: 1}},
                    'R2': {'F': {1: 7, 2: 3}}
                },
                'parent2': {
                    'R1': {'F': {1: 5, 2: 3}, 'Cl': {1: 4, 2: 2}},
                    'R2': {'F': {1: 8, 2: 4}, 'Br': {1: 6, 2: 2}}
                }
            },
            # Independent sources for cross-validation
            'attempts_total': 56,  # Sum of all counts above
            'rule_halogen_counts': {
                'R1': {'F': 15, 'Cl': 11},  # Sum over k
                'R2': {'F': 22, 'Br': 8}
            },
            'rule_totals': {
                'R1': 26,  # 15 + 11
                'R2': 30   # 22 + 8
            },
            'halogen_counts': {
                'F': 37,   # 15 + 22
                'Cl': 11,  # 11
                'Br': 8    # 8
            },
            'k_counts': {
                1: 39,     # 10 + 8 + 15 + 6
                2: 17      # 5 + 3 + 7 + 2
            },
            'rule_counts': {
                'R1': 26,  # Should equal rule_totals
                'R2': 30
            },
            'product_count': 56  # Total products count - should match 3D sum
        }
    
    def test_consistent_data_all_checks_pass(self):
        """Test that consistent data passes all 6 validation checks."""
        result = _validate_rule_halogen_k_consistency(self.consistent_stats)
        
        # Overall consistency should be True
        self.assertTrue(result['consistent'], f"Expected overall consistency but got: {result}")
        
        # All individual checks should pass
        expected_checks = [
            'halogen_vs_flat',
            'rule_halogen_vs_flat', 
            'k_vs_flat',
            'rule_vs_independent',
            'global_product_conservation',
            'by_parent_vs_flat'
        ]
        
        for check_name in expected_checks:
            self.assertIn(check_name, result['checks'], f"Missing check: {check_name}")
            self.assertTrue(result['checks'][check_name]['consistent'], 
                          f"Check {check_name} should pass but failed: {result['checks'][check_name]}")
            self.assertEqual(len(result['checks'][check_name]['mismatches']), 0,
                           f"Check {check_name} should have no mismatches")
        
        # No mismatches
        self.assertEqual(result['mismatch_count'], 0)
        self.assertEqual(len(result['mismatches_sample']), 0)
    
    def test_halogen_dimension_mismatch(self):
        """Test detection of halogen dimension inconsistency."""
        # Corrupt halogen_counts to create mismatch
        corrupt_stats = self.consistent_stats.copy()
        corrupt_stats['halogen_counts'] = corrupt_stats['halogen_counts'].copy()
        corrupt_stats['halogen_counts']['F'] = 30  # Should be 37
        
        result = _validate_rule_halogen_k_consistency(corrupt_stats)
        
        # Overall consistency should fail
        self.assertFalse(result['consistent'])
        
        # Halogen check should fail specifically
        self.assertFalse(result['checks']['halogen_vs_flat']['consistent'])
        self.assertIn('F', result['checks']['halogen_vs_flat']['mismatches'])
        
        # Should detect the specific mismatch
        f_mismatch = result['checks']['halogen_vs_flat']['mismatches']['F']
        self.assertEqual(f_mismatch['calculated'], 37)  # From 3D sum
        self.assertEqual(f_mismatch['actual'], 30)      # Corrupted value
        self.assertEqual(f_mismatch['difference'], 7)
        
        # Should appear in mismatch sample
        self.assertGreater(result['mismatch_count'], 0)
        self.assertTrue(any('H[F]' in msg for msg in result['mismatches_sample']))
    
    def test_rule_halogen_dimension_mismatch(self):
        """Test detection of (rule, halogen) dimension inconsistency."""
        # Corrupt independent rule_halogen_counts
        corrupt_stats = self.consistent_stats.copy()
        corrupt_stats['rule_halogen_counts'] = {
            'R1': {'F': 20, 'Cl': 11},  # F should be 15, not 20
            'R2': {'F': 22, 'Br': 8}
        }
        
        result = _validate_rule_halogen_k_consistency(corrupt_stats)
        
        # Overall consistency should fail
        self.assertFalse(result['consistent'])
        
        # Rule-halogen check should fail
        self.assertFalse(result['checks']['rule_halogen_vs_flat']['consistent'])
        self.assertIn('R1_F', result['checks']['rule_halogen_vs_flat']['mismatches'])
        
        # Should detect the specific mismatch
        r1_f_mismatch = result['checks']['rule_halogen_vs_flat']['mismatches']['R1_F']
        self.assertEqual(r1_f_mismatch['calculated_from_3d'], 15)  # Correct 3D sum
        self.assertEqual(r1_f_mismatch['expected_from_2d'], 20)   # Corrupted independent value
        self.assertEqual(r1_f_mismatch['difference'], -5)
        
        # Should appear in mismatch sample
        self.assertTrue(any('RH[R1,F]' in msg for msg in result['mismatches_sample']))
    
    def test_k_dimension_mismatch(self):
        """Test detection of k dimension inconsistency."""
        # Corrupt k_counts
        corrupt_stats = self.consistent_stats.copy()
        corrupt_stats['k_counts'] = corrupt_stats['k_counts'].copy()
        corrupt_stats['k_counts'][1] = 35  # Should be 39
        
        result = _validate_rule_halogen_k_consistency(corrupt_stats)
        
        # Overall consistency should fail
        self.assertFalse(result['consistent'])
        
        # K dimension check should fail
        self.assertFalse(result['checks']['k_vs_flat']['consistent'])
        self.assertIn('1', result['checks']['k_vs_flat']['mismatches'])
        
        # Should detect the specific mismatch
        k1_mismatch = result['checks']['k_vs_flat']['mismatches']['1']
        self.assertEqual(k1_mismatch['calculated'], 39)  # From 3D sum
        self.assertEqual(k1_mismatch['actual'], 35)      # Corrupted value
        self.assertEqual(k1_mismatch['difference'], 4)
        
        # Should appear in mismatch sample
        self.assertTrue(any('K[1]' in msg for msg in result['mismatches_sample']))
    
    def test_rule_dimension_mismatch(self):
        """Test detection of rule dimension inconsistency using independent rule_totals."""
        # Corrupt independent rule_totals
        corrupt_stats = self.consistent_stats.copy()
        corrupt_stats['rule_totals'] = corrupt_stats['rule_totals'].copy()
        corrupt_stats['rule_totals']['R1'] = 30  # Should be 26
        
        result = _validate_rule_halogen_k_consistency(corrupt_stats)
        
        # Overall consistency should fail
        self.assertFalse(result['consistent'])
        
        # Rule dimension check should fail
        self.assertFalse(result['checks']['rule_vs_independent']['consistent'])
        self.assertIn('R1', result['checks']['rule_vs_independent']['mismatches'])
        
        # Should detect the specific mismatch
        r1_mismatch = result['checks']['rule_vs_independent']['mismatches']['R1']
        self.assertEqual(r1_mismatch['calculated_from_3d'], 26)        # From 3D sum
        self.assertEqual(r1_mismatch['expected_from_independent'], 30) # Corrupted independent
        self.assertEqual(r1_mismatch['difference'], -4)
        
        # Should appear in mismatch sample
        self.assertTrue(any('RULE[R1]' in msg for msg in result['mismatches_sample']))
    
    def test_global_total_mismatch(self):
        """Test detection of global total inconsistency using product_count."""
        # Corrupt product_count
        corrupt_stats = self.consistent_stats.copy()
        corrupt_stats['product_count'] = 50  # Should be 56 (3D sum)
        
        result = _validate_rule_halogen_k_consistency(corrupt_stats)
        
        # Overall consistency should fail
        self.assertFalse(result['consistent'])
        
        # Global check should fail
        self.assertFalse(result['checks']['global_product_conservation']['consistent'])
        self.assertIn('global_products', result['checks']['global_product_conservation']['mismatches'])
        
        # Should detect the specific mismatch
        global_mismatch = result['checks']['global_product_conservation']['mismatches']['global_products']
        self.assertEqual(global_mismatch['calculated_from_3d'], 56)         # From 3D sum
        self.assertEqual(global_mismatch['expected_product_count'], 50)     # Corrupted product_count
        self.assertEqual(global_mismatch['difference'], 6)
        
        # Should appear in mismatch sample
        self.assertTrue(any('PRODUCT_TOTAL:' in msg and 'product_count' in msg for msg in result['mismatches_sample']))
    
    def test_parent_aggregation_mismatch(self):
        """Test detection of parent aggregation inconsistency."""
        # Corrupt parent dimension - remove one count
        corrupt_stats = self.consistent_stats.copy()
        corrupt_stats['rule_halogen_k_by_parent'] = {
            'parent1': {
                'R1': {'F': {1: 5, 2: 2}, 'Cl': {1: 4, 2: 1}},
                'R2': {'F': {1: 7, 2: 3}}
            },
            'parent2': {
                'R1': {'F': {1: 5, 2: 3}, 'Cl': {1: 4, 2: 2}},
                'R2': {'F': {1: 8, 2: 4}, 'Br': {1: 6, 2: 1}}  # Changed 2 to 1 (missing 1 count)
            }
        }
        
        result = _validate_rule_halogen_k_consistency(corrupt_stats)
        
        # Overall consistency should fail
        self.assertFalse(result['consistent'])
        
        # Parent aggregation check should fail
        self.assertFalse(result['checks']['by_parent_vs_flat']['consistent'])
        self.assertIn('R2_Br_2', result['checks']['by_parent_vs_flat']['mismatches'])
        
        # Should detect the specific mismatch
        r2_br_2_mismatch = result['checks']['by_parent_vs_flat']['mismatches']['R2_Br_2']
        self.assertEqual(r2_br_2_mismatch['calculated_from_parents'], 1)  # Corrupted parent sum
        self.assertEqual(r2_br_2_mismatch['actual_flat'], 2)             # Correct 3D value
        self.assertEqual(r2_br_2_mismatch['difference'], -1)
        
        # Should appear in mismatch sample
        self.assertTrue(any('PARENT[R2,Br,2]' in msg for msg in result['mismatches_sample']))
    
    def test_multiple_mismatches_limited_sample(self):
        """Test that multiple mismatches are detected and sample is limited."""
        # Create corrupt stats with multiple issues
        corrupt_stats = self.consistent_stats.copy()
        corrupt_stats['halogen_counts'] = {'F': 30, 'Cl': 15, 'Br': 10}  # All wrong
        corrupt_stats['k_counts'] = {1: 35, 2: 20}                       # All wrong
        corrupt_stats['attempts_total'] = 50                             # Wrong
        
        result = _validate_rule_halogen_k_consistency(corrupt_stats)
        
        # Should fail overall
        self.assertFalse(result['consistent'])
        
        # Multiple checks should fail
        failed_checks = [name for name, check in result['checks'].items() 
                        if not check['consistent']]
        self.assertGreater(len(failed_checks), 1, "Should have multiple failed checks")
        
        # Should have multiple mismatches but sample limited to <=50
        self.assertGreater(result['mismatch_count'], 0)
        self.assertLessEqual(len(result['mismatches_sample']), 50, 
                           "Mismatch sample should be limited to 50 entries")
    
    def test_empty_stats_graceful_handling(self):
        """Test graceful handling of empty/missing statistics."""
        empty_stats = {}
        
        result = _validate_rule_halogen_k_consistency(empty_stats)
        
        # Should complete without error
        self.assertIsInstance(result, dict)
        self.assertIn('consistent', result)
        self.assertIn('checks', result)
        
        # With no data, most checks should pass trivially
        # (empty sets are consistent with each other)
        self.assertTrue(result['consistent'])
    
    def test_incremental_stats_integration(self):
        """Test that incremental stats collection produces consistent results."""
        # Create sample records
        sample_records = [
            {'rule': 'R1', 'halogen': 'F', 'k': 1, 'constraints_ok': True, 'sanitize_ok': True, 'parent_smiles': 'CCO', 'parent_inchikey': 'TEST1'},
            {'rule': 'R1', 'halogen': 'F', 'k': 2, 'constraints_ok': True, 'sanitize_ok': True, 'parent_smiles': 'CCO', 'parent_inchikey': 'TEST1'},
            {'rule': 'R1', 'halogen': 'Cl', 'k': 1, 'constraints_ok': True, 'sanitize_ok': True, 'parent_smiles': 'CCC', 'parent_inchikey': 'TEST2'},
            {'rule': 'R2', 'halogen': 'F', 'k': 1, 'constraints_ok': True, 'sanitize_ok': True, 'parent_smiles': 'CCO', 'parent_inchikey': 'TEST1'},
        ]
        
        # Process using incremental stats
        stats = _init_incremental_stats()
        for record in sample_records:
            _update_incremental_stats(stats, record)
        
        # Add mock QA stats to simulate enumeration results
        mock_qa_stats = {
            'attempts': 4,
            'products': 4,
            'no_product_matches': 0,
            'template_unsupported': 0
        }
        stats['enumeration_qa_stats'] = mock_qa_stats
        stats['attempts_total'] = mock_qa_stats['attempts']
        
        # Convert to serializable format
        from collections import defaultdict, Counter
        stats['k_counts'] = dict(stats['k_counts'])
        stats['rule_counts'] = dict(stats['rule_counts'])
        stats['halogen_counts'] = dict(stats['halogen_counts'])
        stats['rule_totals'] = dict(stats['rule_totals'])
        stats['rule_halogen_counts'] = {r: dict(h) for r, h in stats['rule_halogen_counts'].items()}
        stats['rule_halogen_k_counts'] = {r: {h: dict(k) for h, k in h_dict.items()} for r, h_dict in stats['rule_halogen_k_counts'].items()}
        stats['rule_halogen_k_by_parent'] = {p: {r: {h: dict(k) for h, k in h_dict.items()} for r, h_dict in r_dict.items()} for p, r_dict in stats['rule_halogen_k_by_parent'].items()}
        
        # Validate consistency
        result = _validate_rule_halogen_k_consistency(stats)
        
        # Should be consistent
        self.assertTrue(result['consistent'], 
                       f"Incremental stats should be consistent but got: {result}")
        
        # Verify expected counts
        self.assertEqual(stats['attempts_total'], 4)
        self.assertEqual(stats['rule_totals']['R1'], 3)
        self.assertEqual(stats['rule_totals']['R2'], 1)
        self.assertEqual(stats['halogen_counts']['F'], 3)
        self.assertEqual(stats['halogen_counts']['Cl'], 1)

    def test_three_way_mutex_invariant_violation(self):
        """Test detection of three-way mutex invariant violation."""
        # Create stats with enumeration QA stats that violate the invariant
        stats = {
            'rule_halogen_k_counts': {
                'R1': {'F': {1: 2}, 'Cl': {1: 1}}
            },
            'halogen_counts': {'F': 2, 'Cl': 1},
            'k_counts': {1: 3},
            'rule_counts': {'R1': 3},
            'rule_halogen_k_by_parent': {
                'parent1': {'R1': {'F': {1: 2}}},
                'parent2': {'R1': {'Cl': {1: 1}}}
            },
            'attempts_total': 10,
            'rule_halogen_counts': {'R1': {'F': 2, 'Cl': 1}},
            'rule_totals': {'R1': 3},
            # QA stats that violate the three-way mutex invariant
            'enumeration_qa_stats': {
                'attempts': 10,
                'products': 3,
                'no_product_matches': 2,  # 3 + 2 + 4 = 9 != 10 (violation)
                'template_unsupported': 4
            }
        }
        
        result = _validate_rule_halogen_k_consistency(stats)
        
        # Should detect the three-way mutex violation
        self.assertFalse(result['consistent'],
                        "Should detect three-way mutex invariant violation")
        
        # Should have three-way mutex check
        self.assertIn('three_way_mutex', result['checks'])
        self.assertFalse(result['checks']['three_way_mutex']['consistent'])
        
        # Should have specific mismatch info
        mutex_mismatches = result['checks']['three_way_mutex']['mismatches']
        self.assertEqual(mutex_mismatches['attempts_total'], 10)
        self.assertEqual(mutex_mismatches['products'], 3)
        self.assertEqual(mutex_mismatches['no_product_matches'], 2)
        self.assertEqual(mutex_mismatches['template_unsupported'], 4)
        self.assertEqual(mutex_mismatches['calculated_sum'], 9)
        self.assertEqual(mutex_mismatches['difference'], 1)  # 10 - 9

    def test_three_way_mutex_invariant_passes(self):
        """Test that three-way mutex invariant passes when data is correct."""
        # Create stats with enumeration QA stats that satisfy the invariant
        stats = {
            'rule_halogen_k_counts': {
                'R1': {'F': {1: 2}, 'Cl': {1: 1}}
            },
            'halogen_counts': {'F': 2, 'Cl': 1},
            'k_counts': {1: 3},
            'rule_counts': {'R1': 3},
            'rule_halogen_k_by_parent': {
                'parent1': {'R1': {'F': {1: 2}}},
                'parent2': {'R1': {'Cl': {1: 1}}}
            },
            'attempts_total': 10,
            'rule_halogen_counts': {'R1': {'F': 2, 'Cl': 1}},
            'rule_totals': {'R1': 3},
            # QA stats that satisfy the three-way mutex invariant
            'enumeration_qa_stats': {
                'attempts': 10,
                'products': 3,
                'no_product_matches': 5,  # 3 + 5 + 2 = 10 (correct)
                'template_unsupported': 2
            }
        }
        
        result = _validate_rule_halogen_k_consistency(stats)
        
        # Should have three-way mutex check that passes
        self.assertIn('three_way_mutex', result['checks'])
        self.assertTrue(result['checks']['three_way_mutex']['consistent'])
        
        # Should have empty mismatches for this check
        mutex_mismatches = result['checks']['three_way_mutex']['mismatches']
        self.assertEqual(mutex_mismatches, {})

    def test_global_totals_use_products_not_attempts(self):
        """Test that global totals validation uses product counts, not attempts."""
        # Create stats where attempts > products but product conservation holds
        stats = {
            'rule_halogen_k_counts': {
                'R1': {'F': {1: 2}, 'Cl': {1: 1}}
            },
            'halogen_counts': {'F': 2, 'Cl': 1},  # sum = 3 (matches product_count)
            'k_counts': {1: 3},
            'rule_counts': {'R1': 3},
            'rule_halogen_k_by_parent': {
                'parent1': {'R1': {'F': {1: 2}}},
                'parent2': {'R1': {'Cl': {1: 1}}}
            },
            'product_count': 3,  # Matches sum of halogen_counts
            'attempts_total': 10,  # Much higher than products (has failures)
            'rule_halogen_counts': {'R1': {'F': 2, 'Cl': 1}},
            'rule_totals': {'R1': 3},
            # QA stats with failures (attempts > products)
            'enumeration_qa_stats': {
                'attempts': 10,
                'products': 3,
                'no_product_matches': 5,  
                'template_unsupported': 2  # 3 + 5 + 2 = 10 (mutex correct)
            }
        }
        
        result = _validate_rule_halogen_k_consistency(stats)
        
        # Global product conservation should pass (3D sum = product_count = 3)
        self.assertIn('global_product_conservation', result['checks'])
        global_check = result['checks']['global_product_conservation']
        self.assertTrue(global_check['consistent'], 
                       f"Global product conservation should pass but got: {global_check}")
        
        # Three-way mutex should also pass (separate check)
        self.assertIn('three_way_mutex', result['checks'])
        self.assertTrue(result['checks']['three_way_mutex']['consistent'])
        
        # Overall should be consistent (both checks pass)
        self.assertTrue(result['consistent'], 
                       f"Overall consistency should pass but got: {result}")

    def test_global_product_conservation_sub_buckets_present(self):
        """Ensure global_product_conservation exposes all three mismatch buckets."""
        stats = {
            'rule_halogen_k_counts': {
                'R1': {'F': {1: 2}, 'Cl': {1: 1}}
            },
            'rule_halogen_k_by_parent': {
                'parent1': {'R1': {'F': {1: 2}}},
                'parent2': {'R1': {'Cl': {1: 1}}}
            },
            # Halogen sum will be 5 (2 + 3)
            'halogen_counts': {'F': 2, 'Cl': 3},
            # Product count intentionally 4 (mismatch vs 3D sum=3 and halogen_sum=5)
            'product_count': 4,
            'k_counts': {1: 3},
            'rule_counts': {'R1': 3},
            'rule_halogen_counts': {'R1': {'F': 2, 'Cl': 1}},
            'rule_totals': {'R1': 3},
            # Enumeration stats with matching three-way mutex and products=6 to trigger product_count_vs_enum
            'enumeration_qa_stats': {
                'attempts': 10,
                'products': 6,
                'no_product_matches': 2,
                'template_unsupported': 2
            }
        }

        result = _validate_rule_halogen_k_consistency(stats)
        self.assertIn('global_product_conservation', result['checks'])
        mismatches = result['checks']['global_product_conservation']['mismatches']
        # Expect all three buckets
        self.assertIn('global_products', mismatches)
        self.assertIn('product_count_vs_enum', mismatches)
        self.assertIn('product_vs_halogen_sum', mismatches)


if __name__ == '__main__':
    unittest.main()
