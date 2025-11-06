# -*- coding: ascii -*-
"""Tests for granular QA JSON dimensions functionality."""

import unittest
import tempfile
import os
import json
import shutil
from unittest.mock import patch
from src.halogenator.cli import cmd_enum, _create_granular_qa_summary
from src.halogenator.enumerate_k import EnumConfig
from src.halogenator.io_utils import write_smi


class TestGranularQAJSON(unittest.TestCase):
    """Test that QA JSON includes granular dimensions."""
    
    def setUp(self):
        """Set up test environment."""
        self.test_dir = tempfile.mkdtemp()
        self.config = {
            'io': {
                'smiles_file': os.path.join(self.test_dir, 'parents.smi'),
                'products_table': os.path.join(self.test_dir, 'products.parquet')
            },
            'halogens': ['F', 'Cl'],
            'k_max': 1,
            'constraints': {'per_ring_quota': 2, 'min_graph_distance': 2},
            'standardize': {'do_tautomer': False},
            'qc': {'sanitize_strict': True},
            'pruning': {'enable_symmetry_fold': True}
        }
        
        # Create test parent molecules
        test_parents = [('c1ccc(O)cc1', 'phenol')]
        write_smi(test_parents, self.config['io']['smiles_file'])
    
    def tearDown(self):
        """Clean up test directory."""
        shutil.rmtree(self.test_dir, ignore_errors=True)
    
    def test_granular_qa_summary_creation(self):
        """Test that _create_granular_qa_summary creates proper structure."""
        # Mock QA stats
        total_qa_stats = {
            'isotope_unavailable': 8,
            'isotope_miss': 2,
            'atommap_used': 4,
            'heuristic_used': 1,
            'no_product_matches': 3,
            'template_unsupported': 0,
            'qa_paths': {
                'isotope_unavailable': 8,
                'isotope_miss': 2,
                'atommap_used': 4,
                'heuristic_used': 1
            }
        }
        
        enum_cfg = EnumConfig(
            k_max=1,
            halogens=('F', 'Cl'),
            constraints={'per_ring_quota': 2},
            std_cfg={'do_tautomer': False},
            qc_cfg={'sanitize_strict': True},
            pruning_cfg={'enable_symmetry_fold': True}
        )
        
        granular_summary = _create_granular_qa_summary(total_qa_stats, enum_cfg)
        
        # Verify structure
        self.assertIn('total', granular_summary)
        self.assertIn('by_rule', granular_summary)
        self.assertIn('by_halogen', granular_summary)
        self.assertIn('by_rule_halogen', granular_summary)
        self.assertIn('version', granular_summary)
        
        # Verify total stats are preserved
        self.assertEqual(granular_summary['total'], total_qa_stats)
        
        # Verify by_rule structure
        by_rule = granular_summary['by_rule']
        expected_rules = ['R1', 'R3', 'R4', 'R5']
        self.assertEqual(set(by_rule.keys()), set(expected_rules))
        
        for rule in expected_rules:
            rule_stats = by_rule[rule]
            self.assertIn('isotope_unavailable', rule_stats)
            self.assertIn('isotope_miss', rule_stats)
            self.assertIn('atommap_used', rule_stats)
            self.assertIn('heuristic_used', rule_stats)
            self.assertIn('no_product_matches', rule_stats)
            self.assertIn('template_unsupported', rule_stats)
            
            # Values should be distributed across rules
            self.assertEqual(rule_stats['isotope_unavailable'], 8 // 4)  # 8 / 4 rules = 2
            self.assertEqual(rule_stats['atommap_used'], 4 // 4)  # 4 / 4 rules = 1
        
        # Verify by_halogen structure
        by_halogen = granular_summary['by_halogen']
        expected_halogens = ['F', 'Cl']
        self.assertEqual(set(by_halogen.keys()), set(expected_halogens))
        
        for halogen in expected_halogens:
            halogen_stats = by_halogen[halogen]
            self.assertIn('isotope_unavailable', halogen_stats)
            # Values should be distributed across halogens
            self.assertEqual(halogen_stats['isotope_unavailable'], 8 // 2)  # 8 / 2 halogens = 4
        
        # Verify by_rule_halogen structure
        by_rule_halogen = granular_summary['by_rule_halogen']
        self.assertEqual(set(by_rule_halogen.keys()), set(expected_rules))
        
        for rule in expected_rules:
            self.assertEqual(set(by_rule_halogen[rule].keys()), set(expected_halogens))
            for halogen in expected_halogens:
                combo_stats = by_rule_halogen[rule][halogen]
                self.assertIn('isotope_unavailable', combo_stats)
                # Values should be distributed across rule x halogen combinations
                self.assertEqual(combo_stats['isotope_unavailable'], 8 // 8)  # 8 / (4 rules x 2 halogens) = 1
    
    def test_cli_generates_granular_qa_json(self):
        """Test that CLI enum generates granular QA JSON."""
        class MockArgs:
            def __init__(self, outdir):
                self.k = 1
                self.subset = 'flavonoids'
                self.outdir = outdir
        
        args = MockArgs(self.test_dir)
        
        # Run enumeration
        with patch('builtins.print'):
            cmd_enum(self.config, args)
        
        # Check that QA JSON was generated
        qa_json_path = os.path.join(self.test_dir, 'qa_summary.json')
        self.assertTrue(os.path.exists(qa_json_path))
        
        # Verify JSON structure
        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)
        
        # Should have version and metadata
        self.assertIn('version', qa_data)
        self.assertIn('metadata', qa_data)
        
        # Handle both v1 and v2 formats
        if qa_data.get('version') == '2':
            # v2 format has pivots instead of granular breakdowns
            self.assertIn('pivots', qa_data)
            self.assertIn('attempts', qa_data)  # v2 uses 'attempts' not 'total_attempts'
            # Verify metadata for v2 - basic structure is sufficient
            metadata = qa_data['metadata']
            self.assertIn('description', metadata)
        else:
            # v1 format has granular structure
            self.assertIn('total', qa_data)
            self.assertIn('by_rule', qa_data)
            self.assertIn('by_halogen', qa_data)
            self.assertIn('by_rule_halogen', qa_data)
            
            # Verify metadata for v1
            metadata = qa_data['metadata']
            self.assertIn('description', metadata)
            self.assertIn('semantics', metadata)
            self.assertIn('generated_at', metadata)
            
            # Version should be '1'
            self.assertEqual(qa_data['version'], '1')
    
    def test_granular_dimensions_sum_to_total(self):
        """Test that granular dimension stats sum to total where appropriate."""
        total_qa_stats = {
            'isotope_unavailable': 12,
            'isotope_miss': 8,
            'atommap_used': 6,
            'heuristic_used': 4,
            'no_product_matches': 2,
            'template_unsupported': 0,
            'qa_paths': {'isotope_unavailable': 12, 'isotope_miss': 8, 'atommap_used': 6, 'heuristic_used': 4}
        }
        
        enum_cfg = EnumConfig(
            k_max=1,
            halogens=('F', 'Cl', 'Br'),
            constraints={},
            std_cfg={},
            qc_cfg={},
            pruning_cfg={}
        )
        
        granular_summary = _create_granular_qa_summary(total_qa_stats, enum_cfg)
        
        # Check that by_rule sums approximately match totals (allowing for integer division)
        by_rule = granular_summary['by_rule']
        rule_isotope_sum = sum(rule_stats['isotope_unavailable'] for rule_stats in by_rule.values())
        # Should be close to total (within rounding error)
        self.assertLessEqual(abs(rule_isotope_sum - total_qa_stats['isotope_unavailable']), 4)  # 4 rules max error
        
        # Check that by_halogen sums approximately match totals
        by_halogen = granular_summary['by_halogen']
        halogen_isotope_sum = sum(hal_stats['isotope_unavailable'] for hal_stats in by_halogen.values())
        self.assertLessEqual(abs(halogen_isotope_sum - total_qa_stats['isotope_unavailable']), 3)  # 3 halogens max error
        
        # Check that by_rule_halogen sums approximately match totals
        by_rule_halogen = granular_summary['by_rule_halogen']
        rule_halogen_isotope_sum = 0
        for rule_stats in by_rule_halogen.values():
            for hal_stats in rule_stats.values():
                rule_halogen_isotope_sum += hal_stats['isotope_unavailable']
        self.assertLessEqual(abs(rule_halogen_isotope_sum - total_qa_stats['isotope_unavailable']), 12)  # 12 combinations max error
    
    def test_empty_stats_handled_gracefully(self):
        """Test that empty/zero stats are handled gracefully."""
        empty_qa_stats = {
            'isotope_unavailable': 0,
            'isotope_miss': 0,
            'atommap_used': 0,
            'heuristic_used': 0,
            'no_product_matches': 0,
            'template_unsupported': 0,
            'qa_paths': {}
        }
        
        enum_cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={},
            std_cfg={},
            qc_cfg={},
            pruning_cfg={}
        )
        
        granular_summary = _create_granular_qa_summary(empty_qa_stats, enum_cfg)
        
        # Should not crash and should have proper structure
        self.assertIn('total', granular_summary)
        self.assertIn('by_rule', granular_summary)
        self.assertIn('by_halogen', granular_summary)
        self.assertIn('by_rule_halogen', granular_summary)
        
        # All breakdown values should be 0
        for rule_stats in granular_summary['by_rule'].values():
            for key, value in rule_stats.items():
                self.assertEqual(value, 0)
        
        for hal_stats in granular_summary['by_halogen'].values():
            for key, value in hal_stats.items():
                self.assertEqual(value, 0)
    
    def test_single_halogen_distribution(self):
        """Test granular breakdown with single halogen."""
        total_qa_stats = {
            'isotope_unavailable': 4,
            'isotope_miss': 0,
            'atommap_used': 2,
            'heuristic_used': 1,
            'no_product_matches': 0,
            'template_unsupported': 0,
            'qa_paths': {'isotope_unavailable': 4, 'atommap_used': 2, 'heuristic_used': 1, 'isotope_miss': 0}
        }
        
        enum_cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={},
            std_cfg={},
            qc_cfg={},
            pruning_cfg={}
        )
        
        granular_summary = _create_granular_qa_summary(total_qa_stats, enum_cfg)
        
        # by_halogen should have only F
        by_halogen = granular_summary['by_halogen']
        self.assertEqual(set(by_halogen.keys()), {'F'})
        
        # F should get all the stats
        f_stats = by_halogen['F']
        self.assertEqual(f_stats['isotope_unavailable'], 4)
        self.assertEqual(f_stats['atommap_used'], 2)
        self.assertEqual(f_stats['heuristic_used'], 1)
        
        # by_rule_halogen should have F under each rule
        by_rule_halogen = granular_summary['by_rule_halogen']
        for rule in ['R1', 'R3', 'R4', 'R5']:
            self.assertEqual(set(by_rule_halogen[rule].keys()), {'F'})
            # Each rule x halogen combination gets 1/4 of total
            self.assertEqual(by_rule_halogen[rule]['F']['isotope_unavailable'], 4 // 4)
    
    def test_qa_json_backward_compatibility(self):
        """Test that write_qa_summary_json handles both old and new formats."""
        from src.halogenator.report import write_qa_summary_json
        
        # Test old format (no version)
        old_format_stats = {
            'isotope_unavailable': 5,
            'isotope_miss': 2,
            'no_product_matches': 1
        }
        
        old_json_path = write_qa_summary_json(old_format_stats, self.test_dir)
        
        with open(old_json_path, 'r') as f:
            old_qa_data = json.load(f)
        
        # Should have legacy structure
        self.assertIn('total', old_qa_data)
        self.assertIn('metadata', old_qa_data)
        self.assertNotIn('by_rule', old_qa_data)
        self.assertEqual(old_qa_data['total'], old_format_stats)
        
        # Test new format (with version)
        new_format_stats = {
            'version': '1',
            'total': {'isotope_unavailable': 5, 'isotope_miss': 2},
            'by_rule': {'R3': {'isotope_unavailable': 2}},
            'by_halogen': {'F': {'isotope_unavailable': 3}}
        }
        
        new_json_path = write_qa_summary_json(new_format_stats, self.test_dir)
        
        with open(new_json_path, 'r') as f:
            new_qa_data = json.load(f)
        
        # Should preserve granular structure
        self.assertIn('total', new_qa_data)
        self.assertIn('by_rule', new_qa_data)
        self.assertIn('by_halogen', new_qa_data)
        self.assertIn('metadata', new_qa_data)
        self.assertEqual(new_qa_data['version'], '1')

    def test_v1_rules_default_contains_R2(self):
        """Test that default rules include R2 when metadata.rules is missing."""
        from src.halogenator.report import write_qa_summary_json
        
        # Create v1 input without metadata.rules but with metadata.halogens
        qa_stats_without_rules = {
            'isotope_unavailable': 8,
            'no_product_matches': 2,
            'metadata': {
                'halogens': ['F', 'Cl']
                # Intentionally omit 'rules' to test default
            }
        }
        
        qa_json_path = write_qa_summary_json(qa_stats_without_rules, self.test_dir)
        
        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)
        
        # Should have granular structure with default rules including R2
        self.assertIn('by_rule', qa_data)
        by_rule_keys = set(qa_data['by_rule'].keys())
        expected_rules = {'R1', 'R2', 'R3', 'R4', 'R5'}
        self.assertEqual(by_rule_keys, expected_rules)
        self.assertIn('R2', by_rule_keys)

    def test_v1_distribution_remainder_kept(self):
        """Test that remainder distribution preserves totals exactly."""
        from src.halogenator.report import write_qa_summary_json
        
        # Create input with values that don't divide evenly
        qa_stats_with_remainder = {
            'isotope_unavailable': 5,  # 5 / 2 rules = 2 remainder 1
            'atommap_used': 7,        # 7 / 2 halogens = 3 remainder 1  
            'no_product_matches': 11, # 11 / (2 rules x 2 halogens) = 2 remainder 3
            'metadata': {
                'halogens': ['F', 'Cl'],
                'rules': ['R1', 'R3']
            }
        }
        
        qa_json_path = write_qa_summary_json(qa_stats_with_remainder, self.test_dir)
        
        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)
        
        # Check by_rule sums exactly match totals
        by_rule = qa_data['by_rule']
        rule_isotope_sum = sum(rule_stats['isotope_unavailable'] for rule_stats in by_rule.values())
        self.assertEqual(rule_isotope_sum, 5)
        
        # Check by_halogen sums exactly match totals
        by_halogen = qa_data['by_halogen']
        halogen_atommap_sum = sum(hal_stats['atommap_used'] for hal_stats in by_halogen.values())
        self.assertEqual(halogen_atommap_sum, 7)
        
        # Check by_rule_halogen sums exactly match totals
        by_rule_halogen = qa_data['by_rule_halogen']
        rule_halogen_no_product_sum = 0
        for rule_stats in by_rule_halogen.values():
            for hal_stats in rule_stats.values():
                rule_halogen_no_product_sum += hal_stats['no_product_matches']
        self.assertEqual(rule_halogen_no_product_sum, 11)

    def test_cli_injects_rules_from_config(self):
        """Test that CLI injects actual EnumConfig.rules into metadata."""
        class MockArgs:
            def __init__(self, outdir):
                self.k = 1
                self.subset = 'flavonoids'
                self.outdir = outdir
        
        # Test with custom rules (only R1, R2)
        custom_config = self.config.copy()
        custom_config['rules'] = ['R1', 'R2']  # Custom subset
        
        args = MockArgs(self.test_dir)
        
        # Run enumeration
        with patch('builtins.print'):
            cmd_enum(custom_config, args)
        
        # Check that QA JSON was generated
        qa_json_path = os.path.join(self.test_dir, 'qa_summary.json')
        self.assertTrue(os.path.exists(qa_json_path))
        
        # Verify JSON structure uses custom rules
        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)
        
        if qa_data.get('version') == '1' and 'by_rule' in qa_data:
            # Should only have R1 and R2, not the full set
            by_rule_keys = set(qa_data['by_rule'].keys())
            expected_custom_rules = {'R1', 'R2'}
            self.assertEqual(by_rule_keys, expected_custom_rules)
            
            # by_rule_halogen should also match
            by_rule_halogen_keys = set(qa_data['by_rule_halogen'].keys())
            self.assertEqual(by_rule_halogen_keys, expected_custom_rules)

    def test_v1_granular_passthrough_shape_completion(self):
        """Test that incomplete v1 granular inputs get missing slices completed."""
        from src.halogenator.report import write_qa_summary_json
        
        # Create v1 input with only by_rule (missing by_halogen, by_rule_halogen)
        incomplete_v1_input = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': 2, 'atommap_used': 1},
                'R3': {'isotope_unavailable': 1, 'atommap_used': 2}
            },
            'metadata': {
                'halogens': ['F', 'Cl'],
                'rules': ['R1', 'R3']
            }
        }
        
        qa_json_path = write_qa_summary_json(incomplete_v1_input, self.test_dir)
        
        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)
        
        # Should have completed all granular slices including total
        self.assertIn('by_rule', qa_data)
        self.assertIn('by_halogen', qa_data)
        self.assertIn('by_rule_halogen', qa_data)
        self.assertIn('total', qa_data)  # total should be guaranteed
        
        # Check that total is computed from by_rule (default zero-fill mode)
        total = qa_data['total']
        self.assertEqual(total['isotope_unavailable'], 3)  # 2 + 1 from by_rule
        self.assertEqual(total['atommap_used'], 3)  # 1 + 2 from by_rule

        # In zero-fill mode, 2D is zero-filled, so total != 2D aggregation
        by_rule_halogen = qa_data['by_rule_halogen']
        total_from_2d = sum(by_rule_halogen[rule][halogen]['isotope_unavailable']
                           for rule in ['R1', 'R3'] for halogen in ['F', 'Cl'])
        self.assertEqual(total_from_2d, 0)  # 2D is zero-filled, so should be 0
        
        # Default zero-fill mode: missing dimensions should be zero-filled
        by_halogen = qa_data['by_halogen']
        self.assertEqual(set(by_halogen.keys()), {'F', 'Cl'})
        for halogen in ['F', 'Cl']:
            self.assertIn('isotope_unavailable', by_halogen[halogen])
            self.assertIn('atommap_used', by_halogen[halogen])
            # Should be 0 (default zero-fill behavior for missing structure)
            self.assertEqual(by_halogen[halogen]['isotope_unavailable'], 0)
            self.assertEqual(by_halogen[halogen]['atommap_used'], 0)

        # Check by_rule_halogen structure was properly zero-filled
        by_rule_halogen = qa_data['by_rule_halogen']
        self.assertEqual(set(by_rule_halogen.keys()), {'R1', 'R3'})
        for rule in ['R1', 'R3']:
            self.assertEqual(set(by_rule_halogen[rule].keys()), {'F', 'Cl'})
            for halogen in ['F', 'Cl']:
                self.assertIn('isotope_unavailable', by_rule_halogen[rule][halogen])
                self.assertIn('atommap_used', by_rule_halogen[rule][halogen])
                # Should be 0 (default zero-fill behavior for missing 2D structure)
                self.assertEqual(by_rule_halogen[rule][halogen]['isotope_unavailable'], 0)
                self.assertEqual(by_rule_halogen[rule][halogen]['atommap_used'], 0)

    def test_v1_granular_passthrough_robust_metrics_inference(self):
        """Test that v1 granular passthrough uses union of all rule metrics."""
        from src.halogenator.report import write_qa_summary_json
        
        # Create v1 input where different rules have different metric sets
        incomplete_v1_with_diff_metrics = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': 2, 'atommap_used': 1},
                'R3': {'isotope_unavailable': 1, 'rdkit_error': 3}  # rdkit_error only in R3
            },
            'metadata': {
                'halogens': ['F', 'Cl'],
                'rules': ['R1', 'R3']
            }
        }
        
        qa_json_path = write_qa_summary_json(incomplete_v1_with_diff_metrics, self.test_dir)
        
        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)
        
        # Should have completed all granular slices with union of metrics
        self.assertIn('by_rule', qa_data)
        self.assertIn('by_halogen', qa_data)
        self.assertIn('by_rule_halogen', qa_data)
        
        # Check by_halogen contains all metrics from union
        by_halogen = qa_data['by_halogen']
        for halogen in ['F', 'Cl']:
            self.assertIn('isotope_unavailable', by_halogen[halogen])  # from both rules
            self.assertIn('atommap_used', by_halogen[halogen])         # from R1 only originally  
            self.assertIn('rdkit_error', by_halogen[halogen])          # from R3 only originally
            self.assertIn('no_product_matches', by_halogen[halogen])   # standard top-level metric
            
        # Check by_rule_halogen contains all metrics from union and distributed values
        by_rule_halogen = qa_data['by_rule_halogen']
        for rule in ['R1', 'R3']:
            for halogen in ['F', 'Cl']:
                self.assertIn('isotope_unavailable', by_rule_halogen[rule][halogen])
                self.assertIn('atommap_used', by_rule_halogen[rule][halogen])
                self.assertIn('rdkit_error', by_rule_halogen[rule][halogen])
                self.assertIn('no_product_matches', by_rule_halogen[rule][halogen])
                # no_product_matches should be 0 since missing from input
                self.assertEqual(by_rule_halogen[rule][halogen]['no_product_matches'], 0)

        # Default zero-fill mode: all missing 2D entries should be 0
        # R1 has rdkit_error=0 (missing) -> Cl:0, F:0 (zero-fill)
        # R3 has rdkit_error=3 -> but 2D structure is zero-filled since missing
        self.assertEqual(by_rule_halogen['R1']['Cl']['rdkit_error'], 0)
        self.assertEqual(by_rule_halogen['R1']['F']['rdkit_error'], 0)
        self.assertEqual(by_rule_halogen['R3']['Cl']['rdkit_error'], 0)
        self.assertEqual(by_rule_halogen['R3']['F']['rdkit_error'], 0)

    def test_v1_granular_marginal_consistency(self):
        """Test that by_rule and by_halogen are consistent with by_rule_halogen margins."""
        from src.halogenator.report import write_qa_summary_json
        
        # Create input with values that don't divide evenly to test marginal consistency
        qa_stats_marginal_test = {
            'isotope_unavailable': 17,  # 17 / (2 rules x 3 halogens) = 2 remainder 5
            'atommap_used': 13,         # 13 / (2 rules x 3 halogens) = 2 remainder 1
            'no_product_matches': 11,   # 11 / (2 rules x 3 halogens) = 1 remainder 5
            'metadata': {
                'halogens': ['F', 'Cl', 'Br'],
                'rules': ['R1', 'R3']
            }
        }
        
        qa_json_path = write_qa_summary_json(qa_stats_marginal_test, self.test_dir)
        
        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)
        
        by_rule = qa_data['by_rule']
        by_halogen = qa_data['by_halogen']
        by_rule_halogen = qa_data['by_rule_halogen']
        
        # Test marginal consistency for each metric
        for metric in ['isotope_unavailable', 'atommap_used', 'no_product_matches']:
            # Check rule margins: sum(by_rule_halogen[rule, *]) == by_rule[rule]
            for rule in ['R1', 'R3']:
                rule_sum_from_2d = sum(by_rule_halogen[rule][halogen][metric] for halogen in ['F', 'Cl', 'Br'])
                rule_value_from_margin = by_rule[rule][metric]
                self.assertEqual(rule_sum_from_2d, rule_value_from_margin, 
                               f"Rule {rule}, metric {metric}: 2D sum {rule_sum_from_2d} != margin {rule_value_from_margin}")
            
            # Check halogen margins: sum(by_rule_halogen[*, halogen]) == by_halogen[halogen]  
            for halogen in ['F', 'Cl', 'Br']:
                halogen_sum_from_2d = sum(by_rule_halogen[rule][halogen][metric] for rule in ['R1', 'R3'])
                halogen_value_from_margin = by_halogen[halogen][metric]
                self.assertEqual(halogen_sum_from_2d, halogen_value_from_margin,
                               f"Halogen {halogen}, metric {metric}: 2D sum {halogen_sum_from_2d} != margin {halogen_value_from_margin}")
            
            # Check total consistency: sum(by_rule_halogen[*, *]) == totals
            total_sum_from_2d = sum(by_rule_halogen[rule][halogen][metric] 
                                   for rule in ['R1', 'R3'] for halogen in ['F', 'Cl', 'Br'])
            expected_total = qa_stats_marginal_test[metric]
            self.assertEqual(total_sum_from_2d, expected_total,
                           f"Metric {metric}: 2D total sum {total_sum_from_2d} != expected {expected_total}")

    def test_v1_granular_metadata_transparency(self):
        """Test that v1 granular outputs include distribution semantics metadata."""
        from src.halogenator.report import write_qa_summary_json
        
        # Test constructed v1 granular (from legacy input)
        legacy_input = {
            'isotope_unavailable': 7,
            'no_product_matches': 3,
            'metadata': {
                'halogens': ['F', 'Cl'],
                'rules': ['R1', 'R3']
            }
        }
        
        qa_json_path = write_qa_summary_json(legacy_input, self.test_dir)
        
        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)
        
        # Check that v1 granular has distribution semantics metadata
        self.assertEqual(qa_data['version'], '1')
        self.assertIn('by_rule', qa_data)
        
        metadata = qa_data['metadata']
        self.assertIn('distribution', metadata)
        self.assertIn('marginals_from_2d', metadata)
        self.assertEqual(metadata['distribution'], "equal_with_lexicographic_remainder")
        self.assertEqual(metadata['marginals_from_2d'], True)
        
        # Test passthrough v1 granular (from existing granular input)
        passthrough_input = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': 2},
                'R3': {'isotope_unavailable': 1}
            },
            'metadata': {
                'halogens': ['F', 'Cl'],
                'rules': ['R1', 'R3']
            }
        }
        
        qa_json_path_2 = write_qa_summary_json(passthrough_input, self.test_dir)
        
        with open(qa_json_path_2, 'r') as f:
            qa_data_2 = json.load(f)
        
        # Check that passthrough v1 granular also has distribution semantics metadata
        self.assertEqual(qa_data_2['version'], '1')
        self.assertIn('by_rule', qa_data_2)

        metadata_2 = qa_data_2['metadata']
        self.assertIn('distribution', metadata_2)
        self.assertIn('marginals_from_2d', metadata_2)
        self.assertIn('marginals_source', metadata_2)
        self.assertEqual(metadata_2['distribution'], "equal_with_lexicographic_remainder")
        self.assertEqual(metadata_2['marginals_from_2d'], False)  # Passthrough: margins from input
        self.assertEqual(metadata_2['marginals_source'], "input")

    def test_v1_overview_vs_distributable_separation(self):
        """Test that overview counters and distributable metrics are properly separated."""
        from src.halogenator.report import write_qa_summary_json
        
        # Test input with both distributable metrics and overview counters
        mixed_input = {
            'isotope_unavailable': 6,  # Distributable QA path metric
            'no_product_matches': 4,   # Distributable top-level metric
            'attempts': 100,           # Overview counter (should not be distributed)
            'products': 50,            # Overview counter (should not be distributed)
            'dedup_hits_statesig': 5,  # Overview counter (should not be distributed)
            'metadata': {
                'halogens': ['F', 'Cl'],
                'rules': ['R1', 'R3']
            }
        }
        
        qa_json_path = write_qa_summary_json(mixed_input, self.test_dir)
        
        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)
        
        # Check that v1 granular was constructed
        self.assertEqual(qa_data['version'], '1')
        self.assertIn('by_rule', qa_data)
        self.assertIn('total', qa_data)
        
        total = qa_data['total']
        by_rule = qa_data['by_rule']
        by_halogen = qa_data['by_halogen']
        by_rule_halogen = qa_data['by_rule_halogen']
        
        # Overview counters should be preserved in total
        self.assertIn('attempts', total)
        self.assertIn('products', total)
        self.assertIn('dedup_hits_statesig', total)
        self.assertEqual(total['attempts'], 100)
        self.assertEqual(total['products'], 50)
        self.assertEqual(total['dedup_hits_statesig'], 5)
        
        # Distributable metrics should be present in total and granular dimensions
        self.assertIn('isotope_unavailable', total)
        self.assertIn('no_product_matches', total)
        
        # Overview counters should NOT appear in granular dimensions
        for rule in ['R1', 'R3']:
            self.assertNotIn('attempts', by_rule[rule])
            self.assertNotIn('products', by_rule[rule])
            self.assertNotIn('dedup_hits_statesig', by_rule[rule])
            
            # But distributable metrics should be present
            self.assertIn('isotope_unavailable', by_rule[rule])
            self.assertIn('no_product_matches', by_rule[rule])
        
        for halogen in ['F', 'Cl']:
            self.assertNotIn('attempts', by_halogen[halogen])
            self.assertNotIn('products', by_halogen[halogen])
            self.assertNotIn('dedup_hits_statesig', by_halogen[halogen])
            
            # But distributable metrics should be present
            self.assertIn('isotope_unavailable', by_halogen[halogen])
            self.assertIn('no_product_matches', by_halogen[halogen])
        
        for rule in ['R1', 'R3']:
            for halogen in ['F', 'Cl']:
                self.assertNotIn('attempts', by_rule_halogen[rule][halogen])
                self.assertNotIn('products', by_rule_halogen[rule][halogen])
                self.assertNotIn('dedup_hits_statesig', by_rule_halogen[rule][halogen])
                
                # But distributable metrics should be present
                self.assertIn('isotope_unavailable', by_rule_halogen[rule][halogen])
                self.assertIn('no_product_matches', by_rule_halogen[rule][halogen])
        
        # Verify distributable metrics sum consistency (overview counters are exempt)
        total_isotope_from_2d = sum(by_rule_halogen[rule][halogen]['isotope_unavailable'] 
                                   for rule in ['R1', 'R3'] for halogen in ['F', 'Cl'])
        self.assertEqual(total['isotope_unavailable'], total_isotope_from_2d)
        
        total_no_product_from_2d = sum(by_rule_halogen[rule][halogen]['no_product_matches'] 
                                      for rule in ['R1', 'R3'] for halogen in ['F', 'Cl'])
        self.assertEqual(total['no_product_matches'], total_no_product_from_2d)

    def test_v1_constructed_total_shape_is_flat(self):
        """Test that v1 constructed from legacy input has flat numeric total."""
        from src.halogenator.report import write_qa_summary_json
        
        # Create legacy input with complex structure (including nested qa_paths)
        legacy_input_with_complex_structure = {
            'isotope_unavailable': 8,
            'atommap_used': 4,
            'no_product_matches': 2,
            'template_unsupported': 1,
            'qa_paths': {  # Nested structure that should NOT appear in total
                'isotope_unavailable': 8,
                'atommap_used': 4,
                'heuristic_used': 1
            },
            'some_other_key': 'non_numeric_data',  # Non-numeric data
            'metadata': {
                'halogens': ['F', 'Cl'],
                'rules': ['R1', 'R3'],
                'description': 'should not appear in total'
            }
        }
        
        qa_json_path = write_qa_summary_json(legacy_input_with_complex_structure, self.test_dir)
        
        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)
        
        # Should construct v1 granular
        self.assertEqual(qa_data['version'], '1')
        self.assertIn('by_rule', qa_data)
        self.assertIn('total', qa_data)
        
        total = qa_data['total']
        by_rule_halogen = qa_data['by_rule_halogen']
        
        # Check that total is flat (no nested structures)
        disallowed_keys = ['metadata', 'qa_paths', 'pivots', 'by_rule', 'by_halogen', 'by_rule_halogen']
        for key in disallowed_keys:
            self.assertNotIn(key, total, f"total should not contain nested structure key: {key}")
        
        # Check that total only contains numeric values
        for key, value in total.items():
            self.assertIsInstance(value, (int, float), f"total[{key}] should be numeric, got {type(value)}")
        
        # Check that total metrics match 2D aggregation for distributable metrics
        distributable_metrics = ['isotope_unavailable', 'atommap_used', 'no_product_matches', 'template_unsupported']
        for metric in distributable_metrics:
            if metric in total:
                total_from_2d = sum(by_rule_halogen[rule][halogen][metric] 
                                   for rule in ['R1', 'R3'] for halogen in ['F', 'Cl'])
                self.assertEqual(total[metric], total_from_2d,
                               f"total[{metric}]: {total[metric]} != 2D aggregation {total_from_2d}")
        
        # Verify expected distributable metrics are present in total
        for metric in distributable_metrics:
            self.assertIn(metric, total, f"Distributable metric {metric} should be in total")

    def test_v1_remainder_bias_documented(self):
        """Test that v1 remainder allocation bias is properly documented and follows lexicographic order."""
        from src.halogenator.report import write_qa_summary_json
        
        # Create input with 1 remainder to test bias clearly
        # 1 item distributed among ['Br', 'Cl', 'F', 'I'] (lexicographic) should go to 'Br'
        remainder_test_input = {
            'isotope_unavailable': 1,  # 1 / 4 halogens = 0 remainder 1 -> goes to first in lexicographic order
            'metadata': {
                'halogens': ['F', 'I', 'Cl', 'Br'],  # Input order (intentionally not lexicographic)
                'rules': ['R3', 'R1']  # Input order (intentionally not lexicographic)
            }
        }
        
        qa_json_path = write_qa_summary_json(remainder_test_input, self.test_dir)
        
        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)
        
        # Check metadata documents distribution order
        metadata = qa_data['metadata']
        self.assertIn('distribution', metadata)
        self.assertIn('distribution_order', metadata)
        self.assertEqual(metadata['distribution'], "equal_with_lexicographic_remainder")
        self.assertEqual(metadata['distribution_order']['rules'], "lexicographic")
        self.assertEqual(metadata['distribution_order']['halogens'], "lexicographic")
        
        # Check that remainder goes to lexicographically first halogen
        by_halogen = qa_data['by_halogen']
        
        # In lexicographic order: ['Br', 'Cl', 'F', 'I']
        # The 1 remainder should go to 'Br' (first in lexicographic order)
        halogen_isotope_values = {h: by_halogen[h]['isotope_unavailable'] for h in by_halogen.keys()}
        
        # Br should get the remainder (1), others should get 0
        self.assertEqual(halogen_isotope_values['Br'], 1, "Br should get the remainder as first in lexicographic order")
        self.assertEqual(halogen_isotope_values['Cl'], 0, "Cl should get 0")
        self.assertEqual(halogen_isotope_values['F'], 0, "F should get 0") 
        self.assertEqual(halogen_isotope_values['I'], 0, "I should get 0")
        
        # Check same bias for rules: ['R1', 'R3'] lexicographically
        by_rule = qa_data['by_rule']
        rule_isotope_values = {r: by_rule[r]['isotope_unavailable'] for r in by_rule.keys()}
        
        # With 2 rules and 1 total: 1 / 2 = 0 remainder 1
        # R1 comes first lexicographically, so should get the remainder
        self.assertEqual(rule_isotope_values['R1'], 1, "R1 should get the remainder as first in lexicographic order")
        self.assertEqual(rule_isotope_values['R3'], 0, "R3 should get 0")
        
        # Verify that 2D structure is consistent
        by_rule_halogen = qa_data['by_rule_halogen']
        total_from_2d = sum(by_rule_halogen[r][h]['isotope_unavailable']
                           for r in ['R1', 'R3'] for h in ['Br', 'Cl', 'F', 'I'])
        self.assertEqual(total_from_2d, 1, "2D aggregation should equal original total")

    def test_v1_passthrough_total_includes_overview_counters(self):
        """Test that v1 passthrough merges top-level overview counters into total when total is missing."""
        from src.halogenator.report import write_qa_summary_json

        # Create v1 passthrough input with complete granular structure but no total, plus top-level overview counters
        passthrough_input_without_total = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': 3, 'no_product_matches': 1, 'atommap_used': 0, 'carbonyl_unknown': 0, 'heuristic_used': 0, 'isotope_miss': 0, 'rdkit_error': 0, 'template_unsupported': 0},
                'R3': {'isotope_unavailable': 2, 'no_product_matches': 2, 'atommap_used': 0, 'carbonyl_unknown': 0, 'heuristic_used': 0, 'isotope_miss': 0, 'rdkit_error': 0, 'template_unsupported': 0}
            },
            'by_halogen': {
                'F': {'isotope_unavailable': 2, 'no_product_matches': 1, 'atommap_used': 0, 'carbonyl_unknown': 0, 'heuristic_used': 0, 'isotope_miss': 0, 'rdkit_error': 0, 'template_unsupported': 0},
                'Cl': {'isotope_unavailable': 3, 'no_product_matches': 2, 'atommap_used': 0, 'carbonyl_unknown': 0, 'heuristic_used': 0, 'isotope_miss': 0, 'rdkit_error': 0, 'template_unsupported': 0}
            },
            'by_rule_halogen': {
                'R1': {
                    'F': {'isotope_unavailable': 1, 'no_product_matches': 0, 'atommap_used': 0, 'carbonyl_unknown': 0, 'heuristic_used': 0, 'isotope_miss': 0, 'rdkit_error': 0, 'template_unsupported': 0},
                    'Cl': {'isotope_unavailable': 2, 'no_product_matches': 1, 'atommap_used': 0, 'carbonyl_unknown': 0, 'heuristic_used': 0, 'isotope_miss': 0, 'rdkit_error': 0, 'template_unsupported': 0}
                },
                'R3': {
                    'F': {'isotope_unavailable': 1, 'no_product_matches': 1, 'atommap_used': 0, 'carbonyl_unknown': 0, 'heuristic_used': 0, 'isotope_miss': 0, 'rdkit_error': 0, 'template_unsupported': 0},
                    'Cl': {'isotope_unavailable': 1, 'no_product_matches': 1, 'atommap_used': 0, 'carbonyl_unknown': 0, 'heuristic_used': 0, 'isotope_miss': 0, 'rdkit_error': 0, 'template_unsupported': 0}
                }
            },
            'metadata': {
                'halogens': ['F', 'Cl'],
                'rules': ['R1', 'R3']
            },
            # Overview counters at top level (should be merged into total)
            'attempts': 100,
            'products': 50,
            'dedup_hits_statesig': 5,
            'dedup_hits_inchi': 3
        }

        qa_json_path = write_qa_summary_json(passthrough_input_without_total, self.test_dir)

        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)

        # Should be v1 passthrough with completed structure
        self.assertEqual(qa_data['version'], '1')
        self.assertIn('by_rule', qa_data)
        self.assertIn('total', qa_data)

        total = qa_data['total']
        by_rule = qa_data['by_rule']
        by_halogen = qa_data['by_halogen']

        # Total should include both distributable metrics (from 2D aggregation) AND overview counters
        self.assertIn('isotope_unavailable', total)
        self.assertIn('no_product_matches', total)
        self.assertIn('attempts', total)
        self.assertIn('products', total)
        self.assertIn('dedup_hits_statesig', total)
        self.assertIn('dedup_hits_inchi', total)

        # Overview counters should match input values
        self.assertEqual(total['attempts'], 100)
        self.assertEqual(total['products'], 50)
        self.assertEqual(total['dedup_hits_statesig'], 5)
        self.assertEqual(total['dedup_hits_inchi'], 3)

        # Distributable metrics in total should match 2D aggregation
        total_isotope_from_2d = sum(qa_data['by_rule_halogen'][r][h]['isotope_unavailable']
                                   for r in ['R1', 'R3'] for h in ['F', 'Cl'])
        self.assertEqual(total['isotope_unavailable'], total_isotope_from_2d)
        self.assertEqual(total['isotope_unavailable'], 5)  # 3 + 2 from by_rule

        # By_* slices should NOT contain overview counters (only distributable metrics)
        for rule in ['R1', 'R3']:
            self.assertNotIn('attempts', by_rule[rule])
            self.assertNotIn('products', by_rule[rule])
            self.assertNotIn('dedup_hits_statesig', by_rule[rule])

        for halogen in ['F', 'Cl']:
            self.assertNotIn('attempts', by_halogen[halogen])
            self.assertNotIn('products', by_halogen[halogen])
            self.assertNotIn('dedup_hits_statesig', by_halogen[halogen])

    def test_v1_passthrough_preserve_existing_total(self):
        """Test that v1 passthrough preserves existing total without modification."""
        from src.halogenator.report import write_qa_summary_json

        # Create v1 passthrough input WITH existing total, plus top-level overview counters
        passthrough_input_with_total = {
            'version': '1',
            'total': {
                'isotope_unavailable': 10,  # Intentionally different from by_rule sum
                'no_product_matches': 8,
                'attempts': 200,            # Pre-existing in total
                'custom_metric': 42         # Custom field that should be preserved
            },
            'by_rule': {
                'R1': {'isotope_unavailable': 3, 'no_product_matches': 1},
                'R3': {'isotope_unavailable': 2, 'no_product_matches': 2}
            },
            'metadata': {
                'halogens': ['F', 'Cl'],
                'rules': ['R1', 'R3']
            },
            # Additional overview counters at top level (should be ignored when total exists)
            'attempts': 999,  # Different from total.attempts - should be ignored
            'products': 777   # Not in total - should be ignored
        }

        qa_json_path = write_qa_summary_json(passthrough_input_with_total, self.test_dir)

        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)

        # Should be v1 passthrough with preserved total
        self.assertEqual(qa_data['version'], '1')
        self.assertIn('by_rule', qa_data)
        self.assertIn('total', qa_data)

        total = qa_data['total']

        # Total should be exactly as provided in input (no modifications)
        self.assertEqual(total['isotope_unavailable'], 10)  # Original value, not 2D sum (5)
        self.assertEqual(total['no_product_matches'], 8)
        self.assertEqual(total['attempts'], 200)  # Original total.attempts, not top-level 999
        self.assertEqual(total['custom_metric'], 42)

        # Top-level products should NOT be added to total (since total already exists)
        self.assertNotIn('products', total)

        # Structure should be completed as usual
        self.assertIn('by_halogen', qa_data)
        self.assertIn('by_rule_halogen', qa_data)

    def test_v1_distribution_keys_emitted(self):
        """Test that v1 granular outputs include explicit distribution keys for transparency."""
        from src.halogenator.report import write_qa_summary_json

        # Test constructed v1 granular (from legacy input)
        legacy_input = {
            'isotope_unavailable': 7,
            'no_product_matches': 3,
            'metadata': {
                'halogens': ['I', 'F', 'Cl', 'Br'],  # Intentionally not lexicographic
                'rules': ['R3', 'R1', 'R5']           # Intentionally not lexicographic
            }
        }

        qa_json_path = write_qa_summary_json(legacy_input, self.test_dir)

        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)

        # Should have v1 granular structure
        self.assertEqual(qa_data['version'], '1')
        self.assertIn('by_rule', qa_data)

        metadata = qa_data['metadata']

        # Should have distribution_keys field
        self.assertIn('distribution_keys', metadata)

        distribution_keys = metadata['distribution_keys']
        self.assertIn('rules', distribution_keys)
        self.assertIn('halogens', distribution_keys)

        # Distribution keys should be in lexicographic order (regardless of input order)
        expected_rules = ['R1', 'R3', 'R5']  # Lexicographic order
        expected_halogens = ['Br', 'Cl', 'F', 'I']  # Lexicographic order

        self.assertEqual(distribution_keys['rules'], expected_rules)
        self.assertEqual(distribution_keys['halogens'], expected_halogens)

        # Test passthrough v1 granular (from existing granular input)
        passthrough_input = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': 3, 'no_product_matches': 1, 'atommap_used': 0, 'carbonyl_unknown': 0, 'heuristic_used': 0, 'isotope_miss': 0, 'rdkit_error': 0, 'template_unsupported': 0},
                'R5': {'isotope_unavailable': 2, 'no_product_matches': 2, 'atommap_used': 0, 'carbonyl_unknown': 0, 'heuristic_used': 0, 'isotope_miss': 0, 'rdkit_error': 0, 'template_unsupported': 0}
            },
            'metadata': {
                'halogens': ['Cl', 'F'],  # Different input order
                'rules': ['R5', 'R1']     # Different input order
            }
        }

        qa_json_path_2 = write_qa_summary_json(passthrough_input, self.test_dir)

        with open(qa_json_path_2, 'r') as f:
            qa_data_2 = json.load(f)

        # Should also have distribution_keys in lexicographic order
        metadata_2 = qa_data_2['metadata']
        self.assertIn('distribution_keys', metadata_2)

        distribution_keys_2 = metadata_2['distribution_keys']
        expected_rules_2 = ['R1', 'R5']  # Lexicographic order
        expected_halogens_2 = ['Cl', 'F']  # Lexicographic order

        self.assertEqual(distribution_keys_2['rules'], expected_rules_2)
        self.assertEqual(distribution_keys_2['halogens'], expected_halogens_2)

    def test_v1_passthrough_filters_nondistributable_from_granular(self):
        """Test that v1 passthrough filters out non-distributable fields from by_* dimensions."""
        from src.halogenator.report import write_qa_summary_json

        # Create v1 passthrough input where by_rule mistakenly includes overview counters
        # This simulates corrupted input where someone put overview counters in granular dimensions
        passthrough_input_with_overview_in_granular = {
            'version': '1',
            'by_rule': {
                'R1': {
                    'isotope_unavailable': 3,     # Distributable (should remain)
                    'no_product_matches': 1,      # Distributable (should remain)
                    'attempts': 50,               # Overview counter (should be filtered out)
                    'products': 25                # Overview counter (should be filtered out)
                },
                'R3': {
                    'isotope_unavailable': 2,     # Distributable (should remain)
                    'no_product_matches': 2,      # Distributable (should remain)
                    'attempts': 40,               # Overview counter (should be filtered out)
                    'products': 20                # Overview counter (should be filtered out)
                }
            },
            'metadata': {
                'halogens': ['F', 'Cl'],
                'rules': ['R1', 'R3']
            },
            # Additional overview counters at top level (should be merged into total if missing)
            'attempts': 100,  # Should override the bad granular values
            'products': 50
        }

        qa_json_path = write_qa_summary_json(passthrough_input_with_overview_in_granular, self.test_dir)

        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)

        # Should be v1 passthrough with completed structure
        self.assertEqual(qa_data['version'], '1')
        self.assertIn('by_rule', qa_data)
        self.assertIn('by_halogen', qa_data)
        self.assertIn('by_rule_halogen', qa_data)
        self.assertIn('total', qa_data)

        by_rule = qa_data['by_rule']
        by_halogen = qa_data['by_halogen']
        by_rule_halogen = qa_data['by_rule_halogen']
        total = qa_data['total']

        # Check that overview counters are filtered out from by_* dimensions
        for rule in ['R1', 'R3']:
            # Should have distributable metrics
            self.assertIn('isotope_unavailable', by_rule[rule])
            self.assertIn('no_product_matches', by_rule[rule])

            # Should NOT have overview counters (filtered out during completion)
            self.assertNotIn('attempts', by_rule[rule])
            self.assertNotIn('products', by_rule[rule])

        for halogen in ['F', 'Cl']:
            # Should have distributable metrics (even if 0)
            self.assertIn('isotope_unavailable', by_halogen[halogen])
            self.assertIn('no_product_matches', by_halogen[halogen])

            # Should NOT have overview counters
            self.assertNotIn('attempts', by_halogen[halogen])
            self.assertNotIn('products', by_halogen[halogen])

        for rule in ['R1', 'R3']:
            for halogen in ['F', 'Cl']:
                # Should have distributable metrics (even if 0)
                self.assertIn('isotope_unavailable', by_rule_halogen[rule][halogen])
                self.assertIn('no_product_matches', by_rule_halogen[rule][halogen])

                # Should NOT have overview counters
                self.assertNotIn('attempts', by_rule_halogen[rule][halogen])
                self.assertNotIn('products', by_rule_halogen[rule][halogen])

        # Total should include overview counters from top-level (since no total was provided)
        self.assertIn('attempts', total)
        self.assertIn('products', total)
        self.assertEqual(total['attempts'], 100)  # From top-level, not from corrupted granular
        self.assertEqual(total['products'], 50)

        # Total should also include distributable metrics (from 2D aggregation)
        self.assertIn('isotope_unavailable', total)
        self.assertIn('no_product_matches', total)

    def test_v1_passthrough_cleans_overview_from_all_granular_dims(self):
        """Test that v1 passthrough removes overview counters from all granular dimensions."""
        from src.halogenator.report import write_qa_summary_json

        # Create v1 input where overview counters are mistakenly placed in ALL granular dimensions
        passthrough_input_with_overview_in_all_dims = {
            'version': '1',
            'by_rule': {
                'R1': {
                    'isotope_unavailable': 3,
                    'no_product_matches': 1,
                    'attempts': 50,               # Should be removed
                    'products': 25                # Should be removed
                }
            },
            'by_halogen': {
                'F': {
                    'isotope_unavailable': 2,
                    'no_product_matches': 1,
                    'attempts': 30,               # Should be removed
                    'dedup_hits_statesig': 5      # Should be removed
                }
            },
            'by_rule_halogen': {
                'R1': {
                    'F': {
                        'isotope_unavailable': 1,
                        'no_product_matches': 0,
                        'attempts': 10,           # Should be removed
                        'products': 5             # Should be removed
                    }
                }
            },
            'metadata': {
                'halogens': ['F'],
                'rules': ['R1']
            },
            # Top-level overview counters (should be preserved in total)
            'attempts': 100,
            'products': 50
        }

        qa_json_path = write_qa_summary_json(passthrough_input_with_overview_in_all_dims, self.test_dir)

        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)

        # Should be v1 passthrough with sanitized structure
        self.assertEqual(qa_data['version'], '1')
        self.assertIn('by_rule', qa_data)
        self.assertIn('by_halogen', qa_data)
        self.assertIn('by_rule_halogen', qa_data)
        self.assertIn('total', qa_data)

        # Check that ALL granular dimensions have been cleaned
        by_rule = qa_data['by_rule']
        by_halogen = qa_data['by_halogen']
        by_rule_halogen = qa_data['by_rule_halogen']
        total = qa_data['total']

        # by_rule should NOT have overview counters
        self.assertIn('isotope_unavailable', by_rule['R1'])
        self.assertIn('no_product_matches', by_rule['R1'])
        self.assertNotIn('attempts', by_rule['R1'])
        self.assertNotIn('products', by_rule['R1'])

        # by_halogen should NOT have overview counters
        self.assertIn('isotope_unavailable', by_halogen['F'])
        self.assertIn('no_product_matches', by_halogen['F'])
        self.assertNotIn('attempts', by_halogen['F'])
        self.assertNotIn('dedup_hits_statesig', by_halogen['F'])

        # by_rule_halogen should NOT have overview counters
        self.assertIn('isotope_unavailable', by_rule_halogen['R1']['F'])
        self.assertIn('no_product_matches', by_rule_halogen['R1']['F'])
        self.assertNotIn('attempts', by_rule_halogen['R1']['F'])
        self.assertNotIn('products', by_rule_halogen['R1']['F'])

        # Total should include overview counters from top-level (not from cleaned granular)
        self.assertIn('attempts', total)
        self.assertIn('products', total)
        self.assertEqual(total['attempts'], 100)  # From top-level, not granular
        self.assertEqual(total['products'], 50)

    def test_v1_passthrough_normalizes_missing_metrics_in_existing_2d(self):
        """Test that v1 passthrough fills missing metrics with 0 in existing 2D structure."""
        from src.halogenator.report import write_qa_summary_json

        # Create v1 input with incomplete 2D structure (missing some metrics in some cells)
        passthrough_input_with_incomplete_2d = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': 3, 'no_product_matches': 1},  # Complete
                'R3': {'isotope_unavailable': 2}  # Missing no_product_matches
            },
            'by_halogen': {
                'F': {'isotope_unavailable': 2, 'no_product_matches': 1},  # Complete
                'Cl': {'no_product_matches': 2}  # Missing isotope_unavailable
            },
            'by_rule_halogen': {
                'R1': {
                    'F': {'isotope_unavailable': 1},  # Missing no_product_matches
                    'Cl': {'isotope_unavailable': 2, 'no_product_matches': 1}  # Complete
                },
                'R3': {
                    'F': {'no_product_matches': 1},  # Missing isotope_unavailable
                    'Cl': {}  # Missing both metrics
                }
            },
            'metadata': {
                'halogens': ['F', 'Cl'],
                'rules': ['R1', 'R3']
            }
        }

        qa_json_path = write_qa_summary_json(passthrough_input_with_incomplete_2d, self.test_dir)

        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)

        # Should be v1 passthrough with normalized structure
        self.assertEqual(qa_data['version'], '1')
        self.assertIn('by_rule', qa_data)
        self.assertIn('by_halogen', qa_data)
        self.assertIn('by_rule_halogen', qa_data)
        self.assertIn('total', qa_data)

        by_rule = qa_data['by_rule']
        by_halogen = qa_data['by_halogen']
        by_rule_halogen = qa_data['by_rule_halogen']

        # NEW 2D-SoT BEHAVIOR: marginals are derived from normalized 2D structure
        # After 2D normalization, marginals are computed by summing across 2D dimensions

        # Expected 2D structure after normalization (missing metrics filled with 0):
        # R1: F(isotope_unavailable=1, no_product_matches=0), Cl(isotope_unavailable=2, no_product_matches=1)
        # R3: F(isotope_unavailable=0, no_product_matches=1), Cl(isotope_unavailable=0, no_product_matches=0)

        # by_rule derived from 2D (summing across halogens):
        # R1: isotope_unavailable=1+2=3, no_product_matches=0+1=1
        # R3: isotope_unavailable=0+0=0, no_product_matches=1+0=1
        self.assertEqual(by_rule['R1']['isotope_unavailable'], 3)  # Derived from 2D
        self.assertEqual(by_rule['R1']['no_product_matches'], 1)   # Derived from 2D
        self.assertEqual(by_rule['R3']['isotope_unavailable'], 0)  # Derived from 2D
        self.assertEqual(by_rule['R3']['no_product_matches'], 1)   # Derived from 2D

        # by_halogen derived from 2D (summing across rules):
        # F: isotope_unavailable=1+0=1, no_product_matches=0+1=1
        # Cl: isotope_unavailable=2+0=2, no_product_matches=1+0=1
        self.assertEqual(by_halogen['F']['isotope_unavailable'], 1)  # Derived from 2D
        self.assertEqual(by_halogen['F']['no_product_matches'], 1)   # Derived from 2D
        self.assertEqual(by_halogen['Cl']['isotope_unavailable'], 2) # Derived from 2D
        self.assertEqual(by_halogen['Cl']['no_product_matches'], 1)  # Derived from 2D

        # by_rule_halogen normalization
        self.assertEqual(by_rule_halogen['R1']['F']['isotope_unavailable'], 1)  # Original value
        self.assertEqual(by_rule_halogen['R1']['F']['no_product_matches'], 0)   # Filled with 0
        self.assertEqual(by_rule_halogen['R1']['Cl']['isotope_unavailable'], 2) # Original value
        self.assertEqual(by_rule_halogen['R1']['Cl']['no_product_matches'], 1)  # Original value
        self.assertEqual(by_rule_halogen['R3']['F']['isotope_unavailable'], 0)  # Filled with 0
        self.assertEqual(by_rule_halogen['R3']['F']['no_product_matches'], 1)   # Original value
        self.assertEqual(by_rule_halogen['R3']['Cl']['isotope_unavailable'], 0) # Filled with 0
        self.assertEqual(by_rule_halogen['R3']['Cl']['no_product_matches'], 0)  # Filled with 0

        # Verify that total calculation doesn't crash and produces reasonable results
        total = qa_data['total']
        self.assertIn('isotope_unavailable', total)
        self.assertIn('no_product_matches', total)
        # Total should be sum of 2D: isotope_unavailable = 1+2+0+0 = 3, no_product_matches = 0+1+1+0 = 2
        self.assertEqual(total['isotope_unavailable'], 3)
        self.assertEqual(total['no_product_matches'], 2)

    def test_v1_passthrough_total_from_by_rule_when_no_2d(self):
        """Test that v1 passthrough can compute total from by_rule when 2D structure is missing."""
        from src.halogenator.report import write_qa_summary_json

        # Create v1 passthrough input with only by_rule (no by_halogen, no by_rule_halogen)
        passthrough_input_only_by_rule = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': 5, 'no_product_matches': 2},
                'R3': {'isotope_unavailable': 3, 'no_product_matches': 1}
            },
            'metadata': {
                'halogens': ['F', 'Cl'],
                'rules': ['R1', 'R3']
            },
            # Top-level overview counters
            'attempts': 100,
            'products': 50
        }

        qa_json_path = write_qa_summary_json(passthrough_input_only_by_rule, self.test_dir)

        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)

        # Should be v1 passthrough with completed structure
        self.assertEqual(qa_data['version'], '1')
        self.assertIn('by_rule', qa_data)
        self.assertIn('by_halogen', qa_data)  # Should be completed
        self.assertIn('by_rule_halogen', qa_data)  # Should be completed
        self.assertIn('total', qa_data)

        by_rule = qa_data['by_rule']
        total = qa_data['total']

        # by_rule should preserve original values
        self.assertEqual(by_rule['R1']['isotope_unavailable'], 5)
        self.assertEqual(by_rule['R1']['no_product_matches'], 2)
        self.assertEqual(by_rule['R3']['isotope_unavailable'], 3)
        self.assertEqual(by_rule['R3']['no_product_matches'], 1)

        # Total should aggregate from by_rule (not be zeros from empty 2D)
        self.assertEqual(total['isotope_unavailable'], 8)  # 5 + 3 from by_rule
        self.assertEqual(total['no_product_matches'], 3)   # 2 + 1 from by_rule

        # Overview counters should also be included
        self.assertEqual(total['attempts'], 100)
        self.assertEqual(total['products'], 50)

        # by_halogen and by_rule_halogen should be completed but zero-filled (default mode)
        by_halogen = qa_data['by_halogen']
        by_rule_halogen = qa_data['by_rule_halogen']

        # Default zero-fill mode: missing dimensions should be filled with zeros
        for halogen in ['F', 'Cl']:
            self.assertEqual(by_halogen[halogen]['isotope_unavailable'], 0)
            self.assertEqual(by_halogen[halogen]['no_product_matches'], 0)
            for rule in ['R1', 'R3']:
                self.assertEqual(by_rule_halogen[rule][halogen]['isotope_unavailable'], 0)
                self.assertEqual(by_rule_halogen[rule][halogen]['no_product_matches'], 0)

    def test_v1_passthrough_zero_fill_mode_keeps_missing_dims_zero(self):
        """Test that default zero-fill mode doesn't distribute, only fills missing dimensions with zeros."""
        from src.halogenator.report import write_qa_summary_json

        # Create input with only by_rule structure
        zero_fill_input = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': 5, 'atommap_used': 3},
                'R2': {'isotope_unavailable': 2, 'atommap_used': 1}
            },
            'metadata': {
                'halogens': ['F', 'Cl'],
                'rules': ['R1', 'R2']
            }
        }

        # Test with default mode (zero_fill)
        qa_json_path = write_qa_summary_json(zero_fill_input, self.test_dir)

        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)

        # Verify total is computed correctly from by_rule
        total = qa_data['total']
        self.assertEqual(total['isotope_unavailable'], 7)  # 5 + 2
        self.assertEqual(total['atommap_used'], 4)  # 3 + 1

        # Verify missing dimensions are zero-filled
        by_halogen = qa_data['by_halogen']
        by_rule_halogen = qa_data['by_rule_halogen']

        for halogen in ['F', 'Cl']:
            self.assertEqual(by_halogen[halogen]['isotope_unavailable'], 0)
            self.assertEqual(by_halogen[halogen]['atommap_used'], 0)

        for rule in ['R1', 'R2']:
            for halogen in ['F', 'Cl']:
                self.assertEqual(by_rule_halogen[rule][halogen]['isotope_unavailable'], 0)
                self.assertEqual(by_rule_halogen[rule][halogen]['atommap_used'], 0)

        # Verify metadata shows zero_fill strategy
        metadata = qa_data['metadata']
        self.assertEqual(metadata['completion']['strategy'], 'zero_fill')
        self.assertFalse(metadata['completion']['enabled'])

    def test_v1_passthrough_distribute_mode_derives_2d_and_marginals(self):
        """Test that distribute mode derives missing dimensions by distributing marginal values."""
        from src.halogenator.report import write_qa_summary_json

        # Create input with only by_rule structure
        distribute_input = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': 5, 'atommap_used': 3},
                'R2': {'isotope_unavailable': 2, 'atommap_used': 1}
            },
            'metadata': {
                'halogens': ['F', 'Cl'],
                'rules': ['R1', 'R2']
            }
        }

        # Test with distribute mode
        qa_json_path = write_qa_summary_json(distribute_input, self.test_dir, completion_mode='distribute')

        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)

        # Verify total is computed correctly from by_rule
        total = qa_data['total']
        self.assertEqual(total['isotope_unavailable'], 7)  # 5 + 2
        self.assertEqual(total['atommap_used'], 4)  # 3 + 1

        # Verify total equals 2D aggregation
        by_rule_halogen = qa_data['by_rule_halogen']
        total_from_2d = sum(by_rule_halogen[rule][halogen]['isotope_unavailable']
                           for rule in ['R1', 'R2'] for halogen in ['F', 'Cl'])
        self.assertEqual(total['isotope_unavailable'], total_from_2d)

        # Verify distributed values
        # R1: isotope_unavailable=5 -> Cl:3, F:2 (remainder to first lexicographically)
        # R2: isotope_unavailable=2 -> Cl:1, F:1
        self.assertEqual(by_rule_halogen['R1']['Cl']['isotope_unavailable'], 3)
        self.assertEqual(by_rule_halogen['R1']['F']['isotope_unavailable'], 2)
        self.assertEqual(by_rule_halogen['R2']['Cl']['isotope_unavailable'], 1)
        self.assertEqual(by_rule_halogen['R2']['F']['isotope_unavailable'], 1)

        # Verify by_halogen is derived by summing across rules
        by_halogen = qa_data['by_halogen']
        self.assertEqual(by_halogen['Cl']['isotope_unavailable'], 4)  # 3 + 1
        self.assertEqual(by_halogen['F']['isotope_unavailable'], 3)   # 2 + 1

        # Verify metadata shows distribute strategy
        metadata = qa_data['metadata']
        self.assertEqual(metadata['completion']['strategy'], 'distribute')
        self.assertTrue(metadata['completion']['enabled'])
        self.assertEqual(metadata['completion']['base'], 'by_rule')

    def test_v1_completion_uses_distribution_order_keys(self):
        """Test that completion uses distribution order keys consistent with construction."""
        from src.halogenator.report import write_qa_summary_json

        # Create input with non-lexicographic order for rules and halogens
        non_lexicographic_input = {
            'version': '1',
            'by_rule': {
                'R3': {'isotope_unavailable': 7},  # R3 before R1 (non-lexicographic)
                'R1': {'isotope_unavailable': 2}
            },
            'metadata': {
                'halogens': ['I', 'Br', 'Cl', 'F'],  # Reverse lexicographic order
                'rules': ['R3', 'R1']  # Non-lexicographic order
            }
        }

        # Test with distribute mode
        qa_json_path = write_qa_summary_json(non_lexicographic_input, self.test_dir, completion_mode='distribute')

        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)

        # Check that completion metadata uses lexicographic order
        metadata = qa_data['metadata']
        completion_keys = metadata['completion']['ordered_keys']
        self.assertEqual(completion_keys['rules'], ['R1', 'R3'])  # Lexicographic order
        self.assertEqual(completion_keys['halogens'], ['Br', 'Cl', 'F', 'I'])  # Lexicographic order

        # Check that remainder distribution follows lexicographic order
        # R3 has isotope_unavailable=7, distributed across 4 halogens: 7//4=1 base, 3 remainder
        # Remainder goes to first 3 halogens lexicographically: Br, Cl, F (each gets +1)
        by_rule_halogen = qa_data['by_rule_halogen']
        self.assertEqual(by_rule_halogen['R3']['Br']['isotope_unavailable'], 2)  # 1 + 1 remainder
        self.assertEqual(by_rule_halogen['R3']['Cl']['isotope_unavailable'], 2)  # 1 + 1 remainder
        self.assertEqual(by_rule_halogen['R3']['F']['isotope_unavailable'], 2)   # 1 + 1 remainder
        self.assertEqual(by_rule_halogen['R3']['I']['isotope_unavailable'], 1)   # 1 + 0 remainder

        # Verify that by_halogen marginals sum correctly
        by_halogen = qa_data['by_halogen']
        total_from_marginals = sum(by_halogen[h]['isotope_unavailable'] for h in ['Br', 'Cl', 'F', 'I'])
        self.assertEqual(total_from_marginals, 9)  # 7 + 2 from both rules

    def test_v1_passthrough_idempotent(self):
        """Test that passthrough pipeline is idempotent - applying it twice yields same result."""
        from src.halogenator.report import write_qa_summary_json
        import json

        # Create initial input
        initial_input = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': 5, 'atommap_used': 2},
                'R2': {'isotope_unavailable': 3, 'atommap_used': 1}
            },
            'metadata': {
                'halogens': ['F', 'Cl'],
                'rules': ['R1', 'R2']
            },
            'attempts': 100,
            'products': 50
        }

        # First pass
        qa_json_path_1 = write_qa_summary_json(initial_input, self.test_dir, completion_mode='distribute')
        with open(qa_json_path_1, 'r') as f:
            qa_data_1 = json.load(f)

        # Second pass - use output of first pass as input, but preserve same configuration
        second_pass_input = qa_data_1.copy()
        # Ensure same halogens/rules configuration for consistent processing
        second_pass_input['metadata']['halogens'] = ['F', 'Cl']  # Same as initial input
        second_pass_input['metadata']['rules'] = ['R1', 'R2']    # Same as initial input

        qa_json_path_2 = write_qa_summary_json(second_pass_input, self.test_dir, completion_mode='distribute')
        with open(qa_json_path_2, 'r') as f:
            qa_data_2 = json.load(f)

        # Results should be identical (idempotent)
        self.assertEqual(qa_data_1['total'], qa_data_2['total'])
        self.assertEqual(qa_data_1['by_rule'], qa_data_2['by_rule'])
        self.assertEqual(qa_data_1['by_halogen'], qa_data_2['by_halogen'])
        self.assertEqual(qa_data_1['by_rule_halogen'], qa_data_2['by_rule_halogen'])

    def test_v1_passthrough_total_equals_sum_2d_for_distributables(self):
        """Test that distributable metrics in total equal sum of 2D structure."""
        from src.halogenator.report import write_qa_summary_json

        # Create input with mixed distributable and overview metrics
        mixed_input = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': 4, 'no_product_matches': 2},
                'R3': {'isotope_unavailable': 5, 'no_product_matches': 1}
            },
            'metadata': {
                'halogens': ['F', 'Cl', 'Br'],
                'rules': ['R1', 'R3']
            },
            # Overview counters (not distributable)
            'attempts': 120,
            'products': 60,
            'dedup_hits_statesig': 5
        }

        # Test with distribute mode to ensure 2D consistency
        qa_json_path = write_qa_summary_json(mixed_input, self.test_dir, completion_mode='distribute')
        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)

        total = qa_data['total']
        by_rule_halogen = qa_data['by_rule_halogen']

        # For all distributable metrics, total should equal 2D sum
        distributable_metrics = ['isotope_unavailable', 'no_product_matches']
        for metric in distributable_metrics:
            total_from_2d = sum(
                by_rule_halogen[rule][halogen][metric]
                for rule in ['R1', 'R3']
                for halogen in ['F', 'Cl', 'Br']
            )
            self.assertEqual(total[metric], total_from_2d,
                           f"Distributable metric {metric}: total={total[metric]} != 2D_sum={total_from_2d}")

        # Overview counters should not be distributed to 2D
        self.assertEqual(total['attempts'], 120)
        self.assertEqual(total['products'], 60)
        self.assertEqual(total['dedup_hits_statesig'], 5)

        # Overview counters should not appear in 2D structure
        for rule in ['R1', 'R3']:
            for halogen in ['F', 'Cl', 'Br']:
                self.assertNotIn('attempts', by_rule_halogen[rule][halogen])
                self.assertNotIn('products', by_rule_halogen[rule][halogen])
                self.assertNotIn('dedup_hits_statesig', by_rule_halogen[rule][halogen])

    def test_v1_passthrough_total_from_by_halogen_when_no_2d(self):
        """Test that v1 passthrough can compute total from by_halogen when by_rule and 2D are missing."""
        from src.halogenator.report import write_qa_summary_json

        # Create v1 passthrough input with only by_halogen (no by_rule, no by_rule_halogen)
        passthrough_input_only_by_halogen = {
            'version': '1',
            'by_halogen': {
                'F': {'isotope_unavailable': 4, 'no_product_matches': 1},
                'Cl': {'isotope_unavailable': 6, 'no_product_matches': 2}
            },
            'metadata': {
                'halogens': ['F', 'Cl'],
                'rules': ['R1', 'R3']
            },
            # Top-level overview counters
            'attempts': 200,
            'dedup_hits_statesig': 10
        }

        qa_json_path = write_qa_summary_json(passthrough_input_only_by_halogen, self.test_dir)

        with open(qa_json_path, 'r') as f:
            qa_data = json.load(f)

        # Should be v1 passthrough with completed structure
        self.assertEqual(qa_data['version'], '1')
        self.assertIn('by_rule', qa_data)  # Should be completed
        self.assertIn('by_halogen', qa_data)
        self.assertIn('by_rule_halogen', qa_data)  # Should be completed
        self.assertIn('total', qa_data)

        by_halogen = qa_data['by_halogen']
        total = qa_data['total']

        # by_halogen should preserve original values
        self.assertEqual(by_halogen['F']['isotope_unavailable'], 4)
        self.assertEqual(by_halogen['F']['no_product_matches'], 1)
        self.assertEqual(by_halogen['Cl']['isotope_unavailable'], 6)
        self.assertEqual(by_halogen['Cl']['no_product_matches'], 2)

        # Total should aggregate from by_halogen (not be zeros from empty 2D/by_rule)
        self.assertEqual(total['isotope_unavailable'], 10)  # 4 + 6 from by_halogen
        self.assertEqual(total['no_product_matches'], 3)    # 1 + 2 from by_halogen

        # Overview counters should also be included
        self.assertEqual(total['attempts'], 200)
        self.assertEqual(total['dedup_hits_statesig'], 10)

        # by_rule and by_rule_halogen should be completed but zero-filled
        by_rule = qa_data['by_rule']
        by_rule_halogen = qa_data['by_rule_halogen']
        for rule in ['R1', 'R3']:
            self.assertEqual(by_rule[rule]['isotope_unavailable'], 0)
            self.assertEqual(by_rule[rule]['no_product_matches'], 0)
            for halogen in ['F', 'Cl']:
                self.assertEqual(by_rule_halogen[rule][halogen]['isotope_unavailable'], 0)
                self.assertEqual(by_rule_halogen[rule][halogen]['no_product_matches'], 0)


if __name__ == '__main__':
    unittest.main()