# -*- coding: ascii -*-
"""Test streaming report functionality to avoid OOM."""

import os
import tempfile
import unittest
from src.halogenator.io_utils import write_table, iter_table_records, read_smi, write_smi
from src.halogenator.report import generate_summary_report, _init_incremental_stats, _update_incremental_stats, _finalize_incremental_stats


class TestReportStreaming(unittest.TestCase):
    """Test streaming report processing for large datasets."""
    
    def test_iter_table_records_parquet(self):
        """Test that iter_table_records yields the same data as read_table."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create test data
            test_records = []
            for i in range(100):
                record = {
                    'parent_smiles': 'c1ccccc1',
                    'parent_name': f'compound_{i}',
                    'parent_inchikey': 'UHOVQNZJYSORNB-UHFFFAOYSA-N',
                    'product_smiles': 'Fc1ccccc1',
                    'smiles': 'Fc1ccccc1',
                    'inchikey': 'FNRMXHYUJQADQH-UHFFFAOYSA-N',
                    'k': 1,
                    'rule': 'R1',
                    'halogen': 'F',
                    'substitutions': [{'rule': 'R1', 'site': i % 6, 'halogen': 'F', 'depth': 1}],
                    'constraints_ok': True,
                    'constraints_violations': {}
                }
                test_records.append(record)
            
            # Write to parquet
            parquet_file = os.path.join(temp_dir, 'test.parquet')
            write_table(test_records, parquet_file)
            
            # Compare stream vs batch read
            from src.halogenator.io_utils import read_table
            batch_records = read_table(parquet_file)
            stream_records = list(iter_table_records(parquet_file))
            
            # Should have same count and content
            self.assertEqual(len(batch_records), len(stream_records))
            self.assertEqual(len(stream_records), 100)
            
            # Compare first record structure
            if stream_records:
                self.assertEqual(batch_records[0].keys(), stream_records[0].keys())
                self.assertEqual(batch_records[0]['substitutions'], stream_records[0]['substitutions'])
    
    def test_incremental_stats_equivalence(self):
        """Test that incremental stats produces same results as batch processing."""
        # Create test data
        test_records = []
        for halogen in ['F', 'Cl', 'Br']:
            for k in [1, 2]:
                record = {
                    'parent_smiles': 'c1ccccc1',
                    'parent_name': 'benzene',
                    'parent_inchikey': 'UHOVQNZJYSORNB-UHFFFAOYSA-N',
                    'product_smiles': f'{halogen}c1ccccc1',
                    'k': k,
                    'rule': 'R1',
                    'halogen': halogen,
                    'constraints_ok': True,
                    'constraints_violations': {}
                }
                test_records.append(record)
        
        # Process incrementally
        incremental_stats = _init_incremental_stats()
        for record in test_records:
            _update_incremental_stats(incremental_stats, record)
        
        parent_pairs = [('c1ccccc1', 'benzene')]
        config = {'schema_format': 'P1', 'subset': 'test'}
        _finalize_incremental_stats(incremental_stats, parent_pairs, config)
        
        # Verify key statistics
        self.assertEqual(incremental_stats['product_count'], 6)  # 3 halogens x 2 k values
        self.assertEqual(incremental_stats['parent_count'], 1)
        self.assertEqual(incremental_stats['k_counts'][1], 3)  # F, Cl, Br at k=1
        self.assertEqual(incremental_stats['k_counts'][2], 3)  # F, Cl, Br at k=2
        self.assertEqual(incremental_stats['halogen_counts']['F'], 2)  # k=1 and k=2
        self.assertEqual(incremental_stats['rule_counts']['R1'], 6)
        self.assertEqual(incremental_stats['constraint_stats']['total_passed'], 6)
        self.assertEqual(incremental_stats['constraint_stats']['total_failed'], 0)
    
    def test_streaming_report_with_batches(self):
        """Test full streaming report with multiple parquet files."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create parent file
            parents_file = os.path.join(temp_dir, 'parents.smi')
            write_smi([('c1ccccc1', 'benzene'), ('c1ccc(O)cc1', 'phenol')], parents_file)
            
            # Create multiple batch files
            batch1_records = []
            for i in range(50):
                record = {
                    'parent_smiles': 'c1ccccc1',
                    'parent_name': 'benzene',
                    'parent_inchikey': 'UHOVQNZJYSORNB-UHFFFAOYSA-N',
                    'product_smiles': 'Fc1ccccc1',
                    'smiles': 'Fc1ccccc1',
                    'inchikey': 'FNRMXHYUJQADQH-UHFFFAOYSA-N',
                    'k': 1,
                    'rule': 'R1',
                    'halogen': 'F',
                    'substitutions': [{'rule': 'R1', 'site': 0, 'halogen': 'F', 'depth': 1}],
                    'constraints_ok': True,
                    'constraints_violations': {}
                }
                batch1_records.append(record)
            
            batch2_records = []
            for i in range(30):
                record = {
                    'parent_smiles': 'c1ccc(O)cc1',
                    'parent_name': 'phenol',
                    'parent_inchikey': 'RWXXQIKFVTHVRS-UHFFFAOYSA-N',
                    'product_smiles': 'Fc1ccc(O)cc1',
                    'smiles': 'Fc1ccc(O)cc1',
                    'inchikey': 'VSCGXNMXTNTHQP-UHFFFAOYSA-N',
                    'k': 1,
                    'rule': 'R1',
                    'halogen': 'F',
                    'substitutions': [{'rule': 'R1', 'site': 0, 'halogen': 'F', 'depth': 1}],
                    'constraints_ok': True,
                    'constraints_violations': {}
                }
                batch2_records.append(record)
            
            # Write batch files
            main_file = os.path.join(temp_dir, 'products.parquet')
            part1_file = os.path.join(temp_dir, 'products.part1.parquet')
            
            write_table(batch1_records, main_file)
            write_table(batch2_records, part1_file)
            
            # Configure and run report
            config = {
                'io': {
                    'smiles_file': parents_file,
                    'products_table': main_file,
                    'summary_csv': os.path.join(temp_dir, 'summary.csv')
                }
            }
            
            try:
                generate_summary_report(config, subset='test', k_max=1)
                
                # Verify summary was created
                summary_file = config['io']['summary_csv']
                self.assertTrue(os.path.exists(summary_file), "Summary file should be created")
                
                # Check basic content structure 
                with open(summary_file, 'r') as f:
                    content = f.read()
                    # Should contain expected statistics
                    self.assertIn('total products generated,80', content.lower())
                    self.assertIn('total parent molecules,2', content.lower())
                    self.assertIn('unique parents with products,2', content.lower())
                    
            except Exception as e:
                self.fail(f"Streaming report generation failed: {e}")


if __name__ == '__main__':
    unittest.main()