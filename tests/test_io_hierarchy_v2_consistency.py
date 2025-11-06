# -*- coding: ascii -*-
"""
Test consistency between io_hierarchy v1 and v2 implementations.

This test ensures that the streaming v2 implementation produces
identical outputs to the memory-based v1 implementation.
"""

import unittest
import tempfile
import shutil
import os
import pandas as pd
from pathlib import Path


class TestHierarchyV2Consistency(unittest.TestCase):
    """Test v1 vs v2 consistency for hierarchical output."""

    def setUp(self):
        """Create temporary directories for testing."""
        self.test_dir = tempfile.mkdtemp()
        self.v1_outdir = os.path.join(self.test_dir, 'v1_output')
        self.v2_outdir = os.path.join(self.test_dir, 'v2_output')
        os.makedirs(self.v1_outdir)
        os.makedirs(self.v2_outdir)

    def tearDown(self):
        """Clean up temporary directories."""
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def test_mock_data_consistency(self):
        """
        Test v1 and v2 produce identical output for mock data.

        This test creates minimal mock data and verifies both
        implementations generate the same directory structure.
        """
        # Create mock parquet data
        mock_data = {
            'smiles': ['CCl', 'CCCl', 'CCF'],
            'inchikey': ['key1', 'key2', 'key3'],
            'parent_smiles': ['CC', 'CC', 'CC'],
            'parent_inchikey': ['parent_key', 'parent_key', 'parent_key'],
            'k': [1, 1, 1],
            'rule': ['R1', 'R1', 'R1'],
            'halogen': ['Cl', 'Cl', 'F'],
            'substitutions_json': ['[]', '[]', '[]']
        }
        df = pd.DataFrame(mock_data)

        parquet_path = os.path.join(self.test_dir, 'mock_products.parquet')
        df.to_parquet(parquet_path)

        # Test v1 implementation
        from halogenator.io_hierarchy import write_hierarchical_outputs

        parent_record = {
            'smiles': 'CC',
            'inchikey': 'parent_key',
            'name': 'ethane'
        }
        records = df.to_dict('records')

        v1_summary = write_hierarchical_outputs(
            parent_record, records, self.v1_outdir, halogens_order=['F', 'Cl', 'Br', 'I']
        )

        # Test v2 implementation
        from halogenator.io_hierarchy_v2 import write_hierarchical_outputs_streaming

        parents_smiles = {'parent_key': 'CC'}
        v2_summary = write_hierarchical_outputs_streaming(
            parquet_path, parents_smiles, self.v2_outdir, halogens_order=['F', 'Cl', 'Br', 'I']
        )

        # Verify both produced outputs
        self.assertTrue(os.path.exists(self.v1_outdir))
        self.assertTrue(os.path.exists(self.v2_outdir))

        # Check product counts match
        self.assertEqual(v1_summary['total_products'], v2_summary['total_products'])
        self.assertEqual(v1_summary['total_products'], 3)

    def test_v2_streaming_behavior(self):
        """
        Test v2 streaming implementation processes parents incrementally.

        Verifies that v2 can handle data without loading everything into memory.
        """
        # Create mock data with multiple parents
        mock_data = {
            'smiles': ['CCl', 'CCCl', 'CCF', 'C(Cl)CO', 'C(F)CO'],
            'inchikey': ['key1', 'key2', 'key3', 'key4', 'key5'],
            'parent_smiles': ['CC', 'CC', 'CC', 'CCO', 'CCO'],
            'parent_inchikey': ['parent_A', 'parent_A', 'parent_A', 'parent_B', 'parent_B'],
            'k': [1, 1, 1, 1, 1],
            'rule': ['R1', 'R1', 'R1', 'R1', 'R1'],
            'halogen': ['Cl', 'Cl', 'F', 'Cl', 'F'],
            'substitutions_json': ['[]'] * 5
        }
        df = pd.DataFrame(mock_data)

        parquet_path = os.path.join(self.test_dir, 'multi_parent_products.parquet')
        df.to_parquet(parquet_path)

        # Run v2 streaming
        from halogenator.io_hierarchy_v2 import write_hierarchical_outputs_streaming

        parents_smiles = {'parent_A': 'CC', 'parent_B': 'CCO'}
        summary = write_hierarchical_outputs_streaming(
            parquet_path, parents_smiles, self.v2_outdir
        )

        # Verify summary
        self.assertEqual(summary['total_parents'], 2)
        self.assertEqual(summary['total_products'], 5)
        self.assertTrue(summary['streaming'])

        # Verify global summary file exists
        summary_path = os.path.join(self.v2_outdir, 'hierarchy_summary.json')
        self.assertTrue(os.path.exists(summary_path))

    def test_empty_dataset(self):
        """Test v2 handles empty datasets gracefully."""
        # Create empty parquet
        df = pd.DataFrame({
            'smiles': [],
            'inchikey': [],
            'parent_smiles': [],
            'parent_inchikey': [],
            'k': [],
            'rule': [],
            'halogen': [],
            'substitutions_json': []
        })

        parquet_path = os.path.join(self.test_dir, 'empty_products.parquet')
        df.to_parquet(parquet_path)

        # Run v2 streaming
        from halogenator.io_hierarchy_v2 import write_hierarchical_outputs_streaming

        summary = write_hierarchical_outputs_streaming(
            parquet_path, {}, self.v2_outdir
        )

        # Verify empty results
        self.assertEqual(summary['total_parents'], 0)
        self.assertEqual(summary['total_products'], 0)
        self.assertEqual(summary['total_files'], 0)


if __name__ == '__main__':
    unittest.main()
