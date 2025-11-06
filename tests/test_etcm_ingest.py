# -*- coding: ascii -*-
"""Test ETCM ingest functionality."""

import unittest
import tempfile
import os
import pandas as pd
from scripts.etcm_to_parents import main as etcm_main, ascii_slug
import sys


class TestETCMIngest(unittest.TestCase):
    """Test ETCM ingest conversion."""
    
    def setUp(self):
        """Create synthetic ETCM CSV data for testing."""
        # Create test CSV with Chinese column names (using unicode escapes) and mixed content
        self.test_data = {
            'Ingredient ID': ['TCMIP-I-00001', 'TCMIP-I-00002', 'TCMIP-I-00003'],
            'Ingredient Name in English': ['Quercetin', '5,7-Dihydroxyflavone', 'Apigenin 7-glucoside'],
            '\u6210\u5206\u4e2d\u6587\u540d': ['\u69f2\u76ae\u7d20', '\u4e8c\u7f9f\u57fa\u9ec4\u916e', '\u8299\u83dc\u7d20'],
            '\u6807\u51c6SMILES': [
                'O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12',
                'O=c1cc(-c2ccccc2)oc2cc(O)cc(O)c12',
                'O=c1cc(-c2ccc(O)cc2)oc2cc(O)cc(O)c12'
            ],
            'is_flavonoid': [True, True, True]
        }
        
    def test_ascii_slug_function(self):
        """Test ASCII slug generation."""
        # Test basic conversion
        self.assertEqual(ascii_slug('Quercetin'), 'Quercetin')
        self.assertEqual(ascii_slug('5,7-Dihydroxyflavone'), '5_7-Dihydroxyflavone')
        self.assertEqual(ascii_slug('Apigenin 7-glucoside'), 'Apigenin_7-glucoside')
        
        # Test Chinese characters removal (using unicode escape)
        chinese_text = '\u69f2\u76ae\u7d20'  # Chinese characters
        result = ascii_slug(chinese_text)
        self.assertEqual(result, 'unknown')
        
        # Test None and empty
        self.assertEqual(ascii_slug(None), 'unknown')
        self.assertEqual(ascii_slug(''), 'unknown')
        
        # Test long string truncation
        long_string = 'A' * 100
        result = ascii_slug(long_string)
        self.assertEqual(len(result), 80)
        
    def test_csv_conversion_minimal(self):
        """Test basic CSV conversion with synthetic data."""
        # Create temporary CSV file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False, encoding='utf-8') as tmp_csv:
            df = pd.DataFrame(self.test_data)
            df.to_csv(tmp_csv.name, index=False)
            csv_path = tmp_csv.name
        
        # Create temporary output files
        with tempfile.NamedTemporaryFile(suffix='.smi', delete=False) as tmp_smi:
            smi_path = tmp_smi.name
        with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as tmp_meta:
            meta_path = tmp_meta.name
            
        try:
            # Prepare sys.argv for etcm_main
            old_argv = sys.argv
            sys.argv = [
                'etcm_to_parents.py',
                '-i', csv_path,
                '-o', smi_path,
                '--out-meta', meta_path,
                '--csv'
            ]
            
            # Run conversion
            etcm_main()
            
            # Verify outputs exist
            self.assertTrue(os.path.exists(smi_path), "SMI file should be created")
            self.assertTrue(os.path.exists(meta_path), "Meta CSV should be created")
            
            # Check SMI file content
            with open(smi_path, 'r', encoding='ascii') as f:
                smi_lines = f.readlines()
            
            self.assertEqual(len(smi_lines), 3, "Should have 3 SMILES entries")
            
            # Check first line format: SMILES<TAB>ASCII_NAME
            first_line = smi_lines[0].strip()
            parts = first_line.split('\t')
            self.assertEqual(len(parts), 2, "Each line should have SMILES and name")
            
            smiles, name = parts
            self.assertTrue(smiles.startswith('O='), "Should be valid SMILES")
            self.assertTrue(name.startswith('TCMIP-I-00001__'), "Name should have ID prefix")
            self.assertIn('Quercetin', name, "Name should contain ASCII compound name")
            
            # Check meta CSV content  
            meta_df = pd.read_csv(meta_path)
            self.assertEqual(len(meta_df), 3, "Meta should have 3 rows")
            self.assertIn('_ascii_name', meta_df.columns, "Meta should have ascii_name column")
            
            # Check that ASCII names are actually ASCII-safe
            for ascii_name in meta_df['_ascii_name']:
                # Try to encode as ASCII - should not raise exception
                try:
                    ascii_name.encode('ascii')
                except UnicodeEncodeError:
                    self.fail(f"ASCII name contains non-ASCII characters: {ascii_name}")
            
        finally:
            # Restore sys.argv
            sys.argv = old_argv
            
            # Cleanup temp files
            for path in [csv_path, smi_path, meta_path]:
                if os.path.exists(path):
                    os.unlink(path)
                    
    def test_sample_parameter(self):
        """Test --sample parameter limits output."""
        # Create temporary CSV file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False, encoding='utf-8') as tmp_csv:
            df = pd.DataFrame(self.test_data)
            df.to_csv(tmp_csv.name, index=False)
            csv_path = tmp_csv.name
        
        # Create temporary output files
        with tempfile.NamedTemporaryFile(suffix='.smi', delete=False) as tmp_smi:
            smi_path = tmp_smi.name
        with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as tmp_meta:
            meta_path = tmp_meta.name
            
        try:
            # Prepare sys.argv for etcm_main with --sample 2
            old_argv = sys.argv
            sys.argv = [
                'etcm_to_parents.py',
                '-i', csv_path,
                '-o', smi_path,
                '--out-meta', meta_path,
                '--csv',
                '--sample', '2'
            ]
            
            # Run conversion
            etcm_main()
            
            # Check that only 2 entries were produced
            with open(smi_path, 'r', encoding='ascii') as f:
                smi_lines = f.readlines()
            
            self.assertEqual(len(smi_lines), 2, "Should have 2 SMILES entries with --sample 2")
            
        finally:
            # Restore sys.argv
            sys.argv = old_argv
            
            # Cleanup temp files
            for path in [csv_path, smi_path, meta_path]:
                if os.path.exists(path):
                    os.unlink(path)


if __name__ == '__main__':
    unittest.main()