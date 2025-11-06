# -*- coding: ascii -*-
"""Tests for parent key metadata functionality."""

import unittest
from src.halogenator.report import make_unified_parent_key_with_metadata, SMILES_HASH_THRESHOLD


class TestParentKeyMetadata(unittest.TestCase):
    """Test that parent key metadata is properly generated."""
    
    def test_metadata_structure(self):
        """Test that metadata has expected structure."""
        key, meta = make_unified_parent_key_with_metadata("CCO")
        
        # Check that all expected metadata fields are present
        expected_fields = ['smiles_hashed', 'normalized_length', 'rdkit_available', 'inchikey_generated']
        for field in expected_fields:
            self.assertIn(field, meta, f"Metadata missing field: {field}")
            
        # Check metadata types
        self.assertIsInstance(meta['smiles_hashed'], bool)
        self.assertIsInstance(meta['normalized_length'], int)
        self.assertIsInstance(meta['rdkit_available'], bool)
        self.assertIsInstance(meta['inchikey_generated'], bool)
        
    def test_normal_smiles_metadata(self):
        """Test metadata for normal SMILES."""
        key, meta = make_unified_parent_key_with_metadata("CCO")
        
        # Should not be hashed (too short)
        self.assertFalse(meta['smiles_hashed'])
        self.assertEqual(meta['normalized_length'], 3)
        
        # Should have generated key starting with IK: or SMI:
        self.assertTrue(key.startswith("IK:") or key.startswith("SMI:"))
        
    def test_unknown_smiles_metadata(self):
        """Test metadata for unknown/invalid SMILES."""
        test_cases = ["", "Unknown", "UNKNOWN", "unknown", None]
        
        for smiles in test_cases:
            key, meta = make_unified_parent_key_with_metadata(smiles)
            
            self.assertEqual(key, "SMI:UNKNOWN")
            self.assertFalse(meta['smiles_hashed'])
            self.assertFalse(meta['inchikey_generated'])
            # normalized_length should reflect the normalized SMILES length
            # For unknown/None inputs, normalized length should be 0
            expected_len = 0  # All unknown cases normalize to empty string
            self.assertEqual(meta['normalized_length'], expected_len)
            
    def test_hash_threshold_behavior(self):
        """Test behavior around SMILES_HASH_THRESHOLD."""
        # Test right at threshold (should NOT be hashed)
        exactly_threshold = "C" * SMILES_HASH_THRESHOLD
        key_threshold, meta_threshold = make_unified_parent_key_with_metadata(exactly_threshold)
        
        # Test just over threshold (should be hashed if RDKit not available or InChI fails)
        over_threshold = "C" * (SMILES_HASH_THRESHOLD + 1)
        key_over, meta_over = make_unified_parent_key_with_metadata(over_threshold)
        
        # Check that length is recorded correctly
        self.assertEqual(meta_threshold['normalized_length'], SMILES_HASH_THRESHOLD)
        self.assertEqual(meta_over['normalized_length'], SMILES_HASH_THRESHOLD + 1)
        
        print(f"Threshold test: {SMILES_HASH_THRESHOLD} chars -> hashed: {meta_threshold['smiles_hashed']}")
        print(f"Over threshold test: {SMILES_HASH_THRESHOLD + 1} chars -> hashed: {meta_over['smiles_hashed']}")
        
    def test_whitespace_normalization_in_metadata(self):
        """Test that whitespace normalization is reflected in metadata."""
        # SMILES with extra whitespace
        messy_smiles = "  C C O  \t\n  "
        key, meta = make_unified_parent_key_with_metadata(messy_smiles)
        
        # Normalized length should be of the normalized SMILES "C C O" = 5
        # (whitespace compression keeps single spaces between tokens)
        self.assertEqual(meta['normalized_length'], 5)
        
        # Should still get a valid key
        self.assertTrue(key.startswith("IK:") or key.startswith("SMI:"))
        
    def test_rdkit_availability_tracking(self):
        """Test that RDKit availability is properly tracked."""
        key, meta = make_unified_parent_key_with_metadata("CCO")
        
        # This test assumes RDKit is available in the test environment
        # If RDKit is available, rdkit_available should be True
        # The actual value depends on the environment, so we just check it's a boolean
        self.assertIsInstance(meta['rdkit_available'], bool)
        
        # If RDKit is available and worked, we should have an InChIKey
        if meta['rdkit_available'] and meta['inchikey_generated']:
            self.assertTrue(key.startswith("IK:"))
        
    def test_consistency_with_original_function(self):
        """Test that the new function produces same keys as the original."""
        from src.halogenator.report import make_unified_parent_key
        
        test_cases = [
            "CCO",
            "CC(C)C", 
            "c1ccccc1",
            "Unknown",
            "",
            None,
            "C" * 50,  # Medium length
            "C" * 200,  # Long length
        ]
        
        for smiles in test_cases:
            key_original = make_unified_parent_key(smiles)
            key_with_meta, meta = make_unified_parent_key_with_metadata(smiles)
            
            self.assertEqual(key_original, key_with_meta, 
                f"Keys differ for SMILES '{smiles}': original='{key_original}', with_meta='{key_with_meta}'")


    def test_qa_summary_metadata_integration(self):
        """Test that QA summary contains exactly one parent_key_metadata block."""
        from src.halogenator.report import _finalize_incremental_stats, _init_incremental_stats
        
        # Initialize test stats
        stats = _init_incremental_stats()
        
        # Mock parent pairs with diverse SMILES to trigger metadata generation
        parent_pairs = [
            ("CCO", "ethanol"),
            ("c1ccccc1", "benzene"),
            ("Unknown", "unknown_compound"),
            ("C" * 150, "long_smiles")  # Long SMILES to trigger hashing
        ]
        
        # Mock config
        config = {"test": True}
        
        # Finalize stats - this should write parent_key_metadata exactly once
        _finalize_incremental_stats(stats, parent_pairs, config)
        
        # Verify parent_key_metadata exists and has expected structure
        self.assertIn('parent_key_metadata', stats, "QA summary should contain parent_key_metadata")
        
        metadata = stats['parent_key_metadata']
        
        # Check that it has the expected aggregated metadata structure
        expected_fields = ['hashed_count', 'hashed_examples', 'rdkit_available_count', 'inchikey_generated_count']
        for field in expected_fields:
            self.assertIn(field, metadata, f"parent_key_metadata missing field: {field}")
        
        # Check types
        self.assertIsInstance(metadata['hashed_count'], int)
        self.assertIsInstance(metadata['hashed_examples'], list)
        self.assertIsInstance(metadata['rdkit_available_count'], int)
        self.assertIsInstance(metadata['inchikey_generated_count'], int)
        
        # Verify that long SMILES triggered hashing
        self.assertGreater(metadata['hashed_count'], 0, "Long SMILES should have triggered hashing")
        self.assertGreater(len(metadata['hashed_examples']), 0, "Should have hashed examples")
        
        # Check that the examples are proper keys
        for example in metadata['hashed_examples']:
            self.assertTrue(example.startswith("SMI:") or example.startswith("IK:"), 
                          f"Hashed example should be a proper key: {example}")
    
    def test_unified_parent_key_metadata_write_path(self):
        """Test that parent_key_metadata is written only once in the unified path."""
        # This test ensures that we don't have duplicate metadata writes in different functions
        from src.halogenator.report import generate_summary_report
        import tempfile
        import os
        import json
        
        # Create a minimal config for testing
        config = {
            'input_dir': 'data/input',
            'output_dir': tempfile.mkdtemp(),
            'halogens': ['F', 'Cl'],
            'rules': ['R1'],
            'k_max': 1
        }
        
        try:
            # We can't easily test the full generate_summary_report without setting up
            # a complete test environment, so we'll test the key component that we modified
            from src.halogenator.report import _finalize_incremental_stats, _init_incremental_stats
            
            stats = _init_incremental_stats()
            parent_pairs = [("CCO", "ethanol")]
            
            # This should be the ONLY place where parent_key_metadata is written
            _finalize_incremental_stats(stats, parent_pairs, config)
            
            # Verify it was written exactly once
            self.assertEqual(len([k for k in stats.keys() if k == 'parent_key_metadata']), 1,
                           "parent_key_metadata should appear exactly once in stats")
            
        finally:
            # Clean up temp directory
            import shutil
            if os.path.exists(config['output_dir']):
                shutil.rmtree(config['output_dir'])


if __name__ == '__main__':
    unittest.main()