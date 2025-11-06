# -*- coding: ascii -*-
"""Tests to ensure unified parent key consistency between production and consumption sides."""

import unittest
from src.halogenator.report import make_unified_parent_key, record_to_unified_parent_key


class TestUnifiedParentKeyConsistency(unittest.TestCase):
    """Test that unified parent key generation is consistent and prevents statistical errors."""
    
    def test_make_unified_parent_key_contract(self):
        """Test that make_unified_parent_key follows proper prefixing rules."""
        # Valid SMILES should get IK: prefix (when RDKit available)
        result = make_unified_parent_key("CCO")
        self.assertTrue(result.startswith("IK:") or result.startswith("SMI:"))
        self.assertNotEqual(result, "SMI:UNKNOWN")
        
        # Invalid inputs should return SMI:UNKNOWN
        self.assertEqual(make_unified_parent_key(""), "SMI:UNKNOWN")
        self.assertEqual(make_unified_parent_key(None), "SMI:UNKNOWN")
        self.assertEqual(make_unified_parent_key("Unknown"), "SMI:UNKNOWN")
        self.assertEqual(make_unified_parent_key("UNKNOWN"), "SMI:UNKNOWN")
        self.assertEqual(make_unified_parent_key("unknown"), "SMI:UNKNOWN")
        
    def test_record_to_unified_parent_key_contract(self):
        """Test that record_to_unified_parent_key returns proper types."""
        # Valid InChIKey should return IK: prefixed string
        result = record_to_unified_parent_key({
            'parent_inchikey': 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'
        })
        self.assertIsInstance(result, str)
        self.assertTrue(result.startswith("IK:"))
        
        # Valid SMILES should return string (IK: or SMI: prefix)
        result = record_to_unified_parent_key({
            'parent_smiles': 'CCO'
        })
        self.assertIsInstance(result, str)
        self.assertTrue(result.startswith("IK:") or result.startswith("SMI:"))
        
        # Invalid/missing data should return None
        self.assertIsNone(record_to_unified_parent_key({}))
        self.assertIsNone(record_to_unified_parent_key({
            'parent_inchikey': 'Unknown',
            'parent_smiles': 'Unknown'
        }))
        
    def test_no_incorrect_ik_prefixes(self):
        """Test that SMILES are never incorrectly prefixed with IK:."""
        # This tests the critical fix for P0-A issue
        
        # Even if InChI generation fails, should not get IK:SMILES format
        test_cases = [
            "CCO",
            "CC(C)C", 
            "c1ccccc1",
            "C" * 50,  # Long SMILES
        ]
        
        for smiles in test_cases:
            result = make_unified_parent_key(smiles)
            # Result should either be proper IK:<inchikey> or SMI:<smiles>/<hash>
            if result.startswith("IK:"):
                # If it's IK: prefixed, the part after IK: should look like an InChIKey
                key_part = result[3:]  # Remove "IK:" prefix
                # InChIKeys have a specific format: 14 chars, hyphen, 10 chars, hyphen, 1 char
                self.assertRegex(key_part, r'^[A-Z]{14}-[A-Z]{9}[A-Z]-[A-Z]$', 
                    f"Invalid InChIKey format in result: {result}")
            elif result.startswith("SMI:"):
                # SMI: prefix is valid for SMILES or hash
                self.assertTrue(len(result) > 4)  # Should have content after SMI:
            else:
                self.fail(f"Result should start with IK: or SMI:, got: {result}")
                
    def test_parent_list_filtering_simulation(self):
        """Test that SMI:UNKNOWN keys are properly filtered from parent lists."""
        # Simulate parent list processing
        parent_pairs = [
            ("CCO", "ethanol"),
            ("", "empty_smiles"),
            ("Unknown", "unknown_molecule"), 
            ("CC(C)C", "isobutane"),
            (None, "null_smiles")
        ]
        
        # Simulate the filtering logic that should be applied
        parent_key_to_name = {}
        all_parent_keys = set()
        
        for smiles, name in parent_pairs:
            unified_key = make_unified_parent_key(smiles)
            if unified_key and unified_key != 'SMI:UNKNOWN':
                parent_key_to_name[unified_key] = name
                all_parent_keys.add(unified_key)
        
        # Verify that SMI:UNKNOWN keys are filtered out
        self.assertNotIn("SMI:UNKNOWN", all_parent_keys)
        self.assertNotIn("SMI:UNKNOWN", parent_key_to_name.keys())
        
        # Verify that valid keys are included
        self.assertTrue(len(all_parent_keys) >= 2)  # At least ethanol and isobutane
        
        # Simulate product record processing
        product_records = [
            {"parent_inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"},  # ethanol
            {"parent_smiles": "CC(C)C"},  # isobutane
            {"parent_smiles": "Unknown"},  # should return None and be filtered
        ]
        
        parents_with_products = set()
        for record in product_records:
            parent_key = record_to_unified_parent_key(record)
            if parent_key:  # This filters out None results
                parents_with_products.add(parent_key)
        
        # Verify consistency: parents_with_products should be subset of all_parent_keys
        self.assertTrue(parents_with_products.issubset(all_parent_keys),
            f"parents_with_products {parents_with_products} should be subset of all_parent_keys {all_parent_keys}")
        
        # Calculate parents without products (this should not include SMI:UNKNOWN)
        parents_without_products = all_parent_keys - parents_with_products
        self.assertNotIn("SMI:UNKNOWN", parents_without_products)
        
    def test_case_insensitive_unknown_handling(self):
        """Test that various cases of 'unknown' are handled consistently."""
        unknown_variations = ["Unknown", "UNKNOWN", "unknown", "UnKnOwN"]
        
        for variation in unknown_variations:
            # make_unified_parent_key should return SMI:UNKNOWN
            result = make_unified_parent_key(variation)
            self.assertEqual(result, "SMI:UNKNOWN", f"Failed for variation: {variation}")
            
            # record_to_unified_parent_key should return None
            result = record_to_unified_parent_key({"parent_inchikey": variation})
            self.assertIsNone(result, f"Failed for InChIKey variation: {variation}")
            
            result = record_to_unified_parent_key({"parent_smiles": variation})
            self.assertIsNone(result, f"Failed for SMILES variation: {variation}")
            
    def test_long_smiles_hash_protection(self):
        """Test that very long SMILES are properly hashed."""
        # Create a SMILES longer than 120 characters
        long_smiles = "C" * 150
        result = make_unified_parent_key(long_smiles)
        
        if result.startswith("SMI:"):
            # Should be hashed to 16 characters
            hash_part = result[4:]  # Remove "SMI:" prefix
            self.assertEqual(len(hash_part), 16, "Long SMILES should be hashed to 16 characters")
            # Hash should be hex characters
            self.assertRegex(hash_part, r'^[0-9a-f]{16}$', "Hash should be 16 hex characters")


if __name__ == '__main__':
    unittest.main()