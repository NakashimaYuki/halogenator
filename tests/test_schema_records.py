# -*- coding: ascii -*-
"""Unit tests for schema validation on record collections."""

import unittest
from halogenator.schema import validate_products_records, REQUIRED_PRODUCT_COLUMNS


class TestSchemaRecords(unittest.TestCase):
    """Test validate_products_records functionality."""
    
    def test_validate_products_records_empty_list(self):
        """Test validation passes with empty record list."""
        # Should not raise any exception
        validate_products_records([])
    
    def test_validate_products_records_all_columns_present(self):
        """Test validation passes when all required columns present in all records."""
        records = [
            {col: f'value1_{col}' for col in REQUIRED_PRODUCT_COLUMNS},
            {col: f'value2_{col}' for col in REQUIRED_PRODUCT_COLUMNS}
        ]
        
        # Should not raise any exception
        validate_products_records(records)
    
    def test_validate_products_records_with_extra_columns(self):
        """Test validation passes with extra optional columns."""
        records = [
            {**{col: f'value1_{col}' for col in REQUIRED_PRODUCT_COLUMNS}, 
             'extra_col': 'extra_value', 'sanitize_ok': True},
            {**{col: f'value2_{col}' for col in REQUIRED_PRODUCT_COLUMNS},
             'another_extra': 42}
        ]
        
        # Should not raise any exception
        validate_products_records(records)
    
    def test_validate_products_records_missing_globally(self):
        """Test validation fails when required columns missing from ALL records."""
        # Missing 'rule' and 'halogen' from all records
        incomplete_columns = REQUIRED_PRODUCT_COLUMNS - {'rule', 'halogen'}
        records = [
            {col: f'value1_{col}' for col in incomplete_columns},
            {col: f'value2_{col}' for col in incomplete_columns}
        ]
        
        with self.assertRaises(ValueError) as context:
            validate_products_records(records)
        
        error_msg = str(context.exception).lower()
        self.assertIn("missing required columns", error_msg)
        self.assertIn("not found in any record", error_msg)
        self.assertIn("halogen", str(context.exception))
        self.assertIn("rule", str(context.exception))
    
    def test_validate_products_records_missing_in_some_records(self):
        """Test validation fails when columns missing from some records."""
        # First record complete, second record missing 'rule' and 'parent_type'
        records = [
            {col: f'value1_{col}' for col in REQUIRED_PRODUCT_COLUMNS},
            {col: f'value2_{col}' for col in REQUIRED_PRODUCT_COLUMNS - {'rule', 'parent_type'}}
        ]
        
        with self.assertRaises(ValueError) as context:
            validate_products_records(records)
        
        error_msg = str(context.exception)
        self.assertIn("Required columns missing in some records", error_msg)
        self.assertIn("Record 1:", error_msg)  # 0-indexed, so second record is index 1
        self.assertIn("parent_type", error_msg)
        self.assertIn("rule", error_msg)
    
    def test_validate_products_records_multiple_inconsistent_records(self):
        """Test validation reports multiple inconsistent records with examples."""
        # Multiple records with different missing columns
        records = [
            {col: f'value0_{col}' for col in REQUIRED_PRODUCT_COLUMNS},  # Complete
            {col: f'value1_{col}' for col in REQUIRED_PRODUCT_COLUMNS - {'rule'}},  # Missing rule
            {col: f'value2_{col}' for col in REQUIRED_PRODUCT_COLUMNS - {'halogen', 'depth'}},  # Missing halogen, depth
            {col: f'value3_{col}' for col in REQUIRED_PRODUCT_COLUMNS - {'parent_type'}},  # Missing parent_type
            {col: f'value4_{col}' for col in REQUIRED_PRODUCT_COLUMNS - {'rule'}}   # Another missing rule
        ]
        
        with self.assertRaises(ValueError) as context:
            validate_products_records(records)
        
        error_msg = str(context.exception)
        self.assertIn("Required columns missing in some records", error_msg)
        # Should show first 3 examples
        self.assertIn("Record 1:", error_msg)
        self.assertIn("Record 2:", error_msg) 
        self.assertIn("Record 3:", error_msg)
        # Should indicate more records with issues
        self.assertIn("1 more records", error_msg)
    
    def test_validate_products_records_single_inconsistent_record(self):
        """Test clear error message for single inconsistent record."""
        records = [
            {col: f'value0_{col}' for col in REQUIRED_PRODUCT_COLUMNS},
            {col: f'value1_{col}' for col in REQUIRED_PRODUCT_COLUMNS - {'rule', 'halogen'}}
        ]
        
        with self.assertRaises(ValueError) as context:
            validate_products_records(records)
        
        error_msg = str(context.exception)
        self.assertIn("Record 1: missing ['halogen', 'rule']", error_msg)
        self.assertNotIn("more records", error_msg)  # Only one inconsistent record
    
    def test_validate_products_records_realistic_scenario(self):
        """Test with realistic product record structure."""
        valid_records = [
            {
                'product_smiles': 'CCCl',
                'parent_smiles': 'CCC',
                'parent_inchikey': 'ATUOYWHBWRKTHZ-UHFFFAOYSA-N',
                'rule': 'R1',
                'halogen': 'Cl', 
                'parent_type': 'flavonoid',
                'depth': 1,
                'sanitize_ok': True
            },
            {
                'product_smiles': 'CCF',
                'parent_smiles': 'CCC', 
                'parent_inchikey': 'ATUOYWHBWRKTHZ-UHFFFAOYSA-N',
                'rule': 'R1',
                'halogen': 'F',
                'parent_type': 'flavonoid', 
                'depth': 1
            }
        ]
        
        # Should pass
        validate_products_records(valid_records)
        
        # Now test with missing column in second record
        invalid_records = valid_records.copy()
        del invalid_records[1]['parent_inchikey']
        
        with self.assertRaises(ValueError) as context:
            validate_products_records(invalid_records)
        
        self.assertIn("parent_inchikey", str(context.exception))
    
    def test_validate_products_records_error_message_globally_missing(self):
        """Test error message format for globally missing columns."""
        # Missing 'rule' and 'halogen' from all records
        incomplete_columns = REQUIRED_PRODUCT_COLUMNS - {'rule', 'halogen'}
        records = [
            {col: f'value1_{col}' for col in incomplete_columns},
            {col: f'value2_{col}' for col in incomplete_columns}
        ]
        
        with self.assertRaisesRegex(ValueError, r"Products table.*missing required columns.*not found in any record.*halogen.*rule"):
            validate_products_records(records)
    
    def test_validate_products_records_error_message_record_specific(self):
        """Test error message format for record-specific missing columns."""
        # First record complete, second missing 'parent_smiles'
        records = [
            {col: f'value1_{col}' for col in REQUIRED_PRODUCT_COLUMNS},
            {col: f'value2_{col}' for col in REQUIRED_PRODUCT_COLUMNS - {'parent_smiles'}}
        ]
        
        with self.assertRaisesRegex(ValueError, r"Required columns missing in some records[\s\S]*Record 1[\s\S]*parent_smiles"):
            validate_products_records(records)
    
    def test_validate_products_records_error_message_multiple_records(self):
        """Test error message format shows multiple record examples."""
        # Three records with different issues
        records = [
            {col: f'value0_{col}' for col in REQUIRED_PRODUCT_COLUMNS},  # Complete
            {col: f'value1_{col}' for col in REQUIRED_PRODUCT_COLUMNS - {'rule'}},  # Missing rule
            {col: f'value2_{col}' for col in REQUIRED_PRODUCT_COLUMNS - {'halogen'}},  # Missing halogen
        ]
        
        with self.assertRaisesRegex(ValueError, r"Record 1[\s\S]*rule[\s\S]*Record 2[\s\S]*halogen"):
            validate_products_records(records)


if __name__ == '__main__':
    unittest.main()