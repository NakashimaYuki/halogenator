# -*- coding: ascii -*-
"""Unit tests for schema validation."""

import unittest
import pandas as pd
from halogenator.schema import validate_products_schema, REQUIRED_PRODUCT_COLUMNS


class TestSchema(unittest.TestCase):
    """Test schema validation functionality."""
    
    def test_validate_products_schema_all_columns_present(self):
        """Test validation passes when all required columns present."""
        # All required columns
        columns = REQUIRED_PRODUCT_COLUMNS.copy()
        
        # Should not raise any exception
        validate_products_schema(columns)
    
    def test_validate_products_schema_with_extra_columns(self):
        """Test validation passes with extra optional columns."""
        # Required + some optional columns
        columns = REQUIRED_PRODUCT_COLUMNS.copy()
        columns.update({'sanitize_ok', 'product_inchikey', 'extra_column'})
        
        # Should not raise any exception 
        validate_products_schema(columns)
    
    def test_validate_products_schema_missing_single_column(self):
        """Test validation fails when single required column missing."""
        # Missing 'rule' column
        columns = REQUIRED_PRODUCT_COLUMNS.copy()
        columns.remove('rule')
        
        with self.assertRaises(ValueError) as context:
            validate_products_schema(columns)
        
        self.assertIn("missing required columns", str(context.exception).lower())
        self.assertIn("rule", str(context.exception))
    
    def test_validate_products_schema_missing_multiple_columns(self):
        """Test validation fails when multiple required columns missing."""
        # Missing multiple columns
        columns = REQUIRED_PRODUCT_COLUMNS.copy()
        columns.remove('rule')
        columns.remove('halogen') 
        columns.remove('parent_type')
        
        with self.assertRaises(ValueError) as context:
            validate_products_schema(columns)
        
        error_msg = str(context.exception).lower()
        self.assertIn("missing required columns", error_msg)
        # Should list missing columns in sorted order
        self.assertIn("halogen", str(context.exception))
        self.assertIn("parent_type", str(context.exception))
        self.assertIn("rule", str(context.exception))
    
    def test_validate_products_schema_empty_columns(self):
        """Test validation fails with no columns."""
        columns = set()
        
        with self.assertRaises(ValueError) as context:
            validate_products_schema(columns)
        
        error_msg = str(context.exception).lower()
        self.assertIn("missing required columns", error_msg)
        # Should list all required columns
        for col in REQUIRED_PRODUCT_COLUMNS:
            self.assertIn(col, str(context.exception))
    
    def test_validate_products_schema_dataframe_scenario(self):
        """Test validation with DataFrame columns simulation."""
        # Simulate a valid DataFrame with columns
        valid_df_columns = {
            'product_smiles', 'parent_smiles', 'parent_inchikey',
            'rule', 'halogen', 'parent_type', 'depth', 'sanitize_ok'
        }
        
        # Should pass
        validate_products_schema(valid_df_columns)
        
        # Simulate invalid DataFrame missing key columns
        invalid_df_columns = {
            'product_smiles', 'parent_smiles', 'sanitize_ok'
        }
        
        with self.assertRaises(ValueError):
            validate_products_schema(invalid_df_columns)


if __name__ == '__main__':
    unittest.main()