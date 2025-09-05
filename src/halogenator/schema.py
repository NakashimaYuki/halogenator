# -*- coding: ascii -*-
"""Schema validation for halogenator data structures."""

from typing import Set, List, Dict, Any, Union


# Required columns for products table
REQUIRED_PRODUCT_COLUMNS: Set[str] = {
    'product_smiles',
    'parent_smiles', 
    'parent_inchikey',
    'rule',
    'halogen',
    'parent_type',
    'depth'
}

# Optional columns that may be present  
OPTIONAL_PRODUCT_COLUMNS: Set[str] = {
    'sanitize_ok',
    'parent_name',
    'product_inchikey',
    'rule_description',
    'halogen_position'
}


def validate_products_schema(df_columns: Set[str]) -> None:
    """
    Validate that DataFrame contains all required product columns.
    
    Args:
        df_columns: Set of column names from DataFrame
        
    Raises:
        ValueError: If required columns are missing
    """
    missing = REQUIRED_PRODUCT_COLUMNS - df_columns
    if missing:
        raise ValueError(f"Products table missing required columns: {sorted(missing)}")


def validate_products_records(records: List[Dict[str, Any]]) -> None:
    """
    Validate that all product records contain required columns.
    
    Args:
        records: List of product records (dictionaries)
        
    Raises:
        ValueError: If required columns are missing from any records
    """
    if not records:
        return
    
    # Check all records for column consistency
    all_columns = set()
    record_columns = {}
    
    for i, record in enumerate(records):
        cols = set(record.keys())
        all_columns.update(cols)
        record_columns[i] = cols
    
    # Check for globally missing required columns
    missing_required = REQUIRED_PRODUCT_COLUMNS - all_columns
    if missing_required:
        raise ValueError(f"Products table missing required columns (not found in any record): {sorted(missing_required)}")
    
    # Check for columns that appear in some records but not others
    inconsistent_records = []
    for i, cols in record_columns.items():
        missing_in_record = REQUIRED_PRODUCT_COLUMNS - cols
        if missing_in_record:
            inconsistent_records.append((i, missing_in_record))
    
    if inconsistent_records:
        examples = inconsistent_records[:3]  # Show up to 3 examples
        error_msg = "Required columns missing in some records:\n"
        for record_idx, missing_cols in examples:
            error_msg += f"  Record {record_idx}: missing {sorted(missing_cols)}\n"
        if len(inconsistent_records) > 3:
            error_msg += f"  ... and {len(inconsistent_records) - 3} more records with missing columns"
        raise ValueError(error_msg.strip())


def get_all_known_columns() -> Set[str]:
    """Return all known (required + optional) product columns."""
    return REQUIRED_PRODUCT_COLUMNS | OPTIONAL_PRODUCT_COLUMNS