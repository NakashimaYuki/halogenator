# -*- coding: ascii -*-
"""Schema validation for halogenator data structures."""

from typing import Set, List, Dict, Any, Union, Final, Tuple


# P0 legacy required columns for backwards compatibility
REQUIRED_P0_COLUMNS: Set[str] = {
    'product_smiles',
    'parent_smiles', 
    'parent_inchikey',
    'rule',
    'halogen',
    'parent_type',
    'depth'
}

# P1 required columns for k-dimensional enumeration
REQUIRED_P1_COLUMNS: Set[str] = {
    'smiles',           # Canonical product SMILES
    'inchikey',         # Product InChI key for deduplication
    'parent_smiles',    # Parent SMILES
    'parent_inchikey',  # Parent InChI key
    'k',                # Substitution depth (replaces 'depth')
    'rule',             # Last applied rule (for compatibility)
    'halogen',          # Last applied halogen (for compatibility)
    'substitutions',    # Complete JSON history
    'constraints_ok',   # Boolean constraint validation
    'constraints_violations'  # JSON constraint violation details
}

# Optional columns that may be present in both formats
OPTIONAL_PRODUCT_COLUMNS: Set[str] = {
    'sanitize_ok',
    'parent_name',
    'parent_type',
    'product_smiles',    # P0 compatibility
    'product_inchikey',  # P0 compatibility
    'rule_description',
    'halogen_position',
    'depth',             # P0 compatibility
    'site_index',        # Site information
    'sym_class',         # Symmetry class
    'ring_tag',          # Ring tag information
    'pains_flags',       # PAINS filter results
    'rule_family'        # Rule family for grouping (R2a/R2b -> R2)
}

# Backward compatibility alias
REQUIRED_PRODUCT_COLUMNS = REQUIRED_P0_COLUMNS  # backward compat

# Rule constants for centralized configuration management
ALL_RULES: Final[Tuple[str, ...]] = ('R1', 'R2', 'R3', 'R4', 'R5', 'R6')

# Halogen constants for centralized configuration management
ALL_HALOGENS: Final[Tuple[str, ...]] = ('F', 'Cl', 'Br', 'I')

# Overview counters that appear only in total section (not distributed to granular dimensions)
OVERVIEW_COUNTERS: Final[Tuple[str, ...]] = ('attempts', 'products', 'dedup_hits_statesig', 'dedup_hits_inchi')

# QA path constants - now sourced from qa_utils for single source of truth
# CONSOLIDATION NOTE: This import eliminates dual-source-of-truth between schema.py and qa_utils.py
# All canonical event constants should be defined in qa_utils.py going forward
from .qa_utils import STANDARD_QA_PATHS_KEYS

# Backward compatibility: maintain QA_PATH_KEYS as alias, including schema-specific keys
# Note: qa_utils constants are the authoritative source for new development
QA_PATH_KEYS = tuple(STANDARD_QA_PATHS_KEYS + [
    'carbonyl_unknown',           # Schema-specific legacy key (used in tests)
    'pruned_inchikey_dupe',       # Legacy alias maintained for backward compatibility
    'sugar_proximity_filtered',   # Non-pivot path event (used in tests)
])


def ensure_full_schema_qa_paths(qa_paths: Dict[str, int], emit_legacy_keys: bool = False) -> Dict[str, int]:
    """
    Ensure qa_paths has full schema key set with consistent initialization.

    This function unifies streaming vs snapshot consistency by ensuring both paths
    have identical key sets, eliminating the Gate-1 key set mismatch issue.

    Args:
        qa_paths: Input QA paths dictionary (may have partial keys)
        emit_legacy_keys: If True, emit legacy field aliases for backward compatibility

    Returns:
        QA paths dictionary with full schema keys and compatibility aliases
    """
    # Start with full schema (all keys initialized to 0)
    base = empty_qa_paths()

    # Merge existing counts
    for k, v in (qa_paths or {}).items():
        base[k] = int(v)

    # Apply compatibility layer (legacy key handling)
    return ensure_qa_paths_compatibility(base, emit_legacy_keys=emit_legacy_keys)


def ensure_qa_paths_compatibility(qa_paths: Dict[str, int], emit_legacy_keys: bool = False) -> Dict[str, int]:
    """
    Ensure backward compatibility for QA paths by providing alias mapping.

    Args:
        qa_paths: Original QA paths dictionary
        emit_legacy_keys: If True, emit legacy field aliases for backward compatibility

    Returns:
        QA paths dictionary with compatibility aliases
    """
    # Start with a copy to avoid modifying the original
    result = dict(qa_paths or {})

    # Forward compatibility: if only legacy key exists, copy to new key
    if 'pruned_inchikey_dupe' in result and 'dedup_hits_inchi' not in result:
        result['dedup_hits_inchi'] = result['pruned_inchikey_dupe']

    # Backward compatibility: optionally emit legacy keys when new keys exist
    if emit_legacy_keys:
        # Emit legacy aliases for backward compatibility
        if 'dedup_hits_inchi' in result:
            result['pruned_inchikey_dupe'] = result['dedup_hits_inchi']
    else:
        # Default mode: remove legacy keys if both exist, prefer new keys
        if 'dedup_hits_inchi' in result and 'pruned_inchikey_dupe' in result:
            del result['pruned_inchikey_dupe']

    return result


def empty_qa_paths() -> Dict[str, int]:
    """
    Create empty QA paths dictionary with all known keys initialized to 0.
    
    Note: This factory function focuses only on qa_paths centralization to avoid 
    schema inconsistencies. Full QA stats structures are context-specific and 
    created inline where needed.
    """
    return {k: 0 for k in QA_PATH_KEYS}




def detect_schema_format(columns: Set[str]) -> str:
    """
    Detect whether data uses P0 or P1 schema format.
    
    Args:
        columns: Set of column names
        
    Returns:
        'P0' or 'P1' format identifier
    """
    # P1 format has 'k' and 'substitutions', P0 has 'depth' and 'product_smiles'
    has_p1_markers = 'k' in columns and 'substitutions' in columns
    has_p0_markers = 'depth' in columns and 'product_smiles' in columns
    
    if has_p1_markers:
        return 'P1'
    elif has_p0_markers:
        return 'P0'
    else:
        # Default to P0 for backwards compatibility
        return 'P0'


def validate_products_schema(df_columns: Set[str]) -> None:
    """
    Validate that DataFrame contains all required product columns.
    Supports both P0 (legacy) and P1 (k-dimensional) formats.
    
    Args:
        df_columns: Set of column names from DataFrame
        
    Raises:
        ValueError: If required columns are missing
    """
    schema_format = detect_schema_format(df_columns)
    
    if schema_format == 'P1':
        required_columns = REQUIRED_P1_COLUMNS
    else:
        required_columns = REQUIRED_P0_COLUMNS
    
    missing = required_columns - df_columns
    if missing:
        raise ValueError(f"Products table ({schema_format} format) missing required columns: {sorted(missing)}")


def validate_products_records(records: List[Dict[str, Any]]) -> None:
    """
    Validate that all product records contain required columns.
    Supports both P0 (legacy) and P1 (k-dimensional) formats.
    
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
    
    # Detect schema format
    schema_format = detect_schema_format(all_columns)
    required_columns = REQUIRED_P1_COLUMNS if schema_format == 'P1' else REQUIRED_P0_COLUMNS
    
    # Check for globally missing required columns
    missing_required = required_columns - all_columns
    if missing_required:
        raise ValueError(f"Products table ({schema_format} format) missing required columns (not found in any record): {sorted(missing_required)}")
    
    # Check for columns that appear in some records but not others
    inconsistent_records = []
    for i, cols in record_columns.items():
        missing_in_record = required_columns - cols
        if missing_in_record:
            inconsistent_records.append((i, missing_in_record))
    
    if inconsistent_records:
        examples = inconsistent_records[:3]  # Show up to 3 examples
        error_msg = f"Required columns missing in some records ({schema_format} format):\n"
        for record_idx, missing_cols in examples:
            error_msg += f"  Record {record_idx}: missing {sorted(missing_cols)}\n"
        if len(inconsistent_records) > 3:
            error_msg += f"  ... and {len(inconsistent_records) - 3} more records with missing columns"
        raise ValueError(error_msg.strip())


def get_all_known_columns() -> Set[str]:
    """Return all known (required + optional) product columns from both P0 and P1 formats."""
    return REQUIRED_P0_COLUMNS | REQUIRED_P1_COLUMNS | OPTIONAL_PRODUCT_COLUMNS


def get_required_columns(schema_format: str) -> Set[str]:
    """
    Get required columns for specific schema format.
    
    Args:
        schema_format: 'P0' or 'P1'
        
    Returns:
        Set of required column names
    """
    if schema_format == 'P1':
        return REQUIRED_P1_COLUMNS
    else:
        return REQUIRED_P0_COLUMNS


def validate_p1_constraints_data(record: Dict[str, Any]) -> None:
    """
    Validate P1-specific constraint data integrity.
    
    Args:
        record: Product record dictionary
        
    Raises:
        ValueError: If constraint data is malformed
    """
    if 'constraints_ok' in record:
        constraints_ok = record['constraints_ok']
        violations = record.get('constraints_violations', {})
        
        # Basic type checking
        if not isinstance(constraints_ok, bool):
            raise ValueError("constraints_ok must be boolean")
        
        if violations is not None and not isinstance(violations, dict):
            raise ValueError("constraints_violations must be dict or None")
        
        # Logic consistency check
        if constraints_ok and violations:
            raise ValueError("constraints_ok=True but violations present")
        
        if not constraints_ok and not violations:
            # This is acceptable - constraint failed but no specific violations recorded
            pass


def validate_p1_substitutions_data(record: Dict[str, Any]) -> None:
    """
    Validate P1 substitutions history data.
    
    Args:
        record: Product record dictionary
        
    Raises:
        ValueError: If substitutions data is malformed
    """
    if 'substitutions' in record:
        substitutions = record['substitutions']
        
        if not isinstance(substitutions, list):
            raise ValueError("substitutions must be a list")
        
        k_value = record.get('k', 0)
        if len(substitutions) != k_value:
            raise ValueError(f"substitutions length {len(substitutions)} doesn't match k={k_value}")
        
        # Validate each substitution item
        for i, sub in enumerate(substitutions):
            if not isinstance(sub, dict):
                raise ValueError(f"substitutions[{i}] must be a dict")
            
            required_sub_fields = {'rule', 'halogen', 'depth'}
            missing_sub_fields = required_sub_fields - set(sub.keys())
            if missing_sub_fields:
                raise ValueError(f"substitutions[{i}] missing fields: {sorted(missing_sub_fields)}")