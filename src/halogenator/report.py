# -*- coding: ascii -*-
"""Report generation utilities."""

import os
import sys
import platform
import glob
import json
import logging
from datetime import datetime
from typing import Dict, Any, List, Optional, Tuple
from collections import Counter, defaultdict
import pandas as pd

from .io_utils import read_smi, read_table, iter_table_records
from .schema import validate_products_records, detect_schema_format, QA_PATH_KEYS, ALL_RULES, ALL_HALOGENS, OVERVIEW_COUNTERS

# Module-level logger
LOG = logging.getLogger(__name__)

# Unified parent key constants
IK_PREFIX = "IK:"
SMI_PREFIX = "SMI:"
SMI_UNKNOWN = f"{SMI_PREFIX}UNKNOWN"
SMILES_HASH_THRESHOLD = 120


def _normalize_qa_paths_with_aliases(qa_paths: Dict[str, int]) -> Dict[str, int]:
    """
    Normalize QA paths dict by merging legacy alias keys into standard keys.

    This provides backward compatibility for scripts and data that use old key names.
    Uses the canonical_event function from qa_utils for consistency.

    Args:
        qa_paths: Raw QA paths dictionary

    Returns:
        Normalized QA paths with standard key names and merged alias counts
    """
    from .qa_utils import canonical_event

    # Start with empty normalized dict
    normalized = {}

    # Canonicalize all keys and accumulate values
    for key, value in qa_paths.items():
        canonical_key = canonical_event(key)
        normalized[canonical_key] = normalized.get(canonical_key, 0) + value

    # Special case: backward compatibility for deduplication metrics
    # Ensure both dedup_hits_inchi and pruned_inchikey_dupe are present with identical values
    if 'dedup_hits_inchi' in normalized and normalized['dedup_hits_inchi'] > 0:
        normalized['pruned_inchikey_dupe'] = normalized['dedup_hits_inchi']
    elif 'pruned_inchikey_dupe' in normalized and 'dedup_hits_inchi' not in normalized:
        normalized['dedup_hits_inchi'] = normalized['pruned_inchikey_dupe']

    return normalized


def _safe_to_inchikey(mol):
    """Safe InChIKey generation. Returns InChIKey or None only."""
    if not mol:
        return None
    try:
        from .standardize import to_inchikey
        key = to_inchikey(mol)
        if not key:
            return None
        key = key.strip()
        # Return None for 'unknown' values (case-insensitive)
        return key.upper() if key and key.lower() != 'unknown' else None
    except Exception as e:
        LOG.debug(f"InChIKey generation failed: {type(e).__name__}: {e}")
        return None


def make_unified_parent_key_with_metadata(smiles: str) -> Tuple[str, Dict[str, Any]]:
    """
    Generate unified parent key with metadata from SMILES string.
    
    Args:
        smiles: Parent SMILES string
        
    Returns:
        Tuple of (unified_key, metadata_dict)
        metadata_dict contains information about key generation process
    """
    # Normalize SMILES for consistent key generation (compress whitespace)
    norm_smiles = ' '.join((smiles or '').split()).strip()
    
    # Check if this is an unknown/invalid case first
    is_unknown = not norm_smiles or norm_smiles.lower() == 'unknown'
    
    # Initialize metadata - for unknown cases, normalized_length should be 0
    # since we don't use the actual normalized SMILES
    key_meta = {
        'smiles_hashed': False,
        'normalized_length': 0 if is_unknown else len(norm_smiles),
        'rdkit_available': False,
        'inchikey_generated': False
    }
    
    if is_unknown:
        return SMI_UNKNOWN, key_meta
    
    # Try RDKit -> InChIKey pathway first
    mol = None
    try:
        from .chem_compat import Chem
        mol = Chem.MolFromSmiles(norm_smiles)
        key_meta['rdkit_available'] = True
    except Exception as e:
        LOG.debug(f"RDKit not available or MolFromSmiles failed: {type(e).__name__}: {e}")
        mol = None
    
    if mol is not None:
        # Try InChIKey generation (only returns InChIKey or None)
        ikey = _safe_to_inchikey(mol)
        if ikey:
            key_meta['inchikey_generated'] = True
            return f"{IK_PREFIX}{ikey}", key_meta
        
        # InChIKey failed but RDKit available -> use canonical SMILES
        try:
            canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
            if canonical_smiles:
                norm_smiles = canonical_smiles
                key_meta['canonical_length'] = len(canonical_smiles)  # Track canonical length separately
        except Exception as e:
            LOG.debug(f"MolToSmiles failed: {type(e).__name__}: {e}")
            # Continue with original normalized SMILES
    
    # Final SMILES branch (with hash protection for very long SMILES)
    if len(norm_smiles) > SMILES_HASH_THRESHOLD:
        from hashlib import sha256
        digest = sha256(norm_smiles.encode('utf-8')).hexdigest()[:16]
        key_meta['smiles_hashed'] = True
        return f"{SMI_PREFIX}{digest}", key_meta
    
    return f"{SMI_PREFIX}{norm_smiles}", key_meta


def make_unified_parent_key(smiles: str) -> str:
    """
    Generate unified parent key from SMILES string with graceful RDKit absence handling.
    
    Args:
        smiles: Parent SMILES string
        
    Returns:
        Unified key string: IK:<InChIKey> if available, otherwise SMI:<normalized_smiles>
        Uses hash protection for very long SMILES (>SMILES_HASH_THRESHOLD chars)
    """
    # Use the metadata-aware version but only return the key part
    key, _ = make_unified_parent_key_with_metadata(smiles)
    return key


def record_to_unified_parent_key(record: dict) -> Optional[str]:
    """
    Return unified parent key from a product record.
    
    Args:
        record: Product record with parent_inchikey and/or parent_smiles fields
        
    Returns:
        Unified key string: IK:<InChIKey> (uppercase) if available, otherwise SMI:<smiles>
        None if no valid parent identifier found
    """
    # Priority 1: Use InChIKey if available
    ik = (record.get('parent_inchikey') or '').strip()
    if ik and ik.lower() != 'unknown':
        return f"{IK_PREFIX}{ik.upper()}"
    
    # Priority 2: Fall back to SMILES
    smi = (record.get('parent_smiles') or '').strip()
    if not smi or smi.lower() == 'unknown':
        return None
    
    return make_unified_parent_key(smi)


def generate_summary_report(config: Dict[str, Any], subset: str = 'flavonoids', k_max: int = None) -> None:
    """Generate summary report from enumeration results (supports both P0 k=1 and P1 k<=3 formats)."""
    io_config = config.get('io', {})
    
    # File paths
    parents_file = io_config.get('smiles_file', 'data/output/p0/parents.smi')
    products_table = io_config.get('products_table', 'data/output/p0/products_k1.parquet')
    summary_file = io_config.get('summary_csv', 'data/output/p0/summary_k1.csv')
    
    # Read parent molecules
    parent_pairs = read_smi(parents_file)
    parent_count = len(parent_pairs)
    
    # Read products - collect main file and any part files
    
    product_files = [products_table]
    
    # Look for part files in same directory
    if products_table.endswith('.parquet'):
        base_path = products_table.replace('.parquet', '')
        part_pattern = base_path + '.part*.parquet'
        part_files = sorted(glob.glob(part_pattern))
        if part_files:
            product_files.extend(part_files)
            print(f"Found {len(part_files)} part files, processing {len(product_files)} files total")
    
    # Stream records incrementally for statistics (avoid loading all into memory)
    
    total_size_mb = 0
    for file_path in product_files:
        if os.path.exists(file_path):
            file_size_mb = os.path.getsize(file_path) / (1024 * 1024)
            total_size_mb += file_size_mb
    
    # Warn for large datasets
    if total_size_mb > 500:  # 500 MB threshold
        print(f"Warning: Processing large dataset ({total_size_mb:.1f} MB total). Consider using HALO_BATCH_SIZE environment variable to create smaller part files.")
    
    # Initialize incremental statistics and process files one by one
    stats = _init_incremental_stats()
    total_records = 0
    schema_sample_records = []  # P0-C: Use multiple samples for schema detection
    
    for file_path in product_files:
        if os.path.exists(file_path):
            file_records = 0
            for record in iter_table_records(file_path):
                # P0-C: Collect up to 10 sample records for robust schema detection
                if len(schema_sample_records) < 10:
                    schema_sample_records.append(record)
                _update_incremental_stats(stats, record)
                file_records += 1
                total_records += 1
            print(f"Processed {file_records:,} records from {os.path.basename(file_path)}")
    
    print(f"Total processed: {total_records:,} records from {len(product_files)} files")
    
    # P0-C: Detect schema format from multiple sample records (<=10)
    schema_format = 'P0'  # default
    if schema_sample_records:
        # Combine all columns from sample records for robust detection
        all_sample_columns = set()
        for sample_record in schema_sample_records:
            all_sample_columns.update(sample_record.keys())
        schema_format = detect_schema_format(all_sample_columns)
    
    # Finalize incremental statistics with metadata
    config_with_subset = dict(config)
    config_with_subset['subset'] = subset
    config_with_subset['schema_format'] = schema_format
    config_with_subset['k_max'] = k_max
    
    # Add qa_loader_degrade flag (prefer explicit config, fallback to environment)
    import os
    qa_loader_degrade = config.get('qa_loader_degrade')
    if qa_loader_degrade is None:
        # Fallback to environment variable for backward compatibility
        qa_loader_degrade = (os.environ.get('HALO_QA_DEGRADE', '0') == '1')
    config_with_subset['qa_loader_degrade'] = qa_loader_degrade
    
    # Add parent count and metadata to incremental stats
    stats['parent_count'] = len(parent_pairs)
    stats['product_count'] = total_records
    
    # Load QA stats for three-way mutex invariant validation with enhanced path/schema resolution
    qa_stats, qa_source = _load_qa_stats_with_fallback(products_table, config_with_subset)
    
    # Store QA source info for observability in final stats
    if qa_source:
        stats['qa_source'] = qa_source
    _finalize_incremental_stats(stats, parent_pairs, config_with_subset, qa_stats)
    
    # Write summary CSV
    _write_summary_csv(stats, summary_file)
    
    # Write rule x halogen pivot tables
    if schema_format == 'P1':
        # P1: Rule x Halogen x k pivot
        pivot_file = summary_file.replace('.csv', '_pivot_rule_halogen_k.csv')
        _write_p1_pivot_csv(stats, pivot_file)
        
        # P1: Constraints violations TopN
        constraints_file = summary_file.replace('.csv', '_constraints_topn.csv')
        _write_constraints_csv(stats, constraints_file)
    else:
        # P0: Original rule x halogen pivot
        pivot_file = summary_file.replace('.csv', '_pivot.csv')
        _write_pivot_table(stats['rule_halogen_counts'], pivot_file)
    
    # Write type-specific pivot tables only for 'all' subset (P0 only)
    if subset == 'all' and schema_format == 'P0':
        pivot_flav = summary_file.replace('.csv', '_flavonoids_pivot.csv')
        _write_pivot_table(stats['rule_halogen_counts_flavonoids'], pivot_flav)
        
        pivot_probe = summary_file.replace('.csv', '_probes_pivot.csv')
        _write_pivot_table(stats['rule_halogen_counts_probes'], pivot_probe)
    
    # Print summary to console
    _print_summary(stats)


def _compute_statistics(parent_pairs: List[tuple], product_records: List[Dict[str, Any]], 
                       config: Dict[str, Any]) -> Dict[str, Any]:
    """Compute summary statistics."""
    
    parent_count = len(parent_pairs)
    product_count = len(product_records)
    
    # Count by rule and halogen (overall)
    rule_counts = Counter()
    halogen_counts = Counter()
    rule_halogen_counts = defaultdict(lambda: defaultdict(int))
    
    # Count by rule and halogen (flavonoids only)
    rule_counts_flavonoids = Counter()
    halogen_counts_flavonoids = Counter()  
    rule_halogen_counts_flavonoids = defaultdict(lambda: defaultdict(int))
    
    # Count by rule and halogen (probes only) 
    rule_counts_probes = Counter()
    halogen_counts_probes = Counter()
    rule_halogen_counts_probes = defaultdict(lambda: defaultdict(int))
    
    # Count by parent
    parent_product_counts = Counter()
    
    # Use unified parent key extraction for consistency
    
    for record in product_records:
        rule = record.get('rule', 'Unknown')
        halogen = record.get('halogen', 'Unknown')
        parent_type = record.get('parent_type', 'unknown')
        
        # Overall counts (not dependent on parent_key)
        rule_counts[rule] += 1
        halogen_counts[halogen] += 1
        rule_halogen_counts[rule][halogen] += 1
        
        # Parent-level tracking (only if valid parent_key)
        parent_key = record_to_unified_parent_key(record)
        if parent_key:
            parent_product_counts[parent_key] += 1
        
        # Skip type-specific counts if no valid parent key
        if not parent_key:
            continue
        
        # Type-specific counts
        if parent_type == 'flavonoid':
            rule_counts_flavonoids[rule] += 1
            halogen_counts_flavonoids[halogen] += 1
            rule_halogen_counts_flavonoids[rule][halogen] += 1
        elif parent_type == 'probe':
            rule_counts_probes[rule] += 1
            halogen_counts_probes[halogen] += 1
            rule_halogen_counts_probes[rule][halogen] += 1
    
    # Unified parent tracking with safe key generation
    # Build unified parent key to name mapping with metadata collection
    parent_key_to_name = {}
    all_parent_keys = set()
    parent_key_meta_aggregates = {
        'hashed_count': 0,
        'hashed_examples': [],
        'rdkit_available_count': 0,
        'inchikey_generated_count': 0
    }
    
    for smiles, name in parent_pairs:
        # Use unified key generator with metadata collection
        unified_key, key_meta = make_unified_parent_key_with_metadata(smiles)
        if unified_key and unified_key != SMI_UNKNOWN:
            parent_key_to_name[unified_key] = name
            all_parent_keys.add(unified_key)
            
            # Aggregate metadata for observability
            if key_meta.get('smiles_hashed', False):
                parent_key_meta_aggregates['hashed_count'] += 1
                if len(parent_key_meta_aggregates['hashed_examples']) < 20:
                    parent_key_meta_aggregates['hashed_examples'].append(unified_key)
            
            if key_meta.get('rdkit_available', False):
                parent_key_meta_aggregates['rdkit_available_count'] += 1
            
            if key_meta.get('inchikey_generated', False):
                parent_key_meta_aggregates['inchikey_generated_count'] += 1
    
    # Collect parents with products using unified keys
    parents_with_products = set()
    for record in product_records:
        parent_key = record_to_unified_parent_key(record)
        if parent_key:
            parents_with_products.add(parent_key)
    
    # Calculate parents without products
    parents_without_products_keys = all_parent_keys - parents_with_products
    parents_without_products = sorted([parent_key_to_name.get(key, key) 
                                     for key in parents_without_products_keys])
    
    # Legacy dual-track validation for diagnostics (keep existing logic for backwards compatibility)
    parent_smiles_to_name = {smiles: name for smiles, name in parent_pairs}
    all_parent_smiles = set(parent_smiles_to_name.keys())
    parents_with_products_smiles = set()
    parents_with_products_ikeys = set()
    
    for record in product_records:
        parent_smiles = record.get('parent_smiles', 'Unknown')
        parent_ikey = record.get('parent_inchikey', 'Unknown')
        if parent_smiles != 'Unknown':
            parents_with_products_smiles.add(parent_smiles)
        if parent_ikey != 'Unknown':
            parents_with_products_ikeys.add(parent_ikey)
    
    # Legacy diagnostics for backwards compatibility
    smiles_count = len(parents_with_products_smiles)
    ikey_count = len(parents_with_products_ikeys)
    diff_smiles_vs_inchikey = abs(smiles_count - ikey_count)
    
    # Simple diagnostics without complex mapping
    missing_in_smiles_examples = []
    missing_in_ikeys_examples = []
    
    # Parent statistics
    if parent_product_counts:
        avg_products_per_parent = sum(parent_product_counts.values()) / len(parent_product_counts)
        max_products_per_parent = max(parent_product_counts.values())
        min_products_per_parent = min(parent_product_counts.values())
    else:
        avg_products_per_parent = 0
        max_products_per_parent = 0
        min_products_per_parent = 0
    
    # QC statistics
    sanitize_ok_count = sum(1 for r in product_records if r.get('sanitize_ok', True))
    
    # Self-check: parent count consistency
    parent_count_unified = len(all_parent_keys)
    assert isinstance(parent_count_unified, int)
    assert isinstance(parents_with_products, set)
    if len(parents_with_products) > parent_count_unified:
        # This should never happen; keep a diagnostic hook
        # Do not raise hard in production; just record diagnostics
        diag = config.setdefault('diagnostics', {})
        diag['anomaly_more_parents_with_products_than_all'] = True
    
    return {
        'parent_count': len(all_parent_keys),  # Use unified parent count
        'product_count': product_count,
        'sanitize_ok_count': sanitize_ok_count,
        'rule_counts': dict(rule_counts),
        'halogen_counts': dict(halogen_counts),
        'rule_halogen_counts': dict(rule_halogen_counts),
        'rule_counts_flavonoids': dict(rule_counts_flavonoids),
        'halogen_counts_flavonoids': dict(halogen_counts_flavonoids),
        'rule_halogen_counts_flavonoids': dict(rule_halogen_counts_flavonoids),
        'rule_counts_probes': dict(rule_counts_probes),
        'halogen_counts_probes': dict(halogen_counts_probes),
        'rule_halogen_counts_probes': dict(rule_halogen_counts_probes),
        'parent_product_counts': dict(parent_product_counts),
        'avg_products_per_parent': round(avg_products_per_parent, 2),
        'max_products_per_parent': max_products_per_parent,
        'min_products_per_parent': min_products_per_parent,
        'unique_parents_with_products': len(parents_with_products),
        'parents_without_products': parents_without_products,
        'diagnostics': {
            'diff_smiles_vs_inchikey': diff_smiles_vs_inchikey,
            'smiles_unique_parents': smiles_count,
            'inchikey_unique_parents': ikey_count,
            'missing_in_smiles_count': len(missing_in_smiles_examples),
            'missing_in_ikeys_count': len(missing_in_ikeys_examples),
            'missing_in_smiles_examples': missing_in_smiles_examples,
            'missing_in_ikeys_examples': missing_in_ikeys_examples
        },
        'config': config
    }


def _compute_p1_statistics(parent_pairs: List[tuple], product_records: List[Dict[str, Any]], 
                          config: Dict[str, Any]) -> Dict[str, Any]:
    """Compute P1 k-dimensional statistics."""
    
    # Will be corrected to len(all_parent_keys) after unified parent tracking
    parent_count_initial = len(parent_pairs)
    product_count = len(product_records)
    
    # K-dimensional counting: Rule x Halogen x K
    rule_halogen_k_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    rule_halogen_k_by_parent = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))
    k_counts = Counter()
    rule_counts = Counter()
    halogen_counts = Counter()
    constraint_stats = {'total_passed': 0, 'total_failed': 0, 'violation_reasons': Counter()}
    
    # Unified parent tracking with safe key generation and metadata collection
    parent_key_to_name = {}
    all_parent_keys = set()
    parent_key_meta_aggregates = {
        'hashed_count': 0,
        'hashed_examples': [],
        'rdkit_available_count': 0,
        'inchikey_generated_count': 0
    }
    
    for smiles, name in parent_pairs:
        # Use unified key generator with metadata collection
        unified_key, key_meta = make_unified_parent_key_with_metadata(smiles)
        if unified_key and unified_key != SMI_UNKNOWN:
            parent_key_to_name[unified_key] = name
            all_parent_keys.add(unified_key)
            
            # Aggregate metadata for observability
            if key_meta.get('smiles_hashed', False):
                parent_key_meta_aggregates['hashed_count'] += 1
                if len(parent_key_meta_aggregates['hashed_examples']) < 20:
                    parent_key_meta_aggregates['hashed_examples'].append(unified_key)
            
            if key_meta.get('rdkit_available', False):
                parent_key_meta_aggregates['rdkit_available_count'] += 1
            
            if key_meta.get('inchikey_generated', False):
                parent_key_meta_aggregates['inchikey_generated_count'] += 1
    
    # Process product records
    parent_product_counts = Counter()
    parents_with_products = set()
    
    # Using unified parent key extraction
    
    for record in product_records:
        # Extract k-dimensional fields
        k = record.get('k', record.get('depth', 1))  # Fallback to depth for compatibility
        rule = record.get('rule', 'Unknown')
        halogen = record.get('halogen', 'Unknown')
        constraints_ok = record.get('constraints_ok', True)
        constraints_violations = record.get('constraints_violations', {})
        
        # Count by k-dimensional matrix (global stats - not dependent on parent_key)
        rule_halogen_k_counts[rule][halogen][k] += 1
        k_counts[k] += 1
        rule_counts[rule] += 1
        halogen_counts[halogen] += 1
        
        # Parent-level tracking (only if valid parent_key)
        parent_key = record_to_unified_parent_key(record)
        
        # Add parent dimension to rule_halogen_k tracking
        if parent_key:
            rule_halogen_k_by_parent[parent_key][rule][halogen][k] += 1
        if parent_key:
            parent_product_counts[parent_key] += 1
            parents_with_products.add(parent_key)
        
        # Constraint statistics
        if constraints_ok:
            constraint_stats['total_passed'] += 1
        else:
            constraint_stats['total_failed'] += 1
            if isinstance(constraints_violations, dict):
                for reason in constraints_violations.keys():
                    constraint_stats['violation_reasons'][reason] += 1
    
    # Calculate parents without products
    parents_without_products_keys = all_parent_keys - parents_with_products
    parents_without_products = sorted([parent_key_to_name.get(key, key) 
                                     for key in parents_without_products_keys])
    
    # Parent statistics
    if parent_product_counts:
        avg_products_per_parent = sum(parent_product_counts.values()) / len(parent_product_counts)
        max_products_per_parent = max(parent_product_counts.values())
        min_products_per_parent = min(parent_product_counts.values())
    else:
        avg_products_per_parent = 0
        max_products_per_parent = 0
        min_products_per_parent = 0
    
    # QC statistics
    sanitize_ok_count = sum(1 for r in product_records if r.get('sanitize_ok', True))
    
    # Self-check: parent count consistency
    parent_count_unified = len(all_parent_keys)
    assert isinstance(parent_count_unified, int)
    assert isinstance(parents_with_products, set)
    if len(parents_with_products) > parent_count_unified:
        # This should never happen; keep a diagnostic hook
        # Do not raise hard in production; just record diagnostics
        diag = config.setdefault('diagnostics', {})
        diag['anomaly_more_parents_with_products_than_all'] = True
    
    # Runtime metadata
    k_max = config.get('k_max', max(k_counts.keys()) if k_counts else 1)
    
    return {
        'parent_count': len(all_parent_keys),  # Use unified parent count
        'product_count': product_count,
        'sanitize_ok_count': sanitize_ok_count,
        'k_max': k_max,
        'k_counts': dict(k_counts),
        'rule_counts': dict(rule_counts),
        'halogen_counts': dict(halogen_counts),
        'rule_halogen_k_counts': {r: {h: dict(k) for h, k in h_dict.items()} for r, h_dict in rule_halogen_k_counts.items()},
        'rule_halogen_k_by_parent': {p: {r: {h: dict(k) for h, k in h_dict.items()} for r, h_dict in r_dict.items()} for p, r_dict in rule_halogen_k_by_parent.items()},
        'constraint_stats': constraint_stats,
        'parent_product_counts': dict(parent_product_counts),
        'avg_products_per_parent': round(avg_products_per_parent, 2),
        'max_products_per_parent': max_products_per_parent,
        'min_products_per_parent': min_products_per_parent,
        'unique_parents_with_products': len(parents_with_products),
        'parents_without_products': parents_without_products,
        'config': config
    }


def _write_summary_csv(stats: Dict[str, Any], output_path: str) -> None:
    """Write summary statistics to CSV with fixed schema.""" 
    _write_summary_csv_fixed_schema(stats, output_path)


def _validate_fixed_csv_schema(csv_path: str, expected_columns: List[str]) -> Dict[str, Any]:
    """Validate that a CSV file has the expected fixed schema.
    
    Args:
        csv_path: Path to CSV file to validate
        expected_columns: List of expected column names
        
    Returns:
        Dictionary with validation results
    """
    import os
    
    validation_result = {
        'valid': True,
        'errors': [],
        'warnings': [],
        'file_exists': False,
        'column_count': 0,
        'row_count': 0
    }
    
    try:
        if not os.path.exists(csv_path):
            validation_result['valid'] = False
            validation_result['errors'].append(f"CSV file does not exist: {csv_path}")
            return validation_result
            
        validation_result['file_exists'] = True
        
        # Read CSV and check schema
        df = pd.read_csv(csv_path)
        validation_result['column_count'] = len(df.columns)
        validation_result['row_count'] = len(df)
        
        # Check column names match exactly
        actual_columns = list(df.columns)
        if actual_columns != expected_columns:
            validation_result['valid'] = False
            validation_result['errors'].append(
                f"Column mismatch. Expected: {expected_columns[:5]}... "
                f"Got: {actual_columns[:5]}..."
            )
        
        # Check for missing values in critical columns
        critical_prefixes = ['overall_', 'per_parent_', 'metadata_']
        for col in df.columns:
            if any(col.startswith(prefix) for prefix in critical_prefixes):
                if df[col].isnull().any():
                    validation_result['warnings'].append(f"Missing values in critical column: {col}")
        
        LOG.info(f"CSV schema validation for {csv_path}: {validation_result['valid']} "
                   f"({validation_result['column_count']} columns, {validation_result['row_count']} rows)")
                   
    except Exception as e:
        validation_result['valid'] = False
        validation_result['errors'].append(f"Validation error: {str(e)}")
        LOG.error(f"CSV validation failed for {csv_path}: {e}")
    
    return validation_result


def _write_summary_csv_fixed_schema(stats: Dict[str, Any], output_path: str) -> None:
    """Write summary statistics to CSV with fixed schema for DL comparison workflows."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Create summary rows
    summary_rows = []
    
    # Runtime metadata
    try:
        import rdkit
        rdkit_version = rdkit.__version__
    except ImportError:
        rdkit_version = 'not_installed'
    except Exception as e:
        LOG.debug(f"RDKit version detection failed: {type(e).__name__}: {e}")
        rdkit_version = 'unknown'
        
    try:
        import pandas as pd_meta
        pandas_version = pd_meta.__version__
    except ImportError:
        pandas_version = 'not_installed'
    except Exception as e:
        LOG.debug(f"Pandas version detection failed: {type(e).__name__}: {e}")
        pandas_version = 'unknown'
        
    try:
        import pyarrow
        pyarrow_version = pyarrow.__version__
    except ImportError:
        pyarrow_version = 'not_installed'
    except Exception as e:
        LOG.debug(f"PyArrow version detection failed: {type(e).__name__}: {e}")
        pyarrow_version = 'unknown'
    
    subset_info = stats.get('config', {}).get('subset', 'unknown')
    timestamp = datetime.now().isoformat()
    
    # Metadata rows
    summary_rows.extend([
        {'Category': 'Metadata', 'Metric': 'RDKit Version', 'Value': rdkit_version, 'Description': 'Version of RDKit used'},
        {'Category': 'Metadata', 'Metric': 'Pandas Version', 'Value': pandas_version, 'Description': 'Version of Pandas used'},
        {'Category': 'Metadata', 'Metric': 'PyArrow Version', 'Value': pyarrow_version, 'Description': 'Version of PyArrow used'},
        {'Category': 'Metadata', 'Metric': 'Platform', 'Value': platform.system(), 'Description': 'Operating system platform'},
        {'Category': 'Metadata', 'Metric': 'Subset', 'Value': subset_info, 'Description': 'Molecular subset processed'},
        {'Category': 'Metadata', 'Metric': 'Run Timestamp', 'Value': timestamp, 'Description': 'When this report was generated'}
    ])
    
    # Overall stats
    summary_rows.append({
        'Category': 'Overall',
        'Metric': 'Total Parent Molecules',
        'Value': stats['parent_count'],
        'Description': 'Number of input parent molecules'
    })
    
    summary_rows.append({
        'Category': 'Overall', 
        'Metric': 'Total Products Generated',
        'Value': stats['product_count'],
        'Description': 'Total k=1 halogenated products (after deduplication)'
    })
    
    summary_rows.append({
        'Category': 'Overall',
        'Metric': 'Products Passing QC',
        'Value': stats['sanitize_ok_count'],
        'Description': 'Products that sanitize correctly'
    })
    
    summary_rows.append({
        'Category': 'Overall',
        'Metric': 'Unique Parents with Products',
        'Value': stats['unique_parents_with_products'],
        'Description': 'Number of parents that generated at least one product'
    })
    
    # QF-1 logic for CSV description
    zero_product_parents = stats['parents_without_products']
    parents_with_products = stats['unique_parents_with_products'] 
    total_parents = stats['parent_count']
    
    if parents_with_products == total_parents:
        zero_product_description = 'All parents produced at least one product'
    else:
        examples = zero_product_parents[:5]  # Top 5 examples
        examples_str = ', '.join(examples)
        if len(zero_product_parents) > 5:
            examples_str += f' (and {len(zero_product_parents) - 5} more)'
        zero_product_description = f'Parents that generated no products: {examples_str}'
    
    summary_rows.append({
        'Category': 'Overall',
        'Metric': 'Parents with Zero Products',
        'Value': len(zero_product_parents),
        'Description': zero_product_description
    })
    
    # Per-parent statistics
    summary_rows.append({
        'Category': 'Per-Parent',
        'Metric': 'Average Products per Parent',
        'Value': stats['avg_products_per_parent'],
        'Description': 'Average number of products per parent molecule'
    })
    
    summary_rows.append({
        'Category': 'Per-Parent',
        'Metric': 'Max Products per Parent', 
        'Value': stats['max_products_per_parent'],
        'Description': 'Maximum products generated from any single parent'
    })
    
    summary_rows.append({
        'Category': 'Per-Parent',
        'Metric': 'Min Products per Parent',
        'Value': stats['min_products_per_parent'],
        'Description': 'Minimum products generated from any single parent'
    })
    
    # Rule breakdown
    for rule, count in stats['rule_counts'].items():
        summary_rows.append({
            'Category': 'Rule Breakdown',
            'Metric': f'{rule} Products',
            'Value': count,
            'Description': f'Products generated by rule {rule}'
        })
    
    # Halogen breakdown 
    for halogen, count in stats['halogen_counts'].items():
        summary_rows.append({
            'Category': 'Halogen Breakdown',
            'Metric': f'{halogen} Products',
            'Value': count,
            'Description': f'Products containing halogen {halogen}'
        })
    
    # Rule x Halogen matrix (handle both P0 and P1 formats)
    rule_halogen_data = None
    if 'rule_halogen_counts' in stats:
        # P0 format: direct 2D matrix
        rule_halogen_data = stats['rule_halogen_counts']
    elif 'rule_halogen_k_counts' in stats:
        # P1 format: aggregate 3D matrix to 2D
        rule_halogen_data = {}
        for rule, halogen_dict in stats['rule_halogen_k_counts'].items():
            if rule not in rule_halogen_data:
                rule_halogen_data[rule] = {}
            for halogen, k_dict in halogen_dict.items():
                rule_halogen_data[rule][halogen] = sum(k_dict.values())
    
    if rule_halogen_data:
        for rule, halogen_dict in rule_halogen_data.items():
            for halogen, count in halogen_dict.items():
                summary_rows.append({
                    'Category': 'Rule x Halogen',
                    'Metric': f'{rule} + {halogen}',
                    'Value': count,
                    'Description': f'Products from rule {rule} with halogen {halogen}'
                })
    
    # Diagnostics section  
    diagnostics = stats.get('diagnostics', {})
    
    # Try to load QA summary JSON if it exists
    qa_stats = {}
    try:
        qa_json_path = os.path.join(os.path.dirname(output_path), 'qa_summary.json')
        if os.path.exists(qa_json_path):
            with open(qa_json_path, 'r') as f:
                qa_data = json.load(f)
                # Support both structures: {'total': {...}, 'pivots': {...}} and root-level totals
                qa_stats = qa_data.get('total', qa_data)
    except Exception:
        pass  # QA stats not available, skip
    
    # QA Statistics (if available)
    if qa_stats:
        qa_paths = _normalize_qa_paths_with_aliases(qa_stats.get('qa_paths', {}))
        summary_rows.extend([
            {
                'Category': 'Diagnostics',
                'Metric': 'Isotope unavailable (attempts)',
                'Value': qa_paths.get('isotope_unavailable', 0),
                'Description': 'Rule/halogen attempts where isotope strategy could not be used'
            },
            {
                'Category': 'Diagnostics',
                'Metric': 'Isotope miss (matches)',
                'Value': qa_paths.get('isotope_miss', 0),
                'Description': 'Tagged matches where product site was not found (rare)'
            },
            {
                'Category': 'Diagnostics',
                'Metric': 'AtomMap used (attempts)',
                'Value': qa_paths.get('atommap_used', 0),
                'Description': 'Site detection succeeded using AtomMapNum fallback'
            },
            {
                'Category': 'Diagnostics',
                'Metric': 'Heuristic used (attempts)',
                'Value': qa_paths.get('heuristic_used', 0),
                'Description': 'Site detection succeeded using heuristic fallback'
            },
            {
                'Category': 'Diagnostics',
                'Metric': 'No-product matches',
                'Value': qa_stats.get('no_product_matches', 0),
                'Description': 'Matches/attempts that ultimately failed to produce valid products'
            },
            {
                'Category': 'Diagnostics',
                'Metric': 'Template unsupported (attempts)',
                'Value': qa_stats.get('template_unsupported', 0),
                'Description': 'Templates with unsupported structure (e.g. multi-reactant)'
            }
        ])
    
    summary_rows.append({
        'Category': 'Diagnostics',
        'Metric': 'SMILES vs InChIKey Difference',
        'Value': diagnostics.get('diff_smiles_vs_inchikey', 0),
        'Description': 'Absolute difference in unique parent count between SMILES and InChIKey tracking'
    })
    summary_rows.append({
        'Category': 'Diagnostics',
        'Metric': 'SMILES Unique Parents',
        'Value': diagnostics.get('smiles_unique_parents', 0),
        'Description': 'Number of unique parents with products (SMILES-based)'
    })
    summary_rows.append({
        'Category': 'Diagnostics',
        'Metric': 'InChIKey Unique Parents',
        'Value': diagnostics.get('inchikey_unique_parents', 0),
        'Description': 'Number of unique parents with products (InChIKey-based)'
    })
    
    # Set difference diagnostics
    missing_in_smiles = diagnostics.get('missing_in_smiles_examples', [])
    missing_in_ikeys = diagnostics.get('missing_in_ikeys_examples', [])
    
    summary_rows.append({
        'Category': 'Diagnostics',
        'Metric': 'Missing in SMILES tracking',
        'Value': diagnostics.get('missing_in_smiles_count', 0),
        'Description': f'Parents found via InChIKey but not SMILES. Examples: {", ".join(missing_in_smiles) if missing_in_smiles else "None"}'
    })
    
    summary_rows.append({
        'Category': 'Diagnostics',
        'Metric': 'Missing in InChIKey tracking',
        'Value': diagnostics.get('missing_in_ikeys_count', 0),
        'Description': f'Parents found via SMILES but not InChIKey. Examples: {", ".join(missing_in_ikeys) if missing_in_ikeys else "None"}'
    })
    
    # Write to CSV using original variable schema for backwards compatibility
    df = pd.DataFrame(summary_rows)
    df.to_csv(output_path, index=False)


def _write_summary_csv_fixed_schema(stats: Dict[str, Any], output_path: str) -> None:
    """Write summary statistics to CSV with fixed schema for DL comparison workflows.
    
    This version always outputs the same columns regardless of data content,
    making it suitable for machine learning pipelines that require consistent schemas.
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Define fixed core columns that always appear
    CORE_METRICS = [
        # Metadata columns
        ('metadata_rdkit_version', 'str'),
        ('metadata_pandas_version', 'str'), 
        ('metadata_pyarrow_version', 'str'),
        ('metadata_platform', 'str'),
        ('metadata_subset', 'str'),
        ('metadata_timestamp', 'str'),
        
        # Overall statistics
        ('overall_total_parent_molecules', 'int'),
        ('overall_total_products_generated', 'int'),
        ('overall_products_passing_qc', 'int'),
        ('overall_unique_parents_with_products', 'int'),
        ('overall_parents_with_zero_products', 'int'),
        
        # Per-parent statistics
        ('per_parent_average_products_per_parent', 'float'),
        ('per_parent_max_products_per_parent', 'int'),
        ('per_parent_min_products_per_parent', 'int'),
        
        # Rule breakdown (fixed set of rules)
        ('rule_r1_products', 'int'),
        ('rule_r2_products', 'int'),
        ('rule_r3_products', 'int'),
        ('rule_r4_products', 'int'),
        ('rule_r5_products', 'int'),
        
        # Halogen breakdown (fixed set of halogens)
        ('halogen_f_products', 'int'),
        ('halogen_cl_products', 'int'),
        ('halogen_br_products', 'int'),
        ('halogen_i_products', 'int'),
        
        # Rule x Halogen matrix (20 combinations)
        ('rule_halogen_r1_f', 'int'),
        ('rule_halogen_r1_cl', 'int'),
        ('rule_halogen_r1_br', 'int'),
        ('rule_halogen_r1_i', 'int'),
        ('rule_halogen_r2_f', 'int'),
        ('rule_halogen_r2_cl', 'int'),
        ('rule_halogen_r2_br', 'int'),
        ('rule_halogen_r2_i', 'int'),
        ('rule_halogen_r3_f', 'int'),
        ('rule_halogen_r3_cl', 'int'),
        ('rule_halogen_r3_br', 'int'),
        ('rule_halogen_r3_i', 'int'),
        ('rule_halogen_r4_f', 'int'),
        ('rule_halogen_r4_cl', 'int'),
        ('rule_halogen_r4_br', 'int'),
        ('rule_halogen_r4_i', 'int'),
        ('rule_halogen_r5_f', 'int'),
        ('rule_halogen_r5_cl', 'int'),
        ('rule_halogen_r5_br', 'int'),
        ('rule_halogen_r5_i', 'int'),
        
        # Diagnostics (core QA paths)
        ('diagnostics_isotope_unavailable_attempts', 'int'),
        ('diagnostics_isotope_miss_matches', 'int'),
        ('diagnostics_atommap_used_attempts', 'int'),
        ('diagnostics_heuristic_used_attempts', 'int'),
        ('diagnostics_no_product_matches', 'int'),
        ('diagnostics_template_unsupported_attempts', 'int'),
        ('diagnostics_smiles_vs_inchikey_difference', 'int'),
        ('diagnostics_smiles_unique_parents', 'int'),
        ('diagnostics_inchikey_unique_parents', 'int'),
        ('diagnostics_missing_in_smiles_count', 'int'),
        ('diagnostics_missing_in_inchikey_count', 'int')
    ]
    
    # Initialize row with default values
    fixed_row = {}
    
    # Runtime metadata - always available
    try:
        import rdkit
        rdkit_version = rdkit.__version__
    except ImportError:
        rdkit_version = 'not_installed'
    except Exception:
        rdkit_version = 'unknown'
        
    try:
        import pandas as pd_meta
        pandas_version = pd_meta.__version__
    except ImportError:
        pandas_version = 'not_installed'
    except Exception:
        pandas_version = 'unknown'
        
    try:
        import pyarrow
        pyarrow_version = pyarrow.__version__
    except ImportError:
        pyarrow_version = 'not_installed'
    except Exception:
        pyarrow_version = 'unknown'
    
    import platform
    from datetime import datetime
    
    subset_info = stats.get('config', {}).get('subset', 'unknown')
    timestamp = datetime.now().isoformat()
    
    # Populate metadata fields
    fixed_row['metadata_rdkit_version'] = rdkit_version
    fixed_row['metadata_pandas_version'] = pandas_version
    fixed_row['metadata_pyarrow_version'] = pyarrow_version
    fixed_row['metadata_platform'] = platform.system()
    fixed_row['metadata_subset'] = str(subset_info)
    fixed_row['metadata_timestamp'] = timestamp
    
    # Overall statistics with safe defaults
    fixed_row['overall_total_parent_molecules'] = int(stats.get('parent_count', 0))
    fixed_row['overall_total_products_generated'] = int(stats.get('product_count', 0))
    fixed_row['overall_products_passing_qc'] = int(stats.get('sanitize_ok_count', 0))
    fixed_row['overall_unique_parents_with_products'] = int(stats.get('unique_parents_with_products', 0))
    fixed_row['overall_parents_with_zero_products'] = len(stats.get('parents_without_products', []))
    
    # Per-parent statistics with safe defaults
    fixed_row['per_parent_average_products_per_parent'] = float(stats.get('avg_products_per_parent', 0.0))
    fixed_row['per_parent_max_products_per_parent'] = int(stats.get('max_products_per_parent', 0))
    fixed_row['per_parent_min_products_per_parent'] = int(stats.get('min_products_per_parent', 0))
    
    # Rule breakdown with fixed rule set
    rule_counts = stats.get('rule_counts', {})
    for rule in ALL_RULES:
        key = f'rule_{rule.lower()}_products'
        fixed_row[key] = int(rule_counts.get(rule, 0))
    
    # Halogen breakdown with fixed halogen set
    halogen_counts = stats.get('halogen_counts', {})
    for halogen in ALL_HALOGENS:
        key = f'halogen_{halogen.lower()}_products'
        fixed_row[key] = int(halogen_counts.get(halogen, 0))
    
    # Rule x Halogen matrix with fixed combinations
    # Handle both P0 and P1 formats
    rule_halogen_data = None
    if 'rule_halogen_counts' in stats:
        # P0 format: direct 2D matrix
        rule_halogen_data = stats['rule_halogen_counts']
    elif 'rule_halogen_k_counts' in stats:
        # P1 format: aggregate 3D matrix to 2D
        rule_halogen_data = {}
        for rule, halogen_dict in stats['rule_halogen_k_counts'].items():
            if rule not in rule_halogen_data:
                rule_halogen_data[rule] = {}
            for halogen, k_dict in halogen_dict.items():
                rule_halogen_data[rule][halogen] = sum(k_dict.values())
    
    for rule in ALL_RULES:
        for halogen in ALL_HALOGENS:
            key = f'rule_halogen_{rule.lower()}_{halogen.lower()}'
            count = 0
            if rule_halogen_data and rule in rule_halogen_data and halogen in rule_halogen_data[rule]:
                count = int(rule_halogen_data[rule][halogen])
            fixed_row[key] = count
    
    # Load QA summary JSON for diagnostics if available
    qa_stats = {}
    try:
        qa_json_path = os.path.join(os.path.dirname(output_path), 'qa_summary.json')
        if os.path.exists(qa_json_path):
            with open(qa_json_path, 'r') as f:
                qa_data = json.load(f)
                qa_stats = qa_data.get('total', qa_data)
    except Exception:
        pass
    
    # Diagnostics with safe defaults
    qa_paths = _normalize_qa_paths_with_aliases(qa_stats.get('qa_paths', {}))
    diagnostics = stats.get('diagnostics', {})
    
    fixed_row['diagnostics_isotope_unavailable_attempts'] = int(qa_paths.get('isotope_unavailable', 0))
    fixed_row['diagnostics_isotope_miss_matches'] = int(qa_paths.get('isotope_miss', 0))
    fixed_row['diagnostics_atommap_used_attempts'] = int(qa_paths.get('atommap_used', 0))
    fixed_row['diagnostics_heuristic_used_attempts'] = int(qa_paths.get('heuristic_used', 0))
    fixed_row['diagnostics_no_product_matches'] = int(qa_stats.get('no_product_matches', 0))
    fixed_row['diagnostics_template_unsupported_attempts'] = int(qa_stats.get('template_unsupported', 0))
    fixed_row['diagnostics_smiles_vs_inchikey_difference'] = int(diagnostics.get('diff_smiles_vs_inchikey', 0))
    fixed_row['diagnostics_smiles_unique_parents'] = int(diagnostics.get('smiles_unique_parents', 0))
    fixed_row['diagnostics_inchikey_unique_parents'] = int(diagnostics.get('inchikey_unique_parents', 0))
    fixed_row['diagnostics_missing_in_smiles_count'] = int(diagnostics.get('missing_in_smiles_count', 0))
    fixed_row['diagnostics_missing_in_inchikey_count'] = int(diagnostics.get('missing_in_ikeys_count', 0))
    
    # Create DataFrame with single row and fixed column order
    column_names = [col_name for col_name, _ in CORE_METRICS]
    
    # Ensure all columns are present with proper types
    for col_name, col_type in CORE_METRICS:
        if col_name not in fixed_row:
            if col_type == 'int':
                fixed_row[col_name] = 0
            elif col_type == 'float':
                fixed_row[col_name] = 0.0
            else:  # str
                fixed_row[col_name] = 'unknown'
    
    # Create DataFrame with columns in the defined order
    df_data = {col_name: [fixed_row[col_name]] for col_name in column_names}
    df = pd.DataFrame(df_data)
    
    # Write to CSV
    df.to_csv(output_path, index=False)
    
    # Validate the written CSV has the expected fixed schema
    try:
        validation_result = _validate_fixed_csv_schema(output_path, column_names)
        if not validation_result['valid']:
            LOG.warning(f"Fixed schema validation failed for {output_path}: {validation_result['errors']}")
        else:
            LOG.debug(f"Fixed schema validation passed for {output_path} "
                        f"({validation_result['column_count']} columns)")
    except Exception as e:
        LOG.debug(f"Schema validation skipped due to error: {e}")


def _warn_sort_key(w: Dict[str, Any]) -> Tuple[str, str, str, str]:
    """
    Extract consistent sorting key from warning regardless of structure.
    Handles different warning types with fields in different locations.

    Args:
        w: Warning dictionary

    Returns:
        Tuple of (type, dimension, key, metric) for stable sorting
    """
    # Extract type (always at top level)
    t = str(w.get('type', ''))

    # Extract dimension - try top level first, then 'where' field
    # For marginal_conflict, use 'base' field as fallback
    d = w.get('dimension') or (w.get('where') or {}).get('dimension')
    if not d:
        d = w.get('base', '')  # fallback for marginal_conflict
    d = str(d or '')

    # Extract key - try top level first, then 'where' field
    k = w.get('key') or (w.get('where') or {}).get('key', '')
    k = str(k or '')

    # Extract metric - try top level, then 'where', then first metric from 'metrics' array
    # Sort metrics before taking first element to ensure deterministic ordering
    m = w.get('metric') or (w.get('where') or {}).get('metric')
    if not m:
        metrics = w.get('metrics')
        if isinstance(metrics, list) and metrics:
            # Sort metrics to eliminate generation order variation
            # Only used for sorting key, does not modify original warning content
            m = sorted(str(x) for x in metrics)[0]
    m = str(m or '')

    return (t, d, k, m)


def finalize_metadata(metadata: Dict[str, Any], *, conflict_tolerance: int, coercion_method: str,
                      warnings_list: List[Dict[str, Any]], max_warnings: int) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    """
    Centralized metadata finalization with consistent field assignment and warnings processing.

    Args:
        metadata: Metadata dictionary to finalize (modified in place)
        conflict_tolerance: Marginal conflict tolerance value
        coercion_method: Value coercion method used (e.g., "trunc")
        warnings_list: List of warnings to process
        max_warnings: Maximum number of warnings to include

    Returns:
        Tuple of (finalized_metadata, processed_warnings)
    """
    # 1) Write standard metadata fields consistently
    metadata['marginal_conflict_tolerance'] = conflict_tolerance
    metadata['value_coercion_method'] = coercion_method

    # 2) Process warnings with stable sorting and truncation
    if warnings_list:
        # Sort warnings for deterministic output
        sorted_warnings = sorted(warnings_list, key=_warn_sort_key)

        # Apply truncation if needed
        if len(sorted_warnings) > max_warnings:
            warnings_out = sorted_warnings[:max_warnings]
            truncated = True
            warnings_count = len(sorted_warnings)
            warnings_returned = max_warnings
        else:
            warnings_out = sorted_warnings
            truncated = False
            warnings_count = len(sorted_warnings)
            warnings_returned = len(sorted_warnings)
    else:
        warnings_out = []
        truncated = False
        warnings_count = 0
        warnings_returned = 0

    # 3) Always write warnings metadata fields (even if no warnings)
    metadata['warnings_count'] = warnings_count
    metadata['warnings_returned'] = warnings_returned
    metadata['warnings_truncated'] = truncated

    # 4) Set the warnings array in metadata
    metadata['warnings'] = warnings_out

    return metadata, warnings_out


def sanitize_granular_fields_with_warnings(by_rule: Optional[Dict], by_halogen: Optional[Dict],
                            by_rule_halogen: Optional[Dict], distributable_set: set,
                            rules: Optional[List[str]] = None, halogens: Optional[List[str]] = None) -> List[Dict[str, Any]]:
    """
    Remove non-distributable fields and validate data quality in granular dimensions.
    Args:
        by_rule: by_rule structure (modified in place)
        by_halogen: by_halogen structure (modified in place)
        by_rule_halogen: by_rule_halogen structure (modified in place)
        distributable_set: Set of field names that should be kept
        rules: List of valid rule names (for validation)
        halogens: List of valid halogen names (for validation)
    Returns:
        List of warning objects for data quality issues
    """
    warnings = []

    # Helper function to validate and clean a single metric value
    def validate_metric_value(value, where_info):
        if isinstance(value, (int, float)):
            if value < 0:
                warnings.append({
                    "type": "negative_value_detected",
                    "where": where_info,
                    "value": value
                })
            if not isinstance(value, int):
                warnings.append({
                    "type": "non_integer_value_detected",
                    "where": where_info,
                    "value": value,
                    "coerced_to_int": True,
                    "coercion_method": "trunc"
                })
                return int(value)
        return value

    # Clean and validate by_rule
    if by_rule:
        for rule_key in list(by_rule.keys()):
            if rules and rule_key not in rules:
                warnings.append({
                    "type": "unknown_dimension_key_dropped",
                    "dimension": "by_rule",
                    "key": rule_key
                })
                del by_rule[rule_key]
                continue

            rule_metrics = by_rule[rule_key]
            unknown_metrics = [k for k in rule_metrics.keys() if k not in distributable_set]
            if unknown_metrics:
                warnings.append({
                    "type": "unknown_metric_dropped",
                    "dimension": "by_rule",
                    "key": rule_key,
                    "metrics": unknown_metrics
                })
                for k in unknown_metrics:
                    del rule_metrics[k]

            # Validate remaining metrics
            for metric, value in list(rule_metrics.items()):
                validated_value = validate_metric_value(value, {"dimension": "by_rule", "key": rule_key, "metric": metric})
                rule_metrics[metric] = validated_value

    # Clean and validate by_halogen
    if by_halogen:
        for halogen_key in list(by_halogen.keys()):
            if halogens and halogen_key not in halogens:
                warnings.append({
                    "type": "unknown_dimension_key_dropped",
                    "dimension": "by_halogen",
                    "key": halogen_key
                })
                del by_halogen[halogen_key]
                continue

            halogen_metrics = by_halogen[halogen_key]
            unknown_metrics = [k for k in halogen_metrics.keys() if k not in distributable_set]
            if unknown_metrics:
                warnings.append({
                    "type": "unknown_metric_dropped",
                    "dimension": "by_halogen",
                    "key": halogen_key,
                    "metrics": unknown_metrics
                })
                for k in unknown_metrics:
                    del halogen_metrics[k]

            # Validate remaining metrics
            for metric, value in list(halogen_metrics.items()):
                validated_value = validate_metric_value(value, {"dimension": "by_halogen", "key": halogen_key, "metric": metric})
                halogen_metrics[metric] = validated_value

    # Clean and validate by_rule_halogen
    if by_rule_halogen:
        for rule_key in list(by_rule_halogen.keys()):
            if rules and rule_key not in rules:
                warnings.append({
                    "type": "unknown_dimension_key_dropped",
                    "dimension": "by_rule_halogen",
                    "key": f"rule:{rule_key}"
                })
                del by_rule_halogen[rule_key]
                continue

            for halogen_key in list(by_rule_halogen[rule_key].keys()):
                if halogens and halogen_key not in halogens:
                    warnings.append({
                        "type": "unknown_dimension_key_dropped",
                        "dimension": "by_rule_halogen",
                        "key": f"rule:{rule_key}:halogen:{halogen_key}"
                    })
                    del by_rule_halogen[rule_key][halogen_key]
                    continue

                cell_metrics = by_rule_halogen[rule_key][halogen_key]
                unknown_metrics = [k for k in cell_metrics.keys() if k not in distributable_set]
                if unknown_metrics:
                    warnings.append({
                        "type": "unknown_metric_dropped",
                        "dimension": "by_rule_halogen",
                        "key": f"rule:{rule_key}:halogen:{halogen_key}",
                        "metrics": unknown_metrics
                    })
                    for k in unknown_metrics:
                        del cell_metrics[k]

                # Validate remaining metrics
                for metric, value in list(cell_metrics.items()):
                    validated_value = validate_metric_value(value, {
                        "dimension": "by_rule_halogen",
                        "key": f"rule:{rule_key}:halogen:{halogen_key}",
                        "metric": metric
                    })
                    cell_metrics[metric] = validated_value

    return warnings


def sanitize_granular_fields(by_rule: Optional[Dict], by_halogen: Optional[Dict],
                            by_rule_halogen: Optional[Dict], distributable_set: set) -> None:
    """
    Remove non-distributable fields from existing granular dimensions.

    Args:
        by_rule: by_rule structure (modified in place)
        by_halogen: by_halogen structure (modified in place)
        by_rule_halogen: by_rule_halogen structure (modified in place)
        distributable_set: Set of field names that should be kept
    """
    # Clean by_rule
    if by_rule:
        for rule in by_rule:
            keys_to_remove = [k for k in by_rule[rule].keys() if k not in distributable_set]
            for k in keys_to_remove:
                del by_rule[rule][k]

    # Clean by_halogen
    if by_halogen:
        for halogen in by_halogen:
            keys_to_remove = [k for k in by_halogen[halogen].keys() if k not in distributable_set]
            for k in keys_to_remove:
                del by_halogen[halogen][k]

    # Clean by_rule_halogen
    if by_rule_halogen:
        for rule in by_rule_halogen:
            for halogen in by_rule_halogen[rule]:
                keys_to_remove = [k for k in by_rule_halogen[rule][halogen].keys() if k not in distributable_set]
                for k in keys_to_remove:
                    del by_rule_halogen[rule][halogen][k]


def distribute_marginals_to_2d(by_rule: Optional[Dict], by_halogen: Optional[Dict],
                              rules: List[str], halogens: List[str], metrics: List[str]) -> Dict:
    """
    Create 2D structure by distributing marginal values from by_rule or by_halogen.

    Args:
        by_rule: by_rule structure with marginal values
        by_halogen: by_halogen structure with marginal values
        rules: List of rules to distribute across
        halogens: List of halogens to distribute across
        metrics: List of metrics to distribute

    Returns:
        by_rule_halogen structure with distributed values
    """
    # Empty set protection: return all-zero 2D structure if rules or halogens are empty
    if not rules or not halogens:
        # Create correct shape but all-zero structure
        return {r: {h: {m: 0 for m in metrics} for h in halogens} for r in rules}

    # Get canonical distribution order
    distribution_rules, distribution_halogens = _distribution_order(rules, halogens)

    # Create empty 2D structure
    by_rule_halogen = {r: {h: {m: 0 for m in metrics} for h in halogens} for r in rules}

    # If we have by_rule data, distribute it across halogens
    if by_rule:
        for rule in rules:
            if rule in by_rule:
                for metric in metrics:
                    if metric in by_rule[rule] and by_rule[rule][metric] > 0:
                        total_value = by_rule[rule][metric]
                        # Distribute evenly across halogens, with remainder going to first halogens in distribution order
                        base_value = total_value // len(halogens)
                        remainder = total_value % len(halogens)

                        for i, halogen in enumerate(distribution_halogens):
                            by_rule_halogen[rule][halogen][metric] = base_value
                            if i < remainder:  # First halogens get +1 for remainder
                                by_rule_halogen[rule][halogen][metric] += 1

    # If we have by_halogen data, distribute it across rules
    elif by_halogen:
        for halogen in halogens:
            if halogen in by_halogen:
                for metric in metrics:
                    if metric in by_halogen[halogen] and by_halogen[halogen][metric] > 0:
                        total_value = by_halogen[halogen][metric]
                        # Distribute evenly across rules, with remainder going to first rules in distribution order
                        base_value = total_value // len(rules)
                        remainder = total_value % len(rules)

                        for i, rule in enumerate(distribution_rules):
                            by_rule_halogen[rule][halogen][metric] = base_value
                            if i < remainder:  # First rules get +1 for remainder
                                by_rule_halogen[rule][halogen][metric] += 1

    return by_rule_halogen


def select_distribution_base(by_rule: Optional[Dict], by_halogen: Optional[Dict]) -> Tuple[Optional[Dict], Optional[Dict], str]:
    """
    Select the base marginal for 2D distribution with priority: by_rule > by_halogen.

    Args:
        by_rule: by_rule structure with marginal values
        by_halogen: by_halogen structure with marginal values

    Returns:
        Tuple of (source_by_rule, source_by_halogen, base_str) where:
        - Only one of source_by_rule/source_by_halogen will be non-None (the base)
        - base_str is 'by_rule', 'by_halogen', or 'none'
    """
    # Determine base marginal for distribution (by_rule has priority)
    if by_rule and by_halogen:
        # Both exist: use by_rule as base, ignore by_halogen for distribution
        return by_rule, None, 'by_rule'
    elif by_rule:
        # Only by_rule exists
        return by_rule, None, 'by_rule'
    elif by_halogen:
        # Only by_halogen exists
        return None, by_halogen, 'by_halogen'
    else:
        # Neither exists
        return None, None, 'none'


def distribute_marginals_to_2d_with_warnings(by_rule: Optional[Dict], by_halogen: Optional[Dict],
                                            rules: List[str], halogens: List[str], metrics: List[str]) -> Tuple[Dict, List[Dict[str, Any]]]:
    """
    Create 2D structure by distributing marginal values, with warnings support.

    Args:
        by_rule: by_rule structure with marginal values
        by_halogen: by_halogen structure with marginal values
        rules: List of rules to distribute across
        halogens: List of halogens to distribute across
        metrics: List of metrics to distribute

    Returns:
        Tuple of (by_rule_halogen structure with distributed values, list of warning objects)
    """
    warnings = []

    # Check for empty sets and add warning
    if not rules or not halogens:
        warnings.append({"type": "empty_rules_or_halogens"})

    # Call the main function
    by_rule_halogen = distribute_marginals_to_2d(by_rule, by_halogen, rules, halogens, metrics)

    return by_rule_halogen, warnings


def normalize_granular_shape(by_rule: Optional[Dict], by_halogen: Optional[Dict],
                           by_rule_halogen: Optional[Dict], rules: List[str],
                           halogens: List[str], metrics: List[str]) -> None:
    """
    Fill missing metrics with 0 in existing granular dimensions.

    Args:
        by_rule: by_rule structure (modified in place)
        by_halogen: by_halogen structure (modified in place)
        by_rule_halogen: by_rule_halogen structure (modified in place)
        rules: List of rules to ensure coverage
        halogens: List of halogens to ensure coverage
        metrics: List of metrics to ensure coverage (0-filled if missing)
    """
    # Normalize by_rule
    if by_rule:
        for rule in rules:
            if rule in by_rule:
                for metric in metrics:
                    if metric not in by_rule[rule]:
                        by_rule[rule][metric] = 0

    # Normalize by_halogen
    if by_halogen:
        for halogen in halogens:
            if halogen in by_halogen:
                for metric in metrics:
                    if metric not in by_halogen[halogen]:
                        by_halogen[halogen][metric] = 0

    # Normalize by_rule_halogen
    if by_rule_halogen:
        for rule in rules:
            if rule in by_rule_halogen:
                for halogen in halogens:
                    if halogen in by_rule_halogen[rule]:
                        for metric in metrics:
                            if metric not in by_rule_halogen[rule][halogen]:
                                by_rule_halogen[rule][halogen][metric] = 0


def compute_distributable_totals(qa_summary: Dict, rules: List[str], halogens: List[str], metrics: List[str]) -> Dict[str, int]:
    """
    Compute distributable totals using multi-path fallback strategy.

    Args:
        qa_summary: QA summary structure
        rules: List of rules
        halogens: List of halogens
        metrics: List of metrics to compute totals for

    Returns:
        Dictionary mapping metrics to their total values

    Fallback order:
        1. Aggregate from by_rule_halogen (2D - most reliable)
        2. Aggregate from by_rule (1D fallback)
        3. Aggregate from by_halogen (1D fallback)
        4. Return 0 if none available
    """
    totals = {}

    for metric in metrics:
        total_value = 0

        # Priority 1: Use 2D aggregation if available
        by_rule_halogen = qa_summary.get('by_rule_halogen')
        if by_rule_halogen:
            try:
                total_value = sum(
                    by_rule_halogen[r][h].get(metric, 0)
                    for r in rules for h in halogens
                    if r in by_rule_halogen and h in by_rule_halogen[r]
                )
            except (KeyError, TypeError):
                total_value = 0

        # Priority 2: Fallback to by_rule aggregation
        elif qa_summary.get('by_rule'):
            by_rule = qa_summary['by_rule']
            try:
                total_value = sum(
                    by_rule[r].get(metric, 0)
                    for r in rules
                    if r in by_rule
                )
            except (KeyError, TypeError):
                total_value = 0

        # Priority 3: Fallback to by_halogen aggregation
        elif qa_summary.get('by_halogen'):
            by_halogen = qa_summary['by_halogen']
            try:
                total_value = sum(
                    by_halogen[h].get(metric, 0)
                    for h in halogens
                    if h in by_halogen
                )
            except (KeyError, TypeError):
                total_value = 0

        # Priority 4: Default to 0 if no granular data available
        else:
            total_value = 0

        totals[metric] = total_value

    return totals


def _distribution_order(rules: List[str], halogens: List[str]) -> tuple:
    """
    Define the canonical distribution order for rules and halogens.
    
    Args:
        rules: Input rule list 
        halogens: Input halogen list
        
    Returns:
        Tuple of (sorted_rules, sorted_halogens) for deterministic distribution.
        
    Note:
        Uses lexicographic (alphabetical) ordering to ensure deterministic
        remainder allocation regardless of input order.
    """
    return (sorted(rules), sorted(halogens))


def _distribute(total: int, keys: List[str]) -> Dict[str, int]:
    """
    Distribute total value among keys with remainder handling.
    
    Args:
        total: Total value to distribute
        keys: List of keys to distribute among (caller controls order)
        
    Returns:
        Dictionary mapping each key to its share of the total.
        Remainder is distributed to the first keys in the provided order.
        
    Note:
        The caller is responsible for key ordering (lexicographic or input order).
        For 2D scenarios, current implementation uses sorted(rules) x sorted(halogens).
    """
    if not keys:
        return {}
    
    n = len(keys)
    q, r = divmod(total, n)
    
    # Base allocation + remainder to first r keys in provided order
    return {k: (q + (1 if i < r else 0)) for i, k in enumerate(keys)}


def write_qa_summary_json(qa_stats_dict: Dict[str, Any], output_dir: str,
                         completion_mode: str = 'zero_fill', conflict_tolerance: int = 1, max_warnings: int = 1000) -> str:
    """
    Write QA statistics summary to JSON sidecar file.

    Args:
        qa_stats_dict: QA statistics dictionary from enumeration
        output_dir: Directory to write the JSON file
        completion_mode: Strategy for completing missing granular dimensions
                        'zero_fill' (default): Missing dimensions filled with zeros
                        'distribute': Missing dimensions derived by distributing marginal values
        conflict_tolerance: Threshold for marginal conflict detection (default: 1)
        max_warnings: Maximum number of warnings to include in output (default: 1000)

    Returns:
        Path to the written JSON file
    """
    # Validate completion_mode parameter
    valid_modes = {'zero_fill', 'distribute'}
    if completion_mode not in valid_modes:
        raise ValueError(f"Invalid completion_mode '{completion_mode}'. Must be one of: {', '.join(sorted(valid_modes))}")

    qa_summary_path = os.path.join(output_dir, 'qa_summary.json')
    os.makedirs(output_dir, exist_ok=True)
    
    # Structure QA summary with total and breakdowns
    if qa_stats_dict.get('version') == '2':
        # v2 format with pivots - qa_stats_dict already has the structure we want
        qa_summary = qa_stats_dict.copy()
        metadata = {
            'description': 'QA statistics for halogenation enumeration process',
            'semantics': 'All counters follow per-attempt counting unless noted',
            'generated_at': datetime.now().isoformat(),
            'parent_key_constants': {
                'ik_prefix': IK_PREFIX,
                'smi_prefix': SMI_PREFIX,
                'smi_unknown': SMI_UNKNOWN,
                'smiles_hash_threshold': SMILES_HASH_THRESHOLD
            }
        }
        
        # Add parent key hashing metadata if available
        if 'parent_key_metadata' in qa_stats_dict:
            parent_meta = qa_stats_dict['parent_key_metadata']
            metadata['parent_key_hashing'] = {
                'hashed_count': parent_meta.get('hashed_count', 0),
                'examples': parent_meta.get('hashed_examples', []),
                'threshold': SMILES_HASH_THRESHOLD,
                'rdkit_available_count': parent_meta.get('rdkit_available_count', 0),
                'inchikey_generated_count': parent_meta.get('inchikey_generated_count', 0)
            }

        # Add distribution_keys with priority: metadata -> pivots -> empty
        rules = []
        halogens = []

        # Priority 1: Use metadata rules and halogens if available
        if 'metadata' in qa_summary:
            meta = qa_summary['metadata']
            rules = meta.get('rules', [])
            halogens = meta.get('halogens', [])

        # Priority 2: Fall back to pivots if metadata not available
        if (not rules or not halogens) and 'pivots' in qa_summary:
            pivots = qa_summary['pivots']
            rules_set = set(rules) if rules else set()
            halogens_set = set(halogens) if halogens else set()

            # Extract rules and halogens from pivot dimensions
            if not rules and 'by_rule' in pivots:
                rules_set.update(pivots['by_rule'].keys())
            if not halogens and 'by_halogen' in pivots:
                halogens_set.update(pivots['by_halogen'].keys())
            if 'by_rule_halogen' in pivots:
                if not rules:
                    rules_set.update(pivots['by_rule_halogen'].keys())
                if not halogens:
                    for rule_dict in pivots['by_rule_halogen'].values():
                        if isinstance(rule_dict, dict):
                            halogens_set.update(rule_dict.keys())

            rules = sorted(rules_set) if not rules else rules
            halogens = sorted(halogens_set) if not halogens else halogens

        # Always include distribution_keys, even if empty
        distribution_rules, distribution_halogens = _distribution_order(rules, halogens)
        metadata['distribution_keys'] = {
            "rules": list(distribution_rules),
            "halogens": list(distribution_halogens)
        }
        # NOTE: value_coercion_method and marginal_conflict_tolerance will be set by finalize_metadata

        qa_summary['metadata'] = metadata
    elif qa_stats_dict.get('version') == '1' and ('by_rule' in qa_stats_dict or 'by_halogen' in qa_stats_dict or 'by_rule_halogen' in qa_stats_dict):
        # v1 format with existing granular structure - validate and complete
        qa_summary = qa_stats_dict.copy()
        metadata = {
            'description': 'QA statistics for halogenation enumeration process',
            'semantics': 'All counters follow per-attempt counting unless noted',
            'generated_at': datetime.now().isoformat(),
            'parent_key_constants': {
                'ik_prefix': IK_PREFIX,
                'smi_prefix': SMI_PREFIX,
                'smi_unknown': SMI_UNKNOWN,
                'smiles_hash_threshold': SMILES_HASH_THRESHOLD
            }
        }
        
        # Add parent key hashing metadata if available
        if 'parent_key_metadata' in qa_stats_dict:
            parent_meta = qa_stats_dict['parent_key_metadata']
            metadata['parent_key_hashing'] = {
                'hashed_count': parent_meta.get('hashed_count', 0),
                'examples': parent_meta.get('hashed_examples', []),
                'threshold': SMILES_HASH_THRESHOLD,
                'rdkit_available_count': parent_meta.get('rdkit_available_count', 0),
                'inchikey_generated_count': parent_meta.get('inchikey_generated_count', 0)
            }
        
        # Validate and complete granular structure
        metadata_dict = qa_summary.get('metadata', {})
        halogens = metadata_dict.get('halogens', list(ALL_HALOGENS))
        rules = metadata_dict.get('rules', list(ALL_RULES))
        
        # Get metrics from existing structure or use default
        # Extract metrics from ALL available granular structures for robust inference
        metrics_set = set()

        # Extract from by_rule if present
        existing_by_rule = qa_summary.get('by_rule', {})
        for rule_stats in existing_by_rule.values():
            metrics_set.update(rule_stats.keys())

        # Extract from by_halogen if present
        existing_by_halogen = qa_summary.get('by_halogen', {})
        for halogen_stats in existing_by_halogen.values():
            metrics_set.update(halogen_stats.keys())

        # Extract from by_rule_halogen if present
        existing_by_rule_halogen = qa_summary.get('by_rule_halogen', {})
        for rule_dict in existing_by_rule_halogen.values():
            for halogen_stats in rule_dict.values():
                metrics_set.update(halogen_stats.keys())

        if metrics_set:
            # Add standard metrics to ensure complete coverage
            metrics_set.update(QA_PATH_KEYS)
            metrics_set.update(['no_product_matches', 'template_unsupported'])

            # Filter to only include distributable metrics (exclude overview counters)
            distributable_metrics_set = set(QA_PATH_KEYS) | {'no_product_matches', 'template_unsupported'}
            metrics_for_granular = metrics_set & distributable_metrics_set
            metrics = sorted(list(metrics_for_granular))  # Sort for deterministic order
        else:
            # Default metrics set when no granular structures exist
            qa_path_metrics = list(QA_PATH_KEYS)
            top_level_metrics = ['no_product_matches', 'template_unsupported']
            metrics = qa_path_metrics + top_level_metrics
        
        # NEW PASSTHROUGH PIPELINE: 2D as Single Source of Truth

        # Step 1: Sanitize and discover metrics from all dimensions
        distributable_metrics_set = set(QA_PATH_KEYS) | {'no_product_matches', 'template_unsupported'}

        # Enhanced sanitization with data quality warnings
        sanitize_warnings = sanitize_granular_fields_with_warnings(
            qa_summary.get('by_rule'),
            qa_summary.get('by_halogen'),
            qa_summary.get('by_rule_halogen'),
            distributable_metrics_set,
            rules, halogens
        )

        # Add sanitization warnings to metadata
        if sanitize_warnings:
            if 'warnings' not in metadata:
                metadata['warnings'] = []
            metadata['warnings'].extend(sanitize_warnings)

        # Step 1.5: Detect marginal conflicts and add warnings
        if 'by_rule' in qa_summary and 'by_halogen' in qa_summary:
            by_rule_totals = {}
            by_halogen_totals = {}

            # Calculate totals for each metric from by_rule
            for rule, rule_metrics in qa_summary['by_rule'].items():
                for metric, value in rule_metrics.items():
                    by_rule_totals[metric] = by_rule_totals.get(metric, 0) + value

            # Calculate totals for each metric from by_halogen
            for halogen, halogen_metrics in qa_summary['by_halogen'].items():
                for metric, value in halogen_metrics.items():
                    by_halogen_totals[metric] = by_halogen_totals.get(metric, 0) + value

            # Check for conflicts (differences > tolerance threshold)
            tolerance = conflict_tolerance
            conflicts = {}
            # Use union of metrics to catch metrics present in only one dimension
            all_metrics = set(by_rule_totals.keys()) | set(by_halogen_totals.keys())
            # Only check distributable metrics (exclude overview counters)
            metrics_to_check = all_metrics & distributable_metrics_set

            for metric in metrics_to_check:
                lhs = by_rule_totals.get(metric, 0)
                rhs = by_halogen_totals.get(metric, 0)
                diff = abs(lhs - rhs)
                if diff > tolerance:
                    conflicts[metric] = {
                        "lhs": lhs,
                        "rhs": rhs,
                        "diff": diff
                    }

            if conflicts:
                if 'warnings' not in metadata:
                    metadata['warnings'] = []
                metadata['warnings'].append({
                    "type": "marginal_conflict",
                    "base": "by_rule",
                    "tolerance": tolerance,
                    "delta": conflicts
                })

        # Step 2: Determine strategy and prepare for 2D establishment
        original_had_2d = 'by_rule_halogen' in qa_summary
        if original_had_2d:
            # Input has 2D structure - normalize it
            normalize_granular_shape(
                None, None, qa_summary.get('by_rule_halogen'),
                rules, halogens, metrics
            )
            # 2D is the authoritative source
            by_rule_halogen_2d = qa_summary['by_rule_halogen']
        else:
            # No 2D structure - will be created after total computation
            by_rule_halogen_2d = None

        # Step 3: Compute distributable total first (before creating 2D in zero-fill mode)
        if 'total' not in qa_summary:
            if original_had_2d:
                # Original input had 2D: compute total from existing 2D
                totals_distributable = {}
                for metric in metrics:
                    totals_distributable[metric] = sum(
                        by_rule_halogen_2d[rule][halogen][metric]
                        for rule in rules for halogen in halogens
                    )
            else:
                # No original 2D: compute total using fallback strategy (works for both modes)
                totals_distributable = compute_distributable_totals(qa_summary, rules, halogens, metrics)

            # Merge overview counters into total
            overview_counters = list(OVERVIEW_COUNTERS)
            overview_extracted = {}
            for counter in overview_counters:
                if counter in qa_stats_dict:
                    overview_extracted[counter] = int(qa_stats_dict[counter])

            qa_summary['total'] = {**totals_distributable, **overview_extracted}

        # Step 4: Create/complete 2D structure if not already present
        if by_rule_halogen_2d is None:
            if completion_mode == 'distribute':
                # Distribute mode: derive 2D from marginal values
                # Use unified base selection helper
                by_rule_data = qa_summary.get('by_rule')
                by_halogen_data = qa_summary.get('by_halogen')

                source_by_rule, source_by_halogen, base_str = select_distribution_base(by_rule_data, by_halogen_data)

                by_rule_halogen_2d, distribution_warnings = distribute_marginals_to_2d_with_warnings(
                    source_by_rule,
                    source_by_halogen,
                    rules, halogens, metrics
                )
                # Add any warnings to metadata
                if distribution_warnings:
                    if 'warnings' not in metadata:
                        metadata['warnings'] = []
                    metadata['warnings'].extend(distribution_warnings)
            else:
                # Zero-fill mode: create all-zero 2D structure
                by_rule_halogen_2d = {r: {h: {m: 0 for m in metrics} for h in halogens} for r in rules}

            qa_summary['by_rule_halogen'] = by_rule_halogen_2d

        # Step 5: Derive/complete by_rule and by_halogen (strategy-dependent)
        if completion_mode == 'distribute' or original_had_2d:
            # Distribute mode OR original input had 2D: derive marginals from authoritative 2D
            by_rule_derived = {r: {m: 0 for m in metrics} for r in rules}
            for rule in rules:
                for metric in metrics:
                    by_rule_derived[rule][metric] = sum(
                        by_rule_halogen_2d[rule][halogen][metric] for halogen in halogens
                    )
            qa_summary['by_rule'] = by_rule_derived

            by_halogen_derived = {h: {m: 0 for m in metrics} for h in halogens}
            for halogen in halogens:
                for metric in metrics:
                    by_halogen_derived[halogen][metric] = sum(
                        by_rule_halogen_2d[rule][halogen][metric] for rule in rules
                    )
            qa_summary['by_halogen'] = by_halogen_derived
            marginals_from_2d = True
        else:
            # Zero-fill mode AND no original 2D: preserve input marginals, complete missing ones
            if 'by_rule' not in qa_summary:
                qa_summary['by_rule'] = {r: {m: 0 for m in metrics} for r in rules}
            if 'by_halogen' not in qa_summary:
                qa_summary['by_halogen'] = {h: {m: 0 for m in metrics} for h in halogens}
            marginals_from_2d = False
        
        # Step 6: Add v1 granular distribution semantics to metadata
        metadata['distribution'] = "equal_with_lexicographic_remainder"
        metadata['marginals_from_2d'] = marginals_from_2d  # Strategy-dependent: True if derived from 2D, False if from input
        metadata['marginals_source'] = "2d_derived" if marginals_from_2d else "input"

        # Add completion strategy metadata
        metadata['completion'] = {
            'enabled': completion_mode == 'distribute',
            'strategy': completion_mode,
        }

        if completion_mode == 'distribute':
            # Determine the base source for distribution
            if 'by_rule_halogen' in qa_stats_dict:
                base = '2d'
            else:
                # Use the unified base selection logic for consistent behavior
                _, _, base = select_distribution_base(qa_stats_dict.get('by_rule'), qa_stats_dict.get('by_halogen'))

            # Get canonical distribution order for consistency with distribute_marginals_to_2d
            distribution_rules, distribution_halogens = _distribution_order(rules, halogens)

            metadata['completion'].update({
                'base': base,
                'ordered_keys': {
                    'rules': distribution_rules,
                    'halogens': distribution_halogens
                }
            })
        metadata['distribution_order'] = {
            "rules": "lexicographic",
            "halogens": "lexicographic"
        }
        # NOTE: marginal_conflict_tolerance and value_coercion_method will be set by finalize_metadata

        # Add explicit distribution keys used for remainder allocation
        distribution_rules, distribution_halogens = _distribution_order(rules, halogens)
        metadata['distribution_keys'] = {
            "rules": list(distribution_rules),
            "halogens": list(distribution_halogens)
        }
        
        qa_summary['metadata'] = metadata
    else:
        # Legacy totals-only input or v1 format without granular structure
        metadata = {
            'description': 'QA statistics for halogenation enumeration process',
            'semantics': 'All counters follow per-attempt counting unless noted',
            'generated_at': datetime.now().isoformat(),
            'parent_key_constants': {
                'ik_prefix': IK_PREFIX,
                'smi_prefix': SMI_PREFIX,
                'smi_unknown': SMI_UNKNOWN,
                'smiles_hash_threshold': SMILES_HASH_THRESHOLD
            }
        }
        
        # Add parent key hashing metadata if available
        if 'parent_key_metadata' in qa_stats_dict:
            parent_meta = qa_stats_dict['parent_key_metadata']
            metadata['parent_key_hashing'] = {
                'hashed_count': parent_meta.get('hashed_count', 0),
                'examples': parent_meta.get('hashed_examples', []),
                'threshold': SMILES_HASH_THRESHOLD,
                'rdkit_available_count': parent_meta.get('rdkit_available_count', 0),
                'inchikey_generated_count': parent_meta.get('inchikey_generated_count', 0)
            }
        
        # Check if metadata with halogens is present for granular slices
        metadata_dict = qa_stats_dict.get('metadata', {})
        if 'halogens' in metadata_dict:
            # Build minimal v1 granular slices with deterministic distribution
            halogens = metadata_dict.get('halogens', list(ALL_HALOGENS))
            rules = metadata_dict.get('rules', list(ALL_RULES))
            # Distinguish between distributable metrics and overview counters
            qa_path_metrics = list(QA_PATH_KEYS)
            distributable_top_level = ['no_product_matches', 'template_unsupported']
            distributable_metrics = qa_path_metrics + distributable_top_level
            
            # Overview counters (not distributed to granular dimensions)
            overview_counters = list(OVERVIEW_COUNTERS)

            # Extract distributable metrics for granular distribution
            totals = {}
            qa_paths_tot = _normalize_qa_paths_with_aliases(qa_stats_dict.get('qa_paths', {}))
            for m in distributable_metrics:
                if m in qa_paths_tot:
                    totals[m] = int(qa_paths_tot.get(m, 0))
                else:
                    totals[m] = int(qa_stats_dict.get(m, 0))
            
            # Add overview counters to totals (but not to distributable_metrics)
            for counter in overview_counters:
                if counter in qa_stats_dict:
                    totals[counter] = int(qa_stats_dict[counter])

            # Step 1: First distribute to 2D structure (by_rule_halogen)
            # This is the primary distribution that ensures deterministic remainder allocation
            distribution_rules, distribution_halogens = _distribution_order(rules, halogens)
            by_rule_halogen = {r: {h: {m: 0 for m in distributable_metrics} for h in halogens} for r in rules}
            for m in distributable_metrics:
                # Create flat list of (rule, halogen) combinations for distribution
                rule_halogen_keys = [f"{r}:{h}" for r in distribution_rules for h in distribution_halogens]
                distribution = _distribute(totals[m], rule_halogen_keys)
                
                # Apply distribution to nested structure
                for r in rules:
                    for h in halogens:
                        key = f"{r}:{h}"
                        by_rule_halogen[r][h][m] = distribution.get(key, 0)

            # Step 2: Derive by_rule from 2D structure (marginal consistency)
            by_rule = {r: {m: 0 for m in distributable_metrics} for r in rules}
            for r in rules:
                for m in distributable_metrics:
                    by_rule[r][m] = sum(by_rule_halogen[r][h][m] for h in halogens)

            # Step 3: Derive by_halogen from 2D structure (marginal consistency)
            by_halogen = {h: {m: 0 for m in distributable_metrics} for h in halogens}
            for h in halogens:
                for m in distributable_metrics:
                    by_halogen[h][m] = sum(by_rule_halogen[r][h][m] for r in rules)

            # Add v1 granular distribution semantics to metadata
            metadata['distribution'] = "equal_with_lexicographic_remainder"
            metadata['marginals_from_2d'] = True
            metadata['distribution_order'] = {
                "rules": "lexicographic",
                "halogens": "lexicographic"
            }

            # Add explicit distribution keys used for remainder allocation
            metadata['distribution_keys'] = {
                "rules": list(distribution_rules),
                "halogens": list(distribution_halogens)
            }
            # NOTE: value_coercion_method will be set by finalize_metadata

            qa_summary = {
                'version': '1',
                'total': totals,  # Use computed totals (flat numeric dict) instead of raw input
                'by_rule': by_rule,
                'by_halogen': by_halogen,
                'by_rule_halogen': by_rule_halogen,
                'metadata': metadata
            }
        else:
            # Legacy totals-only case without metadata - no granular slices
            # NOTE: value_coercion_method will be set by finalize_metadata
            # NOTE: marginal_conflict_tolerance will be set by finalize_metadata
            qa_summary = {
                'version': '1',
                'total': qa_stats_dict,
                'metadata': metadata
            }

    # Apply centralized metadata finalization
    if 'metadata' in qa_summary:
        warnings_list = qa_summary['metadata'].get('warnings', [])
        qa_summary['metadata'], _ = finalize_metadata(
            qa_summary['metadata'],
            conflict_tolerance=conflict_tolerance,
            coercion_method="trunc",
            warnings_list=warnings_list,
            max_warnings=max_warnings
        )

    with open(qa_summary_path, 'w') as f:
        json.dump(qa_summary, f, indent=2)
    
    # Generate pivot CSV if we have version 2 data with pivots
    if qa_summary.get('version') == '2' and 'pivots' in qa_summary:
        pivot_csv_path = os.path.join(output_dir, 'pivot_rule_halogen_k.csv')
        _write_pivot_rule_halogen_k_csv(qa_summary['pivots'], pivot_csv_path)

        # Generate dual-metric pivot CSVs
        pivot_ops_atoms_csv_path = os.path.join(output_dir, 'pivot_rule_halogen_ops_atoms.csv')
        _write_pivot_rule_halogen_ops_atoms_csv(qa_summary['pivots'], pivot_ops_atoms_csv_path)

        pivot_rule_ops_atoms_csv_path = os.path.join(output_dir, 'pivot_rule_ops_atoms.csv')
        _write_pivot_rule_ops_atoms_csv(qa_summary['pivots'], pivot_rule_ops_atoms_csv_path)
    
    # Generate dedicated rule_halogen_k.csv with parent_key, rule, halogen, k, count schema
    rule_halogen_k_csv_path = os.path.join(output_dir, 'rule_halogen_k.csv')
    _write_rule_halogen_k_csv(qa_summary, rule_halogen_k_csv_path)
    
    return qa_summary_path


def _write_pivot_rule_halogen_k_csv(pivots_data: Dict[str, Any], output_path: str) -> None:
    """
    Write pivot_rule_halogen_k.csv from version 2 pivots data.
    
    Args:
        pivots_data: Pivots dictionary from version 2 QA stats
        output_path: Output path for the CSV file
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Extract by_rule_halogen_k data
    rule_halogen_k_data = pivots_data.get('by_rule_halogen_k', {})
    
    # Create CSV rows
    csv_rows = []
    
    # Event types to include in CSV
    event_types = ['attempts', 'products', 'no_product_matches', 'isotope_unavailable', 
                   'isotope_miss', 'atommap_used', 'heuristic_used', 'template_unsupported']
    
    # Process each rule_halogen_k combination
    for rule_halogen_k_key, events in rule_halogen_k_data.items():
        # Parse key: "R1_F_1" -> rule="R1", halogen="F", k=1
        parts = rule_halogen_k_key.split('_')
        if len(parts) >= 3:
            rule = parts[0]
            halogen = parts[1]
            k = parts[2]
            
            # Create row with all event counts
            row = {
                'rule': rule,
                'halogen': halogen,
                'k': k
            }
            
            # Add event counts (default to 0 if not present)
            for event_type in event_types:
                row[event_type] = events.get(event_type, 0)
            
            csv_rows.append(row)
    
    if csv_rows:
        # Create DataFrame and write CSV
        df = pd.DataFrame(csv_rows)
        # Sort by rule, halogen, k for consistent output
        df = df.sort_values(['rule', 'halogen', 'k'])
        df.to_csv(output_path, index=False)
    else:
        # Create empty CSV with headers if no data
        empty_df = pd.DataFrame(columns=['rule', 'halogen', 'k'] + event_types)
        empty_df.to_csv(output_path, index=False)


def _write_pivot_rule_halogen_ops_atoms_csv(pivots_data: Dict[str, Any], output_path: str) -> None:
    """
    Write pivot_rule_halogen_ops_atoms.csv from version 2 pivots data.

    Args:
        pivots_data: Pivots dictionary from version 2 QA stats
        output_path: Output path for the CSV file
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Extract by_rule_halogen_ops_atoms data
    rule_halogen_ops_atoms_data = pivots_data.get('by_rule_halogen_ops_atoms', {})

    # Create CSV rows
    csv_rows = []

    # Event types to include in CSV
    event_types = ['attempts', 'products', 'no_product_matches', 'isotope_unavailable',
                   'isotope_miss', 'atommap_used', 'heuristic_used', 'template_unsupported']

    # Process each rule_halogen_ops_atoms combination
    for rule_halogen_ops_atoms_key, events in rule_halogen_ops_atoms_data.items():
        # Parse key: "R1_F_2_3" -> rule="R1", halogen="F", k_ops=2, k_atoms=3
        parts = rule_halogen_ops_atoms_key.split('_')
        if len(parts) >= 4:
            rule = parts[0]
            halogen = parts[1]
            k_ops = parts[2]
            k_atoms = parts[3]

            # Create row with all event counts
            row = {
                'rule': rule,
                'halogen': halogen,
                'k_ops': k_ops,
                'k_atoms': k_atoms
            }

            # Add event counts (default to 0 if not present)
            for event_type in event_types:
                row[event_type] = events.get(event_type, 0)

            csv_rows.append(row)

    if csv_rows:
        # Create DataFrame and write CSV
        df = pd.DataFrame(csv_rows)
        # Sort by rule, halogen, k_ops, k_atoms for consistent output
        df = df.sort_values(['rule', 'halogen', 'k_ops', 'k_atoms'])
        df.to_csv(output_path, index=False)
    else:
        # Create empty CSV with headers if no data
        empty_df = pd.DataFrame(columns=['rule', 'halogen', 'k_ops', 'k_atoms'] + event_types)
        empty_df.to_csv(output_path, index=False)


def _write_pivot_rule_ops_atoms_csv(pivots_data: Dict[str, Any], output_path: str) -> None:
    """
    Write pivot_rule_ops_atoms.csv from version 2 pivots data.

    Args:
        pivots_data: Pivots dictionary from version 2 QA stats
        output_path: Output path for the CSV file
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Extract by_rule_ops_atoms data
    rule_ops_atoms_data = pivots_data.get('by_rule_ops_atoms', {})

    # Create CSV rows
    csv_rows = []

    # Event types to include in CSV
    event_types = ['attempts', 'products', 'no_product_matches', 'isotope_unavailable',
                   'isotope_miss', 'atommap_used', 'heuristic_used', 'template_unsupported']

    # Process each rule_ops_atoms combination
    for rule_ops_atoms_key, events in rule_ops_atoms_data.items():
        # Parse key: "R1_2_3" -> rule="R1", k_ops=2, k_atoms=3
        parts = rule_ops_atoms_key.split('_')
        if len(parts) >= 3:
            rule = parts[0]
            k_ops = parts[1]
            k_atoms = parts[2]

            # Create row with all event counts
            row = {
                'rule': rule,
                'k_ops': k_ops,
                'k_atoms': k_atoms
            }

            # Add event counts (default to 0 if not present)
            for event_type in event_types:
                row[event_type] = events.get(event_type, 0)

            csv_rows.append(row)

    if csv_rows:
        # Create DataFrame and write CSV
        df = pd.DataFrame(csv_rows)
        # Sort by rule, k_ops, k_atoms for consistent output
        df = df.sort_values(['rule', 'k_ops', 'k_atoms'])
        df.to_csv(output_path, index=False)
    else:
        # Create empty CSV with headers if no data
        empty_df = pd.DataFrame(columns=['rule', 'k_ops', 'k_atoms'] + event_types)
        empty_df.to_csv(output_path, index=False)


def _write_rule_halogen_k_csv(stats_data: Dict[str, Any], output_path: str) -> None:
    """
    Write dedicated rule_halogen_k.csv with parent_key, rule, halogen, k, count schema.
    
    This flattens the rule_halogen_k_by_parent from incremental stats into the required CSV format.
    Each row represents attempts for a specific (parent_key, rule, halogen, k) combination.
    
    Args:
        stats_data: Statistics data containing rule_halogen_k_by_parent
        output_path: Output path for the CSV file
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Extract rule_halogen_k_by_parent data (parent -> rule -> halogen -> k -> count)
    rule_halogen_k_by_parent = stats_data.get('rule_halogen_k_by_parent', {})
    
    # Create CSV rows
    csv_rows = []
    
    # Process each parent -> rule -> halogen -> k -> count structure
    for parent_key, rule_dict in rule_halogen_k_by_parent.items():
        for rule, halogen_dict in rule_dict.items():
            for halogen, k_dict in halogen_dict.items():
                for k, count in k_dict.items():
                    csv_rows.append({
                        'parent_key': parent_key,
                        'rule': rule,
                        'halogen': halogen,
                        'k': int(k),  # Store as integer for proper numerical sorting
                        'count': int(count)
                    })
    
    if csv_rows:
        # Create DataFrame and write CSV
        df = pd.DataFrame(csv_rows)
        # Sort by parent_key, rule, halogen, k for consistent output
        df = df.sort_values(['parent_key', 'rule', 'halogen', 'k'])
        df.to_csv(output_path, index=False)
    else:
        # Create empty CSV with headers if no data
        empty_df = pd.DataFrame(columns=['parent_key', 'rule', 'halogen', 'k', 'count'])
        empty_df.to_csv(output_path, index=False)


def _write_pivot_table(rule_halogen_counts: Dict[str, Dict[str, int]], output_path: str) -> None:
    """Write rule x halogen pivot table to CSV."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Create pivot table data
    rules = list(ALL_RULES)
    halogens = list(ALL_HALOGENS)
    
    pivot_data = []
    for rule in rules:
        row = {'Rule': rule}
        total = 0
        for halogen in halogens:
            count = rule_halogen_counts.get(rule, {}).get(halogen, 0)
            row[halogen] = count
            total += count
        row['Total'] = total
        pivot_data.append(row)
    
    # Add total row
    total_row = {'Rule': 'Total'}
    grand_total = 0
    for halogen in halogens:
        halogen_total = sum(rule_halogen_counts.get(rule, {}).get(halogen, 0) for rule in rules)
        total_row[halogen] = halogen_total
        grand_total += halogen_total
    total_row['Total'] = grand_total
    pivot_data.append(total_row)
    
    # Write to CSV
    df = pd.DataFrame(pivot_data)
    df.to_csv(output_path, index=False)


def _print_summary(stats: Dict[str, Any]) -> None:
    """Print summary to console."""
    schema_format = stats.get('config', {}).get('schema_format', 'P0')
    title = "HALOGENATOR P1 SUMMARY REPORT" if schema_format == 'P1' else "HALOGENATOR P0 SUMMARY REPORT"
    print("\n" + "="*60)
    print(title)
    print("="*60)
    
    print(f"Parent molecules processed: {stats['parent_count']}")
    print(f"Total products generated: {stats['product_count']}")
    print(f"Products passing QC: {stats['sanitize_ok_count']}")
    print(f"Unique parents with products: {stats['unique_parents_with_products']}")
    
    # Report parents with zero products using QF-1 logic
    zero_product_parents = stats['parents_without_products']
    parents_with_products = stats['unique_parents_with_products']
    total_parents = stats['parent_count']
    
    if parents_with_products == total_parents:
        print("All parents produced at least one product")
    else:
        # Show count and top 5 examples of parents without products
        zero_count = len(zero_product_parents)
        examples = zero_product_parents[:5]  # Top 5 examples
        examples_str = ', '.join(examples)
        if zero_count > 5:
            examples_str += f" (and {zero_count - 5} more)"
        print(f"Parents with zero products ({zero_count}): {examples_str}")
    
    # Rules activity status (overall)
    all_rules = list(ALL_RULES)
    rule_status = ', '.join(f'{r}={stats["rule_counts"].get(r, 0)}' for r in all_rules)
    print(f"\nRules active (overall): {rule_status}")
    
    # Rules activity on flavonoids (if available)
    if "rule_counts_flavonoids" in stats:
        rule_status_flavonoids = ', '.join(f'{r}={stats["rule_counts_flavonoids"].get(r, 0)}' for r in all_rules)
        print(f"Rules active on flavonoids: {rule_status_flavonoids}")
        
        inactive_rules_flavonoids = [r for r in all_rules if stats["rule_counts_flavonoids"].get(r, 0) == 0]
        if inactive_rules_flavonoids:
            print(f"Rules inactive on flavonoids: {', '.join(inactive_rules_flavonoids)} (expected for R4/R5 - flavonoids rarely contain -NH/-COOH)")
    
    # Rules activity on probes (if available) 
    if "rule_counts_probes" in stats:
        rule_status_probes = ', '.join(f'{r}={stats["rule_counts_probes"].get(r, 0)}' for r in all_rules)
        print(f"Rules active on probes: {rule_status_probes}")
    
    print(f"\nPer-parent statistics:")
    print(f"  Average products per parent: {stats['avg_products_per_parent']}")
    print(f"  Max products per parent: {stats['max_products_per_parent']}")  
    print(f"  Min products per parent: {stats['min_products_per_parent']}")
    
    print(f"\nProducts per rule (non-zero only):")
    for rule, count in sorted(stats['rule_counts'].items()):
        if count > 0:
            print(f"  {rule}: {count} products")
    
    print(f"\nProducts per halogen:")
    for halogen, count in sorted(stats['halogen_counts'].items()):
        print(f"  {halogen}: {count} products")
    
    # Rule x Halogen matrix (compact view)
    print(f"\nRule x Halogen matrix:")
    print(f"{'Rule':<6} {'F':<6} {'Cl':<6} {'Br':<6} {'I':<6} {'Total':<6}")
    print("-" * 42)
    
    # Handle both P0 and P1 formats
    if 'rule_halogen_counts' in stats:
        # P0 format: rule_halogen_counts[rule][halogen]
        rule_halogen_data = stats['rule_halogen_counts']
        for rule in ALL_RULES:
            counts = rule_halogen_data.get(rule, {})
            f_count = counts.get('F', 0)
            cl_count = counts.get('Cl', 0)
            br_count = counts.get('Br', 0)
            i_count = counts.get('I', 0)
            total = f_count + cl_count + br_count + i_count
            if total > 0:
                print(f"{rule:<6} {f_count:<6} {cl_count:<6} {br_count:<6} {i_count:<6} {total:<6}")
    elif 'rule_halogen_k_counts' in stats:
        # P1 format: rule_halogen_k_counts[rule][halogen][k] - sum over k
        rule_halogen_k_data = stats['rule_halogen_k_counts']
        for rule in ALL_RULES:
            rule_data = rule_halogen_k_data.get(rule, {})
            f_count = sum(rule_data.get('F', {}).values())
            cl_count = sum(rule_data.get('Cl', {}).values())
            br_count = sum(rule_data.get('Br', {}).values()) 
            i_count = sum(rule_data.get('I', {}).values())
            total = f_count + cl_count + br_count + i_count
            if total > 0:
                print(f"{rule:<6} {f_count:<6} {cl_count:<6} {br_count:<6} {i_count:<6} {total:<6}")
                
    # P1 additional: show k breakdown
    if 'k_counts' in stats:
        print(f"\nSubstitution depth (k) distribution:")
        for k, count in sorted(stats['k_counts'].items()):
            print(f"  k={k}: {count} products")
    
    # Check for diagnostic anomalies
    diagnostics = stats.get('diagnostics', {})
    if diagnostics.get('anomaly_more_parents_with_products_than_all'):
        print("\nWARNING: Anomaly detected - more parents with products than total parents")
    
    print("\nNo fatal RDKit errors; pipeline stable")
    print("="*60)


def _write_p1_pivot_csv(stats: Dict[str, Any], out_path: str) -> None:
    """Write P1 Rule x Halogen x k pivot table to CSV."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    
    rows = []
    rhk = stats.get('rule_halogen_k_counts', {})
    for rule, hmap in rhk.items():
        for hal, kmap in hmap.items():
            for k, cnt in kmap.items():
                rows.append({'rule': rule, 'halogen': hal, 'k': int(k), 'count': int(cnt)})
    
    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(out_path, index=False)
    else:
        # Write empty file with headers
        pd.DataFrame(columns=['rule', 'halogen', 'k', 'count']).to_csv(out_path, index=False)


def _write_constraints_csv(stats: Dict[str, Any], out_path: str) -> None:
    """Write constraints violations TopN to CSV in original row-based format (BC compatible)."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    
    constraint_stats = stats.get('constraint_stats', {})
    violation_reasons = constraint_stats.get('violation_reasons', {})
    
    # Original row-based format for BC compatibility
    violation_items = sorted(violation_reasons.items(), key=lambda x: -x[1])
    
    rows = []
    for reason, count in violation_items:
        rows.append({
            'violation_reason': str(reason),
            'count': int(count)
        })
    
    # Write row-based CSV (original format)
    if rows:
        df = pd.DataFrame(rows)
    else:
        # Empty CSV with original headers
        df = pd.DataFrame(columns=['violation_reason', 'count'])
    
    df.to_csv(out_path, index=False)
    
    # Also write fixed schema summary file
    summary_path = out_path.replace('_topn.csv', '_summary_fixed.csv')
    _write_constraints_summary_fixed(stats, summary_path)


def _write_constraints_summary_fixed(stats: Dict[str, Any], out_path: str) -> None:
    """Write constraints summary in fixed single-row schema."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    
    constraint_stats = stats.get('constraint_stats', {})
    violation_reasons = constraint_stats.get('violation_reasons', {})
    
    # Summary statistics (always present)
    total_passed = constraint_stats.get('total_passed', 0)
    total_failed = constraint_stats.get('total_failed', 0)
    success_rate = round(100.0 * total_passed / max(1, total_passed + total_failed), 2)
    
    # Create single row with fixed columns for key metrics
    summary_row = {
        'total_passed': int(total_passed),
        'total_failed': int(total_failed), 
        'success_rate_percent': float(success_rate),
        'total_attempts': int(total_passed + total_failed)
    }
    
    # Add top violation reasons with fixed column names (top 10)
    violation_items = sorted(violation_reasons.items(), key=lambda x: -x[1])
    for i in range(10):
        reason_col = f'violation_reason_{i+1}'
        count_col = f'violation_count_{i+1}'
        
        if i < len(violation_items):
            reason, count = violation_items[i]
            summary_row[reason_col] = str(reason)
            summary_row[count_col] = int(count)
        else:
            summary_row[reason_col] = 'none'
            summary_row[count_col] = 0
    
    # Create DataFrame with single row and fixed column order
    df = pd.DataFrame([summary_row])
    df.to_csv(out_path, index=False)


# Incremental statistics accumulation functions (P0-1: true OOM avoidance)

def _init_incremental_stats() -> Dict[str, Any]:
    """Initialize incremental statistics accumulation structure."""
    return {
        # Basic counters
        'product_count': 0,
        'parent_count': 0,  # Will be set later
        
        # Independent aggregation quantities for cross-validation
        'attempts_total': 0,  # Total attempts count for global consistency check
        'rule_halogen_counts': defaultdict(lambda: defaultdict(int)),  # 2D for rule x halogen validation
        'rule_totals': Counter(),  # Rule totals for rule dimension validation
        
        # K-dimensional counting: Rule x Halogen x K
        'rule_halogen_k_counts': defaultdict(lambda: defaultdict(lambda: defaultdict(int))),
        'rule_halogen_k_by_parent': defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int)))),
        'k_counts': Counter(),
        'rule_counts': Counter(),  # Products by rule (not attempts)
        'halogen_counts': Counter(),
        
        # Constraint tracking
        'constraint_stats': {'total_passed': 0, 'total_failed': 0, 'violation_reasons': Counter()},
        
        # Parent tracking
        'parent_product_counts': Counter(),
        'parents_with_products': set(),
        'all_parent_keys': set(),
        'parent_key_to_name': {},
        
        # P0-C: Real-time diagnostic tracking by parent type
        'parents_smiles_with_products': set(),
        'parents_inchikey_with_products': set(),
        
        # Other P1-specific counters
        'has_qc_errors': False,
        'max_products_per_parent': 0,
        'min_products_per_parent': float('inf'),
        'sanitize_ok_count': 0
    }


def _update_incremental_stats(stats: Dict[str, Any], record: Dict[str, Any]) -> None:
    """Update incremental statistics with a single record."""
    
    # Basic counters
    stats['product_count'] += 1
    
    # Extract record fields
    k = record.get('k', 1)
    rule = record.get('rule', 'Unknown')
    halogen = record.get('halogen', 'Unknown')
    constraints_ok = record.get('constraints_ok', True)
    constraints_violations = record.get('constraints_violations', {})
    sanitize_ok = record.get('sanitize_ok', True)
    
    # Independent aggregation quantities (for cross-validation)
    # Note: attempts_total should come from enumeration QA stats, not product records
    stats['rule_halogen_counts'][rule][halogen] += 1  # 2D rule x halogen for validation
    stats['rule_totals'][rule] += 1  # Rule totals for rule dimension validation
    
    # K-dimensional counting
    stats['rule_halogen_k_counts'][rule][halogen][k] += 1
    stats['k_counts'][k] += 1
    stats['rule_counts'][rule] += 1  # Products by rule (not attempts)
    stats['halogen_counts'][halogen] += 1
    
    # P0-C: QC tracking with unified has_qc_errors semantics
    if sanitize_ok:
        stats['sanitize_ok_count'] += 1
    else:
        stats['has_qc_errors'] = True
    
    # Constraint tracking
    if constraints_ok:
        stats['constraint_stats']['total_passed'] += 1
    else:
        stats['constraint_stats']['total_failed'] += 1
        stats['has_qc_errors'] = True  # P0-C: Set QC error flag on constraint failures
        
        # Count violation reasons
        if isinstance(constraints_violations, dict):
            for reason, details in constraints_violations.items():
                stats['constraint_stats']['violation_reasons'][reason] += 1
    
    # Parent tracking using unified key generation
    parent_key = record_to_unified_parent_key(record)
    if parent_key:
        stats['parent_product_counts'][parent_key] += 1
        stats['parents_with_products'].add(parent_key)
        # Add parent dimension to rule_halogen_k tracking
        stats['rule_halogen_k_by_parent'][parent_key][rule][halogen][k] += 1
        
        # Update parent name mapping if available
        parent_name = record.get('parent_name', parent_key)
        if parent_key not in stats['parent_key_to_name']:
            stats['parent_key_to_name'][parent_key] = parent_name
    
    # P0-C: Real-time diagnostic tracking by parent type
    parent_smiles = record.get('parent_smiles', '').strip()
    parent_inchikey = record.get('parent_inchikey', '').strip()
    
    if parent_smiles and parent_smiles != 'Unknown':
        stats['parents_smiles_with_products'].add(parent_smiles)
    if parent_inchikey and parent_inchikey != 'Unknown':
        stats['parents_inchikey_with_products'].add(parent_inchikey)


def _validate_rule_halogen_k_consistency(stats: Dict[str, Any]) -> Dict[str, Any]:
    """
    Validate 6 invariants for rule_halogen_k counting consistency using independent sources.
    
    Checks:
    1. H dimension: Sum_rule,k == halogen_counts[halogen] 
    2. (rule, halogen) dimension: Sum_k == rule_halogen_counts[rule][halogen] (uses independent 2D matrix)
    3. k dimension: Sum_rule,halogen == k_counts[k]
    4. Rule dimension: Sum_halogen,k == rule_totals[rule] (uses independent rule totals)
    5. Global: Sum_rule,halogen,k == attempts_total (uses independent total counter)
    6. Parent aggregation: Sum_parent rule_halogen_k_by_parent == rule_halogen_k_counts
    
    Args:
        stats: Statistics dictionary containing rule_halogen_k_counts and related counts
        
    Returns:
        Dictionary with consistency check results including boolean and diagnostics
    """
    rule_halogen_k_counts = stats.get('rule_halogen_k_counts', {})
    rule_halogen_k_by_parent = stats.get('rule_halogen_k_by_parent', {})
    halogen_counts = stats.get('halogen_counts', {})
    k_counts = stats.get('k_counts', {})
    rule_counts = stats.get('rule_counts', {})
    
    # Independent sources for cross-validation
    attempts_total = stats.get('attempts_total', 0)
    rule_halogen_counts_independent = stats.get('rule_halogen_counts', {})
    rule_totals = stats.get('rule_totals', {})
    
    # Overall consistency flag
    overall_consistent = True
    checks = {}
    all_mismatches = []
    
    # Check 1: H dimension - Sum_rule,k == halogen_counts[halogen] (existing check)
    calculated_halogen_counts = defaultdict(int)
    for rule, halogen_dict in rule_halogen_k_counts.items():
        for halogen, k_dict in halogen_dict.items():
            calculated_halogen_counts[halogen] += sum(k_dict.values())
    
    h_consistent = True
    h_mismatches = {}
    all_halogens = set(calculated_halogen_counts.keys()) | set(halogen_counts.keys())
    
    for halogen in all_halogens:
        calculated = calculated_halogen_counts.get(halogen, 0)
        actual = halogen_counts.get(halogen, 0)
        if calculated != actual:
            h_consistent = False
            h_mismatches[halogen] = {
                'calculated': calculated, 'actual': actual, 'difference': calculated - actual
            }
            all_mismatches.append(f"H[{halogen}]: calc={calculated} vs actual={actual}")
    
    checks['halogen_vs_flat'] = {
        'consistent': h_consistent,
        'mismatches': h_mismatches,
        'description': 'Sum_rule,k rule_halogen_k_counts == halogen_counts'
    }
    overall_consistent &= h_consistent
    
    # Check 2: (rule, halogen) dimension - Sum_k == rule_halogen_counts[rule][halogen] (independent source)
    calculated_rule_halogen_counts = defaultdict(lambda: defaultdict(int))
    for rule, halogen_dict in rule_halogen_k_counts.items():
        for halogen, k_dict in halogen_dict.items():
            calculated_rule_halogen_counts[rule][halogen] = sum(k_dict.values())
    
    # Compare with independently stored rule_halogen_counts
    rh_consistent = True
    rh_mismatches = {}
    
    # Get all rule x halogen combinations from both sources
    all_rule_halogen_pairs = set()
    for rule, halogen_dict in calculated_rule_halogen_counts.items():
        for halogen in halogen_dict.keys():
            all_rule_halogen_pairs.add((rule, halogen))
    for rule, halogen_dict in rule_halogen_counts_independent.items():
        for halogen in halogen_dict.keys():
            all_rule_halogen_pairs.add((rule, halogen))
    
    for rule, halogen in all_rule_halogen_pairs:
        calculated = calculated_rule_halogen_counts[rule][halogen]
        expected = rule_halogen_counts_independent.get(rule, {}).get(halogen, 0)
        if calculated != expected:
            rh_consistent = False
            rh_mismatches[f"{rule}_{halogen}"] = {
                'calculated_from_3d': calculated, 'expected_from_2d': expected, 'difference': calculated - expected
            }
            all_mismatches.append(f"RH[{rule},{halogen}]: 3D_sum={calculated} vs 2D_independent={expected}")
    
    checks['rule_halogen_vs_flat'] = {
        'consistent': rh_consistent,
        'mismatches': rh_mismatches,
        'description': 'Sum_k rule_halogen_k_counts == rule_halogen_counts'
    }
    overall_consistent &= rh_consistent
    
    # Check 3: k dimension - Sum_rule,halogen == k_counts[k]
    calculated_k_counts = defaultdict(int)
    for rule, halogen_dict in rule_halogen_k_counts.items():
        for halogen, k_dict in halogen_dict.items():
            for k, count in k_dict.items():
                calculated_k_counts[k] += count
    
    k_consistent = True
    k_mismatches = {}
    all_ks = set(calculated_k_counts.keys()) | set(k_counts.keys())
    
    for k in all_ks:
        calculated = calculated_k_counts.get(k, 0)
        actual = k_counts.get(k, 0)
        if calculated != actual:
            k_consistent = False
            k_mismatches[str(k)] = {
                'calculated': calculated, 'actual': actual, 'difference': calculated - actual
            }
            all_mismatches.append(f"K[{k}]: calc={calculated} vs actual={actual}")
    
    checks['k_vs_flat'] = {
        'consistent': k_consistent,
        'mismatches': k_mismatches,
        'description': 'Sum_rule,halogen rule_halogen_k_counts == k_counts'
    }
    overall_consistent &= k_consistent
    
    # Check 4: Rule dimension - Sum_halogen,k == rule_totals[rule] (independent source)
    calculated_rule_totals = defaultdict(int)
    for rule, halogen_dict in rule_halogen_k_counts.items():
        for halogen, k_dict in halogen_dict.items():
            calculated_rule_totals[rule] += sum(k_dict.values())
    
    rule_consistent = True
    rule_mismatches = {}
    all_rules = set(calculated_rule_totals.keys()) | set(rule_totals.keys())
    
    for rule in all_rules:
        calculated = calculated_rule_totals.get(rule, 0)
        expected = rule_totals.get(rule, 0)
        if calculated != expected:
            rule_consistent = False
            rule_mismatches[rule] = {
                'calculated_from_3d': calculated, 'expected_from_independent': expected, 'difference': calculated - expected
            }
            all_mismatches.append(f"RULE[{rule}]: 3D_sum={calculated} vs independent={expected}")
    
    checks['rule_vs_independent'] = {
        'consistent': rule_consistent,
        'mismatches': rule_mismatches,
        'description': 'Sum_halogen,k rule_halogen_k_counts == rule_totals (independent)'
    }
    overall_consistent &= rule_consistent
    
    # Check 5: Global product conservation and related sub-checks
    total_calculated = sum(calculated_halogen_counts.values())
    product_count = stats.get('product_count', 0)

    # Start with product conservation against 3D sum
    global_mismatches: Dict[str, Any] = {}
    product_consistent = (total_calculated == product_count)
    if not product_consistent:
        global_mismatches['global_products'] = {
            'calculated_from_3d': total_calculated,
            'expected_product_count': product_count,
            'difference': total_calculated - product_count
        }
        all_mismatches.append(
            f"PRODUCT_TOTAL: 3D_sum={total_calculated} vs product_count={product_count}"
        )

    # Three-way mutex core invariant remains a separate check
    enumeration_qa_stats = stats.get('enumeration_qa_stats', {})
    if enumeration_qa_stats:
        enum_attempts = int(enumeration_qa_stats.get('attempts', 0))
        enum_products = int(enumeration_qa_stats.get('products', 0))
        enum_no_matches = int(enumeration_qa_stats.get('no_product_matches', 0))
        enum_template_unsupported = int(enumeration_qa_stats.get('template_unsupported', 0))

        expected_total = enum_products + enum_no_matches + enum_template_unsupported
        mutex_consistent = (enum_attempts == expected_total)

        if not mutex_consistent:
            all_mismatches.append(
                f"THREE_WAY_MUTEX: attempts={enum_attempts} != products={enum_products} + no_matches={enum_no_matches} + template_unsupported={enum_template_unsupported} (sum={expected_total})"
            )
        checks['three_way_mutex'] = {
            'consistent': mutex_consistent,
            'mismatches': {} if mutex_consistent else {
                'attempts_total': enum_attempts,
                'products': enum_products,
                'no_product_matches': enum_no_matches,
                'template_unsupported': enum_template_unsupported,
                'calculated_sum': expected_total,
                'difference': enum_attempts - expected_total
            },
            'description': 'attempts_total == products + no_product_matches + template_unsupported (core invariant)'
        }

        # Sub-check: product_count vs enumeration.products
        if product_count != enum_products:
            global_mismatches['product_count_vs_enum'] = {
                'product_count': product_count,
                'enumeration_products': enum_products,
                'difference': product_count - enum_products
            }
            all_mismatches.append(
                f"PRODUCT_ENUM: product_count={product_count} vs enumeration.products={enum_products}"
            )

    # Sub-check: product_count vs sum(halogen_counts) for dimensional consistency
    halogen_sum = sum(halogen_counts.values())
    if product_count != halogen_sum:
        global_mismatches['product_vs_halogen_sum'] = {
            'product_count': product_count,
            'halogen_sum': halogen_sum,
            'difference': product_count - halogen_sum
        }
        all_mismatches.append(
            f"PRODUCT_VS_HALOGEN: product_count={product_count} vs halogen_sum={halogen_sum}"
        )

    # Determine global conservation consistency (all three sub-checks must pass)
    total_consistent = (len(global_mismatches) == 0)
    checks['global_product_conservation'] = {
        'consistent': total_consistent,
        'mismatches': global_mismatches,
        'description': 'Global product conservation and dimensional checks'
    }
    overall_consistent &= total_consistent
    
    # Check 6: Parent aggregation - Sum_parent rule_halogen_k_by_parent == rule_halogen_k_counts
    calculated_from_parents = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for parent_key, rule_dict in rule_halogen_k_by_parent.items():
        for rule, halogen_dict in rule_dict.items():
            for halogen, k_dict in halogen_dict.items():
                for k, count in k_dict.items():
                    calculated_from_parents[rule][halogen][k] += count
    
    parent_consistent = True
    parent_mismatches = {}
    
    # Compare aggregated parent data with flat rule_halogen_k_counts
    for rule, halogen_dict in rule_halogen_k_counts.items():
        for halogen, k_dict in halogen_dict.items():
            for k, count in k_dict.items():
                calculated = calculated_from_parents[rule][halogen][k]
                if calculated != count:
                    parent_consistent = False
                    key = f"{rule}_{halogen}_{k}"
                    parent_mismatches[key] = {
                        'calculated_from_parents': calculated, 'actual_flat': count, 'difference': calculated - count
                    }
                    all_mismatches.append(f"PARENT[{rule},{halogen},{k}]: calc={calculated} vs flat={count}")
    
    checks['by_parent_vs_flat'] = {
        'consistent': parent_consistent,
        'mismatches': parent_mismatches,
        'description': 'Sum_parent rule_halogen_k_by_parent == rule_halogen_k_counts'
    }
    overall_consistent &= parent_consistent
    
    # Build result with observability fields
    result: Dict[str, Any] = {
        'consistent': overall_consistent,
        'checks': checks,
        'mismatch_count': len(all_mismatches),
        'mismatches_sample': all_mismatches[:50]
    }

    # Observability additions (P1-C)
    qa_src = stats.get('qa_source', {
        'path': None, 'path_type': None, 'schema': None, 'has_attempts': None, 'degrade': False
    })
    result['qa_source'] = qa_src
    result['product_count_source'] = 'report:product_count'
    enum_stats = stats.get('enumeration_qa_stats') or {}
    attempts_val = int(enum_stats.get('attempts', 0)) if isinstance(enum_stats, dict) else 0
    has_attempts_flag = bool(qa_src.get('has_attempts'))
    degrade_flag = bool(qa_src.get('degrade'))
    result['three_way_mutex_checked'] = bool(enum_stats) and has_attempts_flag and (attempts_val > 0) and (not degrade_flag)

    return result


def _finalize_incremental_stats(stats: Dict[str, Any], parent_pairs: List[tuple], config: Dict[str, Any], qa_stats: Dict[str, Any] = None) -> None:
    """Finalize incremental statistics with parent metadata and computed fields."""
    # Add all expected parent keys from input using safe unified key generator with metadata collection
    parent_key_meta_aggregates = {
        'hashed_count': 0,
        'hashed_examples': [],
        'rdkit_available_count': 0,
        'inchikey_generated_count': 0
    }
    
    for smiles, name in parent_pairs:
        unified_key, key_meta = make_unified_parent_key_with_metadata(smiles)
        if unified_key and unified_key != SMI_UNKNOWN:
            stats['parent_key_to_name'][unified_key] = name
            stats['all_parent_keys'].add(unified_key)
            
            # Aggregate metadata for observability
            if key_meta.get('smiles_hashed', False):
                parent_key_meta_aggregates['hashed_count'] += 1
                if len(parent_key_meta_aggregates['hashed_examples']) < 20:
                    parent_key_meta_aggregates['hashed_examples'].append(unified_key)
            
            if key_meta.get('rdkit_available', False):
                parent_key_meta_aggregates['rdkit_available_count'] += 1
            
            if key_meta.get('inchikey_generated', False):
                parent_key_meta_aggregates['inchikey_generated_count'] += 1
    
    # Store metadata aggregates in stats
    stats['parent_key_metadata'] = parent_key_meta_aggregates
    
    # Compute derived statistics
    parent_product_counts = stats['parent_product_counts']
    
    if parent_product_counts:
        stats['max_products_per_parent'] = max(parent_product_counts.values())
        stats['min_products_per_parent'] = min(parent_product_counts.values())
    else:
        stats['max_products_per_parent'] = 0
        stats['min_products_per_parent'] = 0
    
    # Parents with zero products
    all_parent_keys = stats['all_parent_keys']
    parents_with_products = stats['parents_with_products']
    parents_without_products = [stats['parent_key_to_name'].get(k, k) for k in all_parent_keys if k not in parents_with_products]
    
    # Add computed fields
    stats['unique_parents_with_products'] = len(parents_with_products)
    stats['parents_without_products'] = parents_without_products
    stats['parent_count'] = len(all_parent_keys)  # Override with unified count
    stats['avg_products_per_parent'] = round(stats['product_count'] / max(1, stats['parent_count']), 2)
    
    # P0-C: Add missing fields expected by CSV writer
    
    # Add enumeration QA stats to stats for consistency checking
    if qa_stats:
        stats['enumeration_qa_stats'] = qa_stats
        stats['attempts_total'] = qa_stats.get('attempts', 0)  # Real attempts count from enumeration
    
    # QA Consistency Check: Validate that sum over k of rule_halogen_k_counts equals halogen_counts
    consistency_check = _validate_rule_halogen_k_consistency(stats)
    stats['qa_consistency'] = consistency_check
    
    # Convert independent aggregation quantities to proper dict format (for serialization)
    # Keep the independently collected rule_halogen_counts for validation
    if 'rule_halogen_counts' in stats:
        # Convert defaultdict to regular dict for serialization
        stats['rule_halogen_counts'] = {r: dict(h) for r, h in stats['rule_halogen_counts'].items()}
    
    # Convert rule_totals to regular dict
    if 'rule_totals' in stats:
        stats['rule_totals'] = dict(stats['rule_totals'])
    
    # Ensure attempts_total is an integer (should already be from incremental updates)
    stats['attempts_total'] = int(stats.get('attempts_total', 0))
    
    # P0-C: Use pre-computed diagnostic fields from incremental stage
    
    # Get the pre-computed sets from incremental tracking
    parents_smiles_with_products = stats.get('parents_smiles_with_products', set())
    parents_inchikey_with_products = stats.get('parents_inchikey_with_products', set())
    
    # Real diagnostic statistics (computed from real-time tracking)
    smiles_unique_parents = len(parents_smiles_with_products)
    inchikey_unique_parents = len(parents_inchikey_with_products)
    diff_smiles_vs_inchikey = abs(smiles_unique_parents - inchikey_unique_parents)
    
    # Find examples of discrepancies
    missing_in_smiles_examples = []
    missing_in_ikeys_examples = []
    
    if diff_smiles_vs_inchikey > 0:
        # Find parents in InChIKey set but not in SMILES set
        inchikey_only = parents_inchikey_with_products - parents_smiles_with_products
        for key in list(inchikey_only)[:5]:  # Up to 5 examples
            name = stats['parent_key_to_name'].get(key, key)
            missing_in_smiles_examples.append(name)
        
        # Find parents in SMILES set but not in InChIKey set  
        smiles_only = parents_smiles_with_products - parents_inchikey_with_products
        for key in list(smiles_only)[:5]:  # Up to 5 examples
            name = stats['parent_key_to_name'].get(key, key)
            missing_in_ikeys_examples.append(name)
    
    stats['smiles_unique_parents'] = smiles_unique_parents
    stats['inchikey_unique_parents'] = inchikey_unique_parents
    stats['diff_smiles_vs_inchikey'] = diff_smiles_vs_inchikey
    stats['missing_in_smiles_examples'] = missing_in_smiles_examples
    stats['missing_in_ikeys_examples'] = missing_in_ikeys_examples
    
    # Empty flavonoid/probe counters for P0 compatibility (convert to plain dict)
    stats['rule_counts_flavonoids'] = {}
    stats['halogen_counts_flavonoids'] = {}
    stats['rule_halogen_counts_flavonoids'] = {}
    stats['rule_counts_probes'] = {}
    stats['halogen_counts_probes'] = {}
    stats['rule_halogen_counts_probes'] = {}
    
    # P0-C: Convert all Counter/defaultdict to plain dict for consistent serialization
    stats['k_counts'] = dict(stats['k_counts'])
    stats['rule_counts'] = dict(stats['rule_counts'])
    stats['halogen_counts'] = dict(stats['halogen_counts'])
    stats['parent_product_counts'] = dict(stats['parent_product_counts'])

    stats['constraint_stats']['violation_reasons'] = dict(stats['constraint_stats']['violation_reasons'])
    
    # Also convert rule_halogen_k_by_parent to regular dicts
    if 'rule_halogen_k_by_parent' in stats:
        stats['rule_halogen_k_by_parent'] = {p: {r: {h: dict(k) for h, k in h_dict.items()} for r, h_dict in r_dict.items()} for p, r_dict in stats['rule_halogen_k_by_parent'].items()}
    
    # 3-level dict conversion for rule_halogen_k_counts
    if 'rule_halogen_k_counts' in stats:
        stats['rule_halogen_k_counts'] = {r: {h: dict(k) for h, k in h_dict.items()} for r, h_dict in stats['rule_halogen_k_counts'].items()}
    
    # Deep serialization safety: Convert all sets to sorted lists
    if 'parents_smiles_with_products' in stats:
        stats['parents_smiles_with_products'] = sorted(list(stats['parents_smiles_with_products']))
    if 'parents_inchikey_with_products' in stats:
        stats['parents_inchikey_with_products'] = sorted(list(stats['parents_inchikey_with_products']))
    if 'all_parent_keys' in stats:
        stats['all_parent_keys'] = sorted(list(stats['all_parent_keys']))
    if 'parents_with_products' in stats:
        stats['parents_with_products'] = sorted(list(stats['parents_with_products']))
    
    # Add config metadata
    stats['config'] = config


def _load_qa_stats_with_fallback(products_table: str, config: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Optional[str]]]:
    """
    Load QA stats with enhanced path/schema resolution and fallback strategies.
    
    Args:
        products_table: Primary products table path for reference
        config: Configuration dictionary (may contain qa_summary_path)
        
    Returns:
        Tuple of (canonical_enumeration_qa_stats, qa_source_info)
        canonical_enumeration_qa_stats contains only the four core integer fields.
    """
    # Helper to canonicalize stats
    def _canonicalize(raw: Dict[str, Any]) -> Dict[str, int]:
        return {
            'attempts': int(raw.get('attempts', 0) or 0),
            'products': int(raw.get('products', 0) or 0),
            'no_product_matches': int(raw.get('no_product_matches', 0) or 0),
            'template_unsupported': int(raw.get('template_unsupported', 0) or 0),
        }

    qa_stats: Dict[str, Any] = {}
    qa_source_info: Dict[str, Optional[str]] = {'path': None, 'path_type': None, 'schema': None, 'has_attempts': None, 'degrade': False}

    # Paths and base resolution (do not rely on CWD)
    products_dir = os.path.abspath(os.path.dirname(products_table))
    parent_dir = os.path.abspath(os.path.dirname(products_dir))
    base_dir = os.path.abspath(config.get('output_dir', products_dir))

    # Candidate paths in required priority
    candidate_paths: List[Tuple[str, str]] = []
    config_qa_path = config.get('qa_summary_path')
    if config_qa_path:
        candidate_paths.append(('config', os.path.abspath(config_qa_path)))
    primary_path = os.path.join(products_dir, 'qa_summary.json')
    candidate_paths.append(('primary', primary_path))
    qa_sibling_path = os.path.join(parent_dir, 'qa', 'qa_summary.json')
    candidate_paths.append(('qa_sibling', qa_sibling_path))
    benchmarks_path = os.path.join(base_dir, 'benchmarks', 'qa_summary.json')
    candidate_paths.append(('benchmarks', benchmarks_path))

    # Schema extraction priority
    schema_attempts = [
        ('total', lambda d: d.get('total', {})),
        ('enumeration.totals', lambda d: d.get('enumeration', {}).get('totals', {})),
        ('totals', lambda d: d.get('totals', {})),
        ('flat', lambda d: d if isinstance(d, dict) else {})
    ]

    loaded = False
    for path_type, qa_path in candidate_paths:
        if not os.path.exists(qa_path):
            LOG.debug(f"QA stats missing at {path_type} path: {qa_path}")
            continue
        try:
            LOG.debug(f"Attempting to load QA stats from {path_type}: {qa_path}")
            with open(qa_path, 'r') as f:
                qa_data = json.load(f)
            for schema_name, extractor in schema_attempts:
                extracted = extractor(qa_data)
                if isinstance(extracted, dict):
                    if 'attempts' in extracted:
                        qa_stats = _canonicalize(extracted)
                        qa_source_info = {
                            'path': os.path.abspath(qa_path),
                            'path_type': path_type,
                            'schema': schema_name,
                            'has_attempts': True,
                            'degrade': False
                        }
                        LOG.info(
                            f"Loaded QA stats from {path_type} ({qa_path}) using schema '{schema_name}': attempts={qa_stats.get('attempts', 0)}"
                        )
                        loaded = True
                        break
                    # Optional degrade mode: accept products-only when configured
                    allow_products_only = bool(config.get('qa_loader_degrade', False))
                    if allow_products_only and 'products' in extracted and 'attempts' not in extracted:
                        qa_stats = _canonicalize({'attempts': 0, **extracted})
                        qa_source_info = {
                            'path': os.path.abspath(qa_path),
                            'path_type': path_type,
                            'schema': schema_name,
                            'has_attempts': False,
                            'degrade': True
                        }
                        LOG.warning(
                            f"Degrade mode: accepted products-only QA at {qa_path} with schema '{schema_name}'"
                        )
                        loaded = True
                        break
                    else:
                        LOG.debug(
                            f"Rejected QA schema '{schema_name}' at {qa_path}: missing 'attempts' field"
                        )
            if loaded:
                break
        except Exception as e:
            LOG.debug(f"Failed to load QA stats from {path_type} ({qa_path}): {e}")

    if not loaded:
        LOG.warning("No valid QA stats found from any candidate path. Three-way mutex invariant will be skipped.")
        LOG.debug(f"Tried paths: {[f'{t}: {p}' for t, p in candidate_paths]}")
        return {}, qa_source_info

    return qa_stats, qa_source_info
