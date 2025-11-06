# -*- coding: ascii -*-
"""Enhanced DL comparison module with fixed schema and advanced features."""

from __future__ import annotations

import os
import json
import pandas as pd
from typing import Dict, List, Set, Tuple, Any, Optional, Union
from collections import defaultdict, Counter
import logging
from datetime import datetime
import math

# TYPE_CHECKING import for RDKit (standard approach)
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from rdkit import Chem

LOG = logging.getLogger(__name__)

# Fixed core column schemas
FIXED_CORE_COLUMNS = {
    'matches': [
        'key', 'smiles', 'inchikey', 'rule', 'halogen', 'k', 
        'parent_key_or_smiles', 'count_ours', 'count_dl'
    ],
    'only_ours': [
        'key', 'smiles', 'inchikey', 'rule', 'halogen', 'k',
        'parent_key_or_smiles', 'count'
    ],
    'only_dl': [
        'key', 'smiles', 'inchikey', 'rule', 'halogen', 'k',
        'parent_key_or_smiles', 'count'
    ]
}


class JoinKeyNormalizer:
    """Handle join key normalization with fallback strategies."""
    
    def __init__(self, strategy: str = 'inchikey'):
        """
        Initialize normalizer.
        
        Args:
            strategy: Primary strategy ('inchikey', 'canonical_smiles', 'raw_smiles')
        """
        self.strategy = strategy
        self.rdkit_available = False
        self.normalization_stats = {
            'inchikey_generated': 0,
            'canonical_smiles_used': 0,
            'raw_smiles_used': 0,
            'failures': 0,
            'normalization_path': strategy
        }
        
        # Test RDKit availability
        try:
            from rdkit import Chem
            self.rdkit_available = True
            LOG.debug("RDKit available for join key normalization")
        except ImportError:
            LOG.warning("RDKit not available - using SMILES-only normalization")
    
    def normalize_key(self, smiles: str, inchikey: Optional[str] = None) -> Tuple[str, Dict[str, Any]]:
        """
        Normalize a key using the configured strategy with fallbacks.
        
        Args:
            smiles: SMILES string
            inchikey: Optional InChIKey if available
            
        Returns:
            Tuple of (normalized_key, metadata)
        """
        metadata = {'method_used': 'unknown', 'original_smiles': smiles}
        
        # Clean inputs
        smiles = str(smiles or '').strip()
        inchikey = str(inchikey or '').strip() if inchikey else ''
        
        # Handle empty/invalid inputs
        if not smiles or smiles.lower() == 'unknown':
            metadata['method_used'] = 'empty'
            return '', metadata
        
        # Strategy 1: Use existing InChIKey if valid
        if self.strategy == 'inchikey' and inchikey and inchikey.lower() != 'unknown':
            metadata['method_used'] = 'existing_inchikey'
            return inchikey.upper(), metadata
        
        # Strategy 2: Generate InChIKey from SMILES (if RDKit available)
        if self.strategy == 'inchikey' and self.rdkit_available:
            try:
                from rdkit import Chem
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    generated_inchikey = Chem.MolToInchiKey(mol)
                    if generated_inchikey:
                        self.normalization_stats['inchikey_generated'] += 1
                        metadata['method_used'] = 'generated_inchikey'
                        return generated_inchikey.upper(), metadata
            except Exception as e:
                LOG.debug(f"InChIKey generation failed for {smiles}: {e}")
        
        # Strategy 3: Use canonical SMILES (if RDKit available)
        if self.rdkit_available:
            try:
                from rdkit import Chem
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
                    if canonical_smiles:
                        self.normalization_stats['canonical_smiles_used'] += 1
                        metadata['method_used'] = 'canonical_smiles'
                        return canonical_smiles, metadata
            except Exception as e:
                LOG.debug(f"Canonical SMILES generation failed for {smiles}: {e}")
        
        # Strategy 4: Raw SMILES (fallback)
        self.normalization_stats['raw_smiles_used'] += 1
        metadata['method_used'] = 'raw_smiles'
        return smiles, metadata
    
    def get_stats(self) -> Dict[str, Any]:
        """Get normalization statistics."""
        total_processed = sum([
            self.normalization_stats['inchikey_generated'],
            self.normalization_stats['canonical_smiles_used'], 
            self.normalization_stats['raw_smiles_used'],
            self.normalization_stats['failures']
        ])
        
        return {
            **self.normalization_stats,
            'total_processed': total_processed,
            'rdkit_available': self.rdkit_available,
            'failure_rate': self.normalization_stats['failures'] / max(1, total_processed)
        }


class FixedSchemaDLComparator:
    """Enhanced DL comparator with fixed schema output."""
    
    def __init__(self, ours_csv: str, dl_csv: str, output_dir: str, 
                 key_strategy: str = 'inchikey', low_mem: bool = False):
        """
        Initialize enhanced DL comparator.
        
        Args:
            ours_csv: Path to our system's output CSV
            dl_csv: Path to DL system's output CSV
            output_dir: Directory to write comparison results
            key_strategy: Join key strategy ('inchikey', 'canonical_smiles', 'raw_smiles')
            low_mem: Enable low memory mode for large datasets
        """
        self.ours_csv = ours_csv
        self.dl_csv = dl_csv
        self.output_dir = output_dir
        self.key_strategy = key_strategy
        self.low_mem = low_mem
        
        # Initialize normalizer
        self.normalizer = JoinKeyNormalizer(key_strategy)
        
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Data storage
        self.ours_data = {}
        self.dl_data = {}
        self.validation_report = {
            'ours_validation': {'valid': True, 'errors': [], 'warnings': []},
            'dl_validation': {'valid': True, 'errors': [], 'warnings': []},
            'column_usage': {'ours': {}, 'dl': {}}
        }
    
    def _validate_and_load_csv(self, csv_path: str, source_name: str) -> Tuple[bool, pd.DataFrame, Dict[str, Any]]:
        """
        Validate and load CSV with strict column validation.
        
        Args:
            csv_path: Path to CSV file
            source_name: Name for error reporting ('ours' or 'dl')
            
        Returns:
            Tuple of (success, dataframe, validation_report)
        """
        validation = {'valid': True, 'errors': [], 'warnings': [], 'columns_found': []}
        
        try:
            if not os.path.exists(csv_path):
                validation['valid'] = False
                validation['errors'].append(f"CSV file not found: {csv_path}")
                return False, pd.DataFrame(), validation
            
            df = pd.read_csv(csv_path)
            validation['columns_found'] = list(df.columns)
            LOG.info(f"Loaded {source_name} CSV: {len(df)} rows, {len(df.columns)} columns")
            
            # Check for required columns
            required_columns = ['smiles']
            recommended_columns = ['rule', 'halogen', 'k', 'parent_smiles', 'parent_inchikey', 'inchikey']
            
            missing_required = [col for col in required_columns if col not in df.columns]
            missing_recommended = [col for col in recommended_columns if col not in df.columns]
            
            if missing_required:
                validation['valid'] = False
                validation['errors'].append(f"Missing required columns: {missing_required}")
            
            if missing_recommended:
                validation['warnings'].append(f"Missing recommended columns: {missing_recommended}")
            
            # Check for completely empty required columns
            for col in ['smiles']:
                if col in df.columns and df[col].isna().all():
                    validation['valid'] = False
                    validation['errors'].append(f"Required column '{col}' is completely empty")
            
            return validation['valid'], df, validation
            
        except Exception as e:
            validation['valid'] = False
            validation['errors'].append(f"Error reading CSV: {str(e)}")
            return False, pd.DataFrame(), validation
    
    def _process_dataframe_for_dedup(self, df: pd.DataFrame, source_name: str) -> Dict[str, Any]:
        """
        Process dataframe for deduplication and counting.
        
        Args:
            df: Input dataframe
            source_name: Source name for logging
            
        Returns:
            Dictionary mapping normalized keys to aggregated records
        """
        if self.low_mem:
            return self._process_dataframe_low_mem(df, source_name)
        else:
            return self._process_dataframe_standard(df, source_name)
    
    def _process_dataframe_standard(self, df: pd.DataFrame, source_name: str) -> Dict[str, Any]:
        """Standard memory processing."""
        result = {}
        key_counter = Counter()
        
        for idx, row in df.iterrows():
            # Extract key information
            smiles = str(row.get('smiles', ''))
            inchikey = str(row.get('inchikey', '')) if 'inchikey' in row else None
            
            # Normalize key
            normalized_key, key_metadata = self.normalizer.normalize_key(smiles, inchikey)
            
            if not normalized_key:
                continue  # Skip empty keys
            
            key_counter[normalized_key] += 1
            
            # If this is the first occurrence of this key, store the record
            if normalized_key not in result:
                # Build standardized record
                record = {
                    'key': normalized_key,
                    'smiles': smiles,
                    'inchikey': row.get('inchikey', ''),
                    'rule': row.get('rule', ''),
                    'halogen': row.get('halogen', ''),
                    'k': row.get('k', 0),
                    'parent_key_or_smiles': row.get('parent_inchikey', row.get('parent_smiles', '')),
                    'count': 0,  # Will be set below
                    'key_metadata': key_metadata
                }
                
                # Add source-specific columns with prefix
                prefix = f'{source_name}_'
                for col, value in row.items():
                    if col not in ['smiles', 'inchikey', 'rule', 'halogen', 'k']:  # Skip core columns
                        record[f'{prefix}{col}'] = value
                
                result[normalized_key] = record
        
        # Set final counts
        for key, record in result.items():
            record['count'] = key_counter[key]
        
        LOG.info(f"Processed {source_name}: {len(df)} rows -> {len(result)} unique keys")
        return result
    
    def _process_dataframe_low_mem(self, df: pd.DataFrame, source_name: str) -> Dict[str, Any]:
        """Low memory processing using chunking."""
        # For low memory mode, process in chunks
        chunk_size = min(10000, max(1000, len(df) // 10))
        result = {}
        key_counter = Counter()
        
        for chunk_start in range(0, len(df), chunk_size):
            chunk_end = min(chunk_start + chunk_size, len(df))
            chunk_df = df.iloc[chunk_start:chunk_end]
            
            chunk_result = self._process_dataframe_standard(chunk_df, source_name)
            
            # Merge chunk results
            for key, record in chunk_result.items():
                if key in result:
                    # Merge counts
                    result[key]['count'] += record['count']
                else:
                    result[key] = record
        
        LOG.info(f"Low-mem processed {source_name}: {len(df)} rows -> {len(result)} unique keys")
        return result
    
    def load_and_process_data(self) -> Tuple[bool, str]:
        """
        Load and process both CSV files with validation.
        
        Returns:
            (success: bool, error_message: str)
        """
        try:
            # Load and validate ours CSV
            success, ours_df, ours_validation = self._validate_and_load_csv(self.ours_csv, 'ours')
            self.validation_report['ours_validation'] = ours_validation
            
            if not success:
                return False, f"Our CSV validation failed: {'; '.join(ours_validation['errors'])}"
            
            # Load and validate DL CSV
            success, dl_df, dl_validation = self._validate_and_load_csv(self.dl_csv, 'dl')
            self.validation_report['dl_validation'] = dl_validation
            
            if not success:
                return False, f"DL CSV validation failed: {'; '.join(dl_validation['errors'])}"
            
            # Process dataframes for deduplication
            self.ours_data = self._process_dataframe_for_dedup(ours_df, 'ours')
            self.dl_data = self._process_dataframe_for_dedup(dl_df, 'dl')
            
            # Record column usage
            self.validation_report['column_usage'] = {
                'ours': ours_validation['columns_found'],
                'dl': dl_validation['columns_found']
            }
            
            return True, ""
            
        except Exception as e:
            LOG.error(f"Error in load_and_process_data: {e}")
            return False, str(e)
    
    def compute_comparisons(self) -> Dict[str, Any]:
        """
        Compute comparisons between datasets.
        
        Returns:
            Dictionary with comparison results
        """
        if not self.ours_data or not self.dl_data:
            raise ValueError("Data not loaded. Call load_and_process_data() first.")
        
        ours_keys = set(self.ours_data.keys())
        dl_keys = set(self.dl_data.keys())
        
        # Find matches and differences
        matches = ours_keys & dl_keys
        only_ours = ours_keys - dl_keys
        only_dl = dl_keys - ours_keys
        
        # Prepare structured results
        matches_records = []
        for key in matches:
            ours_record = self.ours_data[key]
            dl_record = self.dl_data[key]
            
            # Create match record with fixed columns
            match_record = {
                'key': key,
                'smiles': ours_record['smiles'],  # Use ours as primary
                'inchikey': ours_record['inchikey'] or dl_record['inchikey'],  # Best available
                'rule': ours_record['rule'] or dl_record['rule'],
                'halogen': ours_record['halogen'] or dl_record['halogen'],
                'k': ours_record['k'] or dl_record['k'],
                'parent_key_or_smiles': ours_record['parent_key_or_smiles'] or dl_record['parent_key_or_smiles'],
                'count_ours': ours_record['count'],
                'count_dl': dl_record['count']
            }
            
            # Add extended columns (ours_* and dl_* prefixed)
            for col, value in ours_record.items():
                if col.startswith('ours_'):
                    match_record[col] = value
            for col, value in dl_record.items():
                if col.startswith('dl_'):
                    match_record[col] = value
            
            matches_records.append(match_record)
        
        only_ours_records = []
        for key in only_ours:
            record = self.ours_data[key].copy()
            # Rename count to match schema
            record['count'] = record.pop('count', 0)
            only_ours_records.append(record)
        
        only_dl_records = []
        for key in only_dl:
            record = self.dl_data[key].copy()
            record['count'] = record.pop('count', 0)
            only_dl_records.append(record)
        
        # Compute summary statistics
        total_ours_count = sum(r['count'] for r in self.ours_data.values())
        total_dl_count = sum(r['count'] for r in self.dl_data.values())
        matched_ours_count = sum(r['count_ours'] for r in matches_records)
        matched_dl_count = sum(r['count_dl'] for r in matches_records)
        
        return {
            'matches': matches_records,
            'only_ours': only_ours_records,
            'only_dl': only_dl_records,
            'summary_stats': {
                'total_ours_unique_keys': len(ours_keys),
                'total_dl_unique_keys': len(dl_keys),
                'total_ours_count': total_ours_count,
                'total_dl_count': total_dl_count,
                'matches_unique_keys': len(matches),
                'matches_ours_count': matched_ours_count,
                'matches_dl_count': matched_dl_count,
                'only_ours_unique_keys': len(only_ours),
                'only_dl_unique_keys': len(only_dl),
                'overlap_percentage': (len(matches) / max(len(ours_keys | dl_keys), 1)) * 100,
                'ours_coverage': (len(matches) / max(len(ours_keys), 1)) * 100,
                'dl_coverage': (len(matches) / max(len(dl_keys), 1)) * 100
            }
        }
    
    def generate_crosstab_analysis(self, matches_records: List[Dict], only_ours_records: List[Dict], 
                                 only_dl_records: List[Dict]) -> Dict[str, Any]:
        """
        Generate dimensional cross-tabulation analysis.
        
        Args:
            matches_records: Match records
            only_ours_records: Only ours records
            only_dl_records: Only DL records
            
        Returns:
            Dictionary with crosstab analysis
        """
        # Collect all dimensional combinations
        dimension_stats = defaultdict(lambda: {
            'ours_count': 0, 'dl_count': 0, 'intersect_count': 0,
            'only_ours_count': 0, 'only_dl_count': 0
        })
        
        # Process matches (intersection)
        for record in matches_records:
            key = (record.get('rule', ''), record.get('halogen', ''), record.get('k', 0))
            dimension_stats[key]['intersect_count'] += 1
            dimension_stats[key]['ours_count'] += record.get('count_ours', 0)
            dimension_stats[key]['dl_count'] += record.get('count_dl', 0)
        
        # Process only_ours
        for record in only_ours_records:
            key = (record.get('rule', ''), record.get('halogen', ''), record.get('k', 0))
            dimension_stats[key]['only_ours_count'] += record.get('count', 0)
            dimension_stats[key]['ours_count'] += record.get('count', 0)
        
        # Process only_dl
        for record in only_dl_records:
            key = (record.get('rule', ''), record.get('halogen', ''), record.get('k', 0))
            dimension_stats[key]['only_dl_count'] += record.get('count', 0)
            dimension_stats[key]['dl_count'] += record.get('count', 0)
        
        # Convert to list format for CSV output
        crosstab_records = []
        for (rule, halogen, k), stats in dimension_stats.items():
            match_flag = stats['only_ours_count'] == 0 and stats['only_dl_count'] == 0
            crosstab_records.append({
                'rule': rule,
                'halogen': halogen, 
                'k': k,
                'ours_count': stats['ours_count'],
                'dl_count': stats['dl_count'],
                'intersect_count': stats['intersect_count'],
                'only_ours_count': stats['only_ours_count'],
                'only_dl_count': stats['only_dl_count'],
                'match_flag': match_flag
            })
        
        # Sort by rule, halogen, k
        crosstab_records.sort(key=lambda x: (x['rule'], x['halogen'], x['k']))
        
        return {
            'crosstab_records': crosstab_records,
            'summary': {
                'total_combinations': len(crosstab_records),
                'perfect_matches': sum(1 for r in crosstab_records if r['match_flag']),
                'partial_matches': sum(1 for r in crosstab_records if r['intersect_count'] > 0 and not r['match_flag']),
                'no_overlap': sum(1 for r in crosstab_records if r['intersect_count'] == 0)
            }
        }
    
    def _write_fixed_schema_csv(self, records: List[Dict], output_path: str, 
                               table_type: str) -> Dict[str, Any]:
        """
        Write CSV with fixed schema columns.
        
        Args:
            records: Records to write
            output_path: Output file path
            table_type: Type of table ('matches', 'only_ours', 'only_dl')
            
        Returns:
            Metadata about written file
        """
        if table_type not in FIXED_CORE_COLUMNS:
            raise ValueError(f"Unknown table type: {table_type}")
        
        core_columns = FIXED_CORE_COLUMNS[table_type]
        
        if not records:
            # Write empty CSV with fixed headers
            empty_df = pd.DataFrame(columns=core_columns)
            empty_df.to_csv(output_path, index=False)
            return {
                'file_path': output_path,
                'core_columns': len(core_columns),
                'total_columns': len(core_columns),
                'rows_written': 0,
                'schema_compliance': True
            }
        
        # Extract core columns data
        core_data = []
        extension_columns = set()
        
        for record in records:
            # Core columns (fixed order)
            core_row = {}
            for col in core_columns:
                core_row[col] = record.get(col, '')
            
            # Extension columns (variable, but prefixed)
            for key, value in record.items():
                if key not in core_columns:
                    extension_columns.add(key)
                    core_row[key] = value
            
            core_data.append(core_row)
        
        # Create DataFrame with fixed column order
        extension_columns = sorted(extension_columns)  # Consistent ordering
        all_columns = core_columns + extension_columns
        
        df = pd.DataFrame(core_data, columns=all_columns)
        df.to_csv(output_path, index=False)
        
        return {
            'file_path': output_path,
            'core_columns': len(core_columns),
            'extension_columns': len(extension_columns),
            'total_columns': len(all_columns),
            'rows_written': len(records),
            'schema_compliance': True,
            'column_order': all_columns
        }
    
    def write_output_files(self, comparison_results: Dict[str, Any]) -> Dict[str, Any]:
        """
        Write all output files with fixed schemas.
        
        Args:
            comparison_results: Results from compute_comparisons()
            
        Returns:
            Dictionary with file metadata
        """
        files_metadata = {}
        
        # 1. Write matches.csv with fixed schema
        matches_path = os.path.join(self.output_dir, 'matches.csv')
        matches_metadata = self._write_fixed_schema_csv(
            comparison_results['matches'], matches_path, 'matches'
        )
        files_metadata['matches'] = matches_metadata
        
        # 2. Write only_ours.csv with fixed schema
        only_ours_path = os.path.join(self.output_dir, 'only_ours.csv')
        only_ours_metadata = self._write_fixed_schema_csv(
            comparison_results['only_ours'], only_ours_path, 'only_ours'
        )
        files_metadata['only_ours'] = only_ours_metadata
        
        # 3. Write only_dl.csv with fixed schema
        only_dl_path = os.path.join(self.output_dir, 'only_dl.csv')
        only_dl_metadata = self._write_fixed_schema_csv(
            comparison_results['only_dl'], only_dl_path, 'only_dl'
        )
        files_metadata['only_dl'] = only_dl_metadata
        
        # 4. Generate and write crosstab analysis
        crosstab_results = self.generate_crosstab_analysis(
            comparison_results['matches'],
            comparison_results['only_ours'],
            comparison_results['only_dl']
        )
        
        # 5. Write crosstab_rule_halogen_k.csv
        crosstab_path = os.path.join(self.output_dir, 'crosstab_rule_halogen_k.csv')
        if crosstab_results['crosstab_records']:
            crosstab_df = pd.DataFrame(crosstab_results['crosstab_records'])
            crosstab_df.to_csv(crosstab_path, index=False)
        else:
            # Empty crosstab with headers
            empty_crosstab = pd.DataFrame(columns=[
                'rule', 'halogen', 'k', 'ours_count', 'dl_count', 'intersect_count',
                'only_ours_count', 'only_dl_count', 'match_flag'
            ])
            empty_crosstab.to_csv(crosstab_path, index=False)
        
        files_metadata['crosstab'] = {
            'file_path': crosstab_path,
            'rows_written': len(crosstab_results['crosstab_records'])
        }
        
        # 6. Write comprehensive summary.json
        summary_path = os.path.join(self.output_dir, 'summary.json')
        summary_data = {
            'comparison_stats': comparison_results['summary_stats'],
            'crosstab_analysis': crosstab_results['summary'],
            'normalization': self.normalizer.get_stats(),
            'validation_report': self.validation_report,
            'configuration': {
                'key_strategy': self.key_strategy,
                'low_mem_mode': self.low_mem,
                'join_key_used': self.key_strategy,
                'output_schema_version': '2.0'
            },
            'files_generated': {
                'matches_csv': matches_path,
                'only_ours_csv': only_ours_path,
                'only_dl_csv': only_dl_path,
                'crosstab_csv': crosstab_path
            },
            'files_metadata': files_metadata,
            'timestamp': datetime.now().isoformat()
        }
        
        with open(summary_path, 'w') as f:
            json.dump(summary_data, f, indent=2)
        
        files_metadata['summary'] = {'file_path': summary_path}
        
        return files_metadata
    
    def run_enhanced_comparison(self) -> Dict[str, Any]:
        """
        Run complete enhanced comparison workflow.
        
        Returns:
            Dictionary with comprehensive results
        """
        try:
            # Load and validate data
            success, error = self.load_and_process_data()
            if not success:
                return {'success': False, 'error': error}
            
            # Compute comparisons
            comparison_results = self.compute_comparisons()
            
            # Write output files
            files_metadata = self.write_output_files(comparison_results)
            
            LOG.info(f"Enhanced comparison completed successfully. Files written to: {self.output_dir}")
            
            return {
                'success': True,
                'summary_stats': comparison_results['summary_stats'],
                'files_metadata': files_metadata,
                'normalization_stats': self.normalizer.get_stats(),
                'validation_report': self.validation_report,
                'output_directory': self.output_dir
            }
            
        except Exception as e:
            LOG.error(f"Enhanced comparison failed: {e}")
            return {'success': False, 'error': str(e)}


def compare_dl_outputs_enhanced(ours_csv: str, dl_csv: str, output_dir: str,
                               key_strategy: str = 'inchikey', low_mem: bool = False) -> Dict[str, Any]:
    """
    Enhanced high-level function to compare DL outputs with fixed schema.
    
    Args:
        ours_csv: Path to our system's CSV output
        dl_csv: Path to DL system's CSV output
        output_dir: Output directory for comparison results
        key_strategy: Join key normalization strategy ('inchikey', 'canonical_smiles', 'raw_smiles')
        low_mem: Enable low memory mode for large datasets
        
    Returns:
        Enhanced comparison results dictionary
    """
    comparator = FixedSchemaDLComparator(ours_csv, dl_csv, output_dir, key_strategy, low_mem)
    return comparator.run_enhanced_comparison()