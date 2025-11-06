# -*- coding: ascii -*-
"""DL comparison module for halogenator vs deep learning outputs."""

import os
import json
import pandas as pd
from typing import Dict, List, Set, Tuple, Any, Optional
from collections import defaultdict
import logging

LOG = logging.getLogger(__name__)


class DLComparator:
    """Compare halogenator outputs with deep learning system outputs."""
    
    def __init__(self, ours_csv: str, dl_csv: str, output_dir: str):
        """
        Initialize DL comparator.
        
        Args:
            ours_csv: Path to our system's output CSV
            dl_csv: Path to DL system's output CSV  
            output_dir: Directory to write comparison results
        """
        self.ours_csv = ours_csv
        self.dl_csv = dl_csv
        self.output_dir = output_dir
        
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Load data
        self.ours_df = None
        self.dl_df = None
        self.matches = []
        self.only_ours = []
        self.only_dl = []
        
    def load_data(self) -> Tuple[bool, str]:
        """
        Load CSV data from both systems.
        
        Returns:
            (success: bool, error_message: str)
        """
        try:
            # Load our system's CSV
            if not os.path.exists(self.ours_csv):
                return False, f"Our system CSV not found: {self.ours_csv}"
            
            self.ours_df = pd.read_csv(self.ours_csv)
            LOG.info(f"Loaded ours CSV: {len(self.ours_df)} rows")
            
            # Load DL system's CSV
            if not os.path.exists(self.dl_csv):
                return False, f"DL system CSV not found: {self.dl_csv}"
            
            self.dl_df = pd.read_csv(self.dl_csv)
            LOG.info(f"Loaded DL CSV: {len(self.dl_df)} rows")
            
            return True, ""
            
        except Exception as e:
            return False, f"Error loading CSV files: {e}"
    
    def _normalize_smiles_key(self, smiles: str) -> str:
        """Normalize SMILES for comparison (basic normalization)."""
        if pd.isna(smiles) or not smiles:
            return ""
        return str(smiles).strip()
    
    def compare_products(self, key_column: str = 'smiles') -> Dict[str, Any]:
        """
        Compare products between the two systems using specified key column.
        
        Args:
            key_column: Column to use for comparison (default: 'smiles')
            
        Returns:
            Dictionary with comparison statistics
        """
        if self.ours_df is None or self.dl_df is None:
            raise ValueError("Data not loaded. Call load_data() first.")
        
        # Check if key column exists
        if key_column not in self.ours_df.columns:
            raise ValueError(f"Key column '{key_column}' not found in our system CSV")
        if key_column not in self.dl_df.columns:
            raise ValueError(f"Key column '{key_column}' not found in DL system CSV")
        
        # Create normalized key sets
        ours_keys = set()
        ours_key_to_row = {}
        
        for idx, row in self.ours_df.iterrows():
            key = self._normalize_smiles_key(row[key_column])
            if key:
                ours_keys.add(key)
                ours_key_to_row[key] = row.to_dict()
        
        dl_keys = set()
        dl_key_to_row = {}
        
        for idx, row in self.dl_df.iterrows():
            key = self._normalize_smiles_key(row[key_column])
            if key:
                dl_keys.add(key)
                dl_key_to_row[key] = row.to_dict()
        
        # Find matches and differences
        matches = ours_keys & dl_keys
        only_ours = ours_keys - dl_keys
        only_dl = dl_keys - ours_keys
        
        # Prepare detailed records
        self.matches = []
        for key in matches:
            match_record = {
                'key': key,
                'ours': ours_key_to_row[key],
                'dl': dl_key_to_row[key]
            }
            self.matches.append(match_record)
        
        self.only_ours = [ours_key_to_row[key] for key in only_ours]
        self.only_dl = [dl_key_to_row[key] for key in only_dl]
        
        # Generate comparison statistics
        stats = {
            'total_ours': len(ours_keys),
            'total_dl': len(dl_keys),
            'matches': len(matches),
            'only_ours': len(only_ours),
            'only_dl': len(only_dl),
            'overlap_percentage': (len(matches) / max(len(ours_keys), len(dl_keys), 1)) * 100,
            'ours_coverage': (len(matches) / max(len(ours_keys), 1)) * 100,
            'dl_coverage': (len(matches) / max(len(dl_keys), 1)) * 100
        }
        
        LOG.info(f"Comparison complete - Matches: {stats['matches']}, "
                   f"Only ours: {stats['only_ours']}, Only DL: {stats['only_dl']}")
        
        return stats
    
    def generate_cross_tabulation(self) -> Dict[str, Any]:
        """
        Generate cross-tabulation by rule, halogen, and k dimensions.
        
        Returns:
            Dictionary with cross-tabulation results
        """
        def extract_dimensions(df, prefix):
            """Extract rule/halogen/k dimensions from dataframe."""
            dims = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
            
            for _, row in df.iterrows():
                # Extract dimensions (with fallbacks)
                rule = row.get('rule', 'unknown')
                halogen = row.get('halogen', 'unknown') 
                k = row.get('k', 0)
                
                dims[rule][halogen][k] += 1
            
            return dims
        
        if self.ours_df is None or self.dl_df is None:
            return {}
        
        # Extract dimensions from both datasets
        ours_dims = extract_dimensions(self.ours_df, 'ours')
        dl_dims = extract_dimensions(self.dl_df, 'dl')
        
        # Create combined dimension space
        all_rules = set(ours_dims.keys()) | set(dl_dims.keys())
        all_halogens = set()
        all_ks = set()
        
        for dims in [ours_dims, dl_dims]:
            for rule_data in dims.values():
                all_halogens.update(rule_data.keys())
                for halogen_data in rule_data.values():
                    all_ks.update(halogen_data.keys())
        
        # Build cross-tabulation
        crosstab = []
        for rule in sorted(all_rules):
            for halogen in sorted(all_halogens):
                for k in sorted(all_ks):
                    ours_count = ours_dims[rule][halogen][k]
                    dl_count = dl_dims[rule][halogen][k]
                    
                    if ours_count > 0 or dl_count > 0:
                        crosstab.append({
                            'rule': rule,
                            'halogen': halogen,
                            'k': k,
                            'ours_count': ours_count,
                            'dl_count': dl_count,
                            'difference': ours_count - dl_count,
                            'match': ours_count == dl_count
                        })
        
        return {
            'crosstab': crosstab,
            'summary': {
                'total_combinations': len(crosstab),
                'matching_combinations': sum(1 for row in crosstab if row['match']),
                'differing_combinations': sum(1 for row in crosstab if not row['match'])
            }
        }
    
    def write_output_files(self, stats: Dict[str, Any], 
                          crosstab: Dict[str, Any]) -> Dict[str, str]:
        """
        Write the three required CSV files and summary JSON.
        
        Returns:
            Dictionary mapping file type to file path
        """
        files_written = {}
        
        # 1. Write matches.csv
        matches_file = os.path.join(self.output_dir, 'matches.csv')
        if self.matches:
            matches_df = pd.DataFrame([
                {
                    **match['ours'],
                    **{f"dl_{k}": v for k, v in match['dl'].items()}
                } for match in self.matches
            ])
            matches_df.to_csv(matches_file, index=False)
        else:
            # Write empty CSV with headers
            pd.DataFrame().to_csv(matches_file, index=False)
        files_written['matches'] = matches_file
        
        # 2. Write only_ours.csv  
        only_ours_file = os.path.join(self.output_dir, 'only_ours.csv')
        if self.only_ours:
            pd.DataFrame(self.only_ours).to_csv(only_ours_file, index=False)
        else:
            pd.DataFrame().to_csv(only_ours_file, index=False)
        files_written['only_ours'] = only_ours_file
        
        # 3. Write only_dl.csv
        only_dl_file = os.path.join(self.output_dir, 'only_dl.csv')
        if self.only_dl:
            pd.DataFrame(self.only_dl).to_csv(only_dl_file, index=False)
        else:
            pd.DataFrame().to_csv(only_dl_file, index=False)
        files_written['only_dl'] = only_dl_file
        
        # 4. Write summary.json
        summary_file = os.path.join(self.output_dir, 'summary.json')
        summary_data = {
            'comparison_stats': stats,
            'cross_tabulation': crosstab,
            'files': {
                'matches_csv': matches_file,
                'only_ours_csv': only_ours_file,
                'only_dl_csv': only_dl_file
            },
            'metadata': {
                'ours_csv_source': self.ours_csv,
                'dl_csv_source': self.dl_csv,
                'total_files_generated': 4
            }
        }
        
        with open(summary_file, 'w') as f:
            json.dump(summary_data, f, indent=2)
        files_written['summary'] = summary_file
        
        return files_written
    
    def run_comparison(self, key_column: str = 'smiles') -> Dict[str, Any]:
        """
        Run complete comparison workflow.
        
        Args:
            key_column: Column to use for comparison
            
        Returns:
            Dictionary with complete results
        """
        # Load data
        success, error = self.load_data()
        if not success:
            return {'success': False, 'error': error}
        
        # Compare products
        try:
            stats = self.compare_products(key_column)
            crosstab = self.generate_cross_tabulation()
            files_written = self.write_output_files(stats, crosstab)
            
            return {
                'success': True,
                'stats': stats,
                'crosstab': crosstab,
                'files': files_written
            }
            
        except Exception as e:
            return {'success': False, 'error': str(e)}


def compare_dl_outputs(ours_csv: str, dl_csv: str, output_dir: str, 
                      key_column: str = 'smiles') -> Dict[str, Any]:
    """
    High-level function to compare DL outputs.
    
    Args:
        ours_csv: Path to our system's CSV output
        dl_csv: Path to DL system's CSV output
        output_dir: Output directory for comparison results
        key_column: Column to use for comparison (default: 'smiles')
        
    Returns:
        Comparison results dictionary
    """
    comparator = DLComparator(ours_csv, dl_csv, output_dir)
    return comparator.run_comparison(key_column)