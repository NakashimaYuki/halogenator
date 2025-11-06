# -*- coding: ascii -*-
"""CSV comparison pipeline for halogenator outputs."""

import os
import pandas as pd
import numpy as np
from typing import Dict, Any, List, Tuple, Optional
from collections import defaultdict
import json
import logging
from datetime import datetime

# Module-level logger
LOG = logging.getLogger(__name__)


class CSVComparison:
    """Handles comparison between two CSV files with detailed analysis."""
    
    def __init__(self, reference_path: str, comparison_path: str, key_columns: List[str] = None):
        """
        Initialize CSV comparison.
        
        Args:
            reference_path: Path to reference CSV file
            comparison_path: Path to comparison CSV file  
            key_columns: List of columns to use as keys for matching rows
        """
        self.reference_path = reference_path
        self.comparison_path = comparison_path
        self.key_columns = key_columns or []
        
        self.reference_df = None
        self.comparison_df = None
        self.comparison_results = {}
        
    def load_data(self):
        """Load CSV data from both files."""
        try:
            self.reference_df = pd.read_csv(self.reference_path)
            LOG.info(f"Loaded reference CSV: {len(self.reference_df)} rows, {len(self.reference_df.columns)} columns")
        except Exception as e:
            LOG.error(f"Failed to load reference CSV {self.reference_path}: {e}")
            raise
            
        try:
            self.comparison_df = pd.read_csv(self.comparison_path)
            LOG.info(f"Loaded comparison CSV: {len(self.comparison_df)} rows, {len(self.comparison_df.columns)} columns")
        except Exception as e:
            LOG.error(f"Failed to load comparison CSV {self.comparison_path}: {e}")
            raise
    
    def compare_schemas(self) -> Dict[str, Any]:
        """Compare schemas (columns) between the two CSV files."""
        ref_cols = set(self.reference_df.columns)
        comp_cols = set(self.comparison_df.columns)
        
        schema_results = {
            'reference_columns': list(ref_cols),
            'comparison_columns': list(comp_cols),
            'common_columns': list(ref_cols & comp_cols),
            'missing_in_comparison': list(ref_cols - comp_cols),
            'extra_in_comparison': list(comp_cols - ref_cols),
            'schema_match': ref_cols == comp_cols
        }
        
        return schema_results
    
    def compare_row_counts(self) -> Dict[str, Any]:
        """Compare row counts and basic statistics."""
        ref_count = len(self.reference_df)
        comp_count = len(self.comparison_df)
        
        count_results = {
            'reference_count': ref_count,
            'comparison_count': comp_count,
            'count_difference': comp_count - ref_count,
            'count_ratio': comp_count / ref_count if ref_count > 0 else float('inf'),
            'counts_match': ref_count == comp_count
        }
        
        return count_results
    
    def compare_numerical_columns(self, common_columns: List[str]) -> Dict[str, Any]:
        """Compare numerical columns with statistical analysis."""
        numerical_results = {}
        
        for col in common_columns:
            if col in self.reference_df.columns and col in self.comparison_df.columns:
                ref_series = self.reference_df[col]
                comp_series = self.comparison_df[col]
                
                # Check if column is numerical
                if pd.api.types.is_numeric_dtype(ref_series) and pd.api.types.is_numeric_dtype(comp_series):
                    col_result = {
                        'column': col,
                        'reference_stats': {
                            'mean': float(ref_series.mean()) if not ref_series.empty else None,
                            'std': float(ref_series.std()) if not ref_series.empty else None,
                            'min': float(ref_series.min()) if not ref_series.empty else None,
                            'max': float(ref_series.max()) if not ref_series.empty else None,
                            'count': int(ref_series.count())
                        },
                        'comparison_stats': {
                            'mean': float(comp_series.mean()) if not comp_series.empty else None,
                            'std': float(comp_series.std()) if not comp_series.empty else None,
                            'min': float(comp_series.min()) if not comp_series.empty else None,
                            'max': float(comp_series.max()) if not comp_series.empty else None,
                            'count': int(comp_series.count())
                        }
                    }
                    
                    # Calculate differences
                    if col_result['reference_stats']['mean'] is not None and col_result['comparison_stats']['mean'] is not None:
                        col_result['mean_difference'] = col_result['comparison_stats']['mean'] - col_result['reference_stats']['mean']
                        if col_result['reference_stats']['mean'] != 0:
                            col_result['mean_relative_change'] = col_result['mean_difference'] / col_result['reference_stats']['mean']
                    
                    numerical_results[col] = col_result
        
        return numerical_results
    
    def compare_categorical_columns(self, common_columns: List[str]) -> Dict[str, Any]:
        """Compare categorical columns with value distribution analysis."""
        categorical_results = {}
        
        for col in common_columns:
            if col in self.reference_df.columns and col in self.comparison_df.columns:
                ref_series = self.reference_df[col]
                comp_series = self.comparison_df[col]
                
                # Check if column is categorical/string
                if not pd.api.types.is_numeric_dtype(ref_series):
                    ref_counts = ref_series.value_counts().to_dict()
                    comp_counts = comp_series.value_counts().to_dict()
                    
                    all_values = set(ref_counts.keys()) | set(comp_counts.keys())
                    
                    col_result = {
                        'column': col,
                        'reference_unique_count': len(ref_counts),
                        'comparison_unique_count': len(comp_counts),
                        'common_values': list(set(ref_counts.keys()) & set(comp_counts.keys())),
                        'missing_in_comparison': list(set(ref_counts.keys()) - set(comp_counts.keys())),
                        'extra_in_comparison': list(set(comp_counts.keys()) - set(ref_counts.keys())),
                        'value_distributions': {}
                    }
                    
                    # Compare distributions
                    for value in all_values:
                        ref_count = ref_counts.get(value, 0)
                        comp_count = comp_counts.get(value, 0)
                        col_result['value_distributions'][str(value)] = {
                            'reference_count': ref_count,
                            'comparison_count': comp_count,
                            'difference': comp_count - ref_count
                        }
                    
                    categorical_results[col] = col_result
        
        return categorical_results
    
    def run_comparison(self) -> Dict[str, Any]:
        """Run complete comparison analysis."""
        if self.reference_df is None or self.comparison_df is None:
            self.load_data()
        
        LOG.info("Starting CSV comparison analysis...")
        
        # Schema comparison
        schema_results = self.compare_schemas()
        common_columns = schema_results['common_columns']
        
        # Row count comparison
        count_results = self.compare_row_counts()
        
        # Numerical columns comparison
        numerical_results = self.compare_numerical_columns(common_columns)
        
        # Categorical columns comparison
        categorical_results = self.compare_categorical_columns(common_columns)
        
        # Overall assessment
        overall_assessment = {
            'schemas_match': schema_results['schema_match'],
            'counts_match': count_results['counts_match'],
            'has_numerical_differences': any(
                abs(result.get('mean_difference', 0)) > 1e-10 
                for result in numerical_results.values()
            ),
            'has_categorical_differences': any(
                len(result['missing_in_comparison']) > 0 or len(result['extra_in_comparison']) > 0
                for result in categorical_results.values()
            )
        }
        
        overall_assessment['files_match'] = (
            overall_assessment['schemas_match'] and 
            overall_assessment['counts_match'] and 
            not overall_assessment['has_numerical_differences'] and
            not overall_assessment['has_categorical_differences']
        )
        
        self.comparison_results = {
            'timestamp': datetime.now().isoformat(),
            'reference_file': self.reference_path,
            'comparison_file': self.comparison_path,
            'schema_comparison': schema_results,
            'count_comparison': count_results,
            'numerical_comparison': numerical_results,
            'categorical_comparison': categorical_results,
            'overall_assessment': overall_assessment
        }
        
        LOG.info(f"Comparison complete. Files match: {overall_assessment['files_match']}")
        return self.comparison_results
    
    def save_report(self, output_path: str):
        """Save comparison report to JSON file."""
        if not self.comparison_results:
            self.run_comparison()
        
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(self.comparison_results, f, indent=2)
        
        LOG.info(f"Comparison report saved to {output_path}")
    
    def save_csv_reports(self, output_dir: str, base_name: str = "comparison"):
        """Save comparison results as 3 CSV files: matches, missing, extra."""
        if self.reference_df is None or self.comparison_df is None:
            self.load_data()
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Use key columns for matching if specified, otherwise use all columns
        if self.key_columns:
            match_columns = self.key_columns
        else:
            # Use common columns between both dataframes
            match_columns = list(set(self.reference_df.columns) & set(self.comparison_df.columns))
            if not match_columns:
                LOG.warning("No common columns found between datasets, using all columns from reference")
                match_columns = list(self.reference_df.columns)
        
        # Create comparison keys
        ref_keys = self.reference_df[match_columns].apply(lambda x: '|'.join(x.astype(str)), axis=1)
        comp_keys = self.comparison_df[match_columns].apply(lambda x: '|'.join(x.astype(str)), axis=1)
        
        ref_key_set = set(ref_keys)
        comp_key_set = set(comp_keys)
        
        # Find matches, missing, and extra records
        matching_keys = ref_key_set & comp_key_set
        missing_keys = ref_key_set - comp_key_set  # In reference but not in comparison
        extra_keys = comp_key_set - ref_key_set    # In comparison but not in reference
        
        # Extract matching records from both datasets
        ref_matches = self.reference_df[ref_keys.isin(matching_keys)].copy()
        comp_matches = self.comparison_df[comp_keys.isin(matching_keys)].copy()
        
        # Add source labels
        ref_matches['source'] = 'reference'
        comp_matches['source'] = 'comparison'
        
        # Combine matches (reference + comparison records for same keys)
        matches_df = pd.concat([ref_matches, comp_matches], ignore_index=True, sort=False)
        
        # Extract missing and extra records
        missing_df = self.reference_df[ref_keys.isin(missing_keys)].copy()
        extra_df = self.comparison_df[comp_keys.isin(extra_keys)].copy()
        
        # Save CSV files
        matches_path = os.path.join(output_dir, f"{base_name}_matches.csv")
        missing_path = os.path.join(output_dir, f"{base_name}_missing.csv") 
        extra_path = os.path.join(output_dir, f"{base_name}_extra.csv")
        
        matches_df.to_csv(matches_path, index=False)
        missing_df.to_csv(missing_path, index=False)
        extra_df.to_csv(extra_path, index=False)
        
        # Log results
        LOG.info(f"CSV comparison results saved:")
        LOG.info(f"  Matches ({len(matches_df)} records): {matches_path}")
        LOG.info(f"  Missing ({len(missing_df)} records): {missing_path}")
        LOG.info(f"  Extra ({len(extra_df)} records): {extra_path}")
        
        return {
            'matches_file': matches_path,
            'missing_file': missing_path,
            'extra_file': extra_path,
            'matches_count': len(matches_df),
            'missing_count': len(missing_df),
            'extra_count': len(extra_df)
        }
    
    def save_business_csv_reports(self, output_dir: str, base_name: str = "business"):
        """Save business-specific comparison results as 3 CSV files for DL analysis."""
        if self.reference_df is None or self.comparison_df is None:
            self.load_data()
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Expected columns for halogenation comparison
        required_cols = ['parent_inchikey', 'rule', 'halogen', 'k', 'product_smiles']
        
        # Check if datasets have the required structure
        ref_cols = set(self.reference_df.columns)
        comp_cols = set(self.comparison_df.columns)
        
        if not all(col in ref_cols for col in required_cols):
            LOG.warning(f"Reference dataset missing required columns: {set(required_cols) - ref_cols}")
            return None
            
        if not all(col in comp_cols for col in required_cols):
            LOG.warning(f"Comparison dataset missing required columns: {set(required_cols) - comp_cols}")
            return None
        
        # Create composite keys for matching with stronger primary key design
        def create_key(df):
            # Use site_id if available, otherwise fallback to product_smiles for uniqueness
            if 'site_id' in df.columns:
                key_cols = ['parent_inchikey', 'rule', 'halogen', 'k', 'site_id'] 
                return df[key_cols].astype(str).agg('|'.join, axis=1)
            else:
                key_cols = ['parent_inchikey', 'rule', 'halogen', 'k', 'product_smiles']
                return df[key_cols].astype(str).agg('|'.join, axis=1)
        
        ref_keys = create_key(self.reference_df)
        comp_keys = create_key(self.comparison_df)
        
        # Deduplicate datasets by key to avoid double-counting
        ref_df_dedup = self.reference_df.copy()
        ref_df_dedup['__key'] = ref_keys
        ref_df_dedup = ref_df_dedup.drop_duplicates(['__key'])
        
        comp_df_dedup = self.comparison_df.copy()
        comp_df_dedup['__key'] = comp_keys
        comp_df_dedup = comp_df_dedup.drop_duplicates(['__key'])
        
        # Create pair-level matches with proper deduplication
        pair_matches_data = []
        all_keys = set(ref_df_dedup['__key']) | set(comp_df_dedup['__key'])
        
        for key in all_keys:
            parts = key.split('|')
            parent_id = parts[0]
            rule = parts[1] 
            halogen = parts[2]
            k = int(parts[3])
            
            ref_match = key in set(ref_df_dedup['__key'])
            comp_match = key in set(comp_df_dedup['__key'])
            
            # Get product info if available (now properly deduplicated)
            if ref_match:
                ref_row = ref_df_dedup[ref_df_dedup['__key'] == key].iloc[0]
                product_id = ref_row.get('product_smiles', 'N/A')
                site_id = ref_row.get('site_id', 'N/A')
            elif comp_match:
                comp_row = comp_df_dedup[comp_df_dedup['__key'] == key].iloc[0]
                product_id = comp_row.get('product_smiles', 'N/A')
                site_id = comp_row.get('site_id', 'N/A')
            else:
                product_id = 'N/A'
                site_id = 'N/A'
            
            pair_matches_data.append({
                'parent_id': parent_id,
                'product_id': product_id,
                'site_id': site_id,
                'rule': rule,
                'halogen': halogen,
                'k': k,
                'matched': 1 if (ref_match and comp_match) else 0,
                'ref_only': 1 if (ref_match and not comp_match) else 0,
                'comp_only': 1 if (not ref_match and comp_match) else 0
            })
        
        pair_matches_df = pd.DataFrame(pair_matches_data)
        
        # Create rule confusion matrix using proper groupby aggregation
        rule_confusion_df = pair_matches_df.groupby(['rule', 'halogen', 'k'])[['matched', 'ref_only', 'comp_only']].sum().reset_index()
        rule_confusion_df = rule_confusion_df.rename(columns={
            'matched': 'tp',
            'ref_only': 'fn', 
            'comp_only': 'fp'
        })
        
        # Create parent summary with correct coverage calculation
        parent_summary_df = pair_matches_df.groupby('parent_id')[['matched', 'ref_only', 'comp_only']].sum().reset_index()
        parent_summary_df = parent_summary_df.rename(columns={
            'matched': 'tp',
            'ref_only': 'fn',
            'comp_only': 'fp'
        })
        # Calculate coverage as tp/(tp+fn) * 100 (note: denominator does not include fp)
        parent_summary_df['coverage_pct'] = parent_summary_df.apply(
            lambda row: (row['tp'] / (row['tp'] + row['fn']) * 100) if (row['tp'] + row['fn']) > 0 else 0, 
            axis=1
        )
        
        # Save CSV files
        pair_matches_path = os.path.join(output_dir, f"{base_name}_pair_matches.csv")
        rule_confusion_path = os.path.join(output_dir, f"{base_name}_rule_confusion.csv")
        parent_summary_path = os.path.join(output_dir, f"{base_name}_parent_summary.csv")
        
        pair_matches_df.to_csv(pair_matches_path, index=False)
        rule_confusion_df.to_csv(rule_confusion_path, index=False)
        parent_summary_df.to_csv(parent_summary_path, index=False)
        
        # Log results
        LOG.info(f"Business CSV comparison results saved:")
        LOG.info(f"  Pair matches ({len(pair_matches_df)} records): {pair_matches_path}")
        LOG.info(f"  Rule confusion ({len(rule_confusion_df)} records): {rule_confusion_path}")
        LOG.info(f"  Parent summary ({len(parent_summary_df)} records): {parent_summary_path}")
        
        return {
            'pair_matches_file': pair_matches_path,
            'rule_confusion_file': rule_confusion_path,
            'parent_summary_file': parent_summary_path,
            'pair_matches_count': len(pair_matches_df),
            'rule_confusion_count': len(rule_confusion_df),
            'parent_summary_count': len(parent_summary_df)
        }


class CSVComparisonPipeline:
    """Pipeline for comparing multiple CSV files and generating comprehensive reports."""
    
    def __init__(self, output_dir: str = "comparison_reports"):
        """
        Initialize comparison pipeline.
        
        Args:
            output_dir: Directory to save comparison reports
        """
        self.output_dir = output_dir
        self.comparisons = []
        
    def add_comparison(self, reference_path: str, comparison_path: str, 
                      name: str = None, key_columns: List[str] = None) -> CSVComparison:
        """
        Add a CSV comparison to the pipeline.
        
        Args:
            reference_path: Path to reference CSV
            comparison_path: Path to comparison CSV
            name: Name for this comparison (default: based on filenames)
            key_columns: Key columns for row matching
            
        Returns:
            CSVComparison instance
        """
        if name is None:
            ref_name = os.path.splitext(os.path.basename(reference_path))[0]
            comp_name = os.path.splitext(os.path.basename(comparison_path))[0]
            name = f"{ref_name}_vs_{comp_name}"
        
        comparison = CSVComparison(reference_path, comparison_path, key_columns)
        self.comparisons.append((name, comparison))
        
        LOG.info(f"Added comparison: {name}")
        return comparison
    
    def run_all_comparisons(self) -> Dict[str, Any]:
        """Run all comparisons in the pipeline."""
        LOG.info(f"Running {len(self.comparisons)} comparisons...")
        
        pipeline_results = {
            'timestamp': datetime.now().isoformat(),
            'total_comparisons': len(self.comparisons),
            'comparisons': {},
            'summary': {
                'all_passed': True,
                'failed_comparisons': [],
                'passed_comparisons': []
            }
        }
        
        for name, comparison in self.comparisons:
            try:
                LOG.info(f"Running comparison: {name}")
                result = comparison.run_comparison()
                pipeline_results['comparisons'][name] = result
                
                if result['overall_assessment']['files_match']:
                    pipeline_results['summary']['passed_comparisons'].append(name)
                else:
                    pipeline_results['summary']['failed_comparisons'].append(name)
                    pipeline_results['summary']['all_passed'] = False
                    
                # Save individual report
                report_path = os.path.join(self.output_dir, f"{name}_comparison.json")
                comparison.save_report(report_path)
                
            except Exception as e:
                LOG.error(f"Failed to run comparison {name}: {e}")
                pipeline_results['summary']['failed_comparisons'].append(name)
                pipeline_results['summary']['all_passed'] = False
        
        # Save pipeline summary
        os.makedirs(self.output_dir, exist_ok=True)
        summary_path = os.path.join(self.output_dir, "pipeline_summary.json")
        with open(summary_path, 'w', encoding='utf-8') as f:
            json.dump(pipeline_results, f, indent=2)
        
        LOG.info(f"Pipeline complete. {len(pipeline_results['summary']['passed_comparisons'])} passed, "
                   f"{len(pipeline_results['summary']['failed_comparisons'])} failed")
        
        return pipeline_results
    
    def compare_summary_files(self, reference_dir: str, comparison_dir: str) -> Dict[str, Any]:
        """
        Compare all summary CSV files between two directories.
        
        Args:
            reference_dir: Directory containing reference CSV files
            comparison_dir: Directory containing comparison CSV files
            
        Returns:
            Pipeline results
        """
        # Find all CSV files in reference directory
        ref_csv_files = []
        for root, dirs, files in os.walk(reference_dir):
            for file in files:
                if file.endswith('.csv'):
                    ref_csv_files.append(os.path.join(root, file))
        
        # Add comparisons for matching files
        for ref_file in ref_csv_files:
            rel_path = os.path.relpath(ref_file, reference_dir)
            comp_file = os.path.join(comparison_dir, rel_path)
            
            if os.path.exists(comp_file):
                name = rel_path.replace(os.sep, '_').replace('.csv', '')
                self.add_comparison(ref_file, comp_file, name)
            else:
                LOG.warning(f"Missing comparison file: {comp_file}")
        
        return self.run_all_comparisons()


def compare_csv_files(reference_path: str, comparison_path: str, 
                     output_path: str = None) -> Dict[str, Any]:
    """
    Convenience function to compare two CSV files.
    
    Args:
        reference_path: Path to reference CSV
        comparison_path: Path to comparison CSV
        output_path: Optional path to save comparison report
        
    Returns:
        Comparison results dictionary
    """
    comparison = CSVComparison(reference_path, comparison_path)
    results = comparison.run_comparison()
    
    if output_path:
        comparison.save_report(output_path)
    
    return results