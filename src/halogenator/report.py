# -*- coding: ascii -*-
"""Report generation utilities."""

import os
import sys
import platform
from datetime import datetime
from typing import Dict, Any, List
from collections import Counter, defaultdict
import pandas as pd

from .io_utils import read_smi, read_table
from .schema import validate_products_records


def generate_summary_report(config: Dict[str, Any], subset: str = 'flavonoids') -> None:
    """Generate summary report from k=1 results."""
    io_config = config.get('io', {})
    
    # File paths
    parents_file = io_config.get('smiles_file', 'data/output/p0/parents.smi')
    products_table = io_config.get('products_table', 'data/output/p0/products_k1.parquet')
    summary_file = io_config.get('summary_csv', 'data/output/p0/summary_k1.csv')
    
    # Read parent molecules
    parent_pairs = read_smi(parents_file)
    parent_count = len(parent_pairs)
    
    # Read products
    product_records = read_table(products_table)
    product_count = len(product_records)
    
    # Validate schema (check all records for consistency)
    if product_records:
        validate_products_records(product_records)
    
    # Generate statistics (pass subset for metadata)
    config_with_subset = dict(config)
    config_with_subset['subset'] = subset
    stats = _compute_statistics(parent_pairs, product_records, config_with_subset)
    
    # Write summary CSV
    _write_summary_csv(stats, summary_file)
    
    # Write rule x halogen pivot tables
    pivot_file = summary_file.replace('.csv', '_pivot.csv')
    _write_pivot_table(stats['rule_halogen_counts'], pivot_file)
    
    # Write type-specific pivot tables only for 'all' subset
    if subset == 'all':
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
    
    for record in product_records:
        rule = record.get('rule', 'Unknown')
        halogen = record.get('halogen', 'Unknown')
        parent_key = record.get('parent_inchikey', record.get('parent_smiles', 'Unknown'))
        parent_type = record.get('parent_type', 'unknown')
        
        # Overall counts
        rule_counts[rule] += 1
        halogen_counts[halogen] += 1
        rule_halogen_counts[rule][halogen] += 1
        parent_product_counts[parent_key] += 1
        
        # Type-specific counts
        if parent_type == 'flavonoid':
            rule_counts_flavonoids[rule] += 1
            halogen_counts_flavonoids[halogen] += 1
            rule_halogen_counts_flavonoids[rule][halogen] += 1
        elif parent_type == 'probe':
            rule_counts_probes[rule] += 1
            halogen_counts_probes[halogen] += 1
            rule_halogen_counts_probes[rule][halogen] += 1
    
    # Dual-track validation: SMILES and InChIKey
    parent_smiles_to_name = {smiles: name for smiles, name in parent_pairs}
    all_parent_smiles = set(parent_smiles_to_name.keys())
    parents_with_products_smiles = set()
    parents_with_products_ikeys = set()
    
    # Build InChIKey to name mapping (one-time for parent pairs)
    parent_smiles_to_ikey = {}
    parent_ikey_to_name = {}
    for smiles, name in parent_pairs:
        # Try to get InChIKey for parent molecules
        from .standardize import to_inchikey
        from rdkit import Chem
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                ikey = to_inchikey(mol)
                parent_smiles_to_ikey[smiles] = ikey
                parent_ikey_to_name[ikey] = name
        except:
            pass
    
    for record in product_records:
        parent_smiles = record.get('parent_smiles', 'Unknown')
        parent_ikey = record.get('parent_inchikey', 'Unknown')
        parents_with_products_smiles.add(parent_smiles)
        parents_with_products_ikeys.add(parent_ikey)
    
    parents_without_products_smiles = all_parent_smiles - parents_with_products_smiles
    # Map back to names for reporting (primary: SMILES-based)
    parents_without_products = sorted([parent_smiles_to_name.get(smiles, smiles) 
                                     for smiles in parents_without_products_smiles])
    
    # Diagnostics: check SMILES vs InChIKey consistency with set differences
    smiles_count = len(parents_with_products_smiles)
    ikey_count = len(parents_with_products_ikeys)
    diff_smiles_vs_inchikey = abs(smiles_count - ikey_count)
    
    # Find set differences for diagnostics using precomputed mapping
    smiles_based_ikeys = set()
    for smiles in parents_with_products_smiles:
        ikey = parent_smiles_to_ikey.get(smiles)
        if ikey:
            smiles_based_ikeys.add(ikey)
    
    missing_in_smiles = parents_with_products_ikeys - smiles_based_ikeys
    missing_in_ikeys = smiles_based_ikeys - parents_with_products_ikeys
    
    # Get examples (up to 3 each)
    missing_in_smiles_examples = sorted([parent_ikey_to_name.get(ikey, ikey) 
                                       for ikey in list(missing_in_smiles)[:3]])
    missing_in_ikeys_examples = sorted([parent_ikey_to_name.get(ikey, ikey) 
                                      for ikey in list(missing_in_ikeys)[:3]])
    
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
    
    return {
        'parent_count': parent_count,
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
        'unique_parents_with_products': len(parent_product_counts),
        'parents_without_products': parents_without_products,
        'diagnostics': {
            'diff_smiles_vs_inchikey': diff_smiles_vs_inchikey,
            'smiles_unique_parents': smiles_count,
            'inchikey_unique_parents': ikey_count,
            'missing_in_smiles_count': len(missing_in_smiles),
            'missing_in_ikeys_count': len(missing_in_ikeys),
            'missing_in_smiles_examples': missing_in_smiles_examples,
            'missing_in_ikeys_examples': missing_in_ikeys_examples
        },
        'config': config
    }


def _write_summary_csv(stats: Dict[str, Any], output_path: str) -> None:
    """Write summary statistics to CSV.""" 
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Create summary rows
    summary_rows = []
    
    # Runtime metadata
    try:
        import rdkit
        rdkit_version = rdkit.__version__
    except:
        rdkit_version = 'unknown'
        
    try:
        import pandas as pd_meta
        pandas_version = pd_meta.__version__
    except:
        pandas_version = 'unknown'
        
    try:
        import pyarrow
        pyarrow_version = pyarrow.__version__
    except:
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
    
    zero_product_list = ', '.join(stats['parents_without_products']) if stats['parents_without_products'] else 'None'
    summary_rows.append({
        'Category': 'Overall',
        'Metric': 'Parents with Zero Products',
        'Value': len(stats['parents_without_products']),
        'Description': f'Parents that generated no products: {zero_product_list}'
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
    
    # Rule x Halogen matrix
    for rule, halogen_dict in stats['rule_halogen_counts'].items():
        for halogen, count in halogen_dict.items():
            summary_rows.append({
                'Category': 'Rule x Halogen',
                'Metric': f'{rule} + {halogen}',
                'Value': count,
                'Description': f'Products from rule {rule} with halogen {halogen}'
            })
    
    # Diagnostics section  
    diagnostics = stats.get('diagnostics', {})
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
    
    # Write to CSV
    df = pd.DataFrame(summary_rows)
    df.to_csv(output_path, index=False)


def _write_pivot_table(rule_halogen_counts: Dict[str, Dict[str, int]], output_path: str) -> None:
    """Write rule x halogen pivot table to CSV."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Create pivot table data
    rules = ['R1', 'R2', 'R3', 'R4', 'R5']
    halogens = ['F', 'Cl', 'Br', 'I']
    
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
    print("\n" + "="*60)
    print("HALOGENATOR P0 SUMMARY REPORT")
    print("="*60)
    
    print(f"Parent molecules processed: {stats['parent_count']}")
    print(f"Total products generated: {stats['product_count']}")
    print(f"Products passing QC: {stats['sanitize_ok_count']}")
    print(f"Unique parents with products: {stats['unique_parents_with_products']}")
    
    # Report parents with zero products for transparency
    zero_product_parents = stats['parents_without_products']
    if zero_product_parents:
        print(f"Parents with zero products ({len(zero_product_parents)}): {', '.join(zero_product_parents)}")
    else:
        print("All parents produced at least one product")
    
    # Rules activity status (overall)
    all_rules = ['R1', 'R2', 'R3', 'R4', 'R5']
    rule_status = ', '.join(f'{r}={stats["rule_counts"].get(r, 0)}' for r in all_rules)
    print(f"\nRules active (overall): {rule_status}")
    
    # Rules activity on flavonoids  
    rule_status_flavonoids = ', '.join(f'{r}={stats["rule_counts_flavonoids"].get(r, 0)}' for r in all_rules)
    print(f"Rules active on flavonoids: {rule_status_flavonoids}")
    
    # Rules activity on probes
    rule_status_probes = ', '.join(f'{r}={stats["rule_counts_probes"].get(r, 0)}' for r in all_rules)
    print(f"Rules active on probes: {rule_status_probes}")
    
    inactive_rules_flavonoids = [r for r in all_rules if stats["rule_counts_flavonoids"].get(r, 0) == 0]
    if inactive_rules_flavonoids:
        print(f"Rules inactive on flavonoids: {', '.join(inactive_rules_flavonoids)} (expected for R4/R5 - flavonoids rarely contain -NH/-COOH)")
    
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
    for rule in ['R1', 'R2', 'R3', 'R4', 'R5']:
        counts = stats['rule_halogen_counts'].get(rule, {})
        f_count = counts.get('F', 0)
        cl_count = counts.get('Cl', 0)
        br_count = counts.get('Br', 0)
        i_count = counts.get('I', 0)
        total = f_count + cl_count + br_count + i_count
        if total > 0:
            print(f"{rule:<6} {f_count:<6} {cl_count:<6} {br_count:<6} {i_count:<6} {total:<6}")
    
    print("\nNo fatal RDKit errors; pipeline stable")
    print("="*60)