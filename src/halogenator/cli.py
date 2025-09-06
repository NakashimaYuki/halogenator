# -*- coding: ascii -*-
"""Command line interface."""

import argparse
import os
import sys
import tempfile
import shutil
import yaml
from typing import Dict, Any

from .io_utils import load_names, resolve_names_to_smiles, write_smi, read_smi, write_sdf_with_props, write_table
from .standardize import std_from_smiles
from .enumerate_k1 import enumerate_k1_halogenation
from .report import generate_summary_report


def load_config(config_path: str) -> Dict[str, Any]:
    """Load configuration from YAML file."""
    try:
        with open(config_path, 'r', encoding='ascii') as f:
            return yaml.safe_load(f)
    except Exception as e:
        print(f"Error loading config {config_path}: {e}")
        sys.exit(1)


def cmd_ingest(config: Dict[str, Any]) -> None:
    """Ingest and standardize parent molecules."""
    io_config = config.get('io', {})
    standardize_config = config.get('standardize', {})
    
    names_file = io_config.get('names_file', 'data/input/flavonoids_p0_names.txt')
    smiles_file = io_config.get('smiles_file', 'data/output/p0/parents.smi')
    do_tautomer = standardize_config.get('do_tautomer', False)
    
    print("Starting ingest process...")
    
    # Check if parents.smi already exists
    existing_pairs = read_smi(smiles_file)
    if existing_pairs:
        print(f"Found existing {smiles_file} with {len(existing_pairs)} entries.")
        print("Using existing file. Delete it to re-resolve names.")
        return
    
    # Load compound names
    try:
        names = load_names(names_file)
        print(f"Loaded {len(names)} compound names from {names_file}")
    except FileNotFoundError:
        print(f"Names file not found: {names_file}")
        print("Please provide names file or create parents.smi manually")
        sys.exit(1)
    
    # Resolve names to SMILES
    try:
        print("Resolving names to SMILES via PubChem...")
        smiles_pairs = resolve_names_to_smiles(names)
        print(f"Successfully resolved {len(smiles_pairs)} compounds")
    except RuntimeError as e:
        print(f"Name resolution failed: {e}")
        sys.exit(1)
    
    # Standardize molecules
    print("Standardizing molecules...")
    standardized_pairs = []
    failed_count = 0
    
    for smiles, name in smiles_pairs:
        mol = std_from_smiles(smiles, do_tautomer)
        if mol is not None:
            # Get canonical SMILES from standardized molecule
            try:
                from rdkit import Chem
                std_smiles = Chem.MolToSmiles(mol, canonical=True)
            except:
                std_smiles = smiles  # Fallback to original
            
            standardized_pairs.append((std_smiles, name))
        else:
            print(f"Failed to standardize: {name}")
            failed_count += 1
    
    print(f"Standardization complete: {len(standardized_pairs)} success, {failed_count} failed")
    
    # Write to file
    write_smi(standardized_pairs, smiles_file)
    print(f"Wrote standardized parents to {smiles_file}")


def cmd_k1(config: Dict[str, Any], args=None) -> None:
    """Run k=1 halogenation enumeration."""
    io_config = config.get('io', {})
    halogens = config.get('halogens', ['F', 'Cl', 'Br', 'I'])
    rules = config.get('rules', ['R1', 'R2', 'R3', 'R4', 'R5'])
    
    # Determine input records based on subset selection
    subset = args.subset if args else 'flavonoids'
    records = []  # list of (smiles, name, parent_type)
    
    if subset == 'flavonoids':
        for smi, nm in read_smi(io_config.get('smiles_file', 'data/output/p0/parents.smi')):
            records.append((smi, nm, 'flavonoid'))
    elif subset == 'probes':
        if not os.path.exists('data/input/rule_probes.smi'):
            print("ERROR: Missing required file: data/input/rule_probes.smi")
            print("\nTo run --subset probes, create this file with sample probes:")
            print("c1ccc(O)cc1\tphenol")
            print("CCN\tethylamine") 
            print("c1ccc(C(=O)O)cc1\tbenzoic_acid")
            print("\nEach line: SMILES<tab>name")
            sys.exit(1)
        for smi, nm in read_smi('data/input/rule_probes.smi'):
            records.append((smi, nm, 'probe'))
    else:  # 'all' - combine both files
        for smi, nm in read_smi(io_config.get('smiles_file', 'data/output/p0/parents.smi')):
            records.append((smi, nm, 'flavonoid'))
        if not os.path.exists('data/input/rule_probes.smi'):
            print("ERROR: Missing required file: data/input/rule_probes.smi")
            print("\nTo run --subset all, create this file with sample probes:")
            print("c1ccc(O)cc1\tphenol")
            print("CCN\tethylamine") 
            print("c1ccc(C(=O)O)cc1\tbenzoic_acid")
            print("\nEach line: SMILES<tab>name")
            sys.exit(1)
        for smi, nm in read_smi('data/input/rule_probes.smi'):
            records.append((smi, nm, 'probe'))
    
    # Update output paths based on subset and outdir
    if subset != 'flavonoids':
        suffix = f"_{subset}"
        products_sdf = io_config.get('products_sdf', 'data/output/p0/products_k1.sdf').replace('.sdf', f'{suffix}.sdf')
        products_table = io_config.get('products_table', 'data/output/p0/products_k1.parquet').replace('.parquet', f'{suffix}.parquet')
    else:
        products_sdf = io_config.get('products_sdf', 'data/output/p0/products_k1.sdf') 
        products_table = io_config.get('products_table', 'data/output/p0/products_k1.parquet')
    
    # Override directory if --outdir provided
    if args and hasattr(args, 'outdir') and args.outdir:
        products_sdf = os.path.join(args.outdir, os.path.basename(products_sdf))
        products_table = os.path.join(args.outdir, os.path.basename(products_table))
        # Ensure output directory exists
        os.makedirs(args.outdir, exist_ok=True)
    
    print(f"Starting k=1 halogenation ({subset} subset)...")
    
    if not records:
        print(f"No parent molecules found for subset {subset}")
        print("Run 'ingest' command first or check input files")
        sys.exit(1)
    
    print(f"Loaded {len(records)} parent molecules")
    print(f"Using rules: {rules}")
    print(f"Using halogens: {halogens}")
    
    # Process each parent
    all_products = []
    all_records = []
    
    from rdkit import Chem
    
    for i, (smiles, name, parent_type) in enumerate(records, 1):
        print(f"Processing {i}/{len(records)}: {name}")
        
        # Parse molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"  Failed to parse SMILES: {smiles}")
            continue
        
        # Enumerate products
        try:
            products = enumerate_k1_halogenation(mol, halogens, rules, config)
            print(f"  Generated {len(products)} products")
            
            # Convert to records and collect molecules
            for product_mol, props in products:
                all_products.append(product_mol)
                
                # Add molecule identifiers to record
                product_smiles = Chem.MolToSmiles(product_mol, canonical=True)
                record = props.copy()
                record['product_smiles'] = product_smiles
                record['parent_name'] = name
                record['parent_type'] = parent_type
                all_records.append(record)
                
        except Exception as e:
            print(f"  Error processing {name}: {e}")
            continue
    
    print(f"\nTotal products generated: {len(all_products)}")
    
    # Write outputs
    if all_products:
        # Set properties on molecules for SDF
        for i, (mol, record) in enumerate(zip(all_products, all_records)):
            for key, value in record.items():
                if isinstance(value, (str, int, float)):
                    mol.SetProp(key, str(value))
        
        write_sdf_with_props(all_products, products_sdf)
        print(f"Wrote SDF to {products_sdf}")
        
        write_table(all_records, products_table) 
        print(f"Wrote table to {products_table}")
    else:
        print("No products generated")


def cmd_report(config: Dict[str, Any], args=None) -> None:
    """Generate summary report."""
    temp_dir_to_cleanup = None  # Initialize for proper scope
    io_conf = dict(config.get('io', {}))
    subset = args.subset if args else 'flavonoids'
    suffix = '' if subset == 'flavonoids' else f'_{subset}'
    
    # Choose products table by subset  
    base_tbl = io_conf.get('products_table', 'data/output/p0/products_k1.parquet')
    io_conf['products_table'] = base_tbl.replace('.parquet', f'{suffix}.parquet')
    
    # Choose summary outputs by subset
    base_sum = io_conf.get('summary_csv', 'data/output/p0/summary_k1.csv') 
    io_conf['summary_csv'] = base_sum.replace('.csv', f'{suffix}.csv')
    
    # Override directory if --outdir provided
    if args and hasattr(args, 'outdir') and args.outdir:
        io_conf['products_table'] = os.path.join(args.outdir, os.path.basename(io_conf['products_table']))
        io_conf['summary_csv'] = os.path.join(args.outdir, os.path.basename(io_conf['summary_csv']))
        # Ensure output directory exists
        os.makedirs(args.outdir, exist_ok=True)
    
    # Choose parent smiles file by subset
    base_smiles = io_conf.get('smiles_file', 'data/output/p0/parents.smi')
    if subset == 'flavonoids':
        io_conf['smiles_file'] = base_smiles
    elif subset == 'probes':
        io_conf['smiles_file'] = 'data/input/rule_probes.smi'
    else:  # 'all'
        # Create a combined parents_all.smi file safely
        temp_dir_to_cleanup = None
        if args and hasattr(args, 'outdir') and args.outdir:
            combined = os.path.join(args.outdir, 'parents_all.smi')
        else:
            # Use temporary directory to avoid polluting default data directory
            temp_dir_to_cleanup = tempfile.mkdtemp()
            combined = os.path.join(temp_dir_to_cleanup, 'parents_all.smi')
        
        # Deduplicate lines to avoid duplicate entries
        seen = set()
        with open(combined, 'w', encoding='utf-8') as w:
            for path in [base_smiles, 'data/input/rule_probes.smi']:
                if os.path.exists(path):
                    with open(path, 'r', encoding='utf-8') as r:
                        for line in r:
                            s = line.strip()
                            if not s:
                                continue
                            # Use original line as deduplication key to avoid parsing complexity
                            if s not in seen:
                                seen.add(s)
                                w.write(s + '\n')
        io_conf['smiles_file'] = combined
    
    # Override into a shallow-copied config object
    cfg2 = dict(config)
    cfg2['io'] = io_conf
    
    print(f"Generating summary report ({subset} subset)...")
    try:
        generate_summary_report(cfg2, subset=subset)
        print("Report generation complete")
    except Exception as e:
        print(f"Error generating report: {e}")
        sys.exit(1)
    finally:
        # Cleanup temporary directory if created
        if subset == 'all' and temp_dir_to_cleanup and os.path.exists(temp_dir_to_cleanup):
            shutil.rmtree(temp_dir_to_cleanup, ignore_errors=True)


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description='Halogenator: k=1 halogen substitution for flavonoids',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Ingest command
    ingest_parser = subparsers.add_parser('ingest', help='Ingest and standardize parent molecules')
    ingest_parser.add_argument('-c', '--config', required=True, help='Configuration file path')
    
    # K1 command
    k1_parser = subparsers.add_parser('k1', help='Run k=1 halogenation enumeration') 
    k1_parser.add_argument('-c', '--config', required=True, help='Configuration file path')
    k1_parser.add_argument('--subset', choices=['flavonoids', 'probes', 'all'], default='flavonoids', help='Which molecular set to process')
    k1_parser.add_argument('--outdir', help='Override output directory for products')
    
    # Report command
    report_parser = subparsers.add_parser('report', help='Generate summary report')
    report_parser.add_argument('-c', '--config', required=True, help='Configuration file path')
    report_parser.add_argument('--subset', choices=['flavonoids', 'probes', 'all'], default='flavonoids', help='Which molecular set to process')
    report_parser.add_argument('--outdir', help='Override output directory for reports')
    
    # Parse arguments
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    # Load configuration
    config = load_config(args.config)
    
    # Dispatch to appropriate command
    if args.command == 'ingest':
        cmd_ingest(config)
    elif args.command == 'k1':
        cmd_k1(config, args)
    elif args.command == 'report':
        cmd_report(config, args)
    else:
        print(f"Unknown command: {args.command}")
        sys.exit(1)


if __name__ == "__main__":
    main()
