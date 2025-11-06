#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simple enumeration runner - processes shard files directly
"""

import os
import sys
import glob
import yaml
import pandas as pd
import logging
from pathlib import Path

# Add src to path
sys.path.insert(0, 'src')

from rdkit import Chem
from halogenator.enumerate_k import enumerate_k_bfs
from halogenator import io_utils
from halogenator.schema import make_history_step

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Load config
with open('configs/flavonoids_k2_prod.yaml', 'r', encoding='utf-8') as f:
    config = yaml.safe_load(f)

# Find shards
shard_files = sorted(glob.glob('data/work/shards/flav_shard_*.smi'))
logger.info(f"Found {len(shard_files)} shard files")

# Extract config parameters
k_max = config.get('k_max', 2)
halogens = config.get('halogens', ['F', 'Cl', 'Br', 'I'])
rules_list = config.get('rules', ['R1', 'R2', 'R3', 'R6_methyl'])
dedup_enabled = config.get('dedup', {}).get('enable', True)
sugar_mode = config.get('sugar_cfg', {}).get('mode', 'heuristic')
r2_fallback = config.get('rules_cfg', {}).get('r2', {}).get('fallback', True)
r6_enable_step = config.get('rules_cfg', {}).get('r6_methyl', {}).get('enable_step', True)
r6_enable_macro = config.get('rules_cfg', {}).get('r6_methyl', {}).get('enable_macro', True)

logger.info(f"Config: k_max={k_max}, halogens={halogens}")
logger.info(f"Rules: {rules_list}, dedup={dedup_enabled}, sugar={sugar_mode}")

# Process each shard
for shard_idx, shard_file in enumerate(shard_files, 1):
    shard_name = Path(shard_file).stem
    outdir = f'data/output/cnpd_flav_k2/{shard_name}'
    os.makedirs(outdir, exist_ok=True)

    logger.info(f"[{shard_idx}/{len(shard_files)}] Processing {shard_name}...")

    # Load parent SMI file
    parents = []
    with open(shard_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                smiles, name = parts[0], parts[1]
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    parents.append({
                        'smiles': smiles,
                        'name': name,
                        'mol': mol,
                        'parent_inchikey': Chem.inchi.MolToInchiKey(mol)
                    })

    logger.info(f"  Loaded {len(parents)} valid parents")

    # Run enumeration for each parent
    all_products = []

    for pidx, parent in enumerate(parents):
        if (pidx + 1) % 100 == 0:
            logger.info(f"    Processing parent {pidx+1}/{len(parents)}...")

        try:
            # Build enum_config for this parent
            enum_cfg = {
                'k_max': k_max,
                'halogens': halogens,
                'rules': rules_list,
                'dedup': {'enable': dedup_enabled},
                'sugar_cfg': {'mode': sugar_mode},
                'rules_cfg': {
                    'r2': {'fallback': r2_fallback},
                    'r6_methyl': {
                        'enable_step': r6_enable_step,
                        'enable_macro': r6_enable_macro,
                        'allow_on_methoxy': True
                    }
                }
            }

            # Run enumeration
            results = enumerate_k_bfs(
                parent['mol'],
                parent_name=parent['name'],
                parent_inchikey=parent['parent_inchikey'],
                enum_cfg=enum_cfg
            )

            # Collect products
            for result in results:
                product_dict = {
                    'smiles': Chem.MolToSmiles(result['mol']),
                    'inchikey': Chem.inchi.MolToInchiKey(result['mol']),
                    'parent_name': parent['name'],
                    'parent_inchikey': parent['parent_inchikey'],
                    'k': result.get('k', 0),
                    'k_ops': result.get('k_ops', 0),
                    'k_atoms': result.get('k_atoms', 0),
                    'rule': result.get('rule', ''),
                    'halogen': result.get('halogen', ''),
                    'macro_label': result.get('macro_label', ''),
                    'substitutions': str(result.get('substitutions', []))
                }
                all_products.append(product_dict)

        except Exception as e:
            logger.error(f"    Error processing {parent['name']}: {e}")
            continue

    # Save products
    if all_products:
        df = pd.DataFrame(all_products)

        # Deduplicate if enabled
        if dedup_enabled:
            initial_count = len(df)
            df = df.drop_duplicates(subset=['inchikey'], keep='first')
            logger.info(f"  Deduplication: {initial_count} → {len(df)} products")

        # Save to parquet
        output_file = os.path.join(outdir, 'products_k2.parquet')
        df.to_parquet(output_file, index=False)
        logger.info(f"  Saved {len(df)} products to {output_file}")

        # Save by-rule CSV
        if config.get('output', {}).get('by_rule_csv', True):
            if 'rule' in df.columns and 'k' in df.columns:
                by_rule = df.groupby(['rule', 'halogen', 'k']).size().reset_index(name='count')
                by_rule.to_csv(os.path.join(outdir, 'by_rule.csv'), index=False)
    else:
        logger.warning(f"  No products generated for {shard_name}")

logger.info("="*80)
logger.info("ENUMERATION COMPLETE")
logger.info("="*80)
