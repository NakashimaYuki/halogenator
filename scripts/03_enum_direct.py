#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Direct Enumeration Script - bypasses CLI issues
"""

import os
import sys
import glob
import yaml
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Load configuration
config_path = "configs/flavonoids_k2_prod.yaml"
with open(config_path, 'r', encoding='utf-8') as f:
    config = yaml.safe_load(f)

# Find shard files
shard_dir = "data/work/shards"
shard_files = sorted(glob.glob(os.path.join(shard_dir, "flav_shard_*.smi")))

logger.info(f"Found {len(shard_files)} shard files")
logger.info(f"Configuration: {config_path}")

# Import halogenator components
from halogenator import io_utils
from halogenator.enumerate_k import enumerate_k_bfs

# Process each shard
for shard_file in shard_files:
    shard_name = Path(shard_file).stem
    outdir = f"data/output/cnpd_flav_k2/{shard_name}"
    os.makedirs(outdir, exist_ok=True)

    logger.info(f"Processing {shard_name}...")
    logger.info(f"  Input: {shard_file}")
    logger.info(f"  Output: {outdir}")

    try:
        # Load parent molecules from .smi file
        parents = []
        with open(shard_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= 2:
                    smiles, name = parts[0], parts[1]
                    parents.append({'smiles': smiles, 'name': name})

        logger.info(f"  Loaded {len(parents)} parent molecules")

        # Run enumeration with config
        k_max = config.get('k_max', 2)
        halogens = config.get('halogens', ['F', 'Cl', 'Br', 'I'])
        rules = config.get('rules', ['R1', 'R2', 'R3', 'R6_methyl'])

        dedup_enable = config.get('dedup', {}).get('enable', True)
        sugar_mode = config.get('sugar_cfg', {}).get('mode', 'heuristic')

        logger.info(f"  k_max={k_max}, halogens={halogens}, rules={rules}")
        logger.info(f"  dedup={dedup_enable}, sugar_mask={sugar_mode}")

        # Create minimal enum config
        enum_config = {
            'k_max': k_max,
            'halogens': halogens,
            'rules': rules,
            'dedup': {'enable': dedup_enable},
            'sugar_cfg': {'mode': sugar_mode},
            'rules_cfg': config.get('rules_cfg', {}),
            'engine_cfg': config.get('engine_cfg', {}),
            'constraints': config.get('constraints', {}),
        }

        # Note: Direct Python API call may need adjustment based on actual halogenator API
        # This is a template that shows the intent

        logger.info(f"  Running enumeration...")

        # For now, log that this needs proper API integration
        logger.warning(f"  Direct Python API integration needed - using CLI workaround instead")

        # Try using subset-based approach with temporary configuration
        import subprocess
        import tempfile

        # Create temporary config for this shard
        shard_config = config.copy()
        shard_config['data_path'] = shard_file  # Try adding data_path

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False, encoding='utf-8') as f:
            yaml.dump(shard_config, f)
            temp_config = f.name

        try:
            # Try running with modified config
            cmd = [
                'python', '-m', 'halogenator.cli',
                '--rdkit-threads', str(config.get('engine_cfg', {}).get('rdkit_threads', 4)),
                'enum',
                '-c', temp_config,
                '--outdir', outdir
            ]

            logger.info(f"  Command: {' '.join(cmd)}")

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

            if result.returncode == 0:
                logger.info(f"  SUCCESS: {shard_name}")
            else:
                logger.error(f"  FAILED: {shard_name}")
                logger.error(f"  stderr: {result.stderr[:500]}")

        finally:
            try:
                os.unlink(temp_config)
            except:
                pass

    except Exception as e:
        logger.error(f"  ERROR processing {shard_name}: {e}")
        import traceback
        traceback.print_exc()

logger.info("Enumeration complete")
