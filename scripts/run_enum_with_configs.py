#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run enumeration by creating per-shard configuration files
"""

import os
import sys
import glob
import yaml
import subprocess
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Load base config
base_config_path = 'configs/flavonoids_k2_prod.yaml'
with open(base_config_path, 'r', encoding='utf-8') as f:
    base_config = yaml.safe_load(f)

# Find shards
shard_files = sorted(glob.glob('data/work/shards/flav_shard_*.smi'))
logger.info(f"Found {len(shard_files)} shard files")

# Create temp config dir
temp_config_dir = 'data/work/temp_configs'
os.makedirs(temp_config_dir, exist_ok=True)

# Process each shard
for shard_idx, shard_file in enumerate(shard_files, 1):
    shard_name = Path(shard_file).stem
    outdir = f'data/output/cnpd_flav_k2/{shard_name}'
    os.makedirs(outdir, exist_ok=True)

    logger.info(f"[{shard_idx}/{len(shard_files)}] Processing {shard_name}...")
    logger.info(f"  Input: {shard_file}")
    logger.info(f"  Output: {outdir}")

    # Create shard-specific config
    shard_config = base_config.copy()
    shard_config['io'] = {'smiles_file': shard_file}
    shard_config['subset'] = 'cnpd_etcm_flavonoids'  # Keep subset name

    # Write temp config
    temp_config_path = os.path.join(temp_config_dir, f'{shard_name}_config.yaml')
    with open(temp_config_path, 'w', encoding='utf-8') as f:
        yaml.dump(shard_config, f)

    # Run enumeration
    rdkit_threads = base_config.get('engine_cfg', {}).get('rdkit_threads', 4)

    cmd = [
        'python', '-m', 'halogenator.cli',
        '--rdkit-threads', str(rdkit_threads),
        'enum',
        '-c', temp_config_path,
        '--outdir', outdir
    ]

    logger.info(f"  Command: {' '.join(cmd)}")

    log_file = f'data/output/cnpd_flav_k2_logs/{shard_name}.log'

    try:
        with open(log_file, 'w') as log_f:
            result = subprocess.run(
                cmd,
                stdout=log_f,
                stderr=subprocess.STDOUT,
                timeout=7200  # 2 hour timeout per shard
            )

        if result.returncode == 0:
            logger.info(f"  ✓ SUCCESS: {shard_name}")
            # Write success marker
            with open(f'data/output/cnpd_flav_k2_logs/completed_shards.txt', 'a') as f:
                f.write(f'{shard_name}\n')
        else:
            logger.error(f"  ✗ FAILED: {shard_name} (see {log_file})")
            # Write failure marker
            with open(f'data/output/cnpd_flav_k2_logs/failed_shards.txt', 'a') as f:
                f.write(f'{shard_name}\n')

    except subprocess.TimeoutExpired:
        logger.error(f"  ✗ TIMEOUT: {shard_name} (exceeded 2 hours)")
        with open(f'data/output/cnpd_flav_k2_logs/failed_shards.txt', 'a') as f:
            f.write(f'{shard_name} (timeout)\n')
    except Exception as e:
        logger.error(f"  ✗ ERROR: {shard_name}: {e}")
        with open(f'data/output/cnpd_flav_k2_logs/failed_shards.txt', 'a') as f:
            f.write(f'{shard_name} (error: {e})\n')

logger.info("="*80)
logger.info("ENUMERATION COMPLETE")
logger.info("="*80)

# Summary
try:
    with open('data/output/cnpd_flav_k2_logs/completed_shards.txt', 'r') as f:
        completed = len(f.readlines())
except:
    completed = 0

try:
    with open('data/output/cnpd_flav_k2_logs/failed_shards.txt', 'r') as f:
        failed = len(f.readlines())
except:
    failed = 0

logger.info(f"Completed: {completed}/{len(shard_files)}")
logger.info(f"Failed: {failed}/{len(shard_files)}")
