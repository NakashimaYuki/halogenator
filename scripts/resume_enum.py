#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Resume enumeration for incomplete shards
Processes only the remainder files and sub-shards
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


def main():
    # Configuration
    base_config_path = 'configs/flavonoids_k2_prod.yaml'
    shard_dir = 'data/work/shards'
    output_root = 'data/output/cnpd_flav_k2_resumed'
    log_dir = 'data/output/cnpd_flav_k2_logs'
    timeout_seconds = 21600  # 6 hours

    # Load base config
    with open(base_config_path, 'r', encoding='utf-8') as f:
        base_config = yaml.safe_load(f)

    rdkit_threads = base_config.get('engine_cfg', {}).get('rdkit_threads', 8)

    # Find shards to process (remainder files and sub-shards)
    shard_patterns = [
        'flav_shard_0001_remainder.smi',
        'flav_shard_0002a.smi',
        'flav_shard_0002b.smi',
        'flav_shard_0002c.smi',
        'flav_shard_0002d.smi',
        'flav_shard_0002e.smi',
        'flav_shard_0002f.smi',
        'flav_shard_0003_remainder.smi',
    ]

    shard_files = []
    for pattern in shard_patterns:
        full_path = os.path.join(shard_dir, pattern)
        if os.path.exists(full_path):
            shard_files.append(full_path)
        else:
            logger.warning(f"Shard file not found: {full_path}")

    if not shard_files:
        logger.error("No shard files found!")
        return 1

    logger.info("="*80)
    logger.info("RESUME ENUMERATION - INCOMPLETE SHARDS")
    logger.info("="*80)
    logger.info(f"Found {len(shard_files)} shard files to process:")
    for sf in shard_files:
        logger.info(f"  - {os.path.basename(sf)}")
    logger.info(f"Config: {base_config_path}")
    logger.info(f"Timeout per shard: {timeout_seconds}s ({timeout_seconds/3600:.1f} hours)")
    logger.info(f"RDKit threads: {rdkit_threads}")
    logger.info("="*80)

    # Create output and log directories
    os.makedirs(output_root, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    # Create temp config dir
    temp_config_dir = 'data/work/temp_configs_resumed'
    os.makedirs(temp_config_dir, exist_ok=True)

    # Track results
    completed = []
    failed = []

    # Process each shard
    for shard_idx, shard_file in enumerate(shard_files, 1):
        shard_name = Path(shard_file).stem
        outdir = os.path.join(output_root, shard_name)
        os.makedirs(outdir, exist_ok=True)

        logger.info("")
        logger.info(f"[{shard_idx}/{len(shard_files)}] Processing {shard_name}...")
        logger.info(f"  Input: {shard_file}")
        logger.info(f"  Output: {outdir}")

        # Create shard-specific config
        shard_config = base_config.copy()
        shard_config['io'] = {'smiles_file': shard_file}
        shard_config['subset'] = 'cnpd_etcm_flavonoids'

        # Write temp config
        temp_config_path = os.path.join(temp_config_dir, f'{shard_name}_config.yaml')
        with open(temp_config_path, 'w', encoding='utf-8') as f:
            yaml.dump(shard_config, f)

        # Run enumeration
        cmd = [
            'python', '-m', 'halogenator.cli',
            '--rdkit-threads', str(rdkit_threads),
            'enum',
            '-c', temp_config_path,
            '--outdir', outdir
        ]

        logger.info(f"  Command: {' '.join(cmd)}")

        log_file = os.path.join(log_dir, f'{shard_name}.log')

        try:
            with open(log_file, 'w', encoding='utf-8') as log_f:
                result = subprocess.run(
                    cmd,
                    stdout=log_f,
                    stderr=subprocess.STDOUT,
                    timeout=timeout_seconds
                )

            if result.returncode == 0:
                logger.info(f"  SUCCESS: {shard_name}")
                completed.append(shard_name)

                # Write success marker
                with open(os.path.join(log_dir, 'completed_shards_resumed.txt'), 'a', encoding='utf-8') as f:
                    f.write(f'{shard_name}\n')
            else:
                logger.error(f"  FAILED: {shard_name} (see {log_file})")
                failed.append(shard_name)

                # Write failure marker
                with open(os.path.join(log_dir, 'failed_shards_resumed.txt'), 'a', encoding='utf-8') as f:
                    f.write(f'{shard_name}\n')

        except subprocess.TimeoutExpired:
            logger.error(f"  TIMEOUT: {shard_name} (exceeded {timeout_seconds/3600:.1f} hours)")
            failed.append(shard_name)

            with open(os.path.join(log_dir, 'failed_shards_resumed.txt'), 'a', encoding='utf-8') as f:
                f.write(f'{shard_name} (timeout)\n')

        except Exception as e:
            logger.error(f"  ERROR: {shard_name}: {e}")
            failed.append(shard_name)

            with open(os.path.join(log_dir, 'failed_shards_resumed.txt'), 'a', encoding='utf-8') as f:
                f.write(f'{shard_name} (error: {e})\n')

    # Summary
    logger.info("")
    logger.info("="*80)
    logger.info("RESUME ENUMERATION COMPLETE")
    logger.info("="*80)
    logger.info(f"Completed: {len(completed)}/{len(shard_files)}")
    logger.info(f"Failed: {len(failed)}/{len(shard_files)}")

    if completed:
        logger.info("\nCompleted shards:")
        for shard in completed:
            logger.info(f"  - {shard}")

    if failed:
        logger.info("\nFailed shards:")
        for shard in failed:
            logger.info(f"  - {shard}")

    logger.info("")
    logger.info(f"Output directory: {output_root}")
    logger.info(f"Log directory: {log_dir}")

    return 0 if len(failed) == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
