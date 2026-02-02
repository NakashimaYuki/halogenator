# Linux Server Setup Guide

## Quick Start

### Option 1: Using Conda (Recommended)

```bash
# 1. Create conda environment from file
conda env create -f environment.yml

# 2. Activate environment
conda activate halogenator

# 3. Verify installation
python -c "from rdkit import Chem; print('RDKit OK')"
python -c "import pyarrow; print('PyArrow OK')"
```

### Option 2: Using pip + venv

```bash
# 1. Create virtual environment
python3.9 -m venv venv

# 2. Activate environment
source venv/bin/activate

# 3. Upgrade pip
pip install --upgrade pip

# 4. Install dependencies
pip install -r requirements.txt

# 5. Verify installation
python -c "from rdkit import Chem; print('RDKit OK')"
python -c "import pyarrow; print('PyArrow OK')"
```

## System Requirements

- **OS**: Ubuntu 20.04+ / CentOS 7+ / RHEL 8+
- **Python**: 3.9 or higher
- **RAM**: 32GB+ recommended for large-scale processing
- **Storage**: 500GB+ for data files
- **CPU**: 16+ cores recommended

## Installing Miniconda (if not already installed)

```bash
# Download Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run installer
bash Miniconda3-latest-Linux-x86_64.sh

# Follow prompts, then restart shell or run:
source ~/.bashrc

# Verify conda installation
conda --version
```

## Data Transfer from Windows

If you have data on Windows that needs to be transferred:

```bash
# From Linux server, use rsync
rsync -avz --progress user@windows-host:/path/to/data/ ./data/

# Or use scp
scp -r user@windows-host:/path/to/data/ ./data/
```

## Testing the Environment

```bash
# Activate environment
conda activate halogenator

# Run a quick test
python -c "
from rdkit import Chem
from rdkit.Chem import Descriptors
import pyarrow.parquet as pq
import pandas as pd
import numpy as np
print('All imports successful!')
"
```

## Processing Shuffled Chunks

The shuffled chunks are ready for processing:

```bash
# Location of shuffled chunks
ls data/output/transforms/polyphenol-2X_SHUFFLED_500K/

# Each chunk should be ~94MB
# 28 chunks total, 500K rows each
# Expected processing: 3-4 hours per chunk
```

## Running Transforms

Example command to process a single chunk:

```bash
python scripts/08_transform_library_v2.py apply \
  --input data/output/transforms/polyphenol-2X_SHUFFLED_500K/chunk_000_input.parquet \
  --outdir data/output/transforms/SHUFFLED_CHUNK_000 \
  --xf-config configs/transforms.yaml \
  --xf-name FG_PHENOL_OH__OH__TO__OMe \
  --workers 16 \
  --max-in-flight 6 \
  --batch-size 50000 \
  --use-bloom-filter \
  --bloom-expected-items 5000000
```

## Performance Optimization for Linux

### CPU Performance

```bash
# Check CPU governor (should be 'performance' for computational work)
cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor

# Set to performance mode (requires root)
sudo cpupower frequency-set -g performance
```

### NUMA Settings

For multi-socket systems:

```bash
# Check NUMA configuration
numactl --hardware

# Run with NUMA awareness (if needed)
numactl --interleave=all python your_script.py
```

### Monitor Resources

```bash
# Monitor CPU, memory, I/O in real-time
htop

# Or use
top

# Monitor I/O specifically
iostat -x 5
```

## Troubleshooting

### RDKit Installation Issues

If conda installation of rdkit fails:

```bash
# Try pip version
pip install rdkit-pypi
```

### Memory Issues

If you encounter memory errors:

```bash
# Check available memory
free -h

# Reduce workers/batch-size:
# --workers 8 --batch-size 25000
```

### Permission Issues

```bash
# Ensure data directories are writable
chmod -R u+w data/

# Check disk space
df -h
```

## Recommended Workflow

1. **Setup Environment** (this guide)
2. **Transfer Data** (if needed)
3. **Test Single Chunk** (validate timing)
4. **Run Batch Processing** (all 28 chunks)
5. **Monitor Progress** (check logs regularly)

## Expected Timeline

- Single chunk (500K rows): ~3-4 hours
- All 28 chunks: ~100 hours total
- Can run multiple chunks in parallel if you have resources

## Support

If you encounter issues:
1. Check logs in output directories
2. Verify environment with test commands above
3. Ensure sufficient disk space and memory
