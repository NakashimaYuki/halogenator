# Halogenator

![CI](https://github.com/NakashimaYuki/halogenator/actions/workflows/ci.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A baseline halogen substitution system for flavonoid natural products.

## Overview

This project implements rule-based halogen substitution (F, Cl, Br, I) for flavonoids, focusing on k=1 single substitutions in P0 phase.

## Rules

- **R1**: Aromatic sp2 carbon with H → replace H with halogen
- **R2**: C ring carbon with exactly 1 H (not carbonyl) → halogenate
- **R3**: Replace -OH with halogen (excluding carboxylic acid)
- **R4**: Replace -NHx with halogen (excluding amide)
- **R5**: Replace -C(=O)OH with halogen (whole carboxyl group)

## Quick Start

```bash
# Setup environment
conda create -n halo-p0 python=3.10 -y
conda activate halo-p0

# Install RDKit via conda (required - not available via pip on all platforms)
conda install -c conda-forge rdkit -y

# Install other dependencies via pip
pip install pandas pyarrow pyyaml

# Install package (pure Python dependencies only)
pip install -e .

# Run P0 demo
bash scripts/run_p0_demo.sh
```

## Usage

```bash
# Ingest flavonoid names and standardize
halogenator ingest -c configs/p0.yaml

# Generate k=1 halogenated products (default: flavonoids only)
halogenator k1 -c configs/p0.yaml

# Generate k=1 products for rule probes (proves R3/R4/R5 functionality)  
halogenator k1 -c configs/p0.yaml --subset probes

# Generate k=1 products for all molecules (flavonoids + probes)
halogenator k1 -c configs/p0.yaml --subset all

# Generate summary report (default: flavonoids)
halogenator report -c configs/p0.yaml

# Generate report for probes subset
halogenator report -c configs/p0.yaml --subset probes

# Generate report for all molecules
halogenator report -c configs/p0.yaml --subset all

# Use custom output directory (optional)
halogenator k1 -c configs/p0.yaml --outdir tmp/run1
halogenator report -c configs/p0.yaml --outdir tmp/run1
```

## Outputs

### Input Files
- `data/output/p0/parents.smi` - Flavonoid molecules (default subset)
- `data/input/rule_probes.smi` - Rule validation probes (phenol, ethylamine, benzoic_acid)
- `examples/input/parents_flavonoids_10.smi` - **Stable reference**: 10 standard flavonoids for reproducible testing

### Product Files
- `data/output/p0/products_k1.sdf` - Flavonoid products (SDF format)
- `data/output/p0/products_k1.parquet` - Flavonoid products (table format)
- `data/output/p0/products_k1_probes.*` - Probe products (with _probes suffix)
- `data/output/p0/products_k1_all.*` - All molecules (with _all suffix)

### Report Files

**Flavonoids subset (no suffix):**
- `data/output/p0/summary_k1.csv` - Summary statistics  
- `data/output/p0/summary_k1_pivot.csv` - Rule x halogen matrix

**Probes subset (_probes suffix):**
- `data/output/p0/summary_k1_probes.csv` - Summary statistics
- `data/output/p0/summary_k1_probes_pivot.csv` - Rule x halogen matrix

**All molecules subset (_all suffix):**
- `data/output/p0/summary_k1_all.csv` - Summary statistics
- `data/output/p0/summary_k1_all_pivot.csv` - Overall rule x halogen matrix
- `data/output/p0/summary_k1_all_flavonoids_pivot.csv` - Flavonoid-specific matrix
- `data/output/p0/summary_k1_all_probes_pivot.csv` - Probe-specific matrix

### Temporary Files
- `parents_all.smi` - Combined parent molecules file (created automatically for 'all' subset)
  - When using `--outdir`, this file is placed in the specified output directory and can be safely deleted after report generation
  - When no `--outdir` is specified, a temporary directory is used and automatically cleaned up

## Testing

### Unit Tests
Run core unit tests (fast, no external dependencies):

```bash
python -m unittest discover -v
```

### Integration Tests  
Run integration tests (slower, requires file I/O):

```bash
HALO_INTEGRATION=1 python -m unittest tests.test_subset_consistency -v
```

### Code Quality
Check ASCII compliance of all source files:

```bash
./scripts/check_ascii.sh
```

## Known Issues and Version Compatibility

### RDKit Warnings
Some RDKit versions may print non-fatal warnings during molecule processing. These warnings do not affect product generation or statistical results. The system is designed to handle RDKit pre-condition violations gracefully with proper sanitization fallbacks.

### Recommended Versions
The system has been tested with the following dependency versions:
- **RDKit**: 2023.09 or 2024.03 (recommended via conda-forge)
- **Python**: 3.10 or 3.11
- **Pandas**: 2.0+ 
- **PyArrow**: 12.0+

### Installation Notes
- Install RDKit via `conda install -c conda-forge rdkit` (recommended)
- Avoid installing RDKit via pip as it may lack platform-specific optimizations
- Different RDKit versions may produce slightly different warning messages but should generate identical molecular outputs

## ASCII Only

All source code files must be pure ASCII. No Unicode characters are allowed.