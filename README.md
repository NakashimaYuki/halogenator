# Halogenator

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
conda install -c conda-forge rdkit -y
pip install pandas pyarrow pyyaml

# Install package
pip install -e .

# Run P0 demo
bash scripts/run_p0_demo.sh
```

## Usage

```bash
# Ingest flavonoid names and standardize
halogenator ingest -c configs/p0.yaml

# Generate k=1 halogenated products
halogenator k1 -c configs/p0.yaml

# Generate summary report
halogenator report -c configs/p0.yaml
```

## Outputs

- `data/output/p0/parents.smi` - Standardized parent molecules
- `data/output/p0/products_k1.sdf` - Halogenated products with properties
- `data/output/p0/products_k1.parquet` - Product table
- `data/output/p0/summary_k1.csv` - Summary statistics

## ASCII Only

All source code files must be pure ASCII. No Unicode characters are allowed.