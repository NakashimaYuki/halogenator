#!/usr/bin/env bash
# -*- coding: ascii -*-
set -euo pipefail
conda activate halo-p0 || true

echo "Starting P0 halogenation demo..."

# 1) 解析名称 -> parents.smi（或用户自备）
echo "Step 1: Ingesting flavonoid names and standardizing..."
python -m halogenator.cli ingest -c configs/p0.yaml

# 2) k=1 单次取代（R1–R5 + QC + 去重）
echo "Step 2: Generating k=1 halogenated products..."
python -m halogenator.cli k1 -c configs/p0.yaml

# 3) 生成摘要报表
echo "Step 3: Generating summary report..."
python -m halogenator.cli report -c configs/p0.yaml

echo "P0 finished. Outputs in data/output/p0/"