# Part B 剩余任务快速执行指南

**本指南提供**: 逐步命令和验收标准，可直接复制执行

---

## 🎯 任务概览

| 任务 | 优先级 | 状态 | 估计时间 |
|------|--------|------|----------|
| B5: 用 v2 重跑四库 | **P1 最高** | 🔲 待执行 | 4-8 小时 |
| B6: 统计+可视化 | P1 | 🔲 待执行 | 2-4 小时 |
| B4: 回归测试 | P1 | 🔲 待执行 | 2-3 小时 |
| B2: Schema 新列 | P2 | 🔲 待执行 | 1-2 小时 |
| B7: 报告整合 | P2 | 🔲 待执行 | 2 小时 |
| B8: VS 导出 | P2 | 🔲 待执行 | 2-3 小时 |

---

## 📦 任务 B5: 用 v2 重跑四个派生库

### 前置检查
```bash
# 1. 检查输入文件是否存在
ls -lh data/output/nplike/Flavone-1X/products.parquet
ls -lh data/output/nplike/Flavone-2X/products.parquet

# 2. 检查 v2 脚本
python scripts/08_transform_library_v2.py --help

# 3. 检查配置文件
cat configs/transforms.yaml | grep -A 10 OH_to_OMe
cat configs/transforms.yaml | grep -A 10 OH_to_NH2
```

### 执行命令

#### Step 1: Flavone-1X-Me
```bash
python scripts/08_transform_library_v2.py apply \
  -i data/output/nplike/Flavone-1X/products.parquet \
  -o data/output/nplike_v2/Flavone-1X-Me/ \
  --xf-config configs/transforms.yaml \
  --xf-name OH_to_OMe \
  --workers 8 \
  --batch-size 100000 \
  2>&1 | tee logs/v2_flavone_1x_me.log
```

**预期输出**:
- 产物数量: ~220k-230k (比 v1 的 ~214k 多)
- 时间: ~1-2 分钟（视机器性能）
- 内存峰值: < 4GB

#### Step 2: Flavone-1X-NH2
```bash
python scripts/08_transform_library_v2.py apply \
  -i data/output/nplike/Flavone-1X/products.parquet \
  -o data/output/nplike_v2/Flavone-1X-NH2/ \
  --xf-config configs/transforms.yaml \
  --xf-name OH_to_NH2 \
  --workers 8 \
  --batch-size 100000 \
  2>&1 | tee logs/v2_flavone_1x_nh2.log
```

#### Step 3: Flavone-2X-Me
```bash
python scripts/08_transform_library_v2.py apply \
  -i data/output/nplike/Flavone-2X/products.parquet \
  -o data/output/nplike_v2/Flavone-2X-Me/ \
  --xf-config configs/transforms.yaml \
  --xf-name OH_to_OMe \
  --workers 8 \
  --batch-size 100000 \
  2>&1 | tee logs/v2_flavone_2x_me.log
```

**预期输出**:
- 产物数量: ~2.5M-3M（比 v1 的 ~2.4M 多）
- 时间: ~30-60 分钟
- 内存峰值: < 6GB

#### Step 4: Flavone-2X-NH2
```bash
python scripts/08_transform_library_v2.py apply \
  -i data/output/nplike/Flavone-2X/products.parquet \
  -o data/output/nplike_v2/Flavone-2X-NH2/ \
  --xf-config configs/transforms.yaml \
  --xf-name OH_to_NH2 \
  --workers 8 \
  --batch-size 100000 \
  2>&1 | tee logs/v2_flavone_2x_nh2.log
```

### 验收检查
```bash
# 检查输出文件
for lib in Flavone-1X-Me Flavone-1X-NH2 Flavone-2X-Me Flavone-2X-NH2; do
  echo "=== $lib ==="
  ls -lh data/output/nplike_v2/$lib/
  python -c "import pandas as pd; df = pd.read_parquet('data/output/nplike_v2/$lib/products.parquet'); print(f'Rows: {len(df):,}')"
  echo ""
done
```

**预期结果**:
```
=== Flavone-1X-Me ===
products.parquet  (~15-20 MB)
dedup.db         (~10-15 MB)
SUMMARY.json
Rows: 220,000-230,000

=== Flavone-1X-NH2 ===
Rows: 220,000-230,000

=== Flavone-2X-Me ===
Rows: 2,500,000-3,000,000

=== Flavone-2X-NH2 ===
Rows: 2,500,000-3,000,000
```

### 故障排查
```bash
# 如果内存不足，降低 batch-size
--batch-size 50000  # 或 20000

# 如果需要续跑（中断后）
--resume

# 查看日志
tail -f logs/v2_flavone_2x_me.log
```

---

## 📊 任务 B6: 生成统计和可视化

### Step 1: 生成统计（05_summaries.py）

```bash
# 确保 05_summaries.py 支持新的 schema（或使用兼容模式）
for lib in Flavone-1X-Me Flavone-1X-NH2 Flavone-2X-Me Flavone-2X-NH2; do
  echo "Generating stats for $lib..."
  python scripts/05_summaries.py \
    -i data/output/nplike_v2/$lib/products.parquet \
    -o data/output/nplike_v2/$lib/ \
    2>&1 | tee logs/stats_$lib.log
done
```

**预期输出文件**:
```
data/output/nplike_v2/Flavone-1X-Me/
  by_rule.csv
  by_rule_family.csv
  halogen_atoms_overall.csv
  halogen_atoms_by_k.csv
  k2_halogen_pairs.csv
  overall_stats.json
```

### Step 2: 批量可视化（10_batch_visualize.py）

**修改配置**:
编辑 `scripts/10_batch_visualize.py`:
```python
# 修改第 33 行
BASE_DIR = Path('data/output/nplike_v2')  # 改为 v2 输出目录
```

**执行**:
```bash
python scripts/10_batch_visualize.py 2>&1 | tee logs/batch_viz_v2.log
```

**预期输出**:
```
data/viz/
  Flavone-1X-Me/
    Flavone-1X-Me_sample_5000.parquet
    Flavone-1X-Me_sample_500.parquet
    Flavone-1X-Me_sample_200.parquet
    Flavone-1X-Me_gallery.html
    Flavone-1X-Me_gallery_thumbs/  (500 PNG files)
    grids/
      page_0001.png
      page_0002.png
    sprite/
      sprite.png
      sprite_index.csv
  Flavone-1X-NH2/
    ...
  Flavone-2X-Me/
    ...
  Flavone-2X-NH2/
    ...

logs/viz/
  Flavone-1X-Me_20251110_*.log
  Flavone-1X-NH2_20251110_*.log
  ...
```

**验收**:
- 打开 HTML 画廊检查是否能正常查看
- 检查 PNG 网格图质量
- 确认采样多样性（观感）

---

## 🧪 任务 B4: 回归测试（v1⊆v2）

### Step 1: 准备测试数据

**选择测试集**（示例：从 Flavone-1X 中采样 1000 个多位点分子）:
```bash
python -c "
import pandas as pd
df = pd.read_parquet('data/output/nplike/Flavone-1X/products.parquet')
# 筛选 k >= 2 的分子（多位点候选）
multi_site = df[df['k'] >= 2].head(1000)
multi_site.to_parquet('data/test/multi_site_sentinel_1k.parquet', index=False)
print(f'Prepared {len(multi_site)} multi-site test molecules')
"
```

### Step 2: 用 v1 和 v2 分别运行

```bash
# v1 (使用原始脚本)
python scripts/08_transform_library.py apply \
  -i data/test/multi_site_sentinel_1k.parquet \
  -o data/test/v1_regression_test/ \
  --xf-config configs/transforms.yaml \
  --xf-name OH_to_OMe

# v2
python scripts/08_transform_library_v2.py apply \
  -i data/test/multi_site_sentinel_1k.parquet \
  -o data/test/v2_regression_test/ \
  --xf-config configs/transforms.yaml \
  --xf-name OH_to_OMe \
  --workers 4 \
  --batch-size 10000
```

### Step 3: 比较结果

```python
import pandas as pd

v1_df = pd.read_parquet('data/test/v1_regression_test/products.parquet')
v2_df = pd.read_parquet('data/test/v2_regression_test/products.parquet')

v1_keys = set(v1_df['inchikey'].dropna())
v2_keys = set(v2_df['inchikey'].dropna())

print(f"V1 products: {len(v1_keys):,}")
print(f"V2 products: {len(v2_keys):,}")
print(f"Common: {len(v1_keys & v2_keys):,}")
print(f"Only in V1: {len(v1_keys - v2_keys):,}")
print(f"Only in V2: {len(v2_keys - v1_keys):,}")

# 断言
assert v1_keys.issubset(v2_keys), "❌ V1 产物未全部包含于 V2！"
print("✅ 回归测试通过: V1 ⊆ V2")
```

### Step 4: 加入 CI（可选）

创建 `tests/test_multi_site_regression.py`:
```python
import pytest
import pandas as pd
from pathlib import Path

def test_v1_v2_regression():
    """Test that V2 contains all V1 products (multi-site molecules)"""
    v1_path = Path('data/test/v1_regression_test/products.parquet')
    v2_path = Path('data/test/v2_regression_test/products.parquet')

    if not v1_path.exists() or not v2_path.exists():
        pytest.skip("Regression test data not found")

    v1_df = pd.read_parquet(v1_path)
    v2_df = pd.read_parquet(v2_path)

    v1_keys = set(v1_df['inchikey'].dropna())
    v2_keys = set(v2_df['inchikey'].dropna())

    # V1 should be subset of V2
    assert v1_keys.issubset(v2_keys), f"V1 has {len(v1_keys - v2_keys)} products not in V2"

    # V2 should have more products (from bug fix)
    assert len(v2_keys) > len(v1_keys), "V2 should have more products than V1"
```

运行测试:
```bash
pytest tests/test_multi_site_regression.py -v
```

---

## 📋 任务 B2: Schema 新列（可选）

**说明**: 此任务为增强功能，非必需。如果时间紧张可跳过。

### 修改位置

1. **修改 output_schema** (`scripts/08_transform_library_v2.py:578-601`):
```python
output_schema = pa.schema([
    # ... 现有列 ...
    ('halogen', pa.string()),
    ('halogen_atom_count', pa.int32()),
    ('halogen_pair', pa.string()),
    # 新增列
    ('halogens_set', pa.string()),           # 如 "F|Cl"
    ('halogen_counts_json', pa.string()),    # 如 '{"F":1,"Cl":1}'
    ('primary_halogen', pa.string()),        # 如 "F" (最后一步引入)
    # ... 其余列 ...
])
```

2. **计算新字段** (在 `TransformationEngineV2.apply_to_molecule` 中):
```python
# 在生成产物时（第 432-445 行附近）
def calculate_halogen_fields(prod_mol, xf_label):
    """计算产物的卤素字段"""
    halogens = {'F': 0, 'Cl': 0, 'Br': 0, 'I': 0}
    for atom in prod_mol.GetAtoms():
        sym = atom.GetSymbol()
        if sym in halogens:
            halogens[sym] += 1

    # 过滤非零卤素
    present_halogens = {k: v for k, v in halogens.items() if v > 0}

    # halogens_set: "F|Cl"
    halogens_set = '|'.join(sorted(present_halogens.keys()))

    # halogen_counts_json: '{"F":1,"Cl":1}'
    halogen_counts_json = json.dumps(present_halogens) if present_halogens else None

    # primary_halogen: 最后一步引入的卤素
    if 'OMe' in xf_label:
        primary_halogen = None  # 不引入卤素
    elif 'NH2' in xf_label:
        primary_halogen = None
    else:
        primary_halogen = list(present_halogens.keys())[0] if present_halogens else None

    return halogens_set, halogen_counts_json, primary_halogen

# 在产物字典中添加
products.append({
    # ... 现有字段 ...
    'halogens_set': halogens_set,
    'halogen_counts_json': halogen_counts_json,
    'primary_halogen': primary_halogen,
    # ...
})
```

---

## 📝 任务 B7: 报告整合

### Step 1: 创建 SCHEMA.json

```bash
cat > SCHEMA.json << 'EOF'
{
  "schema_version": "2.0",
  "description": "Natural Product-Like Library Schema (v2 - Performance Optimized)",
  "last_updated": "2025-11-10",
  "统计口径说明": {
    "产物口径": "每个产物分子按 InChIKey 唯一计数，用于库容量统计",
    "原子口径": "统计每个产物分子中的卤素原子总数（可以是多个同种或不同种卤素），用于卤素分布分析",
    "混合卤素": "如果产物同时含有 F 和 Cl，halogens_set 为 'F|Cl'，halogen_counts_json 为 '{\"F\":1,\"Cl\":1}'"
  },
  "columns": [
    {
      "name": "smiles",
      "type": "string",
      "description": "产物 SMILES（可能非 canonical）"
    },
    {
      "name": "canonical_smiles",
      "type": "string",
      "description": "Canonical isomeric SMILES（用于快速去重）"
    },
    {
      "name": "inchikey",
      "type": "string",
      "description": "InChIKey（唯一标识符）"
    },
    {
      "name": "xf_success",
      "type": "boolean",
      "description": "转化是否成功"
    },
    {
      "name": "xf_label",
      "type": "string",
      "description": "转化标签，如 'OH_to_OMe', 'OH_to_NH2'"
    },
    {
      "name": "xf_site_index",
      "type": "int32",
      "description": "转化位点索引（0-based）"
    },
    {
      "name": "k",
      "type": "int32",
      "description": "卤素化步数（来自基库）"
    },
    {
      "name": "halogen",
      "type": "string",
      "description": "基库卤素类型（单一卤素，如 'F', 'Cl'）"
    },
    {
      "name": "halogens_set",
      "type": "string",
      "description": "产物中的卤素集合，如 'F|Cl'（多种卤素用 | 分隔）"
    },
    {
      "name": "halogen_counts_json",
      "type": "string (JSON)",
      "description": "产物中各卤素的计数，如 '{\"F\":2,\"Cl\":1}'"
    },
    {
      "name": "MW",
      "type": "float64",
      "description": "分子量"
    },
    {
      "name": "TPSA",
      "type": "float64",
      "description": "拓扑极性表面积"
    },
    {
      "name": "cLogP",
      "type": "float64",
      "description": "计算的 logP"
    }
  ]
}
EOF
```

### Step 2: 更新 NPLIKE_LIBRARY_REPORT.md

```bash
# 备份现有报告
cp data/output/nplike/NPLIKE_LIBRARY_REPORT.md \
   data/output/nplike/NPLIKE_LIBRARY_REPORT_v1_backup.md

# 生成新报告（手动或脚本）
```

**新报告结构**:
```markdown
# Natural Product-Like Library Report (v2)

## 概览

| 库名 | 产物数量 (v2) | 产物数量 (v1) | 增幅 | 可视化 |
|------|---------------|---------------|------|--------|
| Flavone-1X-Me | 230,000 | 214,000 | +7.5% | [HTML](../viz/Flavone-1X-Me/Flavone-1X-Me_gallery.html) |
| Flavone-1X-NH2 | 230,000 | 214,000 | +7.5% | [HTML](../viz/Flavone-1X-NH2/Flavone-1X-NH2_gallery.html) |
| Flavone-2X-Me | 2,800,000 | 2,400,000 | +16.7% | [HTML](../viz/Flavone-2X-Me/Flavone-2X-Me_gallery.html) |
| Flavone-2X-NH2 | 2,800,000 | 2,400,000 | +16.7% | [HTML](../viz/Flavone-2X-NH2/Flavone-2X-NH2_gallery.html) |

## 性能对比

| 指标 | v1 | v2 | 提升 |
|------|----|----|------|
| 吞吐量 | 769 mol/s | 2,632 mol/s | +242% |
| 内存占用 (2X) | >10 GB | <6 GB | -40% |
| 产物恢复率 | 基线 | +62% (多位点) | 修复 bug |

## 统计口径说明

见 [SCHEMA.json](../../SCHEMA.json)

## 可视化预览

### Flavone-1X-Me
![Grid Preview](../viz/Flavone-1X-Me/grids/page_0001.png)

[完整 HTML 画廊](../viz/Flavone-1X-Me/Flavone-1X-Me_gallery.html)

...
```

---

## 🔌 任务 B8: VS 导出接口

### 创建导出脚本

```bash
cat > scripts/11_export_for_vs.py << 'EOF'
#!/usr/bin/env python
"""Export molecular libraries for VS (Virtual Screening)"""

import argparse
import pandas as pd
from pathlib import Path

def export_for_vs(input_path, output_dir, format='smi', max_per_file=1_000_000):
    """Export library to VS-compatible format with chunking"""

    df = pd.read_parquet(input_path)
    print(f"Loaded {len(df):,} molecules from {input_path}")

    # Select VS-relevant columns
    vs_cols = ['smiles', 'inchikey', 'MW', 'TPSA', 'HBD', 'HBA', 'cLogP', 'RotB']
    if 'halogens_set' in df.columns:
        vs_cols.append('halogens_set')
    if 'xf_label' in df.columns:
        vs_cols.append('xf_label')

    df_vs = df[vs_cols].dropna(subset=['smiles', 'inchikey'])
    print(f"Filtered to {len(df_vs):,} valid molecules")

    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Split into chunks
    num_chunks = (len(df_vs) + max_per_file - 1) // max_per_file
    print(f"Splitting into {num_chunks} chunks (max {max_per_file:,} per file)")

    for i in range(num_chunks):
        start = i * max_per_file
        end = min((i + 1) * max_per_file, len(df_vs))
        chunk = df_vs.iloc[start:end]

        if format == 'smi':
            output_file = output_dir / f'chunk_{i:04d}.smi'
            chunk.to_csv(output_file, sep='\t', index=False)
        elif format == 'csv':
            output_file = output_dir / f'chunk_{i:04d}.csv'
            chunk.to_csv(output_file, index=False)

        print(f"  Written chunk {i+1}/{num_chunks}: {output_file} ({len(chunk):,} molecules)")

    print(f"\nExport complete: {output_dir}")

def main():
    parser = argparse.ArgumentParser(description='Export libraries for VS')
    parser.add_argument('-i', '--input', required=True, help='Input parquet file')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('--format', choices=['smi', 'csv'], default='smi')
    parser.add_argument('--max-per-file', type=int, default=1_000_000)

    args = parser.parse_args()
    export_for_vs(args.input, args.output, args.format, args.max_per_file)

if __name__ == '__main__':
    main()
EOF

chmod +x scripts/11_export_for_vs.py
```

### 使用示例

```bash
# 导出 Flavone-1X-Me 用于 VS
python scripts/11_export_for_vs.py \
  -i data/output/nplike_v2/Flavone-1X-Me/products.parquet \
  -o data/vs_export/Flavone-1X-Me/ \
  --format smi \
  --max-per-file 100000

# 导出所有库
for lib in Flavone-1X-Me Flavone-1X-NH2 Flavone-2X-Me Flavone-2X-NH2; do
  python scripts/11_export_for_vs.py \
    -i data/output/nplike_v2/$lib/products.parquet \
    -o data/vs_export/$lib/ \
    --format smi \
    --max-per-file 1000000
done
```

### 试跑验证

```bash
# 取第一个 chunk 试跑（假设有 VS 环境）
head -100 data/vs_export/Flavone-1X-Me/chunk_0000.smi > test_sample.smi

# 用 VS 工具验证格式正确性
# （具体命令取决于 VS 软件，如 AutoDock Vina, Glide, etc.）
```

---

## 🎯 总执行时间估算

| 任务 | 时间 | 可并行 |
|------|------|--------|
| B5: 重跑 1X 两库 | 0.5 小时 | ✅ |
| B5: 重跑 2X 两库 | 2 小时 | ✅ |
| B6: 统计 | 0.5 小时 | ✅ |
| B6: 可视化 | 2 小时 | - |
| B4: 回归测试 | 1 小时 | ✅ |
| B2: Schema 新列 | 1 小时 | ✅ |
| B7: 报告 | 1 小时 | - |
| B8: VS 导出 | 0.5 小时 | ✅ |
| **总计** | **~8.5 小时** | **并行后 ~4-5 小时** |

---

## ✅ 最终验收清单

完成后检查以下各项:

### 数据完整性
- [ ] 四库 v2 产物数量 ≥ v1
- [ ] InChIKey 无重复
- [ ] 回归测试通过 (v1 ⊆ v2)

### 性能指标
- [ ] 2X 库内存占用 < 6GB
- [ ] 吞吐量 > 2000 mol/s

### 可视化交付
- [ ] 四库 HTML 画廊可离线查看
- [ ] Grid PNG 清晰
- [ ] Sprite 正常

### 文档完整
- [ ] SCHEMA.json 创建
- [ ] NPLIKE_LIBRARY_REPORT.md 更新
- [ ] VS 导出格式验证

---

**快速指南版本**: 1.0
**最后更新**: 2025-11-10
