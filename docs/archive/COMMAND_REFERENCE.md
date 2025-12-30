# Halogenator 命令参考清单

## 标准枚举命令模板

### 1. Raw 模式 + 层级化输出（推荐用于探索性分析）
```bash
python -m halogenator.cli enum \
  -c configs/<config>.yaml \
  --no-constraints \
  --no-sugar-mask \
  --no-sym-fold \
  --no-dedup \
  --r2-fallback \
  --out-structure hierarchical \
  --group-by family \
  --outdir data/output/<experiment_name>
```

**特点**：
- 无化学约束（产物数最多）
- 无糖环保护（所有位点枚举）
- 无对称性折叠（所有异构体）
- 无去重（保留所有路径）
- R2b 回退启用（flavanone 检测更宽松）
- 层级化输出（易于浏览）
- 家族分组（R2a+R2b → R2）

### 2. Strict 模式 + 层级化输出（推荐用于生产/发布）
```bash
python -m halogenator.cli enum \
  -c configs/<config>.yaml \
  --out-structure hierarchical \
  --group-by family \
  --outdir data/output/<experiment_name>
```

**特点**：
- 完整化学约束（per_ring_quota, min_graph_distance）
- 糖环启发式保护（保留糖苷）
- 对称性折叠（减少冗余）
- InChI 去重（唯一产物）
- R2b 严格检测（默认）
- 层级化输出
- 家族分组

### 3. Flat 输出（快速分析，无 SDF 文件）
```bash
python -m halogenator.cli enum \
  -c configs/<config>.yaml \
  --no-constraints \
  --no-sugar-mask \
  --no-sym-fold \
  --no-dedup \
  --r2-fallback \
  --outdir data/output/<experiment_name>
# 注意：不指定 --out-structure，默认为 flat
```

**特点**：
- 仅生成 products_k2.parquet 和统计 CSV
- 无 SDF 文件（节省磁盘空间）
- 适合纯数据分析

---

## 已验证的实验命令

### Naringenin (Flavanone, R2b Fallback 验证)
```bash
python -m halogenator.cli enum \
  -c configs/naringenin_validation_k2.yaml \
  --no-constraints --no-sugar-mask --no-sym-fold --no-dedup --r2-fallback \
  --out-structure hierarchical --group-by family \
  --outdir data/output/naringenin_k2_raw_hierarchical
```

**预期输出**：
- 992 products (R1=528, R3=348, R2=116)
- R2 全部为 R2b (detection='fallback')
- 132 SDF files

### Flavone (芳香 C 环基准)
```bash
python -m halogenator.cli enum \
  -c configs/flavone_k2.yaml \
  --no-constraints --no-sugar-mask --no-sym-fold --no-dedup --r2-fallback \
  --out-structure hierarchical --group-by family \
  --outdir data/output/flavone_raw_hierarchical
```

**预期输出**：
- 992 products (R1=992, R2=0, R3=0)
- 无 R2 产物（无 sp³ CH₂）
- 132 SDF files

---

## 字段验证命令

### 验证 Parquet 字段完整性
```bash
python scripts/verify_parquet_fields.py <parquet_path>
```

**示例**：
```bash
python scripts/verify_parquet_fields.py \
  data/output/naringenin_k2_raw_hierarchical/products_k2.parquet
```

**预期输出**：
```
==== Parquet Field Verification ====
File: ...
Total products: 992

Field checks:
  Has sub_rule: YES
  Has detection: YES
  Has rule_family: YES

Family counts:
  R1: 528
  R3: 348
  R2: 116

R2 family products: 116
  R2 sub_rule distribution: {'R2b': 116}
  R2 detection distribution: {'fallback': 116}

k == k_ops consistency: YES
```

---

## 测试命令

### 运行 R2b 测试套件
```bash
python tests/mini_r2b_naringenin.py
```

**预期输出**：
```
================================================================================
Minimal R2b Regression Test for Naringenin
================================================================================
[OK] R2b naringenin detection test passed
[OK] R2b fallback test passed
[OK] R2b fallback configuration test passed
[OK] R2b detection field tracking test passed
[OK] By-rule family aggregation test passed
[OK] k>=2 metadata fields test passed

================================================================================
ALL TESTS PASSED
================================================================================
```

---

## 参数说明

### 关键参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--out-structure` | `flat` | 输出结构：`flat`（单 parquet）或 `hierarchical`（层级 SDF） |
| `--group-by` | `rule` | 分组模式：`rule`（R2a, R2b）或 `family`（R2） |
| `--r2-fallback` | 配置决定 | 启用 R2b 回退检测（raw 模式自动启用） |
| `--no-constraints` | False | 禁用所有化学约束 |
| `--no-sugar-mask` | False | 禁用糖环保护 |
| `--no-sym-fold` | False | 禁用对称性折叠 |
| `--no-dedup` | False | 禁用 InChI 去重 |

### Raw 模式开关组合
```bash
--no-constraints --no-sugar-mask --no-sym-fold --no-dedup
```
这4个参数一起使用，进入"完全原始"模式（无任何过滤）。

---

## 输出目录结构

### Flat 输出（`--out-structure flat`）
```
output_dir/
├── products_k2.parquet      # 所有产物
├── by_rule.csv               # 规则统计
├── pivot_*.csv               # 透视表
└── qa_summary.json           # QA 统计
```

### Hierarchical 输出（`--out-structure hierarchical`）
```
output_dir/
├── products_k2.parquet      # Flat 格式备份
├── by_rule.csv               # 家族统计（if --group-by family）
├── pivot_*.csv
├── qa_summary.json
└── <parent_name>/            # 层级化目录
    ├── index.json            # 产物索引
    ├── k1_summary.csv
    ├── k2_summary.csv
    ├── k1/                   # k=1 产物
    │   ├── F/
    │   │   └── <parent>_F.sdf
    │   ├── Cl/
    │   ├── Br/
    │   └── I/
    └── k2/                   # k=2 产物
        ├── F/
        │   ├── <k1_inchikey>/
        │   │   └── <k1_inchikey>_F.sdf
        │   └── ...
        ├── Cl/
        ├── Br/
        └── I/
```

---

## 常见问题

### Q1: 为什么没有生成层级化文件夹？
**A**: 忘记指定 `--out-structure hierarchical`，默认使用 `flat` 模式。

### Q2: by_rule.csv 显示 R2b，但想看 R2 家族统计？
**A**: 添加 `--group-by family` 参数。

### Q3: k=2 的 parquet 缺少 sub_rule 和 detection 字段？
**A**: 确保使用修复后的代码版本（2025-10-28 之后）。

### Q4: R2 产物数为 0，但预期应该有？
**A**: 检查是否启用了 `--r2-fallback`（对 flavanone 结构必需）。

### Q5: 如何验证输出正确性？
**A**: 运行 `python scripts/verify_parquet_fields.py <parquet>` 和 `python tests/mini_r2b_naringenin.py`。

---

*命令参考版本: 1.0*
*最后更新: 2025-10-28*
