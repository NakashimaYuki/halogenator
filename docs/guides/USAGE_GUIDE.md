# 新增脚本使用说明

**版本**: 2025-11-09
**作者**: Task D 性能增强
**适用范围**: 类天然产物库生成与可视化

---

## 目录

1. [09_visualize_library.py - 分子可视化工具](#1-09_visualize_librarypy---分子可视化工具)
2. [08_transform_library_v2.py - 性能优化版转化引擎](#2-08_transform_library_v2py---性能优化版转化引擎)
3. [08a_fill_descriptors.py - 分子描述符填充工具](#3-08a_fill_descriptorspy---分子描述符填充工具)
4. [10_batch_visualize.py - 批量可视化脚本](#4-10_batch_visualizepy---批量可视化脚本)
5. [常见问题与最佳实践](#5-常见问题与最佳实践)

---

## 1. 09_visualize_library.py - 分子可视化工具

### 功能概述

用于大规模分子库的可视化、采样和展示，支持 4 种子命令：
- **sample**: 从大库中抽取代表性子集
- **grid**: 生成 PNG/SVG 网格图像
- **html**: 创建离线可交互的 HTML 画廊
- **sprite**: 生成雪碧图（sprite sheet）

### 系统要求

```bash
pip install rdkit pandas pyarrow pillow numpy
```

---

### 1.1 sample - 采样子命令

#### 功能
从大型分子库中提取代表性样本，支持 3 种采样策略。

#### 基本语法
```bash
python scripts/09_visualize_library.py sample \
  -i <输入文件.parquet> \
  -o <输出文件.parquet> \
  --n <样本数量> \
  --strategy <采样策略> \
  [可选参数]
```

#### 参数说明

| 参数 | 必需 | 说明 | 默认值 |
|------|------|------|--------|
| `-i, --input` | 是 | 输入 Parquet/SMI 文件路径 | - |
| `-o, --output` | 是 | 输出 Parquet/CSV 文件路径 | - |
| `--n` | 否 | 采样数量 | 5000 |
| `--strategy` | 否 | 采样策略：`random`/`stratified`/`diverse` | `random` |
| `--strata-cols` | 否 | 分层采样列（逗号分隔），如 `k,xf_label` | `k` |
| `--diverse-fp` | 否 | 多样性采样指纹类型：`ecfp4`/`mhfp6` | `ecfp4` |
| `--diverse-thresh` | 否 | 多样性采样相似度阈值 (0-1) | 0.4 |
| `--seed` | 否 | 随机种子（可重复性） | 2025 |

#### 采样策略详解

**1. random - 随机采样**
```bash
python scripts/09_visualize_library.py sample \
  -i data/output/nplike/Flavone-2X-Me/products.parquet \
  -o data/viz/Flavone-2X-Me_random_5k.parquet \
  --n 5000 \
  --strategy random \
  --seed 2025
```
- **适用场景**: 快速预览、初步 QC
- **优点**: 速度快（<1 秒）
- **缺点**: 可能遗漏罕见结构类型

**2. stratified - 分层采样**
```bash
python scripts/09_visualize_library.py sample \
  -i data/output/nplike/Flavone-2X-Me/products.parquet \
  -o data/viz/Flavone-2X-Me_stratified_5k.parquet \
  --n 5000 \
  --strategy stratified \
  --strata-cols k,xf_label,halogen_pair \
  --seed 2025
```
- **适用场景**: 需要各类别均衡代表
- **优点**: 保证每个组合（k=1/2, 卤素类型, 转化类型）都有样本
- **缺点**: 若某分层样本数 < n/层数，该层全部抽取

**3. diverse - 多样性采样**
```bash
python scripts/09_visualize_library.py sample \
  -i data/output/nplike/Flavone-2X-Me/products.parquet \
  -o data/viz/Flavone-2X-Me_diverse_5k.parquet \
  --n 5000 \
  --strategy diverse \
  --diverse-fp ecfp4 \
  --diverse-thresh 0.5 \
  --seed 2025
```
- **适用场景**: 最大化化学多样性、结构覆盖
- **原理**: ECFP4 指纹 + Sphere Exclusion（Tanimoto 相似度阈值）
- **优点**: 确保样本化学空间分布广
- **缺点**: 耗时较长（5000 样本约 10-30 秒）

#### 输出示例
```
[23:08:02] INFO - Reading data/output/nplike/Flavone-1X-Me/products.parquet...
[23:08:02] INFO - Loaded 213902 molecules
[23:08:02] INFO - Random sampling 1000 molecules from 213902...
[23:08:02] INFO - Sampled 1000 molecules
[23:08:02] INFO - Saved to data/viz/Flavone-1X-Me_sample_1000.parquet
```

---

### 1.2 grid - 网格图像子命令

#### 功能
生成 PNG 网格图像，每页包含多个分子结构，适合打印或 PDF 导出。

#### 基本语法
```bash
python scripts/09_visualize_library.py grid \
  -i <输入文件.parquet> \
  -o <输出目录/> \
  [可选参数]
```

#### 参数说明

| 参数 | 必需 | 说明 | 默认值 |
|------|------|------|--------|
| `-i, --input` | 是 | 输入 Parquet 文件 | - |
| `-o, --output` | 是 | 输出目录（自动创建） | - |
| `--per-page` | 否 | 每页分子数 | 200 |
| `--cols` | 否 | 每页列数 | 10 |
| `--size` | 否 | 单个分子图像尺寸（宽 高） | 250 250 |
| `--highlight` | 否 | 高亮模式：`none`/`halogens`/`xf_site` | `halogens` |
| `--title` | 否 | 网格标题 | `Molecular Library Grid` |

#### 使用示例

**示例 1: 标准网格（10x10，卤素高亮）**
```bash
python scripts/09_visualize_library.py grid \
  -i data/viz/Flavone-1X-Me_sample_1000.parquet \
  -o data/viz/grids/ \
  --per-page 100 \
  --cols 10 \
  --size 200 200 \
  --highlight halogens \
  --title "Flavone-1X-Me Sample (n=1000)"
```

**示例 2: 高分辨率网格（大图，用于论文插图）**
```bash
python scripts/09_visualize_library.py grid \
  -i data/viz/sample.parquet \
  -o data/viz/grids_hires/ \
  --per-page 50 \
  --cols 5 \
  --size 400 400 \
  --highlight xf_site
```

#### 输出结构
```
data/viz/grids/
├── page_0001.png  (100 分子, 10x10 网格)
├── page_0002.png
├── ...
└── draw_errors.csv  (若有解析失败的分子)
```

#### 图例说明
每个分子下方自动生成图例，包含：
- InChIKey 前 14 位
- k 值（卤代数）
- 转化标签（如 `OH->OMe`）
- 分子量 MW

---

### 1.3 html - HTML 画廊子命令

#### 功能
生成离线可交互的 HTML 画廊，支持搜索、排序、筛选（无需服务器）。

#### 基本语法
```bash
python scripts/09_visualize_library.py html \
  -i <输入文件.parquet> \
  -o <输出文件.html> \
  [可选参数]
```

#### 参数说明

| 参数 | 必需 | 说明 | 默认值 |
|------|------|------|--------|
| `-i, --input` | 是 | 输入 Parquet 文件 | - |
| `-o, --output` | 是 | 输出 HTML 文件路径 | - |
| `--thumb-size` | 否 | 缩略图尺寸（宽 高） | 200 200 |
| `--columns` | 否 | 显示列（逗号分隔） | `inchikey,k,xf_label,MW,TPSA,cLogP` |
| `--highlight` | 否 | 高亮模式 | `halogens` |
| `--title` | 否 | 画廊标题 | `Molecular Library Gallery` |
| `--workers` | 否 | 并行渲染进程数 | 4 |

#### 使用示例

**示例 1: 标准画廊（500 分子）**
```bash
python scripts/09_visualize_library.py html \
  -i data/viz/Flavone-1X-Me_sample_500.parquet \
  -o data/viz/gallery/Flavone-1X-Me.html \
  --thumb-size 180 180 \
  --columns inchikey,k,xf_label,xf_site_index,MW,TPSA,HBD,HBA,cLogP \
  --highlight halogens \
  --title "Flavone-1X-Me Gallery (n=500)" \
  --workers 6
```

**示例 2: 快速预览（200 分子，基础列）**
```bash
python scripts/09_visualize_library.py html \
  -i data/viz/sample_200.parquet \
  -o data/viz/preview.html \
  --columns inchikey,MW,cLogP \
  --workers 4
```

#### 输出特性
- **自包含**: 单个 HTML 文件（SVG 内嵌 base64），可直接分享
- **可交互**: 点击列头排序、搜索框实时筛选
- **轻量化**: 500 分子 ≈ 5-10 MB
- **兼容性**: 任何现代浏览器，无需网络

#### 使用建议
- 推荐样本量: 200-1000（超过 1000 加载变慢）
- 大库可视化: 先用 `sample` 抽取，再生成 HTML
- 分享场景: 发送给导师/合作者审阅、组会展示

---

### 1.4 sprite - 雪碧图子命令

#### 功能
生成单张大图（sprite sheet）+ 索引 CSV，适合 Web 仪表盘或快速加载场景。

#### 基本语法
```bash
python scripts/09_visualize_library.py sprite \
  -i <输入文件.parquet> \
  -o <输出目录/> \
  [可选参数]
```

#### 参数说明

| 参数 | 必需 | 说明 | 默认值 |
|------|------|------|--------|
| `-i, --input` | 是 | 输入 Parquet 文件 | - |
| `-o, --output` | 是 | 输出目录 | - |
| `--thumb` | 否 | 缩略图尺寸（宽 高） | 128 128 |
| `--cols` | 否 | 雪碧图列数 | 64 |
| `--workers` | 否 | 并行渲染进程数 | 4 |

#### 使用示例

```bash
python scripts/09_visualize_library.py sprite \
  -i data/viz/Flavone-1X-Me_sample_200.parquet \
  -o data/viz/sprite/ \
  --thumb 100 100 \
  --cols 20 \
  --workers 6
```

#### 输出结构
```
data/viz/sprite/
├── sprite.png         (2000x1000 像素大图，包含全部分子)
└── sprite_index.csv   (行列坐标 + 元数据映射)
```

#### sprite_index.csv 格式
```csv
index,row,col,x,y,inchikey,smiles
0,0,0,0,0,ABCDEFGHIJKLMN-OPQRST-U,COc1c(F)c(C)cc2oc(C)cc(=O)c12
1,0,1,100,0,BCDEFGHIJKLMNO-PQRSTU-V,COc1c(Cl)c(C)cc2oc(C)cc(=O)c12
...
```

#### 使用场景
- **Web 仪表盘**: 单次加载全部缩略图（减少 HTTP 请求）
- **大屏展示**: 瀑布流/网格快速预览
- **外部工具集成**: 通过 CSV 索引定位分子位置

---

## 2. 08_transform_library_v2.py - 性能优化版转化引擎

### 功能概述

**关键改进**（相对于 v1）:
- ✅ **2.36x 更快**: 10k 分子从 13s → 5.5s
- ✅ **修复关键 bug**: 恢复 62% 被遗漏的产物（多位点转化）
- ✅ **真并行**: ProcessPoolExecutor 多核利用
- ✅ **流式写出**: 恒定内存占用
- ✅ **可恢复执行**: `--resume` 断点续跑

### 系统要求
```bash
pip install rdkit pandas pyarrow pyyaml
```

---

### 2.1 apply - 转化应用子命令

#### 基本语法
```bash
python scripts/08_transform_library_v2.py apply \
  -i <输入库.parquet> \
  -o <输出目录/> \
  --xf-config <转化配置.yaml> \
  --xf-name <转化名称> \
  [可选参数]
```

#### 参数说明

| 参数 | 必需 | 说明 | 默认值 |
|------|------|------|--------|
| `-i, --input` | 是 | 输入基础库（Parquet） | - |
| `-o, --outdir` | 是 | 输出目录 | - |
| `--xf-config` | 是 | 转化配置 YAML 文件 | - |
| `--xf-name` | 是 | 转化规则名称（在 YAML 中定义） | - |
| `--batch-size` | 否 | 批处理大小 | 50000 |
| `--workers` | 否 | 并行进程数（建议 = CPU 核数） | 6 |
| `--resume` | 否 | 断点续跑模式（加载已有去重库） | False |

#### 使用示例

**示例 1: 标准转化（OH → OMe）**
```bash
python scripts/08_transform_library_v2.py apply \
  -i data/output/nplike/Flavone-1X/base.parquet \
  -o data/output/nplike/Flavone-1X-Me_v2 \
  --xf-config configs/transforms.yaml \
  --xf-name OH_to_OMe \
  --batch-size 50000 \
  --workers 6
```

**示例 2: 恢复中断任务**
```bash
# 假设之前任务中断在 60% 处
python scripts/08_transform_library_v2.py apply \
  -i data/output/nplike/Flavone-2X/base.parquet \
  -o data/output/nplike/Flavone-2X-Me_v2 \
  --xf-config configs/transforms.yaml \
  --xf-name OH_to_OMe \
  --workers 8 \
  --resume  # 加载 dedup.db 已有键，跳过重复
```

**示例 3: 大库优化（4.86M 分子）**
```bash
python scripts/08_transform_library_v2.py apply \
  -i data/output/nplike/Flavone-2X/base.parquet \
  -o data/output/nplike/Flavone-2X-NH2_v2 \
  --xf-config configs/transforms.yaml \
  --xf-name OH_to_NH2 \
  --batch-size 100000 \  # 加大批处理（减少 DB 提交）
  --workers 12           # 充分利用多核
```

#### 输出结构
```
data/output/nplike/Flavone-1X-Me_v2/
├── products.parquet     (主产物库)
├── dedup.db            (SQLite 去重数据库)
└── SUMMARY.json        (统计摘要)
```

#### SUMMARY.json 示例
```json
{
  "transform": "OH_to_OMe",
  "source_file": "data/output/nplike/Flavone-1X/base.parquet",
  "source_subset": "Flavone-1X",
  "total_input": 241966,
  "total_processed": 348952,
  "total_products": 348952,
  "unique_products": 213902,
  "elapsed_seconds": 132.5,
  "throughput_mol_per_sec": 2634.2,
  "stats": {
    "success": 280146,
    "fail_no_matching_sites": 68806
  }
}
```

#### 输出日志示例
```
2025-11-09 23:14:52 | INFO     | ================================================================================
2025-11-09 23:14:52 | INFO     | APPLY v2: Performance-Optimized Transformation
2025-11-09 23:14:52 | INFO     | ================================================================================
2025-11-09 23:14:52 | INFO     | [Batch 5/5] Processed: 14,442 | Products: 14,442 | Unique: 11,304 | Rate: 2689.0 mol/s
2025-11-09 23:14:52 | INFO     | TRANSFORMATION COMPLETE
2025-11-09 23:14:52 | INFO     |   Time: 5.5s
2025-11-09 23:14:52 | INFO     |   Throughput: 2632.4 mol/s
2025-11-09 23:14:52 | INFO     |   Uniqueness: 78.3%
```

---

### 2.2 性能调优指南

#### 批处理大小（--batch-size）
- **小库（<100k）**: 10000-20000
- **中库（100k-1M）**: 50000
- **大库（>1M）**: 100000-200000
- **原则**: 批越大，DB 提交次数越少（更快），但内存占用略高

#### 并行进程数（--workers）
- **推荐值**: CPU 核数 - 1（留 1 核给主进程）
- **测试**: `--workers 4/6/8/12` 对比吞吐量
- **瓶颈**: 超过核数收益递减（进程切换开销）

#### 断点续跑（--resume）
- **适用场景**: 大库处理中断、增量更新
- **工作原理**: 读取 `dedup.db` 中已有 InChIKey，跳过重复产物
- **注意**: 需保留原输出目录（含 `dedup.db`）

---

### 2.3 与 v1 对比

| 特性 | v1 (原版) | v2 (优化版) |
|------|-----------|-------------|
| **速度** | 769 mol/s | 2632 mol/s (**+242%**) |
| **产物完整性** | ❌ 仅 site_0 | ✅ 全部位点 (site_0/1/2) |
| **并行** | ❌ 单线程 | ✅ 多进程 |
| **内存** | 线性增长 | 恒定 |
| **可恢复** | ❌ 不支持 | ✅ --resume |
| **日志** | 基础 | 实时吞吐/去重率 |

---

## 3. 08a_fill_descriptors.py - 分子描述符填充工具

### 功能概述

为快速枚举模式（v2）生成的产物库补充完整分子描述符。

**使用场景**:
1. v2 转化仅计算基础属性（MW/TPSA/HBD/HBA/cLogP/RotB）
2. 需要完整描述符（20+ 属性）用于 ADMET 筛选
3. 延迟计算策略：先去重，再算描述符（避免浪费）

### 基本语法
```bash
python scripts/08a_fill_descriptors.py \
  -i <输入库.parquet> \
  -o <输出库.parquet> \
  --mode <quick|full> \
  [可选参数]
```

### 参数说明

| 参数 | 必需 | 说明 | 默认值 |
|------|------|------|--------|
| `-i, --input` | 是 | 输入 Parquet 文件 | - |
| `-o, --output` | 是 | 输出 Parquet 文件 | - |
| `--mode` | 否 | 描述符模式：`quick`/`full` | `quick` |
| `--workers` | 否 | 并行进程数 | 4 |
| `--batch-size` | 否 | 批处理大小 | 10000 |

### 描述符模式

#### quick 模式（6 个基础属性）
- MW (分子量)
- TPSA (拓扑极性表面积)
- HBD (氢键供体)
- HBA (氢键受体)
- cLogP (脂水分配系数)
- RotB (可旋转键)

#### full 模式（20+ 个扩展属性）
在 quick 基础上增加:
- FracCSP3 (sp3 碳占比)
- NumRings / NumAromaticRings / NumAliphaticRings / NumSaturatedRings
- NumHeteroatoms / NumHeavyAtoms / NumBonds / NumAromaticBonds
- MolMR (分子折射率)
- LabuteASA (Labute 可及表面积)
- 等...

### 使用示例

**示例 1: 补充基础属性（推荐）**
```bash
python scripts/08a_fill_descriptors.py \
  -i data/output/nplike/Flavone-1X-Me_v2/products.parquet \
  -o data/output/nplike/Flavone-1X-Me_v2/products_full.parquet \
  --mode quick \
  --workers 6 \
  --batch-size 10000
```

**示例 2: 完整描述符（ADMET 筛选用）**
```bash
python scripts/08a_fill_descriptors.py \
  -i data/output/nplike/Flavone-2X-Me_v2/products.parquet \
  -o data/output/nplike/Flavone-2X-Me_v2/products_admet.parquet \
  --mode full \
  --workers 8
```

### 输出日志示例
```
2025-11-09 23:20:15 | INFO     | ================================================================================
2025-11-09 23:20:15 | INFO     | Fill Descriptors: Post-Processing Pipeline
2025-11-09 23:20:15 | INFO     | ================================================================================
2025-11-09 23:20:15 | INFO     | Column status:
2025-11-09 23:20:15 | INFO     |   inchikey: 0 missing (0.0%)
2025-11-09 23:20:15 | INFO     |   MW: 11304 missing (100.0%)  # 需要填充
2025-11-09 23:20:15 | INFO     |   Processed: 11304/11304 (1852.3 mol/s)
2025-11-09 23:20:15 | INFO     | COMPLETE
```

---

## 4. 10_batch_visualize.py - 批量可视化脚本

### 功能概述

自动化为多个库生成完整可视化报告（采样 + 网格 + HTML + 雪碧图）。

### 基本语法
```bash
python scripts/10_batch_visualize.py
```

### 执行流程

对每个库自动执行：
1. **Diverse 采样** (5000 分子)
2. **HTML 采样** (500 分子，从步骤 1 子集抽取)
3. **生成 HTML 画廊**
4. **网格采样** (200 分子)
5. **生成网格图像**
6. **生成雪碧图**

### 默认处理库列表
```python
LIBRARIES = [
    'Flavone-1X-Me',
    'Flavone-1X-NH2',
    'Flavone-2X-Me',
    'Flavone-2X-NH2'
]
```

### 输出结构
```
data/viz/
├── Flavone-1X-Me/
│   ├── Flavone-1X-Me_sample_5000.parquet
│   ├── Flavone-1X-Me_sample_500.parquet
│   ├── Flavone-1X-Me_sample_200.parquet
│   ├── Flavone-1X-Me_gallery.html
│   ├── grids/
│   │   ├── page_0001.png
│   │   └── page_0002.png
│   └── sprite/
│       ├── sprite.png
│       └── sprite_index.csv
├── Flavone-1X-NH2/
│   └── ...
├── Flavone-2X-Me/
│   └── ...
└── Flavone-2X-NH2/
    └── ...
```

### 自定义修改

编辑 `scripts/10_batch_visualize.py` 中的参数：
```python
# 修改采样数量
sample_5k = lib_viz_dir / f'{lib_name}_sample_5000.parquet'  # 改为 10000
sample_500 = lib_viz_dir / f'{lib_name}_sample_500.parquet'  # 改为 1000

# 修改网格参数
'--per-page', '100',  # 改为 200
'--cols', '10',       # 改为 20

# 修改并行数
'--workers', '6',     # 改为 8
```

---

## 5. 常见问题与最佳实践

### Q1: 哪个脚本处理大库最快？

**A**: 按优先级使用：
1. **v2** (`08_transform_library_v2.py`) - 转化枚举阶段
2. **08a** (`08a_fill_descriptors.py`) - 仅对去重后产物计算描述符
3. **09** (`09_visualize_library.py`) - 仅对小样本可视化

**工作流**:
```bash
# Step 1: 快速转化（v2）
python scripts/08_transform_library_v2.py apply -i base.parquet -o products/ --workers 8

# Step 2: 延迟填充描述符（仅唯一产物）
python scripts/08a_fill_descriptors.py -i products/products.parquet -o products_full.parquet --mode quick

# Step 3: 采样可视化（5000 样本）
python scripts/09_visualize_library.py sample -i products_full.parquet -o sample.parquet --n 5000 --strategy diverse
python scripts/09_visualize_library.py html -i sample.parquet -o gallery.html
```

---

### Q2: 如何选择采样策略？

| 目的 | 推荐策略 | 原因 |
|------|---------|------|
| 快速预览 | `random` | 速度最快（<1s） |
| 均衡代表 | `stratified` | 每个类别都有样本 |
| 最大多样性 | `diverse` | 化学空间覆盖广 |
| 论文插图 | `diverse` | 避免结构冗余 |
| 手动审阅 | `stratified` | 便于按类别检查 |

---

### Q3: 可视化样本数量建议？

| 输出格式 | 推荐样本数 | 性能 |
|----------|-----------|------|
| **grid** | 200-1000 | 快（1-5 秒） |
| **html** | 200-500 | 中（2-10 秒） |
| **sprite** | 100-500 | 快（1-3 秒） |
| **sample** (diverse) | 1000-5000 | 慢（10-60 秒） |

**注意**: HTML 超过 1000 分子加载会变慢（浏览器内存压力）

---

### Q4: v2 转化后产物数量比预期少，怎么排查？

**可能原因**:
1. **糖位点过滤**: 检查 `--xf-config` 中 `exclude_sugar_like: true`
2. **最大位点限制**: 检查 `max_sites_per_mol: 4`
3. **无匹配位点**: 查看 `SUMMARY.json` 中 `fail_no_matching_sites` 数量

**排查步骤**:
```bash
# 1. 查看统计摘要
cat data/output/nplike/Flavone-1X-Me_v2/SUMMARY.json | grep fail

# 2. 检查配置
cat configs/transforms.yaml | grep -A 3 defaults

# 3. 测试小样本（关闭糖过滤）
# 修改 transforms.yaml: exclude_sugar_like: false
python scripts/08_transform_library_v2.py apply -i test_100.parquet -o test_output/ --xf-config configs/transforms.yaml --xf-name OH_to_OMe
```

---

### Q5: 如何并行处理多个转化任务？

**方法 1: 后台运行（推荐）**
```bash
# 启动 4 个转化任务（不同库/不同转化）
nohup python scripts/08_transform_library_v2.py apply -i Flavone-1X/base.parquet -o Flavone-1X-Me/ --xf-config configs/transforms.yaml --xf-name OH_to_OMe --workers 4 > log_1x_me.txt 2>&1 &

nohup python scripts/08_transform_library_v2.py apply -i Flavone-1X/base.parquet -o Flavone-1X-NH2/ --xf-config configs/transforms.yaml --xf-name OH_to_NH2 --workers 4 > log_1x_nh2.txt 2>&1 &

# 监控进度
tail -f log_1x_me.txt
```

**方法 2: GNU Parallel**
```bash
# 创建任务列表
cat > tasks.txt <<EOF
-i Flavone-1X/base.parquet -o Flavone-1X-Me/ --xf-name OH_to_OMe
-i Flavone-1X/base.parquet -o Flavone-1X-NH2/ --xf-name OH_to_NH2
-i Flavone-2X/base.parquet -o Flavone-2X-Me/ --xf-name OH_to_OMe
-i Flavone-2X/base.parquet -o Flavone-2X-NH2/ --xf-name OH_to_NH2
EOF

# 并行执行（每个任务 4 workers，最多同时 2 任务）
parallel -j 2 python scripts/08_transform_library_v2.py apply {} --xf-config configs/transforms.yaml --workers 4 :::: tasks.txt
```

---

### Q6: 输出的 HTML 画廊能否嵌入到网页/PPT？

**A**: 可以，但需注意：

**嵌入网页（iframe）**:
```html
<iframe src="gallery.html" width="100%" height="800px" frameborder="0"></iframe>
```

**嵌入 PPT**:
1. 另存为 PDF（浏览器 → 打印 → 保存为 PDF）
2. 插入 PDF 到 PPT（推荐）

**分享方式**:
- **最佳**: 直接发送 `.html` 文件（单文件，无需服务器）
- **备选**: 上传到内网服务器，共享链接

---

### Q7: 如何从雪碧图中提取单个分子？

**方法 1: 使用索引 CSV + Python**
```python
import pandas as pd
from PIL import Image

# 读取雪碧图和索引
sprite = Image.open('data/viz/sprite/sprite.png')
index = pd.read_csv('data/viz/sprite/sprite_index.csv')

# 提取第 10 个分子（index=9）
row = index.iloc[9]
x, y = row['x'], row['y']
thumb_width, thumb_height = 128, 128

mol_img = sprite.crop((x, y, x + thumb_width, y + thumb_height))
mol_img.save(f"molecule_{row['inchikey'][:14]}.png")
```

**方法 2: 使用 ImageMagick**
```bash
# 提取坐标为 (200, 100) 的 128x128 缩略图
convert sprite.png -crop 128x128+200+100 molecule.png
```

---

### 最佳实践总结

#### ✅ 推荐工作流

```bash
# 1. 转化枚举（v2 - 快速）
python scripts/08_transform_library_v2.py apply \
  -i base.parquet -o products_v2/ \
  --xf-config configs/transforms.yaml --xf-name OH_to_OMe \
  --workers 8 --batch-size 100000

# 2. QC 采样（diverse - 高质量）
python scripts/09_visualize_library.py sample \
  -i products_v2/products.parquet -o products_sample_5k.parquet \
  --n 5000 --strategy diverse --diverse-thresh 0.5

# 3. 快速 HTML 画廊（200 样本）
python scripts/09_visualize_library.py sample \
  -i products_sample_5k.parquet -o products_sample_200.parquet \
  --n 200 --strategy random

python scripts/09_visualize_library.py html \
  -i products_sample_200.parquet -o products_gallery.html \
  --workers 6

# 4. （可选）完整描述符（仅需要时）
python scripts/08a_fill_descriptors.py \
  -i products_v2/products.parquet -o products_full.parquet \
  --mode full --workers 6
```

#### ❌ 避免的操作

- ❌ 直接对 4.86M 库生成 HTML（会卡死浏览器）
- ❌ 在 v2 转化时计算全部 20+ 描述符（浪费）
- ❌ 不做采样直接网格化大库（生成几千张图片）
- ❌ 使用 v1 脚本（慢且有 bug）

---

## 附录：参数速查表

### 09_visualize_library.py

| 子命令 | 关键参数 | 默认值 | 说明 |
|--------|---------|--------|------|
| `sample` | `--n` | 5000 | 样本数量 |
| | `--strategy` | random | random/stratified/diverse |
| | `--diverse-thresh` | 0.4 | 多样性阈值 (0-1) |
| `grid` | `--per-page` | 200 | 每页分子数 |
| | `--cols` | 10 | 列数 |
| | `--size` | 250 250 | 单分子尺寸 |
| `html` | `--thumb-size` | 200 200 | 缩略图尺寸 |
| | `--columns` | inchikey,k,... | 显示列 |
| | `--workers` | 4 | 并行数 |
| `sprite` | `--thumb` | 128 128 | 缩略图尺寸 |
| | `--cols` | 64 | 雪碧图列数 |

### 08_transform_library_v2.py

| 参数 | 默认值 | 推荐范围 | 说明 |
|------|-------|---------|------|
| `--batch-size` | 50000 | 10k-200k | 批处理大小 |
| `--workers` | 6 | CPU 核数 - 1 | 并行进程数 |
| `--resume` | False | - | 断点续跑 |

### 08a_fill_descriptors.py

| 参数 | 默认值 | 选项 | 说明 |
|------|-------|------|------|
| `--mode` | quick | quick/full | 描述符模式 |
| `--workers` | 4 | 2-12 | 并行进程数 |
| `--batch-size` | 10000 | 5k-50k | 批处理大小 |

---

## 技术支持

- **问题反馈**: 查看 `PERFORMANCE_REPORT.md` 第 X 节
- **进阶用例**: 参考 `IMPLEMENTATION_SUMMARY.md`
- **源码注释**: 所有脚本包含详细 docstring

---

**版本历史**:
- 2025-11-09: 初版发布
- 后续更新请查阅 Git 提交记录

**维护者**: Task D 开发团队
