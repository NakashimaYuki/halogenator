# Session Handoff: Batch Pipeline Debugging & Repository Cleanup

**日期:** 2025-12-29
**状态:** 🟡 Pipeline部分运行，需调查深层问题并完成清理
**上一会话:** 参数优化与批处理pipeline实施
**本会话核心任务:**
1. 调查batch pipeline timeout的根本原因
2. 完整清理和整理repository
3. 准备Linux迁移

---

## 执行摘要

### 已完成的重要工作

1. ✅ **参数优化Grid Search** - 完成4个配置测试，识别最优配置
2. ✅ **批处理Pipeline设计与实现** - 创建自动化分批处理系统
3. ✅ **部分数据处理** - 成功完成2/14 chunks

### 当前状态

**Pipeline运行情况:**
- Chunk 0: ✅ 完成 (44分钟, 274MB输出)
- Chunk 1: ✅ 完成 (73分钟, 2.47M products, 190MB)
- Chunk 2: ❌ Timeout失败 (>4小时, 7.31M products, 671MB **部分输出**)
- Chunks 3-13: ⏸️ 待处理

**关键问题:**
- Chunk处理时间差异巨大（44分钟 vs >4小时）
- 可能存在比timeout更深层的系统性问题
- Repository需要整理以便迁移到Linux

---

## 第一部分：已完成任务总结

### 1.1 参数优化Grid Search

**文件:** `E:\Projects\halogenator\optimize_parameters.py`

**测试结果:** (150K rows测试集)

| Config | Workers | Max-In-Flight | Batch | Throughput | Memory | Result |
|--------|---------|---------------|-------|------------|--------|--------|
| Test 1 | 16 | 6 | 50K | 672.5 mol/s | 86.7% | ❌ 内存超限 |
| Test 2 | 16 | 8 | 50K | 661.1 mol/s | 85.7% | ❌ 内存超限 |
| Test 3 | 16 | 6 | 25K | 798.9 mol/s | 76.6% | ⚠️ 接近限制 |
| Test 4 | 12 | 8 | 50K | 657.7 mol/s | 43.8% | ✅ **最安全** |

**关键发现:**
- 150K规模测试全部成功
- 500K+规模测试**全部失败**（workers崩溃）
- 结论：需要分批处理策略

**输出文件:**
- `optimization_results_20251228_041629.csv` - 详细测试数据

### 1.2 批处理Pipeline实现

**文件:** `E:\Projects\halogenator\batch_transform_pipeline.py`

**设计特性:**
- ✅ Checkpoint功能 - 断点续传
- ✅ 自动状态追踪 - `pipeline_state.json`
- ✅ 失败恢复 - 只重跑失败chunks
- ✅ 自动结果合并

**配置参数:**
```python
input_file: polyphenol-2X/products.parquet (13.79M rows)
chunk_size: 1,000,000 rows → 14 chunks
workers: 16
max_in_flight: 6
batch_size: 50,000
timeout: 14400s (4小时) ← 当前设置
```

**输出位置:**
- 主目录: `data/output/transforms/polyphenol-2X_BATCHED/`
- Chunks: `data/output/transforms/polyphenol-2X_BATCHED/chunks/chunk_XXX_output/`
- 状态文件: `data/output/transforms/polyphenol-2X_BATCHED/pipeline_state.json`

### 1.3 实际执行结果

**Chunk 0详情:**
```json
{
  "elapsed_seconds": 2652 (44.2 min),
  "output_file": "274 MB",
  "throughput": 1476.7 mol/s,
  "unique_products": 未记录（JSON解析问题）
}
```

**Chunk 1详情:**
```json
{
  "elapsed_seconds": 4382 (73.0 min),
  "output_file": "190 MB",
  "throughput": 994.7 mol/s,
  "unique_products": 2,465,730
}
```

**Chunk 2详情（关键异常）:**
```
Status: timeout failed
Actual output: 671 MB (最大！)
Flush count: #3650
Products written: 7,307,108 (是chunk 1的3倍！)
Last activity: 06:40 (被4小时timeout终止)
Memory: 45.7% (稳定安全)
Missing: SUMMARY.json (未正常完成)
```

---

## 第二部分：关键问题的深度分析

### 2.1 Timeout问题的表面现象

**现象:**
- Chunk 0: 44分钟完成
- Chunk 1: 73分钟完成（慢65%）
- Chunk 2: >4小时未完成（慢>5倍）

**数据量对比:**
- Chunk 0: 274MB → 推测~3-4M products
- Chunk 1: 190MB → 2.47M products
- Chunk 2: 671MB → 7.31M products ⚠️

**单位时间产出:**
- Chunk 0: 1476 mol/s
- Chunk 1: 995 mol/s (慢33%)
- Chunk 2: <500 mol/s? (推算，未完成)

### 2.2 可能的深层次原因（需验证）

#### 假设1：分子复杂度梯度分布 🔍 **优先调查**

**问题:** polyphenol-2X数据集可能按某种特征排序，导致：
- 前1M (chunk 0): 简单分子，少phenolic OH
- 中1M (chunk 1): 中等复杂度
- 后段 (chunk 2+): 高度复杂，多phenolic OH → 指数级衍生物

**验证方法:**
```python
# 检查各chunk的分子特征分布
import pyarrow.parquet as pq
from rdkit import Chem
from rdkit.Chem import Descriptors

for i in range(3):
    chunk_file = f"data/output/transforms/polyphenol-2X_BATCHED/chunks/chunk_{i:03d}_input.parquet"
    table = pq.read_table(chunk_file)
    smiles_list = table['smiles'].to_pylist()[:1000]  # 采样1000个

    phenolic_oh_counts = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            # 统计phenolic OH数量
            pattern = Chem.MolFromSmarts('[OH][c]')
            matches = mol.GetSubstructMatches(pattern)
            phenolic_oh_counts.append(len(matches))

    print(f"Chunk {i}:")
    print(f"  Avg phenolic OH: {sum(phenolic_oh_counts)/len(phenolic_oh_counts):.2f}")
    print(f"  Max phenolic OH: {max(phenolic_oh_counts)}")
    print(f"  分子量分布: {[Descriptors.MolWt(Chem.MolFromSmiles(s)) for s in smiles_list[:100]]}")
```

**如果验证成立，解决方案:**
- 随机打乱chunks顺序
- 或根据复杂度预分类，平衡分配
- 或动态调整每chunk的row数

#### 假设2：内存碎片化累积 🔍 **次优先**

**问题:**
- RDKit C++对象频繁创建/销毁
- Python GC与C++ memory不同步
- 长时间运行后内存碎片化严重

**验证方法:**
```bash
# 监控chunk处理过程中的内存趋势
grep "memory:" chunk_002_output/transform.log | awk '{print $NF}' | sed 's/%//' > memory_trend.txt

# 分析是否有内存泄漏或碎片化
python -c "
import matplotlib.pyplot as plt
import numpy as np
data = np.loadtxt('memory_trend.txt')
plt.plot(data)
plt.xlabel('Flush Number')
plt.ylabel('Memory %')
plt.title('Chunk 2 Memory Trend')
plt.savefig('chunk2_memory_trend.png')
print(f'Memory start: {data[0]:.1f}%')
print(f'Memory end: {data[-1]:.1f}%')
print(f'Memory drift: {data[-1] - data[0]:.1f}%')
"
```

**如果存在内存漂移，解决方案:**
- 定期重启worker pool
- 减少max-in-flight降低并发压力
- 增加explicit gc.collect()

#### 假设3：Windows ProcessPoolExecutor限制 🔍 **需Linux对比**

**问题:**
- Windows multiprocessing已知有长时间运行稳定性问题
- Process serialization overhead
- File handle泄漏

**验证方法:**
```bash
# 在Linux上运行相同配置的chunk 2测试
# 对比执行时间和稳定性
```

**如果是Windows特定问题:**
- 迁移到Linux服务器（已在计划中）
- 或使用单进程mode (workers=1)

#### 假设4：I/O瓶颈

**问题:**
- 671MB文件写入是否导致I/O阻塞
- Parquet写入性能下降

**验证方法:**
```bash
# 监控I/O wait during chunk processing
# Windows: Resource Monitor -> Disk activity
# Linux: iostat -x 5
```

**解决方案:**
- 使用SSD存储
- 增加flush_interval减少写入频率

### 2.3 诊断计划（下一会话执行）

**Step 1: 数据特征分析（30分钟）**
```bash
cd E:/Projects/halogenator
python scripts/diagnose_chunk_complexity.py  # 需创建此脚本
```

**Step 2: 内存趋势分析（10分钟）**
```bash
python scripts/analyze_memory_trend.py chunk_002_output/transform.log
```

**Step 3: 单chunk深度测试（2小时）**
```bash
# 在高verbosity下重跑chunk 2，监控所有指标
python batch_transform_pipeline.py --single-chunk 2 --verbose
```

**Step 4: Linux对比测试（如果可用）**
```bash
# 在Linux上运行相同测试
```

---

## 第三部分：待完成任务详细说明

### 任务A：深入调查Timeout根本原因 🎯 **最高优先级**

**目标:**
1. 确认chunk 2+处理慢的根本原因
2. 找到稳定可靠的解决方案
3. 避免简单粗暴的"增加timeout"

**实施步骤:**

#### A1. 创建诊断工具

**文件:** `scripts/diagnose_chunk_complexity.py`

```python
#!/usr/bin/env python3
"""
Diagnose chunk complexity distribution to understand timeout issues.
"""

import pyarrow.parquet as pq
from rdkit import Chem
from rdkit.Chem import Descriptors
import numpy as np
from pathlib import Path

def analyze_chunk_complexity(chunk_id, sample_size=1000):
    """Analyze molecular complexity in a chunk."""

    chunk_file = Path(f"data/output/transforms/polyphenol-2X_BATCHED/chunks/chunk_{chunk_id:03d}_input.parquet")

    if not chunk_file.exists():
        print(f"Chunk {chunk_id} input file not found")
        return None

    table = pq.read_table(chunk_file)
    total_rows = len(table)

    # Sample molecules
    sample_indices = np.random.choice(total_rows, min(sample_size, total_rows), replace=False)
    smiles_list = [table['smiles'][i].as_py() for i in sample_indices]

    metrics = {
        'chunk_id': chunk_id,
        'total_rows': total_rows,
        'phenolic_oh_counts': [],
        'mol_weights': [],
        'num_aromatic_rings': [],
        'num_rotatable_bonds': [],
        'complexity_scores': []
    }

    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue

        # Count phenolic OH groups
        pattern = Chem.MolFromSmarts('[OH][c]')
        phenolic_oh = len(mol.GetSubstructMatches(pattern)) if pattern else 0

        # Molecular descriptors
        mw = Descriptors.MolWt(mol)
        aromatic_rings = Descriptors.NumAromaticRings(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)

        # Complexity score (custom metric)
        # Assume each phenolic OH can generate multiple products (per_site mode)
        # Complexity ~ phenolic_oh^k where k depends on transformation
        complexity = phenolic_oh ** 2  # Rough approximation

        metrics['phenolic_oh_counts'].append(phenolic_oh)
        metrics['mol_weights'].append(mw)
        metrics['num_aromatic_rings'].append(aromatic_rings)
        metrics['num_rotatable_bonds'].append(rotatable_bonds)
        metrics['complexity_scores'].append(complexity)

    # Statistics
    results = {
        'chunk_id': chunk_id,
        'total_rows': total_rows,
        'sample_size': len(metrics['phenolic_oh_counts']),
        'phenolic_oh': {
            'mean': np.mean(metrics['phenolic_oh_counts']),
            'median': np.median(metrics['phenolic_oh_counts']),
            'std': np.std(metrics['phenolic_oh_counts']),
            'max': np.max(metrics['phenolic_oh_counts']),
            'min': np.min(metrics['phenolic_oh_counts'])
        },
        'mol_weight': {
            'mean': np.mean(metrics['mol_weights']),
            'std': np.std(metrics['mol_weights'])
        },
        'complexity_score': {
            'mean': np.mean(metrics['complexity_scores']),
            'median': np.median(metrics['complexity_scores']),
            'max': np.max(metrics['complexity_scores'])
        },
        'predicted_products_per_mol': np.mean(metrics['complexity_scores']) * 3  # Rough estimate
    }

    return results

def main():
    print("="*80)
    print("CHUNK COMPLEXITY ANALYSIS")
    print("="*80)

    # Analyze chunks 0-2 (and optionally more)
    chunks_to_analyze = [0, 1, 2, 3, 4]  # First 5 chunks

    all_results = []
    for chunk_id in chunks_to_analyze:
        print(f"\nAnalyzing Chunk {chunk_id}...")
        result = analyze_chunk_complexity(chunk_id)
        if result:
            all_results.append(result)

            print(f"  Total rows: {result['total_rows']:,}")
            print(f"  Phenolic OH (mean): {result['phenolic_oh']['mean']:.2f}")
            print(f"  Phenolic OH (max): {result['phenolic_oh']['max']}")
            print(f"  Complexity score (mean): {result['complexity_score']['mean']:.1f}")
            print(f"  Predicted products/mol: {result['predicted_products_per_mol']:.1f}")

    # Comparison
    print("\n" + "="*80)
    print("COMPARATIVE ANALYSIS")
    print("="*80)

    if len(all_results) >= 2:
        baseline = all_results[0]
        for result in all_results[1:]:
            complexity_ratio = result['complexity_score']['mean'] / baseline['complexity_score']['mean']
            print(f"\nChunk {result['chunk_id']} vs Chunk 0:")
            print(f"  Complexity ratio: {complexity_ratio:.2f}x")
            print(f"  Phenolic OH ratio: {result['phenolic_oh']['mean'] / baseline['phenolic_oh']['mean']:.2f}x")

            if complexity_ratio > 2.0:
                print(f"  ⚠️ WARNING: Chunk {result['chunk_id']} is {complexity_ratio:.1f}x more complex!")
                print(f"     Expected processing time: ~{44 * complexity_ratio:.0f} minutes")

    # Save results
    import json
    with open('chunk_complexity_analysis.json', 'w') as f:
        json.dump(all_results, f, indent=2)

    print(f"\n✓ Results saved to: chunk_complexity_analysis.json")

if __name__ == '__main__':
    main()
```

**运行:**
```bash
python scripts/diagnose_chunk_complexity.py
```

#### A2. 内存趋势分析

**文件:** `scripts/analyze_memory_trend.py`

```python
#!/usr/bin/env python3
"""
Analyze memory usage trends from transform logs.
"""

import re
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def analyze_memory_trend(log_file):
    """Extract and analyze memory trends from log."""

    log_path = Path(log_file)
    if not log_path.exists():
        print(f"Log file not found: {log_file}")
        return

    # Extract memory values
    memory_values = []
    flush_numbers = []

    with open(log_path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            # Match: "Pre-flush memory: XX.X%"
            match = re.search(r'Pre-flush memory: ([0-9.]+)%', line)
            if match:
                memory_values.append(float(match.group(1)))

            # Match: "Flush #XXX:"
            flush_match = re.search(r'Flush #(\d+):', line)
            if flush_match:
                flush_numbers.append(int(flush_match.group(1)))

    if not memory_values:
        print("No memory data found in log")
        return

    memory_array = np.array(memory_values)

    # Statistics
    print("="*80)
    print("MEMORY TREND ANALYSIS")
    print("="*80)
    print(f"Log file: {log_file}")
    print(f"Total samples: {len(memory_array)}")
    print(f"\nMemory Statistics:")
    print(f"  Start: {memory_array[0]:.1f}%")
    print(f"  End: {memory_array[-1]:.1f}%")
    print(f"  Mean: {np.mean(memory_array):.1f}%")
    print(f"  Std: {np.std(memory_array):.2f}%")
    print(f"  Min: {np.min(memory_array):.1f}%")
    print(f"  Max: {np.max(memory_array):.1f}%")
    print(f"  Drift: {memory_array[-1] - memory_array[0]:+.1f}%")

    # Check for memory leak
    if len(memory_array) > 100:
        # Linear regression to detect trend
        x = np.arange(len(memory_array))
        slope, intercept = np.polyfit(x, memory_array, 1)

        print(f"\nTrend Analysis:")
        print(f"  Slope: {slope:.6f}% per flush")
        print(f"  Projected drift (1000 flushes): {slope * 1000:.1f}%")

        if abs(slope) > 0.001:
            print(f"  ⚠️ WARNING: Detected memory trend!")
            if slope > 0:
                print(f"     Memory is gradually increasing (possible leak)")
            else:
                print(f"     Memory is gradually decreasing (unusual)")

    # Plot
    plt.figure(figsize=(12, 6))
    plt.plot(memory_array, alpha=0.7, linewidth=0.5)
    plt.axhline(y=70, color='r', linestyle='--', label='70% limit')
    plt.xlabel('Flush Number')
    plt.ylabel('Memory %')
    plt.title(f'Memory Trend: {log_path.name}')
    plt.legend()
    plt.grid(True, alpha=0.3)

    output_file = log_path.parent / f"{log_path.stem}_memory_trend.png"
    plt.savefig(output_file, dpi=150)
    print(f"\n✓ Plot saved to: {output_file}")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python analyze_memory_trend.py <log_file>")
        print("\nExample:")
        print("  python analyze_memory_trend.py data/output/transforms/polyphenol-2X_BATCHED/chunks/chunk_002_output/transform.log")
        sys.exit(1)

    analyze_memory_trend(sys.argv[1])
```

**运行:**
```bash
python scripts/analyze_memory_trend.py data/output/transforms/polyphenol-2X_BATCHED/chunks/chunk_002_output/transform.log
```

#### A3. 基于分析结果的解决方案

**场景1: 如果是分子复杂度问题**

解决方案：动态chunk大小
```python
# 修改 batch_transform_pipeline.py
# 在 create_chunks() 中添加复杂度感知分割

def create_chunks_adaptive(self):
    """Create chunks with adaptive sizing based on complexity."""

    # First pass: analyze complexity
    complexity_scores = self.estimate_complexity_by_rows(
        self.input_file,
        sample_rate=0.01  # Sample 1%
    )

    # Second pass: create variable-sized chunks
    chunks = []
    current_start = 0
    chunk_id = 0

    target_complexity_per_chunk = np.median(complexity_scores) * 1_000_000

    while current_start < total_rows:
        # Estimate how many rows to include based on complexity
        estimated_rows = self.estimate_rows_for_target_complexity(
            complexity_scores[current_start:],
            target_complexity_per_chunk
        )

        chunk_end = min(current_start + estimated_rows, total_rows)

        # Create chunk...
        chunk_id += 1
        current_start = chunk_end
```

**场景2: 如果是内存问题**

解决方案：定期重启workers
```python
# 在 08_transform_library_v2.py 中添加
# 每处理N个batch后重启worker pool

if batch_idx % 10 == 0 and batch_idx > 0:
    logger.info("Restarting worker pool to clear memory...")
    executor.shutdown(wait=True)
    executor = ProcessPoolExecutor(max_workers=workers)
```

**场景3: 如果是Windows特定问题**

立即迁移到Linux（见任务D）

#### A4. 重新评估Timeout策略

**不要简单增加timeout，而是:**

1. **设置合理timeout基准**
```python
# 基于复杂度动态计算timeout
base_timeout = 3600  # 1 hour baseline
complexity_factor = chunk_complexity_score / baseline_complexity
timeout = base_timeout * complexity_factor * 1.5  # 1.5x safety margin
```

2. **添加进度监控**
```python
# 检测卡死 vs 正常缓慢处理
# 如果长时间无flush输出 → 真的卡死，应该kill
# 如果持续有flush → 正常但慢，应该继续等待
```

---

### 任务B：Repository整理与清理 🗂️ **高优先级**

**目标:**
1. 规范目录结构
2. 删除临时/测试文件
3. 保留重要输出和文档
4. 准备Git提交

**当前目录混乱状态:**

```
E:/Projects/halogenator/
├── *.md                    ← 散落的文档（30+个）
├── *.py                    ← 测试脚本（20+个）
├── *.log                   ← 日志文件（10+个）
├── *.txt                   ← 临时文件
├── configs/                ← 配置文件（部分临时）
├── data/
│   ├── output/
│   │   ├── transforms/
│   │   │   ├── OPT_*       ← Grid search测试输出（需删除）
│   │   │   ├── TEST_*      ← 各种测试输出（需删除）
│   │   │   └── VALIDATION_* ← 验证测试（需删除）
│   │   └── nplike_v2/      ← 正式数据（保留）
│   ├── test/               ← 测试数据（需删除）
│   ├── viz/                ← 可视化临时输出（需删除）
│   └── viz_v2/             ← 可视化临时输出（需删除）
├── scripts/                ← 脚本混乱
└── src/                    ← 源代码（保留）
```

#### B1. 创建目标目录结构

**期望结构:**
```
E:/Projects/halogenator/
├── docs/                   ← 所有.md文档
│   ├── archive/            ← 旧版本文档
│   ├── guides/             ← 使用指南
│   └── reports/            ← 会话报告
├── logs/                   ← 归档日志
├── scripts/
│   ├── production/         ← 生产脚本
│   ├── diagnosis/          ← 诊断工具
│   └── archive/            ← 临时脚本归档
├── configs/
│   ├── production/         ← 生产配置
│   └── archive/            ← 测试配置
├── data/
│   └── output/
│       ├── transforms/
│       │   └── polyphenol-2X_BATCHED/  ← 唯一保留的正式输出
│       └── nplike_v2/      ← 保留
├── src/                    ← 源代码
└── tests/                  ← 测试代码
```

#### B2. 整理脚本

**文件:** `scripts/cleanup_repository.py`

```python
#!/usr/bin/env python3
"""
Repository cleanup and reorganization script.
"""

import shutil
from pathlib import Path
import json

class RepositoryCleanup:
    def __init__(self, repo_root):
        self.repo_root = Path(repo_root)
        self.dry_run = True  # Safety first
        self.actions = []

    def analyze(self):
        """Analyze current state and plan actions."""

        print("="*80)
        print("REPOSITORY CLEANUP ANALYSIS")
        print("="*80)

        # Find all .md files
        md_files = list(self.repo_root.glob("*.md"))
        print(f"\n📄 Markdown files in root: {len(md_files)}")
        for f in md_files:
            self.actions.append({
                'type': 'move',
                'source': str(f),
                'dest': f"docs/{f.name}",
                'reason': 'Documentation to docs/'
            })

        # Find test directories/files to delete
        test_patterns = [
            "data/output/transforms/OPT_*",
            "data/output/transforms/TEST_*",
            "data/output/transforms/VALIDATION_*",
            "data/test/",
            "data/viz/",
            "data/viz_v2/",
            "data/viz_base_libs/",
            "tmp/"
        ]

        for pattern in test_patterns:
            matches = list(self.repo_root.glob(pattern))
            for match in matches:
                if match.exists():
                    size = self.get_size(match)
                    self.actions.append({
                        'type': 'delete',
                        'path': str(match),
                        'size': size,
                        'reason': 'Test/temporary data'
                    })

        # Find .log files
        log_files = list(self.repo_root.glob("*.log"))
        print(f"\n📝 Log files in root: {len(log_files)}")
        for f in log_files:
            self.actions.append({
                'type': 'move',
                'source': str(f),
                'dest': f"logs/{f.name}",
                'reason': 'Log to logs/'
            })

        # Find standalone .py scripts
        py_files = [
            f for f in self.repo_root.glob("*.py")
            if f.name not in ['setup.py', '__init__.py']
        ]
        print(f"\n🐍 Python scripts in root: {len(py_files)}")
        for f in py_files:
            # Categorize
            if 'test' in f.name.lower() or 'validate' in f.name.lower():
                dest = f"scripts/archive/{f.name}"
            elif 'diagnose' in f.name.lower() or 'analyze' in f.name.lower():
                dest = f"scripts/diagnosis/{f.name}"
            elif 'optimize' in f.name.lower() or 'batch' in f.name.lower():
                dest = f"scripts/production/{f.name}"
            else:
                dest = f"scripts/archive/{f.name}"

            self.actions.append({
                'type': 'move',
                'source': str(f),
                'dest': dest,
                'reason': 'Script organization'
            })

        # Summary
        print(f"\n" + "="*80)
        print(f"PLANNED ACTIONS SUMMARY")
        print(f"="*80)

        action_counts = {}
        for action in self.actions:
            action_type = action['type']
            action_counts[action_type] = action_counts.get(action_type, 0) + 1

        for action_type, count in action_counts.items():
            print(f"  {action_type.upper()}: {count} items")

        # Calculate space to be freed
        delete_size = sum(
            action.get('size', 0)
            for action in self.actions
            if action['type'] == 'delete'
        )
        print(f"\n💾 Space to be freed: {delete_size / (1024**3):.2f} GB")

        return self.actions

    def get_size(self, path):
        """Get size of file or directory in bytes."""
        path = Path(path)
        if path.is_file():
            return path.stat().st_size
        elif path.is_dir():
            return sum(f.stat().st_size for f in path.rglob('*') if f.is_file())
        return 0

    def execute(self, dry_run=True):
        """Execute planned actions."""

        self.dry_run = dry_run

        mode_str = "DRY RUN" if dry_run else "EXECUTING"
        print(f"\n{'='*80}")
        print(f"{mode_str} - Repository Cleanup")
        print(f"{'='*80}")

        success_count = 0
        error_count = 0

        for i, action in enumerate(self.actions, 1):
            print(f"\n[{i}/{len(self.actions)}] {action['type'].upper()}: {action.get('source', action.get('path'))}")

            try:
                if action['type'] == 'move':
                    if not dry_run:
                        source = Path(action['source'])
                        dest = self.repo_root / action['dest']
                        dest.parent.mkdir(parents=True, exist_ok=True)
                        shutil.move(str(source), str(dest))
                    print(f"  → {action['dest']}")
                    success_count += 1

                elif action['type'] == 'delete':
                    path = Path(action['path'])
                    size_mb = action.get('size', 0) / (1024**2)
                    if not dry_run:
                        if path.is_dir():
                            shutil.rmtree(path)
                        else:
                            path.unlink()
                    print(f"  ✗ Deleted ({size_mb:.1f} MB)")
                    success_count += 1

            except Exception as e:
                print(f"  ✗ ERROR: {e}")
                error_count += 1

        # Summary
        print(f"\n{'='*80}")
        print(f"CLEANUP {'DRY RUN' if dry_run else 'EXECUTION'} COMPLETE")
        print(f"{'='*80}")
        print(f"  Success: {success_count}/{len(self.actions)}")
        print(f"  Errors: {error_count}")

        if dry_run:
            print(f"\n⚠️  This was a DRY RUN. No files were actually moved or deleted.")
            print(f"   Review the actions above. To execute for real, run:")
            print(f"   python scripts/cleanup_repository.py --execute")
        else:
            print(f"\n✓ Repository cleanup complete!")

            # Save cleanup report
            report_file = self.repo_root / "docs" / "reports" / "cleanup_report.json"
            report_file.parent.mkdir(parents=True, exist_ok=True)
            with open(report_file, 'w') as f:
                json.dump({
                    'actions': self.actions,
                    'success_count': success_count,
                    'error_count': error_count
                }, f, indent=2)
            print(f"   Report saved to: {report_file}")

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Clean up and reorganize repository')
    parser.add_argument('--execute', action='store_true', help='Execute cleanup (default is dry run)')
    parser.add_argument('--repo', default='.', help='Repository root path')
    args = parser.parse_args()

    cleanup = RepositoryCleanup(args.repo)
    cleanup.analyze()
    cleanup.execute(dry_run=not args.execute)

if __name__ == '__main__':
    main()
```

**执行步骤:**

1. **Dry run分析:**
```bash
cd E:/Projects/halogenator
python scripts/cleanup_repository.py
# 仔细检查输出，确认要删除/移动的文件
```

2. **审查删除列表:**
```bash
# 特别注意以下保留项：
# - data/output/nplike_v2/ (原始处理数据)
# - data/output/transforms/polyphenol-2X_BATCHED/ (当前pipeline输出)
# - src/ (源代码)
# - configs/halogen_rules_by_class.yaml (生产配置)
# - configs/transforms.yaml (生产配置)
```

3. **确认后执行:**
```bash
python scripts/cleanup_repository.py --execute
```

#### B3. 文档分类整理

**手动整理建议:**

```bash
mkdir -p docs/{archive,guides,reports}

# 会话报告
mv SESSION_* docs/reports/
mv IMPLEMENTATION_* docs/reports/
mv WORK_REPORT_* docs/reports/

# 技术文档
mv ROOT_CAUSE_ANALYSIS_* docs/reports/
mv SOLUTION_VALIDATED_* docs/reports/

# 使用指南
mv USAGE_GUIDE.md docs/guides/
mv SUGAR_FILTER_GUIDE.md docs/guides/
mv TEST_OBSERVATION_GUIDE.md docs/guides/

# 其他归档
mv *.md docs/archive/  # 剩余文档
```

#### B4. 配置文件整理

```bash
cd configs/

# 保留生产配置
mkdir -p production archive

# 生产配置（不移动）
# - halogen_rules_by_class.yaml
# - transforms.yaml

# 归档测试配置
mv *_k2.yaml archive/
mv *_k3.yaml archive/
mv test_*.yaml archive/
mv macro_verify_*.yaml archive/
```

---

### 任务C：Git提交与Github推送 📤 **中优先级**

**目标:**
1. 提交所有有价值的代码和文档
2. 推送到Github远程仓库
3. 准备在Linux上clone

#### C1. 准备.gitignore

**检查/更新 `.gitignore`:**

```bash
# 查看当前.gitignore
cat .gitignore

# 应该包含：
# Data files
data/output/
*.parquet
*.pkl

# Logs
*.log
logs/

# Temporary files
tmp/
temp/
__pycache__/
*.pyc
*.pyo

# IDE
.vscode/
.idea/

# OS
.DS_Store
Thumbs.db

# Large files
*.gz
*.zip
*.tar

# Test outputs
*_TEST_*
*_OPT_*
optimization_results_*.csv
```

#### C2. Git提交流程

**步骤:**

```bash
cd E:/Projects/halogenator

# 1. 查看当前分支和状态
git status
git branch

# 2. 如果在fix分支，考虑合并或创建新分支
# 当前分支：fix/pr2-contract-and-sugar-accept
# 建议创建新分支用于batch pipeline工作

git checkout -b feature/batch-transform-pipeline

# 3. Stage整理后的文件
git add docs/
git add scripts/production/
git add scripts/diagnosis/
git add batch_transform_pipeline.py
git add optimize_parameters.py
git add .gitignore

# 4. 查看将要提交的内容
git status
git diff --staged

# 5. 创建详细的commit message
cat > commit_message.txt << 'EOF'
feat: Batch transform pipeline with parameter optimization

Major changes:
1. Parameter optimization grid search framework
   - Created optimize_parameters.py for testing worker/batch configs
   - Tested 4 configurations on 150K subset
   - Identified safe config: workers=12, max-in-flight=8

2. Batch processing pipeline implementation
   - Created batch_transform_pipeline.py
   - Automatic chunking of large datasets (1M rows/chunk)
   - Checkpoint and resume capability
   - State tracking with pipeline_state.json
   - Automatic result merging

3. Repository reorganization
   - Moved documentation to docs/
   - Organized scripts into production/diagnosis/archive
   - Cleaned up test outputs and temporary files
   - Improved .gitignore

4. Diagnostic tools
   - diagnose_chunk_complexity.py for complexity analysis
   - analyze_memory_trend.py for memory profiling
   - cleanup_repository.py for automated cleanup

Technical details:
- Workers=16, max-in-flight=6, batch=50K configuration
- 4-hour timeout per chunk (needs investigation)
- Successfully processed 2/14 chunks of polyphenol-2X
- Discovered chunk complexity variance issue

Known issues:
- Timeout issues on complex chunks (needs deeper investigation)
- Chunk 2 timeout after 4 hours (7.3M products generated)
- Possible Windows ProcessPoolExecutor limitations

Next steps:
- Investigate chunk complexity distribution
- Test on Linux for comparison
- Optimize timeout strategy based on complexity

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>
EOF

git commit -F commit_message.txt

# 6. 推送到Github
# 如果是新分支，需要设置upstream
git push -u origin feature/batch-transform-pipeline

# 7. 创建Pull Request (可选)
# 如果需要code review，在Github上创建PR
```

#### C3. 验证远程推送

```bash
# 检查远程仓库
git remote -v

# 如果没有配置，添加Github remote
# git remote add origin https://github.com/YOUR_USERNAME/halogenator.git

# 推送后验证
git log --oneline -5
git branch -r
```

#### C4. 标签重要版本

```bash
# 为batch pipeline实现打标签
git tag -a v0.2.0-batch-pipeline -m "Batch processing pipeline implementation"
git push origin v0.2.0-batch-pipeline
```

---

### 任务D：Linux迁移与兼容性验证 🐧 **高优先级**

**目标:**
1. 在Linux服务器上成功运行
2. 验证性能差异
3. 对比Windows的稳定性问题

#### D1. 环境准备

**Linux服务器要求:**
```
OS: Ubuntu 20.04+ / CentOS 7+
Python: 3.8+
RAM: 32GB+
Storage: 500GB+ (SSD推荐)
CPU: 16+ cores
```

**克隆仓库:**

```bash
# SSH登录Linux服务器
ssh user@linux-server

# 克隆仓库
cd /home/user/projects/  # 或其他工作目录
git clone https://github.com/YOUR_USERNAME/halogenator.git
cd halogenator

# 切换到batch pipeline分支
git checkout feature/batch-transform-pipeline
```

#### D2. 依赖安装

**创建conda环境:**

```bash
# 如果没有conda，先安装Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# 创建环境
conda create -n halo-p0 python=3.9
conda activate halo-p0

# 安装依赖
pip install rdkit-pypi
pip install pyarrow pandas numpy
pip install matplotlib  # for diagnostics
pip install tqdm  # optional, for progress bars

# 验证安装
python -c "from rdkit import Chem; print('RDKit OK')"
python -c "import pyarrow; print('PyArrow OK')"
```

**检查requirements.txt:**

如果没有，创建：
```bash
cat > requirements.txt << EOF
rdkit-pypi>=2022.9.5
pyarrow>=10.0.0
pandas>=1.5.0
numpy>=1.23.0
pyyaml>=6.0
tqdm>=4.64.0
matplotlib>=3.6.0
EOF

pip install -r requirements.txt
```

#### D3. 数据传输

**方案A: rsync (推荐)**

```bash
# 从Windows传输到Linux
# 在Windows上 (使用WSL或Git Bash)
rsync -avz --progress \
  /e/Projects/halogenator/data/output/nplike_v2/ \
  user@linux-server:/home/user/projects/halogenator/data/output/nplike_v2/

# 或只传输polyphenol-2X
rsync -avz --progress \
  /e/Projects/halogenator/data/output/nplike_v2/polyphenol-2X/ \
  user@linux-server:/home/user/projects/halogenator/data/output/nplike_v2/polyphenol-2X/
```

**方案B: 使用scp**

```bash
# 压缩后传输
cd /e/Projects/halogenator/data/output/nplike_v2/
tar -czf polyphenol-2X.tar.gz polyphenol-2X/

scp polyphenol-2X.tar.gz user@linux-server:/home/user/projects/halogenator/data/

# 在Linux上解压
ssh user@linux-server
cd /home/user/projects/halogenator/data/output/nplike_v2/
tar -xzf ../../polyphenol-2X.tar.gz
```

#### D4. 兼容性验证测试

**Test 1: 基础功能测试**

```bash
cd /home/user/projects/halogenator

# 测试transform脚本可以运行
python scripts/08_transform_library_v2.py --help

# 测试batch pipeline
python batch_transform_pipeline.py --help
```

**Test 2: 小规模验证（1小时）**

```bash
# 创建测试chunk (100K rows)
python -c "
import pyarrow.parquet as pq
table = pq.read_table('data/output/nplike_v2/polyphenol-2X/products.parquet')
subset = table.slice(0, 100000)
pq.write_table(subset, 'data/test_100k.parquet')
print('Created test_100k.parquet')
"

# 运行单次transform测试
python scripts/08_transform_library_v2.py apply \
  --input data/test_100k.parquet \
  --outdir data/output/test_linux_100k \
  --xf-config configs/transforms.yaml \
  --xf-name FG_PHENOL_OH__OH__TO__OMe \
  --workers 16 \
  --max-in-flight 6 \
  --batch-size 50000 \
  --use-bloom-filter

# 检查结果
ls -lh data/output/test_linux_100k/
cat data/output/test_linux_100k/SUMMARY.json
```

**Test 3: Chunk 2 对比测试（关键！）**

```bash
# 传输chunk 2输入文件
rsync -avz \
  /e/Projects/halogenator/data/output/transforms/polyphenol-2X_BATCHED/chunks/chunk_002_input.parquet \
  user@linux-server:/home/user/projects/halogenator/data/chunk_002_input.parquet

# 在Linux上运行chunk 2测试
# 使用相同配置，监控是否也会timeout
time python scripts/08_transform_library_v2.py apply \
  --input data/chunk_002_input.parquet \
  --outdir data/output/test_linux_chunk2 \
  --xf-config configs/transforms.yaml \
  --xf-name FG_PHENOL_OH__OH__TO__OMe \
  --workers 16 \
  --max-in-flight 6 \
  --batch-size 50000 \
  --use-bloom-filter \
  --bloom-expected-items 10000000 \
  2>&1 | tee linux_chunk2_test.log

# 关键对比点：
# 1. 是否在4小时内完成？
# 2. 如果完成，用时多少？
# 3. 内存使用是否更稳定？
# 4. 产品数量是否一致？

# 如果Linux能在2-3小时完成chunk 2，说明是Windows特定问题！
```

#### D5. 性能对比分析

**创建对比报告:**

```python
# scripts/compare_windows_linux.py
"""
Compare Windows vs Linux performance for same workload.
"""

import json
from pathlib import Path

def compare_results():
    # Windows chunk 0 results
    windows_c0 = {
        'platform': 'Windows',
        'chunk': 0,
        'time': 2652,  # seconds
        'throughput': 1476.7
    }

    # Linux test results (to be filled)
    linux_test = {
        'platform': 'Linux',
        'chunk': 'test_100k',
        'time': None,  # From actual test
        'throughput': None
    }

    # Analysis...
    print("Platform Performance Comparison")
    print("="*80)
    # ...
```

#### D6. 生产部署决策

**基于Linux测试结果决策:**

**场景A: Linux显著更快/稳定**
→ 立即迁移所有任务到Linux
→ 在Linux上运行完整batch pipeline

**场景B: Linux与Windows相似**
→ 问题不是平台特定
→ 需要更深入调查chunk complexity
→ 可以在任一平台继续

**场景C: Linux更慢**
→ 检查Linux环境配置（CPU governor, NUMA等）
→ 优化Linux系统参数

---

## 第四部分：关键文件索引

### 生产代码

```
src/halogenator/
├── __init__.py
├── cli.py                  # CLI入口
├── enumerate_k.py          # K次枚举核心
├── transform_engine.py     # Transform引擎
├── schema.py              # 数据schema
└── sugar_mask.py          # Sugar过滤

scripts/
├── production/
│   ├── 08_transform_library_v2.py  # Transform主脚本
│   ├── batch_transform_pipeline.py # Batch pipeline
│   └── optimize_parameters.py      # 参数优化
└── diagnosis/
    ├── diagnose_chunk_complexity.py  # 待创建
    └── analyze_memory_trend.py       # 待创建
```

### 配置文件

```
configs/
├── halogen_rules_by_class.yaml  # NP分类halogen规则
└── transforms.yaml              # Transform定义
```

### 当前输出

```
data/output/
├── nplike_v2/
│   └── polyphenol-2X/
│       └── products.parquet     # 13.79M rows, 504MB
└── transforms/
    └── polyphenol-2X_BATCHED/
        ├── pipeline_state.json  # Pipeline状态
        └── chunks/
            ├── chunk_000_output/  # ✅ 完成
            ├── chunk_001_output/  # ✅ 完成
            ├── chunk_002_output/  # ❌ Timeout (部分输出671MB)
            └── chunk_003-013_*/   # ⏸️ 待处理
```

### 文档

```
docs/
├── reports/
│   └── SESSION_HANDOFF_BATCH_PIPELINE_AND_CLEANUP.md  # 本文档
├── guides/
│   └── USAGE_GUIDE.md
└── archive/
    └── *.md
```

---

## 第五部分：技术细节与注意事项

### 5.1 已知的系统性问题

#### 问题1: ProcessPoolExecutor在Windows上的不稳定性

**表现:**
- 150K测试成功，500K+失败
- Worker进程被突然终止（BrokenProcessPool）
- 长时间运行后更容易崩溃

**原因分析:**
- Windows multiprocessing使用spawn而非fork
- Process serialization overhead
- 可能的file handle泄漏

**缓解措施:**
- 使用较小的chunk size
- 定期重启worker pool
- 迁移到Linux

#### 问题2: Chunk复杂度不均匀

**表现:**
- Chunk 0: 44分钟
- Chunk 2: >4小时 (5.4倍差异)
- 产品数量: 2.47M vs 7.31M (3倍差异)

**可能原因:**
- polyphenol-2X可能按分子量或复杂度排序
- 后续chunks包含更多phenolic OH的分子
- Per-site mode导致指数级产品增长

**需要验证:**
- 运行complexity analysis脚本
- 检查原始数据集的排序方式

#### 问题3: 内存使用的非线性增长

**观察:**
- Chunk 0: 未明确记录内存
- Chunk 1: 内存稳定
- Chunk 2: 45.7%稳定（但处理时间异常长）

**可能原因:**
- 不是内存问题（45.7%很安全）
- 而是计算复杂度问题
- 或I/O瓶颈

### 5.2 配置参数的权衡

**Workers数量:**
- 更多workers → 更快，但内存占用高
- Windows建议：12-16
- Linux可能可以更高：24-32

**Max-in-flight:**
- 控制并发batch数量
- 太高 → 内存压力，queue buildup
- 太低 → workers空闲，效率低
- 当前：6-8是平衡点

**Batch size:**
- 50K是测试验证的值
- 更大 → 减少overhead，但增加单batch处理时间
- 更小 → 更灵活，但overhead增加

**Timeout:**
- 当前4小时不够（chunk 2证明）
- 不应简单增加，而应基于复杂度动态计算
- 或改进chunk划分策略

### 5.3 数据完整性检查

**Chunk输出验证:**

```bash
# 检查所有完成的chunks
for chunk_dir in data/output/transforms/polyphenol-2X_BATCHED/chunks/chunk_*_output/; do
    chunk_id=$(basename $chunk_dir | grep -oP '\d+')

    if [ -f "$chunk_dir/SUMMARY.json" ]; then
        products=$(jq '.unique_products' "$chunk_dir/SUMMARY.json")
        echo "Chunk $chunk_id: ✓ $products products"
    else
        if [ -f "$chunk_dir/products.parquet" ]; then
            size=$(du -h "$chunk_dir/products.parquet" | cut -f1)
            echo "Chunk $chunk_id: ⚠ No SUMMARY, but $size products.parquet exists"
        else
            echo "Chunk $chunk_id: ✗ No output"
        fi
    fi
done
```

### 5.4 Resume Pipeline的正确方法

**当前状态文件结构:**

```json
{
  "chunks": [...],
  "completed_chunks": [0, 1],
  "failed_chunks": [2],
  "start_time": "...",
  "end_time": null
}
```

**Resume方法:**

```bash
# Pipeline会自动跳过completed chunks
# 但需要手动处理failed chunks

# 选项1: 清除chunk 2的failed状态，让它重试
python -c "
import json
with open('data/output/transforms/polyphenol-2X_BATCHED/pipeline_state.json', 'r+') as f:
    state = json.load(f)
    # 从failed list移除chunk 2
    state['failed_chunks'] = []
    # 将chunk 2状态改为pending
    state['chunks'][2]['status'] = 'pending'
    f.seek(0)
    json.dump(state, f, indent=2)
    f.truncate()
print('Chunk 2 reset to pending')
"

# 选项2: 跳过chunk 2，先完成其他chunks
# 手动修改pipeline代码跳过特定chunk

# 选项3: 使用chunk 2的部分输出
# 检查products.parquet是否可用
# 可能需要手动创建SUMMARY.json
```

---

## 第六部分：执行优先级与时间估算

### 立即执行（1-2小时）

1. **运行Complexity Analysis** (30分钟)
   ```bash
   python scripts/diagnose_chunk_complexity.py
   ```

2. **运行Memory Trend Analysis** (10分钟)
   ```bash
   python scripts/analyze_memory_trend.py data/.../chunk_002_output/transform.log
   ```

3. **Repository Cleanup Dry Run** (20分钟)
   ```bash
   python scripts/cleanup_repository.py
   # 审查输出
   ```

### 短期执行（半天）

4. **Repository Cleanup Execute** (30分钟)
   ```bash
   python scripts/cleanup_repository.py --execute
   ```

5. **Git提交与推送** (30分钟)
   ```bash
   git add .
   git commit -F commit_message.txt
   git push origin feature/batch-transform-pipeline
   ```

6. **Linux环境搭建** (2小时)
   - SSH登录
   - 安装conda和依赖
   - 克隆仓库
   - 传输测试数据

### 中期执行（1-2天）

7. **Linux兼容性测试** (4小时)
   - 100K测试
   - Chunk 2对比测试
   - 性能对比分析

8. **基于诊断结果的Pipeline优化** (4-8小时)
   - 实现adaptive chunking（如果需要）
   - 或实现dynamic timeout
   - 或其他优化措施

### 长期执行（3-7天）

9. **完整Pipeline运行** (取决于优化结果)
   - 如果每chunk 2-3小时：14 chunks × 3h = 42小时
   - 如果优化后1-2小时：14 chunks × 2h = 28小时
   - 建议在Linux上24/7运行，定期检查

10. **结果验证与合并** (2小时)
    - 验证所有chunks完成
    - 运行merge操作
    - 质量检查

---

## 第七部分：成功标准

### Minimum Viable Product (必须达成)

- ✅ 完成所有14 chunks的处理
- ✅ 每个chunk有valid SUMMARY.json
- ✅ 合并后的产品数量合理（预计40-60M）
- ✅ 无数据损坏或缺失

### Optimal Goals (期望达成)

- ✅ 理解并解决timeout的根本原因
- ✅ 单chunk处理时间可预测（基于复杂度）
- ✅ 在Linux上成功运行并验证更好的性能
- ✅ Repository干净整洁，ready for production

### Stretch Goals (如有时间)

- ✅ 实现adaptive chunking
- ✅ 完整的Linux部署文档
- ✅ 自动化monitoring dashboard
- ✅ 将batch pipeline扩展到terpenoid等其他数据集

---

## 第八部分：紧急联系与回滚

### 如果遇到严重问题

**问题1: Pipeline卡死无响应**

```bash
# 查找并kill进程
ps aux | grep batch_transform_pipeline.py
kill -9 <PID>

# 检查状态文件
cat data/output/transforms/polyphenol-2X_BATCHED/pipeline_state.json

# 恢复：重新运行pipeline（会自动resume）
```

**问题2: 磁盘空间不足**

```bash
# 检查空间
df -h

# 紧急清理临时文件
rm -rf data/output/transforms/OPT_*
rm -rf data/output/transforms/TEST_*
```

**问题3: Git推送失败**

```bash
# 检查大文件
git ls-files | xargs ls -lh | sort -k5 -hr | head -20

# 如果有大文件被误加入
git rm --cached <large_file>
git commit --amend
```

### 回滚到安全状态

```bash
# 如果cleanup破坏了什么
git reset --hard HEAD~1  # 回滚最后一次commit

# 如果需要恢复文件
git checkout HEAD -- <file_path>
```

---

## 第九部分：给下一个会话的具体指令

### 会话开始时的First Steps

1. **读取本文档**（你正在做）✓

2. **确认当前环境**
   ```bash
   pwd  # 应该在 E:/Projects/halogenator
   git branch  # 应该在某个feature分支
   ls data/output/transforms/polyphenol-2X_BATCHED/pipeline_state.json  # 应该存在
   ```

3. **运行诊断脚本**
   ```bash
   python scripts/diagnose_chunk_complexity.py
   ```

4. **基于诊断结果决定下一步**
   - 如果complexity差异巨大 → 实现adaptive chunking
   - 如果无明显差异 → 考虑其他原因（内存、I/O、Windows）
   - 如果不确定 → 先在Linux上测试

### 决策树

```
START
  ├─→ Complexity明显不均？
  │   ├─ YES → 实现adaptive chunking → 重跑pipeline
  │   └─ NO → 继续下一步
  │
  ├─→ 有Linux服务器可用？
  │   ├─ YES → 立即迁移并测试 → 对比结果
  │   └─ NO → 继续Windows优化
  │
  ├─→ 决定继续Windows？
  │   ├─ 增加timeout到8小时
  │   ├─ 或减小chunk size到500K
  │   └─ 重跑pipeline
  │
  └─→ 在Linux测试？
      ├─ 搭建环境 (2小时)
      ├─ 传输数据 (1-2小时)
      ├─ 运行chunk 2测试 (2-4小时)
      └─ 基于结果决定
```

---

## 附录A：完整的命令速查表

### 诊断命令

```bash
# Complexity分析
python scripts/diagnose_chunk_complexity.py

# Memory趋势
python scripts/analyze_memory_trend.py data/output/transforms/polyphenol-2X_BATCHED/chunks/chunk_002_output/transform.log

# 检查pipeline状态
cat data/output/transforms/polyphenol-2X_BATCHED/pipeline_state.json | python -m json.tool

# 检查所有chunks状态
for i in {0..13}; do
    if [ -f "data/output/transforms/polyphenol-2X_BATCHED/chunks/chunk_$(printf '%03d' $i)_output/SUMMARY.json" ]; then
        echo "Chunk $i: DONE"
    else
        echo "Chunk $i: PENDING/FAILED"
    fi
done
```

### Pipeline操作

```bash
# 运行pipeline（自动resume）
python batch_transform_pipeline.py

# 检查正在运行的进程
ps aux | grep python | grep -E "batch_transform|08_transform"

# 查看实时日志
tail -f data/output/transforms/polyphenol-2X_BATCHED/chunks/chunk_XXX_output/transform.log
```

### Repository整理

```bash
# Dry run
python scripts/cleanup_repository.py

# 执行
python scripts/cleanup_repository.py --execute

# 手动整理文档
mkdir -p docs/{reports,guides,archive}
mv SESSION_*.md docs/reports/
mv *_GUIDE.md docs/guides/
```

### Git操作

```bash
# 状态检查
git status
git branch

# 提交
git add .
git commit -F commit_message.txt

# 推送
git push origin feature/batch-transform-pipeline

# 标签
git tag -a v0.2.0-batch-pipeline -m "Batch pipeline implementation"
git push origin v0.2.0-batch-pipeline
```

### Linux迁移

```bash
# 传输数据
rsync -avz --progress /e/Projects/halogenator/data/output/nplike_v2/polyphenol-2X/ user@linux:/path/to/halogenator/data/output/nplike_v2/polyphenol-2X/

# SSH登录
ssh user@linux-server

# 环境搭建
conda create -n halo-p0 python=3.9
conda activate halo-p0
pip install -r requirements.txt

# 测试
python scripts/08_transform_library_v2.py --help
```

---

## 附录B：预期产出文件清单

### 本会话应该产生的文件

```
scripts/
├── diagnosis/
│   ├── diagnose_chunk_complexity.py    ← 新建
│   └── analyze_memory_trend.py         ← 新建
└── production/
    └── cleanup_repository.py           ← 新建

docs/
├── reports/
│   ├── SESSION_HANDOFF_BATCH_PIPELINE_AND_CLEANUP.md  ← 本文档
│   └── chunk_complexity_analysis.json  ← 诊断输出
└── guides/
    └── (整理后的文档)

# Git
.git/
└── (新commits和tags)

# 清理后的目录结构
(见任务B的目标结构)
```

---

## 结语

本会话的核心目标是**调查timeout根本原因**并**完成repository清理准备迁移**。

**关键记忆点：**
1. ⚠️ **不要盲目增加timeout** - 先诊断根本原因
2. 🔍 **Chunk 2的7.3M产品是关键线索** - 复杂度问题
3. 🐧 **Linux测试是critical path** - 可能解决Windows问题
4. 🗂️ **清理repository是迁移前提** - 不要带着垃圾迁移
5. 📊 **运行诊断工具获得数据** - 数据驱动决策

**成功标志：**
- ✅ 理解了为什么chunk 2需要>4小时
- ✅ Repository干净整洁
- ✅ 代码已push到Github
- ✅ 在Linux上成功运行测试
- ✅ 有明确的优化方案

**如果时间有限，优先级：**
1. 运行complexity诊断（必须）
2. Linux测试（必须）
3. Repository清理（重要）
4. Git推送（重要）
5. 实现优化方案（基于1和2的结果）

祝下一个会话顺利！🚀

---

**文档版本:** 1.0
**最后更新:** 2025-12-29
**作者:** Claude Sonnet 4.5
**状态:** Ready for handoff
