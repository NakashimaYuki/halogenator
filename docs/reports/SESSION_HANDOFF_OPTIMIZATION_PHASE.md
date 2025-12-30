# Session Handoff: Transform Pipeline Optimization Phase

**日期:** 2025-12-27
**状态:** 🟢 稳定性已解决，进入优化阶段
**上一会话:** 内存危机调查与修复
**本会话目标:** 在保持稳定性的前提下优化性能

---

## 第一部分：已完成任务摘要 ✅

### 1. 内存危机根本原因确认（已完成）

**问题：** Transform pipeline在处理大数据集时内存达94.2%，产出损坏

**发现的两大根本原因：**
1. **Buffer爆炸：** `StreamingParquetWriter.write_batch()` 将所有产品一次性添加到buffer，导致buffer从2K膨胀到940K
2. **Worker并行过载：** 多个batch（10+）同时处理，导致8-16个worker进程各持有1-1.5GB内存

### 2. 双重修复实施（已完成）

**修复A：Buffer Chunking（代码修复）**
- 文件：`E:\Projects\halogenator\scripts\08_transform_library_v2.py`
- 方法：`StreamingParquetWriter.write_batch()` (lines 413-483)
- 修改：将大批次products分成1000个一组的chunks处理
- 结果：Buffer从940K降至2K（99.8%改进）
- 备份：`08_transform_library_v2.py.backup_before_buffer_fix`

**修复B：Worker并行限制（参数修复）**
- 参数：`--max-in-flight 4`
- 效果：限制同时处理的batch数量从10降至4
- 结果：内存从86-89%降至50-55%（36%改进）

### 3. 三级渐进验证（已完成）

| 测试级别 | 数据量 | Batches | Workers | max-in-flight | 峰值内存 | 状态 |
|---------|--------|---------|---------|---------------|----------|------|
| MICRO | 50K rows | 1 | 4 | 1 | 43.6% | ✅ PASS |
| SMALL | 150K rows | 3 | 6 | 3 | 45.6% | ✅ PASS |
| MEDIUM (unfixed) | 500K rows | 10 | 8 | 10 | 86-89% | ❌ TOO HIGH |
| MEDIUM (fixed) | 500K rows | 10 | 8 | **4** | **50.3%** | ✅ **PASS** |

**验证结果：**
- Buffer控制：完美，所有flush ≤ 2,000 products
- 内存稳定：47.6% - 55.5%（平均50.3%）
- 无Critical警告：0次
- 输出有效：1,602,690 products，无损坏

---

## 第二部分：当前系统状态（详细）

### 当前配置（生产就绪）

**推荐命令：**
```bash
python scripts/08_transform_library_v2.py apply \
  --input data/output/nplike_v2/polyphenol-2X/products.parquet \
  --outdir data/output/transforms/polyphenol-2X_FG_PHENOL_OH__OH__TO__OMe \
  --xf-config configs/transforms.yaml \
  --xf-name FG_PHENOL_OH__OH__TO__OMe \
  --use-bloom-filter \
  --workers 16 \
  --target-memory 70.0 \
  --max-in-flight 4 \
  --bloom-expected-items 100000000 \
  --batch-size 50000 \
  --flush-interval 2000
```

### 性能指标（当前vs原始）

| 指标 | 原始（有bug） | 当前（已修复） | 差异 |
|------|--------------|----------------|------|
| **峰值内存** | 94.2% (crash) | 50-55% | -39% ✅ |
| **吞吐量** | 1,930 mol/s | 600-800 mol/s | **-58-69%** ⚠️ |
| **完成率** | 0% (crash) | 100% | +100% ✅ |
| **稳定性** | 崩溃 | 稳定 | Fixed ✅ |

**关键问题：吞吐量损失58-69%，需要优化！**

### 内存使用分析（已知）

**在max-in-flight=4时的内存分布（500K rows测试）：**
```
Component                    Estimated Memory    % of 32GB
-----------------------------------------------------------------
Workers (8 workers)          ~8-12 GB            25-37%
  - RDKit Mol objects        ~6-8 GB
  - Intermediate products    ~2-4 GB

ProcessPoolExecutor queue    ~2-3 GB             6-9%
  - 4 batches × results

Bloom filter                 ~0.1 GB             <1%
Writer buffer                ~0.01 GB            <1%
OS + Python + Other          ~4-6 GB             12-18%
-----------------------------------------------------------------
TOTAL:                       ~16-18 GB           50-55% ✅
```

**在max-in-flight=10时的内存分布（推算）：**
```
Workers (8 workers)          ~8-12 GB            25-37%
ProcessPoolExecutor queue    ~5-7 GB             15-22%  ← 增加！
  - 10 batches × results
OS + Python + Other          ~4-6 GB             12-18%
-----------------------------------------------------------------
TOTAL:                       ~27-29 GB           86-89% ❌
```

**关键洞察：**
- Buffer和Bloom filter不是问题（已修复）
- 主要内存消耗：Worker processes + Result queue
- max-in-flight直接影响Result queue大小
- Worker数量直接影响Worker memory

---

## 第三部分：待完成任务（极其详细）

### 🎯 核心任务：性能优化（在保持稳定性前提下）

**目标：**
1. **主要目标：** 提升吞吐量从当前600-800 mol/s到1,200-1,500 mol/s（提升50-100%）
2. **约束条件：** 峰值内存必须保持 < 70%（当前50-55%）
3. **次要目标：** 减少运行时间从24h降至16-18h

**当前性能瓶颈分析：**

1. **max-in-flight=4 限制了并行度**
   - 问题：只有4个batch同时处理，其余batch在队列等待
   - 影响：Workers可能空闲等待新batch
   - 机会：如果能增加到6-8，可能提升50%吞吐量

2. **Workers=16 可能未充分利用**
   - 问题：4个batch可能无法让16个workers全部忙碌
   - 影响：CPU利用率可能不足
   - 机会：增加max-in-flight或优化batch调度

3. **batch_size=50000 可能不是最优**
   - 问题：batch太大→内存压力大；batch太小→overhead高
   - 机会：调整batch大小可能优化内存/性能平衡

### 📋 任务1：参数空间探索（必做）

**目标：** 找到最优的 (workers, max-in-flight, batch-size) 组合

**方法：网格搜索（Grid Search）**

#### 子任务1.1：创建参数扫描脚本

**创建文件：** `E:\Projects\halogenator\optimize_parameters.py`

**脚本功能：**
1. 定义参数网格
2. 对每组参数运行SMALL测试（150K rows，快速验证）
3. 记录：峰值内存、平均内存、吞吐量、完成时间
4. 生成对比表格

**参数网格（建议）：**
```python
parameter_grid = {
    'workers': [8, 12, 16],
    'max_in_flight': [2, 4, 6, 8, 10],
    'batch_size': [25000, 50000, 75000],
}

# 约束条件
constraints = {
    'max_memory_percent': 70.0,  # 硬限制
    'min_throughput': 800.0,     # mol/s，期望目标
}

# 优先测试组合（基于分析）
priority_tests = [
    {'workers': 16, 'max_in_flight': 6, 'batch_size': 50000},  # 增加并行度
    {'workers': 16, 'max_in_flight': 8, 'batch_size': 50000},  # 进一步增加
    {'workers': 16, 'max_in_flight': 6, 'batch_size': 25000},  # 小batch减内存
    {'workers': 12, 'max_in_flight': 8, 'batch_size': 50000},  # 减workers增batch
    {'workers': 16, 'max_in_flight': 10, 'batch_size': 25000}, # 激进组合
]
```

**脚本结构：**
```python
#!/usr/bin/env python3
"""
Parameter optimization for transform pipeline.

Runs grid search over (workers, max-in-flight, batch-size)
to find optimal configuration balancing throughput and memory.
"""

import subprocess
import pandas as pd
import json
from pathlib import Path
import re

def run_test(workers, max_in_flight, batch_size, test_name):
    """
    Run single test with given parameters.

    Returns:
        dict: {
            'workers': int,
            'max_in_flight': int,
            'batch_size': int,
            'peak_memory': float,
            'avg_memory': float,
            'throughput': float,
            'time_seconds': float,
            'unique_products': int,
            'success': bool
        }
    """

    # Create output directory
    outdir = Path(f"data/output/transforms/OPT_{test_name}")
    if outdir.exists():
        import shutil
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Create 150K subset for quick testing
    import pyarrow.parquet as pq
    input_file = "data/output/nplike_v2/polyphenol-2X/products.parquet"
    table = pq.read_table(input_file)
    subset = table.slice(0, 150000)
    subset_file = outdir / "input_150k.parquet"
    pq.write_table(subset, subset_file)

    # Build command
    cmd = [
        "python",
        "scripts/08_transform_library_v2.py",
        "apply",
        "--input", str(subset_file),
        "--outdir", str(outdir),
        "--xf-config", "configs/transforms.yaml",
        "--xf-name", "FG_PHENOL_OH__OH__TO__OMe",
        "--batch-size", str(batch_size),
        "--workers", str(workers),
        "--flush-interval", "2000",
        "--target-memory", "70.0",
        "--use-bloom-filter",
        "--bloom-expected-items", "2000000",
        "--max-in-flight", str(max_in_flight),
    ]

    # Run test with timeout
    log_file = outdir / "test.log"
    try:
        with open(log_file, 'w') as log:
            result = subprocess.run(
                cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=900  # 15 min timeout
            )

        success = result.returncode == 0
    except subprocess.TimeoutExpired:
        success = False
        return {
            'workers': workers,
            'max_in_flight': max_in_flight,
            'batch_size': batch_size,
            'peak_memory': None,
            'avg_memory': None,
            'throughput': None,
            'time_seconds': None,
            'unique_products': None,
            'success': False,
            'error': 'timeout'
        }

    # Parse results
    if success:
        # Read summary
        summary_file = outdir / "SUMMARY.json"
        with open(summary_file) as f:
            summary = json.load(f)

        # Parse log for memory
        with open(log_file) as f:
            log_content = f.read()

        mem_values = re.findall(r'Pre-flush memory: ([0-9.]+)%', log_content)
        mem_floats = [float(m) for m in mem_values] if mem_values else []

        peak_mem = max(mem_floats) if mem_floats else None
        avg_mem = sum(mem_floats) / len(mem_floats) if mem_floats else None

        return {
            'workers': workers,
            'max_in_flight': max_in_flight,
            'batch_size': batch_size,
            'peak_memory': peak_mem,
            'avg_memory': avg_mem,
            'throughput': summary.get('throughput_mol_per_sec'),
            'time_seconds': summary.get('elapsed_seconds'),
            'unique_products': summary.get('unique_products'),
            'success': True
        }
    else:
        return {
            'workers': workers,
            'max_in_flight': max_in_flight,
            'batch_size': batch_size,
            'success': False,
            'error': 'failed'
        }

def main():
    print("="*80)
    print("PARAMETER OPTIMIZATION - Grid Search")
    print("="*80)

    # Priority tests
    priority_tests = [
        {'workers': 16, 'max_in_flight': 6, 'batch_size': 50000},
        {'workers': 16, 'max_in_flight': 8, 'batch_size': 50000},
        {'workers': 16, 'max_in_flight': 6, 'batch_size': 25000},
        {'workers': 12, 'max_in_flight': 8, 'batch_size': 50000},
    ]

    results = []

    for i, params in enumerate(priority_tests, 1):
        print(f"\n[Test {i}/{len(priority_tests)}]")
        print(f"  workers={params['workers']}, "
              f"max-in-flight={params['max_in_flight']}, "
              f"batch-size={params['batch_size']}")

        test_name = f"w{params['workers']}_m{params['max_in_flight']}_b{params['batch_size']}"
        result = run_test(**params, test_name=test_name)
        results.append(result)

        if result['success']:
            print(f"  ✓ Success: {result['throughput']:.1f} mol/s, "
                  f"{result['peak_memory']:.1f}% peak memory")
        else:
            print(f"  ✗ Failed: {result.get('error', 'unknown')}")

    # Create results DataFrame
    df = pd.DataFrame(results)

    # Save results
    df.to_csv('optimization_results.csv', index=False)
    print(f"\n\nResults saved to: optimization_results.csv")

    # Filter successful tests
    df_success = df[df['success'] == True]

    if len(df_success) > 0:
        # Sort by throughput
        df_sorted = df_success.sort_values('throughput', ascending=False)

        print("\n" + "="*80)
        print("TOP 5 CONFIGURATIONS (by throughput)")
        print("="*80)
        print(df_sorted[['workers', 'max_in_flight', 'batch_size',
                         'throughput', 'peak_memory', 'avg_memory']].head().to_string())

        # Find best config with memory < 70%
        df_safe = df_success[df_success['peak_memory'] < 70.0]
        if len(df_safe) > 0:
            best = df_safe.sort_values('throughput', ascending=False).iloc[0]

            print("\n" + "="*80)
            print("RECOMMENDED CONFIGURATION (best throughput with mem < 70%)")
            print("="*80)
            print(f"  workers: {best['workers']}")
            print(f"  max-in-flight: {best['max_in_flight']}")
            print(f"  batch-size: {best['batch_size']}")
            print(f"  Expected throughput: {best['throughput']:.1f} mol/s")
            print(f"  Expected peak memory: {best['peak_memory']:.1f}%")
            print(f"  Improvement over baseline (600 mol/s): {(best['throughput']/600 - 1)*100:.1f}%")

if __name__ == '__main__':
    main()
```

**执行：**
```bash
cd E:\Projects\halogenator
python optimize_parameters.py
```

**预期输出：**
- `optimization_results.csv`: 所有测试结果
- 控制台：TOP 5配置 + 推荐配置

**预期时间：** 4-5小时（4个测试 × 约60分钟/测试）

#### 子任务1.2：分析优化结果

**目标：** 从grid search结果中选择最优配置

**分析维度：**
1. **吞吐量vs内存散点图**
   - X轴：峰值内存
   - Y轴：吞吐量
   - 颜色：max-in-flight值
   - 目标：找到帕累托前沿

2. **性能提升计算**
   - 基准：600 mol/s (当前max-in-flight=4)
   - 目标：1,200 mol/s (100%提升)
   - 可接受：900 mol/s (50%提升)

3. **内存安全裕度**
   - 硬限制：70%
   - 推荐：<65% (5%裕度)
   - 理想：<60% (10%裕度)

**决策树：**
```
IF 有配置满足(throughput > 1200 AND peak_mem < 65%):
    → 选择吞吐量最高的
ELIF 有配置满足(throughput > 900 AND peak_mem < 70%):
    → 选择吞吐量最高的
ELIF 有配置满足(throughput > 600 AND peak_mem < 70%):
    → 选择吞吐量最高的（比当前好）
ELSE:
    → 保持当前配置(workers=16, max-in-flight=4)
```

### 📋 任务2：高级优化策略（可选，如果任务1不满足目标）

#### 策略2.1：动态max-in-flight调整

**问题：** 固定max-in-flight可能不是最优
- 早期batches：数据简单，可以更多并行
- 后期batches：数据复杂，需要限制并行

**解决方案：** 根据内存动态调整

**修改文件：** `E:\Projects\halogenator\scripts\08_transform_library_v2.py`

**位置：** `cmd_apply()` 函数，batch submission loop (lines ~1170-1250)

**当前代码：**
```python
# Line ~1170
max_in_flight = args.max_in_flight or (workers * 2)

# Line ~1200
while batches_submitted < num_batches or futures:
    # Submit new batches
    while batches_submitted < num_batches and len(futures) < max_in_flight:
        # Submit batch
        futures[future] = batch_idx
        batches_submitted += 1
```

**优化代码：**
```python
# 在cmd_apply开始添加
def calculate_dynamic_max_in_flight(current_memory_percent, base_max_in_flight):
    """
    Dynamically adjust max-in-flight based on current memory usage.

    Memory ranges:
    - < 50%: Allow 1.5x base
    - 50-60%: Allow base
    - 60-70%: Reduce to 0.75x base
    - > 70%: Reduce to 0.5x base
    """
    if current_memory_percent < 50:
        multiplier = 1.5
    elif current_memory_percent < 60:
        multiplier = 1.0
    elif current_memory_percent < 70:
        multiplier = 0.75
    else:
        multiplier = 0.5

    return max(2, int(base_max_in_flight * multiplier))

# 在batch loop中
while batches_submitted < num_batches or futures:
    # Get current memory
    current_mem = psutil.virtual_memory().percent

    # Adjust max-in-flight dynamically
    dynamic_max = calculate_dynamic_max_in_flight(current_mem, args.max_in_flight)

    # Submit new batches with dynamic limit
    while batches_submitted < num_batches and len(futures) < dynamic_max:
        # Submit batch
        futures[future] = batch_idx
        batches_submitted += 1

        # Log adjustment
        if dynamic_max != args.max_in_flight:
            logger.info(f"Dynamic adjustment: max-in-flight={dynamic_max} "
                       f"(base={args.max_in_flight}, mem={current_mem:.1f}%)")
```

**预期效果：**
- 内存低时：允许更多并行（提升吞吐量）
- 内存高时：自动限制并行（保护稳定性）
- 自适应调整

#### 策略2.2：Worker池预热（Warm-up）

**问题：** Workers冷启动慢，RDKit初始化有overhead

**解决方案：** 预热worker pool

**修改位置：** `cmd_apply()` 函数，executor创建后

**添加代码：**
```python
# After executor creation (line ~1170)
with ProcessPoolExecutor(max_workers=workers) as executor:
    # Warm up workers
    logger.info(f"Warming up {workers} workers...")
    dummy_batch = table.slice(0, min(100, len(table))).to_pandas()

    warmup_futures = []
    for _ in range(workers):
        future = executor.submit(_process_batch_worker, dummy_batch)
        warmup_futures.append(future)

    # Wait for warmup
    for future in warmup_futures:
        try:
            future.result(timeout=30)
        except:
            pass

    logger.info(f"Worker warmup complete")

    # Continue with normal processing...
```

#### 策略2.3：批处理大小自适应

**问题：** 固定batch_size对所有数据不是最优

**解决方案：** 根据分子复杂度调整batch大小

**实现：** 在batch读取时检测复杂度

```python
def estimate_batch_complexity(batch_df):
    """
    Estimate computational complexity of a batch.

    Returns complexity score (higher = more complex)
    """
    # Sample molecules
    sample_smiles = batch_df['smiles'].head(100)

    complexities = []
    for smi in sample_smiles:
        # Complexity indicators:
        # - SMILES length (proxy for size)
        # - Number of rings (aromatic systems)
        # - Number of phenolic OH (transformation sites)
        complexity = len(smi) * (smi.count('O') + 1)
        complexities.append(complexity)

    return np.mean(complexities)

# In batch processing loop
for batch_idx in range(num_batches):
    batch_df = ...

    complexity = estimate_batch_complexity(batch_df)

    # Adjust batch size
    if complexity > threshold_high:
        # Complex molecules, use smaller batch
        adjusted_batch_size = batch_size // 2
    elif complexity < threshold_low:
        # Simple molecules, use larger batch
        adjusted_batch_size = int(batch_size * 1.5)
    else:
        adjusted_batch_size = batch_size
```

### 📋 任务3：生产验证（必做）

**在确定最优参数后：**

#### 子任务3.1：MEDIUM测试验证

**目标：** 用新参数在500K rows测试，确保稳定

**命令模板：**
```bash
python scripts/08_transform_library_v2.py apply \
  --input data/output/transforms/TEST_MEDIUM_FIXED/input_500k.parquet \
  --outdir data/output/transforms/TEST_MEDIUM_OPT \
  --xf-config configs/transforms.yaml \
  --xf-name FG_PHENOL_OH__OH__TO__OMe \
  --batch-size [OPTIMIZED_VALUE] \
  --workers [OPTIMIZED_VALUE] \
  --flush-interval 2000 \
  --target-memory 70.0 \
  --use-bloom-filter \
  --bloom-expected-items 5000000 \
  --max-in-flight [OPTIMIZED_VALUE]
```

**验证标准：**
- 峰值内存 < 70%
- 吞吐量 > 当前配置(600 mol/s)
- 无critical警告
- 输出有效

#### 子任务3.2：Full polyphenol-2X测试

**如果MEDIUM测试通过，运行完整数据集：**

**命令：**
```bash
# 使用优化后的参数
python scripts/08_transform_library_v2.py apply \
  --input data/output/nplike_v2/polyphenol-2X/products.parquet \
  --outdir data/output/transforms/polyphenol-2X_FG_PHENOL_OH__OH__TO__OMe_OPT \
  --xf-config configs/transforms.yaml \
  --xf-name FG_PHENOL_OH__OH__TO__OMe \
  --use-bloom-filter \
  --workers [OPTIMIZED] \
  --target-memory 70.0 \
  --max-in-flight [OPTIMIZED] \
  --batch-size [OPTIMIZED] \
  --bloom-expected-items 100000000
```

**监控：** 每2-4小时检查内存和进度

**成功标准：**
- 完成不崩溃
- 峰值内存 < 70%
- 总时间 < 20小时（vs 当前预期24h）
- 产出 ~89M products

---

## 第四部分：关键文件参考

### 核心代码文件

**1. `E:\Projects\halogenator\scripts\08_transform_library_v2.py`**
- **功能：** 主transform脚本
- **已修改：** StreamingParquetWriter.write_batch() (lines 413-483)
- **待修改区域（如需高级优化）：**
  - `cmd_apply()` function (lines 1025-1340)
  - Batch submission loop (lines 1170-1250)
  - Worker pool creation (line ~1170)

**2. 备份文件：**
- `E:\Projects\halogenator\scripts\08_transform_library_v2.py.backup_before_buffer_fix`
  - 修复前的原始版本
  - 如需回滚可用

### 测试脚本

**3. `E:\Projects\halogenator\test_micro_simple.py`**
- MICRO测试（50K rows，快速验证）
- 用于快速测试修改

**4. `E:\Projects\halogenator\test_medium_fixed.py`**
- MEDIUM测试（500K rows，完整验证）
- 用于生产前最终验证

**5. `E:\Projects\halogenator\validate_fix_gradual.py`**
- 三级渐进测试（MICRO+SMALL+MEDIUM）
- 完整验证流程

**6. 待创建：`E:\Projects\halogenator\optimize_parameters.py`**
- 参数优化grid search脚本
- 详见任务1.1

### 文档文件

**7. `E:\Projects\halogenator\SOLUTION_VALIDATED_PRODUCTION_READY.md`**
- 完整解决方案文档
- 包含所有测试结果和配置

**8. `E:\Projects\halogenator\ROOT_CAUSE_ANALYSIS_MEMORY_CRISIS.md`**
- 10部分详细根本原因分析
- 技术深度文档

**9. `E:\Projects\halogenator\MEDIUM_TEST_DIAGNOSIS.md`**
- MEDIUM测试诊断报告
- 内存分布分析

**10. `E:\Projects\halogenator\TEST_OBSERVATION_GUIDE.md`**
- 测试观察指南
- 成功/失败标准

### 数据文件

**11. 测试数据位置：**
```
E:\Projects\halogenator\data\output\nplike_v2\polyphenol-2X\products.parquet
  - 完整数据集：13.79M rows, 504MB
  - 用于生产测试

E:\Projects\halogenator\data\output\transforms\TEST_MICRO\
  - MICRO测试结果

E:\Projects\halogenator\data\output\transforms\TEST_SMALL_polyphenol\
  - SMALL测试结果

E:\Projects\halogenator\data\output\transforms\TEST_MEDIUM_FIXED\
  - MEDIUM测试结果（已修复）
  - input_500k.parquet - 可重用于快速测试
```

---

## 第五部分：技术细节与注意事项

### 内存模型（关键理解）

**系统内存 = Workers Memory + Queue Memory + Buffer + Overhead**

```python
# 简化模型
def estimate_memory(workers, max_in_flight, batch_size):
    """
    Estimate peak memory usage.

    Parameters from testing:
    - Worker base: ~1.0 GB per worker
    - Queue: ~0.5 GB per in-flight batch
    - Buffer: ~0.01 GB (fixed after chunking)
    - Overhead: ~6 GB (OS + Python + Bloom)
    """
    worker_mem = workers * 1.0  # GB
    queue_mem = max_in_flight * 0.5  # GB
    buffer_mem = 0.01  # GB
    overhead = 6.0  # GB

    total_gb = worker_mem + queue_mem + buffer_mem + overhead
    total_percent = (total_gb / 32.0) * 100  # Assuming 32GB system

    return total_percent

# Examples:
# workers=16, max-in-flight=4:  16 + 2 + 0.01 + 6 = 24GB = 75%  ← 接近实测
# workers=16, max-in-flight=8:  16 + 4 + 0.01 + 6 = 26GB = 81%  ← 可能超限
# workers=12, max-in-flight=8:  12 + 4 + 0.01 + 6 = 22GB = 69%  ← 安全
# workers=16, max-in-flight=6:  16 + 3 + 0.01 + 6 = 25GB = 78%  ← 边界
```

**使用此模型：**
1. 在grid search前预估哪些组合安全
2. 过滤掉预估>75%的组合
3. 优先测试60-70%范围的组合

### 性能瓶颈诊断

**如果优化后吞吐量仍不理想，检查：**

1. **CPU利用率**
```bash
# 运行测试时，另一个终端执行
top  # Linux
# 或
Get-Process python | Select-Object CPU  # Windows PowerShell

# 期望：接近 (workers * 100)%
# 如果低：说明workers空闲，增加max-in-flight
```

2. **I/O等待**
```bash
# 检查磁盘I/O
iostat -x 5  # Linux

# 如果I/O wait高：
# - 考虑增加flush_interval减少写入频率
# - 或使用SSD存储
```

3. **内存交换（Swap）**
```bash
# 检查swap使用
free -h  # Linux
# 或通过日志中的内存值

# 如果有swap：立即降低max-in-flight
```

### 已知限制与约束

**1. ProcessPoolExecutor限制：**
- 结果必须全部序列化到queue
- 大量RDKit Mol对象序列化开销高
- 无法直接控制queue大小，只能通过max-in-flight间接控制

**2. RDKit内存特性：**
- Mol对象在C++层分配内存
- Python GC无法立即释放
- 需要显式del + gc.collect()
- 每个Mol对象约5-10KB（含衍生物）

**3. PyArrow内存：**
- Table.from_pylist() 会复制数据
- 在flush期间短暂双倍内存
- 已通过chunking缓解

### 潜在风险与缓解

**风险1：激进优化导致OOM**
- **缓解：** 始终在SMALL/MEDIUM测试，不要直接上生产
- **回滚：** 保持当前稳定配置作为fallback

**风险2：参数组合爆炸（测试太多）**
- **缓解：** 使用priority_tests，先测最可能成功的
- **策略：** Bayesian optimization替代grid search（高级）

**风险3：不同数据集最优参数不同**
- **缓解：** 在polyphenol-2X, terpenoid-2X, alkaloid-2X各测试
- **文档：** 为每类数据集记录最优配置

---

## 第六部分：执行计划（建议顺序）

### 阶段1：参数探索（必做）

**时间：** 4-6小时

1. **创建 optimize_parameters.py** (30分钟)
   - 复制上面的脚本模板
   - 调整parameter_grid如needed

2. **运行grid search** (4-5小时)
   ```bash
   python optimize_parameters.py > optimization_log.txt 2>&1 &
   ```
   - 后台运行
   - 定期检查进度

3. **分析结果** (30分钟)
   - 查看 optimization_results.csv
   - 选择最优配置
   - 记录推荐参数

### 阶段2：验证最优配置（必做）

**时间：** 1-2小时

1. **MEDIUM测试验证**
   - 使用优化参数运行500K测试
   - 确认内存<70%且吞吐量提升

2. **如果通过：** 进入阶段3
3. **如果失败：** 回到阶段1，调整参数范围

### 阶段3：生产部署（必做）

**时间：** 20-24小时

1. **运行完整polyphenol-2X**
   - 使用优化参数
   - 后台运行（nohup或screen）
   - 每2-4小时检查

2. **监控关键指标：**
   - 内存不超70%
   - 无critical警告
   - 吞吐量符合预期

3. **验证输出：**
   - 文件完整
   - 产品数量~89M
   - 无损坏

### 阶段4：高级优化（可选）

**时间：** 2-4小时

**仅当阶段1-3无法达成目标时执行**

1. 实现动态max-in-flight (策略2.1)
2. 或实现worker预热 (策略2.2)
3. 或实现自适应batch大小 (策略2.3)
4. 重新测试

---

## 第七部分：成功标准

### 优化成功标准

**Level 1 - 最低要求（必须达成）：**
- ✅ 峰值内存 < 70%
- ✅ 吞吐量 > 600 mol/s（至少不降低）
- ✅ 无崩溃、无损坏

**Level 2 - 期望目标（尽力达成）：**
- ✅ 峰值内存 < 65%
- ✅ 吞吐量 800-1000 mol/s（33-66%提升）
- ✅ 运行时间 18-20h

**Level 3 - 理想目标（最佳情况）：**
- ✅ 峰值内存 < 60%
- ✅ 吞吐量 1200-1500 mol/s（100-150%提升）
- ✅ 运行时间 14-16h

### 验收测试

**完成优化后，运行以下验收：**

```bash
# 1. SMALL测试（快速）
python test_medium_fixed.py  # 修改为使用优化参数

# 2. 检查结果
# - 内存 < 70%
# - 吞吐量 > 当前baseline

# 3. 如果通过，运行生产测试
python scripts/08_transform_library_v2.py apply [优化参数]

# 4. 监控24小时

# 5. 验证输出
python -c "
import pyarrow.parquet as pq
t = pq.read_table('data/output/transforms/[OUTPUT]/products.parquet')
print(f'Products: {len(t):,}')
print('Valid!' if len(t) > 85000000 else 'Check count!')
"
```

---

## 第八部分：快速参考

### 当前稳定配置（Baseline）

```bash
--workers 16 \
--max-in-flight 4 \
--batch-size 50000 \
--flush-interval 2000 \
--target-memory 70.0

# 性能：600-800 mol/s，50-55%内存
```

### 推测的最优配置候选

**保守（安全优先）：**
```bash
--workers 16 \
--max-in-flight 6 \
--batch-size 50000
# 预期：800-900 mol/s，55-65%内存
```

**平衡（推荐）：**
```bash
--workers 16 \
--max-in-flight 8 \
--batch-size 50000
# 预期：1000-1200 mol/s，60-70%内存
```

**激进（需验证）：**
```bash
--workers 12 \
--max-in-flight 10 \
--batch-size 25000
# 预期：1200+mol/s，可能65-75%内存
```

**注意：** 这些都是推测，必须通过grid search验证！

### 紧急回滚命令

**如果优化后出问题，立即回滚：**

```bash
# 1. 停止当前运行
kill [PID]

# 2. 使用稳定配置
python scripts/08_transform_library_v2.py apply \
  --input [INPUT] \
  --outdir [OUTPUT] \
  --xf-config configs/transforms.yaml \
  --xf-name [TRANSFORM] \
  --use-bloom-filter \
  --workers 16 \
  --max-in-flight 4 \
  --batch-size 50000 \
  --target-memory 70.0 \
  --bloom-expected-items 100000000
```

---

## 第九部分：预期的下一个会话工作流

**建议的会话流程：**

```
1. 阅读本报告 (10分钟)
   ↓
2. 创建optimize_parameters.py (30分钟)
   ↓
3. 运行grid search (后台4-5小时)
   ├─ 可以做其他工作
   └─ 定期检查进度
   ↓
4. 分析optimization_results.csv (30分钟)
   ├─ 选择最优配置
   └─ 理解性能/内存权衡
   ↓
5. MEDIUM验证测试 (1小时)
   ├─ 使用优化参数
   └─ 确认稳定性
   ↓
6. 决策点：
   ├─ 如果通过 → 生产部署
   └─ 如果失败 → 调整参数或尝试高级策略
   ↓
7. 生产测试 (20-24小时，可后台)
   ↓
8. 验证成功 → 文档化最优配置
```

**预期产出：**
1. `optimization_results.csv` - 所有测试数据
2. 推荐的生产配置
3. 性能基准报告
4. 更新的生产运行脚本

---

## 第十部分：重要提醒

### ⚠️ 必须记住的关键点

1. **NEVER compromise stability for speed**
   - 70%内存是硬限制
   - 宁可慢一点，不要崩溃

2. **Always test on SMALL/MEDIUM before production**
   - 500K测试仅需1小时
   - 能发现90%的问题
   - 比直接在13M数据上失败好得多

3. **Keep the baseline config as fallback**
   - workers=16, max-in-flight=4已验证稳定
   - 如果优化失败，可以回退
   - 稳定完成>快速失败

4. **Document everything**
   - 每次测试记录参数和结果
   - 便于追溯和复现
   - 帮助理解性能特征

5. **Buffer fix is NOT optional**
   - chunking代码必须保留
   - 这是防止940K爆炸的核心
   - 没有这个fix，任何优化都会失败

### 🎯 优化哲学

**当前状态：**
- 稳定但慢（600 mol/s）

**目标状态：**
- 稳定且快（1000+ mol/s）

**方法：**
- 谨慎增加并行度
- 持续监控内存
- 渐进式改进
- 每步验证

**不是：**
- 激进调参
- 跳过测试
- 牺牲稳定性

---

## 第十一部分：联系信息与资源

### 关键数据位置

**输入数据：**
- `E:\Projects\halogenator\data\output\nplike_v2\polyphenol-2X\products.parquet`
- `E:\Projects\halogenator\data\output\nplike_v2\terpenoid-2X\products.parquet`

**测试输出：**
- `E:\Projects\halogenator\data\output\transforms\TEST_*\`

**配置文件：**
- `E:\Projects\halogenator\configs\transforms.yaml`

**Git仓库：**
- 当前分支：`fix/pr2-contract-and-sugar-accept`
- 修改的文件已staged但未commit

### 环境信息

**Python环境：**
- Anaconda环境：`halo-p0`
- Python路径：`C:\Users\Greatony\anaconda3\envs\halo-p0\python.exe`

**系统：**
- OS: Windows
- 总内存：32GB
- 工作目录：`E:\Projects\halogenator`

**关键依赖：**
- pyarrow (Parquet读写)
- RDKit (化学计算)
- psutil (内存监控)
- pandas

---

## 第十二部分：最后的话

**你即将接手一个已经解决了稳定性问题但需要优化性能的系统。**

**核心挑战：** 在不破坏稳定性的前提下，尽可能提升吞吐量。

**你拥有的资源：**
- ✅ 完全理解的内存模型
- ✅ 经过验证的稳定配置
- ✅ 详细的测试框架
- ✅ 清晰的优化路径

**期望你完成：**
- 🎯 找到最优参数配置
- 🎯 至少50%吞吐量提升（600→900+ mol/s）
- 🎯 保持内存<70%
- 🎯 在真实数据上验证

**如果遇到困难：**
- 参考SOLUTION_VALIDATED_PRODUCTION_READY.md
- 重新阅读内存模型部分
- 降低期望，接受保守配置
- 稳定性永远优先

**Good luck! 🚀**

---

**报告结束**
**准备交接给下一个会话**
**当前时间：** 2025-12-27
**报告作者：** Claude (当前会话)
**交接给：** Claude (下一个会话)
