# Transform Pipeline Session Handoff Report
**Date:** 2025-12-20
**Project:** Halogenator Transform Library Generation
**Status:** Pipeline暂停，内存优化待实施
**Next Session Goal:** 修复内存管理问题并完成剩余90个transform jobs

---

## 📋 Executive Summary

### Current Situation
- ✅ **已完成:** 28个transform jobs，生成9.5M products
- ❌ **失败:** 2个大型jobs因内存耗尽僵死
- ⏸️ **暂停:** Pipeline已终止，等待内存优化后重启
- 🎯 **目标:** 完成剩余90个jobs，生成完整transform library

### Critical Issue
**内存占用超过90-95%** 导致系统swap，进程僵死。当前flush机制基于product count，不考虑实际内存使用，在大型jobs上失效。

### User Requirements
1. **内存目标:** 稳定在50%左右（当前90-95%）
2. **速度要求:** Workers=16（不降低）
3. **可接受trade-off:** 用磁盘空间换内存（更频繁flush）

---

## ✅ Part A: 已完成任务总结

### 1. StreamingParquetWriter改进 (完成)
**文件:** `scripts/08_transform_library_v2.py` (lines 301-444)

**已实施改进:**
- ✅ Product buffer + flush_interval机制
- ✅ 可选psutil内存监控
- ✅ Schema enforcement防止null-type错误
- ✅ CLI参数 `--flush-interval` (default: 10000)
- ✅ 超时从1小时增加到24小时 (line 173: `timeout=86400`)

**问题:** flush_interval基于product count，不考虑实际内存/产品大小，导致内存爆炸。

### 2. Python自动化脚本创建 (完成)
**文件:** `scripts/run_transform_pipeline.py`

**功能:**
- 自动解析transforms.yaml构建class-rule映射
- Resume机制（跳过已完成jobs）
- 进度追踪和日志记录
- 错误处理

**问题:** 没有主动内存管理，依赖OS swap（导致僵死）。

### 3. 已完成Transform Jobs (28个)
**总products:** 9,509,456

**Breakdown:**
- **Lipid class (12 jobs):** 315,816 products
  - lipid-1X: 6 jobs, 24,264 products
  - lipid-2X: 6 jobs, 291,552 products

- **Polyphenol-1X (16 jobs):** 9,193,640 products
  - Phenolic OH (4 rules): 9,052,332 products
  - Aromatic Amine (8 rules): 348 products
  - Carboxyl (3 rules): 140,916 products

**存储位置:** `E:/Projects/halogenator/data/output/transforms/`

### 4. 失败Jobs诊断 (已分析)
**失败jobs:** 2个
- Job 17: polyphenol-2X_FG_PHENOL_OH__OH__TO__OMe
- Job 18: polyphenol-2X_FG_PHENOL_OH__OH__TO__NH2

**失败原因:**
1. 输入数据巨大 (13.79M products, 504MB)
2. Dedup数据库膨胀 (16GB+)
3. 内存占用90-95%
4. 系统swap → 进程僵死

**损坏文件位置:**
```
E:/Projects/halogenator/data/output/transforms/polyphenol-2X_FG_PHENOL_OH__OH__TO__OMe/
E:/Projects/halogenator/data/output/transforms/polyphenol-2X_FG_PHENOL_OH__OH__TO__NH2/
```

---

## 🔴 Part B: 当前关键问题详细分析

### 问题1: Flush机制缺陷

**当前实现 (08_transform_library_v2.py, lines 354-378):**
```python
def write_batch(self, products: List[Dict]):
    if not products:
        return

    # Add to buffer
    self.buffer.extend(products)

    # Check flush conditions
    mem_usage = self._check_memory()
    should_flush = (
        len(self.buffer) >= self.flush_interval or  # ← 只看product数量
        (self.memory_monitor_enabled and mem_usage > self.memory_threshold)  # ← 全局内存，不准确
    )
```

**缺陷:**
1. **flush_interval基于product count:**
   - 小product (简单分子): 1KB/product
   - 大product (复杂NP衍生物): 10-50KB/product
   - 同样10,000 products，内存占用可能相差50倍！

2. **全局内存监控不准确:**
   - `psutil.virtual_memory().percent` 包括所有进程
   - 无法区分当前Python进程的实际占用
   - 16个workers同时运行，无法预测峰值

3. **Dedup数据库无限增长:**
   - SQLite dedup.db在job内持续增长
   - Job 17/18的dedup.db达到16GB
   - 不会自动commit/清理，直到job结束

### 问题2: 内存占用计算错误

**观察到的现象:**
- 16 workers × 大型分子 × 10,000 buffer = 内存爆炸
- Polyphenol-2X jobs: 每个product ~5-10KB
- 16 workers × 10,000 products × 8KB = **1.28GB** 仅buffer
- 加上dedup.db (16GB) + 其他开销 = **20GB+**

**用户系统配置推测:**
- 总内存: ~32GB (90% = 28.8GB, 50% = 16GB)
- 期望占用: 16GB
- 当前实际: 28-30GB（触发swap）

### 问题3: Workers并发加剧问题

**当前配置:**
- Workers: 16 (用户要求保持)
- 每个worker独立buffer
- 峰值内存 = 16 × (buffer + overhead)

**问题:**
- 所有workers可能同时flush → 内存峰值
- 没有worker间协调机制
- 没有全局内存预算分配

---

## 🎯 Part C: 待完成任务详细方案

### Task 1: 实现智能内存管理的StreamingParquetWriter (🔴 CRITICAL)

**目标:** 将内存占用稳定在50% (~16GB)，同时保持workers=16

**实施方案:**

#### 1.1 添加进程级内存监控

**修改文件:** `scripts/08_transform_library_v2.py`

**位置:** StreamingParquetWriter类 (lines 301-444)

**新增方法:**
```python
import psutil
import os

def _get_process_memory_mb(self) -> float:
    """获取当前Python进程的内存占用（MB）"""
    try:
        process = psutil.Process(os.getpid())
        return process.memory_info().rss / 1024 / 1024  # Bytes to MB
    except:
        return 0.0

def _get_total_memory_mb(self) -> float:
    """获取系统总内存（MB）"""
    try:
        return psutil.virtual_memory().total / 1024 / 1024
    except:
        return 32000.0  # 假设32GB

def _get_memory_usage_percent(self) -> float:
    """获取当前进程内存占用百分比"""
    process_mem = self._get_process_memory_mb()
    total_mem = self._get_total_memory_mb()
    return (process_mem / total_mem) * 100
```

**插入位置:** 在`_check_memory()`方法之后 (line 352后)

#### 1.2 实现自适应flush策略

**目标:** 根据实际内存使用动态调整flush频率

**修改位置:** `write_batch()` 方法 (lines 354-378)

**新实现:**
```python
def write_batch(self, products: List[Dict]):
    """
    智能内存管理的batch写入

    Strategy:
    1. 监控进程内存占用
    2. 当达到阈值时强制flush
    3. 使用自适应buffer大小
    """
    if not products:
        return

    # Add to buffer
    self.buffer.extend(products)

    # Get current memory status
    process_mem_percent = self._get_memory_usage_percent()
    buffer_size = len(self.buffer)

    # Adaptive flush conditions
    should_flush = False
    flush_reason = ""

    # Condition 1: Memory pressure (PRIORITY)
    if process_mem_percent > 45.0:  # 接近50%目标时开始flush
        should_flush = True
        flush_reason = f"memory_pressure ({process_mem_percent:.1f}%)"

    # Condition 2: Buffer size limit (backup)
    elif buffer_size >= self.flush_interval:
        should_flush = True
        flush_reason = f"buffer_full ({buffer_size:,} products)"

    # Condition 3: Critical memory (URGENT)
    if process_mem_percent > 55.0:  # 紧急flush
        should_flush = True
        flush_reason = f"CRITICAL_MEMORY ({process_mem_percent:.1f}%)"
        self.logger.warning(f"⚠️  Critical memory usage: {process_mem_percent:.1f}%")

    if should_flush:
        self.logger.info(f"Flush triggered: {flush_reason}")
        self._flush()
```

**关键参数:**
- 45%: 开始主动flush（给worker峰值留余量）
- 55%: 紧急flush（防止达到目标50%上限）
- 保留flush_interval作为backup（防止内存监控失效）

#### 1.3 优化_flush()方法

**添加功能:**
1. Flush后立即commit dedup数据库
2. 清理临时对象
3. 可选：显式触发GC

**修改位置:** `_flush()` 方法 (lines 380-429)

**在现有代码最后添加:**
```python
def _flush(self):
    """现有代码保持不变，在最后添加:"""

    # 现有flush逻辑 (lines 387-428) ...

    # Clear buffer
    self.buffer.clear()

    # NEW: 主动内存管理
    import gc
    gc.collect()  # 强制垃圾回收

    # Log memory status after flush
    mem_after = self._get_memory_usage_percent()
    self.logger.info(f"Memory after flush: {mem_after:.1f}%")
```

#### 1.4 修改__init__参数

**位置:** `__init__()` 方法 (lines 312-334)

**修改:**
```python
def __init__(
    self,
    output_path: Path,
    schema: pa.Schema,
    flush_interval: int = 10000,
    memory_threshold: float = 0.50,  # 改为50% (之前是0.8)
    target_memory_percent: float = 45.0  # NEW: 目标内存占用
):
    # 现有代码...
    self.target_memory_percent = target_memory_percent
```

#### 1.5 更新CLI参数

**文件:** `scripts/run_transform_pipeline.py`

**位置:** argparse部分 (搜索 `--flush-interval`)

**添加新参数:**
```python
parser_apply.add_argument(
    '--target-memory',
    type=float,
    default=45.0,
    help='Target process memory usage percent (default: 45.0). '
         'Flush will be triggered when approaching this value.'
)
```

**修改writer初始化 (搜索 `StreamingParquetWriter(`):**
```python
writer = StreamingParquetWriter(
    products_path,
    output_schema,
    flush_interval=args.flush_interval,
    target_memory_percent=args.target_memory  # NEW
)
```

---

### Task 2: 优化Dedup数据库管理 (🟡 IMPORTANT)

**问题:** Dedup.db在大型jobs中膨胀到16GB，占用大量内存

**解决方案1: 定期commit和checkpoint**

**文件:** `scripts/08_transform_library_v2.py`

**位置:** 找到deduplication相关代码（搜索`deduper`或`DedupManager`）

**修改deduplication逻辑:**
```python
# 在处理batch的循环中 (搜索: deduper.filter_new_keys)
# 每N个batch后执行checkpoint

batch_count = 0
CHECKPOINT_INTERVAL = 10  # 每10个batch checkpoint一次

for batch_idx in range(num_batches):
    # 现有处理逻辑...

    # Deduplication
    new_keys = deduper.filter_new_keys(canon_smiles_list)
    # ...
    deduper.commit()

    batch_count += 1

    # NEW: Periodic checkpoint
    if batch_count % CHECKPOINT_INTERVAL == 0:
        # SQLite checkpoint (如果使用SQLite)
        if hasattr(deduper, 'conn'):
            deduper.conn.execute('PRAGMA wal_checkpoint(TRUNCATE)')
            logger.info(f"Dedup DB checkpoint at batch {batch_idx}")
```

**解决方案2: 使用内存限制的dedup策略**

**更激进方案（如果方案1不够）:**
- 定期清空dedup数据库，只保留最近N条记录
- 或使用Bloom filter替代完整dedup（trade-off: 可能有极少重复）

---

### Task 3: 实现Worker协调机制 (🟢 OPTIONAL但推荐)

**目标:** 避免16个workers同时flush导致内存峰值

**方案:** 添加全局flush coordinator

**新建文件:** `scripts/flush_coordinator.py`

```python
import threading
import time

class FlushCoordinator:
    """
    协调多个workers的flush操作，避免同时flush
    使用token bucket算法控制并发flush数量
    """

    def __init__(self, max_concurrent_flush: int = 4):
        self.max_concurrent = max_concurrent_flush
        self.semaphore = threading.Semaphore(max_concurrent_flush)
        self.flush_count = 0
        self.lock = threading.Lock()

    def acquire_flush_permission(self, worker_id: int) -> bool:
        """
        请求flush权限

        Returns:
            True if flush granted, False if should wait
        """
        acquired = self.semaphore.acquire(blocking=False)
        if acquired:
            with self.lock:
                self.flush_count += 1
        return acquired

    def release_flush_permission(self):
        """释放flush权限"""
        self.semaphore.release()
        with self.lock:
            self.flush_count -= 1

    def wait_for_flush_slot(self, timeout: float = 30.0):
        """等待flush槽位"""
        return self.semaphore.acquire(timeout=timeout)
```

**集成到StreamingParquetWriter:**

在`_flush()`方法开始处:
```python
def _flush(self):
    if not self.buffer:
        return

    # NEW: Wait for flush permission (if coordinator available)
    if hasattr(self, 'flush_coordinator') and self.flush_coordinator:
        self.flush_coordinator.wait_for_flush_slot()

    try:
        # 现有flush逻辑...
    finally:
        # Release permission
        if hasattr(self, 'flush_coordinator') and self.flush_coordinator:
            self.flush_coordinator.release_flush_permission()
```

---

### Task 4: 清理失败Jobs并重启Pipeline

#### 4.1 清理损坏文件

**命令:**
```bash
cd /e/Projects/halogenator

# 删除2个失败jobs (~35GB)
rm -rf data/output/transforms/polyphenol-2X_FG_PHENOL_OH__OH__TO__OMe
rm -rf data/output/transforms/polyphenol-2X_FG_PHENOL_OH__OH__TO__NH2

# 验证已删除
ls data/output/transforms/ | wc -l
# 应该返回 28
```

#### 4.2 验证已完成jobs

**Python验证脚本:**
```python
import pyarrow.parquet as pq
from pathlib import Path

transform_dir = Path('E:/Projects/halogenator/data/output/transforms')
completed = []
total_products = 0

for job_dir in sorted(transform_dir.iterdir()):
    if not job_dir.is_dir():
        continue

    parquet_file = job_dir / 'products.parquet'
    if not parquet_file.exists():
        continue

    try:
        table = pq.read_table(parquet_file)
        count = table.num_rows
        completed.append((job_dir.name, count))
        total_products += count
    except Exception as e:
        print(f"ERROR: {job_dir.name} - {e}")

print(f"\n✓ Completed jobs: {len(completed)}")
print(f"✓ Total products: {total_products:,}")
print(f"\nExpected: 28 jobs, 9,509,456 products")
```

#### 4.3 修改后重启Pipeline

**命令:**
```bash
cd /e/Projects/halogenator

# 重启pipeline with新参数
nohup python scripts/run_transform_pipeline.py \
  --workers 16 \
  --flush-interval 5000 \
  --batch-size 50000 \
  --target-memory 45.0 \
  > transform_pipeline_optimized.log 2>&1 &

# 记录PID
echo $!

# 监控启动
tail -f transform_pipeline_optimized.log
```

**关键参数变化:**
- `--workers 16`: 保持不变（用户要求）
- `--flush-interval 5000`: 从10000降低到5000（更频繁flush）
- `--target-memory 45.0`: 新参数，目标45%内存

---

## 📊 Part D: 剩余工作详细清单

### Jobs待执行 (90个)

| Class | K | Rules | Jobs | Est. Products | Est. Time |
|-------|---|-------|------|---------------|-----------|
| Polyphenol-2X | 2X | 16 | 16 | ~250M | 40-60h |
| Terpenoid-1X | 1X | 10 | 10 | ~15M | 10-15h |
| Terpenoid-2X | 2X | 10 | 10 | ~400M | 80-120h |
| Alkaloid-1X | 1X | 12 | 12 | ~3M | 5-8h |
| Alkaloid-2X | 2X | 12 | 12 | ~80M | 20-30h |
| AA_peptide-1X | 1X | 15 | 15 | ~600K | 3-5h |
| AA_peptide-2X | 2X | 15 | 15 | ~8M | 8-12h |
| **Total** | - | - | **90** | **~757M** | **166-250h** |

**预计总时长:** 7-10天（with workers=16, 优化后）

### 特别注意的Jobs

**高风险jobs (需要特别监控):**
1. **Polyphenol-2X phenolic OH (4 jobs):**
   - 输入: 13.79M products each
   - 预计输出: ~60M products each
   - 历史问题: 内存耗尽
   - **监控:** 每2小时检查一次内存

2. **Terpenoid-2X (10 jobs):**
   - 输入: 40.17M products (最大)
   - 预计输出: ~400M products total
   - 风险: 与polyphenol-2X类似
   - **建议:** 可能需要单独处理或分批

**快速jobs (优先执行验证优化效果):**
- AA_peptide-1X (15 jobs, ~3-5小时)
- Alkaloid-1X (12 jobs, ~5-8小时)

**建议执行顺序:**
1. 先跑aa_peptide-1X验证内存优化有效
2. 再跑polyphenol-2X测试大型jobs
3. 最后跑terpenoid-2X

---

## 🗂️ Part E: 关键文件路径速查

### 核心脚本
```
E:/Projects/halogenator/scripts/08_transform_library_v2.py
  - StreamingParquetWriter类: lines 301-444
  - write_batch()方法: lines 354-378
  - _flush()方法: lines 380-429
  - 超时设置: line 173 (timeout=86400)

E:/Projects/halogenator/scripts/run_transform_pipeline.py
  - 主orchestrator脚本
  - 自动resume机制
  - Class-rule映射: build_class_rule_mapping()
```

### 配置文件
```
E:/Projects/halogenator/configs/transforms.yaml
  - 19个transform rules定义
  - scope_classes映射
```

### 数据目录
```
输入 (halogenated library):
E:/Projects/halogenator/data/output/nplike_v2/
  - aa_peptide-1X/products.parquet (30,564 products)
  - aa_peptide-2X/products.parquet (582,730 products)
  - alkaloid-1X/products.parquet (322,280 products)
  - alkaloid-2X/products.parquet (8,461,937 products)
  - lipid-1X/products.parquet (33,084 products)
  - lipid-2X/products.parquet (180,224 products)
  - polyphenol-1X/products.parquet (518,852 products)
  - polyphenol-2X/products.parquet (13,790,820 products) ← 最大
  - terpenoid-1X/products.parquet (1,528,400 products)
  - terpenoid-2X/products.parquet (40,170,816 products) ← 超大
  - polysaccharide-1X/products.parquet (20,524 products)
  - polysaccharide-2X/products.parquet (293,616 products)
  - other-1X/products.parquet (97,012 products)
  - other-2X/products.parquet (1,419,941 products)

输出 (transform library):
E:/Projects/halogenator/data/output/transforms/
  - 当前28个已完成jobs
  - 待生成90个新jobs
```

### 日志文件
```
E:/Projects/halogenator/transform_pipeline_final.log (最近一次运行)
E:/Projects/halogenator/transform_pipeline_optimized.log (优化后新日志)
```

---

## 🔧 Part F: 验证和监控

### 内存监控命令

**实时监控Python进程内存:**
```python
# 创建监控脚本: monitor_memory.py
import psutil
import time

def monitor_python_memory():
    while True:
        total_mem = 0
        for proc in psutil.process_iter(['pid', 'name', 'memory_info']):
            try:
                if 'python' in proc.info['name'].lower():
                    mem_mb = proc.info['memory_info'].rss / 1024 / 1024
                    total_mem += mem_mb
                    print(f"PID {proc.info['pid']}: {mem_mb:.1f} MB")
            except:
                pass

        total_system = psutil.virtual_memory().total / 1024 / 1024
        percent = (total_mem / total_system) * 100
        print(f"Total Python memory: {total_mem:.1f} MB ({percent:.1f}%)")
        print(f"Target: <50% ({total_system * 0.5:.1f} MB)")
        print("-" * 50)
        time.sleep(60)  # 每分钟检查一次

if __name__ == '__main__':
    monitor_python_memory()
```

**使用:**
```bash
python monitor_memory.py &
```

### Pipeline进度监控

**检查完成jobs:**
```bash
ls E:/Projects/halogenator/data/output/transforms/ | wc -l
```

**检查最新日志:**
```bash
tail -50 E:/Projects/halogenator/transform_pipeline_optimized.log
```

**统计总products:**
```python
import pyarrow.parquet as pq
from pathlib import Path

total = sum(
    pq.read_table(f).num_rows
    for f in Path('E:/Projects/halogenator/data/output/transforms').glob('*/products.parquet')
)
print(f'Total products: {total:,}')
```

### 成功标准

**Task 1成功标准:**
- [ ] Python进程内存稳定在40-50%
- [ ] 无内存相关错误/僵死
- [ ] Flush日志显示"memory_pressure"触发
- [ ] 至少完成1个polyphenol-2X phenolic OH job

**Pipeline完成标准:**
- [ ] 118个jobs全部完成
- [ ] 总products > 750M
- [ ] 0个失败jobs
- [ ] 最终报告生成: `transform_pipeline_report.json`

---

## 🚨 Part G: 故障排除

### 问题1: 内存仍然超过50%

**诊断:**
```python
# 检查每个worker的buffer大小
# 在StreamingParquetWriter添加:
def get_buffer_stats(self):
    return {
        'buffer_size': len(self.buffer),
        'memory_mb': self._get_process_memory_mb(),
        'memory_percent': self._get_memory_usage_percent()
    }
```

**解决方案:**
1. 降低`--target-memory`到40或35
2. 降低`--flush-interval`到2000或1000
3. 临时降低workers到8验证

### 问题2: Flush过于频繁影响性能

**症状:** 日志中频繁出现flush消息，throughput降低

**解决方案:**
1. 调整target_memory阈值范围（45-55改为40-60）
2. 添加最小flush间隔（防止连续flush）:
   ```python
   import time

   def _flush(self):
       # Check minimum interval since last flush
       current_time = time.time()
       if hasattr(self, '_last_flush_time'):
           time_since_last = current_time - self._last_flush_time
           if time_since_last < 5.0:  # 至少5秒间隔
               return

       self._last_flush_time = current_time
       # 继续现有flush逻辑...
   ```

### 问题3: Dedup数据库仍然膨胀

**症状:** dedup.db超过20GB

**激进解决方案:**
启用分批dedup（牺牲少量唯一性）:
```python
# 每N个products后重置dedup数据库
DEDUP_RESET_INTERVAL = 1000000  # 100万

if total_processed % DEDUP_RESET_INTERVAL == 0:
    deduper.close()
    os.remove(dedup_db_path)
    deduper = DedupManager(dedup_db_path)  # 重新初始化
    logger.warning(f"Dedup DB reset at {total_processed:,} products")
```

---

## 📝 Part H: Quick Start for Next Session

### Step-by-Step执行清单

**Session开始检查清单:**
```bash
# 1. 确认环境
cd E:/Projects/halogenator
conda activate halo-p0  # 或你的环境
python --version  # 应该是3.8+

# 2. 验证当前状态
ls data/output/transforms/ | wc -l  # 应该是28

# 3. 检查没有遗留进程
ps aux | grep python | grep -v grep  # 应该为空
```

**实施内存优化 (按顺序):**

1. **修改StreamingParquetWriter (CRITICAL):**
   - 打开 `scripts/08_transform_library_v2.py`
   - 找到 `class StreamingParquetWriter` (line 301)
   - 按照 **Part C, Task 1** 的详细方案修改:
     - 添加 `_get_process_memory_mb()` 等3个新方法
     - 修改 `write_batch()` 实现自适应flush
     - 修改 `_flush()` 添加GC
     - 修改 `__init__()` 添加参数
   - 保存文件

2. **修改CLI参数:**
   - 打开 `scripts/run_transform_pipeline.py`
   - 搜索 `--flush-interval`
   - 添加 `--target-memory` 参数
   - 修改writer初始化传参
   - 保存文件

3. **验证修改:**
   ```bash
   # 语法检查
   python -m py_compile scripts/08_transform_library_v2.py
   python -m py_compile scripts/run_transform_pipeline.py

   # 帮助检查
   python scripts/run_transform_pipeline.py apply --help
   # 应该看到 --target-memory 参数
   ```

4. **清理失败jobs:**
   ```bash
   rm -rf data/output/transforms/polyphenol-2X_FG_PHENOL_OH__OH__TO__OMe
   rm -rf data/output/transforms/polyphenol-2X_FG_PHENOL_OH__OH__TO__NH2
   ```

5. **小规模测试 (aa_peptide-1X):**
   ```bash
   # 单个job测试验证内存控制
   python scripts/08_transform_library_v2.py apply \
     --input data/output/nplike_v2/aa_peptide-1X/products.parquet \
     --outdir data/output/transforms/test_aa_peptide-1X_FG_AMINE_AR__NH2__TO__F \
     --xf-config configs/transforms.yaml \
     --xf-name FG_AMINE_AR__NH2__TO__F \
     --workers 16 \
     --flush-interval 5000 \
     --batch-size 50000

   # 监控内存（另一个终端）
   watch -n 5 "ps aux | grep python"

   # 成功标准: 内存<50%, job完成无错误
   ```

6. **启动完整pipeline:**
   ```bash
   nohup python scripts/run_transform_pipeline.py \
     --workers 16 \
     --flush-interval 5000 \
     --batch-size 50000 \
     --target-memory 45.0 \
     > transform_pipeline_optimized.log 2>&1 &

   echo "Pipeline PID: $!"
   ```

7. **持续监控:**
   ```bash
   # 终端1: 日志
   tail -f transform_pipeline_optimized.log

   # 终端2: 内存（如果实现了monitor_memory.py）
   python monitor_memory.py

   # 或手动检查
   watch -n 60 "ps aux | grep python | grep -v grep"
   ```

---

## 🎯 Part I: Success Metrics

### 短期目标 (第一个session)
- [ ] 内存优化实施完成
- [ ] 至少完成10个新jobs（推荐从aa_peptide开始）
- [ ] 内存占用稳定在50%以下
- [ ] 至少1个polyphenol-2X phenolic OH job成功完成

### 中期目标 (2-3个sessions)
- [ ] 完成所有小/中型jobs (aa_peptide, alkaloid, polyphenol-1X剩余)
- [ ] 累计完成60+个jobs
- [ ] 总products > 100M

### 最终目标
- [ ] 118个jobs全部完成
- [ ] 总products > 750M
- [ ] 生成最终报告
- [ ] 所有parquet文件完整可读

---

## 📚 Part J: 参考资源

### 相关文档
- `SESSION_HANDOFF_TRANSFORM_PIPELINE_2025-12-13.md` - 上一个交接报告
- `IMPLEMENTATION_COMPLETE_NP_HALOGEN_PIPELINE.md` - Halogenation pipeline完成报告
- `configs/transforms.yaml` - Transform rules完整定义

### 性能基准
- **Lipid-1X平均:** ~25秒/job, 6 workers
- **Polyphenol-1X平均:**
  - Phenolic OH: ~2,700秒/job (~45分钟)
  - Aromatic Amine: ~80秒/job
  - Carboxyl: ~120秒/job
- **预期polyphenol-2X (优化后):** ~6-8小时/phenolic OH job

### 内存占用估算
- **目标配置:**
  - 16 workers
  - 5,000 flush_interval
  - 45% target memory

- **预期峰值内存:**
  - 16 workers × 5,000 products × 8KB = 640MB (buffer)
  - Dedup DB: ~5-10GB (with checkpointing)
  - Overhead: ~2GB
  - **Total: ~8-13GB (25-40% on 32GB system)** ✓

---

## ⚠️ Part K: Critical Warnings

1. **不要同时运行多个pipeline实例** - 会竞争内存
2. **不要在pipeline运行时关机/重启** - 导致文件损坏
3. **监控磁盘空间** - transform产品需要大量空间（预留500GB+）
4. **定期检查日志** - 至少每4-6小时检查一次
5. **如果内存仍超过60%** - 立即暂停并调整参数

---

## 🏁 Part L: Session Conclusion Checklist

**下一个session结束时，验证:**
- [ ] 修改的代码已保存并通过语法检查
- [ ] 至少1个测试job成功完成并验证内存<50%
- [ ] Pipeline在后台稳定运行
- [ ] 创建了新的session handoff报告（如果未完成）
- [ ] 记录了遇到的任何新问题和解决方案

---

**End of Handoff Report**

**Created:** 2025-12-20
**For Session:** Next Claude session
**Priority:** 🔴 HIGH - Memory optimization critical for pipeline completion
**Estimated Work:** 2-4 hours implementation + 7-10 days pipeline runtime

**Direct Questions to User if Unclear:**
- Exact system RAM amount (assumed 32GB)
- Acceptable max pipeline runtime (assumed 10 days OK)
- Priority: Speed vs Memory safety (assumed: memory safety higher)
