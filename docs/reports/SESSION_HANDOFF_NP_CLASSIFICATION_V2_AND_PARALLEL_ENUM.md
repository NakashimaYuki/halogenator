# Session Handoff: NP Classification v2.0 + Parallel Enumeration

**Date:** 2025-12-08
**Project:** Halogenated Natural Product Library Generation
**Status:** Classification v2.0 Complete, Ready for k=1 and k=2 Enumeration

---

## 本会话已完成的工作

### 1. 并行枚举系统实现 ✅

**问题诊断：**
- 原有枚举系统CPU利用率低（单核运行）
- k=2枚举速度慢（terpenoid预计5-10天）

**解决方案：**
创建了完整的多进程并行枚举系统：

**文件创建/修改：**
- `src/halogenator/parallel_enum.py` - 新建，完整的并行枚举引擎
  - `ParallelEnumerator` 类：管理worker pool和流式写入
  - `_enum_worker()` 函数：worker进程处理单个父分子
  - 周期性刷盘机制：每10,000个产物写入一次（可配置 `--flush-interval`）
  - 内存安全：防止内存爆炸

- `src/halogenator/cli.py` - 修改，集成并行模式
  - 添加 `--workers` 参数（默认1=顺序模式）
  - 添加 `--flush-interval` 参数（默认10000）
  - 检测workers>1时自动切换到并行模式

- `scripts/04_enum_halogen_all_classes.py` - 修改，支持并行参数
  - 添加 `--workers` 和 `--flush-interval` 参数
  - 传递给底层enum-parquet命令

**验证结果：**
- 测试用例：lipid k=2（18个父分子 → 258个产物）
- 顺序模式：~3秒
- 并行模式（4 workers）：~4秒
- **输出完全一致**（258个相同的InChIKey）
- 性能提升：小数据集提升不明显，大数据集预计6-12x加速

**已知问题及修复：**
- Schema不匹配问题：已修复（使用pandas append模式）
- 睡眠模式导致worker死亡：用户已禁用睡眠

---

### 2. NP分类系统v2.0运行完成 ✅

**背景：**
- 旧系统将glycoside设为最高优先级主类
- 导致22,873个分子（33.5%）错误分类为glycoside
- 例如槲皮素-3-葡萄糖苷应该是polyphenol，却被分为glycoside

**新系统架构：**
- **主类（primary class）**：骨架类型（aa_peptide, alkaloid, terpenoid, polyphenol, lipid, polysaccharide, other）
- **标签（tags）**：糖修饰信息（glycoside, saponin, cglycoside_like）

**运行命令：**
```bash
cd E:/Projects/halogenator
python scripts/02_partition_nplike_by_class.py \
  -i data/output/nplike/nplike_merged_clean.parquet \
  -o data/output/nplike_v2 \
  --split \
  --exclude-parquet data/output/nplike/Flavone/base_clean.parquet
```

**分类结果对比：**

| 主类 | 旧系统 | 新系统v2.0 | 变化 |
|------|--------|-----------|------|
| glycoside | 22,873 | **0** | -22,873 ✅ |
| terpenoid | 28,606 | 35,131 | +6,525 |
| polyphenol | 4,047 | 13,168 | **+9,121** |
| alkaloid | 4,287 | 7,871 | +3,584 |
| lipid | 18 | 6,247 | **+6,229** |
| aa_peptide | 3,890 | 1,119 | -2,771 |
| polysaccharide | 265 | 566 | +301 |
| other | 4,262 | 4,146 | -116 |
| **总计** | **68,248** | **68,248** | - |

**关键改进：**
- ✅ Glycoside不再作为主类（0个，完美！）
- ✅ 16,387个分子有glycoside标签（正确保留糖信息）
- ✅ Polyphenol增加225%（黄酮糖苷正确分类）
- ✅ Lipid增加34,000%（之前被严重低估）

**输出文件位置：**
- 合并文件：`data/output/nplike_v2/nplike_with_classes.parquet`
- Per-class文件：`data/output/nplike_v2/{class}/base.parquet`
  - aa_peptide/base.parquet (1,119 molecules)
  - alkaloid/base.parquet (7,871 molecules)
  - lipid/base.parquet (6,247 molecules)
  - polyphenol/base.parquet (13,168 molecules)
  - terpenoid/base.parquet (35,131 molecules)
  - polysaccharide/base.parquet (566 molecules)
  - other/base.parquet (4,146 molecules)

---

## 待完成任务（详细操作指南）

### 任务1：重命名base.parquet为base_clean.parquet

**目标：** 枚举脚本期望文件名为`base_clean.parquet`，但分类脚本生成的是`base.parquet`

**操作：**
```bash
cd E:/Projects/halogenator/data/output/nplike_v2

# 重命名所有per-class文件
for dir in aa_peptide alkaloid lipid polyphenol terpenoid polysaccharide other; do
  mv $dir/base.parquet $dir/base_clean.parquet
done

# 验证
ls -lh */base_clean.parquet
```

**预期结果：**
```
aa_peptide/base_clean.parquet    (418K)
alkaloid/base_clean.parquet      (1.7M)
lipid/base_clean.parquet         (1.3M)
polyphenol/base_clean.parquet    (2.7M)
terpenoid/base_clean.parquet     (6.4M)
polysaccharide/base_clean.parquet (329K)
other/base_clean.parquet         (996K)
```

---

### 任务2：验证分类结果（快速QC）

**目标：** 确认每个class文件内容正确

**操作：**
```bash
cd E:/Projects/halogenator

python -c "
import pandas as pd
from pathlib import Path

base_dir = Path('data/output/nplike_v2')
classes = ['aa_peptide', 'alkaloid', 'lipid', 'polyphenol', 'terpenoid', 'polysaccharide', 'other']

print('Per-Class File Validation:')
print('=' * 70)

for cls in classes:
    file_path = base_dir / cls / 'base_clean.parquet'
    if file_path.exists():
        df = pd.read_parquet(file_path)
        print(f'{cls:15s}: {len(df):6,} molecules, columns: {list(df.columns)[:5]}')
    else:
        print(f'{cls:15s}: FILE NOT FOUND!')

print('=' * 70)
"
```

**预期输出：**
```
aa_peptide     :  1,119 molecules
alkaloid       :  7,871 molecules
lipid          :  6,247 molecules
polyphenol     : 13,168 molecules
terpenoid      : 35,131 molecules
polysaccharide :    566 molecules
other          :  4,146 molecules
```

**成功标准：**
- ✅ 所有7个class文件都存在
- ✅ 分子数量匹配分类结果
- ✅ 总和 = 68,248 molecules

---

### 任务3：运行k=1枚举（所有7个类，16 workers）

**目标：** 对v2.0分类的7个类别分别进行k=1卤代枚举，使用并行模式加速

**重要前提：**
- ✅ 已禁用电脑睡眠模式
- ✅ 有足够磁盘空间（至少50GB可用）
- ✅ parallel_enum.py已创建并集成到CLI

**操作命令：**
```bash
cd E:/Projects/halogenator

# 方式1：逐个运行（推荐，便于监控）
python scripts/04_enum_halogen_all_classes.py \
  --classes aa_peptide \
  --k-values 1 \
  --workers 16 \
  --flush-interval 10000 \
  2>&1 | tee logs/v2_aa_peptide_k1_w16.log

python scripts/04_enum_halogen_all_classes.py \
  --classes alkaloid \
  --k-values 1 \
  --workers 16 \
  --flush-interval 10000 \
  2>&1 | tee logs/v2_alkaloid_k1_w16.log

python scripts/04_enum_halogen_all_classes.py \
  --classes lipid \
  --k-values 1 \
  --workers 16 \
  --flush-interval 10000 \
  2>&1 | tee logs/v2_lipid_k1_w16.log

python scripts/04_enum_halogen_all_classes.py \
  --classes polyphenol \
  --k-values 1 \
  --workers 16 \
  --flush-interval 10000 \
  2>&1 | tee logs/v2_polyphenol_k1_w16.log

python scripts/04_enum_halogen_all_classes.py \
  --classes terpenoid \
  --k-values 1 \
  --workers 16 \
  --flush-interval 10000 \
  2>&1 | tee logs/v2_terpenoid_k1_w16.log

python scripts/04_enum_halogen_all_classes.py \
  --classes polysaccharide \
  --k-values 1 \
  --workers 16 \
  --flush-interval 10000 \
  2>&1 | tee logs/v2_polysaccharide_k1_w16.log

python scripts/04_enum_halogen_all_classes.py \
  --classes other \
  --k-values 1 \
  --workers 16 \
  --flush-interval 10000 \
  2>&1 | tee logs/v2_other_k1_w16.log


# 方式2：一次性运行所有（如果机器稳定）
python scripts/04_enum_halogen_all_classes.py \
  --classes aa_peptide alkaloid lipid polyphenol terpenoid polysaccharide other \
  --k-values 1 \
  --workers 16 \
  --flush-interval 10000 \
  2>&1 | tee logs/v2_all_k1_w16.log
```

**注意事项：**
1. **输入文件路径问题：** 脚本默认查找`data/output/nplike/{class}/base_clean.parquet`
   - 但新分类结果在`data/output/nplike_v2/{class}/base_clean.parquet`
   - **需要修改脚本或创建符号链接！**

**修复方案A - 修改脚本（推荐）：**
```python
# 在 scripts/04_enum_halogen_all_classes.py 第70行左右
# 修改前：
data_dir = Path("E:/Projects/halogenator/data/output/nplike")

# 修改后：
data_dir = Path("E:/Projects/halogenator/data/output/nplike_v2")
```

**修复方案B - 创建符号链接：**
```bash
cd E:/Projects/halogenator/data/output/nplike
rm -rf aa_peptide alkaloid lipid polyphenol terpenoid polysaccharide other

# Windows上使用mklink /D创建目录链接（需要管理员权限）
mklink /D aa_peptide ..\nplike_v2\aa_peptide
mklink /D alkaloid ..\nplike_v2\alkaloid
mklink /D lipid ..\nplike_v2\lipid
mklink /D polyphenol ..\nplike_v2\polyphenol
mklink /D terpenoid ..\nplike_v2\terpenoid
mklink /D polysaccharide ..\nplike_v2\polysaccharide
mklink /D other ..\nplike_v2\other
```

**预期结果（基于旧分类的k=1数据推算）：**

| Class | Parents | Estimated k=1 Products | Estimated Runtime (16 workers) |
|-------|---------|------------------------|--------------------------------|
| aa_peptide | 1,119 | ~28,000 | ~2min |
| alkaloid | 7,871 | ~290,000 | ~15min |
| lipid | 6,247 | ~15,000 | ~5min |
| polyphenol | 13,168 | ~430,000 | ~25min |
| terpenoid | 35,131 | ~1,340,000 | ~90min |
| polysaccharide | 566 | ~50,000 | ~3min |
| other | 4,146 | ~100,000 | ~10min |
| **总计** | **68,248** | **~2,253,000** | **~150min (2.5h)** |

**成功标准：**
- ✅ 所有7个类完成枚举无报错
- ✅ 输出文件：`data/output/nplike_v2/{class}-1X/products.parquet`
- ✅ ALPHA_CARBONYL规则有产物（>1%）
- ✅ RING_SP3规则为主导（~60-70%）

---

### 任务4：验证k=1结果

**目标：** 确认k=1枚举正确，ALPHA_CARBONYL bug已修复

**操作：**
```bash
cd E:/Projects/halogenator

python -c "
import pandas as pd
from pathlib import Path

base_dir = Path('data/output/nplike_v2')
classes = ['aa_peptide', 'alkaloid', 'lipid', 'polyphenol', 'terpenoid', 'polysaccharide', 'other']

print('k=1 Enumeration Validation:')
print('=' * 80)

total_products = 0
for cls in classes:
    k1_file = base_dir / f'{cls}-1X' / 'products.parquet'

    if not k1_file.exists():
        print(f'{cls:15s}: NOT FOUND')
        continue

    df = pd.read_parquet(k1_file)
    parent_file = base_dir / cls / 'base_clean.parquet'
    df_parents = pd.read_parquet(parent_file)

    # Rule distribution
    rule_dist = df['rule'].value_counts()
    alpha_count = rule_dist.get('ALPHA_CARBONYL__CH2__TO__X', 0)
    ring_sp3_count = rule_dist.get('RING_SP3__CH__TO__X', 0)

    print(f'{cls:15s}: {len(df_parents):6,} parents → {len(df):8,} products ({len(df)/len(df_parents):5.1f} prod/parent)')
    print(f'                ALPHA_CARBONYL: {alpha_count:6,} ({alpha_count/len(df)*100:4.1f}%)')
    print(f'                RING_SP3:       {ring_sp3_count:6,} ({ring_sp3_count/len(df)*100:4.1f}%)')

    total_products += len(df)

print('=' * 80)
print(f'TOTAL k=1 PRODUCTS: {total_products:,}')
"
```

**成功标准：**
- ✅ ALPHA_CARBONYL在所有类中都有产物（>1%）
- ✅ RING_SP3为主导规则（50-70%）
- ✅ 总产物数 > 2M
- ✅ 无ERROR SMILES或InChIKeys

---

### 任务5：运行k=2枚举（所有7个类，16 workers）

**目标：** 对v2.0分类的7个类别分别进行k=2卤代枚举

**重要：k=2枚举时间较长，建议分批运行**

**预计时间：**
- **快速批次（lipid, aa_peptide, polysaccharide, other）：** 1-3小时
- **中速批次（alkaloid, polyphenol）：** 8-15小时
- **慢速批次（terpenoid）：** 1-3天

**操作命令：**

**第一批：快速类（建议先运行，验证系统稳定）**
```bash
cd E:/Projects/halogenator

# Lipid (6,247 parents, ~150K products估计, ~1h)
python scripts/04_enum_halogen_all_classes.py \
  --classes lipid \
  --k-values 2 \
  --workers 16 \
  --flush-interval 10000 \
  2>&1 | tee logs/v2_lipid_k2_w16.log &

# AA_peptide (1,119 parents, ~560K products估计, ~2h)
python scripts/04_enum_halogen_all_classes.py \
  --classes aa_peptide \
  --k-values 2 \
  --workers 16 \
  --flush-interval 10000 \
  2>&1 | tee logs/v2_aa_peptide_k2_w16.log &

# Polysaccharide (566 parents, ~1M products估计, ~3h)
python scripts/04_enum_halogen_all_classes.py \
  --classes polysaccharide \
  --k-values 2 \
  --workers 16 \
  --flush-interval 10000 \
  2>&1 | tee logs/v2_polysaccharide_k2_w16.log &

# Other (4,146 parents, ~2M products估计, ~5h)
python scripts/04_enum_halogen_all_classes.py \
  --classes other \
  --k-values 2 \
  --workers 16 \
  --flush-interval 10000 \
  2>&1 | tee logs/v2_other_k2_w16.log &
```

**第二批：中速类（快速批次成功后运行）**
```bash
# Polyphenol (13,168 parents, ~9M products估计, ~12h)
python scripts/04_enum_halogen_all_classes.py \
  --classes polyphenol \
  --k-values 2 \
  --workers 16 \
  --flush-interval 10000 \
  2>&1 | tee logs/v2_polyphenol_k2_w16.log &

# Alkaloid (7,871 parents, ~6M products估计, ~10h)
python scripts/04_enum_halogen_all_classes.py \
  --classes alkaloid \
  --k-values 2 \
  --workers 16 \
  --flush-interval 10000 \
  2>&1 | tee logs/v2_alkaloid_k2_w16.log &
```

**第三批：慢速类（最后运行，耗时最长）**
```bash
# Terpenoid (35,131 parents, ~28M products估计, ~48-72h即2-3天)
python scripts/04_enum_halogen_all_classes.py \
  --classes terpenoid \
  --k-values 2 \
  --workers 16 \
  --flush-interval 5000 \
  2>&1 | tee logs/v2_terpenoid_k2_w16.log &

# 注意：terpenoid使用更小的flush-interval (5000) 以防内存问题
```

**监控进度：**
```bash
# 查看实时日志
tail -f logs/v2_terpenoid_k2_w16.log

# 检查产物文件大小（判断进度）
watch -n 60 'ls -lh data/output/nplike_v2/*-2X/products.parquet'

# 检查CPU占用
top -p $(pgrep -f "enum_halogen_all_classes")
```

**预期结果（基于21.9x扩增因子）：**

| Class | k=1 Products | k=2 Products (估计) | Runtime (16w) |
|-------|--------------|---------------------|---------------|
| lipid | ~15K | ~330K | ~1h |
| aa_peptide | ~28K | ~613K | ~2h |
| polysaccharide | ~50K | ~1.1M | ~3h |
| other | ~100K | ~2.2M | ~5h |
| alkaloid | ~290K | ~6.4M | ~10h |
| polyphenol | ~430K | ~9.4M | ~12h |
| terpenoid | ~1,340K | ~29.4M | **2-3天** |
| **总计** | **~2.25M** | **~49.4M** | **3-4天** |

**成功标准：**
- ✅ 所有类完成k=2枚举
- ✅ 输出文件：`data/output/nplike_v2/{class}-2X/products.parquet`
- ✅ 无schema错误（已修复）
- ✅ extreme_site hard skip < 5%

---

### 任务6：生成最终统计报告

**目标：** 汇总k=1和k=2的全部结果，生成完整报告

**操作：**
```bash
cd E:/Projects/halogenator

python -c "
import pandas as pd
from pathlib import Path
import json

base_dir = Path('data/output/nplike_v2')
classes = ['aa_peptide', 'alkaloid', 'lipid', 'polyphenol', 'terpenoid', 'polysaccharide', 'other']

print('=' * 100)
print('NP CLASSIFICATION V2.0 + HALOGENATION ENUMERATION - FINAL REPORT')
print('=' * 100)
print()

# Classification summary
print('CLASSIFICATION V2.0 SUMMARY:')
print('-' * 100)
df_all = pd.read_parquet(base_dir / 'nplike_with_classes.parquet')
print(df_all['np_primary_class'].value_counts())
print()

# k=1 summary
print('K=1 ENUMERATION SUMMARY:')
print('-' * 100)
print(f'{'Class':15s} {'Parents':>8s} {'k=1 Products':>12s} {'Prod/Parent':>12s} {'ALPHA_CARBONYL%':>15s}')
print('-' * 100)

k1_total = 0
for cls in classes:
    parent_file = base_dir / cls / 'base_clean.parquet'
    k1_file = base_dir / f'{cls}-1X' / 'products.parquet'

    if not k1_file.exists():
        continue

    df_parents = pd.read_parquet(parent_file)
    df_k1 = pd.read_parquet(k1_file)

    alpha_pct = (df_k1['rule'] == 'ALPHA_CARBONYL__CH2__TO__X').sum() / len(df_k1) * 100

    print(f'{cls:15s} {len(df_parents):8,} {len(df_k1):12,} {len(df_k1)/len(df_parents):12.1f} {alpha_pct:14.1f}%')
    k1_total += len(df_k1)

print('-' * 100)
print(f'{'TOTAL':15s} {len(df_all):8,} {k1_total:12,}')
print()

# k=2 summary
print('K=2 ENUMERATION SUMMARY:')
print('-' * 100)
print(f'{'Class':15s} {'k=1 Products':>12s} {'k=2 Products':>12s} {'Expansion':>10s} {'File Size':>10s}')
print('-' * 100)

k2_total = 0
for cls in classes:
    k1_file = base_dir / f'{cls}-1X' / 'products.parquet'
    k2_file = base_dir / f'{cls}-2X' / 'products.parquet'

    if not k2_file.exists():
        continue

    df_k1 = pd.read_parquet(k1_file)
    df_k2 = pd.read_parquet(k2_file)

    expansion = len(df_k2) / len(df_k1)
    file_size_mb = k2_file.stat().st_size / 1024 / 1024

    print(f'{cls:15s} {len(df_k1):12,} {len(df_k2):12,} {expansion:10.1f}x {file_size_mb:9.1f}MB')
    k2_total += len(df_k2)

print('-' * 100)
print(f'{'TOTAL':15s} {k1_total:12,} {k2_total:12,}')
print()
print('=' * 100)
print(f'GRAND TOTAL (k=1 + k=2 unique): Requires deduplication')
print('=' * 100)
"
```

**报告保存：**
将上述输出重定向到文件：
```bash
python [上述脚本] > FINAL_ENUMERATION_REPORT_V2.txt
```

---

## 技术要点总结

### 并行枚举关键参数

**--workers N (推荐值)：**
- 小数据集（<1000 parents）：4-8 workers
- 中数据集（1000-10000 parents）：8-12 workers
- 大数据集（>10000 parents）：12-16 workers
- 不建议超过CPU核心数

**--flush-interval N (推荐值)：**
- 默认：10000（平衡性能和内存）
- 内存紧张：5000（更频繁写入，更安全）
- 内存充足：20000（减少I/O，更快）

### 内存估算

**k=2枚举内存需求：**
- 每个产物：~1KB
- Flush buffer：10,000 × 1KB = ~10MB
- Worker并发：16 workers × 800 products/parent × 1KB = ~13MB
- **总峰值：~50-100MB**（非常安全）

### 故障排查

**如果枚举卡死：**
1. 检查CPU占用：`top -p $(pgrep -f python)`
2. 检查日志最后输出：`tail -f logs/v2_*.log`
3. 检查磁盘空间：`df -h`
4. 如果真的卡死，杀进程重启：`pkill -f enum_halogen_all_classes`

**如果schema错误：**
- 已修复！使用pandas append模式
- 如果仍出现，检查parallel_enum.py是否最新版本

**如果内存爆炸：**
- 降低flush-interval到5000或更低
- 减少workers数量

---

## 文件路径速查

**分类输出：**
- 主文件：`data/output/nplike_v2/nplike_with_classes.parquet`
- Per-class：`data/output/nplike_v2/{class}/base_clean.parquet`

**枚举输出：**
- k=1：`data/output/nplike_v2/{class}-1X/products.parquet`
- k=2：`data/output/nplike_v2/{class}-2X/products.parquet`

**日志：**
- 分类：`logs/classification_v2_clean.log`
- k=1：`logs/v2_{class}_k1_w16.log`
- k=2：`logs/v2_{class}_k2_w16.log`

**代码：**
- 并行引擎：`src/halogenator/parallel_enum.py`
- CLI入口：`src/halogenator/cli.py`
- 分类脚本：`scripts/02_partition_nplike_by_class.py`
- 枚举脚本：`scripts/04_enum_halogen_all_classes.py`

---

## 下一步建议

1. **首先**：完成任务1（重命名base.parquet）和任务2（验证）
2. **然后**：运行k=1枚举（全部7类，预计2.5小时）
3. **验证k=1**：确认ALPHA_CARBONYL正常工作
4. **分批k=2**：先快速类，再中速类，最后terpenoid
5. **最终报告**：所有枚举完成后生成统计

**预计总时间：**
- k=1：2.5小时
- k=2快速批：12小时
- k=2中速批：24小时
- k=2慢速批：2-3天
- **总计：约4天完成全部枚举**

---

## 重要提醒

⚠️ **运行terpenoid k=2时：**
- 确保电脑不会睡眠
- 确保磁盘有50GB+空间
- 使用flush-interval=5000（更安全）
- 可能需要2-3天连续运行
- 定期检查日志确认正常运行

✅ **成功标志：**
- 分类v2.0：glycoside主类 = 0
- k=1：ALPHA_CARBONYL有产物（>1%）
- k=2：所有类完成，总产物约50M

---

**会话结束时间：** 2025-12-08 21:05
**下一个Claude会话应从任务1开始执行**
