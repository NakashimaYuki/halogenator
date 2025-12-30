# NP卤代流水线关键问题与解决方案

**日期：** 2025-12-01
**版本：** v2.0
**状态：** 待执行

---

## 执行摘要

当前NP卤代流水线k=1枚举已完成，但发现三个关键问题需要解决：

1. **规则覆盖不足**：仅使用R1和R5，大量合理规则未启用
2. **黄酮重复问题**：黄酮类已单独处理，但未从NP合并库中排除
3. **多重分类冲突**：部分分子可能同时满足多个NP类别标准

本文档提供详细的问题分析和解决方案。

---

## 问题1：规则覆盖不足

### 1.1 现状分析

**当前k=1枚举结果显示：**
- 仅2个规则生效：R1 (芳环卤代) 和 R5 (羧基卤代)
- R1占90-98%，R5占2-10%
- 其他规则（R2, R3, R4, RING_SP3, ALPHA_CARBONYL, PRIMARY_OH）**完全未产生产物**

**根本原因：**

查看 `configs/halogen_rules_by_class.yaml`，发现：

1. **include_rules列表过于保守**
   ```yaml
   # 示例：polyphenol k1配置
   k1:
     include_rules:
       - RING_SP2__CH__TO__X     # R1 - 唯一芳环规则
       - RING_SP3__CH__TO__X     # 未生效？
       - COOH__TO__CX            # R5
   ```

2. **缺失关键规则**
   - **R3 (OH卤代)**: 完全未包含在任何类的include_rules中
   - **R4 (醛基卤代)**: 完全未包含
   - **R2 (黄酮C环α位)**: 仅在黄酮特定流水线使用，未包含在通用NP配置中

3. **RING_SP3规则未生效**
   - 虽然在include_rules中，但k=1结果显示无产物
   - 可能原因：
     - 该规则定义的SMARTS模式过于严格
     - 或者在NP分子中匹配的位点极少

### 1.2 可用规则清单

根据 `configs/transforms.yaml` 和 `scripts/list_halogen_rules.py` 输出：

| 规则ID | 语义名称 | 化学描述 | 当前状态 |
|--------|----------|----------|----------|
| **R1** | RING_SP2__CH__TO__X | 芳环/杂芳环C-H卤代 | ✅ 已用 |
| **R2a/R2b** | FLAVONE_C_RING | 黄酮C环特定位点 | ❌ 未用（黄酮特定） |
| **R3** | OH__TO__X | 羟基卤代 (OH → X) | ❌ **未包含** |
| **R4** | CHO__TO__CX | 醛基卤代 (CHO → CX) | ❌ **未包含** |
| **R5** | COOH__TO__CX | 羧基卤代 (COOH → COX) | ✅ 已用 |
| **R6_methyl** | METHYL_HALOGEN | 甲基多卤代 (CH3 → CX3) | ⚠️ 配置中disabled |
| RING_SP3__CH__TO__X | - | 环sp3 C-H卤代 | ⚠️ 包含但无产物 |
| ALPHA_CARBONYL__CH2__TO__X | - | 羰基α位CH2卤代 | ⚠️ 仅k2包含 |
| PRIMARY_OH__CH2OH__TO__X | - | 伯醇CH2OH卤代 | ❌ 多数类disabled |

### 1.3 解决方案

#### 方案A：激进扩展（推荐）

**目标：** 最大化化学空间，启用所有合理规则

**配置修改：**

```yaml
# 以terpenoid为例（其他类类似）
terpenoid:
  k1:
    include_rules:
      - RING_SP2__CH__TO__X           # R1 芳环
      - RING_SP3__CH__TO__X           # 脂环C-H
      - COOH__TO__CX                  # R5 羧基
      - R3                            # NEW: OH卤代
      - R4                            # NEW: 醛基卤代
      - ALPHA_CARBONYL__CH2__TO__X    # α-羰基
      - PRIMARY_OH__CH2OH__TO__X      # 伯醇

  k2:
    include_rules:
      - RING_SP2__CH__TO__X
      - RING_SP3__CH__TO__X
      - COOH__TO__CX
      - R3                            # NEW
      - R4                            # NEW
      - ALPHA_CARBONYL__CH2__TO__X
      - PRIMARY_OH__CH2OH__TO__X

  per_rule_overrides:
    # 移除所有max_sites_per_parent限制
    # 仅保留极端情况的safety nets
    R3:
      max_sites_per_parent: 8         # OH通常较多，设较高限制
    R4:
      max_sites_per_parent: 3         # 醛基较少
    PRIMARY_OH__CH2OH__TO__X:
      enabled: true                   # 重新启用
      max_sites_per_parent: 5
```

**各类推荐规则组合：**

| NP类 | 核心规则 | 可选规则 | 禁用规则 |
|------|----------|----------|----------|
| **Terpenoid** | R1, R3, R4, R5, RING_SP3, ALPHA_CARBONYL | PRIMARY_OH | R6_methyl |
| **Alkaloid** | R1, R3, R4, R5, RING_SP3 | ALPHA_CARBONYL | PRIMARY_OH, R6_methyl |
| **Polyphenol** | R1, R3, R5, RING_SP3 | R4, PRIMARY_OH | R6_methyl |
| **Glycoside** | R1, R3, R5 (+ sugar_mask) | RING_SP3 | PRIMARY_OH, R4 |
| **Lipid** | R5, ALPHA_CARBONYL | PRIMARY_OH | R1, R3, R4 |
| **AA_peptide** | R1, R5, ALPHA_CARBONYL | R3 | R4, PRIMARY_OH |

#### 方案B：保守扩展

仅添加R3和R4，其他规则保持当前配置。

**优点：** 风险低，增量可控
**缺点：** 化学空间扩展有限

### 1.4 实施步骤

1. **备份当前配置**
   ```bash
   cp configs/halogen_rules_by_class.yaml configs/halogen_rules_by_class.yaml.backup_k1only
   ```

2. **更新配置文件**
   - 为每个类的k1和k2添加R3, R4
   - 移除所有max_sites_per_parent（全局和规则级）
   - 仅保留极端情况的safety nets (8-12)

3. **重新运行POC测试**（强烈推荐）
   ```bash
   # 测试扩展规则集不会导致组合爆炸
   python scripts/03_enum_halogen_poc.py --class terpenoid --k 1 --max-parents 1000
   python scripts/03_enum_halogen_poc.py --class polyphenol --k 1 --max-parents 1000
   ```

4. **如果POC通过，重新运行k=1全量**
   ```bash
   # 删除旧的1X输出
   rm -rf data/output/nplike/*-1X

   # 重新枚举
   python scripts/04_enum_halogen_all_classes.py --classes all --k-values 1
   ```

---

## 问题2：黄酮重复问题

### 2.1 问题描述

**现状：**
- 黄酮类已单独处理：`data/output/nplike/Flavone/base.parquet`
- CNPD-ETCM合并库被分类为6个NP类（包括polyphenol）
- **问题：** 黄酮类（Flavonoids）应该属于polyphenol的子类，但：
  - 黄酮没有从合并库中排除
  - 导致黄酮分子被重复枚举（一次作为Flavone，一次作为polyphenol）

**数据量估算：**
- Polyphenol base: 4,047 molecules
- Flavone base: ??? molecules (需要检查)
- **重叠数量：** 未知

### 2.2 解决方案

#### 方案A：排他性分类（推荐）

**原则：** 每个分子仅属于一个主类别，优先级排序

**分类优先级（从高到低）：**
1. **Flavone**（黄酮 - 最具体）
2. **Polyphenol**（多酚 - 中等具体）
3. **Terpenoid**（萜类）
4. **Alkaloid**（生物碱）
5. **Glycoside**（糖苷）
6. **AA_peptide**（肽类）
7. **Lipid**（脂类）
8. **Other**（其他）

**实施流程：**

```python
# 伪代码
def assign_unique_class(molecule_inchikey, all_classes_assigned):
    """
    为每个分子分配唯一的主类别

    Args:
        molecule_inchikey: 分子标识
        all_classes_assigned: {inchikey: [class1, class2, ...]}

    Returns:
        primary_class: str
    """
    classes = all_classes_assigned[molecule_inchikey]

    # 优先级排序
    priority = ['flavone', 'polyphenol', 'terpenoid', 'alkaloid',
                'glycoside', 'aa_peptide', 'lipid', 'other']

    for cls in priority:
        if cls in classes:
            return cls

    return 'other'  # fallback
```

**具体执行：**

```bash
# Step 1: 检查黄酮与polyphenol的重叠
python -c "
import pandas as pd

flavone_df = pd.read_parquet('data/output/nplike/Flavone/base.parquet')
polyphenol_df = pd.read_parquet('data/output/nplike/polyphenol/base.parquet')

flavone_keys = set(flavone_df['inchikey'])
polyphenol_keys = set(polyphenol_df['inchikey'])

overlap = flavone_keys & polyphenol_keys

print(f'Flavone molecules: {len(flavone_keys)}')
print(f'Polyphenol molecules: {len(polyphenol_keys)}')
print(f'Overlap: {len(overlap)}')
print(f'Overlap ratio: {len(overlap)/len(polyphenol_keys)*100:.1f}%')
"

# Step 2: 从polyphenol中排除黄酮
python -c "
import pandas as pd

flavone_df = pd.read_parquet('data/output/nplike/Flavone/base.parquet')
polyphenol_df = pd.read_parquet('data/output/nplike/polyphenol/base.parquet')

flavone_keys = set(flavone_df['inchikey'])

# 排除黄酮后的polyphenol
polyphenol_noflavone = polyphenol_df[~polyphenol_df['inchikey'].isin(flavone_keys)]

print(f'Original polyphenol: {len(polyphenol_df)}')
print(f'After excluding flavones: {len(polyphenol_noflavone)}')

# 保存
polyphenol_noflavone.to_parquet('data/output/nplike/polyphenol/base_no_flavone.parquet')
print('Saved to: data/output/nplike/polyphenol/base_no_flavone.parquet')
"
```

#### 方案B：分层分类（备选）

**原则：** 黄酮作为polyphenol的特殊子类

**优点：** 保持分类层次性
**缺点：** 需要修改数据结构，支持parent_class字段

#### 方案C：多标签+后处理去重（不推荐）

允许重复，在下游分析时通过InChIKey去重。

**缺点：**
- 浪费计算资源
- 下游分析复杂
- 数据冗余

### 2.3 推荐执行方案

**采用方案A：排他性分类**

1. **立即行动：**
   - 检查Flavone/base.parquet的分子数量
   - 统计与polyphenol的重叠率
   - 从polyphenol/base.parquet中排除黄酮分子
   - 生成 `polyphenol/base_no_flavone.parquet`

2. **重新枚举polyphenol:**
   ```bash
   # 使用排除黄酮后的base
   halogenator enum-parquet \
     --input-parquet data/output/nplike/polyphenol/base_no_flavone.parquet \
     --outdir data/output/nplike/polyphenol-1X-clean \
     --k 1 --np-class polyphenol
   ```

3. **全局去重检查：**
   ```python
   # 检查所有类之间是否还有重叠
   all_classes = ['flavone', 'polyphenol', 'terpenoid', 'alkaloid',
                  'glycoside', 'aa_peptide', 'lipid']

   for i, cls1 in enumerate(all_classes):
       for cls2 in all_classes[i+1:]:
           # 检查InChIKey重叠
           # 报告冲突
   ```

---

## 问题3：多重分类冲突

### 3.1 问题分析

**潜在冲突场景：**

| 分子类型 | 可能的多重分类 | 示例 |
|---------|----------------|------|
| **黄酮苷** | Flavone + Glycoside | Quercetin-3-O-glucoside |
| **萜类苷** | Terpenoid + Glycoside | Ginsenosides |
| **芳香萜** | Terpenoid + Polyphenol | Thymol, Carvacrol |
| **生物碱苷** | Alkaloid + Glycoside | Solanine |
| **含氮多酚** | Polyphenol + Alkaloid | 某些吲哚衍生物 |

**冲突率预估：**
- 黄酮苷：可能占glycoside的20-40%
- 萜类苷：可能占glycoside的10-20%
- 其他冲突：<5%

### 3.2 解决方案

#### 方案A：优先级分类（推荐）

**分类决策树：**

```
if has_sugar_moiety:
    if is_flavonoid_aglycone:
        return "Flavone-Glycoside" (单独类)
    elif is_terpenoid_aglycone:
        return "Terpenoid-Glycoside" (or just Glycoside)
    else:
        return "Glycoside"
elif is_flavonoid:
    return "Flavone"
elif is_terpenoid:
    return "Terpenoid"
elif is_alkaloid:
    return "Alkaloid"
elif is_polyphenol:
    return "Polyphenol"
else:
    return "Other"
```

**优先级规则：**
1. **糖苷优先**：含糖基团的分子优先归为Glycoside（因为sugar_mask是特殊处理）
2. **具体性优先**：Flavone > Polyphenol（黄酮是多酚的子类）
3. **结构主导优先**：Alkaloid > Polyphenol（含氮杂环是主导结构特征）

#### 方案B：多标签保留+metadata标记

**原则：** 保留多重分类信息，但在枚举时只执行一次

**实施：**
```python
# 为每个分子保留所有匹配的类别
molecule_metadata = {
    'inchikey': 'XXX',
    'primary_class': 'glycoside',        # 用于枚举
    'all_classes': ['flavone', 'glycoside'],  # 用于分析
    'classification_priority': 'sugar_moiety'  # 分类依据
}
```

**优点：** 保留完整分类信息
**缺点：** 实施复杂度高

### 3.3 推荐方案

**采用方案A + 增强报告：**

1. **执行优先级分类**（如上决策树）

2. **生成分类冲突报告：**
   ```python
   # 输出示例
   Classification Conflicts Report
   ==============================
   Total molecules: 65,000
   Single-class: 58,000 (89.2%)
   Multi-class: 7,000 (10.8%)

   Conflict breakdown:
   - Flavone + Glycoside: 3,500 → Assigned to Glycoside
   - Terpenoid + Glycoside: 2,000 → Assigned to Glycoside
   - Terpenoid + Polyphenol: 1,000 → Assigned to Terpenoid
   - Other conflicts: 500
   ```

3. **为下游分析保留元数据：**
   - 在枚举输出中添加 `secondary_classes` 字段
   - 便于后续分析时按"真实化学分类"聚类

---

## 综合实施计划

### Phase 1: 配置修复与扩展（1-2天）

1. **备份当前状态**
   ```bash
   git commit -m "checkpoint: k=1 enumeration complete with R1+R5 only"
   cp configs/halogen_rules_by_class.yaml configs/halogen_rules_by_class.yaml.backup
   ```

2. **扩展规则集**
   - 修改 `halogen_rules_by_class.yaml`
   - 为所有类添加R3, R4
   - 移除max_sites_per_parent限制
   - 重新启用合理的规则（PRIMARY_OH等）

3. **处理黄酮重复**
   - 统计Flavone与Polyphenol重叠
   - 生成 `polyphenol/base_no_flavone.parquet`

### Phase 2: POC验证（1天）

1. **小规模测试扩展规则集**
   ```bash
   python scripts/03_enum_halogen_poc.py --class terpenoid --k 1 --max-parents 1000
   python scripts/03_enum_halogen_poc.py --class polyphenol --k 1 --max-parents 1000
   ```

2. **检查产物规模**
   - 确认新规则（R3, R4）是否产生产物
   - 确认无组合爆炸
   - 记录products/parent比例变化

### Phase 3: 全量重新枚举（3-5天）

1. **清理旧输出**
   ```bash
   rm -rf data/output/nplike/*-1X
   ```

2. **重新运行k=1**（使用扩展规则集）
   ```bash
   python scripts/04_enum_halogen_all_classes.py --classes all --k-values 1
   ```

3. **QC检查**
   - 验证所有规则都产生了产物
   - 统计规则分布
   - 对比扩展前后的产物数量变化

### Phase 4: 多重分类处理（2-3天）

1. **全局冲突检测**
   - 编写脚本检查所有类之间的InChIKey重叠
   - 生成冲突矩阵

2. **应用优先级分类**
   - 根据决策树重新分配冲突分子
   - 生成干净的base.parquet（每个分子唯一类别）

3. **重新枚举（如需要）**

### Phase 5: k=2枚举（时间待定）

使用验证通过的配置运行k=2全量枚举。

---

## 附录A：完整配置示例（Terpenoid）

```yaml
terpenoid:
  description: "Isoprene-derived polycyclic natural products"
  molecule_count: 28606
  enabled: true

  k1:
    include_rules:
      - RING_SP2__CH__TO__X           # R1 芳环
      - RING_SP3__CH__TO__X           # 脂环
      - COOH__TO__CX                  # R5 羧基
      - R3                            # NEW: OH卤代
      - R4                            # NEW: 醛基卤代
      - ALPHA_CARBONYL__CH2__TO__X    # α-羰基
      - PRIMARY_OH__CH2OH__TO__X      # 伯醇
    sugar_mask: false

  k2:
    include_rules:
      - RING_SP2__CH__TO__X
      - RING_SP3__CH__TO__X
      - COOH__TO__CX
      - R3
      - R4
      - ALPHA_CARBONYL__CH2__TO__X
      - PRIMARY_OH__CH2OH__TO__X
    sugar_mask: false

  per_rule_overrides:
    # 仅保留极端safety nets
    R3:
      max_sites_per_parent: 10        # OH较多
    R4:
      max_sites_per_parent: 3         # 醛基较少
    RING_SP3__CH__TO__X:
      max_sites_per_parent: 8         # 环较多
    PRIMARY_OH__CH2OH__TO__X:
      enabled: true
      max_sites_per_parent: 5
```

---

## 附录B：脚本工具

### B.1 检查黄酮重叠

```python
# scripts/check_flavone_overlap.py
import pandas as pd

flavone = pd.read_parquet('data/output/nplike/Flavone/base.parquet')
polyphenol = pd.read_parquet('data/output/nplike/polyphenol/base.parquet')

flavone_keys = set(flavone['inchikey'])
polyphenol_keys = set(polyphenol['inchikey'])
overlap = flavone_keys & polyphenol_keys

print(f'Flavone: {len(flavone_keys):,}')
print(f'Polyphenol: {len(polyphenol_keys):,}')
print(f'Overlap: {len(overlap):,} ({len(overlap)/len(polyphenol_keys)*100:.1f}%)')

# 排除黄酮
polyphenol_clean = polyphenol[~polyphenol['inchikey'].isin(flavone_keys)]
polyphenol_clean.to_parquet('data/output/nplike/polyphenol/base_no_flavone.parquet')
print(f'Clean polyphenol: {len(polyphenol_clean):,}')
```

### B.2 全局冲突检测

```python
# scripts/detect_classification_conflicts.py
import pandas as pd
from itertools import combinations

classes = ['flavone', 'polyphenol', 'terpenoid', 'alkaloid',
           'glycoside', 'aa_peptide', 'lipid']

# 加载所有base.parquet
class_molecules = {}
for cls in classes:
    path = f'data/output/nplike/{cls}/base.parquet'
    try:
        df = pd.read_parquet(path)
        class_molecules[cls] = set(df['inchikey'])
    except:
        print(f'Warning: {cls} not found')

# 检测冲突
print('Classification Conflicts Matrix')
print('='*60)

for cls1, cls2 in combinations(classes, 2):
    if cls1 in class_molecules and cls2 in class_molecules:
        overlap = class_molecules[cls1] & class_molecules[cls2]
        if overlap:
            pct1 = len(overlap) / len(class_molecules[cls1]) * 100
            pct2 = len(overlap) / len(class_molecules[cls2]) * 100
            print(f'{cls1:<15} vs {cls2:<15}: {len(overlap):>5} '
                  f'({pct1:>5.1f}% / {pct2:>5.1f}%)')
```

---

## 后续会话衔接要点

1. **当前状态：** k=1枚举完成，但仅使用R1+R5
2. **待执行任务：**
   - [ ] 扩展规则集（添加R3, R4等）
   - [ ] 处理黄酮重复问题
   - [ ] 解决多重分类冲突
   - [ ] 重新运行k=1枚举
   - [ ] 继续k=2枚举
3. **关键文件：**
   - 配置：`configs/halogen_rules_by_class.yaml`
   - 脚本：`scripts/04_enum_halogen_all_classes.py`
   - 数据：`data/output/nplike/*/base.parquet`
4. **Git状态：** commit ebd0967

---

**文档结束**
