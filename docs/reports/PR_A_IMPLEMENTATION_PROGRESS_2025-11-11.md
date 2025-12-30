# PR-A Implementation Progress Report
**Date**: 2025-11-11
**Status**: Hotfixes完成，PR-A核心实现待继续

---

## ✅ 已完成的任务

### Hotfix 1-5：可视化层稳定性修复

#### ✅ Hotfix 1: cmd_html中diag_df访问的NameError修复
**文件**: `scripts/09_visualize_library.py` (lines 1340-1360)

**修改**:
```python
# 修复前：diagnostics为空时不会创建diag_df，后续访问会NameError
if diagnostics and SITE_FINDER_AVAILABLE:
    diag_df = pd.DataFrame(diagnostics)
    ...

# 修复后：统一创建diag_df，即使为空DataFrame
diag_df = pd.DataFrame(diagnostics) if (SITE_FINDER_AVAILABLE and diagnostics) else pd.DataFrame()

if not diag_df.empty:
    # 安全访问
    ...
```

#### ✅ Hotfix 2: Grid模式添加site_finder缺失时的兜底
**文件**: `scripts/09_visualize_library.py` (lines 1100-1115)

**修改**:
```python
# 添加了elif分支处理SITE_FINDER_AVAILABLE=False的情况
if highlight_sites and SITE_FINDER_AVAILABLE:
    # 智能高亮
    ...
elif highlight_sites:
    # 打印警告（每页只打印一次）
    if idx == start_idx:
        logger.warning("site_finder module not available, skipping site highlighting")
```

#### ✅ Hotfix 3: 统一清理脏列名（空格）
**文件**: `scripts/09_visualize_library.py`

**修改1 - 扩展normalize_columns函数** (lines 141-182):
```python
def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    增强：先清理所有列名中的空格（如 'x f_site_atoms' → 'xf_site_atoms'）
    然后应用规范化映射
    """
    rename_map = {}

    # Hotfix 3: 首先移除所有列名中的空格
    for col in df.columns:
        if ' ' in col:
            cleaned = col.replace(' ', '')
            rename_map[col] = cleaned

    # 然后应用COLUMN_CANON映射
    ...
```

**修改2 - 在三个命令中调用normalize_columns**:
- cmd_grid (line 1070): `df = normalize_columns(df)`
- cmd_html (line 1224): `df = normalize_columns(df)`
- cmd_sprite (line 1616): `df = normalize_columns(df)`

#### ✅ Hotfix 4: 清理重复的diagnostic透传逻辑
**文件**: `scripts/09_visualize_library.py` (line 867)

**修改**:
```python
# 删除前：在render_params中存储_diagnostic（多余）
render_params['_diagnostic'] = diagnostic

# 删除后：只在result字典中存储diagnostic（line 877）
result['diagnostic'] = diagnostic
```

#### ✅ Hotfix 5: Windows控制台编码问题
**状态**: 已在之前测试脚本中修复（特殊符号如✓改为[OK]），logger输出无特殊字符。

---

## 🔄 正在进行/待完成的任务

### PR-A核心实现：产后回贴777映射号

**背景**：RDKit的`RunReactants()`不保留产物端的原子映射号。需要在反应后手动标记需要高亮的原子。

#### 📝 待完成任务清单

##### **Step 1: 更新transforms.yaml添加highlight_mapnums**
**文件**: `configs/transforms.yaml`

**待添加内容**:
```yaml
transforms:
  - name: OH_to_OMe
    label: "OH->OMe"
    query_smarts: "[c:1][OX2H:2]"
    smirks: "[c:1][OX2H:2]>>[c:1][O:2]C"
    # 新增字段：
    highlight_mapnums: [2]              # 产品模板中映射号2的原子（氧）
    product_highlight_symbol: "O"       # 验证用：期望的原子符号

  - name: OH_to_NH2
    label: "OH->NH2"
    query_smarts: "[c:1][OX2H:2]"
    smirks: "[c:1][OX2H:2]>>[c:1][N:2]([H])[H]"
    # 新增字段：
    highlight_mapnums: [2]              # 产品模板中映射号2的原子（氮）
    product_highlight_symbol: "N"       # 验证用：期望的原子符号
```

**说明**：
- `highlight_mapnums`: 告诉下游"产品模板中哪个映射号代表需要高亮的位点"
- `product_highlight_symbol`: 可选，用于验证（确保标记的是正确的原子类型）

---

##### **Step 2: 实现annotate_product_sites_with_mapnum函数**
**文件**: `scripts/08_transform_library_v2.py`

**待添加函数**（插入到文件顶部，在TransformationEngineV2类之前）:

```python
def annotate_product_sites_with_mapnum(
    rxn: AllChem.ChemicalReaction,
    prod_mol: Chem.Mol,
    highlight_mapnums: List[int],
    target_symbol: Optional[str] = None,
    stamp_mapnum: int = 777
) -> List[int]:
    """
    用产品模板在产物上定位highlight_mapnums对应原子，标记为stamp_mapnum。

    原理：
    1. 获取反应的产品模板（保留了映射号）
    2. 用产品模板在实际产物上做子结构匹配
    3. 找到模板中映射号==highlight_mapnums的原子
    4. 将对应的产物原子标记为stamp_mapnum (777)

    Args:
        rxn: RDKit反应对象
        prod_mol: 产物分子（来自rxn.RunReactants()）
        highlight_mapnums: 要高亮的映射号列表（如[2]）
        target_symbol: 可选，验证原子符号（如'O'或'N'）
        stamp_mapnum: 要标记的映射号（默认777）

    Returns:
        被标记的原子索引列表
    """
    marked_atoms = set()

    # 遍历所有产品模板（通常只有一个）
    for pidx in range(rxn.GetNumProductTemplates()):
        prod_template = rxn.GetProductTemplate(pidx)

        # 在模板中找到目标映射号的原子索引
        target_template_idxs = [
            atom.GetIdx()
            for atom in prod_template.GetAtoms()
            if atom.GetAtomMapNum() in highlight_mapnums
        ]

        if not target_template_idxs:
            continue

        # 用模板匹配产物
        matches = prod_mol.GetSubstructMatches(
            prod_template,
            useChirality=True,
            maxMatches=64,
            uniquify=True
        )

        for match in matches:
            for tmpl_idx in target_template_idxs:
                prod_idx = match[tmpl_idx]
                atom = prod_mol.GetAtomWithIdx(prod_idx)

                # 验证原子符号（如果指定）
                if target_symbol and atom.GetSymbol() != target_symbol:
                    continue

                # 标记原子
                atom.SetAtomMapNum(stamp_mapnum)
                marked_atoms.add(prod_idx)

    return sorted(marked_atoms)
```

**关键说明**:
- **产品模板**: `rxn.GetProductTemplate(0)` 保留了SMIRKS中的映射号
- **子结构匹配**: 将模板映射到实际产物分子
- **映射号转移**: 找到模板中`:2`的原子，在产物中对应的原子上标记`:777`

---

##### **Step 3: 修改08_transform取消位点上限**
**文件**: `scripts/08_transform_library_v2.py`

**位置1 - TransformationEngineV2.__init__()** (约line 386):
```python
# 修改前：
self.max_sites = self.defaults.get('max_sites_per_mol', 4)

# 修改后：
self.max_sites = self.defaults.get('max_sites_per_mol', -1)  # -1表示不限制
```

**位置2 - find_sites()方法** (约line 410):
```python
# 修改前：
if len(matches) > self.max_sites:
    matches = matches[:self.max_sites]

# 修改后：
if self.max_sites > 0 and len(matches) > self.max_sites:
    matches = matches[:self.max_sites]
```

**位置3 - 读取配置** (约line 356-392):
```python
# 在__init__中读取highlight_mapnums
self.highlight_mapnums = self.transform.get('highlight_mapnums', [])
self.product_highlight_symbol = self.transform.get('product_highlight_symbol')
```

---

##### **Step 4: 集成产后回贴777到apply_to_molecule**
**文件**: `scripts/08_transform_library_v2.py`
**方法**: `TransformationEngineV2.apply_to_molecule()`

**修改位置** (约line 477-511，在产物循环中):

```python
try:
    prod_mol = rxn_products[site_idx][0]

    # Sanitize
    Chem.SanitizeMol(prod_mol)

    # ===== PR-A 核心修改：产后回贴777 =====
    # 使用产品模板匹配，标记需要高亮的原子
    marked_atoms = []
    if self.highlight_mapnums:
        marked_atoms = annotate_product_sites_with_mapnum(
            self.reaction,
            prod_mol,
            highlight_mapnums=self.highlight_mapnums,
            target_symbol=self.product_highlight_symbol,
            stamp_mapnum=777
        )

    # Get canonical SMILES (fast)
    prod_smiles_canon = canonical_smiles(prod_mol)
    if not prod_smiles_canon:
        continue

    # ===== NEW: Get mapped SMILES (preserve :777) =====
    prod_smiles_mapped = Chem.MolToSmiles(prod_mol, canonical=False)

    # InChIKey generation
    prod_inchikey = get_inchikey(prod_mol)

    # Quick properties only
    props = calculate_properties_quick(prod_mol)

    # Calculate halogen fields
    halogen_fields = calculate_halogen_fields(prod_mol, self.transform['label'])

    products.append({
        'smiles': prod_smiles_canon,
        'canonical_smiles': prod_smiles_canon,
        'inchikey': prod_inchikey,
        'xf_success': True,
        'xf_reason': None,
        'xf_label': self.transform['label'],
        'xf_rule_id': self.transform['name'],
        'xf_site_index': site_idx,
        'xf_site_atoms': str(match),
        # ===== NEW FIELDS for PR-A =====
        'product_mapped_smiles': prod_smiles_mapped,         # 包含:777的SMILES
        'xf_site_mapnums': json.dumps([777] * len(marked_atoms)),  # [777]或[777, 777]
        'source_smiles': smiles,
        **source_data,
        **props,
        **halogen_fields
    })

except Exception as e:
    logger.debug(f"Site {site_idx} failed: {e}")
    continue
```

**关键说明**:
- `product_mapped_smiles`: 用`canonical=False`保留原子映射号
- `xf_site_mapnums`: JSON列表，记录使用的映射号（通常是[777]）

---

##### **Step 5: 小规模测试(10分子)验证映射**

**创建测试脚本**: `scripts/test_pr_a_mapping.py`

```python
#!/usr/bin/env python3
"""测试PR-A产后回贴777功能"""
import sys
sys.path.insert(0, 'scripts')

from rdkit import Chem
from rdkit.Chem import AllChem
import yaml

# 加载配置
with open('configs/transforms.yaml', 'r', encoding='utf-8') as f:
    config = yaml.safe_load(f)

# 找到OH->OMe转化
transform = None
for xf in config['transforms']:
    if xf['name'] == 'OH_to_OMe':
        transform = xf
        break

print(f"Transform: {transform['label']}")
print(f"SMIRKS: {transform['smirks']}")
print(f"Highlight mapnums: {transform.get('highlight_mapnums', [])}")

# 创建反应
rxn = AllChem.ReactionFromSmarts(transform['smirks'])

# 测试分子：苯酚
test_mol = Chem.MolFromSmiles('c1ccccc1O')
print(f"\nTest molecule: c1ccccc1O (phenol)")

# 运行反应
products = rxn.RunReactants((test_mol,))
print(f"Reaction produced {len(products)} product(s)")

if products:
    prod_mol = products[0][0]
    Chem.SanitizeMol(prod_mol)

    # 手动标记（模拟annotate_product_sites_with_mapnum）
    prod_template = rxn.GetProductTemplate(0)
    print(f"\nProduct template atoms:")
    for atom in prod_template.GetAtoms():
        print(f"  Atom {atom.GetIdx()}: {atom.GetSymbol()} mapnum={atom.GetAtomMapNum()}")

    # 匹配
    matches = prod_mol.GetSubstructMatches(prod_template, useChirality=True)
    print(f"\nMatches: {len(matches)}")

    if matches:
        match = matches[0]
        print(f"Match: {match}")

        # 找到映射号2的原子
        highlight_mapnums = transform.get('highlight_mapnums', [])
        for tmpl_idx, tmpl_atom in enumerate(prod_template.GetAtoms()):
            if tmpl_atom.GetAtomMapNum() in highlight_mapnums:
                prod_idx = match[tmpl_idx]
                print(f"\nMarking product atom {prod_idx} (was template atom {tmpl_idx} with mapnum {tmpl_atom.GetAtomMapNum()})")
                prod_mol.GetAtomWithIdx(prod_idx).SetAtomMapNum(777)

        # 检查结果
        mapped_smiles = Chem.MolToSmiles(prod_mol, canonical=False)
        print(f"\nMapped SMILES: {mapped_smiles}")

        if ':777' in mapped_smiles:
            print("[SUCCESS] :777 found in mapped SMILES!")
        else:
            print("[FAIL] :777 NOT found in mapped SMILES")

        # 验证可以提取
        for atom in prod_mol.GetAtoms():
            if atom.GetAtomMapNum() == 777:
                print(f"[OK] Atom {atom.GetIdx()} ({atom.GetSymbol()}) has mapnum=777")
```

**运行**: `python scripts/test_pr_a_mapping.py`

**预期输出**:
```
[SUCCESS] :777 found in mapped SMILES!
[OK] Atom X (O) has mapnum=777
```

---

##### **Step 6-10: 产库任务（4个派生库）**

**待执行命令**（需在Step 1-5完成后）:

```bat
REM 确保基底库存在
REM 基底库位置：E:\Projects\halogenator\data\base\Flavone-1X.parquet
REM                E:\Projects\halogenator\data\base\Flavone-2X.parquet

REM 1) Flavone-1X + OH->OMe
python scripts\08_transform_library_v2.py apply ^
  -c configs\transforms.yaml ^
  -i E:\Projects\halogenator\data\base\Flavone-1X.parquet ^
  -o E:\Projects\halogenator\data\derived\Flavone-1X-OMe_full.parquet ^
  --transform OH_to_OMe ^
  --workers 4

REM 2) Flavone-1X + OH->NH2
python scripts\08_transform_library_v2.py apply ^
  -c configs\transforms.yaml ^
  -i E:\Projects\halogenator\data\base\Flavone-1X.parquet ^
  -o E:\Projects\halogenator\data\derived\Flavone-1X-NH2_full.parquet ^
  --transform OH_to_NH2 ^
  --workers 4

REM 3) Flavone-2X + OH->OMe
python scripts\08_transform_library_v2.py apply ^
  -c configs\transforms.yaml ^
  -i E:\Projects\halogenator\data\base\Flavone-2X.parquet ^
  -o E:\Projects\halogenator\data\derived\Flavone-2X-OMe_full.parquet ^
  --transform OH_to_OMe ^
  --workers 4

REM 4) Flavone-2X + OH->NH2
python scripts\08_transform_library_v2.py apply ^
  -c configs\transforms.yaml ^
  -i E:\Projects\halogenator\data\base\Flavone-2X.parquet ^
  -o E:\Projects\halogenator\data\derived\Flavone-2X-NH2_full.parquet ^
  --transform OH_to_NH2 ^
  --workers 4
```

**验证**:
- 检查输出parquet文件包含 `product_mapped_smiles` 和 `xf_site_mapnums` 列
- 随机抽查几行，确保 `product_mapped_smiles` 包含 `:777`
- 确保 `xf_site_mapnums` 为 `[777]`

---

##### **Step 11: 可视化验收**

**命令**（对4个派生库分别执行）:

```bat
REM 1) Flavone-1X-OMe 可视化
python scripts\09_visualize_library.py html ^
  -i E:\Projects\halogenator\data\derived\Flavone-1X-OMe_full.parquet ^
  -o E:\Projects\halogenator\data\viz_pr_a\Flavone-1X-OMe_gallery.html ^
  --highlight-sites --preset hq --workers 4

REM 2) Flavone-1X-NH2 可视化
python scripts\09_visualize_library.py html ^
  -i E:\Projects\halogenator\data\derived\Flavone-1X-NH2_full.parquet ^
  -o E:\Projects\halogenator\data\viz_pr_a\Flavone-1X-NH2_gallery.html ^
  --highlight-sites --preset hq --workers 4

REM 3) Flavone-2X-OMe 可视化
python scripts\09_visualize_library.py html ^
  -i E:\Projects\halogenator\data\derived\Flavone-2X-OMe_full.parquet ^
  -o E:\Projects\halogenator\data\viz_pr_a\Flavone-2X-OMe_gallery.html ^
  --highlight-sites --preset hq --workers 4

REM 4) Flavone-2X-NH2 可视化
python scripts\09_visualize_library.py html ^
  -i E:\Projects\halogenator\data\derived\Flavone-2X-NH2_full.parquet ^
  -o E:\Projects\halogenator\data\viz_pr_a\Flavone-2X-NH2_gallery.html ^
  --highlight-sites --preset hq --workers 4
```

**检查点**:
- 查看 `viz_diagnostics.csv` 中的 `method_used` 列
  - 预期：大部分为 `'mapnum'`（使用原子映射）
  - 高置信度应达到 ~98%
- 打开HTML gallery，随机检查几个分子
  - 确认高亮仅在新引入的甲氧基/氨基上
  - 不误高亮桥氧、羰基氧、或基底中已有的OMe/NH2

---

##### **Step 12: 最终验收**

**8个问题样本回归测试**:

这8个样本在新数据中应该都获得High confidence + method_used='mapnum':
- mol_000012
- mol_000027
- mol_000029
- mol_000006
- mol_000007
- mol_000008
- mol_000017
- mol_000010

**验收标准**:
1. ✅ **映射优先率**: method_used='mapnum' ≥ 95%
2. ✅ **高置信度率**: confidence='high' ≥ 98%
3. ✅ **零失败率**: confidence='none' = 0%
4. ✅ **视觉准确**: 手工检查20个样本，高亮位置100%正确
5. ✅ **一致性**: HTML/Grid/Sprite三端风格一致

---

## 📂 关键文件位置

### 修改的文件（Hotfixes已完成）:
- ✅ `scripts/09_visualize_library.py`
  - Lines 141-182: normalize_columns() 增强
  - Lines 1070, 1224, 1616: 调用normalize_columns()
  - Line 1100-1115: Grid模式兜底
  - Line 1340-1360: diag_df安全处理
  - Line 864: 删除重复diagnostic

### 待修改的文件（PR-A核心）:
- ⏳ `configs/transforms.yaml` - 添加highlight_mapnums字段
- ⏳ `scripts/08_transform_library_v2.py` - 实现产后回贴777
  - 新增函数: annotate_product_sites_with_mapnum()
  - 修改方法: TransformationEngineV2.__init__()
  - 修改方法: find_sites()
  - 修改方法: apply_to_molecule()

### 新增的文件:
- ✅ `scripts/site_finder.py` (632行) - 已完成
- ✅ `scripts/test_site_finder.py` - 已完成
- ⏳ `scripts/test_pr_a_mapping.py` - 待创建

---

## 🎯 下一步行动

### 立即执行（顺序很重要）:
1. **Step 1**: 修改 `configs/transforms.yaml`
2. **Step 2**: 在 `08_transform_library_v2.py` 中实现 `annotate_product_sites_with_mapnum()`
3. **Step 3**: 修改 `08_transform_library_v2.py` 取消位点上限
4. **Step 4**: 修改 `apply_to_molecule()` 集成产后回贴
5. **Step 5**: 运行小规模测试验证
6. **Step 6-9**: 执行4个产库命令（~1-2小时）
7. **Step 10**: 可视化验收
8. **Step 11**: 最终验收

### 预计时间:
- Step 1-5 (代码实现): ~30-60分钟
- Step 6-9 (产库运行): ~1-2小时（取决于数据量）
- Step 10-11 (验收): ~30分钟

### 总计: ~2-4小时可完成整个PR-A实现

---

## 📊 当前状态总结

| 阶段 | 状态 | 完成度 |
|------|------|--------|
| PR-B~E Hotfixes | ✅ 完成 | 100% |
| PR-A设计 | ✅ 完成 | 100% |
| PR-A代码实现 | ⏳ 待执行 | 0% |
| 产库执行 | ⏳ 待执行 | 0% |
| 可视化验收 | ⏳ 待执行 | 0% |
| 最终验收 | ⏳ 待执行 | 0% |

**Overall Progress**: ~35% (Hotfixes + 设计完成)

---

## 📝 后续对话需要的信息

如果上下文中断，下一个Claude实例需要知道：

1. **Hotfixes已全部完成**（1-5），无需重做
2. **PR-A核心逻辑已设计完成**，按本文档的Step 1-12执行即可
3. **关键决策**：
   - 使用产品模板匹配 + 手动SetAtomMapNum(777)
   - 映射号777作为标准高亮标记
   - 新字段：product_mapped_smiles, xf_site_mapnums
4. **执行顺序严格**：必须先完成Step 1-5（代码），再执行Step 6-9（产库）
5. **验收标准**：映射优先率≥95%，高置信度≥98%，零失败

---

**报告生成时间**: 2025-11-11 15:40
**预计完成时间**: 2025-11-11 18:00 (如连续工作)
