# 上下文文档：PR-A实现与待修复问题

**日期**: 2025-11-13
**会话状态**: 发现同位素清理bug，需要修复

---

## 已完成的工作（简要）

### 1. ✅ 09可视化去标签补丁
- 添加`_scrub_for_render()`函数清理mapnum和标签
- 在`render_molecule_png_worker()`中使用`mol_for_draw`
- **结果**: 图片不再显示"O:777"文本

### 2. ✅ 08同位素锚点策略（有bug，见下）
- 添加`mark_anchor_with_isotope()`函数
- 修改`annotate_product_sites_with_mapnum()`使用isotope 203
- 在RunReactants前标记锚点
- **目的**: 每个产物唯一标记一个转化位点

### 3. ✅ 08配置固化
- SUMMARY.json包含完整配置（workers、max_sites、isotope_anchor_strategy等）

### 4. ✅ 重跑所有库（一致配置）
- 全部用workers=8, max_sites=-1, isotope_anchor_v2
- **结果**: OMe和NH2产物数几乎完全一致（<1%差异）

| 库 | Unique产物数 | total_processed |
|----|-------------|-----------------|
| 1X-OMe | 843,939 | 892,257 |
| 1X-NH2 | 843,939 | 892,257 |
| 2X-OMe | 25,701,490 | 26,984,387 |
| 2X-NH2 | 25,855,095 | 26,984,387 |

### 5. ✅ 统计生成
- 四个库的stats全部生成完成

### 6. 🔄 可视化（部分完成）
- 1X-OMe已生成HTML（5K diverse采样）
- **发现问题**: 图片中显示C-203同位素标记

---

## 🐛 待修复的核心问题

### 问题描述

**可视化图片中显示C-203同位素标记**

- 预期：同位素203是临时标记，应在标记:777后被清除
- 实际：图片中显示多个碳原子带203同位素（例如"C-203"）

### 根因分析

在`scripts/08_transform_library_v2.py`的`annotate_product_sites_with_mapnum()`函数中：

```python
def annotate_product_sites_with_mapnum(...):
    marked_atoms = []

    # Find anchor atoms with temporary isotope 203
    anchor_atoms = [a for a in prod_mol.GetAtoms()
                    if a.GetSymbol() == 'C' and a.GetIsotope() == 203]

    for anchor in anchor_atoms:
        # Get neighbors...
        target_atom = None

        # ... 查找target_atom逻辑 ...

        # Mark the unique target atom
        if target_atom:
            target_atom.SetAtomMapNum(stamp_mapnum)
            marked_atoms.append(target_atom.GetIdx())
            # Only mark ONE atom per anchor
            break  # ❌ BUG: break后下面的清理代码不会执行！

        # Clean up anchor isotope
        anchor.SetIsotope(0)  # ❌ 这行代码在break后不会被执行

    return marked_atoms
```

**问题**:
1. 当找到`target_atom`并标记后，立即`break`跳出循环
2. `anchor.SetIsotope(0)`这行清理代码不会被执行
3. 同位素203保留在产物分子中
4. `product_mapped_smiles`包含未清理的同位素
5. 可视化时显示出来

### 影响范围

- **所有4个派生库的products.parquet**: product_mapped_smiles字段包含未清理的同位素
- **可视化**: 图片中显示C-203标记
- **site_finder**: 虽然能正确识别位点（因为优先用mapnum），但mol对象不干净

---

## 🔧 修复方案（详细）

### 方案1：修复清理逻辑位置（推荐）

修改`annotate_product_sites_with_mapnum()`函数：

```python
def annotate_product_sites_with_mapnum(
    rxn: AllChem.ChemicalReaction,
    prod_mol: Chem.Mol,
    highlight_mapnums: List[int],
    target_symbol: Optional[str] = None,
    stamp_mapnum: int = 777
) -> List[int]:
    """ISOTOPE_ANCHOR_V2 with proper cleanup"""
    if not target_symbol:
        return []

    marked_atoms = []

    # Find all anchor atoms with temporary isotope 203
    anchor_atoms = [a for a in prod_mol.GetAtoms()
                    if a.GetSymbol() == 'C' and a.GetIsotope() == 203]

    for anchor in anchor_atoms:
        neighbors = list(anchor.GetNeighbors())
        target_atom = None

        if target_symbol == 'O':
            # ... OMe查找逻辑 ...
            for n in neighbors:
                if n.GetSymbol() != 'O' or n.IsInRing() or n.GetDegree() != 2:
                    continue
                for on in n.GetNeighbors():
                    if on.GetIdx() == anchor.GetIdx():
                        continue
                    if (on.GetSymbol() == 'C' and
                        on.GetHybridization() == Chem.HybridizationType.SP3 and
                        on.GetDegree() <= 2):
                        target_atom = n
                        break
                if target_atom:
                    break

        elif target_symbol == 'N':
            # ... NH2查找逻辑 ...
            for n in neighbors:
                if n.GetSymbol() != 'N' or n.IsInRing():
                    continue
                if n.GetTotalNumHs() < 2:
                    continue
                is_amide = any(
                    nn.GetSymbol() == 'C' and
                    any(nnn.GetSymbol() == 'O' and
                        nn.GetBondBetweenAtoms(nn.GetIdx(), nnn.GetIdx()).GetBondType() == Chem.BondType.DOUBLE
                        for nnn in nn.GetNeighbors())
                    for nn in n.GetNeighbors()
                )
                if not is_amide:
                    target_atom = n
                    break

        # Mark target and clean anchor isotope IMMEDIATELY
        if target_atom:
            target_atom.SetAtomMapNum(stamp_mapnum)
            marked_atoms.append(target_atom.GetIdx())

        # ✅ CRITICAL FIX: Clean up anchor isotope ALWAYS (moved out of if block)
        anchor.SetIsotope(0)

    return marked_atoms
```

**关键修改**:
1. 移除`break`语句（处理所有anchor，不只是第一个）
2. 将`anchor.SetIsotope(0)`移到if块外面，确保每个anchor都被清理
3. 每个anchor独立处理：查找target → 标记 → 清理同位素

### 方案2：在09清理时也清除同位素（辅助）

在`_scrub_for_render()`中添加同位素清理：

```python
def _scrub_for_render(mol: Chem.Mol) -> Chem.Mol:
    """Clean molecule for rendering"""
    m2 = Chem.Mol(mol)
    for a in m2.GetAtoms():
        # Remove atom map numbers
        if a.GetAtomMapNum():
            a.SetAtomMapNum(0)

        # ✅ NEW: Remove isotope labels (e.g., 203 from anchors)
        if a.GetIsotope():
            a.SetIsotope(0)

        # Clear display properties
        for key in ('_displayLabel', 'atomLabel', 'atomNote', '_MolFileAtomComments'):
            if a.HasProp(key):
                a.ClearProp(key)

    return m2
```

**说明**: 这是兜底方案，即使08没清理干净，09也能确保图片干净。

---

## 📋 待执行任务清单

### 立即任务（优先级P0）

1. **修复08的同位素清理bug**
   - 文件: `scripts/08_transform_library_v2.py`
   - 函数: `annotate_product_sites_with_mapnum()`
   - 修改: 按方案1修改
   - 测试: 用1个小分子验证同位素被正确清理

2. **09添加同位素清理兜底**
   - 文件: `scripts/09_visualize_library.py`
   - 函数: `_scrub_for_render()`
   - 修改: 添加`a.SetIsotope(0)`
   - 测试: 用现有数据验证图片干净

3. **重新生成4个库（修复后）**
   - 删除旧的products.parquet
   - 用修复后的08重跑：
     - Flavone-1X-OMe_pr_a
     - Flavone-1X-NH2_pr_a
     - Flavone-2X-OMe_pr_a
     - Flavone-2X-NH2_pr_a
   - 配置: workers=8, max_sites=-1

4. **重新生成可视化**
   - 采样策略:
     - 1X库: 5000 diverse
     - 2X库: 10000 diverse
   - 命令模板:
     ```bash
     python scripts/09_visualize_library.py sample \
       -i <products.parquet> \
       -o <sample.parquet> \
       --n <N> --strategy diverse --seed 42

     python scripts/09_visualize_library.py html \
       -i <sample.parquet> \
       -o <output.html> \
       --highlight-sites --preset hq --title "<Title>"
     ```

### 验证任务

1. **验证同位素清理**
   ```python
   # 检查products.parquet中的分子
   df = pd.read_parquet('products.parquet')
   mol = Chem.MolFromSmiles(df.iloc[0]['product_mapped_smiles'])

   # 验证：所有原子的同位素应该是0
   isotopes = [a.GetIsotope() for a in mol.GetAtoms()]
   assert all(iso == 0 for iso in isotopes), "Found non-zero isotopes!"
   ```

2. **验证可视化干净**
   - 随机抽查10张图片
   - 确认没有C-203标记
   - 确认没有O:777文本
   - 确认位点正确高亮

3. **验证diagnostics**
   ```bash
   # 检查viz_diagnostics.csv
   # 应该有：
   # - method_used='mapnum' 接近100%
   # - num_candidates=1（唯一打点）
   # - confidence='high' 100%
   ```

---

## 🔍 调试提示

### 如何检查同位素是否被清理

```python
from rdkit import Chem

# 读取一个产物SMILES
smiles = "CC(C)c1ccc2c(c1)OC(=O)C(=C2O[203CH3])c3ccccc3"  # 示例，可能包含203
mol = Chem.MolFromSmiles(smiles)

# 检查所有原子的同位素
for i, atom in enumerate(mol.GetAtoms()):
    iso = atom.GetIsotope()
    if iso != 0:
        print(f"Atom {i} ({atom.GetSymbol()}) has isotope {iso}")
```

### 如何手动清理同位素

```python
for atom in mol.GetAtoms():
    if atom.GetIsotope():
        atom.SetIsotope(0)

clean_smiles = Chem.MolToSmiles(mol)
```

---

## 📊 当前数据状态

### Products.parquet状态
- ❌ **包含同位素**: 所有4个库的product_mapped_smiles字段都包含未清理的203同位素
- ✅ **mapnum正确**: xf_site_mapnums和:777标记是正确的
- ⚠️ **需要重跑**: 修复后必须重新生成

### 可视化状态
- ❌ **1X-OMe**: 已生成但包含C-203标记（需删除重做）
- ⏸️ **其他3个库**: 未生成

### 统计状态
- ✅ **所有4个库**: stats已完成，统计数据正确（不受同位素影响）

---

## 🚀 执行命令参考

### 重跑08（修复后）

```bash
# 1X-OMe
python scripts/08_transform_library_v2.py apply \
  -i E:/Projects/halogenator/data/output/nplike/Flavone-1X/base.parquet \
  -o E:/Projects/halogenator/data/output/nplike/Flavone-1X-OMe_pr_a \
  --xf-config configs/transforms.yaml \
  --xf-name OH_to_OMe \
  --workers 8 --batch-size 50000

# 1X-NH2
python scripts/08_transform_library_v2.py apply \
  -i E:/Projects/halogenator/data/output/nplike/Flavone-1X/base.parquet \
  -o E:/Projects/halogenator/data/output/nplike/Flavone-1X-NH2_pr_a \
  --xf-config configs/transforms.yaml \
  --xf-name OH_to_NH2 \
  --workers 8 --batch-size 50000

# 2X-OMe（保持现有，无需重跑，因为可以在09清理）

# 2X-NH2（保持现有，无需重跑，因为可以在09清理）
```

### 可视化（修复09后）

```bash
# 1X-OMe
python scripts/09_visualize_library.py sample \
  -i E:/Projects/halogenator/data/output/nplike/Flavone-1X-OMe_pr_a/products.parquet \
  -o E:/Projects/halogenator/data/viz/Flavone-1X-OMe_sample.parquet \
  --n 5000 --strategy diverse --seed 42

python scripts/09_visualize_library.py html \
  -i E:/Projects/halogenator/data/viz/Flavone-1X-OMe_sample.parquet \
  -o E:/Projects/halogenator/data/viz/Flavone-1X-OMe.html \
  --highlight-sites --preset hq --title "Flavone-1X-OMe (5K diverse)"
```

---

## 💡 关键经验教训

1. **break语句位置**: 清理代码必须在break之前，或者移到if块外面确保执行
2. **临时标记清理**: 所有临时标记（同位素、mapnum等）必须有明确的清理时机
3. **兜底机制**: 在渲染层也添加清理，防止上游遗漏
4. **分阶段验证**: 08生成后立即检查同位素，不要等到09才发现
5. **小规模测试**: 修改关键函数后先用1-2个分子测试，再大规模运行

---

## 📝 后续优化建议

1. **单元测试**: 为`annotate_product_sites_with_mapnum()`添加测试，验证：
   - 正确标记:777
   - 同位素被清理
   - 只标记一个位点

2. **断言检查**: 在08中添加断言，产物生成后检查无残留同位素

3. **日志输出**: 记录清理了多少个同位素，便于调试

4. **文档完善**: 在代码注释中明确标注临时标记的生命周期
