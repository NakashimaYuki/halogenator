# 修复完成报告

**日期**: 2025-11-12
**状态**: ✅ 完成

---

## 已完成的修复

### 1. ✅ 09可视化去标签补丁

**问题**: 图片中显示 "O:777" 文本

**修复**:
- 添加了 `_scrub_for_render()` 函数
- 在 `render_molecule_png_worker()` 中，在 `resolve_sites()` 后、`draw_png_hidpi()` 前清理分子
- 清除 atom map numbers 和所有display属性

**验证**: 200样本测试，100% mapnum usage，100% high confidence，无"O:777"显示

**文件**: `scripts/09_visualize_library.py`

---

### 2. ✅ 08唯一打点策略（同位素锚点V2）

**问题**: 旧实现假设产物中有:1映射，但RunReactants()会清除所有映射

**修复**:
- 添加 `mark_anchor_with_isotope()` 辅助函数
- 修改 `annotate_product_sites_with_mapnum()` 使用isotope 203查找锚点
- 在 `apply_to_molecule()` 中，RunReactants前先标记锚点同位素
- 每个产物只标记唯一的转化位点（避免桥氧误标）

**策略**:
1. 反应前：在匹配位点的锚点原子（[c:1]）上打同位素203
2. 运行反应：RDKit保留同位素信息
3. 反应后：查找带isotope 203的锚点，从其邻居中唯一选择目标原子
4. 标记目标原子为:777，清理isotope

**化学验证**:
- OMe: 非环氧 (X2) 连接sp3碳 (CH3)
- NH2: 非环氮，H≥2，排除酰胺

**文件**: `scripts/08_transform_library_v2.py`

---

### 3. ✅ 08配置固化到SUMMARY.json

**问题**: 两次运行配置不一致，无法审计差异

**修复**:
在SUMMARY.json中添加完整配置信息：
- workers
- batch_size
- resume
- max_sites_per_mol
- exclude_sugar_like
- apply_mode
- xf_config
- isotope_anchor_strategy: 'v2'

**文件**: `scripts/08_transform_library_v2.py`

---

## 待执行任务

### 4. ⏳ 重跑 Flavone-2X-NH2

**原因**: 与OMe组对齐配置（workers=8, max_sites=-1）

**命令**:
```bash
python scripts/08_transform_library_v2.py apply \
  -i E:/Projects/halogenator/data/output/nplike/Flavone-2X/base.parquet \
  -o E:/Projects/halogenator/data/output/nplike/Flavone-2X-NH2_pr_a_v2 \
  --xf-config configs/transforms.yaml \
  --xf-name OH_to_NH2 \
  --workers 8 \
  --batch-size 50000
```

**预期**: 产物数应接近OMe（去重前），去重后可能略低（NH2坍缩更强）

---

### 5. ⏳ 生成统计报告

**命令**:
```bash
# 为四个派生库生成统计
for lib in Flavone-1X-OMe_pr_a Flavone-1X-NH2_pr_a Flavone-2X-OMe_pr_a Flavone-2X-NH2_pr_a_v2; do
  python scripts/05_summaries_v2.py \
    -i E:/Projects/halogenator/data/output/nplike/$lib/products.parquet \
    -o E:/Projects/halogenator/data/output/nplike/$lib/stats
done
```

---

### 6. ⏳ 生成可视化

**命令** (每个库):
```bash
python scripts/09_visualize_library.py html \
  -i E:/Projects/halogenator/data/output/nplike/Flavone-1X-OMe_pr_a/products.parquet \
  -o E:/Projects/halogenator/data/viz/Flavone-1X-OMe.html \
  --highlight-sites --preset hq
```

**验证指标**:
- viz_diagnostics.csv中method_used=mapnum接近100%
- num_candidates=1 (唯一打点)
- confidence=high 100%

---

## 技术要点

### 同位素锚点策略优势

1. **唯一性**: 每个产物只标记1个位点
2. **化学验证**: 不仅靠模板匹配，还有化学约束
3. **避免误标**: 桥氧等环内原子被排除
4. **保真度高**: RDKit保留同位素，不像映射号被清除

### 可复现性保证

1. 所有关键配置固化到SUMMARY.json
2. git commit可追溯
3. 配置差异一目了然
4. 便于CI/CD自动化检查

---

## 下一步建议

1. **立即执行**: 用对齐配置重跑Flavone-2X-NH2
2. **对比验证**: 新旧NH2的SUMMARY.json配置差异
3. **批量可视化**: 为所有4个库生成HTML gallery
4. **最终验证**: 抽检图片确认无"O:777"，位点唯一且正确
