# 糖环筛选脚本使用指南

## 概述

`filter_sugar_rings.py` 是一个用于筛选分子中糖环数量的工具脚本，基于halogenator项目中的`sugar_mask`模块开发。该脚本可以读取SDF、CSV或XLSX格式的文件，过滤掉糖环数量超过指定阈值的分子。

## 功能特性

- **多格式支持**: SDF、CSV、XLSX文件格式
- **自动糖环检测**: 基于halogenator的启发式糖环识别算法
- **灵活的阈值设置**: 可自定义最大糖环数量
- **详细的输出**: 为每个分子添加糖环计数信息
- **完善的日志**: 提供详细的处理过程和统计信息

## 糖环识别标准

脚本使用以下标准识别糖环（基于halogenator的sugar_mask模块）：

1. **环大小**: 5或6元环
2. **芳香性**: 非芳香环
3. **环氧原子**: 正好包含1个环氧原子
4. **羰基**: 环内无羰基（C=O，其中C在环内）
5. **SP3碳比例**:
   - 5元环: ≥40% SP3碳
   - 6元环: ≥50% SP3碳
6. **外环氧取代基**: ≥2个单键连接的外环氧，或存在C-糖苷证据
7. **评分系统**: 使用证据评分系统（阈值≥8.0）综合判断

## 安装要求

```bash
# 必需的依赖
pip install rdkit pandas openpyxl
```

## 使用方法

### 基本用法

```bash
# SDF文件筛选
python scripts/filter_sugar_rings.py input.sdf -o output.sdf --max-rings 3

# CSV文件筛选
python scripts/filter_sugar_rings.py input.csv -o output.csv --max-rings 3

# XLSX文件筛选
python scripts/filter_sugar_rings.py input.xlsx -o output.xlsx --max-rings 3
```

### 命令行参数

```
必需参数:
  input                 输入文件路径 (SDF/CSV/XLSX)
  -o, --output         输出文件路径

可选参数:
  --max-rings N        最大糖环数量阈值 (默认: 3)
  --smiles-column COL  SMILES列名，用于CSV/XLSX (默认: SMILES)
  --sheet SHEET        工作表名称，仅用于XLSX (默认: 第一个工作表)
  --verbose, -v        显示详细日志
  -h, --help           显示帮助信息
```

### 使用示例

```bash
# 示例1: 过滤SDF文件，只保留糖环数≤2的分子
python scripts/filter_sugar_rings.py compounds.sdf -o filtered.sdf --max-rings 2

# 示例2: 处理CSV文件，SMILES列名为"smiles"
python scripts/filter_sugar_rings.py data.csv -o filtered.csv --smiles-column smiles

# 示例3: 处理XLSX文件的特定工作表，显示详细日志
python scripts/filter_sugar_rings.py data.xlsx -o filtered.xlsx --sheet "Sheet2" --verbose

# 示例4: 默认阈值(3个糖环)
python scripts/filter_sugar_rings.py input.csv -o output.csv
```

## 测试验证

### 测试数据集

脚本已经过全面测试，测试数据包括：

| 分子名称 | SMILES | 糖环数 |
|---------|--------|--------|
| Benzene | C1=CC=CC=C1 | 0 |
| Glucose | OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O | 1 |
| Aspirin | CC(=O)Oc1ccccc1C(=O)O | 0 |
| Quercetin-3-O-glucoside | O=C1C=C(O[C@@H]2...)c2O1 | 1 |
| Lactose | O[C@H]1[C@H](O)[C@@H](O)... | 2 |
| Quercetin-diglucoside | O=C1C=C(O[C@@H]2O[C@H]...)c2O1 | 2 |

### 测试结果

#### 测试1: max-rings=3（默认值）
```
输入: 10个分子
通过: 10个 (100%)
过滤: 0个 (0%)
```
所有分子都通过（因为最大糖环数为2，≤3）

#### 测试2: max-rings=1
```
输入: 10个分子
通过: 8个 (80%)
过滤: 2个 (20%)
```
正确过滤掉了Lactose和Quercetin-diglucoside（各2个糖环）

#### 测试3: 多格式验证
- ✅ SDF格式: 通过 (8/10 分子)
- ✅ CSV格式: 通过 (8/10 分子)
- ✅ XLSX格式: 通过 (8/10 分子)

所有格式的筛选结果一致

## 输出说明

### CSV/XLSX输出

输出文件会在原有列的基础上添加`sugar_ring_count`列，记录每个分子的糖环数量：

```csv
SMILES,Name,sugar_ring_count
C1=CC=CC=C1,Benzene,0
OC[C@H]1O[C@H](O)...,Glucose,1
O[C@H]1[C@H](O)...,Lactose,2
```

### SDF输出

SDF文件会为每个分子添加`sugar_ring_count`属性：

```
> <sugar_ring_count>
1
```

### 终端输出

脚本运行时会显示处理摘要：

```
============================================================
筛选摘要:
  输入文件: test_sugar_filter_input.csv
  输出文件: test_sugar_filter_output.csv
  最大糖环数: 1
  总分子数: 10
  通过筛选: 8 (80.0%)
  被过滤: 2 (20.0%)
============================================================
```

## 常见问题

### Q1: 如何查看详细的处理日志？
使用`--verbose`或`-v`参数：
```bash
python scripts/filter_sugar_rings.py input.csv -o output.csv -v
```

### Q2: CSV文件的SMILES列名不是"SMILES"怎么办？
使用`--smiles-column`参数指定列名：
```bash
python scripts/filter_sugar_rings.py input.csv -o output.csv --smiles-column "molecule_smiles"
```

### Q3: XLSX文件有多个工作表，如何指定？
使用`--sheet`参数：
```bash
python scripts/filter_sugar_rings.py input.xlsx -o output.xlsx --sheet "Data"
```

### Q4: 什么是"糖环"？
糖环是指符合糖类结构特征的环状结构，通常包含：
- 5或6元环
- 1个环氧原子
- 多个羟基（-OH）取代基
- 较高的SP3碳比例
- 非芳香性

典型的糖环包括葡萄糖环、半乳糖环、核糖环等。

### Q5: 为什么myo-Inositol不被识别为糖环？
尽管myo-Inositol是环己醇，但它不包含环氧原子，因此不符合糖环的定义（必须有环氧）。

### Q6: 脚本处理大文件需要多长时间？
处理速度取决于：
- 文件大小
- 分子复杂度
- 硬件性能

典型速度: ~100-1000 分子/秒

## 技术细节

### 糖环检测算法

脚本使用halogenator的`sugar_mask._find_sugar_rings()`函数，该函数采用证据评分系统：

**评分因子**:
- 环大小 (5元环: +1.0, 6元环: +2.0)
- 正好1个环氧: +3.0
- 无内部羰基: +2.0
- SP3比例评分: 最高+3.0
- 外环氧数量: 每个+1.0 (最高+4.0)
- C-糖苷证据: +2.0

**接受阈值**: 总分≥8.0 且有强证据（≥2个外环氧或C-糖苷证据）

### 依赖的模块

- `halogenator.sugar_mask._find_sugar_rings()`: 糖环检测核心函数
- `halogenator.sugar_mask._get_default_sugar_cfg()`: 默认配置

## 脚本位置

```
E:\Projects\halogenator\scripts\filter_sugar_rings.py
```

## 测试文件

测试输入文件:
```
E:\Projects\halogenator\test_sugar_filter_input.csv
E:\Projects\halogenator\test_sugar_filter_input.sdf
E:\Projects\halogenator\test_sugar_filter_input.xlsx
```

测试输出文件:
```
E:\Projects\halogenator\test_sugar_filter_output.csv
E:\Projects\halogenator\test_sugar_filter_output.sdf
E:\Projects\halogenator\test_sugar_filter_output.xlsx
E:\Projects\halogenator\test_sugar_filter_output_max1.csv
```

## 已知限制

1. **编码问题**: 在某些Windows终端中，中文日志可能显示为乱码（不影响功能）
2. **立体化学**: 脚本依赖SMILES中的立体化学信息，无立体化学标记可能影响准确性
3. **特殊结构**: 某些非常规糖结构可能无法被正确识别

## 更新日志

### v1.0 (2025-12-08)
- ✅ 初始版本发布
- ✅ 支持SDF、CSV、XLSX格式
- ✅ 基于halogenator的sugar_mask模块
- ✅ 完整的测试覆盖
- ✅ 详细的日志和统计信息

## 作者与许可

基于halogenator项目的sugar_mask模块开发。

## 联系方式

如有问题或建议，请在halogenator项目中提issue。
