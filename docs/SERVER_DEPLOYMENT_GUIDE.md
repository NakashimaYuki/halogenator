# Halogenator 服务器部署指南

本指南帮助你将 Halogenator 项目从 Windows 部署到 Linux 服务器。

---

## 概述

由于输入数据文件（约2.5GB）不在Git版本控制中，需要分两步部署：
1. **Git clone** - 获取代码
2. **数据传输** - 传输输入数据文件

---

## 第一步：Windows 端准备

### 1.1 打包数据

在 Windows 的 PowerShell 或 Git Bash 中运行：

```bash
cd E:\Projects\halogenator
python scripts/pack_data_for_server.py
```

这会生成一个压缩包：`nplike_v2_data_YYYYMMDD.tar.gz`（约 1-2GB）

### 1.2 确保代码已推送到 GitHub

```bash
cd E:\Projects\halogenator

# 查看状态
git status

# 如果有需要提交的更改
git add configs/transforms.yaml
git commit -m "Add transforms config for server deployment"
git push
```

---

## 第二步：Linux 服务器端

### 2.1 克隆代码

```bash
# SSH 到服务器
ssh user@your-server

# 进入工作目录
cd /path/to/projects

# 克隆仓库
git clone https://github.com/YOUR_USERNAME/halogenator.git
cd halogenator

# 切换到正确的分支
git checkout fix/pr2-contract-and-sugar-accept
```

### 2.2 创建 Python 环境

```bash
# 创建虚拟环境
python3 -m venv venv
source venv/bin/activate

# 安装依赖
pip install -e .

# 或者如果有 requirements.txt
pip install -r requirements.txt
```

### 2.3 创建数据目录结构

```bash
mkdir -p data/output
```

---

## 第三步：传输数据

### 方法A：使用 scp（简单直接）

在 **Windows** 端运行（Git Bash 或 PowerShell）：

```bash
# 替换 user@server 和路径
scp nplike_v2_data_20260203.tar.gz user@server:/path/to/halogenator/
```

### 方法B：使用 rsync（支持断点续传）

```bash
# 替换 user@server 和路径
rsync -avz --progress nplike_v2_data_20260203.tar.gz user@server:/path/to/halogenator/
```

---

## 第四步：服务器端解压数据

SSH 到服务器后：

```bash
cd /path/to/halogenator

# 解压数据
tar -xzvf nplike_v2_data_20260203.tar.gz

# 验证数据
ls -la data/output/nplike_v2/

# 清理压缩包（可选）
rm nplike_v2_data_20260203.tar.gz
```

---

## 第五步：验证部署

```bash
# 激活环境
source venv/bin/activate

# 验证数据文件存在
python -c "
import pyarrow.parquet as pq
from pathlib import Path

nplike = Path('data/output/nplike_v2')
for p in sorted(nplike.glob('**/products.parquet')):
    pf = pq.ParquetFile(p)
    print(f'{p.parent.name}: {pf.metadata.num_rows:,} rows')
"

# 验证 transform 配置
python -c "
import yaml
with open('configs/transforms.yaml') as f:
    cfg = yaml.safe_load(f)
print(f'Loaded {len(cfg.get(\"transforms\", []))} transform rules')
"
```

预期输出类似：
```
aa_peptide-1X: 43,456 rows
aa_peptide-2X: 789,123 rows
alkaloid-1X: 234,567 rows
...
Loaded 24 transform rules
```

---

## 第六步：运行 Transform 任务

示例：运行 alkaloid-2X 的变换任务

```bash
# 激活环境
source venv/bin/activate

# 创建输出目录
mkdir -p data/output/transforms

# 运行变换
python scripts/08_transform_library_v2.py apply \
  --input data/output/nplike_v2/alkaloid-2X/products.parquet \
  --outdir data/output/transforms/alkaloid-2X_FG_AMINE_AR__NH2__TO__F \
  --xf-config configs/transforms.yaml \
  --xf-name FG_AMINE_AR__NH2__TO__F \
  --workers 16 \
  --max-in-flight 6 \
  --batch-size 50000 \
  --use-bloom-filter \
  --bloom-expected-items 5000000
```

对于大数据集（>10M分子），参考 `SESSION_HANDOFF_2026_01_22.md` 中的 shuffle-rechunking 策略。

---

## 快速命令汇总

### Windows 端
```bash
# 1. 打包数据
python scripts/pack_data_for_server.py

# 2. 推送代码
git add configs/transforms.yaml
git commit -m "Add transforms config"
git push

# 3. 传输数据
scp nplike_v2_data_*.tar.gz user@server:/path/to/halogenator/
```

### Linux 端
```bash
# 1. 克隆代码
git clone https://github.com/YOUR_USERNAME/halogenator.git
cd halogenator
git checkout fix/pr2-contract-and-sugar-accept

# 2. 设置环境
python3 -m venv venv
source venv/bin/activate
pip install -e .

# 3. 解压数据
mkdir -p data/output
tar -xzvf nplike_v2_data_*.tar.gz

# 4. 验证
ls data/output/nplike_v2/
```

---

## 常见问题

### Q: scp 传输中断怎么办？
使用 rsync，它支持断点续传：
```bash
rsync -avz --progress --partial nplike_v2_data_*.tar.gz user@server:/path/
```

### Q: 解压时磁盘空间不足？
压缩包约 1-2GB，解压后约 2.5GB，确保有至少 5GB 可用空间。

### Q: 如何只传输特定数据集？
可以直接传输单个 parquet 文件：
```bash
scp data/output/nplike_v2/alkaloid-2X/products.parquet user@server:/path/to/halogenator/data/output/nplike_v2/alkaloid-2X/
```

---

## 文件清单

需要在服务器上存在的关键文件：

| 文件 | 来源 | 说明 |
|------|------|------|
| `configs/transforms.yaml` | Git | Transform 规则定义 |
| `configs/halogen_rules_by_class.yaml` | Git | 卤化规则配置 |
| `scripts/08_transform_library_v2.py` | Git | Transform 引擎 |
| `data/output/nplike_v2/*/products.parquet` | 数据包 | 输入数据 |

---

**文档版本:** 1.0
**最后更新:** 2026-02-03
