# Pipeline恢复指南

## 当前状态（2025-12-29）

### 已完成
✅ **Chunk 0**: 44分钟, 274MB, 已保存
✅ **Chunk 1**: 73分钟, 190MB, 2.5M products, 已保存

### 失败（需重试）
❌ **Chunk 2**: timeout（4小时不够），但产生了671MB部分数据

### 待处理
⏸️ **Chunks 3-13**: 11个chunks待处理

### 进度
**2/14 完成 = 14.3%**

---

## 关机前状态确认

✅ 无进程运行（pipeline已安全停止）
✅ 状态已保存到：`data/output/transforms/polyphenol-2X_BATCHED/pipeline_state.json`
✅ 已完成chunks的数据已保存

---

## 重启方法（开机后）

### 步骤1：恢复pipeline（自动跳过已完成的chunks）

```bash
cd E:\Projects\halogenator
python batch_transform_pipeline.py
```

Pipeline会自动：
- ✅ 跳过chunk 0和1（已完成）
- 🔄 重试chunk 2（从头开始，因为没有checkpoint）
- ▶️ 继续处理chunks 3-13

### 步骤2：监控进度

```bash
# 查看整体进度
cat data/output/transforms/polyphenol-2X_BATCHED/pipeline_state.json

# 查看当前运行的chunk
ps aux | grep "08_transform_library_v2.py"

# 查看最新chunk日志
tail -f data/output/transforms/polyphenol-2X_BATCHED/chunks/chunk_*/transform.log
```

---

## 优化建议（重启前考虑）

### 当前配置
```python
workers = 16
max_in_flight = 6
batch_size = 50000
timeout = 14400  # 4小时
chunk_size = 1000000  # 1M rows
```

### 建议修改（解决timeout问题）

**选项A：增加timeout到8-10小时**
```python
# 在batch_transform_pipeline.py第157行
timeout=28800  # 8小时
```

**选项B：减小chunk大小到500K**
```python
# 在batch_transform_pipeline.py第347行或运行时参数
rows_per_chunk=500000  # 500K rows
```

推荐：**选项A**（增加timeout），因为：
- 不改变已完成chunks的结构
- 更简单，只改一个数字
- 能处理复杂chunk

---

## 文件位置

**主脚本：** `E:\Projects\halogenator\batch_transform_pipeline.py`
**状态文件：** `data/output/transforms/polyphenol-2X_BATCHED/pipeline_state.json`
**输出目录：** `data/output/transforms/polyphenol-2X_BATCHED/chunks/`

**已完成数据：**
- Chunk 0: `chunks/chunk_000_output/products.parquet` (274MB)
- Chunk 1: `chunks/chunk_001_output/products.parquet` (190MB)

---

## 预计完成时间

假设每个chunk平均需要2-4小时：
- 剩余12个chunks × 3小时平均 = **36小时**
- 如果优化后每个2小时 = **24小时**

建议让它在后台运行过夜。

---

## 注意事项

⚠️ **重启前检查：**
1. 确认没有旧的python进程残留（任务管理器）
2. 考虑是否增加timeout（避免chunk 2再次失败）
3. 确保磁盘空间充足（预计最终需要~10-15GB）

📝 **Pipeline特性：**
- 自动断点续传（跳过已完成chunks）
- 失败的chunks会自动重试
- 每个chunk独立运行，互不影响

---

生成时间：2025-12-29
当前进度：2/14 chunks (14.3%)
