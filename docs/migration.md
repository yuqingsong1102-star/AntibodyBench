
# Migration Notes

本文档记录 AntibodyBench 从旧脚本结构迁移到新 MVP 架构后的入口变化与兼容状态。

---

## 1. Why migration

旧版本项目存在以下问题：

- 数据准备、模型运行、候选收集、评估逻辑分散在多个脚本中
- 路径和参数存在较多硬编码
- benchmark 逻辑与 evaluation 逻辑边界不清晰
- 新增模型或修改评估流程的维护成本较高

因此项目重构为：

- `configs/`：统一配置层
- `antibench/`：统一入口和核心逻辑层
- `evaluation/`：保留原统一评估流水线
- `analysis/`：保留结果分析脚本
- `scripts/`：暂时保留为 legacy bridge

---

## 2. New entry points

当前推荐的主入口：

```bash
python -m antibench collect
python -m antibench benchmark --part runtime
python -m antibench benchmark --part b2 --k 10
python -m antibench evaluate
python -m antibench benchmark --part hits
python -m antibench analyze