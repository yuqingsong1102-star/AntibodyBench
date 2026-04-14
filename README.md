# AntibodyBench

AntibodyBench 用于在**统一靶点集**和**统一后处理评估**下，比较不同 AI 模型在 **VHH de novo design** 任务中的**质量—效率权衡（quality-efficiency trade-off）**。

当前项目支持的核心流程包括：

1. 数据准备（data preparation）
2. 模型运行（model running）
3. 候选收集（candidate collection）
4. benchmark 构建（default benchmark / common-K benchmark）
5. 统一 external evaluation
6. 结果分析与可视化（analysis）

---

## 1. Project overview

本项目的研究目标不是简单比较“哪个模型绝对更强”，而是：

> 在统一靶点集和统一后处理评估下，比较不同模型在真实使用场景中的质量、效率与鲁棒性表现。

当前支持的模型包括：

- BoltzGen
- mBER-open
- Germinal
- RFantibody

---

## 2. Benchmark design

本项目包含两层 benchmark：

### Benchmark 1: default-setting benchmark

各模型按**默认参数**运行，不强制统一生成候选数量。  
这一层主要反映模型在真实使用场景中的表现，关注：

- 总 GPU 时间
- 总候选数
- hit 数
- time-to-first-hit
- GPU hours per hit

对应输出：

- `outputs/benchmark/benchmark_1_runtime.csv`
- `outputs/benchmark/benchmark_1_hits.csv`

### Benchmark 2: common-K benchmark

从每个模型在每个靶点上的输出中，按该模型内部排序截取相同数量的候选（例如 top-10），再进入统一 external evaluation。  
这一层用于增强横向可比性，关注：

- iPTM
- pDockQ2
- BSA
- clash
- sequence developability

对应输出：

- `outputs/benchmark/benchmark_2_k{K}.csv`

> 注意：模型内部得分只用于该模型内部候选截取，不用于跨模型直接比较。  
> 跨模型优劣判断以统一 external evaluation 为准。

---

## 3. Project structure

```text
AntibodyBench/
├── configs/                    # 配置文件
│   ├── dataset.yaml
│   ├── benchmark.yaml
│   ├── evaluation.yaml
│   └── models/
│       ├── boltzgen.yaml
│       ├── mber_open.yaml
│       ├── germinal.yaml
│       └── rfantibody.yaml
│
├── antibench/                  # 新的统一入口和核心逻辑
│   ├── __init__.py
│   ├── __main__.py
│   ├── cli.py
│   ├── dataset.py
│   ├── candidate.py
│   ├── runner.py
│   ├── collection.py
│   ├── benchmark.py
│   └── utils.py
│
├── evaluation/                 # 统一外部评估流水线
│   ├── run_evaluation.py
│   ├── schema.py
│   ├── aggregate.py
│   ├── structure_inference.py
│   └── steps/
│
├── analysis/
│   └── analyze_results.py
│
├── scripts/                    # 旧脚本层（当前仍作为兼容执行层保留）
│   ├── run.sh
│   ├── data_prep/
│   └── native_runners/
│
├── data/                       # 数据集、原始结构、准备后数据
├── native_inputs/              # 各模型输入文件
└── outputs/                    # 所有运行产物

模型运行结果必须统一放在：
outputs/runs/<model>/<target_id>/

当前旧的 pilot_* 目录仅作为历史结果目录保留。
如需接入新 pipeline，应通过软链接或迁移到 outputs/runs/ 下。

当前 evaluation 为多环境运行：Step 1 使用 chai 环境，Step 2–4 使用主环境。

---

## 4. 快速上手（Quick Start）

> 所有命令均在项目根目录 `AntibodyBench/` 下执行。  
> 主 Python 环境为 `/home/yqsong/enter/bin/python`（base conda）。

### 4.1 前置检查

```bash
# 确认依赖
/home/yqsong/enter/bin/python -c "import yaml, biopython, matplotlib; print('OK')"

# 确认 chai conda 环境（用于 Step 1）
conda run -n chai python -c "import chai_lab; print('chai OK')"
```

### 4.2 一键完整流程

```
模型输出 → collect → benchmark_2 → evaluate → analyze
```

**Step 1：收集候选**（扫描 `outputs/runs/<model>/<target_id>/`）

```bash
python -m antibench collect
# 输出：outputs/candidates/candidates_manifest.csv
```

**Step 2：构建 common-K benchmark**

```bash
python -m antibench benchmark --part b2
# 输出：outputs/benchmark/benchmark_2_k10.csv（k 由 configs/benchmark.yaml 控制）
```

**Step 3：Evaluation**（多步骤，自动多环境调度）

```bash
python -m antibench evaluate
# Step 1 (Chai-1) 自动用 chai conda 环境
# Step 2–4 自动用主环境
# 输出：outputs/evaluation/eval_final.csv
```

**Step 4：分析与可视化**

```bash
python -m antibench analyze
# 输出：outputs/analysis/
```

---

### 4.3 按需运行（部分流程）

只跑 evaluation Step 2–4（跳过 Chai-1，适合断点续跑）：

```bash
python -m antibench evaluate --steps 2 3 4
```

只跑 benchmark 1 运行时统计（不需要 eval 结果）：

```bash
python -m antibench benchmark --part b1_runtime
```

跑完 eval 后生成命中统计：

```bash
python -m antibench benchmark --part b1_hits
```

指定自定义候选 CSV 和输出目录：

```bash
python -m antibench evaluate \
  --candidates-csv outputs/benchmark/benchmark_2_k10.csv \
  --out-dir outputs/evaluation_v2
```

多目录对比分析：

```bash
python -m antibench analyze \
  --eval-dirs outputs/evaluation_v1 outputs/evaluation_v2 \
  --model-labels run_v1 run_v2 \
  --rank-by chai1_iptm
```

---

### 4.4 关键配置文件

| 文件 | 用途 |
|---|---|
| `configs/dataset.yaml` | 靶点集路径、数据集过滤条件 |
| `configs/benchmark.yaml` | k_common、hit 阈值、模型列表 |
| `configs/evaluation.yaml` | 评估步骤、设备、Python 环境路径 |
| `configs/models/*.yaml` | 各模型运行参数、输出文件模式、得分字段 |

修改后**无需重启**，所有命令实时读取配置。

---

### 4.5 输出目录结构

```
outputs/
├── runs/                          # 模型运行原始结果（必须）
│   └── <model>/<target_id>/
├── candidates/
│   └── candidates_manifest.csv   # collect 产物
├── benchmark/
│   ├── benchmark_1_runtime.csv
│   ├── benchmark_1_hits.csv
│   └── benchmark_2_k10.csv
├── evaluation/
│   ├── chai1_runs/                # Step 1 结构预测中间结果
│   └── eval_final.csv            # 全量评估结果
└── analysis/                     # 图表和汇总表格
```

---

### 4.6 常见问题

**Q: `collect` 找到 0 个候选？**  
检查模型输出是否放在 `outputs/runs/<model>/<target_id>/`。  
旧的 `outputs/pilot_*` 目录需迁移或用软链接指向新路径：
```bash
ln -s ../../pilot_germinal/germinal/8dtn_A_B outputs/runs/germinal/8dtn_A_B
```

**Q: `evaluate` Step 1 报 `ModuleNotFoundError: chai_lab`？**  
`configs/evaluation.yaml` 中确认 `step1_python: "conda run -n chai python"`，且 `chai` conda 环境已安装 `chai_lab`。

**Q: `benchmark --part b1_hits` 输出全 0 命中？**  
说明 `eval_final.csv` 中 `chai1_iptm` / `pdockq2` 字段为空，需先成功完成 Step 1 evaluation（Chai-1 推断）。

**Q: 想新增一个模型？**  
1. 在 `configs/models/` 新建 `<model>.yaml`，填写 `output_patterns` 和 `internal_score_field`  
2. 在 `configs/benchmark.yaml` 的 `models:` 列表中加入模型名  
3. 将模型输出放到 `outputs/runs/<model>/`，然后正常跑 collect 即可
