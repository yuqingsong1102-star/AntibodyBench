# AntibodyBench

用于“从头抗体模型评测”的统一框架，支持多模型接入、自动化执行与结果分析。

## 目录结构

```text
AntibodyBench/
├── data/                        # raw 源数据与 prepared 派生基准数据
├── native_inputs/              # 四个模型各自的原生输入树
├── outputs/                     # native 推理结果与评估结果（生成物）
├── scripts/
│   ├── run.sh                  # 原生推理编排入口
│   ├── data_prep/              # 数据准备脚本
│   ├── native_runners/         # 四个模型的原生 runner
│   └── tools/                  # 可选工具（如 TMscore / DockQ 包装，评估流水线当前不调用）
├── evaluation/                  # ingest / judge / metrics / report / pipeline
└── docs/                        # 项目地图、问题记录、后续计划
```

说明：四个模型代码仓位于仓库外的 `models/` 目录（RFantibody / germinal / BindCraft / boltzgen），`AntibodyBench/` 只负责统一数据、运行和评估。

## 快速开始

```bash
cd AntibodyBench
bash scripts/run.sh --model boltzgen --sample-id 8q3j_B_A --out-root outputs/native_predictions_run2
```

## 第一步：准备数据索引与原生输入

```bash
python scripts/data_prep/prepare_native_inputs.py
```

现在主流程只强调两个索引概念：

- `data/raw/dataset_index.csv`：原始源索引
- `data/prepared/dataset_index_ready.csv`：最终用于构建四模型输入、运行与评估的 ready 子集

派生产物现在统一落在两个位置：`data/prepared/` 放 ready 索引与 epitope 缓存，`native_inputs/` 放四模型原生输入树。

## 第二步：运行原生推理（默认 smoke_1）

```bash
# 默认只跑 1 个样本（--max-samples 默认值为 1）
bash scripts/run.sh --model RFantibody
bash scripts/run.sh --model germinal
bash scripts/run.sh --model BindCraft
bash scripts/run.sh --model boltzgen
```

常用参数：

- `--sample-id <id>`：只跑一个指定样本
- `--max-samples <N>`：限制样本数（`0` 表示不限制）
- `--continue-on-error <0|1>`：失败是否继续（默认 `1`）
- `--input-root` / `--out-root`：覆盖输入输出根目录

## 输出契约（统一）

每样本目录：

`outputs/native_predictions_<run_tag>/<model>/<sample_id>/`

最小交付文件：

- `candidate_manifest.csv`
- `run_meta.json`
- `native_run/`
- `run_stdout.log`
- `run_stderr.log`

模型级汇总：

- `outputs/native_predictions_<run_tag>/<model>/manifest.csv`

## 第三步：运行统一评估（Phase1）

在仓库根目录 `AntibodyBench/` 下执行：

```bash
python -m evaluation.pipeline.run_eval_pipeline \
  --native-root outputs/native_predictions_run2 \
  --index-csv data/prepared/dataset_index_ready.csv \
  --out-dir outputs/evaluation/external_binding_benchmark_phase1 \
  --fallback-to-source-structure
```

无外部 judge manifest 时需加 `--fallback-to-source-structure`（用模型 `top1` 结构占位跑通评分，非正式 judge 结论）。默认会优先使用 `outputs/native_predictions_run2`，否则回退到 `native_predictions_real` 与 `native_predictions_smoke`。

详见 `evaluation/README.md`。
