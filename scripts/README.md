# Scripts

## 主入口

```bash
bash scripts/run.sh --model <boltzgen|mber-open|germinal|RFantibody> [options]
```

| 参数 | 说明 | 默认值 |
|---|---|---|
| `--model` | 模型名（必填） | — |
| `--sample-list` | 样本列表文件 | — |
| `--sample-id` | 仅运行指定样本 | — |
| `--input-root` | 原生输入根目录 | `native_inputs` |
| `--out-root` | 输出根目录 | `outputs/native_predictions` |
| `--max-samples` | 最多处理样本数（0=全部） | `1` |
| `--continue-on-error` | 出错是否继续 | `1` |

## 数据准备

```bash
python scripts/data_prep/prepare_native_inputs.py
```

只需运行一次。生成 `native_inputs/<model>/<sample_id>/` 目录树。

## 模型运行器

| 文件 | 模型 |
|---|---|
| `native_runners/boltzgen.sh` | boltzgen |
| `native_runners/mber_open.sh` | mBER-open |
| `native_runners/germinal.sh` | germinal |
| `native_runners/rfantibody.sh` | RFantibody |
| `native_runners/collect_candidates.py` | 候选收集器（通用） |

## 输出结构

每个样本的输出：

```
outputs/<run_name>/<model>/<sample_id>/
├── candidate_manifest.csv
├── run_meta.json
├── native_run/
├── run_stdout.log
└── run_stderr.log
```

模型级汇总：

- outputs/native_predictions/<model>/manifest.csv

核心字段包括：

- candidate_id
- sample_id
- model
- candidate_rank
- candidate_name
- status
- sequence_path
- structure_path
- meta_path
- duration_sec
- error_summary

说明：

- 如果模型原始输出只有 CSV 序列列，collector 会生成标准化 FASTA sidecar
- structure_path 是可选字段，没有结构也允许保留候选
- manifest.csv 每次运行都会覆盖重写，不会向旧结果追加

- `scripts/data_prep/build_model_inputs_native.py`
	- 直接读取 `data/prepared/dataset_index_ready.csv` 与 `data/prepared/epitopes/<sample_id>.json`，优先使用 `design_structure_path`，一次生成四个模型的 native 输入。
	- 默认输出到 `native_inputs/<model>/<sample_id>/`，并刷新 `native_inputs/_native_manifest.csv`。

## Ops

- `scripts/build_apptainer_images.sh`

## Legacy / utility

- `scripts/legacy/collect_top1.py`

## Evaluation

一键评估（Phase1，见 `evaluation/README.md`）：

```bash
python -m evaluation.pipeline.run_eval_pipeline --fallback-to-source-structure
```

默认产物目录：`outputs/evaluation/external_binding_benchmark/`（含 `external_eval_*.csv`、`report.md`、`figures/`、`judge_inputs/`）。

