# AntibodyBench Project Layout

This document is a practical map for daily work.

## Top-Level Roles

- `models/`（仓库外）: 四个原生模型代码仓（RFantibody / germinal / BindCraft / boltzgen）
- `data/raw/`: 原始源数据，包含源索引、目标结构与 full complex
- `data/prepared/`: 派生基准数据，包含 ready 索引与 epitope 缓存
- `native_inputs/`: 四个模型各自可直接运行的原生输入树
- `outputs/`: 运行与评估产物
- `evaluation/`: ingest / judge / metrics / report / pipeline
- `scripts/`: 入口脚本与 runner
- `docs/`: 快速地图、问题记录、执行清单

## Recommended Daily Entry Points

Use these entry points first:

- `scripts/data_prep/prepare_native_inputs.py`
- `scripts/run.sh`
- `evaluation/pipeline/run_eval_pipeline.py`
- `scripts/build_apptainer_images.sh`

Notes:

- `scripts/data_prep/*` is the canonical data-prep namespace.
- `prepare_native_inputs.py` is the only recommended day-to-day data-prep entry.
- Evaluation commands should use `evaluation/*` paths.

See:

- `scripts/README.md`
- `data/prepared/README.md`
- `outputs/native_predictions_run2/`
- `outputs/evaluation/external_binding_benchmark_phase1/`

## Key Data Files

- `data/raw/dataset_index.csv`: 原始源索引
- `data/prepared/dataset_index_ready.csv`: 过滤后的 ready 子集
- `native_inputs/_native_manifest.csv`: 当前 native 输入清单

## Keep vs Ignore

- Keep: `data/raw/`, `data/prepared/`, `native_inputs/`, `scripts/`, `evaluation/`, `docs/`
- Generated only: `outputs/`
- Optional/local-only: local dependency folders, caches, pycache files

当前常用输出入口：`outputs/native_predictions_run2/` 与 `outputs/evaluation/external_binding_benchmark_phase1/`。

