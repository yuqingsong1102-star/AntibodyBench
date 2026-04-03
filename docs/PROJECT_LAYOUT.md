# AntibodyBench Project Layout

This document is a practical map for daily work.

## Top-Level Roles

- `models/` (仓库外): 各原生模型代码仓（RFantibody/germinal/BindCraft/boltzgen）
- `data/`: raw/reference data (do not edit casually)
- `inputs/`: benchmark index and model input metadata
- `outputs/`: generated artifacts (input/prediction/evaluation)
- `evaluation/`: metrics implementation and reporting pipeline
- `scripts/`: runnable entry scripts
- `docs/`: quick project maps and operation notes

## Recommended Daily Entry Points

Use these entry points first:

- `scripts/run.sh`
- `scripts/data_prep/generate_dataset_index.py`
- `scripts/data_prep/fetch_reference_complexes.py`
- `scripts/data_prep/fill_cdr_h3_from_anarci.py`
- `scripts/build_apptainer_images.sh`
- `evaluation/pipeline/run_eval_pipeline.py`

Notes:

- `scripts/data_prep/*` is the canonical data-prep entry namespace.
- Data-prep commands should always use `scripts/data_prep/*` paths.
- Evaluation commands should use `evaluation/*` paths.

See:

- `scripts/README.md`
- `inputs/README.md`
- `outputs/native_predictions_run2/`
- `outputs/evaluation/all_models/`

## Key CSV Files

- `inputs/dataset_index.csv`: canonical index
- `inputs/dataset_index_h3_annotated.csv`: H3 annotation output
- `inputs/dataset_index_ready.csv`: filtered ready subset

## Keep vs Ignore

- Keep: `inputs/`, `scripts/`, `evaluation/`, `data/raw/reference_complexes/`, `docs/`
- Generated only: `outputs/`
- Optional/local-only: local dependency folders, caches, pycache files

当前常用输出入口：`outputs/native_predictions_run2/` 与 `outputs/evaluation/all_models/`。

