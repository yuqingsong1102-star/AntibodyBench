# Scripts Index

## Native inference entry

`scripts/run.sh` 已改为**原生推理编排入口**，不再执行 preprocess/postprocess/aggregate 的三阶段评测流程。

### 用法

```bash
bash scripts/run.sh --model <RFantibody|germinal|BindCraft|boltzgen> [options]
```

参数说明：

- `--model`：必填，模型名。
- `--input-root`：原生输入根目录，默认 `data/model_inputs_native`。
- `--out-root`：输出根目录，默认 `outputs/native_predictions`。
- `--max-samples`：最多处理样本数，默认 `1`（即默认 smoke_1）。
- `--sample-id`：仅运行一个样本（会忽略 `--max-samples`）。
- `--continue-on-error`：`0|1`，默认 `1`，即 `sequential_continue`。

### 输入来源

仅从以下目录读取样本输入，不再依赖适配层 `outputs/input`：

`data/model_inputs_native/<model>/<sample_id>/`

### 运行器

模型执行由 `scripts/native_runners/` 下的薄封装脚本负责：

- `rfantibody.sh`
- `germinal.sh`
- `bindcraft.sh`
- `boltzgen.sh`
- `collect_top1.py`（统一 Top-1 抽取）

### 输出契约

每个样本输出目录：

`outputs/native_predictions/<model>/<sample_id>/`

最小交付文件：

- `top1_structure.(pdb|cif)`
- `top1_sequence.fasta`
- `top1_meta.json`
- `run_stdout.log`
- `run_stderr.log`

模型级汇总：

- `outputs/native_predictions/<model>/manifest.csv`

`manifest.csv` 字段：

- `sample_id`
- `status`
- `structure_path`
- `sequence_path`
- `meta_path`
- `duration_sec`
- `error_summary`

### 失败与可追溯

- 单样本失败默认不中断后续样本（`--continue-on-error 1`）。
- 外部命令标准输出/错误分别写入 `run_stdout.log` / `run_stderr.log`。
- 即使运行环境不完整，也会产出 `top1_meta.json` 和 `manifest.csv` 失败记录。

### 最小验证（每个模型 1 个样本）

```bash
SAMPLE_ID=8q3j_B_A
OUT_ROOT=outputs/native_predictions_smoke

bash scripts/run.sh --model RFantibody --sample-id "${SAMPLE_ID}" --out-root "${OUT_ROOT}"
bash scripts/run.sh --model germinal --sample-id "${SAMPLE_ID}" --out-root "${OUT_ROOT}"
bash scripts/run.sh --model BindCraft --sample-id "${SAMPLE_ID}" --out-root "${OUT_ROOT}"
bash scripts/run.sh --model boltzgen --sample-id "${SAMPLE_ID}" --out-root "${OUT_ROOT}"
```

检查这两个文件是否存在：

- `outputs/native_predictions_smoke/<model>/manifest.csv`
- `outputs/native_predictions_smoke/<model>/<sample_id>/top1_meta.json`

### 每个模型怎么跑（环境变量 + 命令模板）

下面 4 个 runner 都支持覆盖工作目录和命令模板：

- RFantibody：`RFANTIBODY_WORKDIR`、`RFANTIBODY_CMD`
- germinal：`GERMINAL_WORKDIR`、`GERMINAL_CONDA_ENV`、`GERMINAL_USE_CONDA_RUN`、`GERMINAL_CMD`
- BindCraft：`BINDCRAFT_WORKDIR`、`BINDCRAFT_CONDA_ENV`、`BINDCRAFT_USE_CONDA_RUN`、`BINDCRAFT_CMD`
- boltzgen：`BOLTZGEN_WORKDIR`、`BOLTZGEN_CONDA_ENV`、`BOLTZGEN_USE_CONDA_RUN`、`BOLTZGEN_CMD`

#### RFantibody

必填：`RFANTIBODY_CMD`（未设置会直接失败并写入失败元数据）。

可用占位符：

- `__SAMPLE_ID__`
- `__INPUT_DIR__`
- `__OUTPUT_DIR__`
- `__CONFIG_FILE__`
- `__TARGET_PDB__`
- `__FRAMEWORK_PDB__`
- `__HOTSPOTS__`
- `__DESIGN_LOOPS__`
- `__NUM_DESIGNS__`
- `__NUM_SEQS__`
- `__NUM_RECYCLES__`

示例：

```bash
export RFANTIBODY_WORKDIR=/path/to/models/RFantibody
export RFANTIBODY_CMD='python run_pipeline.py --target "__TARGET_PDB__" --out "__OUTPUT_DIR__" --sample "__SAMPLE_ID__"'
bash scripts/run.sh --model RFantibody --sample-id 8q3j_B_A
```

#### germinal

默认会执行内置命令；你也可以用 `GERMINAL_CMD` 覆盖。

可用占位符：

- `__SAMPLE_ID__`
- `__INPUT_DIR__`
- `__OUTPUT_DIR__`
- `__TARGET_YAML__`
- `__OVERRIDES_FILE__`
- `__OVERRIDES__`（`run_overrides.txt` 合并为一行后替换）

示例：

```bash
export GERMINAL_WORKDIR=/path/to/models/germinal
export GERMINAL_CONDA_ENV=germinal
export GERMINAL_USE_CONDA_RUN=1
export GERMINAL_CMD='python run_germinal.py --config "__TARGET_YAML__" results_dir="__OUTPUT_DIR__" experiment_name="__SAMPLE_ID__"'
bash scripts/run.sh --model germinal --sample-id 8q3j_B_A
```

#### BindCraft

默认会执行内置命令；你也可以用 `BINDCRAFT_CMD` 覆盖。

可用占位符：

- `__SAMPLE_ID__`
- `__INPUT_DIR__`
- `__OUTPUT_DIR__`
- `__SETTINGS_FILE__`
- `__FILTERS_FILE__`
- `__ADVANCED_FILE__`

示例：

```bash
export BINDCRAFT_WORKDIR=/path/to/models/BindCraft
export BINDCRAFT_CONDA_ENV=BindCraft
export BINDCRAFT_USE_CONDA_RUN=1
export BINDCRAFT_CMD='python bindcraft.py --settings "__SETTINGS_FILE__" --filters "__FILTERS_FILE__" --advanced "__ADVANCED_FILE__"'
bash scripts/run.sh --model BindCraft --sample-id 8q3j_B_A
```

#### boltzgen

默认会执行内置命令；你也可以用 `BOLTZGEN_CMD` 覆盖。

可用占位符：

- `__SAMPLE_ID__`
- `__INPUT_DIR__`
- `__OUTPUT_DIR__`
- `__SPEC_FILE__`
- `__PROTOCOL__`
- `__NUM_DESIGNS__`
- `__BUDGET__`

示例：

```bash
export BOLTZGEN_WORKDIR=/path/to/models/boltzgen
export BOLTZGEN_CONDA_ENV=boltzgen
export BOLTZGEN_USE_CONDA_RUN=1
export BOLTZGEN_CMD='boltzgen run "__SPEC_FILE__" --output "__OUTPUT_DIR__" --protocol "__PROTOCOL__" --num_designs __NUM_DESIGNS__ --budget __BUDGET__'
bash scripts/run.sh --model boltzgen --sample-id 8q3j_B_A
```

### 当前机器实测注意事项

- `germinal`：可以启动，但当前环境会在 IgLM 加载阶段报 `torch<2.6`（需要升级 torch）。
- `RFantibody`：能启动到模型权重加载，但会报 CUDA OOM（需要更大显存/减少并发）。
- `boltzgen`：当前可用环境名是 `bg`，不是 `boltzgen`（可设置 `BOLTZGEN_CONDA_ENV=bg`）。
- `BindCraft`：真实运行会比较久，单样本可能持续很长时间才结束。

## Data preparation

- `scripts/data_prep/generate_dataset_index.py`
- `scripts/data_prep/fetch_reference_complexes.py`
- `scripts/data_prep/fill_cdr_h3_from_anarci.py`
- `scripts/data_prep/build_model_inputs.py`
- `scripts/data_prep/build_model_inputs_native.py`
- `scripts/data_prep/extract_epitopes_from_complexes.py`
- `scripts/data_prep/apply_epitopes_to_native_inputs.py`

## Ops

- `scripts/ops/build_apptainer_images.sh`

## Legacy / utility

- `scripts/prepare_antibody_dataset.py`

## Evaluation (Track B)

评估分析主入口已迁移到 `evaluation/` 目录，建议使用一键命令：

```bash
python -m evaluation.pipeline.run_eval_pipeline
```

默认产物：

- `outputs/evaluation/all_models/evaluation_long.csv`
- `outputs/evaluation/all_models/summary_by_model.csv`
- `outputs/evaluation/all_models/summary_by_sample.csv`
- `outputs/evaluation/all_models/figures/*`
- `outputs/evaluation/all_models/report.md`

详细说明见：`evaluation/README.md`

