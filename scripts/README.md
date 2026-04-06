# Scripts Index

## Native inference entry

`scripts/run.sh` 已改为**原生推理编排入口**，不再执行 preprocess/postprocess/aggregate 的三阶段评测流程。
新的输出契约不再压成 Top-1，而是保留**候选清单**。

### 用法

```bash
bash scripts/run.sh --model <RFantibody|germinal|BindCraft|boltzgen> [options]
```

参数说明：

- `--model`：必填，模型名。
- `--input-root`：原生输入根目录，默认 `native_inputs`。
- `--out-root`：输出根目录，默认 `outputs/native_predictions`。
- `--max-samples`：最多处理样本数，默认 `1`（即默认 smoke_1）。
- `--sample-id`：仅运行一个样本（会忽略 `--max-samples`）。
- `--continue-on-error`：`0|1`，默认 `1`，即 `sequential_continue`。

### 输入来源

仅从以下目录读取样本输入，不再依赖适配层 `outputs/input`：

`native_inputs/<model>/<sample_id>/`

### 运行器

模型执行由 `scripts/native_runners/` 下的薄封装脚本负责：

- `rfantibody.sh`
- `germinal.sh`
- `bindcraft.sh`
- `boltzgen.sh`
- `collect_candidates.py`（统一候选清单抽取）

旧的 `collect_top1.py` 已移到 `scripts/legacy/collect_top1.py`，仅作为旧流程备份，不再属于正式主流程。

### 输出契约

每个样本输出目录：

`outputs/native_predictions/<model>/<sample_id>/`

最小交付文件：

- `candidate_manifest.csv`
- `run_meta.json`
- `native_run/`（原始模型输出，不做 Top-1 抽取）
- `run_stdout.log`
- `run_stderr.log`

模型级汇总：

- `outputs/native_predictions/<model>/manifest.csv`

`candidate_manifest.csv` 与模型级 `manifest.csv` 使用同一字段集合，核心字段如下：

- `candidate_id`
- `sample_id`
- `model`
- `candidate_rank`
- `candidate_name`
- `status`
- `sequence_path`
- `structure_path`
- `meta_path`
- `duration_sec`
- `error_summary`
- `source_stage`
- `native_score_name`
- `native_score_value`

说明：

- `sequence_path` 指向每个候选自己的序列文件；如果模型原始输出只有 CSV 序列列，抽取器会生成标准化 FASTA sidecar。
- `structure_path` 是可选字段；没有结构也允许保留候选。
- `meta_path` 指向每个候选的元数据 JSON，而不是样本级 Top-1 元数据。
- `manifest.csv` 每次运行都会**覆盖重写**，不再向旧结果追加。

### 失败与可追溯

- 单样本失败默认不中断后续样本（`--continue-on-error 1`）。
- 外部命令标准输出/错误分别写入 `run_stdout.log` / `run_stderr.log`。
- 即使运行环境不完整，也会产出 `candidate_manifest.csv`、`run_meta.json` 和模型级 `manifest.csv` 的失败记录。
- 不再在 scripts 层做 Top-1 选择，也不再从结构反推序列补全候选。

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
- `outputs/native_predictions_smoke/<model>/<sample_id>/candidate_manifest.csv`

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

默认 runner 会把 `settings_target.json` 中的 `design_path` 在运行时重写到 `__OUTPUT_DIR__/designs` 对应的真实 sample 输出目录，避免候选结构落回输入目录。

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

当前 data_prep 主流程只保留 3 个脚本：

- `scripts/data_prep/prepare_native_inputs.py`
	- 一键入口：从 `data/raw/dataset_index.csv` 出发，补齐 full complex、生成 ready 索引、抽取 epitope，并产出四模型原生输入。
	- 默认写入 `data/prepared/dataset_index_ready.csv`、`data/prepared/epitopes/` 和 `native_inputs/`。
- `scripts/data_prep/extract_epitopes_from_complexes.py`
	- 读取 `data/prepared/dataset_index_ready.csv`，为每个样本生成 `data/prepared/epitopes/<sample_id>.json`。
	- 当 full complex 链名与数据集链名不一致时，会按接触关系、链类型与 CDR3 风格自动回退到更合理的链映射。
	- 可配合 `--repair-target-pdbs`，把缺失或空壳的 `reference_structure_path` 目标链从 full complex 反抽重建。
- `scripts/data_prep/build_model_inputs_native.py`
	- 直接读取 `data/prepared/dataset_index_ready.csv` 与 `data/prepared/epitopes/<sample_id>.json`，一次生成四个模型的 native 输入。
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

