# AntibodyBench

用于“从头抗体模型评测”的统一框架，支持多模型接入、自动化执行与结果分析。

## 目录结构

```text
AntibodyBench/
├── algorithms/                   # 模型适配代码（每个模型 4 个核心文件）
├── inputs/
│   ├── antibody_datasets/        # 训练/测试划分与索引
│   └── alphafold3_inputs.json    # 抗体输入（AF3 风格）
├── outputs/
│   └── native_predictions/
├── scripts/
│   ├── run.sh                        # native inference orchestrator
│   ├── data_prep/
│   └── native_runners/
└── evaluation_tools/
```

## 快速开始

```bash
cd AntibodyBench
bash scripts/run.sh --model boltzgen --max-samples 1
```

## 第一步：准备数据索引与原生输入

```bash
python scripts/data_prep/generate_dataset_index.py
python scripts/data_prep/fetch_reference_complexes.py --input-csv inputs/antibody_datasets/dataset_index.csv --output-csv inputs/antibody_datasets/dataset_index.csv
python scripts/data_prep/fill_cdr_h3_from_anarci.py --input-csv inputs/antibody_datasets/dataset_index.csv --output-csv inputs/antibody_datasets/dataset_index_h3_annotated.csv
python scripts/data_prep/build_model_inputs_native.py
```

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

`outputs/native_predictions/<model>/<sample_id>/`

最小交付文件：

- `top1_structure.(pdb|cif)`
- `top1_sequence.fasta`
- `top1_meta.json`
- `run_stdout.log`
- `run_stderr.log`

模型级汇总：

- `outputs/native_predictions/<model>/manifest.csv`
