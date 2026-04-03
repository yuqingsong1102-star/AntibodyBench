# AntibodyBench

用于“从头抗体模型评测”的统一框架，支持多模型接入、自动化执行与结果分析。

## 目录结构

```text
AntibodyBench/
├── data/                        # 参考结构、ready 样本对应原始数据、模型输入
├── inputs/                      # 样本索引 CSV（base / annotated / ready）
├── outputs/                     # native 推理结果与评估结果（生成物）
├── scripts/
│   ├── run.sh                  # 原生推理编排入口
│   ├── data_prep/              # 数据准备脚本
│   ├── native_runners/         # 四个模型的原生 runner
│   └── tools/                  # 项目自带 TMscore / DockQ 工具
├── evaluation/                  # ingest / metrics / report / pipeline 四层评估体系
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
python scripts/data_prep/generate_dataset_index.py
python scripts/data_prep/fetch_reference_complexes.py --input-csv inputs/dataset_index.csv --output-csv inputs/dataset_index.csv
python scripts/data_prep/fill_cdr_h3_from_anarci.py --input-csv inputs/dataset_index.csv --output-csv inputs/dataset_index_h3_annotated.csv
python scripts/data_prep/build_model_inputs_native.py
```

三个索引文件的含义：

- `inputs/dataset_index.csv`：原始样本索引
- `inputs/dataset_index_h3_annotated.csv`：加入 CDR-H3 标注状态后的审核表
- `inputs/dataset_index_ready.csv`：过滤后的稳定可运行样本子集

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

- `top1_structure.(pdb|cif)`
- `top1_sequence.fasta`
- `top1_meta.json`
- `run_stdout.log`
- `run_stderr.log`

模型级汇总：

- `outputs/native_predictions_<run_tag>/<model>/manifest.csv`

## 第三步：运行统一评估

```bash
python -m evaluation.pipeline.run_eval_pipeline \
	--native-root outputs/native_predictions_run2 \
	--index-csv inputs/dataset_index_ready.csv \
	--out-dir outputs/evaluation/all_models
```

默认情况下，评估管线会优先读取 `outputs/native_predictions_run2`；若不存在，再回退到 `real + smoke` 历史目录。

## 评估工具

- `TMscore` 默认使用项目自带二进制：`scripts/tools/bin/TMscore`
- `DockQ` 默认使用项目内包装脚本：`scripts/tools/bin/DockQ.py`
- 两者都不依赖系统 PATH 中的同名命令
