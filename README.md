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
│   ├── input/
│   ├── prediction/
│   └── evaluation/
├── scripts/
│   ├── run.sh
│   └── aggregate_metrics.py
└── evaluation_tools/
```

## 快速开始

```bash
cd AntibodyBench
bash scripts/run.sh --model boltzgen
```

## 第一步：仅输入准备（不跑模型）

先完成输入准备里程碑，验证每个模型的输入是否齐全可追踪：

```bash
# 三个模型分别执行（不会触发推理）
bash scripts/run.sh --model germinal --prepare-only
bash scripts/run.sh --model boltzgen --prepare-only
bash scripts/run.sh --model BindCraft --prepare-only
```

v1 默认只依赖一个索引文件：

- `inputs/antibody_datasets/dataset_index.csv`

说明：

- 入口脚本 `scripts/run.sh` 默认就读取 `dataset_index.csv`。
- `alphafold3_inputs.json` 为可选增强信息；不存在时预处理会自动回退，不阻塞输入准备。

输入准备产物：

- `outputs/input/<model>/<sample_id>.json`
- `outputs/input/<model>/_manifest.json`

`_manifest.json` 包含：

- `total_samples`
- `ready_samples`
- `skipped_samples`
- `skip_reasons`
- `records`（逐样本状态与原因）

支持模型：

- `germinal`
- `boltzgen`
- `BindCraft`

## 模型接入说明

默认可用 `*_DRY_RUN=1` 快速打通流程；真实推理时可直接使用内置默认命令模板（也可手动覆盖）。

每个模型默认会通过 `conda run` 使用独立环境（推荐）：

- germinal: `GERMINAL_CONDA_ENV=germinal`
- boltzgen: `BOLTZGEN_CONDA_ENV=boltzgen`
- BindCraft: `BINDCRAFT_CONDA_ENV=BindCraft`

如需关闭 `conda run`（不推荐）可设置 `*_USE_CONDA_RUN=0`。

### Germinal torch 修复

若运行 `germinal` 时出现 `transformers` 要求 `torch>=2.6` 的错误，可执行：

```bash
bash scripts/fix_germinal_torch.sh
```

脚本会自动打印升级前后版本并校验 `torch >= 2.6.0`。

### 轻量调试建议

```bash
# boltzgen
BOLTZGEN_MAX_SAMPLES=1 \
BOLTZGEN_NUM_DESIGNS=2 \
BOLTZGEN_BUDGET=1 \
bash scripts/run.sh --model boltzgen

# BindCraft
BINDCRAFT_MAX_SAMPLES=1 \
BINDCRAFT_MAX_TRAJ=1 \
BINDCRAFT_NUM_FINAL=1 \
bash scripts/run.sh --model BindCraft
```

## 输入准备阶段验收清单（v1）

- 三个模型在 `--prepare-only` 下均返回退出码 0。
- 每个模型目录都生成 `_manifest.json`，且包含统计字段。
- 抽查 5 个 `ready` 样本：
  - 对应 `<sample_id>.json` 文件存在；
  - `records[*].input_json` 路径可访问。
- `skip` 样本必须有明确 `reason`，不能静默丢失。
