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
