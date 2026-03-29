# algorithms 接入说明

每个模型目录都应包含以下 4 个核心文件：

- `preprocess.py`
- `make_predictions.sh`
- `postprocess.py`
- `container.def`

当前已接入：

- `germinal`
- `boltzgen`
- `BindCraft`

三者统一支持：

- `*_DRY_RUN=1`：用输入PDB占位输出，先打通流程
- `*_CMD`：真实推理命令模板（可选；不设置时使用内置默认模板）
- `*_WORKDIR`：执行命令时的工作目录

默认会在各自 conda 环境执行（可覆盖）：

- germinal：`GERMINAL_CONDA_ENV`（默认 `germinal`）/ `GERMINAL_USE_CONDA_RUN`（默认 `1`）
- boltzgen：`BOLTZGEN_CONDA_ENV`（默认 `boltzgen`）/ `BOLTZGEN_USE_CONDA_RUN`（默认 `1`）
- BindCraft：`BINDCRAFT_CONDA_ENV`（默认 `BindCraft`）/ `BINDCRAFT_USE_CONDA_RUN`（默认 `1`）

germinal 还支持：

- `GERMINAL_MAX_SAMPLES`：限制本次处理样本数（0 表示不限制）
- `GERMINAL_MAX_TRAJ` / `GERMINAL_MAX_HALLUCINATED` / `GERMINAL_MAX_PASSING`

boltzgen 还支持：

- `BOLTZGEN_MAX_SAMPLES`
- `BOLTZGEN_NUM_DESIGNS` / `BOLTZGEN_BUDGET` / `BOLTZGEN_PROTOCOL` / `BOLTZGEN_BINDER_LEN`

BindCraft 还支持：

- `BINDCRAFT_MAX_SAMPLES`
- `BINDCRAFT_MAX_TRAJ` / `BINDCRAFT_NUM_FINAL`
- `BINDCRAFT_LENGTH_MIN` / `BINDCRAFT_LENGTH_MAX`
- `BINDCRAFT_FILTERS` / `BINDCRAFT_ADVANCED_TEMPLATE`
