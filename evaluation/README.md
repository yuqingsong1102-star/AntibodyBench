# Evaluation (Track B)

本目录用于“先分析、后补全结果”的半自动评估流程。  
即使当前只有部分模型结果，也能产出统一总表与可视化报告。

## 目录分层

- `ingest/`：入库层，把 native outputs 汇总为统一长表
- `metrics/`：指标层，按行计算 sequence / RMSD / TM-score / DockQ / iRMSD
- `report/`：报告层，生成 summary、figure、report
- `pipeline/`：编排层，一键串起整个评估流程

## 一键运行

在仓库根目录执行：

```bash
python -m evaluation.pipeline.run_eval_pipeline
```

默认读取：

- `outputs/native_predictions_run2`（若存在，优先使用）
- 否则回退到 `outputs/native_predictions_real` + `outputs/native_predictions_smoke`
- `inputs/dataset_index_h3_annotated.csv`

## 输出目录

默认输出到：`outputs/evaluation/all_models`

核心文件：

- `evaluation_long.csv`：拼接后的统一长表（run 级）
- `evaluation_long_seq.csv`：补充 Level-0/1 序列指标
- `evaluation_long_metrics.csv`：补充结构指标（当前含 `cdr_h3_rmsd`、`tm_score`、`dockq_score`、`irmsd`）
- `summary_by_model.csv`：模型级汇总
- `summary_by_sample.csv`：样本级汇总
- `report.md`：一页报告
- `figures/`：5 张固定图（成功率、RMSD、耗时、热力图、Pareto）

## 分步执行（可选）

```bash
python -m evaluation.ingest.build_evaluation_long
python -m evaluation.metrics.compute_sequence_metrics
python -m evaluation.metrics.compute_cdr_h3_rmsd
python -m evaluation.metrics.compute_tm_score
python -m evaluation.metrics.compute_dockq
python -m evaluation.report.build_summaries
python -m evaluation.report.visualize_metrics
python -m evaluation.report.write_report
```

## 当前实现范围

- 已实现：
  - 长表拼接（`manifest + top1 + dataset_index`）
  - 完整性指标：`success_rate`、`has_structure_rate`、`has_sequence_rate`、`median_runtime_sec`
  - 序列指标：`pred_seq_len`、`aa_valid_ratio`、`cys_count`、`motif_h3_len_in_8_20`
  - 结构指标：`cdr_h3_rmsd`、`tm_score`、`dockq_score`、`irmsd`
  - 模型/样本汇总：`summary_by_model.csv`、`summary_by_sample.csv`
  - 5 张标准图 + `report.md`
- 未实现（后续可扩展）：
  - `cdr_l3_rmsd`、全局 RMSD、更多 interface/confidence 指标等

## 工具依赖

- `TMscore` 默认使用项目自带二进制：`scripts/tools/bin/TMscore`
- `DockQ` 默认使用项目内包装脚本：`scripts/tools/bin/DockQ.py`
- `DockQ.py` 会优先从 `scripts/tools/vendor/` 加载项目本地安装的 DockQ 包，而不是依赖系统 PATH

## 兼容说明

若环境里 `matplotlib` 与 `numpy` 不兼容，脚本会降级为占位图并继续输出 CSV 与 `report.md`，不会中断评估流程。

