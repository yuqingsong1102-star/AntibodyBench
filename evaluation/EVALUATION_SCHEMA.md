# Evaluation Long Table Schema

文件：`outputs/evaluation/all_models/evaluation_long.csv`

## 字段说明

- `run_id`：唯一运行行标识，格式 `{run_tag}:{model}:{sample_id}:{manifest_row_idx}`。
- `run_tag`：运行批次标签，默认取 native 输出目录名（如 `native_predictions_real`）。
- `model`：模型名（`RFantibody` / `germinal` / `BindCraft` / `boltzgen`）。
- `manifest_row_idx`：该样本在 `manifest.csv` 的行序号（从 1 开始）。
- `sample_id`：样本 ID（如 `8q3j_B_A`）。
- `status`：样本状态（通常为 `ok` 或 `failed`）。
- `error_summary`：失败摘要。
- `duration_sec`：推理耗时（秒）。
- `meta_path` / `sequence_path` / `structure_path`：对应产物路径。
- `pred_sequence`：从 `top1_sequence.fasta` 提取的序列字符串。
- `pred_seq_len`：序列长度。
- `aa_valid_ratio`：序列中合法氨基酸比例。
- `cys_count`：序列中 Cys 数量。
- `motif_h3_len_in_8_20`：`pred_seq_len` 是否落在 `[8,20]`（`1/0`）。
- `has_structure` / `has_sequence`：结构/序列文件是否存在（`1/0`）。
- `reference_complex_path` / `reference_complex_status`：参考复合物信息（来自数据索引）。
- `antigen_chain_id` / `antibody_chain_id`：链 ID（来自数据索引）。
- `cdr_h3_start` / `cdr_h3_end` / `cdr_h3_status`：H3 区间与标注状态（来自数据索引）。
- `cdr_h3_rmsd` / `cdr_h3_atom_count`：结构评估结果。
- `tm_score`：TM-score 结构相似性指标。
- `dockq_score`：DockQ 复合物质量指标。
- `irmsd`：interface RMSD（由 DockQ 输出解析）。
- `struct_eval_eligible`：是否满足结构评估条件（`1/0`）。
- `metric_error_code`：指标计算失败的结构化错误码。
- `metric_error`：指标计算错误摘要。

## 产物分层

- `evaluation_long.csv`：入库后的统一长表
- `evaluation_long_seq.csv`：加 sequence 指标
- `evaluation_long_metrics.csv`：加 structure / docking 指标
- `summary_by_model.csv`：按模型聚合
- `summary_by_sample.csv`：按样本聚合
- `report.md`：Markdown 报告

## 缺失值规则

- 统一使用空字符串 `""` 表示缺失（包括数值列和文本列）。
- 布尔列统一用字符串 `1` / `0`。
- 即使样本失败也必须保留该行，禁止静默丢弃。

