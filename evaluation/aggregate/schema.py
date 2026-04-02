from __future__ import annotations

"""
统一评估长表 schema。

约定：
- 缺失值统一写空字符串 ""（CSV 友好）。
- 布尔值统一写 "0"/"1"（便于 shell/csv 工具处理）。
- 数值字段如缺失同样写空字符串。
"""

EVALUATION_LONG_FIELDS = [
  "run_id",
  "run_tag",
  "model",
  "manifest_row_idx",
  "sample_id",
  "status",
  "error_summary",
  "duration_sec",
  "meta_path",
  "sequence_path",
  "structure_path",
  "pred_sequence",
  "pred_seq_len",
  "aa_valid_ratio",
  "cys_count",
  "motif_h3_len_in_8_20",
  "has_structure",
  "has_sequence",
  "reference_complex_path",
  "reference_complex_status",
  "antigen_chain_id",
  "antibody_chain_id",
  "cdr_h3_start",
  "cdr_h3_end",
  "cdr_h3_status",
  "cdr_h3_rmsd",
  "cdr_h3_atom_count",
  "struct_eval_eligible",
  "metric_error",
]


SUMMARY_BY_MODEL_FIELDS = [
  "model",
  "n_total",
  "n_success",
  "success_rate",
  "has_structure_rate",
  "has_sequence_rate",
  "median_runtime_sec",
  "median_seq_len",
  "median_cdr_h3_rmsd",
  "p25_cdr_h3_rmsd",
  "p75_cdr_h3_rmsd",
]


SUMMARY_BY_SAMPLE_FIELDS = [
  "sample_id",
  "n_models_attempted",
  "n_models_success",
  "best_model_by_cdr_h3_rmsd",
  "best_cdr_h3_rmsd",
]

