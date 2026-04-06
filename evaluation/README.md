# Evaluation（Unified External）

本目录按 `UNIFIED_EXTERNAL_EVALUATION_PLAN.md` 重建，目标是统一外部 judge 评估并输出多层评分卡，而不是 native-similarity 排名。

## 目录分层

- `ingest/build_external_eval_inputs.py`：构建 `external_eval_candidates.csv` 与 judge 任务模板
- `judge/run_primary_judge.py`：合并主 judge 结果或 fallback 占位
- `judge/run_secondary_judge.py`：可选 secondary judge 一致性导入
- `metrics/compute_binding_plausibility.py`：Viability + Binding plausibility
- `metrics/compute_affinity_proxy.py`：Affinity proxy 层
- `metrics/compute_developability.py`：Developability / novelty 字段
- `metrics/compute_robustness.py`：seed 稳健性字段
- `report/build_scorecards.py`：生成 `external_eval_sequence_level.csv`
- `report/build_model_comparison.py`：生成模型/样本汇总
- `pipeline/run_eval_pipeline.py`：一键编排

## 一键运行

```bash
python -m evaluation.pipeline.run_eval_pipeline --fallback-to-source-structure
```

常用参数：

- `--native-root <path>`：可重复指定；默认优先 `outputs/native_predictions_run2`
- `--index-csv`：默认 `data/prepared/dataset_index_ready.csv`
- `--out-dir`：默认 `outputs/evaluation/external_binding_benchmark`
- `--judge-manifest`：primary judge 结果 CSV（可选）
- `--secondary-judge-manifest`：secondary judge 结果 CSV（可选）
- `--judge-name`：默认 `chai1_primary`
- `--n-requested-seeds`：默认 `3`
- `--max-candidates-per-model-sample`：默认 `10`，先做 top-K 唯一序列筛选再进入评估
- `--judge-cmd`：直接执行外部 primary judge 的命令模板（可选）

## 默认输出

- `external_eval_candidates.csv`
- `external_eval_candidates_selection_audit.csv`
- `external_eval_seed_level.csv`
- `external_eval_sequence_level.csv`
- `external_eval_model_summary.csv`
- `external_eval_sample_summary.csv`
- `report.md`
- `figures/`
- `judge_inputs/`
- `primary_judge_runs/`

## 候选筛选规则

- 官方候选集只保留每个 `model + sample` 下的 top-K 唯一序列
- 去掉空序列、非法氨基酸序列、重复序列、以及超出 top-K 的尾部候选
- 所有保留/丢弃原因会写进 `external_eval_candidates_selection_audit.csv`

## 直接执行 Primary Judge

如果你不想先手工准备 `judge_manifest`，可以直接传 `--judge-cmd`。命令模板支持这些占位符：

- `__JOB_JSON__`
- `__JUDGE_SEED__`
- `__OUTPUT_DIR__`
- `__RESULT_JSON__`
- `__CANDIDATE_ID__`
- `__SAMPLE_ID__`
- `__MODEL__`
- `__PRED_SEQUENCE__`
- `__TARGET_STRUCTURE_PATH__`
- `__TARGET_CHAIN_ID__`
- `__DESIGNED_ANTIBODY_CHAIN_ID__`
- `__HOTSPOT_STRING__`

外部命令建议在 `__RESULT_JSON__` 写回一个结果 JSON，至少包含：

- `status` 或 `external_status`
- `structure_path` 或 `external_structure_path`

可选字段：

- `confidence_path` 或 `external_confidence_path`
- `iptm` / `external_iptm`
- `ipae` / `external_ipae`
- `plddt_binder` / `external_plddt_binder`
- `error_summary` / `external_error_summary`
- `target_chain_id`
- `designed_antibody_chain_id`

如果结果 JSON 没写结构路径，评估器会尝试在对应 seed 输出目录里自动寻找 `.cif/.mmcif/.pdb` 文件。

## 兼容入口

- `ingest/build_external_eval_candidates.py`
- `metrics/compute_phase1_scorecard.py`
- `report/build_external_summaries.py`

这些脚本会转调新实现，便于旧命令平滑过渡。
