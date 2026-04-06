# Prepared Benchmark Data

`data/prepared/` 用来存放从 `data/raw/` 派生出来、但仍属于 benchmark 准备阶段的中间成果。

当前主流程只保留两类内容：

- `dataset_index_ready.csv`
  - 由 `python scripts/data_prep/prepare_native_inputs.py` 生成。
  - 来源是 `data/raw/dataset_index.csv`，并会补齐 `sample_id`、`reference_structure_path`、`reference_complex_path` 等字段。
  - 这是后续 native 输入构建、原生运行和统一评估默认使用的索引。
- `epitopes/<sample_id>.json`
  - 由 `extract_epitopes_from_complexes.py` 生成。
  - 保存每个样本的抗原热点、实际使用链 ID、修复状态等信息。

推荐工作流：

```bash
python scripts/data_prep/prepare_native_inputs.py
```

四个模型的实际运行输入不再放在 `data/` 里，而是统一写到顶层 `native_inputs/`。