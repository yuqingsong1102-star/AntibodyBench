# Evaluation v2 使用说明（详细版）

本目录提供 AntibodyBench 的统一评估流水线，用来对不同生成模型产出的候选进行**可复现、可拆分、可增量**的后处理评估。

核心原则：

- 以 `candidate_id` 作为统一主键，逐步补齐指标；
- 各 step 可独立执行，也可通过统一入口串行执行；
- 指标以连续值为主，不在本层做强制二分类筛选；
- 支持缓存/断点续跑（尤其是 Step 1 Chai-1）。

---

## 1. 流水线总览（当前实现）

当前 `run_evaluation.py` 实际执行为 **4 步 + 2 次合并**：

```text
模型候选（native outputs） + index.csv
            │
            ├─ Step 1: Chai-1 Judge（GPU）
            │    输出: chai1_results.csv, chai1_runs/
            │    关键指标: chai1_iptm, chai1_ptm, chai1_plddt_binder
            │
            ├─ Step 2: Interface Metrics（CPU，可选 PyRosetta）
            │    输出: interface_results.csv
            │    关键指标: interface_bsa, pdockq2, clash/contact, approx_dG, (rosetta_*)
            │
            ├─ Step 3: Sequence Metrics（CPU）
            │    输出: sequence_results.csv
            │    关键指标: liabilities, charge, entropy, CDR 长度, VHH hallmarks
            │
            ├─ Aggregate
            │    输出: eval_merged.csv
            │    作用: 按 candidate_id 合并 1/2/3 步
            │
            ├─ Step 4: Diversity + Epitope（CPU + 可选 MMseqs2）
            │    输出: diversity_epitope_results.csv
            │    关键指标: 聚类多样性、表位覆盖率
            │
            └─ Final Merge
                 输出: eval_final.csv
                 作用: eval_merged + diversity_epitope
```

---

## 2. 目录与输入契约

### 2.1 `--native-root` 目录约定

代码扫描逻辑默认你给的每个 `native_root` 下是如下结构：

```text
<native_root>/
  <model_name>/
    <sample_id>/
      candidate_manifest.csv
```

其中 `candidate_manifest.csv` 至少应包含（代码实际读取）：

- `candidate_id`
- `status`（仅 `ok` 会进入评估）
- `sequence_path`（可选）
- `structure_path`（可选，但很多流程会依赖）

### 2.2 `--index-csv` 约定

`index_csv` 用于补充 target 元信息。代码中实际会用到的列包括：

- `sample_id`
- `reference_complex_path`（target 结构路径）
- `antigen_chain`
- `antibody_chain`
- `hotspot_string`（Step2 hotspot 联系计数）

如果部分列缺失，对应指标会降级/留空，但流程通常仍可继续。

---

## 3. 快速开始

### 3.1 跑完整流程（推荐）

```bash
python -m evaluation.run_evaluation \
  --native-root outputs/final_boltzgen \
  --native-root outputs/native_predictions_cropped_formal/RFantibody \
  --index-csv data/prepared/dataset_index_ready.csv \
  --out-dir outputs/evaluation/results
```

### 3.2 只跑 Step 1（需要 GPU）

```bash
python -m evaluation.run_evaluation \
  --step 1 \
  --native-root outputs/final_boltzgen \
  --device cuda:0 \
  --out-dir outputs/evaluation/results
```

### 3.3 只跑 Step 2 + Step 3（CPU）

```bash
python -m evaluation.run_evaluation \
  --step 2 --step 3 \
  --native-root outputs/final_boltzgen \
  --out-dir outputs/evaluation/results
```

### 3.4 跑 Step 4（多样性 + 表位覆盖）

```bash
python -m evaluation.run_evaluation \
  --step 4 \
  --native-root outputs/final_boltzgen \
  --out-dir outputs/evaluation/results \
  --epitope-dir data/epitopes
```

### 3.5 启用 PyRosetta（Step 2 增强）

```bash
python -m evaluation.run_evaluation \
  --step 2 \
  --use-pyrosetta \
  --native-root outputs/final_boltzgen \
  --out-dir outputs/evaluation/results
```

---

## 4. 统一入口 CLI 说明（`run_evaluation.py`）

| 参数 | 类型 | 默认值 | 说明 |
|---|---|---|---|
| `--native-root` | 可重复 Path | 无（必填） | 模型输出根目录，可传多个 |
| `--index-csv` | Path | `data/prepared/dataset_index_ready.csv` | 数据索引 |
| `--out-dir` | Path | `outputs/evaluation/results` | 评估输出目录 |
| `--step` | 可重复 int | 全部（1,2,3,4） | 指定执行步骤 |
| `--device` | str | `cuda:0` | Step1 的 Chai-1 设备 |
| `--model` | 可重复 str | None | 仅评估指定模型名 |
| `--max-candidates` | int | `0` | Step1 候选上限（0=不限） |
| `--use-pyrosetta` | flag | False | Step2 追加 Rosetta 指标 |
| `--epitope-dir` | Path | None | Step4 表位 JSON 目录 |

说明：

- `--step` 不传时会跑 1/2/3/4 全部步骤；
- 任何步骤执行后都会做一次 `eval_merged.csv` 聚合；
- 仅当包含 Step4 时才会额外生成 `eval_final.csv`。

---

## 5. Step 1：Chai-1 Judge（`steps/chai1_judge.py`）

### 5.1 功能

- 收集候选（来自 `candidate_manifest.csv`）；
- 准备 target+binder FASTA；
- 调用 `chai_lab.chai1.run_inference`；
- 提取并保存置信度指标；
- 支持两级缓存（`out_csv` 和 `run_dir/<cid>/result.json`）。

### 5.2 输出字段（`chai1_results.csv`）

- 基本：`candidate_id`, `sample_id`, `model`, `pred_sequence`
- 状态：`chai1_status`, `chai1_error`
- 置信度：`chai1_iptm`, `chai1_ptm`, `chai1_plddt_binder`
- 路径：`chai1_structure_path`, `chai1_confidence_path`
- target 元数据：`target_structure_path`, `target_chain_id`, `binder_chain_id`

### 5.3 关键实现细节

- 仅处理 `manifest` 中 `status=ok` 的候选；
- `pred_sequence` 获取优先级：
  1) `sequence_path` 文件；
  2) 从 `structure_path` 解析（优先提示链，否则取最短链）；
- Chai-1 输出按 `sorted()` 取 top 结构；
- binder pLDDT 取 `per_chain_plddt` 的最后一条（因为 FASTA 中 binder 后写入）；
- 失败会写 `chai1_status=failed` + 错误类型，便于排错和续跑。

### 5.4 环境变量（Step1 生效）

- `CHAI_NUM_DIFFN_SAMPLES`（默认 `5`）
- `CHAI_NUM_TRUNK_RECYCLES`（默认 `3`）
- `CHAI_NUM_DIFFN_TIMESTEPS`（默认 `200`）
- `CHAI_USE_MSA_SERVER`（`1/true` 开启）
- `CHAI_LOW_MEMORY`（默认开启，设 `0/false` 可关闭）

---

## 6. Step 2：Interface Metrics（`steps/interface_metrics.py`）

### 6.1 功能

从复合物结构中计算界面几何与近似能量指标，默认基于 BioPython；可选调用 PyRosetta 增强。

### 6.2 结构来源策略

- 若提供且可用 `chai1_results.csv`：默认优先 `chai1_structure_path`；
- 原生结构仅在 `prefer_chai1_structure=False`（CLI 对应 `--prefer-native`）时作为补充；
- 若结构缺失或链推断失败，对应候选标记为 `interface_status=failed`。

### 6.3 输出字段（`interface_results.csv`）

- 状态：`interface_status`, `interface_error`
- 接触/碰撞：`interface_residue_count`, `antigen_contact_residue_count`, `interface_contact_pair_count`, `interface_clash_count`, `interface_clash_ratio`
- 界面面积：`interface_bsa`, `binder_bsa_fraction`
- 质量：`pdockq2`, `avg_interface_plddt`, `binder_mean_plddt`
- 表面属性：`surface_hydrophobicity`, `interface_fragmentation_index`
- hotspot：`hotspot_contact_count`
- 近似能量：`approx_dG`, `approx_hbond_count`, `approx_saltbridge_count`
- Rosetta（可选）：`rosetta_ddG`, `rosetta_sc`, `rosetta_dSASA`, `rosetta_hbonds_total`
- 来源：`structure_source`, `structure_path`

### 6.4 关键阈值与公式（实现口径）

- 接触阈值：`5.0 Å`（`CONTACT_CUTOFF`）
- 碰撞阈值：`2.0 Å`（`CLASH_CUTOFF`）
- 表面残基阈值：`SASA >= 15.0`
- pDockQ2：
  - 以 `CB`（无则 `CA`）构建跨链 8Å 接触；
  - `x = avg_plddt * log10(n_contacts)`；
  - `pdockq2 = 0.724 / (1 + exp(-0.052*(x-152.611))) + 0.018`
- `approx_dG` 为基于 BSA/contact/clash 的经验近似（非 Rosetta 正式能量）。

---

## 7. Step 3：Sequence Metrics（`steps/sequence_metrics.py`）

### 7.1 功能

纯序列指标分析，不依赖 GPU，不做结构重预测。

### 7.2 数据来源

- 优先读取 `chai1_results.csv` 的 `pred_sequence`；
- 缺失时回退 `manifest` 的 `sequence_path`；
- 再缺失时尝试从结构提取链序列（链提示失败则取最短链）。

### 7.3 输出字段（`sequence_results.csv`）

- 状态/基础：`sequence_status`, `pred_sequence`, `seq_length`
- 理化：`net_charge`, `isoelectric_point_approx`, `molecular_weight_approx`, `gravy`, `sequence_entropy`
- liability：`free_cys_count`, `glycosylation_motif_count`, `deamidation_motif_count`, `isomerization_motif_count`, `total_liability_count`
- 组成：`aromatic_fraction`, `charged_fraction`, `hydrophobic_fraction`
- VHH 特征：`has_vhh_hallmarks`, `framework_gapped_positions`
- CDR：`cdr1_length`, `cdr2_length`, `cdr3_length`, `cdr3_sequence`

### 7.4 说明

- pI 是粗略近似，不是精确滴定模型；
- CDR 长度是启发式估计（非严格编号体系）；
- `has_vhh_hallmarks` 使用 FR2 区段正则近似判定。

---

## 8. 聚合：`aggregate.py`

### 8.1 `eval_merged.csv`

- 按 `candidate_id` 对 Step1/2/3 做并集；
- 先放通用列：`candidate_id`, `sample_id`, `model`；
- 再拼接各 step 字段（避免重复列）；
- 若某 step 缺失，字段留空字符串。

### 8.2 `eval_final.csv`

- 先运行 Step4 生成 `diversity_epitope_results.csv`；
- 再把 Step4 新字段按 `candidate_id` 合并到 `eval_merged.csv`；
- 输出为 `eval_final.csv`。

---

## 9. Step 4：多样性 + 表位覆盖（`steps/diversity_epitope.py`）

### 9.1 多样性（Diversity）

以 `(sample_id, model)` 分组，组内按序列聚类：

- 优先使用 MMseqs2（`mmseqs cluster`）；
- 参数：`--min-seq-id 0.9`, `-c 0.8`, `--cov-mode 0`；
- 若找不到 `mmseqs` 或组内序列不足 2 条，则该组聚类指标降级（`cluster_id=0` 等）。

输出字段：

- `cluster_id`, `cluster_size`
- `group_cluster_count`, `group_total_seqs`
- `group_largest_cluster_frac`
- `group_singleton_frac`

### 9.2 表位覆盖（Epitope Coverage）

- 读取 `epitope_dir/<sample_id>.json` 的 hotspot 残基集合；
- 在结构中统计 binder 与 hotspot（antigen 残基）是否 < 5Å 接触；
- 输出：
  - `epitope_contact_count`
  - `epitope_total_hotspots`
  - `epitope_coverage = contact_count / total_hotspots`

支持多种 JSON key（实现中）：

- `hotspot_residues`
- `residues`
- `epitope_residues`
- `hotspot_indices`
- 或直接 JSON list

---

## 10. 结果文件清单

默认输出到 `--out-dir`：

- `chai1_results.csv`
- `interface_results.csv`
- `sequence_results.csv`
- `eval_merged.csv`
- `diversity_epitope_results.csv`（执行 Step4 时）
- `eval_final.csv`（执行 Step4 时）
- `chai1_runs/`（Step1 中间产物与缓存）

---

## 11. 断点续跑与缓存机制

### Step1 缓存

- 若 `chai1_results.csv` 已存在且某 `candidate_id` 状态是 `ok`，会直接复用；
- 若 `chai1_runs/<cid>/result.json` 记录 `status=ok`，也会复用；
- 每处理一个候选会增量写 CSV，崩溃后可继续。

### Step2/Step3/Step4

- 当前实现默认每次全量重算并覆盖目标 CSV；
- 如果只想“补跑”，建议先用 `--step` 限定步骤，或上游裁剪输入候选。

---

## 12. 依赖与运行要求

基础依赖（按代码导入）：

- Python 3.10+（建议）
- `numpy`
- `biopython`

可选依赖：

- `chai_lab`（Step1 必需）
- `pyrosetta`（Step2 `--use-pyrosetta` 时）
- `mmseqs2` 可执行文件（Step4 多样性聚类）

设备要求：

- Step1 需要 GPU（`--device`）
- Step2/3/4 默认 CPU 即可

---

## 13. 常见问题与排错

### 13.1 Step1 全是 `no_sequence` / `no_target_structure`

- 检查 `candidate_manifest.csv` 中 `sequence_path` / `structure_path`；
- 检查 `index_csv` 的 `reference_complex_path` 是否存在；
- 确认 `sample_id` 能在 `index_csv` 中匹配到。

### 13.2 Step2 链解析失败（`cannot_resolve_chains`）

- 优先保证 `index_csv` 中 `antigen_chain` / `antibody_chain` 正确；
- 确保结构文件里确实有对应链；
- 复杂链重命名场景建议先用 Chai-1 结构作为输入。

### 13.3 Step4 显示 `mmseqs2 unavailable`

- 安装 `mmseqs` 并确保在 PATH；
- 或放置在 `~/.local/bin/mmseqs`。

### 13.4 Step4 没有 epitope 覆盖结果

- 确认 `--epitope-dir` 指向的目录存在；
- 确认文件命名为 `<sample_id>.json`；
- JSON 中 hotspot key 是否在支持列表内。

---

## 14. 代码结构

```text
evaluation/
├── README.md
├── __init__.py
├── run_evaluation.py          # 统一入口（Step1~4 + 聚合）
├── aggregate.py               # eval_merged / eval_final 合并逻辑
├── structure_inference.py     # 结构解析与链推断工具
├── schema.py                  # 外部统一 schema 常量（供其它模块使用）
└── steps/
    ├── chai1_judge.py         # Step1
    ├── interface_metrics.py   # Step2
    ├── sequence_metrics.py    # Step3
    └── diversity_epitope.py   # Step4
```

---

## 15. 与 `analysis/` 的分工

- `evaluation/`：生产标准化指标 CSV（数据生产层）
- `analysis/`：消费 CSV 做统计图、对比图、论文图（数据分析层）

建议流程：先稳定 `evaluation` 输出口径，再在 `analysis` 中迭代图表和报告逻辑。
