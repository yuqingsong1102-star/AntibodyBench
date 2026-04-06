# Unified External Evaluation Plan

## 1. 背景与结论

当前评估体系默认把“与原抗体复合物是否相似”当成主问题，因此主指标集中在 `cdr_h3_rmsd`、`tm_score`、`dockq_score`、`irmsd` 这一类 reference-based 指标上。

这套思路不再作为正式评估体系，原因有三点：

- 新设计与原抗体越像，不等于结合越好。
- 新设计与原抗体完全不同，也不等于不能结合同一靶点。
- 用原复合物相似性做主指标，会系统性偏向“复刻 native binder”，而不是评估“是否形成可信的新 binder”。

因此，本项目正式改用统一的外部评估体系，主问题改为：

1. 该模型生成的抗体，是否在统一外部 judge 下表现出可信的结合可能性。
2. 如果可能结合，其界面质量和结合强度代理是否足够好。
3. 该结论是否稳健，而不是某一次预测偶然得到的好看结构。

本方案直接删除“与原抗体像不像”这一块在正式评分卡中的地位。原有 reference-based 指标不再进入官方主分面板，不再用于模型排序，不再用于核心结论。

## 2. 评估范围

- 任务范围：仅覆盖 antibody-antigen。
- 模型范围：RFantibody、BoltzGen、Germinal、BindCraft。
- 输入对象：模型生成的抗体序列为主，模型自带结构只作为辅助手段，不作为最终 binding 主评判依据。
- 输出目标：统一的外部评分卡、模型级比较表、样本级比较表、通过率与稳健性报告。

当前阶段明确不做以下事情：

- 不用原抗体复合物相似度作为主评价目标。
- 不在第一阶段依赖实验 KD、BLI、SPR 数据。
- 不在第一阶段用单一总分替代多维评分卡。
- 不扩展到 protein-protein、protein-ligand、protein-RNA、protein-DNA 等其他任务。

## 3. 核心原则

新的统一外部评估体系遵守以下原则：

- 所有模型必须由同一套外部 judge 重评，不能直接使用各模型内部自评结果做横向胜负判断。
- 主评价对象是“能否形成可信结合界面”，不是“是否重现 native binder”。
- 所有高层结论必须以外部 judge 的统一输出为准。
- 模型自带分数只能保留为 metadata，用于解释模型行为，不能作为官方比较主依据。
- 先做多维评分卡，再决定是否需要单一综合分；第一阶段不允许直接压成总分。
- 所有阈值都必须通过内部校准集确定，禁止拍脑袋给 final cutoff。

## 4. 统一外部 Judge 设计

### 4.1 官方 Judge

第一阶段官方 judge 定义如下：

- 主 judge：外部共折叠预测器，统一输入 `target structure + designed antibody sequence`，输出预测复合物结构与置信度。
- 推荐实现：`Chai-1` 作为默认主 judge。
- 可选第二 judge：`AF3`，若资源允许则作为稳健性复核 judge，而不是第一阶段硬依赖。

选择 Chai-1 作为默认主 judge 的原因：

- 可重复。
- 可以统一重评四个模型的输出。
- 不直接继承四个模型各自的内部偏置。
- 相比把各模型自带结构直接拿来评分，更适合作为外部评委。

### 4.2 Judge 输入标准化

对四个模型一律统一成如下输入：

- 靶点结构：仅保留 antigen 链。
- 抗体输入：模型生成的 top-K 抗体序列。
- 链定义：统一使用 `target_chain_id` 与 `designed_antibody_chain_id`。
- 每条序列在外部 judge 下运行固定数量的 seeds。

建议默认设置：

- 每个样本每个模型取 `K=10` 条候选序列进入外部评估。
- 每条序列用主 judge 重评 `S=3` 次不同 seed。
- 若启用第二 judge，则第二 judge 只重评主 judge 最好的 `top-M` 条候选，建议 `M=3`。

### 4.3 Judge 输出

统一保存以下产物：

- 外部 judge 预测复合物结构。
- 外部 judge 置信度原始输出。
- 结构几何与界面统计表。
- 能量代理与 developability 统计表。
- 每条序列的 seed 级结果。
- seed 聚合后的 sequence 级结果。

## 5. 正式评分卡结构

新体系不使用单个总分，而使用五层评分卡。

### A. Viability Gate

这层只回答“该候选是否有资格进入结合评估”。

建议字段：

- `has_valid_sequence`
- `has_external_structure`
- `sequence_valid_ratio`
- `sequence_length`
- `cys_count`
- `num_missing_residues`
- `chain_parse_ok`
- `external_prediction_ok`
- `severe_clash_flag`
- `basic_liability_flag`

判定逻辑：

- 任何候选若未通过 viability gate，不再进入 binding 主评分。
- 失败原因必须明确落盘，不能被静默过滤。

### B. Binding Plausibility

这层是正式主层，回答“它像不像一个真的在结合的复合物”。

建议主指标：

- `external_iptm`
- `external_ipae`
- `external_plddt_binder`
- `pdockq2`
- `interface_residue_count`
- `interface_contact_pair_count`
- `interface_bsa`
- `interface_binder_fraction`
- `interface_clash_count`
- `interface_clash_ratio`
- `interface_shape_complementarity`
- `interface_packstat`
- `hotspot_contact_count`
- `antigen_contact_residue_count`
- `interface_fragmentation_index`

这层的主旨是：

- 必须先确认存在像样的跨链接口。
- 必须确认接口不是靠严重碰撞或虚假接触支撑起来的。
- 必须确认接口对靶标有足够接触，而不是偶然擦边。

### C. Affinity Proxy

这层回答“如果它在结合，那这个界面从代理能量上看强不强”。

建议指标：

- `interface_dG`
- `interface_dSASA`
- `interface_dG_dSASA_ratio`
- `interface_hbond_count`
- `interface_saltbridge_count`
- `interface_unsat_hbond_count`
- `interface_hydrophobicity`
- `surface_hydrophobicity`
- `binder_score`
- `sap_score`

说明：

- 这一层不是亲和力真值，只是亲和力代理层。
- 只在候选已通过 binding plausibility 后才解释其高低。
- 这一层不应单独决定“是否能结合”，只能决定“优先级更高还是更低”。

### D. Robustness

这层回答“这个好结果稳不稳”。

建议指标：

- `judge_seed_pass_rate`
- `judge_seed_metric_std`
- `best_seed_vs_median_gap`
- `secondary_judge_agreement`
- `interface_residue_jaccard_across_seeds`
- `epitope_consistency_across_seeds`

解释规则：

- 单次 seed 很好、但多次 seed 波动很大，不应视为高可信 binder。
- 主 judge 与第二 judge 一致时，可信度显著提高。

### E. Developability And Novelty

这层不回答“能不能结合”，而回答“值不值得做”。

建议指标：

- `framework_mutation_count`
- `surface_hydrophobicity`
- `free_cys_flag`
- `glycosylation_motif_flag`
- `deamidation_risk_flag`
- `isomerization_risk_flag`
- `liability_motif_count`
- `net_charge`
- `sequence_entropy`
- `accepted_set_diversity`

说明：

- developability 只能影响优先级，不应推翻明确的 binding 证据。
- novelty 只作为辅助分析，不再用“是否像原抗体”作为好坏标准。

## 6. 原有相似性指标的处理原则

以下指标从正式主评分体系中删除：

- `cdr_h3_rmsd`
- `tm_score`
- `dockq_score`
- `irmsd`
- 任意以原抗体复合物为中心的全局/局部 RMSD
- 任意以原抗体表位重合度为主排序依据的指标

处理规则：

- 正式评分卡、summary、figure、report 首页一律不展示这些指标。
- 模型间比较一律不使用这些指标。
- 如研究者确实需要分析“是否复刻了原 binder”，可在单独的 exploratory analysis 中保留，但必须从主流程剥离。

## 7. 与四个模型内部评分的关系

四个模型内部自带的评分只能作为 metadata，不得直接进入官方比较主表。

统一口径如下：

- RFantibody：内部 `pLDDT / PAE / optional ddG` 只代表其自带复核，不代表官方 binding 结论。
- BoltzGen：内部 `design_to_target_iptm / design_ptm / PAE / refolded interface metrics` 只代表其内置筛选。
- Germinal：内部 `pDockQ2 / external iPTM / interface filters / AbMPNN score` 只代表其本地 accept-reject 系统。
- BindCraft：内部 `Average_i_pTM + PyRosetta interface metrics` 只代表其本地 ranking 逻辑。

统一原则：

- 各模型内部评分可保留到 `model_native_metadata` 表中。
- 官方比较只能看统一外部 judge 的结果。

## 8. 模型级比较的正式定义

后续模型比较不再问“谁更像 native”，而问以下问题：

1. 在相同样本、相同候选预算下，谁更容易产出可评估候选。
2. 在统一外部 judge 下，谁更容易产出通过 binding plausibility 的候选。
3. 在通过 binding plausibility 的候选中，谁的 affinity proxy 更强。
4. 在多 seed 或双 judge 下，谁的结果更稳健。
5. 在同等通过率下，谁能提供更多样性、更好的 developability。

建议的模型级主指标：

- `viability_rate`
- `bindability_pass_at_1`
- `bindability_pass_at_5`
- `bindability_pass_at_10`
- `strong_binder_proxy_rate`
- `robust_pass_rate`
- `median_binding_plausibility_metrics`
- `median_affinity_proxy_metrics`
- `accepted_diversity`
- `runtime_per_accepted_binder`

## 9. 阈值策略

阈值不直接硬编码为最终结论，而是分两步确定。

### 9.1 初始占位阈值

第一轮可先使用保守占位阈值，便于系统落地：

- `external_iptm >= 0.60`
- `external_ipae <= 10`
- `pdockq2 > 0.23`
- `interface_bsa >= 400 A^2`
- `interface_contact_pair_count >= 15`
- `interface_clash_ratio <= 0.05`
- `interface_hbond_count >= 2`
- `hotspot_contact_count >= 1`（若样本定义了热点）

这些阈值只作为系统上线前的起始门槛，不能直接作为正式论文结论。

### 9.2 校准阈值

正式阈值必须通过内部校准集获得：

- 正对照：真实 reference complex 的抗原-抗体配对。
- 近正对照：对 reference complex 做轻微扰动。
- 负对照：错位平移、随机旋转、打散接口、跨样本错配抗体序列。

校准目标：

- 主指标能稳定把“明显合理接口”和“明显坏接口”分开。
- 通过率对不同模型不存在明显结构偏置。
- 不因为某一模型输出风格不同而系统性吃亏。

## 10. 外部评估流水线设计

建议在 `evaluation/` 下按以下分层重建：

- `ingest/`
- `judge/`
- `metrics/`
- `report/`
- `pipeline/`

建议新增模块：

- `evaluation/ingest/build_external_eval_inputs.py`
- `evaluation/judge/run_primary_judge.py`
- `evaluation/judge/run_secondary_judge.py`
- `evaluation/metrics/compute_binding_plausibility.py`
- `evaluation/metrics/compute_affinity_proxy.py`
- `evaluation/metrics/compute_robustness.py`
- `evaluation/metrics/compute_developability.py`
- `evaluation/report/build_scorecards.py`
- `evaluation/report/build_model_comparison.py`
- `evaluation/pipeline/run_eval_pipeline.py`

### 10.1 新长表分层

建议新的官方长表分为四层：

- `external_eval_candidates.csv`
- `external_eval_seed_level.csv`
- `external_eval_sequence_level.csv`
- `external_eval_model_summary.csv`

### 10.2 推荐字段

`external_eval_candidates.csv` 建议包含：

- `candidate_id`
- `model`
- `sample_id`
- `candidate_rank_in_model`
- `source_sequence_path`
- `source_structure_path`
- `pred_sequence`
- `target_structure_path`
- `target_chain_id`
- `judge_name`
- `judge_seed`

`external_eval_seed_level.csv` 建议包含：

- `candidate_id`
- `judge_name`
- `judge_seed`
- `viability_pass`
- `binding_plausibility_pass`
- `affinity_proxy_pass`
- 全部 seed 级主指标字段

`external_eval_sequence_level.csv` 建议包含：

- `candidate_id`
- `viability_gate_status`
- `binding_pass_rate`
- `affinity_proxy_pass_rate`
- `robust_pass`
- 各主指标中位数
- developability 汇总字段

`external_eval_model_summary.csv` 建议包含：

- `model`
- `n_candidates`
- `viability_rate`
- `bindability_pass_at_1`
- `bindability_pass_at_5`
- `bindability_pass_at_10`
- `robust_pass_rate`
- `median_pdockq2`
- `median_external_iptm`
- `median_interface_bsa`
- `median_interface_dG`
- `median_surface_hydrophobicity`
- `accepted_diversity`

## 11. 报告层要求

报告首页必须按以下顺序展示：

1. Viability
2. Binding plausibility
3. Affinity proxy
4. Robustness
5. Developability

报告首页明确禁止出现：

- `TM-score ranking`
- `DockQ ranking`
- `CDR-H3 RMSD ranking`
- “与原抗体越像越好”类描述

建议图表：

- 各模型 `bindability_pass@K`
- 各模型 `robust_pass_rate`
- `external_iptm` 与 `pdockq2` 分布图
- `interface_bsa` 与 `interface_dG` 联合散点图
- `surface_hydrophobicity` 与 `binder_score` 分布图
- accepted candidate diversity heatmap

## 12. 迁移策略

旧体系到新体系按以下步骤迁移：

1. 先冻结现有 reference-based pipeline，不再扩展新图表和新 summary。
2. 在文档层明确标注旧体系为 legacy。
3. 新建统一外部评估流水线，不复用 `cdr_h3_rmsd / tm_score / dockq_score / irmsd` 作为主字段。
4. 老结果目录与新结果目录彻底分开，避免混用。

建议目录：

- 旧目录保留为：`outputs/evaluation/legacy_native_similarity/`
- 新目录使用：`outputs/evaluation/external_binding_benchmark/`

## 13. 分阶段实施计划

### Phase 1

- 用主 judge 重评 top-K 抗体序列。
- 完成 viability、binding plausibility、developability 三层。
- 完成模型级 pass@K 报告。

### Phase 2

- 增加 Rosetta 或等价能量层。
- 完成 affinity proxy。
- 增加 seed 级稳健性统计。

### Phase 3

- 接入第二 judge。
- 计算 judge agreement。
- 完成双 judge 稳健性报告。

### Phase 4

- 若后续获得实验数据，再做 external benchmark 与实验真值对齐。
- 再决定是否需要训练校准器或构建单一综合分。

## 14. 正式决策

本项目从现在开始采用如下正式决策：

- 删除“与原抗体像不像”在官方评估体系中的主地位。
- 不再用 `RMSD / TM-score / DockQ / iRMSD` 作为模型主排序依据。
- 统一使用外部 judge 作为四个模型的共同评委。
- 正式结果以多维评分卡输出，而不是单一总分。
- 模型内部自带评分一律降级为 metadata。

这份文档是后续重写 `evaluation/` 模块的正式设计依据。