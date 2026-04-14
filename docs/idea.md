# 抗体设计模型 Benchmark 实验设计思路

## 1. 核心问题：四个模型的产出逻辑完全不同

| 模型 | 架构 | 产出逻辑 | 内部过滤 | 默认产出 |
|------|------|---------|---------|---------|
| BoltzGen | Boltzmann 生成模型（单步 forward） | 生成 N 个 → 排序取 top BUDGET | 排序（无硬阈值） | 30 |
| mBER-open | AF + ColabDesign 迭代优化 | 逐条 trajectory → 实时双过滤(iPTM+pLDDT) → 凑够 num_accepted 停 | iPTM≥0.75, pLDDT≥0.7 | ≤100 |
| RFantibody | RFdiffusion → ProteinMPNN → RF2 | NUM_DESIGNS × NUM_SEQS = 全部产出，无内置过滤 | 无 | 10 |
| Germinal | AF2 hallucination → AbMPNN redesign | 逐条 trajectory → 早期质控 → 初始过滤 → AbMPNN 重设计(×4) → 14项严格终筛 → 凑够 max_passing 停 | 极严（14 项指标，含 Sc>0.6, pDockQ2>0.23, external_iPTM>0.74 等） | ≤100 |

### 关键差异

- **BoltzGen/RFantibody**：出得多但没筛选，候选质量参差不齐
- **mBER-open**：每步都含 AF 优化和评估，有中等过滤（iPTM+pLDDT），产出质量有保证
- **Germinal**：过滤最严格（14 项指标），出得少但质量最有保证

## 2. 四个模型默认参数（源码）

| 模型 | 关键参数 | 默认值 | 源码位置 |
|------|---------|--------|---------|
| **BoltzGen** | `--num_designs` | 10000 | `boltzgen/src/boltzgen/cli/boltzgen.py` |
| | `--budget`（最终保留数） | 30 | 同上 |
| | `top_budget`（先选 top 质量集） | 10 | `boltzgen/src/boltzgen/resources/config/filtering.yaml` |
| | `--inverse_fold_num_sequences` | 1 | `boltzgen/src/boltzgen/cli/boltzgen.py` |
| | `alpha`（多样性权重） | 0.001（非 peptide 协议） | `boltzgen/src/boltzgen/resources/config/filtering.yaml` |
| **mBER-open** | `num_accepted` | 100 | `mber-open/.../VHH_binder_design/cli.py` |
| | `max_trajectories` | 10000 | 同上 |
| | `min_iptm` | 0.75 | 同上 |
| | `min_plddt` | 0.70 | 同上 |
| **RFantibody** | `--num-designs`（RFdiffusion） | 10 | `RFantibody/src/rfantibody/cli/inference.py` |
| | `inference.num_designs`（底层配置） | 10 | `RFantibody/scripts/config/inference/base.yaml` |
| | `--seqs-per-struct`（ProteinMPNN） | 1 | `RFantibody/src/rfantibody/cli/inference.py` |
| | `--num-recycles`（RF2） | 10 | `RFantibody/src/rfantibody/cli/inference.py` |
| **Germinal** | `max_trajectories` | 10000 | `germinal/configs/config.yaml` |
| | `max_hallucinated_trajectories` | 1000 | 同上 |
| | `max_passing_designs`（接受上限） | 100 | 同上 |
| | `num_recycles_design`（VHH 配置） | 3 | `germinal/configs/run/vhh.yaml` |
| | `learning_rate`（VHH 配置） | 0.1 | 同上 |

**默认情况下**：BoltzGen 出 30，RFantibody 出 10，mBER/Germinal 最多出 100，各模型计算量差 100-1000 倍。默认参数没有直接可比性。

## 3. 公平性问题分析

### 问题

如果简单地让每个模型输出 100 个候选，然后比较平均质量：
- BoltzGen/RFantibody：100 个**未过滤**的候选
- mBER/Germinal：100 个**已通过严格内部筛选**的候选

比平均值不公平——内部有过滤的模型天然占优。

### 关键认知：内部过滤是模型设计的一部分

mBER 的 AF 迭代优化 + 实时过滤、Germinal 的 14 项终筛，这些不是"作弊"，而是模型 pipeline 的核心组成部分。去掉过滤就不是这个模型了。同理，BoltzGen/RFantibody 的快速大量生成也是它们的设计哲学。

**不要为了"公平"去关闭 mBER/Germinal 的内部筛选**——那就不是这个模型了。也**不要因为 BoltzGen/RFantibody 没有筛选就觉得不公平**——它们的设计理念就是快速大量生成。

## 4. 推荐的公平比较策略

### 核心原则

让每个模型按照作者推荐的方式运行，然后用**统一的外部评估**来打分。

### 实操方案

```
模型A 出了 100 个未过滤候选 → 统一外部评估 → 取 top-K
模型B 出了 50 个内部已过滤候选 → 统一外部评估 → 取 top-K
                                                ↓
                                      比较 top-K 的质量
```

### 推荐运行参数

| 模型 | 推荐参数 | 实际产出 | 每靶耗时 (1×GPU) | 理由 |
|---|---|---|---|---|
| BoltzGen | `NUM_DESIGNS=500, BUDGET=100` | 100 (排序取) | ~1.5h | 生成足够多再选 top100 |
| mBER-open | `num_accepted=100, max_traj=500` | ≤100 | ~4-40h | 用默认阈值，限制 trajectory 上限控时间 |
| RFantibody | `NUM_DESIGNS=25, NUM_SEQS=4` | 100 (全部) | ~4-8h | 25 结构×4 序列=100 |
| Germinal | 默认 `max_passing=100` | ≤100 | 数小时-数天 | 自带严格筛，通过率低 |

### 为什么这样公平

1. **外部评估是盲审**——评估 pipeline（Chai/AF3 预测 → DockQ 打分）不关心候选是怎么来的，只看最终结构质量
2. **top-K 归一化了产出量差异**——不管模型出了 100 还是 500 个，最终都只取 top-K（如 top-10）来比
3. **反映真实使用场景**——用户实际使用时就是：跑模型 → 选最好的几个 → 送去实验验证
4. **内部过滤的代价已体现在耗时中**——mBER 跑 2-4 天才出 100 个，BoltzGen 1.5 小时出 100 个，效率本身就是评估维度

## 5. 推荐的评估指标体系

### 5.1 质量指标（top-K 候选的外部评估）

- **最佳候选质量**：每个靶点 top-1 候选的 DockQ / iPTM / pLDDT
- **top-K 平均质量**：top-10 候选的平均 DockQ
- **通过率**：top-10 中 DockQ > 0.23 的比例

### 5.2 效率指标

- **每个靶点的 GPU 小时数**
- **每个通过质量阈值的候选的平均 GPU 耗时**（性价比）

### 5.3 鲁棒性指标

- **成功率**：在多少个靶点上至少有 1 个候选通过质量阈值
- **epitope 覆盖率**：设计的抗体是否真正覆盖了目标 epitope
- **序列多样性**：top-K 候选之间的序列聚类数

## 6. Benchmark 本质

> **给定一个靶点，各模型用合理配置跑完后，最终拿到的最好候选质量如何？**

- 每个模型用推荐参数跑出各自的候选
- 统一外部评估 + top-K 归一化
- 同时报告质量、效率、鲁棒性

模型有没有内部过滤只是 pipeline 特征之一，不影响可比性。过滤的代价（时间成本）和收益（质量提升）都会在最终结果中体现。

## 7. 代码层面：evaluation/ 与 analysis/ 的分工

### evaluation/ — 数据生产层（不需要改）

`evaluation/` 对**每个候选**逐一计算统一指标，产出 `eval_final.csv`。它不关心候选来自哪个模型、有没有经过内部筛选——这正是它应有的设计。

```
evaluation/ 的职责:
  Step 1: Chai-1 Judge (GPU)     → chai1_iptm, chai1_ptm, chai1_plddt_binder
  Step 2: Interface Metrics (CPU) → interface_bsa, pdockq2, clash/contact, approx_dG
  Step 3: Sequence Metrics (CPU)  → liabilities, charge, entropy, CDR 长度
  Aggregate                       → eval_merged.csv
  Step 4: Diversity + Epitope     → 聚类多样性、表位覆盖率
  Final Merge                     → eval_final.csv
```

### analysis/ — 数据分析层（公平性在这里解决）

`analysis/analyze_results.py` 消费 `eval_final.csv`，做统计和可视化。**公平性问题在这一层通过 top-K 归一化解决**：

```bash
# 不做 top-K（使用所有候选，可能不公平）
python analysis/analyze_results.py \
  --eval-dirs outputs/evaluation/results \
  --out-dir outputs/analysis/all_candidates

# 做 top-K 归一化（推荐：公平比较）
python analysis/analyze_results.py \
  --eval-dirs outputs/evaluation/results \
  --top-k 10 --rank-by chai1_iptm \
  --out-dir outputs/analysis/top10_comparison
```

### top-K 选择的工作方式

```
eval_final.csv (所有候选，每个候选一行)
      │
      ▼
按 (model, sample_id) 分组
      │
      ▼
组内按 --rank-by 指标降序排列
      │
      ▼
每组取前 K 个候选
      │
      ▼
在 top-K 子集上做所有图表和统计
```

### 支持的排序指标

| `--rank-by` 值 | 含义 |
|---|---|
| `chai1_iptm` | Chai-1 预测的 interface TM-score（默认） |
| `pdockq2` | 预测的 pDockQ2 对接质量 |
| `composite` | 0.5×iPTM + 0.5×pDockQ2 综合排序 |

### 建议的分析策略

1. **先跑一次 `--top-k 0`**（不限制），观察各模型的产出数量和原始分布
2. **再跑 `--top-k 10`**，用统一的 top-10 做公平比较
3. 两组结果放在一起看，既能体现模型间的公平质量差异，也能看到产出量和过滤策略的影响
