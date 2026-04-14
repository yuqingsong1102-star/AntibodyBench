# Evaluation — 四步外部评估架构

## 设计思路

### 定位

本模块是 AntibodyBench 的**统一外部评估层**——对四个 VHH 设计模型（BoltzGen、RFantibody、mBER-open、Germinal）的候选序列/结构做模型无关、靶点无关的质量打分。

核心原则：
- **盲审**：评估 pipeline 不知道候选来自哪个模型，仅根据序列和结构计算指标
- **连续值**：所有指标输出为连续数值，不在此层做 pass/fail 二元判定（硬阈值留给下游 `analysis/`）
- **可拆分**：四步独立运行，支持断点续跑，按步骤灵活调度 GPU/CPU 资源
- **可复现**：固定 seed、缓存中间结果、增量写入，崩溃后可无损恢复

### 四步架构与评估维度

```
候选序列/结构（来自 outputs/<model>/<sample_id>/candidate_manifest.csv）
     │
     ├── Step 1: Chai-1 Judge ─── AI 置信度维度（GPU, ~30-60s/候选）
     │     输入: 候选序列 + 抗原结构 → Chai-1 共折叠
     │     输出: chai1_results.csv
     │     指标: iPTM, pTM, pLDDT (binder)
     │     意义: 独立结构预测模型判断该序列是否能与靶点形成稳定复合物
     │
     ├── Step 2: Interface Metrics ─── 物理/几何维度（CPU, 可选 PyRosetta）
     │     输入: Chai-1 预测的复合物 CIF（或模型原生 PDB）
     │     输出: interface_results.csv
     │     指标: BSA, pDockQ2, contacts, clash, approx_dG, (rosetta_ddG, rosetta_sc)
     │     意义: 从三维结构中提取可量化的界面质量特征
     │
     ├── Step 3: Sequence Analysis ─── 可开发性维度（CPU, 纯序列）
     │     输入: 候选序列
     │     输出: sequence_results.csv
     │     指标: net_charge, liability motifs, CDR 长度, VHH hallmarks, GRAVY
     │     意义: 评估序列的可制造性和类药性（与结构无关的序列内在属性）
     │
     ├── Step 4: Diversity + Epitope ─── 多样性与表位特异性维度（CPU + MMseqs2）
     │     输入: eval_merged.csv + epitope JSON
     │     输出: diversity_epitope_results.csv
     │     指标: 聚类多样性（90% identity）、表位覆盖率
     │     意义: 衡量模型生成的候选是否多样、是否真的结合到目标表位
     │
     └── Aggregate
           输出: eval_merged.csv → eval_final.csv
           作用: 以 candidate_id 为主键，将四步指标合并为一行一候选的宽表
```

### 为什么这样设计

| 设计决策 | 理由 |
|---------|------|
| 选 Chai-1 而非 AlphaFold3 作为 judge | Chai-1 开源可本地运行，结果可复现；AF3 需 API 且有配额限制 |
| Step 2 默认 BioPython，可选 PyRosetta | BioPython 零依赖、跨平台；PyRosetta 需额外安装但提供 ddG/Sc 物理量 |
| 序列分析独立成步 | 不依赖结构预测，可在 GPU 步骤之前/之后任意运行，开销极低 |
| 聚类用 MMseqs2 | 比 CD-HIT 更快更准确，90% identity 是标准多样性阈值 |
| 全程输出连续值 | Top-K 选择和阈值判定留给 `analysis/` 层，evaluation 只负责"测量" |
| 以 candidate_id 为主键 | 所有步骤、所有 CSV 用同一主键，可任意 join |

### 三个评估维度的互补关系

```
                AI 置信度（Chai-1 iPTM/pLDDT）
               ╱                               ╲
              ╱   高 iPTM ≠ 物理上合理           ╲
             ╱    （AI 可能"幻觉"结合）           ╲
  物理/几何指标 ────────────────────────———————— 可开发性
  （BSA, pDockQ2, ddG）              （liabilities, charge）
        ↕                                    ↕
  结构上看起来结合 ≠ 能制造           序列"干净" ≠ 能结合
```

三个维度互相验证：
1. iPTM 高 + pDockQ2 高 + liability 少 → 最可靠的候选
2. iPTM 高但 pDockQ2 低 → Chai-1 结构预测可能不可靠
3. 物理指标好但序列 liability 多 → 结构可能好但不可制造

## 用法

```bash
# 全部四步
python -m evaluation.run_evaluation \
  --native-root outputs/pilot_boltzgen \
  --native-root outputs/pilot_rfantibody \
  --native-root outputs/pilot_mber \
  --native-root outputs/pilot_germinal \
  --index-csv data/prepared/dataset_index_ready.csv \
  --out-dir outputs/evaluation/results

# 仅 Step 1: Chai-1（需要 GPU，最耗时）
python -m evaluation.run_evaluation --step 1 \
  --native-root outputs/pilot_boltzgen \
  --device cuda:0

# 仅 Step 2 + Step 3（不需要 GPU，分钟级）
python -m evaluation.run_evaluation --step 2 --step 3 \
  --native-root outputs/pilot_boltzgen

# 启用 PyRosetta 增强 Step 2
python -m evaluation.run_evaluation --step 2 --use-pyrosetta \
  --native-root outputs/pilot_boltzgen

# 仅 Step 4: 多样性 + 表位（需先跑过 1-3）
python -m evaluation.run_evaluation --step 4 \
  --native-root outputs/pilot_boltzgen \
  --epitope-dir data/epitopes
```

## CLI 参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--native-root` | 模型输出目录（可重复，每个对应一个模型） | 必填 |
| `--index-csv` | 数据集索引 CSV（含 target 结构路径和链 ID） | `data/prepared/dataset_index_ready.csv` |
| `--out-dir` | 结果输出目录 | `outputs/evaluation/results` |
| `--step` | 只运行指定步骤（1/2/3/4，可重复） | 全部 |
| `--device` | Step 1 的 GPU 设备 | `cuda:0` |
| `--model` | 仅评估指定模型名（可重复） | 全部 |
| `--max-candidates` | 候选数上限（0=不限） | 0 |
| `--use-pyrosetta` | 启用 PyRosetta ddG/Sc 指标 | 关 |
| `--epitope-dir` | Step 4 表位 JSON 文件目录 | 自动检测 |

## 输出文件

| 文件 | Step | 内容 |
|------|------|------|
| `chai1_results.csv` | 1 | 每个候选的 Chai-1 iPTM, pTM, pLDDT |
| `interface_results.csv` | 2 | 界面 BSA, pDockQ2, contacts, approx_dG, (rosetta_*) |
| `sequence_results.csv` | 3 | 序列 charge, liability, CDR, GRAVY |
| `eval_merged.csv` | 1→3 | 合并表：一行一候选，包含 Step 1-3 全部指标 |
| `diversity_epitope_results.csv` | 4 | 聚类多样性、表位覆盖率 |
| `eval_final.csv` | 全部 | 最终宽表：eval_merged + Step 4 指标 |
| `chai1_runs/` | 1 | Chai-1 预测的 CIF 结构文件（可断点续跑） |

## 输入约定

### `--native-root` 目录结构

```
<native_root>/
  <model_name>/           # boltzgen / RFantibody / mber-open / germinal
    <sample_id>/          # 8dtn_A_B / 8hxq_A_B / ...
      candidate_manifest.csv   # 必须: candidate_id, status, structure_path, sequence_path
```

### `--index-csv` 必需列

| 列名 | 用途 |
|------|------|
| `sample_id` | 关联 native_root 中的样本目录 |
| `reference_complex_path` | target 抗原结构路径（Step 1 用） |
| `antigen_chain` | 抗原链 ID |
| `antibody_chain` | 抗体链 ID |
| `hotspot_string` | 表位残基列表（Step 2 hotspot 计数） |

## 各步骤关键指标速查

### Step 1 — AI 置信度

| 指标 | 含义 | 好坏方向 | 参考值 |
|------|------|---------|--------|
| `chai1_iptm` | 链间预测 TM-score | ↑ 越高越好 | >0.6 较好，>0.8 高置信 |
| `chai1_ptm` | 整体预测 TM-score | ↑ | >0.5 |
| `chai1_plddt_binder` | binder 链 pLDDT | ↑ | >0.7 |

### Step 2 — 物理/几何

| 指标 | 含义 | 好坏方向 |
|------|------|---------|
| `pdockq2` | 基于 pLDDT 和接触数的对接质量 | ↑ >0.23 acceptable |
| `interface_bsa` | 埋入表面积 (Å²) | ↑ >600 较好 |
| `interface_clash_ratio` | 碰撞比例 | ↓ 越低越好 |
| `approx_dG` | 近似结合能 (kcal/mol) | ↓ 越负越好 |
| `hotspot_contact_count` | 接触的表位残基数 | ↑ |
| `rosetta_ddG` | Rosetta ddG（可选） | ↓ |

### Step 3 — 可开发性

| 指标 | 含义 | 好坏方向 |
|------|------|---------|
| `total_liability_count` | 序列风险位点总数 | ↓ |
| `net_charge` | 净电荷 | 接近 0 |
| `has_vhh_hallmarks` | 是否具有 VHH 标志残基 | 1 |
| `gravy` | 疏水性指数 | 适中 (-0.5~0.5) |

### Step 4 — 多样性与表位

| 指标 | 含义 | 好坏方向 |
|------|------|---------|
| `group_cluster_count` | 聚类数 (90% identity) | ↑ 越多越多样 |
| `group_singleton_frac` | 单独聚类比例 | ↑ 越高越多样 |
| `epitope_coverage` | 覆盖的表位残基比例 | ↑ |

## 代码结构

```
evaluation/
├── __init__.py
├── README.md              ← 本文件（设计思路 + 用法）
├── README_v2.md           ← 详细版文档（实现细节 + 排错）
├── run_evaluation.py      — 统一入口（CLI + 四步调度）
├── aggregate.py           — 合并多步 CSV（eval_merged / eval_final）
├── schema.py              — 外部评估 CSV 字段常量定义
├── structure_inference.py — BioPython 结构解析 + 链推断工具
└── steps/
    ├── __init__.py
    ├── chai1_judge.py      — Step 1: Chai-1 共折叠评估
    ├── interface_metrics.py — Step 2: 界面几何 + 可选 PyRosetta
    ├── sequence_metrics.py  — Step 3: 序列理化 + 可开发性
    └── diversity_epitope.py — Step 4: MMseqs2 聚类 + 表位覆盖

~2160 行代码（不含 schema）
```

## 与 analysis/ 的分工

| 层 | 职责 | 输入 | 输出 |
|----|------|------|------|
| `evaluation/` | 计算指标（数据生产） | 模型候选 + index CSV | eval_final.csv |
| `analysis/` | 统计图表（数据消费） | eval_final.csv | 论文图片 + 汇总表 |

`evaluation/` 只负责"测量"，不做排序和比较判断；`analysis/` 负责 top-K 选择、跨模型对比、论文可视化。
