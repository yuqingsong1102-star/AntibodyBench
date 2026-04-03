# Native Run Issues Log

## 当前状态总览 (2026-04-03)

| 模型 | 状态 | 说明 |
|---|---|---|
| boltzgen | ✅ 成功 | 6 步全部完成，生成 2 个 VHH 设计 |
| BindCraft | ⏳ 待运行 | 配置已修复，等待 GPU |
| germinal | ❌ GPU OOM | 需 ~11 GB 连续显存，当前无可用 GPU |
| RFantibody | ⏳ 待运行 | 配置已优化，等待 GPU |

---

## 2026-04-02 run2 实际运行 (`sample_id=8q3j_B_A`)

Output root: `outputs/native_predictions_run2`

### 1) boltzgen — ✅ 成功

- **Status:** 成功完成全部 6 步 pipeline (design → inverse_folding → folding → design_folding → analysis → filtering)
- **结果:** 2 个 VHH 设计
  - Rank 1 (`design_spec_0`): 123 aa, quality_score=1.0, designfolding-rmsd=0.98Å
  - Rank 2 (`design_spec_1`): 137 aa, quality_score=0.0, designfolding-rmsd=2.22Å
- **输出位置:** `outputs/native_predictions_run2/boltzgen/8q3j_B_A/native_run/final_ranked_designs/`
- **已修复问题:**
  - env 名称 `boltzgen` → `bg`
  - 强制单 GPU 模式（`CUDA_VISIBLE_DEVICES=0`），避免 distributed TCPStore 失败
  - NATIVE_OUT_DIR 改为绝对路径（修复 runner cd 后相对路径错位问题）
- **注意:** collector (`collect_top1.py`) 需适配 boltzgen 的 CIF 输出格式

### 2) germinal — ❌ GPU OOM

- **Status:** 失败，JAX/XLA 在 AF design step 时 OOM
- **Key error:** `RESOURCE_EXHAUSTED: Out of memory while trying to allocate 11121629096 bytes` (~10.36 GB)
- **Evidence:** `outputs/native_predictions_run2/germinal/8q3j_B_A/run_stderr.log`
- **已修复问题:**
  - torch 升级到 2.6.0+cu124（CUDA 可用）
  - IgLM tokenizer 兼容性问题已修（transformers 5.3.0 + WordLevel tokenizer patch）
  - af_params_dir 指向 `models/BindCraft/params`
  - NATIVE_OUT_DIR 改为绝对路径
- **阻塞:** 需要一张 ≥11 GB 连续空闲显存的 GPU。当前最大空闲 GPU 4 仅 ~10.7 GB
- **Next action:** 等待 GPU 空闲，或考虑降低 AF 模型规模（fewer models / fewer recycles）

### 3) BindCraft — ⏳ 待运行

- **Status:** 配置已修复，尚未重新运行
- **已修复问题:**
  - `settings_target.json`: `number_of_final_designs: 10 → 1`
  - 创建 sample-local `advanced_settings.json`: `max_trajectories: 20`, `soft_iterations: 50`
  - `bindcraft.sh` 优先读取 sample 目录下的 `advanced_settings.json`
  - NATIVE_OUT_DIR 改为绝对路径
- **Next action:** 在 GPU 有空余时运行，预计需 ~6-10 分钟/trajectory

### 4) RFantibody — ⏳ 待运行

- **Status:** 配置已优化，尚未尝试运行
- **已修复问题:**
  - runner 自动选择空闲显存最多的 GPU
  - `rfantibody_config.env`: `NUM_DESIGNS=1, NUM_SEQS=1, NUM_RECYCLES=3`（降低显存需求）
  - NATIVE_OUT_DIR 改为绝对路径
- **Next action:** 在 GPU 有空余时运行

---

## 通用修复 (2026-04-02)

- **输出路径修复:** 所有 4 个 runner 的 `NATIVE_OUT_DIR` 从相对路径改为绝对路径解析。原因：runner 会 `cd` 到模型目录，导致相对路径写到错误位置
- **GPU 环境:** 服务器 8×RTX 4090 (24GB)，全部被其他用户进程占用大部分显存

---

## 历史记录

### 2026-04-01 首次 real run

Output root: `outputs/native_predictions_real`

所有 4 个模型均失败：
- **RFantibody:** CUDA OOM（未优化配置）
- **germinal:** torch 版本过低 → 升级后 torchvision 不匹配 → tokenizer 不兼容
- **BindCraft:** 运行时间过长，无限制 trajectory，未产出设计
- **boltzgen:** env 名称错误 + 分布式通信失败

详细修复过程见上方各模型"已修复问题"。
