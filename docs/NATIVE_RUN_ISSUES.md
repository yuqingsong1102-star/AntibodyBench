# Native Run Issues Log

## 2026-04-01 one-sample real run (`sample_id=8q3j_B_A`)

Output root: `outputs/native_predictions_real`

### 1) RFantibody

- **Status:** failed
- **Observed:** runner starts, model loads, then no `top1_structure`
- **Key error:** CUDA OOM while loading/running RFdiffusion checkpoint
- **Evidence:** `outputs/native_predictions_real/RFantibody/8q3j_B_A/run_stderr.log`
- **Next action:**
  - run on freer GPU (`CUDA_VISIBLE_DEVICES=6` or `7`)
  - reduce diffusion load (lower `--diffuser-t`, smaller design count)
  - avoid running with other heavy GPU jobs in parallel

### 2) germinal

- **Status:** failed
- **Observed:** pipeline starts, then exits before accepted structures
- **Key error:** transformers blocks `torch.load` because torch < 2.6
- **Evidence:** `outputs/native_predictions_real/germinal/8q3j_B_A/run_stderr.log`
- **Next action:**
  - upgrade torch in `germinal` env to `>=2.6`
  - if network unstable, retry download or use local wheel/cache mirror
  - proxy test note (2026-04-01): `http://1user:pass1@10.70.106.147:7893` returned `Connection refused` from this host
  - after torch upgrade to 2.6.0, new blocker: `torchvision::nms does not exist` (torch/torchvision mismatch), need aligned torchvision build
  - after torchvision fix, new blocker is GPU OOM in JAX step; needs lower GPU pressure and explicit GPU selection
  - after cleaning duplicate torch metadata, `transformers` version gate is fixed (`is_torch_greater_or_equal('2.6') == True`)
  - latest blocker now: `ValueError: Unrecognized amino acid token: A` in `colabdesign/iglm/model.py` during `CustomIgLM` init

### 3) BindCraft

- **Status:** failed (no final structure produced in this run)
- **Observed:** process runs very long; no final accepted structure detected by collector
- **Evidence:** `outputs/native_predictions_real/BindCraft/8q3j_B_A/top1_meta.json`
- **Next action:**
  - run a quick config first (`number_of_final_designs=1`, lower iterations)
  - confirm expected output files in `native_run/Accepted` or stats CSV
  - if run is intentionally long, set longer wait window before collect

### 4) boltzgen

- **Status:** failed (collector did not find final structure in this run window)
- **Observed:** run creates config files and `design_spec.cif`, but no ranked outputs found before stop
- **Evidence:** `outputs/native_predictions_real/boltzgen/8q3j_B_A/top1_meta.json`
- **Next action:**
  - use existing env `bg` (`BOLTZGEN_CONDA_ENV=bg`)
  - allow full run to finish (longer wait) before collecting Top-1
  - verify `native_run/final_ranked_designs/` appears

## Notes

- Current GPU nodes are heavily occupied by other workloads.
- For stable validation, prefer running models one-by-one during lower GPU contention windows.

---

## 2026-04-02 Fixes Applied

### 1) RFantibody
- **Fix:** runner now auto-selects GPU with most free memory (`nvidia-smi` query) when `CUDA_VISIBLE_DEVICES` not set
- **Fix:** reduced `rfantibody_config.env`: `NUM_DESIGNS=1`, `NUM_SEQS=1`, `NUM_RECYCLES=3` to lower VRAM footprint
- **Remaining risk:** still needs a GPU with ≥4 GB free; all 8 GPUs currently under heavy load

### 2) germinal
- **Fix:** reinstalled `torch==2.6.0+cu124` with CUDA support (`torch.cuda.is_available() == True` confirmed)
- **Fix:** patched `iglm/model/IgLM.py` — replaced `BertTokenizerFast` with `PreTrainedTokenizerFast` built from `tokenizers.models.WordLevel`, because `transformers 5.3.0` silently drops non-special tokens from simple vocab files
- **Verified:** all 20 standard amino acids now correctly tokenized (A→5, R→19, ... V→22)
- **Next risk:** GPU OOM in JAX step may still appear under GPU contention

### 3) BindCraft
- **Fix:** `settings_target.json` changed `number_of_final_designs: 10 → 1`
- **Fix:** created sample-local `advanced_settings.json` with `max_trajectories: 20` (was `false`/unlimited) and `soft_iterations: 50` (was 75)
- **Fix:** `bindcraft.sh` now checks for `${SAMPLE_INPUT_DIR}/advanced_settings.json` before falling back to global default

### 4) boltzgen
- **Fix:** `boltzgen.sh` default `BOLTZGEN_CONDA_ENV` changed from `boltzgen` to `bg` (matching actual env name)
- **Fix:** runner now exports `CUDA_VISIBLE_DEVICES=0` by default to force single-GPU mode, avoiding distributed TCPStore failures
