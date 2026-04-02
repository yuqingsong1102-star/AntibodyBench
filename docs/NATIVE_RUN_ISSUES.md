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
