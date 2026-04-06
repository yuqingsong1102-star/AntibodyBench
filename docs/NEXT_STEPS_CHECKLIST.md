# Next Steps Checklist

Use this as a low-stress execution list.
Only do 1-2 items per session.

## Day 1: Environment + Data Readiness

- [ ] Activate your working conda env (recommended: `abbench`)
- [ ] Run one-step data prep:
  - [ ] `python scripts/data_prep/prepare_native_inputs.py`

## Day 2: Build Ready Subset + First Metrics

- [ ] Confirm `data/prepared/dataset_index_ready.csv` was refreshed
- [ ] Run a small native smoke check (single sample):
  - [ ] `bash scripts/run.sh --model RFantibody --sample-id <sample_id> --max-samples 1`
- [ ] Verify unified output contract:
  - [ ] `outputs/native_predictions/<model>/<sample_id>/run_meta.json` exists
  - [ ] `outputs/native_predictions/<model>/<sample_id>/candidate_manifest.csv` exists
  - [ ] `outputs/native_predictions/<model>/manifest.csv` has one new row

## Day 3: Multi-Model Consistency

- [ ] Repeat Day 2 checks for:
  - [ ] Germinal
  - [ ] BoltzGen
  - [ ] BindCraft
- [ ] Confirm each model has:
  - [ ] per-sample files in `outputs/native_predictions/<model>/<sample_id>/`
  - [ ] rolling summary in `outputs/native_predictions/<model>/manifest.csv`

## Day 4: Review + Failure Triage

- [ ] Review model manifests and count failures:
  - [ ] check `status` / `error_summary` columns in `outputs/native_predictions/<model>/manifest.csv`
- [ ] Spot-check 3-5 failed samples:
  - [ ] open `run_stdout.log` / `run_stderr.log`
  - [ ] open `run_meta.json` for failure reason

## Troubleshooting Quick Notes

- If `cdr_h3_rmsd` is empty:
  - Check `cdr_h3_start/cdr_h3_end` in source CSV
  - Check chain IDs in `dataset_index` match reference complex structure
  - Check `reference_complex_path` points to full complex, not antigen-only structure

- If native input preparation fails:
  - Check `data/raw/dataset_index.csv` formatting and required columns
  - Check `data/raw/raw_pdbs/` and `data/raw/reference_complexes/` are readable
  - Check `data/prepared/epitopes/` is writable

## Definition of Done (v1)

- [ ] At least one model has valid non-empty `cdr_h3_rmsd` for a non-trivial sample set
- [ ] Same workflow runs reproducibly from `data/prepared/dataset_index_ready.csv`
- [ ] Project docs (`docs/`, `scripts/README.md`, `data/prepared/README.md`) are enough for rerun without memory

