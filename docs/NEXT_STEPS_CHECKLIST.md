# Next Steps Checklist

Use this as a low-stress execution list.
Only do 1-2 items per session.

## Day 1: Environment + Data Readiness

- [ ] Activate your working conda env (recommended: `abbench`)
- [ ] Ensure ANARCI dependency chain is ready:
  - [ ] `python -c "import anarci; print('anarci ok')"`
  - [ ] `which hmmscan`
- [ ] Refresh full complex paths:
  - [ ] `python scripts/data_prep/fetch_reference_complexes.py --input-csv inputs/antibody_datasets/dataset_index.csv --output-csv inputs/antibody_datasets/dataset_index.csv`
- [ ] Re-run H3 annotation:
  - [ ] `python scripts/data_prep/fill_cdr_h3_from_anarci.py --input-csv inputs/antibody_datasets/dataset_index.csv --output-csv inputs/antibody_datasets/dataset_index_h3_annotated.csv`

## Day 2: Build Ready Subset + First Metrics

- [ ] Create/update ready subset (`dataset_index_ready.csv`) from rows with `cdr_h3_status=ok`
- [ ] Run a small native smoke check (single sample):
  - [ ] `bash scripts/run.sh --model RFantibody --sample-id <sample_id> --max-samples 1`
- [ ] Verify unified output contract:
  - [ ] `outputs/native_predictions/<model>/<sample_id>/top1_meta.json` exists
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
  - [ ] open `top1_meta.json` for selection rule and failure reason

## Troubleshooting Quick Notes

- If `cdr_h3_rmsd` is empty:
  - Check `cdr_h3_start/cdr_h3_end` in source CSV
  - Check chain IDs in `dataset_index` match reference complex structure
  - Check `reference_complex_path` points to full complex, not antigen-only structure

- If ANARCI fails:
  - Confirm `hmmscan` exists in PATH
  - Re-run with a tiny CSV subset first

## Definition of Done (v1)

- [ ] At least one model has valid non-empty `cdr_h3_rmsd` for a non-trivial sample set
- [ ] Same workflow runs reproducibly from `dataset_index_ready.csv`
- [ ] Project docs (`docs/`, `scripts/README.md`, `inputs/README.md`, `outputs/README.md`) are enough for rerun without memory

