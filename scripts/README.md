# Scripts Index

Original scripts are kept at `scripts/` root for backward compatibility.

Use categorized entry points below to reduce confusion:

## Core pipeline

- `scripts/core/run.sh`
- `scripts/core/aggregate_metrics.py`

## Data preparation

- `scripts/data_prep/generate_dataset_index.py`
- `scripts/data_prep/fetch_reference_complexes.py`
- `scripts/data_prep/fill_cdr_h3_from_anarci.py`
- `scripts/data_prep/build_model_inputs.py`
- `scripts/data_prep/build_model_inputs_native.py` (templates aligned with original `models/*` entry formats)
- `scripts/data_prep/extract_epitopes_from_complexes.py`
- `scripts/data_prep/apply_epitopes_to_native_inputs.py`

## Ops

- `scripts/ops/build_apptainer_images.sh`

## Legacy / utility (still available)

- `scripts/prepare_antibody_dataset.py`
- `scripts/visualize_metrics.py`

