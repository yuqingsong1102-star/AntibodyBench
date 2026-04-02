# AntibodyBench Project Layout

This document is a practical map for daily work.

## Top-Level Roles

- `algorithms/`: model adapters (preprocess, prediction, postprocess per model)
- `data/`: raw/reference data (do not edit casually)
- `inputs/`: benchmark index and model input metadata
- `outputs/`: generated artifacts (input/prediction/evaluation)
- `evaluation_tools/`: metrics implementation
- `scripts/`: runnable entry scripts
- `docs/`: quick project maps and operation notes

## Recommended Daily Entry Points

Use these entry points first:

- `scripts/run.sh`
- `scripts/data_prep/generate_dataset_index.py`
- `scripts/data_prep/fetch_reference_complexes.py`
- `scripts/data_prep/fill_cdr_h3_from_anarci.py`
- `scripts/ops/build_apptainer_images.sh`

Notes:

- `scripts/data_prep/*` is the canonical data-prep entry namespace.
- Data-prep commands should always use `scripts/data_prep/*` paths.

See:

- `scripts/README.md`
- `inputs/README.md`
- `inputs/antibody_datasets/README.md`
- `outputs/README.md`

## Key CSV Files

- `inputs/antibody_datasets/dataset_index.csv`: canonical index
- `inputs/antibody_datasets/dataset_index_h3_annotated.csv`: H3 annotation output
- `inputs/antibody_datasets/dataset_index_ready.csv`: filtered ready subset

## Keep vs Ignore

- Keep: `inputs/`, `scripts/`, `algorithms/`, `evaluation_tools/`, `data/raw/reference_complexes/`
- Generated only: `outputs/`
- Optional/local-only: local dependency folders, caches, pycache files

See also: `outputs/README.md` for result navigation.

