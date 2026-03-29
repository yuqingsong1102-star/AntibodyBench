#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ALG_DIR="${ROOT_DIR}/algorithms"

for d in "${ALG_DIR}"/*; do
  [[ -d "${d}" ]] || continue
  def_file="${d}/container.def"
  sif_file="${d}/container.sif"
  if [[ ! -f "${def_file}" ]]; then
    continue
  fi
  echo "[INFO] build: ${def_file}"
  apptainer build "${sif_file}" "${def_file}"
done
