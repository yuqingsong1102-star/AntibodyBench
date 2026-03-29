#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR=""
OUTPUT_DIR=""
MODEL_NAME="RFantibody"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir) INPUT_DIR="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --model-name) MODEL_NAME="$2"; shift 2 ;;
    *) echo "[ERROR] 未知参数: $1"; exit 2 ;;
  esac
done

if [[ -z "${INPUT_DIR}" || -z "${OUTPUT_DIR}" ]]; then
  echo "[ERROR] 缺少 --input-dir 或 --output-dir"
  exit 2
fi

mkdir -p "${OUTPUT_DIR}"
MODEL_WORKDIR="${RFANTIBODY_WORKDIR:-/home/yqsong/projects/antibody_benchmark/models/RFantibody}"
CMD_TEMPLATE="${RFANTIBODY_CMD:-}"
DRY_RUN="${RFANTIBODY_DRY_RUN:-0}"

pick_first_structure() {
  local search_dir="$1"
  SEARCH_DIR="${search_dir}" python3 - <<'PY'
import os
from pathlib import Path
root = Path(os.environ["SEARCH_DIR"])
files = sorted(root.rglob("*.pdb")) + sorted(root.rglob("*.cif"))
print(str(files[0]) if files else "", end="")
PY
}

shopt -s nullglob
for sample_json in "${INPUT_DIR}"/*.json; do
  sample_id="$(basename "${sample_json}" .json)"
  sample_out="${OUTPUT_DIR}/${sample_id}"
  mkdir -p "${sample_out}"
  input_pdb="$(
    SAMPLE_JSON="${sample_json}" python3 - <<'PY'
import json, os
with open(os.environ["SAMPLE_JSON"], "r", encoding="utf-8") as f:
  d = json.load(f)
print((d.get("reference_structure_path") or "").strip(), end="")
PY
  )"
  if [[ -z "${input_pdb}" || ! -f "${input_pdb}" ]]; then
    echo "[WARN] ${sample_id}: 输入PDB缺失，跳过"
    continue
  fi

  if [[ "${DRY_RUN}" == "1" || -z "${CMD_TEMPLATE}" ]]; then
    cp "${input_pdb}" "${OUTPUT_DIR}/${sample_id}.pdb"
    echo "[OK] ${MODEL_NAME} (dry-run) 产物: ${OUTPUT_DIR}/${sample_id}.pdb"
    continue
  fi

  cmd="${CMD_TEMPLATE//__INPUT__/${input_pdb}}"
  cmd="${cmd//__OUTDIR__/${sample_out}}"
  cmd="${cmd//__SAMPLE__/${sample_id}}"
  if (cd "${MODEL_WORKDIR}" && bash -lc "${cmd}") >"${sample_out}/stdout.log" 2>"${sample_out}/stderr.log"; then
    pred_file="$(pick_first_structure "${sample_out}")"
    if [[ -n "${pred_file}" ]]; then
      ext="${pred_file##*.}"
      cp "${pred_file}" "${OUTPUT_DIR}/${sample_id}.${ext}"
      echo "[OK] ${MODEL_NAME} 产物: ${OUTPUT_DIR}/${sample_id}.${ext}"
    else
      echo "[WARN] ${sample_id}: 命令执行成功但未找到pdb/cif产物"
    fi
  else
    echo "[WARN] ${sample_id}: 执行失败，见 ${sample_out}/stderr.log"
  fi
done
