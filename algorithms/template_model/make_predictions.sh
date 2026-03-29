#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR=""
OUTPUT_DIR=""
MODEL_NAME="template_model"

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
shopt -s nullglob
for sample_json in "${INPUT_DIR}"/*.json; do
  sample_id="$(basename "${sample_json}" .json)"
  input_pdb="$(
    SAMPLE_JSON="${sample_json}" python3 - <<'PY'
import json, os
with open(os.environ["SAMPLE_JSON"], "r", encoding="utf-8") as f:
  d = json.load(f)
print((d.get("reference_structure_path") or "").strip(), end="")
PY
  )"
  if [[ -n "${input_pdb}" && -f "${input_pdb}" ]]; then
    cp "${input_pdb}" "${OUTPUT_DIR}/${sample_id}.pdb"
    echo "[OK] ${MODEL_NAME} 产物: ${OUTPUT_DIR}/${sample_id}.pdb"
  else
    echo "[WARN] ${sample_id}: 输入PDB缺失，跳过"
  fi
done
