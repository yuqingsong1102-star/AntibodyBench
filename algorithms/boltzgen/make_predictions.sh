#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR=""
OUTPUT_DIR=""
MODEL_NAME="boltzgen"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir)
      INPUT_DIR="$2"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --model-name)
      MODEL_NAME="$2"
      shift 2
      ;;
    *)
      echo "[ERROR] 未知参数: $1"
      exit 2
      ;;
  esac
done

if [[ -z "${INPUT_DIR}" || -z "${OUTPUT_DIR}" ]]; then
  echo "[ERROR] 缺少 --input-dir 或 --output-dir"
  exit 2
fi

mkdir -p "${OUTPUT_DIR}"
CMD_TEMPLATE="${BOLTZGEN_CMD:-}"
MODEL_WORKDIR="${BOLTZGEN_WORKDIR:-/home/yqsong/projects/antibody_benchmark/models/boltzgen}"
DRY_RUN="${BOLTZGEN_DRY_RUN:-0}"
MAX_SAMPLES="${BOLTZGEN_MAX_SAMPLES:-0}"
NUM_DESIGNS="${BOLTZGEN_NUM_DESIGNS:-2}"
BUDGET="${BOLTZGEN_BUDGET:-1}"
PROTOCOL="${BOLTZGEN_PROTOCOL:-protein-anything}"
BINDER_LEN="${BOLTZGEN_BINDER_LEN:-60..90}"
USE_CONDA_RUN="${BOLTZGEN_USE_CONDA_RUN:-1}"
CONDA_ENV="${BOLTZGEN_CONDA_ENV:-boltzgen}"
HAS_CONDA=0

if [[ "${USE_CONDA_RUN}" == "1" ]]; then
  if command -v conda >/dev/null 2>&1; then
    HAS_CONDA=1
    echo "[INFO] boltzgen 将使用 conda 环境: ${CONDA_ENV}"
  else
    echo "[WARN] 未找到 conda，boltzgen 将在当前环境执行。"
  fi
fi

run_model_cmd() {
  local cmd="$1"
  if [[ "${USE_CONDA_RUN}" == "1" && "${HAS_CONDA}" == "1" ]]; then
    (cd "${MODEL_WORKDIR}" && conda run -n "${CONDA_ENV}" bash -lc "${cmd}")
  else
    (cd "${MODEL_WORKDIR}" && bash -lc "${cmd}")
  fi
}

if [[ "${DRY_RUN}" != "1" && -z "${CMD_TEMPLATE}" ]]; then
  CMD_TEMPLATE="boltzgen run '__SPEC__' --output '__OUTDIR__' --protocol '__PROTOCOL__' --num_designs __NUM_DESIGNS__ --budget __BUDGET__"
  echo "[INFO] 未设置 BOLTZGEN_CMD，使用内置默认命令模板。"
fi

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
processed_count=0
for sample_json in "${INPUT_DIR}"/*.json; do
  if [[ "${MAX_SAMPLES}" != "0" && "${processed_count}" -ge "${MAX_SAMPLES}" ]]; then
    echo "[INFO] 达到 BOLTZGEN_MAX_SAMPLES=${MAX_SAMPLES}，停止继续样本。"
    break
  fi
  processed_count=$((processed_count + 1))
  sample_id="$(basename "${sample_json}" .json)"
  sample_out="${OUTPUT_DIR}/${sample_id}"
  mkdir -p "${sample_out}"

  input_pdb="$(
    SAMPLE_JSON="${sample_json}" python3 - <<'PY'
import json
import os
with open(os.environ["SAMPLE_JSON"], "r", encoding="utf-8") as f:
  d = json.load(f)
print(str(d.get("reference_structure_path", "")).strip(), end="")
PY
  )"

  if [[ -z "${input_pdb}" || ! -f "${input_pdb}" ]]; then
    echo "[WARN] ${sample_id}: 输入PDB缺失，跳过"
    continue
  fi

  target_chain="$(
    SAMPLE_JSON="${sample_json}" python3 - <<'PY'
import json
import os
with open(os.environ["SAMPLE_JSON"], "r", encoding="utf-8") as f:
  d = json.load(f)
chain = str(d.get("antigen_chain", "")).strip()
if "," in chain:
  chain = chain.split(",", 1)[0].strip()
print(chain or "A", end="")
PY
  )"

  design_yaml="${sample_out}/design_spec.yaml"
  INPUT_PDB="${input_pdb}" TARGET_CHAIN="${target_chain}" BINDER_LEN="${BINDER_LEN}" DESIGN_YAML="${design_yaml}" python3 - <<'PY'
import os
from pathlib import Path

input_pdb = Path(os.environ["INPUT_PDB"]).resolve()
target_chain = os.environ["TARGET_CHAIN"] or "A"
binder_len = os.environ["BINDER_LEN"] or "60..90"
design_yaml = Path(os.environ["DESIGN_YAML"])
design_yaml.write_text(
  "\n".join(
    [
      "entities:",
      "  - protein:",
      "      id: B",
      f"      sequence: {binder_len}",
      "  - file:",
      f"      path: {input_pdb}",
      "      include:",
      "        - chain:",
      f"            id: {target_chain}",
      "",
    ]
  ),
  encoding="utf-8",
)
PY

  if [[ "${DRY_RUN}" == "1" || -z "${CMD_TEMPLATE}" ]]; then
    cp "${input_pdb}" "${OUTPUT_DIR}/${sample_id}.pdb"
    echo "[OK] ${MODEL_NAME} (dry-run) 产物: ${OUTPUT_DIR}/${sample_id}.pdb"
    continue
  fi

  cmd="${CMD_TEMPLATE//__INPUT__/${input_pdb}}"
  cmd="${cmd//__OUTDIR__/${sample_out}}"
  cmd="${cmd//__SAMPLE__/${sample_id}}"
  cmd="${cmd//__SPEC__/${design_yaml}}"
  cmd="${cmd//__PROTOCOL__/${PROTOCOL}}"
  cmd="${cmd//__NUM_DESIGNS__/${NUM_DESIGNS}}"
  cmd="${cmd//__BUDGET__/${BUDGET}}"

  if run_model_cmd "${cmd}" >"${sample_out}/stdout.log" 2>"${sample_out}/stderr.log"; then
    pred_file="$(pick_first_structure "${sample_out}")"
    if [[ -z "${pred_file}" ]]; then
      echo "[WARN] ${sample_id}: 命令执行成功但未找到pdb/cif产物"
      continue
    fi
    ext="${pred_file##*.}"
    cp "${pred_file}" "${OUTPUT_DIR}/${sample_id}.${ext}"
    echo "[OK] ${MODEL_NAME} 产物: ${OUTPUT_DIR}/${sample_id}.${ext}"
  else
    echo "[WARN] ${sample_id}: 执行失败，见 ${sample_out}/stderr.log"
  fi
done
