#!/usr/bin/env bash
set -euo pipefail

SAMPLE_ID=""
SAMPLE_INPUT_DIR=""
SAMPLE_OUTPUT_DIR=""
PROJECT_ROOT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --sample-input-dir) SAMPLE_INPUT_DIR="$2"; shift 2 ;;
    --sample-output-dir) SAMPLE_OUTPUT_DIR="$2"; shift 2 ;;
    --project-root) PROJECT_ROOT="$2"; shift 2 ;;
    *) echo "[ERROR] 未知参数: $1"; exit 2 ;;
  esac
done

if [[ -z "${SAMPLE_ID}" || -z "${SAMPLE_INPUT_DIR}" || -z "${SAMPLE_OUTPUT_DIR}" ]]; then
  echo "[ERROR] 缺少必填参数"
  exit 2
fi

SETTINGS_FILE="${SAMPLE_INPUT_DIR}/settings_target.json"
if [[ ! -f "${SETTINGS_FILE}" ]]; then
  echo "[ERROR] 缺少输入文件: ${SETTINGS_FILE}"
  exit 3
fi

MODEL_WORKDIR="${BINDCRAFT_WORKDIR:-/home/yqsong/projects/antibody_benchmark/models/BindCraft}"
CONDA_ENV="${BINDCRAFT_CONDA_ENV:-BindCraft}"
USE_CONDA_RUN="${BINDCRAFT_USE_CONDA_RUN:-1}"
FILTERS_JSON="${BINDCRAFT_FILTERS:-${MODEL_WORKDIR}/settings_filters/default_filters.json}"
# Prefer sample-local advanced settings if present
if [[ -f "${SAMPLE_INPUT_DIR}/advanced_settings.json" ]]; then
  ADVANCED_JSON="${SAMPLE_INPUT_DIR}/advanced_settings.json"
else
  ADVANCED_JSON="${BINDCRAFT_ADVANCED:-${MODEL_WORKDIR}/settings_advanced/default_4stage_multimer.json}"
fi
FINAL_DESIGNS_OVERRIDE="${BINDCRAFT_NUM_FINAL_DESIGNS:-}"
MAX_TRAJECTORIES_OVERRIDE="${BINDCRAFT_MAX_TRAJECTORIES:-}"
NATIVE_OUT_DIR="$(cd "$(dirname "${SAMPLE_OUTPUT_DIR}")" && pwd)/$(basename "${SAMPLE_OUTPUT_DIR}")/native_run"
mkdir -p "${NATIVE_OUT_DIR}"
DESIGN_ROOT="${BINDCRAFT_DESIGN_PATH:-${NATIVE_OUT_DIR}/designs}"
RUNTIME_INPUT_DIR="${NATIVE_OUT_DIR}/runtime_inputs"
RUNTIME_SETTINGS_FILE="${RUNTIME_INPUT_DIR}/settings_target.runtime.json"
RUNTIME_ADVANCED_FILE="${RUNTIME_INPUT_DIR}/advanced_settings.runtime.json"
mkdir -p "${RUNTIME_INPUT_DIR}"

python3 - <<'PY' "${SETTINGS_FILE}" "${RUNTIME_SETTINGS_FILE}" "${DESIGN_ROOT}" "${FINAL_DESIGNS_OVERRIDE}"
import json
import sys
from pathlib import Path

src = Path(sys.argv[1])
dst = Path(sys.argv[2])
design_root = sys.argv[3]
final_designs_raw = sys.argv[4].strip()
payload = json.loads(src.read_text(encoding="utf-8"))
payload["design_path"] = design_root
if final_designs_raw:
    payload["number_of_final_designs"] = int(final_designs_raw)
dst.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
PY

python3 - <<'PY' "${ADVANCED_JSON}" "${RUNTIME_ADVANCED_FILE}" "${MAX_TRAJECTORIES_OVERRIDE}"
import json
import sys
from pathlib import Path

src = Path(sys.argv[1])
dst = Path(sys.argv[2])
max_trajectories_raw = sys.argv[3].strip()
payload = json.loads(src.read_text(encoding="utf-8"))
if max_trajectories_raw:
    payload["max_trajectories"] = int(max_trajectories_raw)
dst.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
PY

default_cmd="python bindcraft.py --settings \"${RUNTIME_SETTINGS_FILE}\" --filters \"${FILTERS_JSON}\" --advanced \"${RUNTIME_ADVANCED_FILE}\""
CMD_TEMPLATE="${BINDCRAFT_CMD:-${default_cmd}}"
cmd="${CMD_TEMPLATE//__SAMPLE_ID__/${SAMPLE_ID}}"
cmd="${cmd//__INPUT_DIR__/${SAMPLE_INPUT_DIR}}"
cmd="${cmd//__OUTPUT_DIR__/${NATIVE_OUT_DIR}}"
cmd="${cmd//__SETTINGS_FILE__/${RUNTIME_SETTINGS_FILE}}"
cmd="${cmd//__FILTERS_FILE__/${FILTERS_JSON}}"
cmd="${cmd//__ADVANCED_FILE__/${RUNTIME_ADVANCED_FILE}}"

if [[ ! -d "${MODEL_WORKDIR}" ]]; then
  echo "[ERROR] BindCraft 工作目录不存在: ${MODEL_WORKDIR}"
  exit 4
fi

if [[ "${USE_CONDA_RUN}" == "1" ]] && command -v conda >/dev/null 2>&1; then
  (
    cd "${MODEL_WORKDIR}"
    conda run -n "${CONDA_ENV}" bash -lc "${cmd}"
  )
else
  (
    cd "${MODEL_WORKDIR}"
    bash -lc "${cmd}"
  )
fi

echo "[OK] BindCraft native runner 完成: ${SAMPLE_ID}"
