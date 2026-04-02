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
ADVANCED_JSON="${BINDCRAFT_ADVANCED:-${MODEL_WORKDIR}/settings_advanced/default_4stage_multimer.json}"
NATIVE_OUT_DIR="${SAMPLE_OUTPUT_DIR}/native_run"
mkdir -p "${NATIVE_OUT_DIR}"

default_cmd="python bindcraft.py --settings \"${SETTINGS_FILE}\" --filters \"${FILTERS_JSON}\" --advanced \"${ADVANCED_JSON}\""
CMD_TEMPLATE="${BINDCRAFT_CMD:-${default_cmd}}"
cmd="${CMD_TEMPLATE//__SAMPLE_ID__/${SAMPLE_ID}}"
cmd="${cmd//__INPUT_DIR__/${SAMPLE_INPUT_DIR}}"
cmd="${cmd//__OUTPUT_DIR__/${NATIVE_OUT_DIR}}"
cmd="${cmd//__SETTINGS_FILE__/${SETTINGS_FILE}}"
cmd="${cmd//__FILTERS_FILE__/${FILTERS_JSON}}"
cmd="${cmd//__ADVANCED_FILE__/${ADVANCED_JSON}}"

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
