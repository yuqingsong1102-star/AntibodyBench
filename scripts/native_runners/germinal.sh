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

TARGET_YAML="${SAMPLE_INPUT_DIR}/target.yaml"
OVERRIDES_FILE="${SAMPLE_INPUT_DIR}/run_overrides.txt"
if [[ ! -f "${TARGET_YAML}" ]]; then
  echo "[ERROR] 缺少输入文件: ${TARGET_YAML}"
  exit 3
fi
if [[ ! -f "${OVERRIDES_FILE}" ]]; then
  echo "[ERROR] 缺少输入文件: ${OVERRIDES_FILE}"
  exit 3
fi

MODEL_WORKDIR="${GERMINAL_WORKDIR:-/home/yqsong/projects/antibody_benchmark/models/germinal}"
CONDA_ENV="${GERMINAL_CONDA_ENV:-germinal}"
USE_CONDA_RUN="${GERMINAL_USE_CONDA_RUN:-1}"
NATIVE_OUT_DIR="$(cd "$(dirname "${SAMPLE_OUTPUT_DIR}")" && pwd)/$(basename "${SAMPLE_OUTPUT_DIR}")/native_run"
mkdir -p "${NATIVE_OUT_DIR}"

default_cmd="python run_germinal.py __OVERRIDES__ results_dir=\"__OUTPUT_DIR__\" experiment_name=\"__SAMPLE_ID__\""
CMD_TEMPLATE="${GERMINAL_CMD:-${default_cmd}}"

if [[ "${CMD_TEMPLATE}" == *"__OVERRIDES__"* ]]; then
  OVERRIDES_CONTENT="$(tr '\n' ' ' < "${OVERRIDES_FILE}")"
  # Keep quotes around hotspot strings for Hydra parsing.
  OVERRIDES_CONTENT="${OVERRIDES_CONTENT//\'/\\\'}"
  CMD_TEMPLATE="${CMD_TEMPLATE//__OVERRIDES__/${OVERRIDES_CONTENT}}"
fi
cmd="${CMD_TEMPLATE//__SAMPLE_ID__/${SAMPLE_ID}}"
cmd="${cmd//__INPUT_DIR__/${SAMPLE_INPUT_DIR}}"
cmd="${cmd//__OUTPUT_DIR__/${NATIVE_OUT_DIR}}"
cmd="${cmd//__TARGET_YAML__/${TARGET_YAML}}"
cmd="${cmd//__OVERRIDES_FILE__/${OVERRIDES_FILE}}"

if [[ ! -d "${MODEL_WORKDIR}" ]]; then
  echo "[ERROR] germinal 工作目录不存在: ${MODEL_WORKDIR}"
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

echo "[OK] germinal native runner 完成: ${SAMPLE_ID}"
