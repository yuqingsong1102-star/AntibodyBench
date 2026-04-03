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

SPEC_FILE="${SAMPLE_INPUT_DIR}/design_spec.yaml"
if [[ ! -f "${SPEC_FILE}" ]]; then
  echo "[ERROR] 缺少输入文件: ${SPEC_FILE}"
  exit 3
fi

MODEL_WORKDIR="${BOLTZGEN_WORKDIR:-/home/yqsong/projects/antibody_benchmark/models/boltzgen}"
CONDA_ENV="${BOLTZGEN_CONDA_ENV:-bg}"
USE_CONDA_RUN="${BOLTZGEN_USE_CONDA_RUN:-1}"
# Default to single GPU to avoid distributed communication failures
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
PROTOCOL="${BOLTZGEN_PROTOCOL:-protein-anything}"
NUM_DESIGNS="${BOLTZGEN_NUM_DESIGNS:-2}"
BUDGET="${BOLTZGEN_BUDGET:-1}"
NATIVE_OUT_DIR="$(cd "$(dirname "${SAMPLE_OUTPUT_DIR}")" && pwd)/$(basename "${SAMPLE_OUTPUT_DIR}")/native_run"
mkdir -p "${NATIVE_OUT_DIR}"

default_cmd="boltzgen run \"${SPEC_FILE}\" --output \"${NATIVE_OUT_DIR}\" --protocol \"${PROTOCOL}\" --num_designs ${NUM_DESIGNS} --budget ${BUDGET}"
CMD_TEMPLATE="${BOLTZGEN_CMD:-${default_cmd}}"
cmd="${CMD_TEMPLATE//__SAMPLE_ID__/${SAMPLE_ID}}"
cmd="${cmd//__INPUT_DIR__/${SAMPLE_INPUT_DIR}}"
cmd="${cmd//__OUTPUT_DIR__/${NATIVE_OUT_DIR}}"
cmd="${cmd//__SPEC_FILE__/${SPEC_FILE}}"
cmd="${cmd//__PROTOCOL__/${PROTOCOL}}"
cmd="${cmd//__NUM_DESIGNS__/${NUM_DESIGNS}}"
cmd="${cmd//__BUDGET__/${BUDGET}}"

if [[ ! -d "${MODEL_WORKDIR}" ]]; then
  echo "[ERROR] boltzgen 工作目录不存在: ${MODEL_WORKDIR}"
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

echo "[OK] boltzgen native runner 完成: ${SAMPLE_ID}"
