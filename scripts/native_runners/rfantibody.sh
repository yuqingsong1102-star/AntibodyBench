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

CONFIG_FILE="${SAMPLE_INPUT_DIR}/rfantibody_config.env"
if [[ ! -f "${CONFIG_FILE}" ]]; then
  echo "[ERROR] 缺少输入配置: ${CONFIG_FILE}"
  exit 3
fi

# shellcheck disable=SC1090
source "${CONFIG_FILE}"

MODEL_WORKDIR="${RFANTIBODY_WORKDIR:-/home/yqsong/projects/antibody_benchmark/models/RFantibody}"
NATIVE_OUT_DIR="$(cd "$(dirname "${SAMPLE_OUTPUT_DIR}")" && pwd)/$(basename "${SAMPLE_OUTPUT_DIR}")/native_run"
mkdir -p "${NATIVE_OUT_DIR}"

# Auto-select GPU with most free memory if CUDA_VISIBLE_DEVICES not set
if [[ -z "${CUDA_VISIBLE_DEVICES:-}" ]] && command -v nvidia-smi >/dev/null 2>&1; then
  BEST_GPU=$(nvidia-smi --query-gpu=index,memory.free --format=csv,noheader,nounits | sort -t',' -k2 -nr | head -1 | cut -d',' -f1 | tr -d ' ')
  export CUDA_VISIBLE_DEVICES="${BEST_GPU}"
  echo "[INFO] Auto-selected GPU ${BEST_GPU} (most free memory)"
fi

if [[ -z "${RFANTIBODY_CMD:-}" ]]; then
  echo "[ERROR] 未设置 RFANTIBODY_CMD。"
  echo "示例: export RFANTIBODY_CMD='python run_pipeline.py --input __TARGET_PDB__ --out __OUTPUT_DIR__ --sample __SAMPLE_ID__'"
  echo "可用占位符: __SAMPLE_ID__ __INPUT_DIR__ __OUTPUT_DIR__ __CONFIG_FILE__ __TARGET_PDB__ __FRAMEWORK_PDB__ __HOTSPOTS__ __DESIGN_LOOPS__ __NUM_DESIGNS__ __NUM_SEQS__ __NUM_RECYCLES__"
  exit 11
fi

cmd="${RFANTIBODY_CMD}"
cmd="${cmd//__SAMPLE_ID__/${SAMPLE_ID}}"
cmd="${cmd//__INPUT_DIR__/${SAMPLE_INPUT_DIR}}"
cmd="${cmd//__OUTPUT_DIR__/${NATIVE_OUT_DIR}}"
cmd="${cmd//__CONFIG_FILE__/${CONFIG_FILE}}"
cmd="${cmd//__TARGET_PDB__/${TARGET_PDB:-}}"
cmd="${cmd//__FRAMEWORK_PDB__/${FRAMEWORK_PDB:-}}"
cmd="${cmd//__HOTSPOTS__/${HOTSPOTS:-}}"
cmd="${cmd//__DESIGN_LOOPS__/${DESIGN_LOOPS:-}}"
cmd="${cmd//__NUM_DESIGNS__/${NUM_DESIGNS:-10}}"
cmd="${cmd//__NUM_SEQS__/${NUM_SEQS:-4}}"
cmd="${cmd//__NUM_RECYCLES__/${NUM_RECYCLES:-10}}"

if [[ ! -d "${MODEL_WORKDIR}" ]]; then
  echo "[ERROR] RFantibody 工作目录不存在: ${MODEL_WORKDIR}"
  exit 4
fi

(
  cd "${MODEL_WORKDIR}"
  bash -lc "${cmd}"
)

echo "[OK] RFantibody native runner 完成: ${SAMPLE_ID}"
