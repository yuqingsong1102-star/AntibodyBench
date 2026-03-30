#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MODEL_NAME="template_model"
DATA_SPLIT_CSV="${ROOT_DIR}/inputs/antibody_datasets/dataset_index.csv"
PREPARE_ONLY=0

usage() {
  cat <<'EOF'
用法:
  bash scripts/run.sh --model <model_name> [--split-csv <path>] [--prepare-only]
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --model)
      MODEL_NAME="$2"
      shift 2
      ;;
    --split-csv)
      DATA_SPLIT_CSV="$2"
      shift 2
      ;;
    --prepare-only)
      PREPARE_ONLY=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "[ERROR] 未知参数: $1"
      usage
      exit 2
      ;;
  esac
done

MODEL_DIR="${ROOT_DIR}/algorithms/${MODEL_NAME}"
INPUT_DIR="${ROOT_DIR}/outputs/input/${MODEL_NAME}"
PRED_DIR="${ROOT_DIR}/outputs/prediction/${MODEL_NAME}"
EVAL_DIR="${ROOT_DIR}/outputs/evaluation/${MODEL_NAME}"

mkdir -p "${INPUT_DIR}" "${PRED_DIR}" "${EVAL_DIR}"

if [[ ! -d "${MODEL_DIR}" ]]; then
  echo "[ERROR] 模型目录不存在: ${MODEL_DIR}"
  exit 1
fi
if [[ ! -f "${DATA_SPLIT_CSV}" ]]; then
  echo "[ERROR] 数据索引不存在: ${DATA_SPLIT_CSV}"
  echo "请先准备 inputs/antibody_datasets/dataset_index.csv"
  exit 1
fi

echo "[INFO] 阶段1/3: 预处理"
python "${MODEL_DIR}/preprocess.py" \
  --af3-input "${ROOT_DIR}/inputs/alphafold3_inputs.json" \
  --dataset-csv "${DATA_SPLIT_CSV}" \
  --out-dir "${INPUT_DIR}"

if [[ "${PREPARE_ONLY}" == "1" ]]; then
  echo "[INFO] 已启用 --prepare-only，仅完成输入准备。"
  echo "[INFO] 产物目录: ${INPUT_DIR}"
  exit 0
fi

echo "[INFO] 阶段2/3: 推理"
bash "${MODEL_DIR}/make_predictions.sh" \
  --input-dir "${INPUT_DIR}" \
  --output-dir "${PRED_DIR}" \
  --model-name "${MODEL_NAME}"

echo "[INFO] 阶段3/3: 后处理"
python "${MODEL_DIR}/postprocess.py" \
  --dataset-csv "${DATA_SPLIT_CSV}" \
  --prediction-dir "${PRED_DIR}" \
  --out-dir "${EVAL_DIR}"

echo "[INFO] 指标汇总"
python "${ROOT_DIR}/scripts/aggregate_metrics.py" --model "${MODEL_NAME}"

echo "[INFO] 运行完成: ${MODEL_NAME}"
