#!/usr/bin/env bash
set -uo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

MODEL_NAME=""
INPUT_ROOT="${ROOT_DIR}/native_inputs"
OUT_ROOT="${ROOT_DIR}/outputs/native_predictions"
MAX_SAMPLES=1
SAMPLE_ID=""
CONTINUE_ON_ERROR=1

usage() {
  cat <<'EOF'
用法:
  bash scripts/run.sh --model <RFantibody|germinal|BindCraft|boltzgen> [选项]

选项:
  --model <name>            必填，模型名（大小写按目录名）
  --input-root <path>       原生输入根目录（默认: native_inputs）
  --out-root <path>         输出根目录（默认: outputs/native_predictions）
  --max-samples <N>         最多处理样本数（默认: 1，0 表示不限制）
  --sample-id <id>          仅运行指定样本（会忽略 --max-samples）
  --continue-on-error <0|1> 单样本失败是否继续（默认: 1）
  -h, --help                显示帮助
EOF
}

normalize_model() {
  local raw="$1"
  case "${raw}" in
    RFantibody|rfantibody|RFANTIBODY) echo "RFantibody" ;;
    germinal|GERMINAL) echo "germinal" ;;
    BindCraft|bindcraft|BINDCRAFT) echo "BindCraft" ;;
    boltzgen|BOLTZGEN|BoltzGen) echo "boltzgen" ;;
    *) return 1 ;;
  esac
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --model)
      MODEL_NAME="$(normalize_model "$2")" || {
        echo "[ERROR] 不支持的模型: $2"
        exit 2
      }
      shift 2
      ;;
    --input-root)
      INPUT_ROOT="$2"
      shift 2
      ;;
    --out-root)
      OUT_ROOT="$2"
      shift 2
      ;;
    --max-samples)
      MAX_SAMPLES="$2"
      shift 2
      ;;
    --sample-id)
      SAMPLE_ID="$2"
      shift 2
      ;;
    --continue-on-error)
      CONTINUE_ON_ERROR="$2"
      shift 2
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

if [[ -z "${MODEL_NAME}" ]]; then
  echo "[ERROR] 缺少必填参数 --model"
  usage
  exit 2
fi
if [[ ! "${MAX_SAMPLES}" =~ ^[0-9]+$ ]]; then
  echo "[ERROR] --max-samples 必须为非负整数，当前: ${MAX_SAMPLES}"
  exit 2
fi
if [[ "${CONTINUE_ON_ERROR}" != "0" && "${CONTINUE_ON_ERROR}" != "1" ]]; then
  echo "[ERROR] --continue-on-error 仅支持 0 或 1"
  exit 2
fi

MODEL_INPUT_DIR="${INPUT_ROOT}/${MODEL_NAME}"
MODEL_OUT_DIR="${OUT_ROOT}/${MODEL_NAME}"
RUNNER_DIR="${ROOT_DIR}/scripts/native_runners"
COLLECTOR="${RUNNER_DIR}/collect_candidates.py"

case "${MODEL_NAME}" in
  RFantibody) RUNNER="${RUNNER_DIR}/rfantibody.sh" ;;
  germinal) RUNNER="${RUNNER_DIR}/germinal.sh" ;;
  BindCraft) RUNNER="${RUNNER_DIR}/bindcraft.sh" ;;
  boltzgen) RUNNER="${RUNNER_DIR}/boltzgen.sh" ;;
esac

if [[ ! -d "${MODEL_INPUT_DIR}" ]]; then
  echo "[ERROR] 原生输入目录不存在: ${MODEL_INPUT_DIR}"
  exit 1
fi
if [[ ! -x "${RUNNER}" ]]; then
  echo "[ERROR] runner 不存在或不可执行: ${RUNNER}"
  exit 1
fi
if [[ ! -f "${COLLECTOR}" ]]; then
  echo "[ERROR] collector 不存在: ${COLLECTOR}"
  exit 1
fi

mkdir -p "${MODEL_OUT_DIR}"
MANIFEST_PATH="${MODEL_OUT_DIR}/manifest.csv"
python3 "${COLLECTOR}" --print-fieldnames > "${MANIFEST_PATH}"

collect_field() {
  local meta_path="$1"
  local field="$2"
  META_PATH="${meta_path}" FIELD_NAME="${field}" python3 - <<'PY'
import json
import os
from pathlib import Path

meta_path = Path(os.environ["META_PATH"])
field = os.environ["FIELD_NAME"]
if not meta_path.exists():
  print("", end="")
  raise SystemExit(0)
try:
  data = json.loads(meta_path.read_text(encoding="utf-8"))
except Exception:
  print("", end="")
  raise SystemExit(0)
value = data.get(field, "")
if value is None:
  value = ""
print(str(value), end="")
PY
}

declare -a SAMPLE_IDS
if [[ -n "${SAMPLE_ID}" ]]; then
  if [[ ! -d "${MODEL_INPUT_DIR}/${SAMPLE_ID}" ]]; then
    echo "[ERROR] 指定样本不存在: ${MODEL_INPUT_DIR}/${SAMPLE_ID}"
    exit 1
  fi
  SAMPLE_IDS=("${SAMPLE_ID}")
else
  while IFS= read -r sid; do
    SAMPLE_IDS+=("${sid}")
  done < <(python3 - <<'PY' "${MODEL_INPUT_DIR}"
import sys
from pathlib import Path

root = Path(sys.argv[1])
for p in sorted(root.iterdir()):
  if p.is_dir():
    print(p.name)
PY
)
fi

if [[ "${#SAMPLE_IDS[@]}" -eq 0 ]]; then
  echo "[WARN] 没有可运行样本: ${MODEL_INPUT_DIR}"
  exit 0
fi

if [[ "${MAX_SAMPLES}" != "0" && -z "${SAMPLE_ID}" && "${#SAMPLE_IDS[@]}" -gt "${MAX_SAMPLES}" ]]; then
  SAMPLE_IDS=("${SAMPLE_IDS[@]:0:${MAX_SAMPLES}}")
fi

echo "[INFO] 原生推理模式启动"
echo "[INFO] model=${MODEL_NAME}"
echo "[INFO] samples=${#SAMPLE_IDS[@]}"
echo "[INFO] input_root=${INPUT_ROOT}"
echo "[INFO] out_root=${OUT_ROOT}"

overall_exit=0
for sid in "${SAMPLE_IDS[@]}"; do
  sample_input_dir="${MODEL_INPUT_DIR}/${sid}"
  sample_out_dir="${MODEL_OUT_DIR}/${sid}"
  rm -rf "${sample_out_dir}"
  mkdir -p "${sample_out_dir}"
  stdout_log="${sample_out_dir}/run_stdout.log"
  stderr_log="${sample_out_dir}/run_stderr.log"
  run_meta_path="${sample_out_dir}/run_meta.json"
  candidate_manifest_path="${sample_out_dir}/candidate_manifest.csv"

  echo "[INFO] 运行样本: ${sid}"
  start_epoch="$(date +%s)"
  if "${RUNNER}" \
    --sample-id "${sid}" \
    --sample-input-dir "${sample_input_dir}" \
    --sample-output-dir "${sample_out_dir}" \
    --project-root "${ROOT_DIR}" \
    > "${stdout_log}" 2> "${stderr_log}"; then
    runner_exit=0
  else
    runner_exit=$?
  fi

  if ! python3 "${COLLECTOR}" \
    --model "${MODEL_NAME}" \
    --sample-id "${sid}" \
    --sample-input-dir "${sample_input_dir}" \
    --sample-output-dir "${sample_out_dir}" \
    --project-root "${ROOT_DIR}" \
    --runner-exit-code "${runner_exit}" \
    --duration-sec "$(( $(date +%s) - start_epoch ))"; then
    echo "[WARN] collect_candidates 异常: ${sid}"
  fi

  end_epoch="$(date +%s)"
  duration="$((end_epoch - start_epoch))"
  status="$(collect_field "${run_meta_path}" "status")"
  error_summary="$(collect_field "${run_meta_path}" "error_summary")"
  n_candidates="$(collect_field "${run_meta_path}" "n_candidates")"
  if [[ -z "${status}" ]]; then
    status="failed"
  fi

  if [[ -f "${candidate_manifest_path}" ]]; then
    tail -n +2 "${candidate_manifest_path}" >> "${MANIFEST_PATH}"
  fi

  if [[ -z "${n_candidates}" ]]; then
    n_candidates="0"
  fi

  if [[ "${status}" == "failed" || "${status}" == "partial" ]]; then
    overall_exit=1
    echo "[WARN] 样本未完全成功: ${sid} (candidates=${n_candidates}; ${error_summary:-unknown_error})"
    if [[ "${CONTINUE_ON_ERROR}" == "0" ]]; then
      echo "[ERROR] 因 --continue-on-error=0，提前终止。"
      break
    fi
  else
    echo "[OK] 样本完成: ${sid} (candidates=${n_candidates})"
  fi
done

echo "[INFO] manifest: ${MANIFEST_PATH}"
exit "${overall_exit}"
