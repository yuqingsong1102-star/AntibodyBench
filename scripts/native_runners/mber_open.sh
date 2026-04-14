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

SETTINGS_FILE="${SAMPLE_INPUT_DIR}/mber_vhh_settings.json"
if [[ ! -f "${SETTINGS_FILE}" ]]; then
  echo "[ERROR] 缺少输入配置: ${SETTINGS_FILE}"
  exit 3
fi

MODEL_WORKDIR="${MBER_OPEN_WORKDIR:-/home/yqsong/projects/antibody_benchmark/models/mber-open}"
CONDA_ENV="${MBER_OPEN_CONDA_ENV:-mber}"
USE_CONDA_RUN="${MBER_OPEN_USE_CONDA_RUN:-1}"
NATIVE_OUT_DIR="$(cd "$(dirname "${SAMPLE_OUTPUT_DIR}")" && pwd)/$(basename "${SAMPLE_OUTPUT_DIR}")/native_run"
RUNTIME_INPUT_DIR="${NATIVE_OUT_DIR}/runtime_inputs"
RUNTIME_SETTINGS_FILE="${RUNTIME_INPUT_DIR}/mber_vhh_settings.runtime.json"
mkdir -p "${RUNTIME_INPUT_DIR}"

python3 - <<'PY' "${SETTINGS_FILE}" "${RUNTIME_SETTINGS_FILE}" "${NATIVE_OUT_DIR}" "${SAMPLE_ID}" \
  "${MBER_OPEN_NUM_ACCEPTED:-}" "${MBER_OPEN_MAX_TRAJECTORIES:-}" "${MBER_OPEN_MIN_IPTM:-}" "${MBER_OPEN_MIN_PLDDT:-}" \
  "${MBER_OPEN_MASKED_SEQUENCE:-}" "${MBER_OPEN_SKIP_ANIMATIONS:-}" "${MBER_OPEN_SKIP_PICKLE:-}" "${MBER_OPEN_SKIP_PNG:-}"
import json
import sys
from pathlib import Path

src = Path(sys.argv[1])
dst = Path(sys.argv[2])
output_dir = sys.argv[3]
sample_id = sys.argv[4]
num_accepted_raw = sys.argv[5].strip()
max_trajectories_raw = sys.argv[6].strip()
min_iptm_raw = sys.argv[7].strip()
min_plddt_raw = sys.argv[8].strip()
masked_sequence_raw = sys.argv[9]
skip_animations_raw = sys.argv[10].strip()
skip_pickle_raw = sys.argv[11].strip()
skip_png_raw = sys.argv[12].strip()

payload = json.loads(src.read_text(encoding="utf-8"))
payload.setdefault("output", {})["dir"] = output_dir
payload.setdefault("target", {})["name"] = payload.get("target", {}).get("name") or sample_id

if num_accepted_raw:
  payload.setdefault("stopping", {})["num_accepted"] = int(num_accepted_raw)
if max_trajectories_raw:
  payload.setdefault("stopping", {})["max_trajectories"] = int(max_trajectories_raw)
if min_iptm_raw:
  payload.setdefault("filters", {})["min_iptm"] = float(min_iptm_raw)
if min_plddt_raw:
  payload.setdefault("filters", {})["min_plddt"] = float(min_plddt_raw)
if masked_sequence_raw:
  payload.setdefault("binder", {})["masked_sequence"] = masked_sequence_raw

def parse_bool(raw: str):
  if not raw:
    return None
  text = raw.lower()
  if text in {"1", "true", "yes", "on"}:
    return True
  if text in {"0", "false", "no", "off"}:
    return False
  raise SystemExit(f"[ERROR] 无法解析布尔值: {raw}")

skip_animations = parse_bool(skip_animations_raw)
skip_pickle = parse_bool(skip_pickle_raw)
skip_png = parse_bool(skip_png_raw)
if skip_animations is not None:
  payload.setdefault("output", {})["skip_animations"] = skip_animations
if skip_pickle is not None:
  payload.setdefault("output", {})["skip_pickle"] = skip_pickle
if skip_png is not None:
  payload.setdefault("output", {})["skip_png"] = skip_png

dst.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
PY

default_cmd='mber-vhh --settings "__SETTINGS_FILE__"'
CMD_TEMPLATE="${MBER_OPEN_CMD:-${default_cmd}}"
cmd="${CMD_TEMPLATE//__SAMPLE_ID__/${SAMPLE_ID}}"
cmd="${cmd//__INPUT_DIR__/${SAMPLE_INPUT_DIR}}"
cmd="${cmd//__OUTPUT_DIR__/${NATIVE_OUT_DIR}}"
cmd="${cmd//__SETTINGS_FILE__/${RUNTIME_SETTINGS_FILE}}"
cmd="${cmd//__CONFIG_FILE__/${RUNTIME_SETTINGS_FILE}}"

if [[ ! -d "${MODEL_WORKDIR}" ]]; then
  echo "[ERROR] mBER-open 工作目录不存在: ${MODEL_WORKDIR}"
  exit 4
fi

HF_HOME_VALUE="${MBER_OPEN_HF_HOME:-${HOME}/.mber/huggingface}"
HF_ENDPOINT_VALUE="${MBER_OPEN_HF_ENDPOINT:-${HF_ENDPOINT:-https://hf-mirror.com}}"
HF_HUB_ENABLE_HF_TRANSFER_VALUE="${MBER_OPEN_HF_HUB_ENABLE_HF_TRANSFER:-1}"
HF_HUB_OFFLINE_VALUE="${MBER_OPEN_HF_HUB_OFFLINE:-1}"

if [[ "${USE_CONDA_RUN}" == "1" ]] && command -v conda >/dev/null 2>&1; then
  (
    cd "${MODEL_WORKDIR}"
    conda run -n "${CONDA_ENV}" bash -lc "export HF_HOME='${HF_HOME_VALUE}' HF_ENDPOINT='${HF_ENDPOINT_VALUE}' HF_HUB_ENABLE_HF_TRANSFER='${HF_HUB_ENABLE_HF_TRANSFER_VALUE}' HF_HUB_OFFLINE='${HF_HUB_OFFLINE_VALUE}' TRANSFORMERS_OFFLINE='${HF_HUB_OFFLINE_VALUE}' && ${cmd}"
  )
else
  (
    cd "${MODEL_WORKDIR}"
    bash -lc "export HF_HOME='${HF_HOME_VALUE}' HF_ENDPOINT='${HF_ENDPOINT_VALUE}' HF_HUB_ENABLE_HF_TRANSFER='${HF_HUB_ENABLE_HF_TRANSFER_VALUE}' HF_HUB_OFFLINE='${HF_HUB_OFFLINE_VALUE}' TRANSFORMERS_OFFLINE='${HF_HUB_OFFLINE_VALUE}' && ${cmd}"
  )
fi

echo "[OK] mBER-open native runner 完成: ${SAMPLE_ID}"