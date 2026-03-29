#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR=""
OUTPUT_DIR=""
MODEL_NAME="BindCraft"

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
CMD_TEMPLATE="${BINDCRAFT_CMD:-}"
MODEL_WORKDIR="${BINDCRAFT_WORKDIR:-/home/yqsong/projects/antibody_benchmark/models/BindCraft}"
DRY_RUN="${BINDCRAFT_DRY_RUN:-0}"
MAX_SAMPLES="${BINDCRAFT_MAX_SAMPLES:-0}"
LENGTH_MIN="${BINDCRAFT_LENGTH_MIN:-60}"
LENGTH_MAX="${BINDCRAFT_LENGTH_MAX:-90}"
NUM_FINAL="${BINDCRAFT_NUM_FINAL:-1}"
MAX_TRAJ="${BINDCRAFT_MAX_TRAJ:-1}"
FILTERS_JSON="${BINDCRAFT_FILTERS:-${MODEL_WORKDIR}/settings_filters/default_filters.json}"
ADVANCED_TEMPLATE="${BINDCRAFT_ADVANCED_TEMPLATE:-${MODEL_WORKDIR}/settings_advanced/default_4stage_multimer.json}"
USE_CONDA_RUN="${BINDCRAFT_USE_CONDA_RUN:-1}"
CONDA_ENV="${BINDCRAFT_CONDA_ENV:-BindCraft}"
HAS_CONDA=0

if [[ "${USE_CONDA_RUN}" == "1" ]]; then
  if command -v conda >/dev/null 2>&1; then
    HAS_CONDA=1
    echo "[INFO] BindCraft 将使用 conda 环境: ${CONDA_ENV}"
  else
    echo "[WARN] 未找到 conda，BindCraft 将在当前环境执行。"
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
  CMD_TEMPLATE="python -u ./bindcraft.py --settings '__SETTINGS__' --filters '__FILTERS__' --advanced '__ADVANCED__'"
  echo "[INFO] 未设置 BINDCRAFT_CMD，使用内置默认命令模板。"
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
    echo "[INFO] 达到 BINDCRAFT_MAX_SAMPLES=${MAX_SAMPLES}，停止继续样本。"
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
    SAMPLE_JSON="${sample_json}" INPUT_PDB="${input_pdb}" python3 - <<'PY'
import json
import os
from Bio.PDB import PDBParser

with open(os.environ["SAMPLE_JSON"], "r", encoding="utf-8") as f:
  d = json.load(f)
chain = str(d.get("antigen_chain", "")).strip()
if "," in chain:
  chain = chain.split(",", 1)[0].strip()
if chain:
  print(chain, end="")
  raise SystemExit(0)

parser = PDBParser(QUIET=True)
structure = parser.get_structure("target", os.environ["INPUT_PDB"])
for model in structure:
  chains = [c.id for c in model]
  if chains:
    print(chains[0], end="")
    break
  break
PY
  )"
  if [[ -z "${target_chain}" ]]; then
    target_chain="A"
  fi

  settings_json="${sample_out}/bindcraft_settings.json"
  advanced_json="${sample_out}/bindcraft_advanced.json"
  INPUT_PDB="${input_pdb}" SETTINGS_JSON="${settings_json}" ADVANCED_JSON="${advanced_json}" ADVANCED_TEMPLATE="${ADVANCED_TEMPLATE}" SAMPLE_ID="${sample_id}" SAMPLE_OUT="${sample_out}" TARGET_CHAIN="${target_chain}" LENGTH_MIN="${LENGTH_MIN}" LENGTH_MAX="${LENGTH_MAX}" NUM_FINAL="${NUM_FINAL}" MAX_TRAJ="${MAX_TRAJ}" python3 - <<'PY'
import json
import os
from pathlib import Path

settings_path = Path(os.environ["SETTINGS_JSON"])
advanced_path = Path(os.environ["ADVANCED_JSON"])
advanced_template = Path(os.environ["ADVANCED_TEMPLATE"])

target_settings = {
  "design_path": str(Path(os.environ["SAMPLE_OUT"]).resolve()),
  "binder_name": os.environ["SAMPLE_ID"],
  "starting_pdb": str(Path(os.environ["INPUT_PDB"]).resolve()),
  "chains": os.environ["TARGET_CHAIN"] or "A",
  "target_hotspot_residues": None,
  "lengths": [int(os.environ["LENGTH_MIN"]), int(os.environ["LENGTH_MAX"])],
  "number_of_final_designs": int(os.environ["NUM_FINAL"]),
}
settings_path.write_text(json.dumps(target_settings, ensure_ascii=False, indent=2), encoding="utf-8")

advanced = {}
if advanced_template.exists():
  with advanced_template.open("r", encoding="utf-8") as f:
    advanced = json.load(f)
advanced["max_trajectories"] = int(os.environ["MAX_TRAJ"])
advanced_path.write_text(json.dumps(advanced, ensure_ascii=False, indent=2), encoding="utf-8")
PY

  if [[ "${DRY_RUN}" == "1" || -z "${CMD_TEMPLATE}" ]]; then
    cp "${input_pdb}" "${OUTPUT_DIR}/${sample_id}.pdb"
    echo "[OK] ${MODEL_NAME} (dry-run) 产物: ${OUTPUT_DIR}/${sample_id}.pdb"
    continue
  fi

  cmd="${CMD_TEMPLATE//__INPUT__/${input_pdb}}"
  cmd="${cmd//__OUTDIR__/${sample_out}}"
  cmd="${cmd//__SAMPLE__/${sample_id}}"
  cmd="${cmd//__SETTINGS__/${settings_json}}"
  cmd="${cmd//__FILTERS__/${FILTERS_JSON}}"
  cmd="${cmd//__ADVANCED__/${advanced_json}}"

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
