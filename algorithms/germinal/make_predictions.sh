#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR=""
OUTPUT_DIR=""
MODEL_NAME="germinal"

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
CMD_TEMPLATE="${GERMINAL_CMD:-}"
MODEL_WORKDIR="${GERMINAL_WORKDIR:-/home/yqsong/projects/antibody_benchmark/models/germinal}"
DRY_RUN="${GERMINAL_DRY_RUN:-0}"
MAX_TRAJ="${GERMINAL_MAX_TRAJ:-1}"
MAX_PASSING="${GERMINAL_MAX_PASSING:-1}"
MAX_HAL="${GERMINAL_MAX_HALLUCINATED:-1}"
MAX_SAMPLES="${GERMINAL_MAX_SAMPLES:-0}"
USE_CONDA_RUN="${GERMINAL_USE_CONDA_RUN:-1}"
CONDA_ENV="${GERMINAL_CONDA_ENV:-germinal}"
HAS_CONDA=0
SKIP_TORCH_CHECK="${GERMINAL_SKIP_TORCH_CHECK:-0}"

if [[ "${USE_CONDA_RUN}" == "1" ]]; then
  if command -v conda >/dev/null 2>&1; then
    HAS_CONDA=1
    echo "[INFO] germinal 将使用 conda 环境: ${CONDA_ENV}"
  else
    echo "[WARN] 未找到 conda，germinal 将在当前环境执行。"
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

check_torch_version() {
  run_model_cmd "python -c \"import re,sys,torch; m=re.match(r'(\\d+\\.\\d+\\.\\d+)', torch.__version__); v=tuple(map(int, m.group(1).split('.'))) if m else (0,0,0); sys.exit(0 if v >= (2,6,0) else 3)\""
}

if [[ "${SKIP_TORCH_CHECK}" != "1" ]]; then
  if ! check_torch_version >/dev/null 2>&1; then
    echo "[ERROR] germinal 环境的 torch 版本低于 2.6，当前 transformers 将拒绝 torch.load。"
    echo "[ERROR] 请先升级（示例命令）:"
    echo "  /home/yqsong/.conda/envs/${CONDA_ENV}/bin/python -m pip install --upgrade 'torch>=2.6,<2.7' --index-url https://download.pytorch.org/whl/cu121 --extra-index-url https://pypi.org/simple"
    echo "[ERROR] 如需临时跳过检测，可设置 GERMINAL_SKIP_TORCH_CHECK=1"
    exit 1
  fi
fi

if [[ "${DRY_RUN}" != "1" && -z "${CMD_TEMPLATE}" ]]; then
  if run_model_cmd "python -c 'import pyrosetta'" >/dev/null 2>&1; then
    CMD_TEMPLATE="python run_germinal.py target.target_pdb_path='__TARGET__' target.target_chain='A' target.binder_chain='B' target.target_hotspots=__HOTSPOTS_ESC__ results_dir='__OUTDIR__' experiment_name='__SAMPLE__' max_trajectories=__MAX_TRAJ__ max_hallucinated_trajectories=__MAX_HAL__ max_passing_designs=__MAX_PASSING__"
    echo "[INFO] 未设置 GERMINAL_CMD，使用内置默认命令模板。"
  else
    echo "[WARN] 未检测到 pyrosetta，自动切换 germinal dry-run。"
    echo "       若需真实运行，请确认 conda 环境 ${CONDA_ENV} 已正确安装依赖。"
    DRY_RUN="1"
  fi
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

prepare_target_pdb() {
  local sample_json="$1"
  local input_pdb="$2"
  local output_pdb="$3"
  SAMPLE_JSON="${sample_json}" INPUT_PDB="${input_pdb}" OUTPUT_PDB="${output_pdb}" python3 - <<'PY'
import json
import os
from Bio.PDB import PDBIO, PDBParser, Select

sample_json = os.environ["SAMPLE_JSON"]
input_pdb = os.environ["INPUT_PDB"]
output_pdb = os.environ["OUTPUT_PDB"]

with open(sample_json, "r", encoding="utf-8") as f:
  d = json.load(f)
antigen_chain = str(d.get("antigen_chain", "")).strip()

parser = PDBParser(QUIET=True)
structure = parser.get_structure("target", input_pdb)

chains = []
for model in structure:
  chains = [c.id for c in model]
  break

if antigen_chain and antigen_chain not in chains:
  if len(chains) == 1:
    antigen_chain = chains[0]
  else:
    print("ERR_MULTICHAIN_MISMATCH", end="")
    raise SystemExit(3)
elif not antigen_chain:
  if len(chains) == 1:
    antigen_chain = chains[0]
  else:
    antigen_chain = chains[0]

class ChainSelector(Select):
  def accept_chain(self, chain):
    return 1 if chain.id == antigen_chain else 0

io = PDBIO()
io.set_structure(structure)
io.save(output_pdb, ChainSelector())

structure_a = parser.get_structure("target_a", output_pdb)
for model in structure_a:
  for chain in model:
    chain.id = "A"
io2 = PDBIO()
io2.set_structure(structure_a)
io2.save(output_pdb)
print("OK", end="")
PY
}

guess_hotspots() {
  local target_pdb="$1"
  TARGET_PDB="${target_pdb}" python3 - <<'PY'
import os
nums = []
with open(os.environ["TARGET_PDB"], "r", encoding="utf-8", errors="ignore") as f:
  for line in f:
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
      continue
    if len(line) < 27:
      continue
    if line[21] != "A":
      continue
    s = line[22:26].strip()
    if not s:
      continue
    try:
      nums.append(int(s))
    except ValueError:
      pass
nums = sorted(set(nums))
if not nums:
  print("A10,A20,A30", end="")
else:
  a, b, c = nums[len(nums)//4], nums[len(nums)//2], nums[3*len(nums)//4]
  print(f"A{a},A{b},A{c}", end="")
PY
}

shopt -s nullglob
processed_count=0
for sample_json in "${INPUT_DIR}"/*.json; do
  if [[ "${MAX_SAMPLES}" != "0" && "${processed_count}" -ge "${MAX_SAMPLES}" ]]; then
    echo "[INFO] 达到 GERMINAL_MAX_SAMPLES=${MAX_SAMPLES}，停止继续样本。"
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

  target_pdb="${sample_out}/target_chain_A.pdb"
  if ! prep_status="$(prepare_target_pdb "${sample_json}" "${input_pdb}" "${target_pdb}" 2>/dev/null)"; then
    echo "[WARN] ${sample_id}: 目标链准备失败，跳过"
    continue
  fi
  if [[ "${prep_status}" != "OK" ]]; then
    echo "[WARN] ${sample_id}: 目标链准备失败(${prep_status})，跳过"
    continue
  fi
  hotspots="$(guess_hotspots "${target_pdb}")"
  hotspots_esc="\\'${hotspots}\\'"

  if [[ "${DRY_RUN}" == "1" || -z "${CMD_TEMPLATE}" ]]; then
    cp "${target_pdb}" "${OUTPUT_DIR}/${sample_id}.pdb"
    echo "[OK] ${MODEL_NAME} (dry-run) 产物: ${OUTPUT_DIR}/${sample_id}.pdb"
    continue
  fi

  cmd="${CMD_TEMPLATE//__INPUT__/${input_pdb}}"
  cmd="${cmd//__OUTDIR__/${sample_out}}"
  cmd="${cmd//__SAMPLE__/${sample_id}}"
  cmd="${cmd//__TARGET__/${target_pdb}}"
  cmd="${cmd//__HOTSPOTS__/${hotspots}}"
  cmd="${cmd//__HOTSPOTS_ESC__/${hotspots_esc}}"
  cmd="${cmd//__MAX_TRAJ__/${MAX_TRAJ}}"
  cmd="${cmd//__MAX_HAL__/${MAX_HAL}}"
  cmd="${cmd//__MAX_PASSING__/${MAX_PASSING}}"

  if run_model_cmd "${cmd}" >"${sample_out}/stdout.log" 2>"${sample_out}/stderr.log"; then
    pred_file="$(pick_first_structure "${sample_out}")"
    if [[ -z "${pred_file}" ]]; then
      pred_file="$(pick_first_structure "${MODEL_WORKDIR}/results")"
    fi
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
