"""antibench/collection.py — 候选标准化收集。

从各模型的 outputs/runs/<model>/<target_id>/ 中读取出候选序列和结构，
转成统一的 Candidate 对象列表，写入 outputs/candidates/candidates_manifest.csv。

设计策略：
- 优先读取现有 collect_candidates.py 产出的 candidate_manifest.csv（复用）
- 若不存在，直接扫描 outputs/runs/<model>/<target_id>/ 中的 FASTA/PDB 文件
- 按模型内部得分（internal_score）排序，写入 rank_in_model
- 支持追加 native_reference 候选（source_type="native_reference"）
"""
from __future__ import annotations

import re
import csv
from pathlib import Path
from typing import Any

from antibench.candidate import (
    Candidate,
    CANDIDATE_CSV_FIELDS,
    load_candidates,
    save_candidates,
)
from antibench.dataset import Target
from antibench.utils import (
    PROJECT_ROOT,
    PathManager,
    get_logger,
    load_model_config,
    read_csv,
    read_json,
    safe_str,
    to_project_rel,
)

logger = get_logger("collection")

FASTA_EXTENSIONS = frozenset({".fa", ".faa", ".fasta", ".fas", ".seq"})
STRUCTURE_EXTENSIONS = frozenset({".pdb", ".cif", ".mmcif"})


# ─────────────────────────────────────────────────────────────────────────────
# 序列工具
# ─────────────────────────────────────────────────────────────────────────────

def _normalize_seq(seq: str) -> str:
    return "".join((seq or "").split()).upper()


def _parse_fasta(path: Path) -> list[tuple[str, str]]:
    """解析 FASTA 文件，返回 (header, sequence) 列表。"""
    records = []
    header = ""
    seq_lines: list[str] = []
    try:
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            line = line.strip()
            if line.startswith(">"):
                if header:
                    records.append((header, _normalize_seq("".join(seq_lines))))
                header = line[1:].strip()
                seq_lines = []
            elif line:
                seq_lines.append(line)
        if header:
            records.append((header, _normalize_seq("".join(seq_lines))))
    except Exception as exc:
        logger.warning(f"解析 FASTA 失败 {path}: {exc}")
    return records


# PDB → 单字母序列的氨基酸映射
_AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def _extract_seq_from_pdb(pdb_path: Path) -> str:
    """从 PDB 文件的 ATOM 记录提取单字母序列（取第一条链最长序列）。

    不依赖 biopython 以减少依赖，纯文本解析 ATOM 行。
    返回空字符串表示失败。
    """
    chains: dict[str, dict[int, str]] = {}  # chain_id → {res_seq: one_letter}
    try:
        for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
            if not line.startswith("ATOM"):
                continue
            if len(line) < 26:
                continue
            chain_id = line[21].strip()
            try:
                res_seq = int(line[22:26].strip())
            except ValueError:
                continue
            res_name = line[17:20].strip()
            aa = _AA3TO1.get(res_name)
            if aa is None:
                continue
            if chain_id not in chains:
                chains[chain_id] = {}
            chains[chain_id].setdefault(res_seq, aa)
    except Exception as exc:
        logger.warning(f"解析 PDB 失败 {pdb_path}: {exc}")
        return ""

    if not chains:
        return ""
    # 取残基数最多的链
    best_chain = max(chains, key=lambda c: len(chains[c]))
    seq = "".join(chains[best_chain][k] for k in sorted(chains[best_chain]))
    return seq


def _extract_score_from_metrics_csv(
    run_dir: Path,
    pdb_stem: str,
    meta_patterns: list[str],
    score_field: str,
) -> float | None:
    """从 germinal 风格的 `<stem>_metrics.csv` 中读取最后一行得分。

    Args:
        run_dir: 目标 run 目录（扫描 meta_patterns 的根目录）。
        pdb_stem: PDB 文件 stem（如 "8dtn_A_B_nb_s108166"），
            用于匹配 `<stem>_metrics.csv`。
        meta_patterns: 来自 model_cfg["output_patterns"]["meta"] 的 glob 列表。
        score_field: 要读取的列名（如 "ablm_ll"）。

    Returns:
        最后一行的 score_field 值，失败返回 None。
    """
    # 优先精确匹配该 stem 对应的 metrics CSV
    candidate_name = f"{pdb_stem}_metrics.csv"
    for pattern in meta_patterns:
        for path in sorted(run_dir.glob(pattern)):
            if path.name == candidate_name:
                return _read_last_csv_field(path, score_field)
    # 如果没有精确匹配，尝试模糊（只在 meta_patterns 结果里搜）
    return None


def _read_last_csv_field(csv_path: Path, field: str) -> float | None:
    """读取 CSV 最后一行的指定字段，返回 float 或 None。"""
    try:
        import csv as _csv
        rows = []
        with csv_path.open(encoding="utf-8", errors="replace") as f:
            reader = _csv.DictReader(f)
            for row in reader:
                rows.append(row)
        if not rows or field not in rows[-1]:
            return None
        return float(rows[-1][field])
    except Exception:
        return None


def _slugify(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return text.strip("._") or "candidate"


# ─────────────────────────────────────────────────────────────────────────────
# 模型输出解析（从 run 目录直接扫描）
# ─────────────────────────────────────────────────────────────────────────────

def _glob_first(root: Path, patterns: list[str]) -> Path | None:
    for pattern in patterns:
        matches = sorted(root.glob(pattern))
        if matches:
            return matches[0]
    return None


def _glob_all(root: Path, patterns: list[str]) -> list[Path]:
    seen: set[str] = set()
    out: list[Path] = []
    for pattern in patterns:
        for m in sorted(root.glob(pattern)):
            key = str(m.resolve())
            if key not in seen:
                seen.add(key)
                out.append(m)
    return out


def _extract_score_from_json(meta_path: Path, score_field: str) -> float | None:
    data = read_json(meta_path)
    if not data:
        return None
    # 递归查找第一层 key
    value = data.get(score_field)
    if value is not None:
        try:
            return float(value)
        except (TypeError, ValueError):
            pass
    return None


# ─────────────────────────────────────────────────────────────────────────────
# 从现有 collect_candidates.py 输出 CSV 中读取候选
# ─────────────────────────────────────────────────────────────────────────────

def _candidates_from_legacy_csv(
    legacy_csv: Path,
    target: Target,
    model_name: str,
    model_cfg: dict[str, Any],
) -> list[Candidate]:
    """读取现有 collect_candidates.py 输出的 candidate_manifest.csv，
    转成新 Candidate 对象。"""
    rows = read_csv(legacy_csv)
    if not rows:
        return []

    score_field = model_cfg.get("internal_score_field", "")
    ascending = model_cfg.get("internal_score_ascending", False)

    # 按 native_score_value 排序（复用旧字段）
    def _score(row: dict[str, str]) -> float:
        v = row.get("native_score_value") or row.get("score") or ""
        try:
            return float(v)
        except (ValueError, TypeError):
            return float("-inf") if not ascending else float("inf")

    rows_sorted = sorted(rows, key=_score, reverse=not ascending)

    candidates: list[Candidate] = []
    for rank, row in enumerate(rows_sorted, start=1):
        status = safe_str(row.get("status") or "ok")
        seq = _normalize_seq(row.get("pred_sequence") or row.get("sequence") or "")
        structure_str = safe_str(
            row.get("source_structure_path")
            or row.get("structure_path")
            or ""
        )
        structure_path = Path(structure_str) if structure_str else None
        if structure_path and not structure_path.is_absolute():
            structure_path = PROJECT_ROOT / structure_path

        score_value: float | None = None
        score_str = safe_str(row.get("native_score_value") or "")
        if score_str:
            try:
                score_value = float(score_str)
            except ValueError:
                pass

        score_name = safe_str(row.get("native_score_name") or score_field)
        duration_str = safe_str(row.get("source_duration_sec") or row.get("duration_sec") or "")
        duration: float | None = None
        try:
            duration = float(duration_str) if duration_str else None
        except ValueError:
            pass

        source_file_str = safe_str(
            row.get("source_meta_path")
            or row.get("source_csv_path")
            or row.get("source_file")
            or ""
        )
        source_file = Path(source_file_str) if source_file_str else legacy_csv

        candidate_id = (
            safe_str(row.get("candidate_id") or "")
            or f"{model_name}__{target.target_id}__{rank:04d}"
        )

        c = Candidate(
            candidate_id=candidate_id,
            target_id=target.target_id,
            model_name=model_name,
            source_type="generated",
            sequence=seq,
            structure_path=structure_path,
            internal_score_name=score_name,
            internal_score_value=score_value,
            rank_in_model=rank,
            source_file=source_file,
            duration_sec=duration,
            status=status if status in ("ok", "failed", "empty_seq") else "ok",
            error_summary=safe_str(row.get("source_error_summary") or ""),
            crop_rule_version=target.crop_rule_version,
        )
        candidates.append(c)

    return candidates


# ─────────────────────────────────────────────────────────────────────────────
# 从 run 目录直接扫描（无 legacy CSV 时的后备路径）
# ─────────────────────────────────────────────────────────────────────────────

def _candidates_from_run_dir(
    run_dir: Path,
    target: Target,
    model_name: str,
    model_cfg: dict[str, Any],
) -> list[Candidate]:
    """直接扫描 outputs/runs/<model>/<target_id>/ 提取候选。

    策略：
    1. 优先从 FASTA 文件提取序列（通用方式）
    2. 若 FASTA 模式为空或无匹配，则 fallback 到结构文件（PDB/CIF）提取序列
       （适用于 germinal 等只输出 PDB 而无 FASTA 的模型）
    3. 得分来源：meta JSON (result.json) 或 metrics CSV (<stem>_metrics.csv)
    """
    patterns = model_cfg.get("output_patterns", {})
    seq_patterns = patterns.get("sequence", ["**/*.fasta", "**/*.fa"])
    struct_patterns = patterns.get("structure", ["**/*.pdb", "**/*.cif"])
    meta_patterns = patterns.get("meta", ["**/result.json"])

    score_field = model_cfg.get("internal_score_field", "")
    ascending = model_cfg.get("internal_score_ascending", False)

    candidates_raw: list[tuple[float | None, Candidate]] = []
    seen_seqs: set[str] = set()

    # ── 路径 1: FASTA ──────────────────────────────────────────────────────
    fasta_files = _glob_all(run_dir, seq_patterns) if seq_patterns else []
    for fasta_path in fasta_files:
        records = _parse_fasta(fasta_path)
        for header, seq in records:
            if not seq or seq in seen_seqs:
                continue
            seen_seqs.add(seq)

            # 尝试找同名 PDB
            structure_path: Path | None = None
            stem = fasta_path.stem
            for ext in (".pdb", ".cif"):
                candidate_struct = fasta_path.parent / f"{stem}{ext}"
                if candidate_struct.exists():
                    structure_path = candidate_struct
                    break

            # 尝试从 meta JSON 读内部得分
            score: float | None = None
            meta_path = _glob_first(fasta_path.parent, meta_patterns)
            if meta_path and score_field:
                score = _extract_score_from_json(meta_path, score_field)

            placeholder_id = f"{model_name}__{target.target_id}__tmp"
            c = Candidate(
                candidate_id=placeholder_id,
                target_id=target.target_id,
                model_name=model_name,
                sequence=seq,
                structure_path=structure_path,
                internal_score_name=score_field,
                internal_score_value=score,
                rank_in_model=0,
                source_file=fasta_path,
                crop_rule_version=target.crop_rule_version,
            )
            candidates_raw.append((score, c))

    # ── 路径 2: 结构文件（PDB/CIF）fallback ────────────────────────────────
    # 当 FASTA 列表为空时（如 germinal 模型）才启用
    if not candidates_raw and struct_patterns:
        struct_files = _glob_all(run_dir, struct_patterns)
        if not struct_files:
            logger.debug(f"[collect/{model_name}] {run_dir} 中未找到 FASTA/PDB，跳过")
        for struct_path in struct_files:
            if struct_path.suffix.lower() not in STRUCTURE_EXTENSIONS:
                continue
            seq = _extract_seq_from_pdb(struct_path)
            if not seq or seq in seen_seqs:
                continue
            seen_seqs.add(seq)

            # 从同 stem 的 _metrics.csv 读得分
            pdb_stem = struct_path.stem
            score_val: float | None = None
            if score_field and meta_patterns:
                score_val = _extract_score_from_metrics_csv(
                    run_dir, pdb_stem, meta_patterns, score_field
                )

            placeholder_id = f"{model_name}__{target.target_id}__tmp"
            c = Candidate(
                candidate_id=placeholder_id,
                target_id=target.target_id,
                model_name=model_name,
                sequence=seq,
                structure_path=struct_path,
                internal_score_name=score_field,
                internal_score_value=score_val,
                rank_in_model=0,
                source_file=struct_path,
                crop_rule_version=target.crop_rule_version,
            )
            candidates_raw.append((score_val, c))

        if candidates_raw:
            logger.info(
                f"[collect/{model_name}] {run_dir.name}: "
                f"从 {len(candidates_raw)} 个 PDB 文件提取序列（无 FASTA 模式）"
            )

    # ── 排序 + 写 rank_in_model ─────────────────────────────────────────────
    def _sort_key(item: tuple[float | None, Candidate]) -> tuple[int, float]:
        s, _ = item
        if s is None:
            return (1, 0.0)
        return (0, s if ascending else -s)

    candidates_raw.sort(key=_sort_key)

    candidates: list[Candidate] = []
    for rank, (_, c) in enumerate(candidates_raw, start=1):
        c.rank_in_model = rank
        c.candidate_id = f"{model_name}__{target.target_id}__{rank:04d}"
        candidates.append(c)

    return candidates


# ─────────────────────────────────────────────────────────────────────────────
# 公开接口
# ─────────────────────────────────────────────────────────────────────────────

def collect_model_candidates(
    target: Target,
    model_name: str,
    path_manager: PathManager | None = None,
    model_cfg: dict[str, Any] | None = None,
) -> list[Candidate]:
    """收集单个模型针对单个靶点的所有候选。

    查找顺序：
    1. outputs/runs/<model>/<target_id>/candidate_manifest.csv（旧 collect_candidates.py 输出）
    2. outputs/runs/<model>/<target_id>/ 直接扫描
    """
    pm = path_manager or PathManager()
    cfg = model_cfg or load_model_config(model_name)
    run_dir = pm.run_dir(model_name, target.target_id)

    if not run_dir.exists():
        logger.debug(f"[{model_name}] run_dir 不存在: {run_dir}")
        return []

    # 路径 1：旧 collect_candidates.py 输出
    legacy_csv = run_dir / "candidate_manifest.csv"
    if legacy_csv.exists():
        candidates = _candidates_from_legacy_csv(legacy_csv, target, model_name, cfg)
        if candidates:
            logger.debug(
                f"[{model_name}/{target.target_id}] 从 legacy CSV 读取 {len(candidates)} 个候选"
            )
            return candidates

    # 路径 2：直接扫描
    candidates = _candidates_from_run_dir(run_dir, target, model_name, cfg)
    logger.debug(
        f"[{model_name}/{target.target_id}] 扫描 run_dir 得到 {len(candidates)} 个候选"
    )
    return candidates


def collect_all(
    targets: list[Target],
    model_names: list[str],
    path_manager: PathManager | None = None,
    out_path: Path | None = None,
    include_failed: bool = False,
) -> list[Candidate]:
    """收集所有模型、所有靶点的候选，写入 candidates_manifest.csv。

    Args:
        targets: Target 列表（来自 dataset.load_targets()）。
        model_names: 要收集的模型名列表。
        path_manager: PathManager 实例。
        out_path: 输出 CSV 路径，None 则用 pm.candidates_manifest。
        include_failed: 是否包含状态为 failed 的候选。

    Returns:
        全部 Candidate 对象列表。
    """
    pm = path_manager or PathManager()
    out_path = out_path or pm.candidates_manifest

    all_candidates: list[Candidate] = []

    for model_name in model_names:
        try:
            model_cfg = load_model_config(model_name)
        except FileNotFoundError:
            logger.warning(f"模型配置不存在，跳过: {model_name}")
            continue

        for target in targets:
            candidates = collect_model_candidates(
                target, model_name, pm, model_cfg
            )
            if not include_failed:
                candidates = [c for c in candidates if c.status != "failed"]
            all_candidates.extend(candidates)

    logger.info(f"[collection] 共收集 {len(all_candidates)} 个候选（{len(model_names)} 个模型）")
    save_candidates(all_candidates, out_path)
    return all_candidates
