"""antibench/candidate.py — Candidate 与 RunMetadata 数据对象。

Candidate：一个模型针对一个靶点生成的单个候选（序列 + 可选结构）。
RunMetadata：记录一次模型运行的时间、GPU、状态等元信息。

设计原则：
- source_type / is_reference 明确区分生成候选与天然参考
- internal_score 只用于该模型内部 top-K 截取，不跨模型比较
- rank_in_model 按 internal_score 排序后确定，由 collection.py 写入
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field, fields
from datetime import datetime
from pathlib import Path
from typing import Any

from antibench.utils import safe_str, to_project_rel, PROJECT_ROOT


# ─────────────────────────────────────────────────────────────────────────────
# RunMetadata — 运行元信息（pre-eval，来自 runner.py）
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class RunMetadata:
    """记录一次模型运行（单靶点）的完整元信息。"""
    model_name: str
    target_id: str
    start_time: datetime | None = None
    end_time: datetime | None = None
    wall_clock_sec: float | None = None   # end - start（秒）
    gpu_count: int = 1
    exit_code: int | None = None
    log_path: Path | None = None
    run_dir: Path | None = None           # outputs/runs/<model>/<target_id>/
    status: str = "unknown"               # "ok" | "failed" | "timeout" | "skipped"
    error_summary: str = ""
    extra: dict[str, Any] = field(default_factory=dict)

    @property
    def gpu_hours(self) -> float | None:
        if self.wall_clock_sec is None:
            return None
        return self.wall_clock_sec * self.gpu_count / 3600.0

    def to_dict(self) -> dict[str, str]:
        return {
            "model_name": self.model_name,
            "target_id": self.target_id,
            "start_time": self.start_time.isoformat() if self.start_time else "",
            "end_time": self.end_time.isoformat() if self.end_time else "",
            "wall_clock_sec": f"{self.wall_clock_sec:.2f}" if self.wall_clock_sec is not None else "",
            "gpu_count": str(self.gpu_count),
            "gpu_hours": f"{self.gpu_hours:.4f}" if self.gpu_hours is not None else "",
            "exit_code": str(self.exit_code) if self.exit_code is not None else "",
            "log_path": to_project_rel(self.log_path),
            "run_dir": to_project_rel(self.run_dir),
            "status": self.status,
            "error_summary": self.error_summary,
        }

    @classmethod
    def from_dict(cls, row: dict[str, str]) -> "RunMetadata":
        def _dt(s: str) -> datetime | None:
            return datetime.fromisoformat(s) if s else None

        def _float(s: str) -> float | None:
            try:
                return float(s) if s else None
            except ValueError:
                return None

        def _int(s: str) -> int | None:
            try:
                return int(s) if s else None
            except ValueError:
                return None

        return cls(
            model_name=row.get("model_name", ""),
            target_id=row.get("target_id", ""),
            start_time=_dt(row.get("start_time", "")),
            end_time=_dt(row.get("end_time", "")),
            wall_clock_sec=_float(row.get("wall_clock_sec", "")),
            gpu_count=int(row.get("gpu_count") or 1),
            exit_code=_int(row.get("exit_code", "")),
            log_path=Path(row["log_path"]) if row.get("log_path") else None,
            run_dir=Path(row["run_dir"]) if row.get("run_dir") else None,
            status=row.get("status", "unknown"),
            error_summary=row.get("error_summary", ""),
        )


# ─────────────────────────────────────────────────────────────────────────────
# Candidate — 单个设计候选
# ─────────────────────────────────────────────────────────────────────────────

# candidate_id 命名约定：<model_name>__<target_id>__<rank:04d>
# 例如：boltzgen__8dtn_A_B__0001
# 对于 native_reference：native__<target_id>__ref

@dataclass
class Candidate:
    """一个 AI 模型针对给定靶点生成的单个候选。

    source_type / is_reference 用于区分：
        - source_type="generated"  → is_reference=False  → 模型生成候选
        - source_type="native_reference" → is_reference=True → 天然参考（用于 baseline）
    """
    # ── 标识 ─────────────────────────────────────────────────────────────
    candidate_id: str                       # 全局唯一 ID
    target_id: str
    model_name: str

    # ── 来源类型 ─────────────────────────────────────────────────────────
    source_type: str = "generated"          # "generated" | "native_reference"

    @property
    def is_reference(self) -> bool:
        return self.source_type == "native_reference"

    # ── 序列与结构 ───────────────────────────────────────────────────────
    sequence: str = ""                      # 氨基酸单字母序列（大写）
    structure_path: Path | None = None      # .pdb 或 .cif，可能为 None

    # ── 模型内部得分（只用于该模型内部排序，不跨模型比较） ───────────────
    internal_score_name: str = ""           # e.g. "boltz_confidence"
    internal_score_value: float | None = None
    rank_in_model: int = 0                  # 按 internal_score 降序/升序排列后的序号（1-based）

    # ── 来源追踪 ─────────────────────────────────────────────────────────
    source_file: Path | None = None         # 该候选来自哪个文件
    duration_sec: float | None = None       # 该候选的生成时间（如模型提供）
    status: str = "ok"                      # "ok" | "failed" | "empty_seq"
    error_summary: str = ""

    # ── 裁剪版本（从 Target 透传，确保可追溯） ──────────────────────────
    crop_rule_version: str = "none"

    # ── 其他透传字段 ─────────────────────────────────────────────────────
    extra: dict[str, Any] = field(default_factory=dict)

    # ── CSV 字段定义（collection.py 写 CSV 时按此顺序） ─────────────────
    CSV_FIELDS: list[str] = field(default_factory=lambda: [
        "candidate_id", "target_id", "model_name",
        "source_type", "is_reference",
        "sequence", "seq_len",
        "structure_path",
        "internal_score_name", "internal_score_value",
        "rank_in_model",
        "source_file",
        "duration_sec", "status", "error_summary",
        "crop_rule_version",
    ], repr=False, compare=False)

    @property
    def seq_len(self) -> int:
        return len(self.sequence)

    def to_dict(self, root: Path | None = None) -> dict[str, str]:
        """序列化为 CSV 行（路径转相对路径字符串）。"""
        root = root or PROJECT_ROOT
        return {
            "candidate_id": self.candidate_id,
            "target_id": self.target_id,
            "model_name": self.model_name,
            "source_type": self.source_type,
            "is_reference": "1" if self.is_reference else "0",
            "sequence": self.sequence,
            "seq_len": str(self.seq_len),
            "structure_path": to_project_rel(self.structure_path, root),
            "internal_score_name": self.internal_score_name,
            "internal_score_value": (
                f"{self.internal_score_value:.6f}"
                if self.internal_score_value is not None
                else ""
            ),
            "rank_in_model": str(self.rank_in_model),
            "source_file": to_project_rel(self.source_file, root),
            "duration_sec": f"{self.duration_sec:.2f}" if self.duration_sec is not None else "",
            "status": self.status,
            "error_summary": self.error_summary,
            "crop_rule_version": self.crop_rule_version,
        }

    @classmethod
    def from_dict(cls, row: dict[str, str], root: Path | None = None) -> "Candidate":
        """从 CSV 行反序列化 Candidate。"""
        root = root or PROJECT_ROOT

        def _path(s: str) -> Path | None:
            if not s:
                return None
            p = Path(s)
            return (root / p) if not p.is_absolute() else p

        def _float(s: str) -> float | None:
            try:
                return float(s) if s else None
            except ValueError:
                return None

        return cls(
            candidate_id=row.get("candidate_id", ""),
            target_id=row.get("target_id", ""),
            model_name=row.get("model_name", ""),
            source_type=row.get("source_type", "generated"),
            sequence=row.get("sequence", ""),
            structure_path=_path(row.get("structure_path", "")),
            internal_score_name=row.get("internal_score_name", ""),
            internal_score_value=_float(row.get("internal_score_value", "")),
            rank_in_model=int(row.get("rank_in_model") or 0),
            source_file=_path(row.get("source_file", "")),
            duration_sec=_float(row.get("duration_sec", "")),
            status=row.get("status", "ok"),
            error_summary=row.get("error_summary", ""),
            crop_rule_version=row.get("crop_rule_version", "none"),
        )


# ─────────────────────────────────────────────────────────────────────────────
# 批量 IO
# ─────────────────────────────────────────────────────────────────────────────

CANDIDATE_CSV_FIELDS = [
    "candidate_id", "target_id", "model_name",
    "source_type", "is_reference",
    "sequence", "seq_len",
    "structure_path",
    "internal_score_name", "internal_score_value",
    "rank_in_model",
    "source_file",
    "duration_sec", "status", "error_summary",
    "crop_rule_version",
]


def save_candidates(candidates: list[Candidate], path: Path, root: Path | None = None) -> None:
    """将 Candidate 列表写入 CSV。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [c.to_dict(root) for c in candidates]
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CANDIDATE_CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"[candidate] 已保存 {len(candidates)} 个候选到 {path}")


def load_candidates(path: Path, root: Path | None = None) -> list[Candidate]:
    """从 CSV 文件加载 Candidate 列表。"""
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    return [Candidate.from_dict(r, root) for r in rows]
