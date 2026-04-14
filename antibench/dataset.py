"""antibench/dataset.py — Target 数据对象与靶点集加载。

Target 是整个实验数据流的起点，记录一个靶点的所有元信息：
- 原始结构位置（raw）
- 裁剪后结构位置（cropped）
- 表位信息
- 裁剪规则版本（保证所有模型使用同一份裁剪结果）
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from antibench.utils import (
    PROJECT_ROOT,
    PathManager,
    load_config,
    read_json,
    safe_str,
)


# ─────────────────────────────────────────────────────────────────────────────
# 数据对象
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class Epitope:
    """单个靶点的表位描述。"""
    hotspot_string: str                 # e.g. "A:45,A:67"（链:残基号）
    hotspot_count: int
    json_path: Path | None = None       # data/epitopes/<target_id>.json 原始文件


@dataclass
class Target:
    """一个 VHH 设计靶点的完整元信息。

    target_id 命名约定：<pdb_lower>_<antibody_chain>_<antigen_chain>
    例如：8dtn_A_B
    """
    # ── 基本标识 ──────────────────────────────────────────────────────────
    target_id: str                      # e.g. "8dtn_A_B"
    pdb_id: str                         # e.g. "8dtn"（小写，不含 assembly 后缀）
    antibody_chain: str                 # e.g. "A"
    antigen_chain: str                  # e.g. "B"
    antigen_name: str
    antigen_type: str
    antigen_length: int

    # ── 结构路径 ──────────────────────────────────────────────────────────
    # 原始 PDB（未裁剪），可能是复合物或仅抗原链
    raw_complex_path: Path | None = None
    # 裁剪后用于所有模型的标准化抗原结构
    # 若 needs_cropping=False，此路径与 raw_complex_path 相同
    cropped_complex_path: Path | None = None
    # 裁剪规则版本号（写入后续 Candidate / EvalRecord，确保可追溯）
    crop_rule_version: str = "none"

    # ── 裁剪标志 ──────────────────────────────────────────────────────────
    needs_cropping: bool = False
    crop_note: str = ""                 # 原始 dataset_index.csv 中的 note 字段

    # ── 表位 ──────────────────────────────────────────────────────────────
    epitope: Epitope | None = None

    # ── 其他原始字段（透传保存，供 runner 使用） ───────────────────────────
    extra: dict[str, Any] = field(default_factory=dict)

    @property
    def has_epitope(self) -> bool:
        return (
            self.epitope is not None
            and self.epitope.hotspot_count > 0
            and bool(self.epitope.hotspot_string)
        )

    @property
    def effective_complex_path(self) -> Path | None:
        """返回实际用于模型输入的结构路径（裁剪后优先）。"""
        return self.cropped_complex_path or self.raw_complex_path

    def to_dict(self) -> dict[str, str]:
        """序列化为 CSV 友好的 dict（路径转字符串）。"""
        return {
            "target_id": self.target_id,
            "pdb_id": self.pdb_id,
            "antibody_chain": self.antibody_chain,
            "antigen_chain": self.antigen_chain,
            "antigen_name": self.antigen_name,
            "antigen_type": self.antigen_type,
            "antigen_length": str(self.antigen_length),
            "raw_complex_path": str(self.raw_complex_path or ""),
            "cropped_complex_path": str(self.cropped_complex_path or ""),
            "crop_rule_version": self.crop_rule_version,
            "needs_cropping": "1" if self.needs_cropping else "0",
            "crop_note": self.crop_note,
            "epitope_hotspot_string": self.epitope.hotspot_string if self.epitope else "",
            "epitope_hotspot_count": str(self.epitope.hotspot_count) if self.epitope else "0",
            "epitope_json_path": str(self.epitope.json_path or "") if self.epitope else "",
        }


# ─────────────────────────────────────────────────────────────────────────────
# 加载函数
# ─────────────────────────────────────────────────────────────────────────────

def _base_pdb_id(raw: str) -> str:
    """从 '8dtn-assembly1' 提取 '8dtn'。"""
    return safe_str(raw).split("-")[0].lower()


def _load_epitope(target_id: str, epitope_dir: Path) -> Epitope | None:
    """
    尝试从 epitope_dir/<target_id>.json 加载表位信息。
    也检查 <pdb_lower>_<chain1>_<chain2>.json 命名格式。
    """
    candidates = [
        epitope_dir / f"{target_id}.json",
        epitope_dir / f"{target_id.lower()}.json",
    ]
    for ep_path in candidates:
        if ep_path.exists():
            data = read_json(ep_path)
            if not data:
                continue
            return Epitope(
                hotspot_string=safe_str(data.get("hotspot_string") or data.get("hotspots", "")),
                hotspot_count=int(data.get("hotspot_count") or data.get("count") or 0),
                json_path=ep_path,
            )
    return None


def load_targets(
    *,
    dataset_config: dict | None = None,
    path_manager: PathManager | None = None,
    filter_ids: list[str] | None = None,
) -> list[Target]:
    """从 dataset_index.csv 加载所有靶点，构造 Target 对象列表。

    Args:
        dataset_config: configs/dataset.yaml 的内容，None 时自动加载。
        path_manager: PathManager 实例，None 时使用默认项目根。
        filter_ids: 只加载指定 target_id 列表，None 表示全部加载。

    Returns:
        list[Target]，按 dataset_index.csv 行序排列。
    """
    cfg = dataset_config or load_config("dataset")
    pm = path_manager or PathManager()

    index_csv = pm.root / cfg["dataset_index"]
    pdb_dir = pm.root / cfg.get("pdb_dir", "data/pdbs")
    epitope_dir = pm.root / cfg.get("epitope_dir", "data/epitopes")
    crop_cfg = cfg.get("cropping", {})
    crop_threshold = int(crop_cfg.get("length_threshold", 300))
    crop_rule_version = safe_str(crop_cfg.get("rule_version", "v1_radius12A"))

    if not index_csv.exists():
        raise FileNotFoundError(f"dataset_index.csv 不存在: {index_csv}")

    targets: list[Target] = []
    seen_ids: set[str] = set()

    with index_csv.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # 兼容中英文列名
            raw_pdb_id = safe_str(row.get("pdb_id") or row.get("PDB ID") or "")
            antibody_chain = safe_str(row.get("antibody_chain") or "")
            antigen_chain = safe_str(row.get("antigen_chain") or "")

            if not raw_pdb_id or not antibody_chain or not antigen_chain:
                continue

            pdb_id = _base_pdb_id(raw_pdb_id)
            target_id = f"{pdb_id}_{antibody_chain}_{antigen_chain}"

            # 去重（dataset_index.csv 中同一靶点可能有多行）
            if target_id in seen_ids:
                continue
            seen_ids.add(target_id)

            # 按 filter_ids 过滤
            if filter_ids and target_id not in filter_ids:
                continue

            antigen_length = int(safe_str(row.get("antigen_length") or "0") or 0)
            needs_cropping = antigen_length > crop_threshold

            # ── 结构路径 ────────────────────────────────────────────────
            raw_pdb = pdb_dir / f"{pdb_id}.pdb"
            raw_complex_path = raw_pdb if raw_pdb.exists() else None

            # 裁剪后路径（data_prep 阶段写入，此时可能还不存在）
            cropped_dir_for_target = pm.cropped_dir / target_id
            cropped_pdb = cropped_dir_for_target / "antigen_cropped.pdb"
            cropped_complex_path = cropped_pdb if cropped_pdb.exists() else None

            # ── 表位 ────────────────────────────────────────────────────
            epitope = _load_epitope(target_id, epitope_dir)

            target = Target(
                target_id=target_id,
                pdb_id=pdb_id,
                antibody_chain=antibody_chain,
                antigen_chain=antigen_chain,
                antigen_name=safe_str(row.get("antigen_name") or ""),
                antigen_type=safe_str(row.get("antigen_type") or ""),
                antigen_length=antigen_length,
                raw_complex_path=raw_complex_path,
                cropped_complex_path=cropped_complex_path,
                crop_rule_version=crop_rule_version if needs_cropping else "none",
                needs_cropping=needs_cropping,
                crop_note=safe_str(row.get("note") or ""),
                epitope=epitope,
                extra={k: safe_str(v) for k, v in row.items()},
            )
            targets.append(target)

    return targets


def save_prepared_index(targets: list[Target], out_path: Path) -> None:
    """将 Target 列表序列化为 CSV（data_prep 阶段完成后调用）。"""
    if not targets:
        return
    rows = [t.to_dict() for t in targets]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"[dataset] 已保存 {len(targets)} 个靶点到 {out_path}")
