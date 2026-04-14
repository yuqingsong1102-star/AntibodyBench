#!/usr/bin/env python3
"""
AntibodyBench — 论文级结果分析与可视化

用法:
  python scripts/analyze_results.py \
    --eval-dirs outputs/evaluation/chai_boltzgen_eval outputs/evaluation/chai_mber_eval \
    --model-labels boltzgen mber-open \
    --out-dir outputs/analysis/thesis_figures

支持单模型（只传一个 --eval-dirs）或多模型对比。
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

# ── 论文级画图风格 ────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# 色盲友好色板 (Okabe-Ito)
MODEL_COLORS = {
    "boltzgen":  "#0072B2",
    "mber-open": "#D55E00",
    "germinal":  "#009E73",
    "RFantibody":"#CC79A7",
}
DEFAULT_COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#F0E442", "#56B4E9"]

# ── 数据加载 ─────────────────────────────────────────────────

def load_eval_merged(eval_dir: Path) -> list[dict]:
    """Load eval_merged.csv (new 3-step format) or fall back to external_eval_seed_level.csv (legacy)."""
    for name in ("eval_final.csv", "eval_merged.csv", "external_eval_seed_level.csv"):
        csv_path = eval_dir / name
        if csv_path.exists():
            with open(csv_path) as f:
                rows = list(csv.DictReader(f))
            # Normalise field names: map new names to legacy aliases used by figures
            normalised = []
            for r in rows:
                r.setdefault("external_iptm", r.get("chai1_iptm", ""))
                r.setdefault("external_ptm", r.get("chai1_ptm", ""))
                r.setdefault("external_plddt_binder", r.get("chai1_plddt_binder", ""))
                normalised.append(r)
            return normalised
    print(f"[WARN] No eval CSV found in {eval_dir}, skipping", file=sys.stderr)
    return []


def safe_float(val: str | None, default: float = float("nan")) -> float:
    if not val or val.strip() == "":
        return default
    try:
        return float(val)
    except ValueError:
        return default


# ── Top-K 归一化选择 ──────────────────────────────────────────

TOP_K_RANKING_METRICS = {
    "chai1_iptm": "external_iptm",
    "pdockq2": "pdockq2",
    "composite": None,  # 特殊: 加权组合
}


def select_topk(
    all_data: dict[str, list[dict]],
    top_k: int,
    rank_by: str = "chai1_iptm",
) -> dict[str, list[dict]]:
    """按 (model, sample_id) 分组，组内排序取 top-K 候选。

    Args:
        all_data: {model_label: [row, ...]}
        top_k: 每个 (model, target) 保留的候选数。0 表示不限。
        rank_by: 排序指标名 (chai1_iptm / pdockq2 / composite)

    Returns:
        与 all_data 同结构，但每组只保留 top-K。
    """
    if top_k <= 0:
        return all_data

    result: dict[str, list[dict]] = {}
    for model, rows in all_data.items():
        # 按 sample_id 分组
        groups: dict[str, list[dict]] = defaultdict(list)
        for r in rows:
            groups[r["sample_id"]].append(r)

        selected: list[dict] = []
        for sid, group in groups.items():
            # 排序: 按指定指标降序 (数值越大越好)
            if rank_by == "composite":
                # 组合排序: 0.5*iPTM + 0.5*pDockQ2 (归一化到 [0,1])
                def _composite_score(r: dict) -> float:
                    iptm = safe_float(r.get("external_iptm", r.get("chai1_iptm", "")))
                    pdockq = safe_float(r.get("pdockq2", ""))
                    s = 0.0
                    n = 0
                    if not np.isnan(iptm):
                        s += iptm
                        n += 1
                    if not np.isnan(pdockq):
                        s += pdockq
                        n += 1
                    return s / n if n > 0 else -1.0
                group.sort(key=_composite_score, reverse=True)
            else:
                field = TOP_K_RANKING_METRICS.get(rank_by, rank_by)
                if field is None:
                    field = "external_iptm"
                group.sort(
                    key=lambda r: safe_float(r.get(field, r.get(rank_by, ""))),
                    reverse=True,
                )

            selected.extend(group[:top_k])

        result[model] = selected

    # 打印 top-K 选择摘要
    print(f"\n[Top-K] rank_by={rank_by}, K={top_k}")
    for model, rows in result.items():
        targets = set(r["sample_id"] for r in rows)
        print(f"  {model}: {len(rows)} candidates across {len(targets)} targets "
              f"(avg {len(rows)/len(targets):.1f}/target)" if targets else f"  {model}: 0 candidates")

    return result


# ── 序列分析工具 ──────────────────────────────────────────────

# VHH 典型 CDR 长度范围 (IMGT numbering approximation)
VHH_CDR_PATTERNS = {
    # 简化版: 基于框架保守基序定位 CDR
    "FR1_end":  r"SCAAS",       # FR1 终止于 Cys-Ala-Ala-Ser 附近
    "CDR1_end": r"W[FYVR]RQAPGK",  # CDR1 终止于 Trp 后接 FR2
    "CDR2_start": r"[AVILM][A-Z]{1,3}[A-Z]G[A-Z]{0,3}[STDN]",
    "FR3_start": r"[RK]F[TI]ISR",   # FR3 开始 RF[TI]ISR
    "CDR3_start": r"YYC",       # CDR3 前序 Tyr-Tyr-Cys
    "FR4_start": r"WG[QK]G[TLMS]",  # FR4 开始 WG.G
}


def extract_cdr3(seq: str) -> str | None:
    """粗略提取 CDR3 区域 (YYC...WGxG 之间)"""
    m1 = re.search(r"YYC(.+?)WG[QK]G", seq)
    if m1:
        return m1.group(1)
    return None


def pairwise_identity(seqs: list[str]) -> float:
    """计算一组序列的平均 pairwise identity"""
    if len(seqs) < 2:
        return 1.0
    total, count = 0.0, 0
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            s1, s2 = seqs[i], seqs[j]
            minlen = min(len(s1), len(s2))
            if minlen == 0:
                continue
            matches = sum(a == b for a, b in zip(s1[:minlen], s2[:minlen]))
            total += matches / minlen
            count += 1
    return total / count if count else 1.0


# ── 图表生成 ─────────────────────────────────────────────────

def fig_iptm_by_target(all_data: dict[str, list[dict]], out_dir: Path):
    """每个靶点的 Chai-1 iPTM 分值，按模型分组"""
    models = list(all_data.keys())
    # 收集所有 sample_id
    sample_ids = sorted({r["sample_id"] for rows in all_data.values() for r in rows})

    if len(models) == 1:
        # 单模型: 水平条形图，按 iPTM 排序
        model = models[0]
        data = {r["sample_id"]: safe_float(r["external_iptm"]) for r in all_data[model]}
        sorted_items = sorted(data.items(), key=lambda x: x[1])
        names, vals = zip(*sorted_items)

        fig, ax = plt.subplots(figsize=(7, max(5, len(names) * 0.35)))
        colors = [plt.cm.RdYlGn(v / max(0.8, max(vals))) for v in vals]
        bars = ax.barh(range(len(names)), vals, color=colors, edgecolor="white", linewidth=0.5)
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels([n.replace("_", " ") for n in names], fontsize=9)
        ax.set_xlabel("Chai-1 iPTM")
        ax.set_title(f"Interface Predicted TM-score by Target ({model})")
        ax.axvline(x=0.6, color="#888", linestyle="--", linewidth=0.8, alpha=0.6)
        ax.text(0.61, len(names) - 0.5, "iPTM=0.6", fontsize=8, color="#666")
        ax.set_xlim(0, min(1.0, max(vals) * 1.2))
        for i, v in enumerate(vals):
            ax.text(v + 0.01, i, f"{v:.2f}", va="center", fontsize=8, color="#333")
    else:
        # 多模型: 分组条形图
        fig, ax = plt.subplots(figsize=(max(8, len(sample_ids) * 0.6), 5))
        x = np.arange(len(sample_ids))
        width = 0.8 / len(models)
        for i, model in enumerate(models):
            data = {r["sample_id"]: safe_float(r["external_iptm"]) for r in all_data[model]}
            vals = [data.get(sid, float("nan")) for sid in sample_ids]
            color = MODEL_COLORS.get(model, DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
            ax.bar(x + i * width - 0.4 + width / 2, vals, width * 0.9, label=model, color=color, alpha=0.85)
        ax.set_xticks(x)
        ax.set_xticklabels([s.replace("_", "\n") for s in sample_ids], fontsize=8, rotation=45, ha="right")
        ax.set_ylabel("Chai-1 iPTM")
        ax.set_title("Interface Predicted TM-score by Target")
        ax.axhline(y=0.6, color="#888", linestyle="--", linewidth=0.8, alpha=0.6)
        ax.legend(loc="upper left")

    fig.tight_layout()
    fig.savefig(out_dir / "fig01_iptm_by_target.png")
    fig.savefig(out_dir / "fig01_iptm_by_target.pdf")
    plt.close(fig)
    print(f"  [OK] fig01_iptm_by_target")


def fig_pdockq2_by_target(all_data: dict[str, list[dict]], out_dir: Path):
    """每个靶点的 pDockQ2"""
    models = list(all_data.keys())
    sample_ids = sorted({r["sample_id"] for rows in all_data.values() for r in rows})

    if len(models) == 1:
        model = models[0]
        data = {r["sample_id"]: safe_float(r["pdockq2"]) for r in all_data[model]}
        sorted_items = sorted(data.items(), key=lambda x: x[1])
        names, vals = zip(*sorted_items)

        fig, ax = plt.subplots(figsize=(7, max(5, len(names) * 0.35)))
        colors = [plt.cm.RdYlGn(v / max(0.6, max(vals))) for v in vals]
        bars = ax.barh(range(len(names)), vals, color=colors, edgecolor="white", linewidth=0.5)
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels([n.replace("_", " ") for n in names], fontsize=9)
        ax.set_xlabel("pDockQ2")
        ax.set_title(f"Predicted Interface Quality (pDockQ2) by Target ({model})")
        ax.axvline(x=0.23, color="#888", linestyle="--", linewidth=0.8, alpha=0.6)
        ax.text(0.24, len(names) - 0.5, "pDockQ2=0.23", fontsize=8, color="#666")
        ax.set_xlim(0, min(1.0, max(vals) * 1.2))
        for i, v in enumerate(vals):
            ax.text(v + 0.01, i, f"{v:.2f}", va="center", fontsize=8, color="#333")
    else:
        fig, ax = plt.subplots(figsize=(max(8, len(sample_ids) * 0.6), 5))
        x = np.arange(len(sample_ids))
        width = 0.8 / len(models)
        for i, model in enumerate(models):
            data = {r["sample_id"]: safe_float(r["pdockq2"]) for r in all_data[model]}
            vals = [data.get(sid, float("nan")) for sid in sample_ids]
            color = MODEL_COLORS.get(model, DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
            ax.bar(x + i * width - 0.4 + width / 2, vals, width * 0.9, label=model, color=color, alpha=0.85)
        ax.set_xticks(x)
        ax.set_xticklabels([s.replace("_", "\n") for s in sample_ids], fontsize=8, rotation=45, ha="right")
        ax.set_ylabel("pDockQ2")
        ax.set_title("Predicted Interface Quality (pDockQ2) by Target")
        ax.axhline(y=0.23, color="#888", linestyle="--", linewidth=0.8, alpha=0.6)
        ax.legend(loc="upper left")

    fig.tight_layout()
    fig.savefig(out_dir / "fig02_pdockq2_by_target.png")
    fig.savefig(out_dir / "fig02_pdockq2_by_target.pdf")
    plt.close(fig)
    print(f"  [OK] fig02_pdockq2_by_target")


def fig_metric_distributions(all_data: dict[str, list[dict]], out_dir: Path):
    """多指标分布箱线图/小提琴图"""
    metrics = [
        ("external_iptm", "Chai-1 iPTM", (0, 1)),
        ("pdockq2", "pDockQ2", (0, 0.8)),
        ("interface_bsa", "Interface BSA (Å²)", None),
        ("avg_interface_plddt", "Avg Interface pLDDT", (40, 100)),
        ("external_plddt_binder", "Binder pLDDT", (50, 100)),
        ("interface_contact_pair_count", "Contact Pairs", None),
    ]
    models = list(all_data.keys())
    n_metrics = len(metrics)
    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    axes = axes.flatten()

    for idx, (col, label, ylim) in enumerate(metrics):
        ax = axes[idx]
        data_per_model = []
        labels = []
        for model in models:
            vals = [safe_float(r[col]) for r in all_data[model]]
            vals = [v for v in vals if not np.isnan(v)]
            data_per_model.append(vals)
            labels.append(model)

        if len(models) == 1:
            parts = ax.violinplot(data_per_model, showmeans=True, showmedians=True)
            for pc in parts["bodies"]:
                color = MODEL_COLORS.get(models[0], DEFAULT_COLORS[0])
                pc.set_facecolor(color)
                pc.set_alpha(0.6)
            ax.set_xticks([1])
            ax.set_xticklabels(labels)
        else:
            positions = range(1, len(models) + 1)
            parts = ax.violinplot(data_per_model, positions=positions, showmeans=True, showmedians=True)
            for i, pc in enumerate(parts["bodies"]):
                color = MODEL_COLORS.get(models[i], DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
                pc.set_facecolor(color)
                pc.set_alpha(0.6)
            ax.set_xticks(positions)
            ax.set_xticklabels(labels)

        ax.set_title(label)
        if ylim:
            ax.set_ylim(ylim)

        # 标注中位数
        for i, vals in enumerate(data_per_model):
            if vals:
                med = np.median(vals)
                ax.text(i + 1, med, f" {med:.2f}", fontsize=8, va="bottom", color="#333")

    fig.suptitle("Metric Distributions", fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(out_dir / "fig03_metric_distributions.png")
    fig.savefig(out_dir / "fig03_metric_distributions.pdf")
    plt.close(fig)
    print(f"  [OK] fig03_metric_distributions")


def fig_iptm_vs_pdockq2(all_data: dict[str, list[dict]], out_dir: Path):
    """iPTM vs pDockQ2 散点图"""
    models = list(all_data.keys())
    fig, ax = plt.subplots(figsize=(7, 6))

    for i, model in enumerate(models):
        color = MODEL_COLORS.get(model, DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
        for r in all_data[model]:
            iptm = safe_float(r["external_iptm"])
            pdockq = safe_float(r["pdockq2"])
            bsa = safe_float(r["interface_bsa"], 400)
            if np.isnan(iptm) or np.isnan(pdockq):
                continue
            size = max(20, min(200, bsa / 10))
            ax.scatter(iptm, pdockq, s=size, color=color, alpha=0.7, edgecolors="white", linewidth=0.5)
            # 标注样本名
            sid = r["sample_id"].split("_")[0]  # PDB code only
            ax.annotate(sid, (iptm, pdockq), fontsize=6, alpha=0.7,
                        xytext=(3, 3), textcoords="offset points")

    # 参考线
    ax.axvline(x=0.6, color="#ccc", linestyle="--", linewidth=0.8)
    ax.axhline(y=0.23, color="#ccc", linestyle="--", linewidth=0.8)

    # 象限标注
    ax.text(0.8, 0.6, "Best", fontsize=9, color="#2a9d2a", ha="center", alpha=0.6)
    ax.text(0.2, 0.05, "Worst", fontsize=9, color="#cc3333", ha="center", alpha=0.6)

    ax.set_xlabel("Chai-1 iPTM")
    ax.set_ylabel("pDockQ2")
    ax.set_title("Interface Quality: iPTM vs pDockQ2")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, max(0.6, max(safe_float(r["pdockq2"]) for rows in all_data.values() for r in rows) * 1.1))

    # 图例：模型颜色 + 大小图例
    if len(models) > 1:
        handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=MODEL_COLORS.get(m, DEFAULT_COLORS[i]),
                              markersize=8, label=m) for i, m in enumerate(models)]
        ax.legend(handles=handles, loc="upper left")

    # 大小图例
    for bsa_val, label in [(400, "BSA=400"), (800, "800"), (1200, "1200")]:
        ax.scatter([], [], s=bsa_val / 10, c="gray", alpha=0.4, label=f"BSA≈{bsa_val}")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.8)

    fig.tight_layout()
    fig.savefig(out_dir / "fig04_iptm_vs_pdockq2.png")
    fig.savefig(out_dir / "fig04_iptm_vs_pdockq2.pdf")
    plt.close(fig)
    print(f"  [OK] fig04_iptm_vs_pdockq2")


def fig_success_rate_curve(all_data: dict[str, list[dict]], out_dir: Path):
    """不同 iPTM 阈值下的"成功率"曲线 (cumulative distribution)"""
    models = list(all_data.keys())
    thresholds = np.arange(0.05, 0.85, 0.025)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    for i, model in enumerate(models):
        color = MODEL_COLORS.get(model, DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
        iptm_vals = sorted([safe_float(r["external_iptm"]) for r in all_data[model]
                            if not np.isnan(safe_float(r["external_iptm"]))])
        pdockq_vals = sorted([safe_float(r["pdockq2"]) for r in all_data[model]
                              if not np.isnan(safe_float(r["pdockq2"]))])
        n_iptm = len(iptm_vals)
        n_pdockq = len(pdockq_vals)

        # iPTM success curve
        rates_iptm = [sum(1 for v in iptm_vals if v >= t) / n_iptm for t in thresholds] if n_iptm else [0] * len(thresholds)
        ax1.plot(thresholds, rates_iptm, color=color, linewidth=2, label=model, marker=".", markersize=3)

        # pDockQ2 success curve
        rates_pdockq = [sum(1 for v in pdockq_vals if v >= t) / n_pdockq for t in thresholds] if n_pdockq else [0] * len(thresholds)
        ax2.plot(thresholds, rates_pdockq, color=color, linewidth=2, label=model, marker=".", markersize=3)

    ax1.set_xlabel("iPTM Threshold")
    ax1.set_ylabel("Fraction of Targets ≥ Threshold")
    ax1.set_title("Cumulative Success Rate (iPTM)")
    ax1.axvline(x=0.6, color="#aaa", linestyle=":", linewidth=0.8)
    ax1.set_xlim(0.05, 0.8)
    ax1.set_ylim(0, 1.05)
    ax1.yaxis.set_major_formatter(mticker.PercentFormatter(1.0))
    ax1.legend()

    ax2.set_xlabel("pDockQ2 Threshold")
    ax2.set_ylabel("Fraction of Targets ≥ Threshold")
    ax2.set_title("Cumulative Success Rate (pDockQ2)")
    ax2.axvline(x=0.23, color="#aaa", linestyle=":", linewidth=0.8)
    ax2.set_xlim(0.05, 0.6)
    ax2.set_ylim(0, 1.05)
    ax2.yaxis.set_major_formatter(mticker.PercentFormatter(1.0))
    ax2.legend()

    fig.tight_layout()
    fig.savefig(out_dir / "fig05_success_rate_curve.png")
    fig.savefig(out_dir / "fig05_success_rate_curve.pdf")
    plt.close(fig)
    print(f"  [OK] fig05_success_rate_curve")


def fig_interface_heatmap(all_data: dict[str, list[dict]], out_dir: Path):
    """靶点 × 指标 热力图"""
    metrics = [
        ("external_iptm", "iPTM"),
        ("pdockq2", "pDockQ2"),
        ("avg_interface_plddt", "Interface\npLDDT"),
        ("interface_bsa", "Interface\nBSA"),
        ("interface_contact_pair_count", "Contact\nPairs"),
        ("hotspot_contact_count", "Hotspot\nContacts"),
    ]
    models = list(all_data.keys())

    for model in models:
        rows = all_data[model]
        sample_ids = sorted({r["sample_id"] for r in rows})
        data_matrix = np.zeros((len(sample_ids), len(metrics)))

        for si, sid in enumerate(sample_ids):
            row = next((r for r in rows if r["sample_id"] == sid), None)
            if not row:
                continue
            for mi, (col, _) in enumerate(metrics):
                data_matrix[si, mi] = safe_float(row[col])

        # 按列归一化 (min-max)
        norm_matrix = np.zeros_like(data_matrix)
        for mi in range(data_matrix.shape[1]):
            col_data = data_matrix[:, mi]
            vmin, vmax = np.nanmin(col_data), np.nanmax(col_data)
            if vmax > vmin:
                norm_matrix[:, mi] = (col_data - vmin) / (vmax - vmin)
            else:
                norm_matrix[:, mi] = 0.5

        fig, ax = plt.subplots(figsize=(8, max(5, len(sample_ids) * 0.35)))
        im = ax.imshow(norm_matrix, cmap="RdYlGn", aspect="auto", vmin=0, vmax=1)

        ax.set_xticks(range(len(metrics)))
        ax.set_xticklabels([m[1] for m in metrics], fontsize=9)
        ax.set_yticks(range(len(sample_ids)))
        ax.set_yticklabels([s.replace("_", " ") for s in sample_ids], fontsize=9)

        # 标注原始数值
        for si in range(len(sample_ids)):
            for mi in range(len(metrics)):
                val = data_matrix[si, mi]
                txt = f"{val:.0f}" if val > 10 else f"{val:.2f}"
                color = "white" if norm_matrix[si, mi] < 0.3 or norm_matrix[si, mi] > 0.85 else "black"
                ax.text(mi, si, txt, ha="center", va="center", fontsize=7, color=color)

        suffix = f" ({model})" if len(models) > 1 else ""
        ax.set_title(f"Per-Target Interface Metrics{suffix}")
        fig.colorbar(im, ax=ax, label="Normalized Score", shrink=0.6)

        fig.tight_layout()
        tag = f"_{model}" if len(models) > 1 else ""
        fig.savefig(out_dir / f"fig06_interface_heatmap{tag}.png")
        fig.savefig(out_dir / f"fig06_interface_heatmap{tag}.pdf")
        plt.close(fig)
        print(f"  [OK] fig06_interface_heatmap{tag}")


def fig_sequence_analysis(all_data: dict[str, list[dict]], out_dir: Path):
    """序列级分析: CDR3 长度、序列长度分布"""
    models = list(all_data.keys())

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))

    # Panel 1: 序列长度分布
    ax = axes[0]
    for i, model in enumerate(models):
        lens = [int(r.get("seq_length") or r.get("pred_seq_len") or 0) for r in all_data[model]
                if (r.get("seq_length") or r.get("pred_seq_len", ""))]
        if not lens:
            continue
        color = MODEL_COLORS.get(model, DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
        ax.hist(lens, bins=range(min(lens) - 2, max(lens) + 3), alpha=0.6, color=color, label=model, edgecolor="white")
    ax.set_xlabel("VHH Sequence Length")
    ax.set_ylabel("Count")
    ax.set_title("VHH Length Distribution")
    ax.legend()

    # Panel 2: CDR3 长度分布
    ax = axes[1]
    for i, model in enumerate(models):
        cdr3_lens = []
        for r in all_data[model]:
            cdr3 = extract_cdr3(r["pred_sequence"])
            if cdr3:
                cdr3_lens.append(len(cdr3))
        color = MODEL_COLORS.get(model, DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
        if cdr3_lens:
            ax.hist(cdr3_lens, bins=range(min(cdr3_lens) - 1, max(cdr3_lens) + 2),
                    alpha=0.6, color=color, label=model, edgecolor="white")
    ax.set_xlabel("CDR3 Length (residues)")
    ax.set_ylabel("Count")
    ax.set_title("CDR3 Length Distribution")
    ax.legend()

    # Panel 3: Developability 指标
    ax = axes[2]
    # New schema uses count fields; old schema used flag fields
    dev_metrics = ["free_cys_count", "deamidation_motif_count", "isomerization_motif_count", "glycosylation_motif_count"]
    dev_fallback = ["free_cys_flag", "deamidation_risk_flag", "isomerization_risk_flag", "glycosylation_motif_flag"]
    dev_labels = ["Free Cys", "Deamidation", "Isomerization", "Glycosylation"]
    x = np.arange(len(dev_metrics))
    width = 0.8 / max(1, len(models))
    for i, model in enumerate(models):
        rates = []
        for dm, fb in zip(dev_metrics, dev_fallback):
            vals = []
            for r in all_data[model]:
                v = r.get(dm, r.get(fb, ""))
                if v != "":
                    vals.append(1 if float(v) > 0 else 0)
            rates.append(sum(vals) / len(vals) if vals else 0)
        color = MODEL_COLORS.get(model, DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
        ax.bar(x + i * width - 0.4 + width / 2, rates, width * 0.9, label=model, color=color, alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(dev_labels, fontsize=9)
    ax.set_ylabel("Fraction with Liability")
    ax.set_title("Sequence Liability Flags")
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(1.0))
    ax.legend()

    fig.tight_layout()
    fig.savefig(out_dir / "fig07_sequence_analysis.png")
    fig.savefig(out_dir / "fig07_sequence_analysis.pdf")
    plt.close(fig)
    print(f"  [OK] fig07_sequence_analysis")


def fig_model_comparison_radar(all_data: dict[str, list[dict]], out_dir: Path):
    """多模型雷达图对比 (仅当有 ≥2 个模型时)"""
    models = list(all_data.keys())
    if len(models) < 2:
        print(f"  [SKIP] fig08_radar (need ≥2 models, have {len(models)})")
        return

    # 定义雷达图维度 (越高越好的指标)
    dims = [
        ("external_iptm", "iPTM"),
        ("pdockq2", "pDockQ2"),
        ("avg_interface_plddt", "Interface pLDDT"),
        ("interface_bsa", "Interface BSA"),
        ("external_plddt_binder", "Binder pLDDT"),
    ]

    # 计算每个模型各维度的中位数
    raw_vals = {}
    for model in models:
        raw_vals[model] = []
        for col, _ in dims:
            vals = [safe_float(r[col]) for r in all_data[model]]
            vals = [v for v in vals if not np.isnan(v)]
            raw_vals[model].append(np.median(vals) if vals else 0)

    # 归一化 (相对于所有模型的 max)
    n_dims = len(dims)
    norm_vals = {}
    for model in models:
        norm_vals[model] = []
    for di in range(n_dims):
        max_val = max(raw_vals[m][di] for m in models)
        for model in models:
            norm_vals[model].append(raw_vals[model][di] / max_val if max_val > 0 else 0)

    # 画图
    angles = np.linspace(0, 2 * np.pi, n_dims, endpoint=False).tolist()
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=(7, 7), subplot_kw=dict(polar=True))
    for i, model in enumerate(models):
        vals = norm_vals[model] + norm_vals[model][:1]
        color = MODEL_COLORS.get(model, DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
        ax.plot(angles, vals, 'o-', linewidth=2, label=model, color=color)
        ax.fill(angles, vals, alpha=0.15, color=color)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels([d[1] for d in dims], fontsize=10)
    ax.set_ylim(0, 1.1)
    ax.set_title("Multi-Metric Model Comparison", y=1.1)
    ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1))

    fig.tight_layout()
    fig.savefig(out_dir / "fig08_model_comparison_radar.png")
    fig.savefig(out_dir / "fig08_model_comparison_radar.pdf")
    plt.close(fig)
    print(f"  [OK] fig08_model_comparison_radar")


def fig_head_to_head(all_data: dict[str, list[dict]], out_dir: Path):
    """模型对比: 每个靶点谁赢 (仅 ≥2 个模型时)"""
    models = list(all_data.keys())
    if len(models) < 2:
        print(f"  [SKIP] fig09_head_to_head (need ≥2 models)")
        return

    sample_ids = sorted({r["sample_id"] for rows in all_data.values() for r in rows})

    # 每个 sample 每个 model 的 best iPTM
    model_iptm = {}
    for model in models:
        model_iptm[model] = {}
        for r in all_data[model]:
            sid = r["sample_id"]
            val = safe_float(r["external_iptm"])
            if sid not in model_iptm[model] or val > model_iptm[model][sid]:
                model_iptm[model][sid] = val

    # 统计胜率
    wins = {m: 0 for m in models}
    ties = 0
    for sid in sample_ids:
        best_val = -1
        best_model = None
        for m in models:
            v = model_iptm[m].get(sid, float("-inf"))
            if v > best_val:
                best_val = v
                best_model = m
        if best_model:
            wins[best_model] += 1

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw={"width_ratios": [3, 1]})

    # Left: grouped bar chart
    x = np.arange(len(sample_ids))
    width = 0.8 / len(models)
    for i, model in enumerate(models):
        vals = [model_iptm[model].get(sid, 0) for sid in sample_ids]
        color = MODEL_COLORS.get(model, DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
        ax1.bar(x + i * width - 0.4 + width / 2, vals, width * 0.9, label=model, color=color, alpha=0.85)
    ax1.set_xticks(x)
    ax1.set_xticklabels([s.split("_")[0] for s in sample_ids], fontsize=8, rotation=60, ha="right")
    ax1.set_ylabel("Best Chai-1 iPTM")
    ax1.set_title("Head-to-Head: iPTM per Target")
    ax1.legend()

    # Right: win rate pie
    colors = [MODEL_COLORS.get(m, DEFAULT_COLORS[i]) for i, m in enumerate(models)]
    ax2.pie([wins[m] for m in models], labels=[f"{m}\n({wins[m]}/{len(sample_ids)})" for m in models],
            colors=colors, autopct="%1.0f%%", startangle=90)
    ax2.set_title("Win Rate")

    fig.tight_layout()
    fig.savefig(out_dir / "fig09_head_to_head.png")
    fig.savefig(out_dir / "fig09_head_to_head.pdf")
    plt.close(fig)
    print(f"  [OK] fig09_head_to_head")


# ── 新图表: 帕累托 Hit Rate / 多样性 / 表位覆盖 ─────────────

# 帕累托判据阈值 (可调)
PARETO_CRITERIA = {
    "A_iptm":  ("external_iptm",    0.85, ">="),
    "B_ddG":   ("approx_dG",       -30.0, "<="),
    "C_bsa":   ("interface_bsa",   800.0, ">="),  # BSA 作为 Sc 的替代
}


def _pareto_pass(row: dict, criterion: str, criteria: dict = PARETO_CRITERIA) -> bool:
    """Check if a row passes one Pareto criterion."""
    field, threshold, direction = criteria[criterion]
    val = safe_float(row.get(field, ""))
    if np.isnan(val):
        return False
    return val >= threshold if direction == ">=" else val <= threshold


def fig_pareto_hit_rate(all_data: dict[str, list[dict]], out_dir: Path):
    """帕累托命中率: A(iPTM>=0.85), B(ddG<=-30), C(BSA>=800), Hit Rate = |A∪B∪C|/N"""
    models = list(all_data.keys())

    # Per model: compute hit rates by target
    model_stats = {}
    for model in models:
        per_target: dict[str, dict] = defaultdict(lambda: {"total": 0, "A": 0, "B": 0, "C": 0, "union": 0})
        for r in all_data[model]:
            sid = r["sample_id"]
            per_target[sid]["total"] += 1
            a = _pareto_pass(r, "A_iptm")
            b = _pareto_pass(r, "B_ddG")
            c = _pareto_pass(r, "C_bsa")
            if a: per_target[sid]["A"] += 1
            if b: per_target[sid]["B"] += 1
            if c: per_target[sid]["C"] += 1
            if a or b or c: per_target[sid]["union"] += 1
        model_stats[model] = per_target

    all_targets = sorted({sid for st in model_stats.values() for sid in st})

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Left: grouped bar chart — hit rate per target
    ax = axes[0]
    x = np.arange(len(all_targets))
    n_models = len(models)
    width = 0.8 / max(1, n_models)
    for i, model in enumerate(models):
        rates = []
        for sid in all_targets:
            st = model_stats[model].get(sid, {"total": 0, "union": 0})
            rates.append(st["union"] / st["total"] * 100 if st["total"] > 0 else 0)
        color = MODEL_COLORS.get(model, DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
        ax.bar(x + i * width - 0.4 + width / 2, rates, width * 0.9, label=model, color=color, alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels([s.split("_")[0] for s in all_targets], fontsize=8, rotation=60, ha="right")
    ax.set_ylabel("Hit Rate (%)")
    ax.set_title("Pareto Hit Rate per Target\n(A: iPTM≥0.85 ∪ B: ΔG≤−30 ∪ C: BSA≥800)")
    ax.legend(fontsize=9)

    # Right: stacked bar for overall breakdown per model
    ax = axes[1]
    labels_crit = ["A: iPTM≥0.85", "B: ΔG≤−30", "C: BSA≥800", "Hit (A∪B∪C)"]
    x2 = np.arange(len(models))
    bar_data = {lbl: [] for lbl in labels_crit}
    for model in models:
        total = sum(st["total"] for st in model_stats[model].values())
        for key, lbl in [("A", labels_crit[0]), ("B", labels_crit[1]), ("C", labels_crit[2]), ("union", labels_crit[3])]:
            count = sum(st[key] for st in model_stats[model].values())
            bar_data[lbl].append(count / total * 100 if total > 0 else 0)
    bar_colors = ["#0072B2", "#D55E00", "#009E73", "#333333"]
    for j, lbl in enumerate(labels_crit):
        ax.bar(x2 + j * 0.18 - 0.27, bar_data[lbl], 0.16, label=lbl, color=bar_colors[j], alpha=0.85)
    ax.set_xticks(x2)
    ax.set_xticklabels(models, fontsize=10)
    ax.set_ylabel("Candidates (%)")
    ax.set_title("Overall Pareto Criteria Pass Rates")
    ax.legend(fontsize=8, loc="upper right")

    fig.tight_layout()
    fig.savefig(out_dir / "fig10_pareto_hit_rate.png")
    fig.savefig(out_dir / "fig10_pareto_hit_rate.pdf")
    plt.close(fig)
    print(f"  [OK] fig10_pareto_hit_rate")


def fig_diversity(all_data: dict[str, list[dict]], out_dir: Path):
    """序列多样性: cluster 数量 & 最大簇占比"""
    models = list(all_data.keys())
    # Check if diversity fields exist
    has_div = any(r.get("group_cluster_count", "") for rows in all_data.values() for r in rows)
    if not has_div:
        print(f"  [SKIP] fig11_diversity (no clustering data)")
        return

    all_targets = sorted({r["sample_id"] for rows in all_data.values() for r in rows})

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: cluster count per target
    ax = axes[0]
    x = np.arange(len(all_targets))
    width = 0.8 / max(1, len(models))
    for i, model in enumerate(models):
        target_clusters = {}
        for r in all_data[model]:
            sid = r["sample_id"]
            c = safe_float(r.get("group_cluster_count", ""))
            if not np.isnan(c) and c > 0:
                target_clusters[sid] = int(c)
        vals = [target_clusters.get(sid, 0) for sid in all_targets]
        color = MODEL_COLORS.get(model, DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
        ax.bar(x + i * width - 0.4 + width / 2, vals, width * 0.9, label=model, color=color, alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels([s.split("_")[0] for s in all_targets], fontsize=8, rotation=60, ha="right")
    ax.set_ylabel("Cluster Count (90% identity)")
    ax.set_title("Sequence Diversity: Cluster Count")
    ax.legend(fontsize=9)

    # Right: largest cluster fraction per model (boxplot)
    ax = axes[1]
    box_data = []
    box_labels = []
    for model in models:
        # One value per target (deduplicated)
        target_fracs = {}
        for r in all_data[model]:
            sid = r["sample_id"]
            f = safe_float(r.get("group_largest_cluster_frac", ""))
            if not np.isnan(f):
                target_fracs[sid] = f
        if target_fracs:
            box_data.append(list(target_fracs.values()))
            box_labels.append(model)
    if box_data:
        bp = ax.boxplot(box_data, tick_labels=box_labels, patch_artist=True, widths=0.5)
        for j, patch in enumerate(bp["boxes"]):
            color = MODEL_COLORS.get(box_labels[j], DEFAULT_COLORS[j % len(DEFAULT_COLORS)])
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
    ax.set_ylabel("Largest Cluster Fraction")
    ax.set_title("Sequence Redundancy (lower = more diverse)")
    ax.set_ylim(0, 1.05)

    fig.tight_layout()
    fig.savefig(out_dir / "fig11_diversity.png")
    fig.savefig(out_dir / "fig11_diversity.pdf")
    plt.close(fig)
    print(f"  [OK] fig11_diversity")


def fig_epitope_coverage(all_data: dict[str, list[dict]], out_dir: Path):
    """表位覆盖率: 设计的纳米抗体是否结合在预期表位上"""
    models = list(all_data.keys())
    has_epi = any(r.get("epitope_coverage", "") for rows in all_data.values() for r in rows)
    if not has_epi:
        print(f"  [SKIP] fig12_epitope_coverage (no epitope data)")
        return

    all_targets = sorted({r["sample_id"] for rows in all_data.values() for r in rows})

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: coverage per target
    ax = axes[0]
    x = np.arange(len(all_targets))
    width = 0.8 / max(1, len(models))
    for i, model in enumerate(models):
        target_cov = {}
        for r in all_data[model]:
            sid = r["sample_id"]
            c = safe_float(r.get("epitope_coverage", ""))
            if not np.isnan(c):
                target_cov.setdefault(sid, []).append(c)
        # Mean coverage per target
        vals = [np.mean(target_cov.get(sid, [0])) for sid in all_targets]
        color = MODEL_COLORS.get(model, DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
        ax.bar(x + i * width - 0.4 + width / 2, vals, width * 0.9, label=model, color=color, alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels([s.split("_")[0] for s in all_targets], fontsize=8, rotation=60, ha="right")
    ax.set_ylabel("Mean Epitope Coverage")
    ax.set_title("Epitope Coverage per Target")
    ax.axhline(y=0.7, color="#888", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.text(0.1, 0.72, "coverage=0.7", fontsize=8, color="#666")
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=9)

    # Right: distribution (violin or box)
    ax = axes[1]
    box_data = []
    box_labels = []
    for model in models:
        covs = [safe_float(r.get("epitope_coverage", "")) for r in all_data[model]]
        covs = [c for c in covs if not np.isnan(c)]
        if covs:
            box_data.append(covs)
            box_labels.append(model)
    if box_data:
        bp = ax.boxplot(box_data, tick_labels=box_labels, patch_artist=True, widths=0.5)
        for j, patch in enumerate(bp["boxes"]):
            color = MODEL_COLORS.get(box_labels[j], DEFAULT_COLORS[j % len(DEFAULT_COLORS)])
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
    ax.set_ylabel("Epitope Coverage")
    ax.set_title("Epitope Coverage Distribution")
    ax.axhline(y=0.7, color="#888", linestyle="--", linewidth=0.8)
    ax.set_ylim(0, 1.05)

    fig.tight_layout()
    fig.savefig(out_dir / "fig12_epitope_coverage.png")
    fig.savefig(out_dir / "fig12_epitope_coverage.pdf")
    plt.close(fig)
    print(f"  [OK] fig12_epitope_coverage")


def generate_summary_table(all_data: dict[str, list[dict]], out_dir: Path):
    """生成汇总统计表 (CSV + 打印)"""
    models = list(all_data.keys())
    summary_rows = []

    for model in models:
        rows = all_data[model]
        n = len(rows)
        iptms = [safe_float(r["external_iptm"]) for r in rows]
        iptms_clean = [v for v in iptms if not np.isnan(v)]
        pdockqs = [safe_float(r["pdockq2"]) for r in rows]
        pdockqs_clean = [v for v in pdockqs if not np.isnan(v)]
        bsas = [safe_float(r["interface_bsa"]) for r in rows]
        bsas_clean = [v for v in bsas if not np.isnan(v)]
        plddts = [safe_float(r["external_plddt_binder"]) for r in rows]
        plddts_clean = [v for v in plddts if not np.isnan(v)]

        cdr3_lens = [len(c) for r in rows if (c := extract_cdr3(r["pred_sequence"]))]

        summary_rows.append({
            "model": model,
            "n_targets": n,
            "median_iptm": f"{np.median(iptms_clean):.3f}" if iptms_clean else "",
            "mean_iptm": f"{np.mean(iptms_clean):.3f}" if iptms_clean else "",
            "max_iptm": f"{max(iptms_clean):.3f}" if iptms_clean else "",
            "iptm_ge_0.4": f"{sum(1 for v in iptms_clean if v >= 0.4)}/{n}",
            "iptm_ge_0.5": f"{sum(1 for v in iptms_clean if v >= 0.5)}/{n}",
            "median_pdockq2": f"{np.median(pdockqs_clean):.3f}" if pdockqs_clean else "",
            "mean_pdockq2": f"{np.mean(pdockqs_clean):.3f}" if pdockqs_clean else "",
            "pdockq2_ge_0.23": f"{sum(1 for v in pdockqs_clean if v > 0.23)}/{n}",
            "pdockq2_ge_0.4": f"{sum(1 for v in pdockqs_clean if v >= 0.4)}/{n}",
            "median_bsa": f"{np.median(bsas_clean):.0f}" if bsas_clean else "",
            "median_binder_plddt": f"{np.median(plddts_clean):.1f}" if plddts_clean else "",
            "median_cdr3_len": f"{np.median(cdr3_lens):.0f}" if cdr3_lens else "",
        })

    # 写 CSV
    csv_path = out_dir / "summary_table.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=summary_rows[0].keys())
        w.writeheader()
        w.writerows(summary_rows)

    # 打印
    print("\n" + "=" * 80)
    print("Model Comparison Summary")
    print("=" * 80)
    for sr in summary_rows:
        print(f"\n  Model: {sr['model']} ({sr['n_targets']} targets)")
        print(f"    iPTM:    median={sr['median_iptm']}  mean={sr['mean_iptm']}  max={sr['max_iptm']}")
        print(f"    iPTM:    ≥0.4: {sr['iptm_ge_0.4']}  ≥0.5: {sr['iptm_ge_0.5']}")
        print(f"    pDockQ2: median={sr['median_pdockq2']}  mean={sr['mean_pdockq2']}")
        print(f"    pDockQ2: >0.23: {sr['pdockq2_ge_0.23']}  ≥0.4: {sr['pdockq2_ge_0.4']}")
        print(f"    BSA:     median={sr['median_bsa']} Å²")
        print(f"    Binder pLDDT: median={sr['median_binder_plddt']}")
        print(f"    CDR3 length:  median={sr['median_cdr3_len']}")
    print()

    return summary_rows


# ── 主入口 ───────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="AntibodyBench thesis-level result analysis")
    parser.add_argument("--eval-dirs", nargs="+", required=True,
                        help="One or more evaluation output directories (each containing eval_merged.csv)")
    parser.add_argument("--model-labels", nargs="+", default=None,
                        help="Labels for each eval dir (default: auto-detect from CSV)")
    parser.add_argument("--out-dir", type=str, default="outputs/analysis/thesis_figures",
                        help="Output directory for figures and summary")
    parser.add_argument("--top-k", type=int, default=0,
                        help="Per (model, target) keep only top-K candidates by --rank-by metric. "
                             "0 = no filtering (use all candidates). Recommended: 10")
    parser.add_argument("--rank-by", type=str, default="chai1_iptm",
                        choices=["chai1_iptm", "pdockq2", "composite"],
                        help="Metric for top-K ranking (default: chai1_iptm)")
    args = parser.parse_args()

    # 加载数据
    all_data: dict[str, list[dict]] = {}
    for i, eval_dir_str in enumerate(args.eval_dirs):
        eval_dir = Path(eval_dir_str)
        rows = load_eval_merged(eval_dir)
        if not rows:
            continue
        if args.model_labels and i < len(args.model_labels):
            label = args.model_labels[i]
        else:
            # 自动检测模型名
            models_in_csv = set(r["model"] for r in rows)
            label = list(models_in_csv)[0] if len(models_in_csv) == 1 else eval_dir.name
        # 如果同一个 label 已有数据，合并
        if label in all_data:
            all_data[label].extend(rows)
        else:
            all_data[label] = rows

    if not all_data:
        print("[ERROR] No data loaded. Check --eval-dirs paths.", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {sum(len(v) for v in all_data.values())} candidates from {len(all_data)} model(s)")
    for model, rows in all_data.items():
        print(f"  {model}: {len(rows)} candidates across {len(set(r['sample_id'] for r in rows))} targets")

    # ── Top-K 归一化选择（公平比较核心步骤）──
    if args.top_k > 0:
        all_data = select_topk(all_data, top_k=args.top_k, rank_by=args.rank_by)
        suffix = f"_top{args.top_k}"
    else:
        suffix = "_all"

    # 输出目录加上 top-K 标识
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nGenerating figures in {out_dir}/\n")

    # 生成所有图表
    fig_iptm_by_target(all_data, out_dir)
    fig_pdockq2_by_target(all_data, out_dir)
    fig_metric_distributions(all_data, out_dir)
    fig_iptm_vs_pdockq2(all_data, out_dir)
    fig_success_rate_curve(all_data, out_dir)
    fig_interface_heatmap(all_data, out_dir)
    fig_sequence_analysis(all_data, out_dir)
    fig_model_comparison_radar(all_data, out_dir)
    fig_head_to_head(all_data, out_dir)
    fig_pareto_hit_rate(all_data, out_dir)
    fig_diversity(all_data, out_dir)
    fig_epitope_coverage(all_data, out_dir)

    # 汇总表
    generate_summary_table(all_data, out_dir)

    print(f"All figures saved to {out_dir}/")
    print("  PNG (300 dpi) for slides, PDF for LaTeX")


if __name__ == "__main__":
    main()
