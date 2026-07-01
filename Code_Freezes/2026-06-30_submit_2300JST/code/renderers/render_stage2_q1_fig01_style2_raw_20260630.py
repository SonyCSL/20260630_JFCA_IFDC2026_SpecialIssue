#!/usr/bin/env python3
"""Downstream-only Stage2 q1 Fig01 style2 redraw for raw outputs.

This renderer reads existing downstream combo payload ``outer.csv`` files and
does not recompute Stage 1/Stage 2 estimates.
"""

from __future__ import annotations

import argparse
import gzip
import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


DOMAINS = ["Formula", "Brite", "Pathway"]
SUBSETS = ["All", "PM", "SM"]
MODES = [
    "Combo6yr",
    "Period1_2015_2017",
    "Period2_2018_2020",
    "2015",
    "2016",
    "2017",
    "2018",
    "2019",
    "2020",
]
FIXED_ABS_R_ORDER = list(range(100, -1, -1))

TAIL_ORDER = ["High-tail cutoff r<0", "Low-tail cutoff r>0", "No cutoff r=0"]
TAIL_RENAME = {
    "Low tail r<0": "High-tail cutoff r<0",
    "High tail r>0": "Low-tail cutoff r>0",
    "No cutoff": "No cutoff r=0",
    "No cutoff r=0": "No cutoff r=0",
}
BOX_PALETTE = {
    "High-tail cutoff r<0": "#FDBA74",
    "Low-tail cutoff r>0": "#93C5FD",
    "No cutoff r=0": "#D1D5DB",
}
MEAN_MARKERS = {
    "High-tail cutoff r<0": "#F97316",
    "Low-tail cutoff r>0": "#2563EB",
    "No cutoff r=0": "#555555",
}
TAIL_OFFSETS = {
    "High-tail cutoff r<0": -0.27,
    "Low-tail cutoff r>0": 0.0,
    "No cutoff r=0": 0.27,
}
CLIFF_GUIDES = [
    (0.147, ":", "Small-or-larger: |delta|>=0.147"),
    (0.330, "--", "Medium-or-larger: |delta|>=0.330"),
    (0.474, "-.", "Large-or-larger: |delta|>=0.474"),
]
GUIDE_COLOR = "#D62728"
EFFECT_BANDS = [
    (0.330, 0.474, "#BFF6FA", 0.22),
    (0.474, 1.08, "#55DDE7", 0.18),
    (-0.474, -0.330, "#F9C7EF", 0.22),
    (-1.08, -0.474, "#F472D0", 0.18),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--domains", default=",".join(DOMAINS))
    parser.add_argument("--subsets", default=",".join(SUBSETS))
    parser.add_argument("--modes", default=",".join(MODES))
    parser.add_argument("--q-label", default="q1")
    parser.add_argument("--estimate", default="Stage2_iNEXT_TD_m_est")
    parser.add_argument("--metric", default="CliffsDelta")
    parser.add_argument("--sweep-sc", type=float, default=0.96)
    parser.add_argument("--dpi", type=int, default=240)
    parser.add_argument("--skip-main", action="store_true")
    return parser.parse_args()


def split_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def display_tick_label(abs_r: int) -> str:
    return str(abs_r) if abs_r % 5 == 0 else ""


def add_effect_background(ax) -> None:
    for ymin, ymax, color, alpha in EFFECT_BANDS:
        ax.axhspan(ymin, ymax, facecolor=color, alpha=alpha, edgecolor="none", zorder=0)


def add_cliff_guides(ax) -> None:
    for value, linestyle, _label in CLIFF_GUIDES:
        for sign in (-1, 1):
            ax.axhline(sign * value, color=GUIDE_COLOR, linestyle=linestyle, linewidth=1.05, alpha=0.88, zorder=1)


def add_boxplot_layer(ax, plot_df: pd.DataFrame, order_int: list[int]) -> None:
    for tail in TAIL_ORDER:
        data: list[np.ndarray] = []
        positions: list[float] = []
        for idx, abs_r in enumerate(order_int):
            values = plot_df[
                (plot_df["Tail_Label"].eq(tail)) & (plot_df["Abs_R_Pct"].eq(abs_r))
            ]["Value"].dropna().to_numpy(dtype=float)
            if len(values) == 0:
                continue
            data.append(values)
            positions.append(idx + TAIL_OFFSETS[tail])
        if not data:
            continue
        bp = ax.boxplot(
            data,
            positions=positions,
            widths=0.22,
            patch_artist=True,
            showfliers=False,
            manage_ticks=False,
            medianprops={"color": "#333333", "linewidth": 1.0},
            whiskerprops={"color": "#555555", "linewidth": 0.82},
            capprops={"color": "#555555", "linewidth": 0.82},
            boxprops={"edgecolor": "#555555", "linewidth": 0.95},
        )
        for box in bp["boxes"]:
            box.set_facecolor(BOX_PALETTE[tail])
            box.set_alpha(1.0)


def add_mean_markers(ax, plot_df: pd.DataFrame, order_int: list[int]) -> None:
    grouped = (
        plot_df.groupby(["Abs_R_Pct", "Tail_Label"], observed=True, as_index=False)["Value"]
        .mean()
        .dropna()
    )
    x_pos = {abs_r: idx for idx, abs_r in enumerate(order_int)}
    for _, row in grouped.iterrows():
        tail = str(row["Tail_Label"])
        abs_r = int(row["Abs_R_Pct"])
        if abs_r not in x_pos:
            continue
        ax.plot(
            x_pos[abs_r] + TAIL_OFFSETS.get(tail, 0.0),
            float(row["Value"]),
            marker="x",
            linestyle="None",
            color=MEAN_MARKERS.get(tail, "#333333"),
            markersize=5.0,
            markeredgewidth=1.5,
            zorder=5,
        )


def tail_legend_handles() -> list:
    return [
        Patch(facecolor=BOX_PALETTE["High-tail cutoff r<0"], edgecolor="#555555", label="High-tail cutoff r<0"),
        Patch(facecolor=BOX_PALETTE["Low-tail cutoff r>0"], edgecolor="#555555", label="Low-tail cutoff r>0"),
        Patch(facecolor=BOX_PALETTE["No cutoff r=0"], edgecolor="#555555", label="No cutoff r=0"),
        Line2D([0], [0], marker="x", linestyle="None", color=MEAN_MARKERS["High-tail cutoff r<0"], label="Mean r<0"),
        Line2D([0], [0], marker="x", linestyle="None", color=MEAN_MARKERS["Low-tail cutoff r>0"], label="Mean r>0"),
        Patch(facecolor="#55DDE7", alpha=0.18, edgecolor="none", label="Positive band"),
        Patch(facecolor="#F472D0", alpha=0.18, edgecolor="none", label="Negative band"),
    ]


def threshold_handles() -> list:
    return [
        Line2D([0], [0], color=GUIDE_COLOR, linestyle=linestyle, linewidth=1.3, label=label)
        for _value, linestyle, label in CLIFF_GUIDES
    ]


def add_outside_legends(fig, ax) -> None:
    first = ax.legend(
        handles=tail_legend_handles(),
        loc="upper left",
        bbox_to_anchor=(1.012, 1.0),
        borderaxespad=0.0,
        frameon=True,
        fontsize=8.0,
    )
    ax.add_artist(first)
    ax.legend(
        handles=threshold_handles(),
        title="Effect-size thresholds",
        loc="lower left",
        bbox_to_anchor=(1.012, 0.0),
        borderaxespad=0.0,
        frameon=True,
        fontsize=8.0,
        title_fontsize=8.5,
    )
    fig.subplots_adjust(right=0.765)


def render_legend_only(out_path: Path, *, dpi: int) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12.0, 0.82))
    ax.axis("off")
    handles = tail_legend_handles() + threshold_handles()
    ax.legend(
        handles=handles,
        loc="center",
        ncol=len(handles),
        frameon=True,
        fontsize=8.5,
        columnspacing=1.25,
        handlelength=2.0,
        borderpad=0.35,
    )
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", transparent=True)
    plt.close(fig)


def render_one(plot_df: pd.DataFrame, out_path: Path, *, dpi: int, include_legend: bool) -> None:
    if plot_df.empty:
        raise RuntimeError(f"Empty plot dataframe for {out_path}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    order_int = FIXED_ABS_R_ORDER
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.edgecolor": "#CFCFCF",
            "axes.linewidth": 0.8,
            "grid.color": "#D8D8D8",
            "grid.linewidth": 0.85,
            "grid.alpha": 1.0,
        }
    )
    fig_width = 20.4 if include_legend else 17.6
    fig, ax = plt.subplots(figsize=(fig_width, 6.6))
    ax.set_facecolor("#FFFFFF")
    add_effect_background(ax)
    add_boxplot_layer(ax, plot_df, order_int)
    add_cliff_guides(ax)
    add_mean_markers(ax, plot_df, order_int)
    ax.set_xticks(range(len(order_int)))
    ax.set_xticklabels([display_tick_label(v) for v in order_int], rotation=90, fontsize=8)
    ax.set_xlim(-0.5, len(order_int) - 0.5)
    ax.set_xlabel("Intensity cutoff |r|", fontsize=11)
    ax.set_ylabel("Effect size delta", fontsize=11)
    ax.set_ylim(-1.08, 1.08)
    ax.grid(axis="y", zorder=0)
    ax.grid(axis="x", color="#ECECEC", linewidth=0.4, alpha=0.65, zorder=0)
    if include_legend:
        add_outside_legends(fig, ax)
    else:
        fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def normalize_tail(row: pd.Series) -> str:
    side = str(row.get("Cutoff_Side", "")).lower()
    if side == "low":
        return "Low tail r<0"
    if side == "high":
        return "High tail r>0"
    return "No cutoff"


def read_filtered_outer(args: argparse.Namespace, domains: list[str], subsets: list[str], modes: list[str]) -> pd.DataFrame:
    payload_root = args.source_root / "csv" / "combo_payloads"
    paths = sorted(payload_root.glob("combo_*/outer.csv"))
    if not paths:
        raise FileNotFoundError(f"No outer.csv files found under {payload_root}")
    usecols = [
        "Analysis_Tier",
        "Domain",
        "Subset",
        "Mode",
        "Q_Label",
        "Estimate_Definition",
        "Retain_Pct",
        "Cutoff_Side",
        "Reference_Target_SC",
        "Sweep_SC",
        "Sweep_SC_Label",
        "OuterRep",
        "Metric",
        "Value",
    ]
    parts: list[pd.DataFrame] = []
    for path in paths:
        df = pd.read_csv(path, usecols=lambda c: c in usecols, low_memory=False)
        df["Mode"] = df["Mode"].astype(str)
        df["Domain"] = df["Domain"].astype(str)
        df["Subset"] = df["Subset"].astype(str)
        df["Q_Label"] = df["Q_Label"].astype(str)
        df["Estimate_Definition"] = df["Estimate_Definition"].astype(str)
        df["Metric"] = df["Metric"].astype(str)
        df = df[
            df["Analysis_Tier"].eq("main_q12")
            & df["Domain"].isin(domains)
            & df["Subset"].isin(subsets)
            & df["Mode"].isin(modes)
            & df["Q_Label"].eq(args.q_label)
            & df["Estimate_Definition"].eq(args.estimate)
            & df["Metric"].eq(args.metric)
        ].copy()
        if df.empty:
            continue
        df["Sweep_SC"] = pd.to_numeric(df["Sweep_SC"], errors="coerce")
        df = df[np.isclose(df["Sweep_SC"], args.sweep_sc, atol=1e-9)].copy()
        if not df.empty:
            parts.append(df)
    if not parts:
        raise RuntimeError("No rows matched the requested q1 Fig01 filters.")
    out = pd.concat(parts, ignore_index=True)
    out["Retain_Pct"] = pd.to_numeric(out["Retain_Pct"], errors="coerce")
    out["Abs_R_Pct"] = (101 - out["Retain_Pct"]).round().astype("Int64")
    out.loc[out["Cutoff_Side"].astype(str).str.lower().eq("none"), "Abs_R_Pct"] = 0
    out["Abs_R_Pct"] = out["Abs_R_Pct"].clip(lower=0, upper=100).astype(int)
    out["Abs_R_F"] = out["Abs_R_Pct"].astype(str)
    out["Tail_Label"] = out.apply(normalize_tail, axis=1).map(TAIL_RENAME)
    out["Value"] = pd.to_numeric(out["Value"], errors="coerce")
    out["OuterRep"] = pd.to_numeric(out["OuterRep"], errors="coerce")
    out = out[np.isfinite(out["Value"])].copy()
    return out[
        [
            "OuterRep",
            "Domain",
            "Subset",
            "Mode",
            "Retain_Pct",
            "Abs_R_Pct",
            "Abs_R_F",
            "Tail_Label",
            "Reference_Target_SC",
            "Sweep_SC",
            "Sweep_SC_Label",
            "Value",
        ]
    ]


def validate_expected(df: pd.DataFrame, domains: list[str], subsets: list[str], modes: list[str]) -> pd.DataFrame:
    rows: list[dict] = []
    for domain in domains:
        for subset in subsets:
            for mode in modes:
                one = df[(df["Domain"].eq(domain)) & (df["Subset"].eq(subset)) & (df["Mode"].eq(mode))]
                rows.append(
                    {
                        "Transform": "raw",
                        "Domain": domain,
                        "Subset": subset,
                        "Mode": mode,
                        "Rows": int(len(one)),
                        "OuterRep_N": int(one["OuterRep"].nunique()) if not one.empty else 0,
                        "Abs_R_N": int(one["Abs_R_Pct"].nunique()) if not one.empty else 0,
                        "Tail_N": int(one["Tail_Label"].nunique()) if not one.empty else 0,
                    }
                )
    summary = pd.DataFrame(rows)
    missing = summary[summary["Rows"].eq(0)]
    if not missing.empty:
        raise RuntimeError(f"Missing required plot rows:\n{missing.to_string(index=False)}")
    return summary


def render_all(df: pd.DataFrame, args: argparse.Namespace, domains: list[str], subsets: list[str], modes: list[str]) -> list[dict]:
    manifest_rows: list[dict] = []
    for domain in domains:
        for subset in subsets:
            for mode in modes:
                one = df[(df["Domain"].eq(domain)) & (df["Subset"].eq(subset)) & (df["Mode"].eq(mode))].copy()
                base_dir = (
                    args.out_root
                    / "raw"
                    / "figures"
                    / "main_q12"
                    / domain
                    / subset
                    / mode
                    / args.q_label
                    / args.estimate
                    / args.metric
                )
                main_path = base_dir / "01_reference_bridge_boxplots_paired.png"
                nolegend_path = base_dir / "01_reference_bridge_boxplots_paired_nolegend.png"
                print(f"[PLOT] raw {domain} {subset} {mode} rows={len(one)}", flush=True)
                if not args.skip_main:
                    render_one(one, main_path, dpi=args.dpi, include_legend=True)
                render_one(one, nolegend_path, dpi=args.dpi, include_legend=False)
                manifest_rows.append(
                    {
                        "Transform": "raw",
                        "Transform_Title": "raw",
                        "Domain": domain,
                        "Subset": subset,
                        "Mode": mode,
                        "Q_Label": args.q_label,
                        "Estimate_Definition": args.estimate,
                        "Metric": args.metric,
                        "Rows": int(len(one)),
                        "Figure_File": str(main_path),
                        "NoLegend_Figure_File": str(nolegend_path),
                    }
                )
    return manifest_rows


def main() -> int:
    args = parse_args()
    started = time.time()
    args.out_root.mkdir(parents=True, exist_ok=True)
    domains = split_csv(args.domains)
    subsets = split_csv(args.subsets)
    modes = split_csv(args.modes)

    df = read_filtered_outer(args, domains, subsets, modes)
    filtered_path = args.out_root / f"filtered_outer_raw_{args.q_label}_stage2_inext_sweep096.tsv.gz"
    with gzip.open(filtered_path, "wt", encoding="utf-8") as handle:
        df.to_csv(handle, sep="\t", index=False)

    summary = validate_expected(df, domains, subsets, modes)
    summary_path = args.out_root / "plot_row_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    manifest = pd.DataFrame(render_all(df, args, domains, subsets, modes))
    manifest_path = args.out_root / "figure01_boxplot_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)
    legend_path = args.out_root / "figure01_boxplot_legend_style2_horizontal_single_row.png"
    render_legend_only(legend_path, dpi=args.dpi)

    note = {
        "elapsed_sec": round(time.time() - started, 3),
        "source_root": str(args.source_root),
        "out_root": str(args.out_root),
        "domains": domains,
        "subsets": subsets,
        "modes": modes,
        "q_label": args.q_label,
        "estimate": args.estimate,
        "metric": args.metric,
        "sweep_sc": args.sweep_sc,
        "filtered_tsv": str(filtered_path),
        "summary_tsv": str(summary_path),
        "manifest_tsv": str(manifest_path),
        "legend_png": str(legend_path),
        "main_figure_count": int(len(manifest)),
        "nolegend_figure_count": int(len(manifest)),
        "upstream_recomputed": False,
    }
    (args.out_root / "run_note.json").write_text(json.dumps(note, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    print(f"[DONE] wrote {len(manifest)} main figures under {args.out_root}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
