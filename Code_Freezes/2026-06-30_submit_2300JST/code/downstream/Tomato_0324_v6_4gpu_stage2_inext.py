#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# Python port of Tomato_0324_v2.R adapted for v36.9 Stage2 iNEXT-bootstrap outputs
#
# 日本語メモ:
#   2026-06-21 の v6_4gpu_stage2_inext は v5_4gpu からの薄い互換分岐。
#   新 upstream `v36_9_tomato_long_q012_patch_4gpu_v3_stage2_inext.py` の
#   `Full_v36_9_stage2_inext_4gpu_*` / `PoC_v36_9_stage2_inext_4gpu_*`
#   を auto-detect し、BRITE / Disease などの domain alias を展開する。
#
#   2026-04-08 の v5_4gpu では、A 図に差なし基準 0.5 の破線を追加した。
#   あわせて、upstream v2 が出力する Welch の -log10(p) も metric として読めるようにした。
#   BM を主に読む方針は維持し、Welch は補助線として扱う。
#
# Purpose
#   Read the v36.9 Stage2 iNEXT-bootstrap long CSVs and generate signed-cutoff plots where:
#     Stage 1) reference pseudo-samples are fixed at Reference_Target_SC = 0.965
#     Stage 2) evaluation sweeps across Sweep_SC = 1%..100% (100% = reference endpoint)
#
# Inputs
#   1) comparison_cells_long_all_q012_v36_9.csv
#   2) comparison_outer_v36_9_q012.csv
#   3) group_outer_v36_9_q012.csv
#
# Main recommendations implemented here
#   - keep the stage-1 reference target SC distinct from the stage-2 Sweep_SC axis
#   - show paired high-tail vs low-tail cutoff views
#   - make Sweep_SC the main x-axis for rarefaction / extrapolation interpretation
#   - keep effort-ratio support as supplemental only
#   - NegLogP / WelchNegLogP get p=0.05 / p=0.01 guides
#   - A gets a no-difference guide at 0.5
# ------------------------------------------------------------------------------

from __future__ import annotations

import math
import os
import re
import subprocess
import sys
import time
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", message="All-NaN slice encountered", category=RuntimeWarning)
try:
    import seaborn as sns
except ModuleNotFoundError:
    sns = None
try:
    import torch
except ModuleNotFoundError:
    torch = None


DEFAULT_ROOT_DIR = Path("/home/nonaka/work/nonaka/Chao1_Intensity")
LOCAL_FALLBACK_ROOT = Path(__file__).resolve().parent
ROOT_DIR = Path(os.getenv("TOMATO_ROOT_DIR", str(DEFAULT_ROOT_DIR)))
if not ROOT_DIR.exists():
    ROOT_DIR = LOCAL_FALLBACK_ROOT


def parse_int_env(name: str, default: int) -> int:
    raw = os.getenv(name, "")
    try:
        return int(raw)
    except (TypeError, ValueError):
        return default


def parse_float_env(name: str, default: float) -> float:
    raw = os.getenv(name, "")
    try:
        val = float(raw)
    except (TypeError, ValueError):
        return default
    return val if math.isfinite(val) else default


def parse_env_list(name: str, fallback: Optional[List[str]] = None) -> Optional[List[str]]:
    raw = os.getenv(name, "").strip()
    if not raw:
        return fallback
    vals = [item.strip() for item in raw.split(",")]
    vals = [item for item in vals if item]
    return vals if vals else fallback


DISEASE_DOMAIN_SLUGS = ["Disease_NE"]
DOMAIN_FILTER_ALIASES = {
    "Formula": ["Formula"],
    "formula": ["Formula"],
    "BRITE": ["Brite"],
    "Brite": ["Brite"],
    "brite": ["Brite"],
    "Network": ["Network"],
    "network": ["Network"],
    "Disease": DISEASE_DOMAIN_SLUGS,
    "disease": DISEASE_DOMAIN_SLUGS,
    "Disease_All": ["Disease_NE", "Disease_ICD11", "Disease_PathCL"],
    "Disease_NE": ["Disease_NE"],
    "Disease_ICD11": ["Disease_ICD11"],
    "Disease_PathCL": ["Disease_PathCL"],
}

MODE_FILTER_ALIASES = {
    "Annual": [str(year) for year in range(2015, 2021)],
    "annual": [str(year) for year in range(2015, 2021)],
    "SingleYear": [str(year) for year in range(2015, 2021)],
    "singleyear": [str(year) for year in range(2015, 2021)],
    "First3": ["Period1_2015_2017"],
    "first3": ["Period1_2015_2017"],
    "Last3": ["Period2_2018_2020"],
    "last3": ["Period2_2018_2020"],
    "Period1": ["Period1_2015_2017"],
    "Period2": ["Period2_2018_2020"],
    "Combo": ["Combo6yr"],
    "combo": ["Combo6yr"],
}


def expand_filter_aliases(values: Optional[List[str]], aliases: Dict[str, List[str]]) -> Optional[List[str]]:
    if values is None:
        return None
    out: List[str] = []
    for value in values:
        expanded = aliases.get(value, [value])
        for item in expanded:
            if item not in out:
                out.append(item)
    return out


POC_MODE = parse_int_env("TOMATO_POC_MODE", 0)
PRIMARY_REFERENCE_TARGET_SC = parse_float_env("TOMATO_REFERENCE_TARGET_SC", 0.965)
REFERENCE_SWEEP_OVERRIDE = parse_float_env("TOMATO_REFERENCE_SWEEP_SC", float("nan"))
PRIMARY_SC_TOL = 1e-9
REPRESENTATIVE_RETAIN_PCTS = [1, 5, 10, 20, 40, 60, 80, 95, 100]
REPRESENTATIVE_SWEEP_TARGETS = [0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.965, 1.0]
TAIL_PALETTE = {"High-tail": "#62C7D0", "Low-tail": "#E36A9B", "Baseline": "#7A7A7A"}
FIG01_ABS_R_ORDER = list(range(100, -1, -1))
FIG01_TAIL_ORDER = ["High-tail cutoff r<0", "Low-tail cutoff r>0", "No cutoff r=0"]
FIG01_TAIL_PALETTE = {
    "High-tail cutoff r<0": "#FDBA74",
    "Low-tail cutoff r>0": "#93C5FD",
    "No cutoff r=0": "#D1D5DB",
}
FIG01_MEAN_MARKERS = {
    "High-tail cutoff r<0": "#F97316",
    "Low-tail cutoff r>0": "#2563EB",
    "No cutoff r=0": "#555555",
}
FIG01_TAIL_OFFSETS = {
    "High-tail cutoff r<0": -0.27,
    "Low-tail cutoff r>0": 0.0,
    "No cutoff r=0": 0.27,
}
CLIFF_SMALL_THRESHOLD = 0.147
CLIFF_HEATMAP_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "cliff_cyan_white_magenta",
    ["#C000B8", "#FFFFFF", "#00BFD8"],
)
CLIFF_HEATMAP_CMAP.set_bad("#E6E6E6")
CLIFF_MEDIUM_THRESHOLD = 0.33
CLIFF_LARGE_THRESHOLD = 0.474
CLIFF_MASK_CMAP = mcolors.ListedColormap(["#D73027", "#FFFFFF", "#2C7BB6"])
CLIFF_MASK_CMAP.set_bad("#D9D9D9")
CLIFF_MASK_NORM = mcolors.BoundaryNorm([-1.5, -0.5, 0.5, 1.5], CLIFF_MASK_CMAP.N)
SIGNED_CUTOFF_DISPLAY_COL = "Signed_Cutoff_Display_Pct"
SIGNED_CUTOFF_AXIS_LABEL = "Intensity cutoff $r$ [%]"
SIGNED_CUTOFF_DISPLAY_VALUES = list(range(-99, 100))
SIGNED_CUTOFF_TICK_TARGETS = [-95, -75, -50, -25, 0, 25, 50, 75, 95]
POSITIVE_SURFACE_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "cliff_positive_surface",
    ["#FFFFFF", "#B2E2E2", "#0571B0"],
)
POSITIVE_SURFACE_CMAP.set_bad("#D9D9D9")
NEGATIVE_SURFACE_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "cliff_negative_surface",
    ["#FFFFFF", "#F1B6DA", "#C51B7D"],
)
NEGATIVE_SURFACE_CMAP.set_bad("#D9D9D9")
POSITIVE_MASK_CMAP = mcolors.ListedColormap(["#FFFFFF", "#B2E2E2", "#0571B0"])
POSITIVE_MASK_CMAP.set_bad("#D9D9D9")
NEGATIVE_MASK_CMAP = mcolors.ListedColormap(["#FFFFFF", "#F1B6DA", "#C51B7D"])
NEGATIVE_MASK_CMAP.set_bad("#D9D9D9")
SIGNED_MASK_NORM = mcolors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5], POSITIVE_MASK_CMAP.N)
EMBED_FIGURE_CAPTIONS = os.getenv("TOMATO_EMBED_FIGURE_CAPTIONS", "0").strip() == "1"
SWEEP_SHAPES = {
    "both_rarefaction": "o",
    "mixed": "^",
    "both_extrapolation": "s",
    "both_asymptotic": "D",
    "unknown": "x",
}

MAX_COMBOS_TO_EXPORT = parse_int_env("TOMATO_MAX_COMBOS", 0)
if MAX_COMBOS_TO_EXPORT <= 0:
    MAX_COMBOS_TO_EXPORT = math.inf
WORKER_MODE = parse_int_env("TOMATO_4GPU_WORKER", 0)
WORKER_GPU_ID = parse_int_env("TOMATO_4GPU_GPU_ID", 0)
REQUESTED_NUM_WORKERS = max(1, parse_int_env("TOMATO_4GPU_NUM_WORKERS", 4))
COMBO_START = parse_int_env("TOMATO_4GPU_COMBO_START", 0)
COMBO_END = parse_int_env("TOMATO_4GPU_COMBO_END", 0)
SKIP_SUMMARY_WRITE = parse_int_env("TOMATO_4GPU_SKIP_SUMMARY_WRITE", 0) == 1
SKIP_FIGURE_EXPORT = parse_int_env("TOMATO_4GPU_SKIP_FIGURE_EXPORT", 0) == 1
OUT_DIR_OVERRIDE = os.getenv("TOMATO_OUT_DIR", "").strip()
GPU_ID_RAW = parse_env_list("TOMATO_4GPU_GPU_IDS")
FIGURE_FILTER_RAW = parse_env_list("TOMATO_FIGURE_FILTER")
EXPORT_FORMATS_RAW = parse_env_list("TOMATO_EXPORT_FORMATS")
WRITE_COMBO_PAYLOADS = parse_int_env("TOMATO_WRITE_COMBO_PAYLOADS", 1) == 1


def normalize_figure_filter(raw: Optional[List[str]]) -> Optional[set[str]]:
    if raw is None:
        return None
    out: set[str] = set()
    for item in raw:
        token = str(item).strip().lower()
        if token in {"06_2", "6_2", "062", "06-2", "6-2"}:
            out.add("06_2")
            continue
        match = re.search(r"(\d+)", token)
        if not match:
            continue
        out.add(f"{int(match.group(1)):02d}")
    return out or None


FIGURE_FILTER = normalize_figure_filter(FIGURE_FILTER_RAW)


def normalize_export_formats(raw: Optional[List[str]]) -> set[str]:
    if raw is None:
        formats = {"pdf", "png"}
    else:
        formats = {str(item).strip().lower().lstrip(".") for item in raw if str(item).strip()}
        formats = {fmt for fmt in formats if fmt in {"pdf", "png"}}
        if not formats:
            formats = {"pdf", "png"}
    export_pdf_raw = os.getenv("TOMATO_EXPORT_PDF", "").strip().lower()
    if export_pdf_raw in {"0", "false", "no", "off"}:
        formats.discard("pdf")
    export_png_raw = os.getenv("TOMATO_EXPORT_PNG", "").strip().lower()
    if export_png_raw in {"0", "false", "no", "off"}:
        formats.discard("png")
    if not formats:
        raise RuntimeError("No figure export formats remain after TOMATO_EXPORT_FORMATS / TOMATO_EXPORT_* filters.")
    return formats


EXPORT_FORMATS = normalize_export_formats(EXPORT_FORMATS_RAW)


def should_export_figure(figure_id: str) -> bool:
    return FIGURE_FILTER is None or figure_id in FIGURE_FILTER


ONLY_FIGURE_06 = FIGURE_FILTER == {"06"}
ONLY_FIGURE_06_2 = FIGURE_FILTER == {"06_2"}
ONLY_FIGURE_01 = FIGURE_FILTER == {"01"}
ONLY_FIGURE_01_06 = FIGURE_FILTER == {"01", "06"}
LIGHTWEIGHT_FIGURE_FILTER = (
    FIGURE_FILTER is not None
    and len(FIGURE_FILTER) > 0
    and FIGURE_FILTER.issubset({"01", "06", "06_2"})
)


def detect_latest_run(parent_dir: Path, prefix: str) -> Path:
    runs = [
        path
        for path in parent_dir.iterdir()
        if path.is_dir() and re.match(rf"^{re.escape(prefix)}_\d{{8}}_\d{{6}}$", path.name)
    ]
    if not runs:
        raise RuntimeError(f"No run directory found for prefix: {prefix}")
    runs.sort(key=lambda path: path.stat().st_mtime, reverse=True)
    return runs[0]


RUN_DIR_SOURCE = "env:TOMATO_RUN_DIR"
run_dir_env = os.getenv("TOMATO_RUN_DIR", "").strip()
if run_dir_env:
    RUN_DIR = Path(run_dir_env)
else:
    RUN_DIR = detect_latest_run(
        ROOT_DIR / "out_sensitivity",
        "PoC_v36_9_stage2_inext_4gpu" if POC_MODE == 1 else "Full_v36_9_stage2_inext_4gpu",
    )
    RUN_DIR_SOURCE = "auto_detected_latest"

CELLS_CSV = RUN_DIR / "comparison_cells_long_all_q012_v36_9.csv"
OUTER_CSV = RUN_DIR / "comparison_outer_v36_9_q012.csv"
GROUP_CSV = RUN_DIR / "group_outer_v36_9_q012.csv"

TS_TAG = time.strftime("%Y%m%d_%H%M%S")
OUT_DIR = (
    Path(OUT_DIR_OVERRIDE)
    if OUT_DIR_OVERRIDE
    else ROOT_DIR / "out_heatmap_cutoff_quantification_from_csv" / (
        f"{'v36_9_stage2_inext_poc_recommended_4gpu_' if POC_MODE == 1 else 'v36_9_stage2_inext_recommended_4gpu_'}{TS_TAG}"
    )
)
DIR_CSV = OUT_DIR / "csv"
DIR_FIG = OUT_DIR / "figures"
for path in (OUT_DIR, DIR_CSV, DIR_FIG):
    path.mkdir(parents=True, exist_ok=True)
DIR_COMBO_PAYLOAD = DIR_CSV / "combo_payloads"

DEVICE = (
    f"cuda:{WORKER_GPU_ID}"
    if torch is not None and torch.cuda.is_available() and os.getenv("TOMATO_FORCE_CPU", "0") != "1"
    else "cpu"
)
DOMAIN_FILTER = None
SUBSET_FILTER = None
MODE_FILTER = None
Q_FILTER = None
ESTIMATE_FILTER = None
METRIC_FILTER = None


def log(msg: str) -> None:
    print(msg, flush=True)


def clean_id(x: object) -> str:
    out = re.sub(r"[^A-Za-z0-9._-]+", "_", str(x))
    out = re.sub(r"_+", "_", out)
    out = re.sub(r"^_|_$", "", out)
    return out


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def combo_payload_dir(combo_id: int) -> Path:
    return DIR_COMBO_PAYLOAD / f"combo_{int(combo_id):04d}"


def must_read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise RuntimeError(f"Required CSV is missing: {path}")
    return pd.read_csv(path)


def maybe_filter_in(df: pd.DataFrame, col: str, vals: Optional[Sequence[str]]) -> pd.DataFrame:
    if vals is None:
        return df
    return df[df[col].isin(vals)].copy()


def metric_pretty(x: str) -> str:
    return {
        "NegLogP": "-log10(BM p)",
        "WelchNegLogP": "-log10(Welch p)",
        "A": "A = P(Syneco > Conv)",
        "CliffsDelta": "Cliff's delta",
    }.get(x, str(x))


def metric_value_col(x: str) -> Optional[str]:
    return {
        "NegLogP": "NegLogP",
        "WelchNegLogP": "WelchNegLogP",
        "A": "A",
        "CliffsDelta": "CliffsDelta",
    }.get(x)


def estimate_pretty(x: str) -> str:
    return {
        "Observed_Standardized": "Observed standardized",
        "Chao1_asymptotic": "Chao1 asymptotic",
        "Empirical_Hill": "Empirical Hill",
        "Stage2_iNEXT_Bootstrap_Hill": "Stage2 iNEXT-bootstrap Hill",
        "Stage2_iNEXT_TD_m_est": "Stage2 iNEXT TD.m.est",
    }.get(x, str(x))


def analysis_tier(q_label: str) -> str:
    return "main_q12" if q_label in {"q1", "q2"} else "reference_q0"


def q0_caution_text(q_label: str, estimate_definition: str) -> Optional[str]:
    if q_label == "q0" and estimate_definition == "Chao1_asymptotic":
        return "q0 Chao1 asymptotic should be read as a cautious lower-bound-oriented reference."
    if q_label == "q0" and estimate_definition == "Observed_Standardized":
        return "q0 observed standardized is empirical, not an asymptotic endpoint."
    return None


def tail_display(cutoff_side: str) -> str:
    if cutoff_side == "high":
        return "High-tail"
    if cutoff_side == "low":
        return "Low-tail"
    return "Baseline"


def format_sweep_label(x: float) -> str:
    return "Asymptotic" if abs(x - 1.0) < PRIMARY_SC_TOL else f"{x:.3f}"


def format_sweep_pct(x: float) -> str:
    if abs(x - 1.0) < PRIMARY_SC_TOL:
        return "Asymptotic"
    return f"{int(round(x * 100))}%"


def format_sweep_pct_number(x: float) -> str:
    if abs(x - 1.0) < PRIMARY_SC_TOL:
        return "Asymptotic"
    return f"{int(round(x * 100))}"


def retain_breaks(values: Iterable[int]) -> List[str]:
    vals = sorted({int(v) for v in values})
    keep = [v for v in vals if v % 5 == 0 or v in {1, 100}]
    return [str(v) for v in keep]


def pick_representative_retain_pcts(values: Iterable[int]) -> List[int]:
    avail = sorted({int(v) for v in values})
    return [v for v in REPRESENTATIVE_RETAIN_PCTS if v in avail]


def pick_nearest_available(available: Iterable[float], targets: Iterable[float]) -> List[float]:
    avail = sorted({float(v) for v in available if math.isfinite(float(v))})
    if not avail:
        return []
    out = []
    for target in targets:
        out.append(min(avail, key=lambda v: abs(v - float(target))))
    dedup: List[float] = []
    for value in out:
        if all(abs(value - existing) > 1e-12 for existing in dedup):
            dedup.append(value)
    return dedup


def pick_reference_sweep(available: Iterable[float], reference_target: float, override: float = float("nan")) -> float:
    target = override if math.isfinite(override) else reference_target
    out = pick_nearest_available(available, [target])
    return out[0]


def summarise_distribution(values: Sequence[float]) -> Dict[str, float]:
    arr = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return {
            "n": 0,
            "mean_value": np.nan,
            "sd_value": np.nan,
            "median_value": np.nan,
            "q1": np.nan,
            "q3": np.nan,
            "iqr": np.nan,
            "p05": np.nan,
            "p95": np.nan,
            "min_value": np.nan,
            "max_value": np.nan,
        }
    return {
        "n": int(arr.size),
        "mean_value": float(np.mean(arr)),
        "sd_value": float(np.std(arr, ddof=1)) if arr.size >= 2 else np.nan,
        "median_value": float(np.median(arr)),
        "q1": float(np.quantile(arr, 0.25)),
        "q3": float(np.quantile(arr, 0.75)),
        "iqr": float(np.quantile(arr, 0.75) - np.quantile(arr, 0.25)),
        "p05": float(np.quantile(arr, 0.05)),
        "p95": float(np.quantile(arr, 0.95)),
        "min_value": float(np.min(arr)),
        "max_value": float(np.max(arr)),
    }


def compute_norms_vs_zero(values: Sequence[float]) -> Dict[str, float]:
    arr = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return {"frobenius": np.nan, "rmse": np.nan}
    return {"frobenius": float(np.sqrt(np.sum(arr ** 2))), "rmse": float(np.sqrt(np.mean(arr ** 2)))}


def empty_outer_metric_summary(group_cols: Sequence[str]) -> pd.DataFrame:
    summary_cols = list(summarise_distribution([]).keys()) + list(compute_norms_vs_zero([]).keys())
    extra_cols = [
        "Metric",
        "Metric_Pretty",
        "Sweep_Status",
        "Primary_Eligible",
        "Max_Sweep_M_over_N",
    ]
    return pd.DataFrame(columns=list(group_cols) + extra_cols + summary_cols)


def safe_max_or_na(values: Sequence[float]) -> float:
    arr = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(arr.max()) if arr.size else np.nan


def collapse_primary_eligible(values: Sequence[object]) -> object:
    vals = pd.Series(values).dropna().astype(bool)
    if vals.empty:
        return np.nan
    return bool(vals.all())


def collapse_sweep_status(status_vec: Sequence[object], sweep_is_asymptotic_vec: Sequence[object]) -> str:
    asym = pd.Series(sweep_is_asymptotic_vec).fillna(False).astype(bool)
    if bool(asym.any()):
        return "both_asymptotic"
    statuses = pd.Series(status_vec).dropna().astype(str).tolist()
    statuses = list(dict.fromkeys(statuses))
    if not statuses:
        return "unknown"
    if len(statuses) == 1:
        return statuses[0]
    return "mixed"


def empty_plot(title_txt: str, subtitle_txt: str = "No data available after filtering."):
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.text(0.5, 0.5, subtitle_txt, ha="center", va="center", fontsize=11)
    ax.set_axis_off()
    fig.suptitle(title_txt)
    return fig


def save_plot_pdf_png(fig, pdf_path: Path, png_path: Path, width: float, height: float) -> None:
    if fig is None:
        return
    fig.set_size_inches(width, height)
    if not getattr(fig, "_tomato_skip_tight_layout", False):
        fig.tight_layout(rect=(0, 0.03, 1, 0.98))
    if "pdf" in EXPORT_FORMATS:
        fig.savefig(pdf_path, bbox_inches="tight")
    if "png" in EXPORT_FORMATS:
        fig.savefig(png_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def make_combo_title(df_one: pd.DataFrame, suffix: str) -> str:
    row = df_one.iloc[0]
    return (
        f"{row['Domain']} | {row['Subset']} | {row['Mode']} | {row['Q_Label']} | "
        f"{row['Estimate_Pretty']} | {row['Metric_Pretty']} : {suffix}"
    )


def make_group_title(df_one: pd.DataFrame, suffix: str) -> str:
    row = df_one.iloc[0]
    return f"{row['Domain']} | {row['Subset']} | {row['Mode']} : {suffix}"


def add_neglogp_guides(ax, metric_key: str) -> None:
    if metric_key not in {"NegLogP", "WelchNegLogP"}:
        return
    ax.axhline(-math.log10(0.05), linestyle="--", color="#6E6E6E", linewidth=0.8)
    ax.axhline(-math.log10(0.01), linestyle=":", color="#444444", linewidth=0.8)


def add_a_guides(ax, metric_key: str, axis: str = "y") -> None:
    if metric_key != "A":
        return
    if axis == "x":
        ax.axvline(0.5, linestyle="--", color="#6E6E6E", linewidth=0.8)
    else:
        ax.axhline(0.5, linestyle="--", color="#6E6E6E", linewidth=0.8)


def summarise_outer_metric(df: pd.DataFrame, value_col: str, metric_key: str) -> pd.DataFrame:
    group_cols = [
        "Analysis_Tier",
        "Domain",
        "Subset",
        "Mode",
        "Diversity_Order",
        "Q_Label",
        "Estimate_Definition",
        "Estimate_Pretty",
        "Q0_Caution",
        "Cutoff_Key",
        "Cutoff_Side",
        "Tail_Display",
        "Retain_Pct",
        "Retain_Ratio",
        "Cutoff_Signed_Pct",
        "Cutoff_Label",
        "Reference_Target_SC",
        "Reference_Sweep_Selected",
        "Reference_Sweep_Label",
        "Is_Reference_Sweep_Selected",
        "Sweep_SC",
        "Sweep_SC_Label",
        "Sweep_SC_Pct",
        "Sweep_Is_Asymptotic",
    ]
    if df.empty:
        return empty_outer_metric_summary(group_cols)
    rows: List[Dict[str, object]] = []
    for keys, grp in df.groupby(group_cols, dropna=False, sort=False):
        row = dict(zip(group_cols, keys if isinstance(keys, tuple) else (keys,)))
        row["Metric"] = metric_key
        row["Metric_Pretty"] = metric_pretty(metric_key)
        row["Sweep_Status"] = collapse_sweep_status(grp["Sweep_Status"], grp["Sweep_Is_Asymptotic"])
        row["Primary_Eligible"] = collapse_primary_eligible(grp["Primary_Eligible"])
        row["Max_Sweep_M_over_N"] = safe_max_or_na(grp["Max_Sweep_M_over_N"])
        row.update(summarise_distribution(grp[value_col]))
        row.update(compute_norms_vs_zero(grp[value_col]))
        rows.append(row)
    return pd.DataFrame(rows)


def filter_by_combo_keys(df: pd.DataFrame, one: pd.Series, key_cols: Sequence[str]) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    missing = [col for col in key_cols if col not in df.columns]
    if missing:
        return df.iloc[0:0].copy()
    mask = pd.Series(True, index=df.index)
    for col in key_cols:
        mask &= df[col] == one[col]
    return df[mask].copy()


def materialize_outer_metric(df: pd.DataFrame, metric_key: str) -> pd.DataFrame:
    value_col = metric_value_col(metric_key)
    if value_col is None:
        raise RuntimeError(f"Unsupported metric key: {metric_key}")
    out = df[
        [
            "Analysis_Tier",
            "Domain",
            "Subset",
            "Mode",
            "Diversity_Order",
            "Q_Label",
            "Estimate_Definition",
            "Estimate_Pretty",
            "Q0_Caution",
            "Cutoff_Key",
            "Cutoff_Side",
            "Tail_Display",
            "Retain_Pct",
            "Retain_Ratio",
            "Cutoff_Signed_Pct",
            "Cutoff_Label",
            "Reference_Target_SC",
            "Reference_Sweep_Selected",
            "Reference_Sweep_Label",
            "Is_Reference_Sweep_Selected",
            "Sweep_SC",
            "Sweep_SC_Label",
            "Sweep_SC_Pct",
            "Sweep_Is_Asymptotic",
            "Sweep_Status",
            "Primary_Eligible",
            "Max_Sweep_M_over_N",
            "OuterRep",
        ]
    ].copy()
    out["Metric"] = metric_key
    out["Metric_Pretty"] = metric_pretty(metric_key)
    out["Value"] = pd.to_numeric(df[value_col], errors="coerce")
    return out


def summarise_group_qc(df: pd.DataFrame) -> pd.DataFrame:
    group_cols = [
        "Domain",
        "Subset",
        "Mode",
        "Group",
        "Cutoff_Key",
        "Cutoff_Side",
        "Tail_Display",
        "Retain_Pct",
        "Retain_Ratio",
        "Cutoff_Signed_Pct",
        "Cutoff_Label",
        "Reference_Target_SC",
        "Reference_Sweep_Selected",
        "Reference_Sweep_Label",
        "Is_Reference_Sweep_Selected",
        "Sweep_SC",
        "Sweep_SC_Label",
        "Sweep_SC_Pct",
        "Sweep_Is_Asymptotic",
    ]
    rows: List[Dict[str, object]] = []
    for keys, grp in df.groupby(group_cols, dropna=False, sort=False):
        row = dict(zip(group_cols, keys if isinstance(keys, tuple) else (keys,)))
        row["Sweep_Status"] = collapse_sweep_status(grp["Sweep_Status"], grp["Sweep_Is_Asymptotic"])
        row["Primary_Eligible"] = collapse_primary_eligible(grp["Primary_Eligible"])
        row["Reference_N"] = float(np.nanmedian(grp["Reference_N"].to_numpy(dtype=float)))
        row["Sweep_Planned_M"] = float(np.nanmedian(grp["Sweep_Planned_M"].to_numpy(dtype=float)))
        row["Sweep_M_over_N"] = float(np.nanmedian(grp["Sweep_M_over_N"].to_numpy(dtype=float)))
        row["Uses_Extrapolation"] = bool(pd.Series(grp["Uses_Extrapolation"]).fillna(False).astype(bool).any())
        row["Realized_SC_Model_Median"] = float(np.nanmedian(grp["Realized_SC_Model"].to_numpy(dtype=float)))
        row["Realized_SC_Model_Q025"] = float(np.nanquantile(grp["Realized_SC_Model"].to_numpy(dtype=float), 0.025))
        row["Realized_SC_Model_Q975"] = float(np.nanquantile(grp["Realized_SC_Model"].to_numpy(dtype=float), 0.975))
        row["Realized_SC_Empirical_Median"] = float(np.nanmedian(grp["Realized_SC_Empirical"].to_numpy(dtype=float)))
        row["Mean_f1_Median"] = float(np.nanmedian(grp["Mean_f1"].to_numpy(dtype=float)))
        row["Mean_Realized_N_Median"] = float(np.nanmedian(grp["Mean_Realized_N"].to_numpy(dtype=float)))
        rows.append(row)
    return pd.DataFrame(rows)


def sweep_axis_breaks(values: Iterable[float]) -> List[float]:
    return pick_nearest_available(values, REPRESENTATIVE_SWEEP_TARGETS)


def add_signed_cutoff_display_pct(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    side = out["Cutoff_Side"].astype(str).str.lower()
    retain = pd.to_numeric(out["Retain_Pct"], errors="coerce")
    signed = pd.Series(np.nan, index=out.index, dtype=float)
    signed.loc[side == "none"] = 0.0
    signed.loc[side == "low"] = -(100.0 - retain.loc[side == "low"])
    signed.loc[side == "high"] = 100.0 - retain.loc[side == "high"]
    out[SIGNED_CUTOFF_DISPLAY_COL] = signed
    return out


def signed_cutoff_tick_positions(values: Sequence[float]) -> tuple[List[int], List[str]]:
    vals = [float(v) for v in values]
    positions: List[int] = []
    labels: List[str] = []
    for target in SIGNED_CUTOFF_TICK_TARGETS:
        if not vals:
            continue
        idx = min(range(len(vals)), key=lambda i: abs(vals[i] - target))
        if idx in positions:
            continue
        positions.append(idx)
        labels.append(f"{vals[idx]:.0f}")
    return positions, labels


def fixed_signed_cutoff_display_values(values: Optional[Iterable[float]] = None) -> List[float]:
    """Use a constant display domain so missing cells do not shrink the x-axis."""
    fixed = [float(v) for v in SIGNED_CUTOFF_DISPLAY_VALUES]
    if values is None:
        return fixed
    extras = []
    for value in values:
        try:
            x = float(value)
        except (TypeError, ValueError):
            continue
        if not math.isfinite(x):
            continue
        rounded = round(x)
        if abs(x - rounded) < 1e-9 and int(rounded) in SIGNED_CUTOFF_DISPLAY_VALUES:
            continue
        extras.append(x)
    if extras:
        return sorted(set(fixed + extras))
    return fixed


def apply_signed_cutoff_axes(ax, x_cols: Sequence[float], y_rows: Sequence[float], show_y_labels: bool = True) -> None:
    x_positions, x_labels = signed_cutoff_tick_positions(x_cols)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, rotation=90, fontsize=8)
    if x_cols:
        zero_idx = min(range(len(x_cols)), key=lambda i: abs(float(x_cols[i])))
        ax.axvline(zero_idx, color="#777777", linewidth=0.8)
    y_break_vals = sweep_axis_breaks(y_rows)
    y_positions = [min(range(len(y_rows)), key=lambda i: abs(y_rows[i] - target)) for target in y_break_vals]
    ax.set_yticks(y_positions)
    ax.set_yticklabels([format_sweep_pct_number(y_rows[idx]) for idx in y_positions] if show_y_labels else [], fontsize=8)
    ax.set_xlabel(SIGNED_CUTOFF_AXIS_LABEL)
    ax.set_ylabel("Sample coverage $s$ [%]" if show_y_labels else "")
    for spine in ax.spines.values():
        spine.set_visible(False)


def fig01_display_tick_label(abs_r: int) -> str:
    return str(abs_r) if abs_r % 5 == 0 else ""


def add_fig01_effect_background(ax) -> None:
    bands = [
        (CLIFF_MEDIUM_THRESHOLD, CLIFF_LARGE_THRESHOLD, "#BFF6FA", 0.22),
        (CLIFF_LARGE_THRESHOLD, 1.08, "#55DDE7", 0.18),
        (-CLIFF_LARGE_THRESHOLD, -CLIFF_MEDIUM_THRESHOLD, "#F9C7EF", 0.22),
        (-1.08, -CLIFF_LARGE_THRESHOLD, "#F472D0", 0.18),
    ]
    for ymin, ymax, color, alpha in bands:
        ax.axhspan(ymin, ymax, facecolor=color, alpha=alpha, edgecolor="none", zorder=0)


def add_fig01_cliff_guides(ax) -> None:
    for value, linestyle in [
        (CLIFF_SMALL_THRESHOLD, ":"),
        (CLIFF_MEDIUM_THRESHOLD, "--"),
        (CLIFF_LARGE_THRESHOLD, "-."),
    ]:
        for sign in (-1, 1):
            ax.axhline(sign * value, color="#D62728", linestyle=linestyle, linewidth=1.05, alpha=0.88, zorder=1)


def fig01_tail_label(cutoff_side: object) -> str:
    side = str(cutoff_side).lower()
    if side == "low":
        return "High-tail cutoff r<0"
    if side == "high":
        return "Low-tail cutoff r>0"
    return "No cutoff r=0"


def fig01_abs_r_pct(row: pd.Series) -> int:
    side = str(row.get("Cutoff_Side", "")).lower()
    if side == "none":
        return 0
    retain = pd.to_numeric(pd.Series([row.get("Retain_Pct")]), errors="coerce").iloc[0]
    if not math.isfinite(float(retain)):
        return 0
    return int(max(0, min(100, round(101 - float(retain)))))


def add_fig01_boxplot_layer(ax, plot_df: pd.DataFrame, order_int: Sequence[int]) -> None:
    for tail in FIG01_TAIL_ORDER:
        data: List[np.ndarray] = []
        positions: List[float] = []
        for idx, abs_r in enumerate(order_int):
            values = plot_df[
                (plot_df["Fig01_Tail_Label"].eq(tail)) & (plot_df["Fig01_Abs_R_Pct"].eq(abs_r))
            ]["Value"].dropna().to_numpy(dtype=float)
            if len(values) == 0:
                continue
            data.append(values)
            positions.append(idx + FIG01_TAIL_OFFSETS[tail])
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
            box.set_facecolor(FIG01_TAIL_PALETTE[tail])
            box.set_alpha(1.0)


def add_fig01_mean_markers(ax, plot_df: pd.DataFrame, order_int: Sequence[int]) -> None:
    grouped = (
        plot_df.groupby(["Fig01_Abs_R_Pct", "Fig01_Tail_Label"], observed=True, as_index=False)["Value"]
        .mean()
        .dropna()
    )
    x_pos = {abs_r: idx for idx, abs_r in enumerate(order_int)}
    for _, row in grouped.iterrows():
        tail = str(row["Fig01_Tail_Label"])
        abs_r = int(row["Fig01_Abs_R_Pct"])
        if abs_r not in x_pos:
            continue
        ax.plot(
            x_pos[abs_r] + FIG01_TAIL_OFFSETS.get(tail, 0.0),
            float(row["Value"]),
            marker="x",
            linestyle="None",
            color=FIG01_MEAN_MARKERS.get(tail, "#333333"),
            markersize=5.0,
            markeredgewidth=1.5,
            zorder=5,
        )


def fig01_tail_legend_handles() -> List[object]:
    return [
        Patch(facecolor=FIG01_TAIL_PALETTE["High-tail cutoff r<0"], edgecolor="#555555", label="High-tail cutoff r<0"),
        Patch(facecolor=FIG01_TAIL_PALETTE["Low-tail cutoff r>0"], edgecolor="#555555", label="Low-tail cutoff r>0"),
        Patch(facecolor=FIG01_TAIL_PALETTE["No cutoff r=0"], edgecolor="#555555", label="No cutoff r=0"),
        Line2D([0], [0], marker="x", linestyle="None", color=FIG01_MEAN_MARKERS["High-tail cutoff r<0"], label="Mean r<0"),
        Line2D([0], [0], marker="x", linestyle="None", color=FIG01_MEAN_MARKERS["Low-tail cutoff r>0"], label="Mean r>0"),
        Patch(facecolor="#55DDE7", alpha=0.18, edgecolor="none", label="Positive band"),
        Patch(facecolor="#F472D0", alpha=0.18, edgecolor="none", label="Negative band"),
    ]


def fig01_threshold_handles() -> List[object]:
    return [
        Line2D([0], [0], color="#D62728", linestyle=":", linewidth=1.3, label="Small-or-larger: |delta|>=0.147"),
        Line2D([0], [0], color="#D62728", linestyle="--", linewidth=1.3, label="Medium-or-larger: |delta|>=0.330"),
        Line2D([0], [0], color="#D62728", linestyle="-.", linewidth=1.3, label="Large-or-larger: |delta|>=0.474"),
    ]


def add_fig01_outside_legends(fig, ax) -> None:
    first = ax.legend(
        handles=fig01_tail_legend_handles(),
        loc="upper left",
        bbox_to_anchor=(1.012, 1.0),
        borderaxespad=0.0,
        frameon=False,
        fontsize=8.0,
    )
    ax.add_artist(first)
    ax.legend(
        handles=fig01_threshold_handles(),
        title="Effect-size thresholds",
        loc="lower left",
        bbox_to_anchor=(1.012, 0.0),
        borderaxespad=0.0,
        frameon=False,
        fontsize=8.0,
        title_fontsize=8.5,
    )
    fig.subplots_adjust(right=0.765)


def render_fig01_legend_only(path: Path, dpi: int = 220) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12.0, 0.82))
    ax.axis("off")
    ax.legend(
        handles=fig01_tail_legend_handles() + fig01_threshold_handles(),
        loc="center",
        ncol=len(fig01_tail_legend_handles() + fig01_threshold_handles()),
        frameon=False,
        fontsize=8.5,
        columnspacing=1.25,
        handlelength=2.0,
        borderpad=0.35,
    )
    fig.savefig(path, dpi=dpi, bbox_inches="tight", transparent=True)
    plt.close(fig)


def attach_common_style(fig, title: str, caption: Optional[str] = None) -> None:
    # Stage2 mask review decks use slide-level titles; keep panel PNGs title-free.
    if caption and EMBED_FIGURE_CAPTIONS:
        fig.text(0.01, 0.01, caption, ha="left", va="bottom", fontsize=9)


def plot_reference_bridge_boxplots(df_one: pd.DataFrame, include_legend: bool = True):
    plot_df = df_one[
        (np.abs(df_one["Reference_Target_SC"] - PRIMARY_REFERENCE_TARGET_SC) < PRIMARY_SC_TOL)
        & df_one["Is_Reference_Sweep_Selected"].astype(bool)
        & np.isfinite(df_one["Value"].astype(float))
    ].copy()
    title_txt = make_combo_title(df_one, "bridge boxplots at reference-like sweep")
    if plot_df.empty:
        return empty_plot(title_txt)

    plot_df["Fig01_Abs_R_Pct"] = plot_df.apply(fig01_abs_r_pct, axis=1)
    plot_df["Fig01_Tail_Label"] = plot_df["Cutoff_Side"].map(fig01_tail_label)
    plot_df["Value"] = pd.to_numeric(plot_df["Value"], errors="coerce")
    plot_df = plot_df[np.isfinite(plot_df["Value"])].copy()
    if plot_df.empty:
        return empty_plot(title_txt)

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
    fig._tomato_skip_tight_layout = True
    ax.set_facecolor("#FFFFFF")
    add_fig01_effect_background(ax)
    add_fig01_boxplot_layer(ax, plot_df, FIG01_ABS_R_ORDER)
    add_fig01_cliff_guides(ax)
    add_fig01_mean_markers(ax, plot_df, FIG01_ABS_R_ORDER)
    ax.set_xticks(range(len(FIG01_ABS_R_ORDER)))
    ax.set_xticklabels([fig01_display_tick_label(v) for v in FIG01_ABS_R_ORDER], rotation=90, fontsize=8)
    ax.set_xlim(-0.5, len(FIG01_ABS_R_ORDER) - 0.5)
    ax.set_xlabel("Intensity cutoff |r|", fontsize=11)
    ax.set_ylabel("Effect size δ" if plot_df["Metric"].iloc[0] == "CliffsDelta" else plot_df["Metric_Pretty"].iloc[0], fontsize=11)
    if plot_df["Metric"].iloc[0] == "CliffsDelta":
        ax.set_ylim(-1.08, 1.08)
    else:
        add_neglogp_guides(ax, plot_df["Metric"].iloc[0])
        add_a_guides(ax, plot_df["Metric"].iloc[0], axis="y")
    ax.grid(axis="y", zorder=0)
    ax.grid(axis="x", color="#ECECEC", linewidth=0.4, alpha=0.65, zorder=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    if include_legend:
        add_fig01_outside_legends(fig, ax)
    else:
        fig.tight_layout()
    return fig


def plot_reference_bridge_quantile_lines(summary_one: pd.DataFrame):
    plot_df = summary_one[summary_one["Cutoff_Side"] != "none"].copy()
    ribbon_df = plot_df[np.isfinite(plot_df["q1"]) & np.isfinite(plot_df["q3"])].copy()
    line_df = plot_df[np.isfinite(plot_df["median_value"]) | np.isfinite(plot_df["mean_value"])].copy()
    base_df = summary_one[(summary_one["Cutoff_Side"] == "none") & np.isfinite(summary_one["median_value"])].copy()
    title_txt = make_combo_title(summary_one, "bridge quantile lines at reference-like sweep")
    if line_df.empty:
        return empty_plot(title_txt)

    fig, ax = plt.subplots(figsize=(12, 6))
    for tail_name in ["High-tail", "Low-tail"]:
        tail_rib = ribbon_df[ribbon_df["Tail_Display"] == tail_name].sort_values("Retain_Pct")
        tail_line = line_df[line_df["Tail_Display"] == tail_name].sort_values("Retain_Pct")
        if not tail_rib.empty:
            ax.fill_between(tail_rib["Retain_Pct"], tail_rib["q1"], tail_rib["q3"], color=TAIL_PALETTE[tail_name], alpha=0.18)
        if not tail_line.empty:
            ax.plot(tail_line["Retain_Pct"], tail_line["median_value"], color=TAIL_PALETTE[tail_name], linewidth=1.5, label=tail_name)
            ax.plot(tail_line["Retain_Pct"], tail_line["mean_value"], color=TAIL_PALETTE[tail_name], linewidth=1.0, linestyle="--")
    if not base_df.empty:
        ax.scatter(base_df["Retain_Pct"], base_df["median_value"], color=TAIL_PALETTE["Baseline"], s=24, label="Baseline")
    ax.set_xlabel("Intensity retained (%)")
    ax.set_ylabel(summary_one["Metric_Pretty"].iloc[0])
    ax.set_xticks([1] + list(range(5, 101, 5)))
    ax.legend(loc="upper right")
    add_neglogp_guides(ax, summary_one["Metric"].iloc[0])
    add_a_guides(ax, summary_one["Metric"].iloc[0], axis="y")
    attach_common_style(
        fig,
        title_txt,
        f"Bridge slice evaluated at Sweep_SC = {summary_one['Reference_Sweep_Label'].dropna().iloc[0]}. Solid = median; dashed = mean; shaded ribbon = Q1-Q3.",
    )
    return fig


def plot_reference_bridge_histograms(df_one: pd.DataFrame):
    chosen = pick_representative_retain_pcts(df_one["Retain_Pct"])
    plot_df = df_one[
        (np.abs(df_one["Reference_Target_SC"] - PRIMARY_REFERENCE_TARGET_SC) < PRIMARY_SC_TOL)
        & df_one["Is_Reference_Sweep_Selected"].astype(bool)
        & df_one["Retain_Pct"].isin(chosen)
        & np.isfinite(df_one["Value"].astype(float))
    ].copy()
    title_txt = make_combo_title(df_one, "bridge histograms at reference-like sweep")
    if plot_df.empty:
        return empty_plot(title_txt)

    plot_df["Facet_Label"] = np.where(
        plot_df["Cutoff_Side"] == "high",
        "High " + plot_df["Retain_Pct"].astype(int).astype(str).str.zfill(2) + "%",
        np.where(
            plot_df["Cutoff_Side"] == "low",
            "Low " + plot_df["Retain_Pct"].astype(int).astype(str).str.zfill(2) + "%",
            "Baseline 100%",
        ),
    )
    facets = list(dict.fromkeys(plot_df["Facet_Label"].tolist()))
    ncols = 3
    nrows = math.ceil(len(facets) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(13, 10), squeeze=False)
    for ax, facet in zip(axes.flatten(), facets):
        sub = plot_df[plot_df["Facet_Label"] == facet]
        for tail_name in ["High-tail", "Low-tail", "Baseline"]:
            tail_sub = sub[sub["Tail_Display"] == tail_name]
            if tail_sub.empty:
                continue
            ax.hist(tail_sub["Value"], bins=14, alpha=0.82, color=TAIL_PALETTE[tail_name], edgecolor="white", label=tail_name)
        add_a_guides(ax, plot_df["Metric"].iloc[0], axis="x")
        ax.set_title(facet, fontsize=9)
        ax.set_xlabel(plot_df["Metric_Pretty"].iloc[0])
        ax.set_ylabel("Count")
    for ax in axes.flatten()[len(facets):]:
        ax.set_axis_off()
    handles, labels = axes.flatten()[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center", ncol=3)
    attach_common_style(
        fig,
        title_txt,
        "Representative bar histograms at the bridge slice. This is the stage-2 distribution at Sweep_SC closest to the fixed reference target.",
    )
    return fig


def plot_representative_sweep_boxplots(df_one: pd.DataFrame):
    if sns is None:
        raise RuntimeError("seaborn is required for boxplot export. Use TOMATO_FIGURE_FILTER=06 for heatmap-only export.")
    chosen = pick_nearest_available(df_one["Sweep_SC"], REPRESENTATIVE_SWEEP_TARGETS)
    plot_df = df_one[df_one["Sweep_SC"].isin(chosen) & np.isfinite(df_one["Value"].astype(float))].copy()
    title_txt = make_combo_title(df_one, "paired boxplots across representative Sweep_SC slices")
    if plot_df.empty:
        return empty_plot(title_txt)

    facets = plot_df[["Sweep_SC", "Sweep_SC_Label"]].drop_duplicates().sort_values("Sweep_SC")["Sweep_SC_Label"].tolist()
    ncols = 3
    nrows = math.ceil(len(facets) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(18, 12), squeeze=False)
    for ax, facet in zip(axes.flatten(), facets):
        sub = plot_df[plot_df["Sweep_SC_Label"] == facet].copy()
        sub["Retain_F"] = sub["Retain_Pct"].astype(int).astype(str)
        sns.boxplot(
            data=sub,
            x="Retain_F",
            y="Value",
            hue="Tail_Display",
            order=[str(v) for v in sorted(sub["Retain_Pct"].unique())],
            hue_order=["High-tail", "Low-tail", "Baseline"],
            palette=TAIL_PALETTE,
            dodge=True,
            fliersize=0,
            ax=ax,
        )
        keep_ticks = set(retain_breaks(sub["Retain_Pct"]))
        labels = [tick.get_text() for tick in ax.get_xticklabels()]
        ax.set_xticklabels([lab if lab in keep_ticks else "" for lab in labels], rotation=90, fontsize=7)
        ax.set_title(facet, fontsize=9)
        ax.set_xlabel("Intensity retained (%)")
        ax.set_ylabel(sub["Metric_Pretty"].iloc[0])
        add_neglogp_guides(ax, sub["Metric"].iloc[0])
        add_a_guides(ax, sub["Metric"].iloc[0], axis="y")
        if ax.get_legend() is not None:
            ax.get_legend().remove()
    for ax in axes.flatten()[len(facets):]:
        ax.set_axis_off()
    handles, labels = axes.flatten()[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center", ncol=3)
    attach_common_style(fig, title_txt, "Panels are representative stage-2 sample-coverage slices across the full Sweep_SC axis.")
    return fig


def plot_sweep_sc_support_lines(summary_one: pd.DataFrame):
    chosen = pick_representative_retain_pcts(summary_one["Retain_Pct"])
    plot_df = summary_one[
        summary_one["Retain_Pct"].isin(chosen)
        & (summary_one["Cutoff_Side"] != "none")
        & np.isfinite(summary_one["median_value"])
    ].copy()
    title_txt = make_combo_title(summary_one, "full Sweep_SC support lines")
    if plot_df.empty:
        return empty_plot(title_txt)

    plot_df["Facet_Label"] = plot_df["Retain_Pct"].astype(int).map(lambda x: f"Retain {x:03d}%")
    facets = list(dict.fromkeys(plot_df["Facet_Label"].tolist()))
    base_template = pd.DataFrame({"Facet_Label": facets})
    base_df = summary_one[(summary_one["Cutoff_Side"] == "none") & np.isfinite(summary_one["median_value"])][["Sweep_SC", "median_value"]].copy()
    base_df = base_template.assign(_key=1).merge(base_df.assign(_key=1), on="_key").drop(columns="_key")

    ncols = 3
    nrows = math.ceil(len(facets) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 10), squeeze=False)
    for ax, facet in zip(axes.flatten(), facets):
        sub = plot_df[plot_df["Facet_Label"] == facet].sort_values("Sweep_SC")
        for tail_name in ["High-tail", "Low-tail"]:
            tail_sub = sub[sub["Tail_Display"] == tail_name]
            if tail_sub.empty:
                continue
            ax.plot(tail_sub["Sweep_SC"], tail_sub["median_value"], color=TAIL_PALETTE[tail_name], linewidth=1.2, label=tail_name)
            for status, mark in SWEEP_SHAPES.items():
                status_sub = tail_sub[tail_sub["Sweep_Status"] == status]
                if status_sub.empty:
                    continue
                ax.scatter(status_sub["Sweep_SC"], status_sub["median_value"], color=TAIL_PALETTE[tail_name], marker=mark, s=18)
        base_sub = base_df[base_df["Facet_Label"] == facet]
        if not base_sub.empty:
            ax.plot(base_sub["Sweep_SC"], base_sub["median_value"], color=TAIL_PALETTE["Baseline"], linestyle="--", linewidth=1.0)
        ax.set_title(facet, fontsize=9)
        ax.set_xlabel("Sample coverage $s$")
        ax.set_ylabel(sub["Metric_Pretty"].iloc[0])
        breaks = sweep_axis_breaks(summary_one["Sweep_SC"])
        ax.set_xticks(breaks)
        ax.set_xticklabels([format_sweep_pct(v) for v in breaks], rotation=90)
        add_neglogp_guides(ax, summary_one["Metric"].iloc[0])
        add_a_guides(ax, summary_one["Metric"].iloc[0], axis="y")
    for ax in axes.flatten()[len(facets):]:
        ax.set_axis_off()
    attach_common_style(
        fig,
        title_txt,
        "This is the main stage-2 view. Dashed grey is the no-cutoff baseline. The Asymptotic tick is the separate endpoint.",
    )
    return fig


def plot_sweep_sc_heatmap(cell_one: pd.DataFrame):
    plot_df = cell_one[np.isfinite(cell_one["Value"].astype(float))].copy()
    title_txt = make_combo_title(cell_one, "sample coverage x intensity cutoff heatmap")
    if plot_df.empty:
        return empty_plot(title_txt)

    plot_df = add_signed_cutoff_display_pct(plot_df)
    plot_df = plot_df[np.isfinite(plot_df[SIGNED_CUTOFF_DISPLAY_COL].astype(float))].copy()
    if plot_df.empty:
        return empty_plot(title_txt, "No rows remain after intensity cutoff r display-axis conversion.")

    # Keep the two mask panels close to the Period/Annual overlay geometry used
    # in review boards; otherwise Combo6yr becomes much wider and visibly
    # smaller after board placement.
    fig = plt.figure(figsize=(12.8, 8.2))
    fig._tomato_skip_tight_layout = True
    gs = fig.add_gridspec(
        nrows=1,
        ncols=4,
        width_ratios=[1.0, 0.035, 1.0, 0.035],
        wspace=0.10,
    )
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 2])]
    cbar_axes = [fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 3])]
    metric_key = str(plot_df["Metric"].dropna().iloc[0]) if "Metric" in plot_df.columns and not plot_df["Metric"].dropna().empty else ""
    piv = plot_df.pivot_table(
        index="Sweep_SC",
        columns=SIGNED_CUTOFF_DISPLAY_COL,
        values="Value",
        aggfunc="mean",
    ).sort_index()
    piv = piv.reindex(fixed_signed_cutoff_display_values(piv.columns), axis=1)
    x_cols = [float(x) for x in piv.columns]
    y_rows = [float(y) for y in piv.index]
    panels = [
        ("Positive effect surface", lambda values: np.where(np.isfinite(values), np.maximum(values, 0.0), np.nan), POSITIVE_SURFACE_CMAP),
        ("Negative effect surface", lambda values: np.where(np.isfinite(values), np.maximum(-values, 0.0), np.nan), NEGATIVE_SURFACE_CMAP),
    ]
    for panel_idx, (ax, cbar_ax, (panel_title, mapper, cmap)) in enumerate(zip(axes, cbar_axes, panels)):
        arr = piv.to_numpy(dtype=float)
        plot_arr = mapper(arr)
        im = ax.imshow(plot_arr, aspect="auto", origin="lower", vmin=0.0, vmax=1.0, cmap=cmap)
        ax.set_title("")
        apply_signed_cutoff_axes(ax, x_cols, y_rows, show_y_labels=(panel_idx == 0))
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.ax.tick_params(labelsize=8)
    attach_common_style(
        fig,
        title_txt,
        "Full two-stage map. Intensity cutoff r uses 0% for no cutoff, negative for high-intensity cutoff, and positive for low-intensity cutoff.",
    )
    fig.subplots_adjust(left=0.08, right=0.93, top=0.86, bottom=0.18, wspace=0.08)
    return fig


def plot_sweep_sc_large_effect_mask(cell_one: pd.DataFrame):
    plot_df = cell_one[np.isfinite(cell_one["Value"].astype(float))].copy()
    title_txt = make_combo_title(cell_one, "medium/large mask over sample coverage x intensity cutoff")
    if plot_df.empty:
        return empty_plot(title_txt)
    metric_key = str(plot_df["Metric"].dropna().iloc[0]) if "Metric" in plot_df.columns and not plot_df["Metric"].dropna().empty else ""
    if metric_key != "CliffsDelta":
        return empty_plot(f"{title_txt}\n06_2 medium/large mask is defined only for Cliff's delta.")

    plot_df = add_signed_cutoff_display_pct(plot_df)
    plot_df = plot_df[np.isfinite(plot_df[SIGNED_CUTOFF_DISPLAY_COL].astype(float))].copy()
    if plot_df.empty:
        return empty_plot(title_txt, "No rows remain after intensity cutoff r display-axis conversion.")

    fig = plt.figure(figsize=(15, 6.4))
    fig._tomato_skip_tight_layout = True
    gs = fig.add_gridspec(nrows=1, ncols=2, wspace=0.10)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]

    piv = plot_df.pivot_table(
        index="Sweep_SC",
        columns=SIGNED_CUTOFF_DISPLAY_COL,
        values="Value",
        aggfunc="mean",
    ).sort_index()
    piv = piv.reindex(fixed_signed_cutoff_display_values(piv.columns), axis=1)
    values = piv.to_numpy(dtype=float)
    x_cols = [float(x) for x in piv.columns]
    y_rows = [float(y) for y in piv.index]

    panel_specs = [
        (r"Effect size $\delta > 0$", POSITIVE_MASK_CMAP, lambda arr: (arr >= 0, arr >= CLIFF_MEDIUM_THRESHOLD, arr >= CLIFF_LARGE_THRESHOLD)),
        (r"Effect size $\delta < 0$", NEGATIVE_MASK_CMAP, lambda arr: (arr <= 0, arr <= -CLIFF_MEDIUM_THRESHOLD, arr <= -CLIFF_LARGE_THRESHOLD)),
    ]
    for panel_idx, (ax, (panel_title, cmap, classifier)) in enumerate(zip(axes, panel_specs)):
        sign_ok, medium_ok, large_ok = classifier(values)
        mask_arr = np.full(values.shape, np.nan, dtype=float)
        finite = np.isfinite(values)
        mask_arr[finite] = 0.0
        mask_arr[finite & sign_ok & medium_ok] = 1.0
        mask_arr[finite & sign_ok & large_ok] = 2.0
        ax.imshow(mask_arr, aspect="auto", origin="lower", cmap=cmap, norm=SIGNED_MASK_NORM)
        ax.set_title(panel_title, fontsize=9)
        apply_signed_cutoff_axes(ax, x_cols, y_rows, show_y_labels=(panel_idx == 0))

    attach_common_style(fig, title_txt)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.84, bottom=0.18, wspace=0.10)
    return fig


def export_heatmap_only_from_cells(cells_df: pd.DataFrame) -> None:
    if cells_df.empty:
        raise RuntimeError("No rows remain for heatmap-only export.")
    cells_df = cells_df.copy()
    if "Metric_Pretty" not in cells_df.columns:
        cells_df["Metric_Pretty"] = cells_df["Metric"].map(metric_pretty)
    if "Estimate_Pretty" not in cells_df.columns:
        cells_df["Estimate_Pretty"] = cells_df["Estimate_Definition"].map(estimate_pretty)
    if "Analysis_Tier" not in cells_df.columns:
        cells_df["Analysis_Tier"] = cells_df["Q_Label"].map(analysis_tier)
    if "Q0_Caution" not in cells_df.columns:
        cells_df["Q0_Caution"] = [q0_caution_text(q, est) for q, est in zip(cells_df["Q_Label"], cells_df["Estimate_Definition"])]
    if "Tail_Display" not in cells_df.columns:
        cells_df["Tail_Display"] = cells_df["Cutoff_Side"].map(tail_display)

    manifest_cols = [
        "Analysis_Tier",
        "Domain",
        "Subset",
        "Mode",
        "Diversity_Order",
        "Q_Label",
        "Estimate_Definition",
        "Estimate_Pretty",
        "Q0_Caution",
        "Metric",
        "Metric_Pretty",
        "Reference_Target_SC",
        "Reference_Sweep_Selected",
        "Reference_Sweep_Label",
    ]
    combo_manifest = (
        cells_df[manifest_cols]
        .drop_duplicates()
        .sort_values(["Analysis_Tier", "Domain", "Subset", "Mode", "Diversity_Order", "Metric"])
        .reset_index(drop=True)
    )
    combo_manifest["combo_id"] = np.arange(1, len(combo_manifest) + 1)
    combo_manifest["output_rel_dir"] = combo_manifest.apply(
        lambda row: str(
            Path("figures")
            / clean_id(row["Analysis_Tier"])
            / clean_id(row["Domain"])
            / clean_id(row["Subset"])
            / clean_id(row["Mode"])
            / clean_id(row["Q_Label"])
            / clean_id(row["Estimate_Definition"])
            / clean_id(row["Metric"])
        ),
        axis=1,
    )
    if COMBO_START > 0:
        combo_end = COMBO_END if COMBO_END > 0 else int(combo_manifest["combo_id"].max())
        combo_manifest["export_selected"] = (
            (combo_manifest["combo_id"] >= COMBO_START)
            & (combo_manifest["combo_id"] <= combo_end)
        )
    else:
        combo_manifest["export_selected"] = (
            True if math.isinf(MAX_COMBOS_TO_EXPORT) else combo_manifest["combo_id"] <= MAX_COMBOS_TO_EXPORT
        )

    combo_manifest.to_csv(DIR_CSV / "combo_manifest.csv", index=False)
    cells_df.to_csv(DIR_CSV / "comparison_cells_filtered.csv", index=False)
    q0_note = [
        "q = 0 caution for v36.9 recommended outputs",
        "- q0 / Chao1_asymptotic is retained, but should be interpreted as a cautious lower-bound-oriented reference.",
        "- q0 / Observed_Standardized is empirical and is not an asymptotic endpoint.",
        "- Main recommended figures remain q1 and q2; q0 is exported alongside them with this note.",
    ]
    (DIR_CSV / "q0_caution_note.txt").write_text("\n".join(q0_note) + "\n", encoding="utf-8")

    if WRITE_COMBO_PAYLOADS:
        for _, one in combo_manifest[combo_manifest["export_selected"]].iterrows():
            payload_dir = ensure_dir(combo_payload_dir(int(one["combo_id"])))
            key_mask_cells = pd.Series(True, index=cells_df.index)
            key_mask_outer = pd.Series(True, index=outer_wide.index)
            for col in combo_key_cols:
                key_mask_cells &= cells_df[col] == one[col]
                key_mask_outer &= outer_wide[col] == one[col]
            cells_df[key_mask_cells].to_csv(payload_dir / "cells.csv", index=False)
            materialize_outer_metric(outer_wide[key_mask_outer].copy(), one["Metric"]).to_csv(
                payload_dir / "outer.csv",
                index=False,
            )

    combos_to_export = combo_manifest[combo_manifest["export_selected"]].copy()
    if SKIP_FIGURE_EXPORT:
        combos_to_export = combos_to_export.iloc[0:0].copy()
    if should_export_figure("01") and not combos_to_export.empty:
        render_fig01_legend_only(DIR_FIG / "figure01_boxplot_legend_style2_horizontal_single_row.png")
    for _, one in combos_to_export.iterrows():
        mask = (
            (cells_df["Analysis_Tier"] == one["Analysis_Tier"])
            & (cells_df["Domain"] == one["Domain"])
            & (cells_df["Subset"] == one["Subset"])
            & (cells_df["Mode"] == one["Mode"])
            & (cells_df["Diversity_Order"] == one["Diversity_Order"])
            & (cells_df["Q_Label"] == one["Q_Label"])
            & (cells_df["Estimate_Definition"] == one["Estimate_Definition"])
            & (cells_df["Metric"] == one["Metric"])
        )
        cell_one = cells_df[mask].copy()
        outdir_one = ensure_dir(
            DIR_FIG
            / clean_id(one["Analysis_Tier"])
            / clean_id(one["Domain"])
            / clean_id(one["Subset"])
            / clean_id(one["Mode"])
            / clean_id(one["Q_Label"])
            / clean_id(one["Estimate_Definition"])
            / clean_id(one["Metric"])
        )
        p6 = plot_sweep_sc_heatmap(cell_one)
        save_plot_pdf_png(p6, outdir_one / "06_sweep_sc_heatmap.pdf", outdir_one / "06_sweep_sc_heatmap.png", 15, 6.4)

    run_note = [
        f"RUN_DIR: {RUN_DIR}",
        f"RUN_DIR_SOURCE: {RUN_DIR_SOURCE}",
        f"CELLS_CSV: {CELLS_CSV}",
        f"OUT_DIR: {OUT_DIR}",
        f"ROOT_DIR: {ROOT_DIR}",
        f"DEVICE: {DEVICE}",
        "FIGURE_FILTER: 06",
        "",
        "Interpretation notes:",
        "- Heatmap-only downstream export used TOMATO_FIGURE_FILTER=06; outer and group CSVs were not read.",
        f"- PRIMARY_REFERENCE_TARGET_SC = {PRIMARY_REFERENCE_TARGET_SC}",
        "- Sample coverage s is the main interpretation axis. The Asymptotic row is a separate endpoint.",
        "- 06_sweep_sc_heatmap uses intensity cutoff r as a signed display coordinate: 0% is no cutoff, negative values cut high-intensity mass, and positive values cut low-intensity mass.",
        "- Positive and negative Cliff's delta surfaces are shown in separate panels to align with the overlay convention.",
        "- A = P(Syneco > Conv) on the probability scale; Cliff's delta = 2A - 1.",
        (
            f"- Exported combos: {len(combos_to_export)} / {len(combo_manifest)}"
            if math.isinf(MAX_COMBOS_TO_EXPORT)
            else f"- Exported combos: {len(combos_to_export)} / {len(combo_manifest)} (limited by TOMATO_MAX_COMBOS)"
        ),
    ]
    (OUT_DIR / "run_note.txt").write_text("\n".join(run_note) + "\n", encoding="utf-8")
    log(f"[SUCCESS] Heatmap-only output written to: {OUT_DIR}")


def export_mask_only_from_cells(cells_df: pd.DataFrame) -> None:
    if cells_df.empty:
        raise RuntimeError("No rows remain for 06_2-only export.")
    manifest_cols = [
        "Analysis_Tier",
        "Domain",
        "Subset",
        "Mode",
        "Diversity_Order",
        "Q_Label",
        "Estimate_Definition",
        "Estimate_Pretty",
        "Q0_Caution",
        "Metric",
        "Metric_Pretty",
        "Reference_Target_SC",
        "Reference_Sweep_Selected",
        "Reference_Sweep_Label",
    ]
    combo_manifest = (
        cells_df[manifest_cols]
        .drop_duplicates()
        .sort_values(["Analysis_Tier", "Domain", "Subset", "Mode", "Diversity_Order", "Metric"])
        .reset_index(drop=True)
    )
    combo_manifest["combo_id"] = np.arange(1, len(combo_manifest) + 1)
    combo_manifest["output_rel_dir"] = combo_manifest.apply(
        lambda row: str(
            Path("figures")
            / clean_id(row["Analysis_Tier"])
            / clean_id(row["Domain"])
            / clean_id(row["Subset"])
            / clean_id(row["Mode"])
            / clean_id(row["Q_Label"])
            / clean_id(row["Estimate_Definition"])
            / clean_id(row["Metric"])
        ),
        axis=1,
    )
    if COMBO_START > 0:
        combo_end = COMBO_END if COMBO_END > 0 else int(combo_manifest["combo_id"].max())
        combo_manifest["export_selected"] = (
            (combo_manifest["combo_id"] >= COMBO_START)
            & (combo_manifest["combo_id"] <= combo_end)
        )
    else:
        combo_manifest["export_selected"] = (
            True if math.isinf(MAX_COMBOS_TO_EXPORT) else combo_manifest["combo_id"] <= MAX_COMBOS_TO_EXPORT
        )

    combo_manifest.to_csv(DIR_CSV / "combo_manifest.csv", index=False)
    cells_df.to_csv(DIR_CSV / "comparison_cells_filtered.csv", index=False)
    combos_to_export = combo_manifest[combo_manifest["export_selected"]].copy()
    if SKIP_FIGURE_EXPORT:
        combos_to_export = combos_to_export.iloc[0:0].copy()
    combo_key_cols = [
        "Analysis_Tier",
        "Domain",
        "Subset",
        "Mode",
        "Diversity_Order",
        "Q_Label",
        "Estimate_Definition",
        "Metric",
    ]
    for _, one in combos_to_export.iterrows():
        cell_one = filter_by_combo_keys(cells_df, one, combo_key_cols)
        outdir_one = ensure_dir(
            DIR_FIG
            / clean_id(one["Analysis_Tier"])
            / clean_id(one["Domain"])
            / clean_id(one["Subset"])
            / clean_id(one["Mode"])
            / clean_id(one["Q_Label"])
            / clean_id(one["Estimate_Definition"])
            / clean_id(one["Metric"])
        )
        p6_2 = plot_sweep_sc_large_effect_mask(cell_one)
        save_plot_pdf_png(
            p6_2,
            outdir_one / "06_2_sweep_sc_large_effect_mask.pdf",
            outdir_one / "06_2_sweep_sc_large_effect_mask.png",
            15,
            6.4,
        )

    run_note = [
        f"RUN_DIR: {RUN_DIR}",
        f"RUN_DIR_SOURCE: {RUN_DIR_SOURCE}",
        f"CELLS_CSV: {CELLS_CSV}",
        f"OUT_DIR: {OUT_DIR}",
        f"ROOT_DIR: {ROOT_DIR}",
        f"DEVICE: {DEVICE}",
        "FIGURE_FILTER: 06_2",
        "",
        "Interpretation notes:",
        "- 06_2-only downstream export used comparison_cells_filtered.csv and did not read outer or group CSVs.",
        "- This fast path is intended for 9-panel board production; 01 and 06 can be rendered later from retained CSVs.",
        "- Sample coverage s is the main interpretation axis. The Asymptotic row is a separate endpoint.",
        "- 06_2_sweep_sc_large_effect_mask uses intensity cutoff r as a signed display coordinate.",
        (
            f"- Exported combos: {len(combos_to_export)} / {len(combo_manifest)}"
            if math.isinf(MAX_COMBOS_TO_EXPORT)
            else f"- Exported combos: {len(combos_to_export)} / {len(combo_manifest)} (limited by TOMATO_MAX_COMBOS)"
        ),
    ]
    (OUT_DIR / "run_note.txt").write_text("\n".join(run_note) + "\n", encoding="utf-8")
    log(f"[SUCCESS] 06_2-only output written to: {OUT_DIR}")


def export_reference_and_heatmap_from_cells_outer(cells_df: pd.DataFrame, outer_wide: pd.DataFrame) -> None:
    if cells_df.empty or outer_wide.empty:
        raise RuntimeError("No rows remain for 01/06-only export.")
    manifest_cols = [
        "Analysis_Tier",
        "Domain",
        "Subset",
        "Mode",
        "Diversity_Order",
        "Q_Label",
        "Estimate_Definition",
        "Estimate_Pretty",
        "Q0_Caution",
        "Metric",
        "Metric_Pretty",
        "Reference_Target_SC",
        "Reference_Sweep_Selected",
        "Reference_Sweep_Label",
    ]
    combo_manifest = (
        cells_df[manifest_cols]
        .drop_duplicates()
        .sort_values(["Analysis_Tier", "Domain", "Subset", "Mode", "Diversity_Order", "Metric"])
        .reset_index(drop=True)
    )
    combo_manifest["combo_id"] = np.arange(1, len(combo_manifest) + 1)
    combo_manifest["output_rel_dir"] = combo_manifest.apply(
        lambda row: str(
            Path("figures")
            / clean_id(row["Analysis_Tier"])
            / clean_id(row["Domain"])
            / clean_id(row["Subset"])
            / clean_id(row["Mode"])
            / clean_id(row["Q_Label"])
            / clean_id(row["Estimate_Definition"])
            / clean_id(row["Metric"])
        ),
        axis=1,
    )
    if COMBO_START > 0:
        combo_end = COMBO_END if COMBO_END > 0 else int(combo_manifest["combo_id"].max())
        combo_manifest["export_selected"] = (
            (combo_manifest["combo_id"] >= COMBO_START)
            & (combo_manifest["combo_id"] <= combo_end)
        )
    else:
        combo_manifest["export_selected"] = (
            True if math.isinf(MAX_COMBOS_TO_EXPORT) else combo_manifest["combo_id"] <= MAX_COMBOS_TO_EXPORT
        )

    combo_key_cols = [
        "Analysis_Tier",
        "Domain",
        "Subset",
        "Mode",
        "Diversity_Order",
        "Q_Label",
        "Estimate_Definition",
    ]
    combo_manifest.to_csv(DIR_CSV / "combo_manifest.csv", index=False)
    cells_df.to_csv(DIR_CSV / "comparison_cells_filtered.csv", index=False)
    q0_note = [
        "q = 0 caution for v36.9 recommended outputs",
        "- q0 / Chao1_asymptotic is retained, but should be interpreted as a cautious lower-bound-oriented reference.",
        "- q0 / Observed_Standardized is empirical and is not an asymptotic endpoint.",
        "- Main recommended figures remain q1 and q2; q0 is exported alongside them with this note.",
    ]
    (DIR_CSV / "q0_caution_note.txt").write_text("\n".join(q0_note) + "\n", encoding="utf-8")

    if WRITE_COMBO_PAYLOADS:
        for _, one in combo_manifest[combo_manifest["export_selected"]].iterrows():
            payload_dir = ensure_dir(combo_payload_dir(int(one["combo_id"])))
            key_mask_cells = pd.Series(True, index=cells_df.index)
            key_mask_outer = pd.Series(True, index=outer_wide.index)
            for col in combo_key_cols:
                key_mask_cells &= cells_df[col] == one[col]
                key_mask_outer &= outer_wide[col] == one[col]
            cells_df[key_mask_cells].to_csv(payload_dir / "cells.csv", index=False)
            materialize_outer_metric(outer_wide[key_mask_outer].copy(), one["Metric"]).to_csv(
                payload_dir / "outer.csv",
                index=False,
            )

    combos_to_export = combo_manifest[combo_manifest["export_selected"]].copy()
    if SKIP_FIGURE_EXPORT:
        combos_to_export = combos_to_export.iloc[0:0].copy()
    if should_export_figure("01") and not combos_to_export.empty:
        render_fig01_legend_only(DIR_FIG / "figure01_boxplot_legend_style2_horizontal_single_row.png")
    for _, one in combos_to_export.iterrows():
        key_mask_cells = pd.Series(True, index=cells_df.index)
        key_mask_outer = pd.Series(True, index=outer_wide.index)
        for col in combo_key_cols:
            key_mask_cells &= cells_df[col] == one[col]
            key_mask_outer &= outer_wide[col] == one[col]
        cell_one = cells_df[key_mask_cells].copy()
        outer_one = materialize_outer_metric(outer_wide[key_mask_outer].copy(), one["Metric"])
        outdir_one = ensure_dir(
            DIR_FIG
            / clean_id(one["Analysis_Tier"])
            / clean_id(one["Domain"])
            / clean_id(one["Subset"])
            / clean_id(one["Mode"])
            / clean_id(one["Q_Label"])
            / clean_id(one["Estimate_Definition"])
            / clean_id(one["Metric"])
        )
        if should_export_figure("01"):
            p1 = plot_reference_bridge_boxplots(outer_one, include_legend=True)
            save_plot_pdf_png(
                p1,
                outdir_one / "01_reference_bridge_boxplots_paired.pdf",
                outdir_one / "01_reference_bridge_boxplots_paired.png",
                20.4,
                6.6,
            )
            p1_nolegend = plot_reference_bridge_boxplots(outer_one, include_legend=False)
            save_plot_pdf_png(
                p1_nolegend,
                outdir_one / "01_reference_bridge_boxplots_paired_nolegend.pdf",
                outdir_one / "01_reference_bridge_boxplots_paired_nolegend.png",
                17.6,
                6.6,
            )
        if should_export_figure("06"):
            p6 = plot_sweep_sc_heatmap(cell_one)
            save_plot_pdf_png(p6, outdir_one / "06_sweep_sc_heatmap.pdf", outdir_one / "06_sweep_sc_heatmap.png", 15, 6.4)
        if should_export_figure("06_2"):
            p6_2 = plot_sweep_sc_large_effect_mask(cell_one)
            save_plot_pdf_png(p6_2, outdir_one / "06_2_sweep_sc_large_effect_mask.pdf", outdir_one / "06_2_sweep_sc_large_effect_mask.png", 15, 6.4)

    run_note = [
        f"RUN_DIR: {RUN_DIR}",
        f"RUN_DIR_SOURCE: {RUN_DIR_SOURCE}",
        f"CELLS_CSV: {CELLS_CSV}",
        f"OUTER_CSV: {OUTER_CSV}",
        f"OUT_DIR: {OUT_DIR}",
        f"ROOT_DIR: {ROOT_DIR}",
        f"DEVICE: {DEVICE}",
        f"FIGURE_FILTER: {','.join(sorted(FIGURE_FILTER))}",
        "",
        "Interpretation notes:",
        "- Lightweight 01/06/06_2 downstream export path was used; group and sweep support CSVs were not read.",
        f"- PRIMARY_REFERENCE_TARGET_SC = {PRIMARY_REFERENCE_TARGET_SC}",
        "- Sample coverage s is the main interpretation axis. The Asymptotic row is a separate endpoint.",
        "- 06_sweep_sc_heatmap and 06_2_sweep_sc_large_effect_mask use intensity cutoff r as a signed display coordinate: 0% is no cutoff, negative values cut high-intensity mass, and positive values cut low-intensity mass.",
        "- Positive and negative Cliff's delta surfaces/masks are shown in separate panels to align with the overlay convention.",
        "- A = P(Syneco > Conv) on the probability scale; Cliff's delta = 2A - 1.",
        (
            f"- Exported combos: {len(combos_to_export)} / {len(combo_manifest)}"
            if math.isinf(MAX_COMBOS_TO_EXPORT)
            else f"- Exported combos: {len(combos_to_export)} / {len(combo_manifest)} (limited by TOMATO_MAX_COMBOS)"
        ),
    ]
    (OUT_DIR / "run_note.txt").write_text("\n".join(run_note) + "\n", encoding="utf-8")
    log(f"[SUCCESS] 01/06/06_2 lightweight output written to: {OUT_DIR}")


def plot_sweep_m_over_n_lines(group_one: pd.DataFrame):
    chosen = pick_representative_retain_pcts(group_one["Retain_Pct"])
    plot_df = group_one[
        group_one["Retain_Pct"].isin(chosen)
        & (group_one["Tail_Display"] != "Baseline")
        & np.isfinite(group_one["Sweep_M_over_N"])
    ].copy()
    title_txt = make_group_title(group_one, "supplemental m / n support lines")
    if plot_df.empty:
        return empty_plot(title_txt)

    plot_df["Facet_Label"] = plot_df["Retain_Pct"].astype(int).map(lambda x: f"Retain {x:03d}%")
    facets = list(dict.fromkeys(plot_df["Facet_Label"].tolist()))
    ncols = 3
    nrows = math.ceil(len(facets) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 9), squeeze=False)
    group_colors = {"Conv": "#4C72B0", "Syneco": "#55A868"}
    tail_styles = {"High-tail": "-", "Low-tail": "--"}
    for ax, facet in zip(axes.flatten(), facets):
        sub = plot_df[plot_df["Facet_Label"] == facet].sort_values("Sweep_SC")
        ax.axhline(1.0, linestyle="--", color="#6E6E6E", linewidth=0.8)
        ax.axhline(2.0, linestyle=":", color="#AA3333", linewidth=0.8)
        for group in ["Conv", "Syneco"]:
            for tail_name in ["High-tail", "Low-tail"]:
                grp = sub[(sub["Group"] == group) & (sub["Tail_Display"] == tail_name)]
                if grp.empty:
                    continue
                ax.plot(grp["Sweep_SC"], grp["Sweep_M_over_N"], color=group_colors[group], linestyle=tail_styles[tail_name], linewidth=1.2, label=f"{group} / {tail_name}")
                ax.scatter(grp["Sweep_SC"], grp["Sweep_M_over_N"], color=group_colors[group], s=18)
        breaks = sweep_axis_breaks(plot_df["Sweep_SC"])
        ax.set_xticks(breaks)
        ax.set_xticklabels([format_sweep_pct(v) for v in breaks], rotation=90)
        ax.set_title(facet, fontsize=9)
        ax.set_xlabel("Sample coverage $s$")
        ax.set_ylabel("Median planned m / reference n")
    for ax in axes.flatten()[len(facets):]:
        ax.set_axis_off()
    handles, labels = axes.flatten()[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center", ncol=2)
    attach_common_style(fig, title_txt, "Supplemental effort-based support. Horizontal lines mark 1x and 2x.")
    return fig


def plot_realized_sc_qc(group_one: pd.DataFrame):
    chosen = pick_representative_retain_pcts(group_one["Retain_Pct"])
    plot_df = group_one[
        group_one["Retain_Pct"].isin(chosen)
        & (group_one["Tail_Display"] != "Baseline")
        & np.isfinite(group_one["Realized_SC_Model_Median"])
    ].copy()
    title_txt = make_group_title(group_one, "realized SC QC across Sweep_SC")
    if plot_df.empty:
        return empty_plot(title_txt)

    plot_df["Facet_Label"] = plot_df["Retain_Pct"].astype(int).map(lambda x: f"Retain {x:03d}%")
    facets = list(dict.fromkeys(plot_df["Facet_Label"].tolist()))
    ncols = 3
    nrows = math.ceil(len(facets) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 9), squeeze=False)
    group_colors = {"Conv": "#4C72B0", "Syneco": "#55A868"}
    tail_styles = {"High-tail": "-", "Low-tail": "--"}
    for ax, facet in zip(axes.flatten(), facets):
        sub = plot_df[plot_df["Facet_Label"] == facet].sort_values("Sweep_SC")
        ax.plot([0, 1], [0, 1], linestyle="--", color="#9A9A9A", linewidth=0.8)
        for group in ["Conv", "Syneco"]:
            for tail_name in ["High-tail", "Low-tail"]:
                grp = sub[(sub["Group"] == group) & (sub["Tail_Display"] == tail_name)]
                if grp.empty:
                    continue
                ax.plot(grp["Sweep_SC"], grp["Realized_SC_Model_Median"], color=group_colors[group], linestyle=tail_styles[tail_name], linewidth=1.2, label=f"{group} / {tail_name}")
                ax.scatter(grp["Sweep_SC"], grp["Realized_SC_Model_Median"], color=group_colors[group], s=18)
        breaks = sweep_axis_breaks(plot_df["Sweep_SC"])
        ax.set_xticks(breaks)
        ax.set_xticklabels([format_sweep_pct(v) for v in breaks], rotation=90)
        ax.set_yticks(breaks)
        ax.set_yticklabels([format_sweep_pct(v) for v in breaks])
        ax.set_title(facet, fontsize=9)
        ax.set_xlabel("Sample coverage $s$")
        ax.set_ylabel("Realized sample coverage (median over OuterRep)")
    for ax in axes.flatten()[len(facets):]:
        ax.set_axis_off()
    handles, labels = axes.flatten()[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center", ncol=2)
    attach_common_style(
        fig,
        title_txt,
        "Diagonal line marks perfect target matching. Deviations show where the empirical realized coverage drifts from the requested Stage2 sample coverage.",
    )
    return fig


def main() -> None:
    if not SKIP_SUMMARY_WRITE and not (CELLS_CSV.exists() and OUTER_CSV.exists() and GROUP_CSV.exists()):
        raise RuntimeError("Expected v36.9 upstream CSVs are missing.")

    global DOMAIN_FILTER, SUBSET_FILTER, MODE_FILTER, Q_FILTER, ESTIMATE_FILTER, METRIC_FILTER
    DOMAIN_FILTER = expand_filter_aliases(parse_env_list("TOMATO_DOMAIN_FILTER", DOMAIN_FILTER), DOMAIN_FILTER_ALIASES)
    SUBSET_FILTER = parse_env_list("TOMATO_SUBSET_FILTER", SUBSET_FILTER)
    MODE_FILTER = expand_filter_aliases(parse_env_list("TOMATO_MODE_FILTER", MODE_FILTER), MODE_FILTER_ALIASES)
    Q_FILTER = parse_env_list("TOMATO_Q_FILTER", Q_FILTER)
    ESTIMATE_FILTER = parse_env_list("TOMATO_ESTIMATE_FILTER", ESTIMATE_FILTER)
    METRIC_FILTER = parse_env_list("TOMATO_METRIC_FILTER", METRIC_FILTER)

    log(f"[INFO] ROOT_DIR = {ROOT_DIR}")
    log(f"[INFO] DEVICE = {DEVICE}")
    log(f"[INFO] RUN_DIR = {RUN_DIR}")
    log(f"[INFO] RUN_DIR_SOURCE = {RUN_DIR_SOURCE}")
    log(f"[INFO] MAX_COMBOS_TO_EXPORT = {'ALL' if math.isinf(MAX_COMBOS_TO_EXPORT) else MAX_COMBOS_TO_EXPORT}")
    log(
        f"[INFO] WORKER_MODE = {WORKER_MODE}, WORKER_GPU_ID = {WORKER_GPU_ID}, "
        f"COMBO_RANGE = {COMBO_START}-{COMBO_END if COMBO_END > 0 else 'ALL'}, "
        f"SKIP_SUMMARY_WRITE = {SKIP_SUMMARY_WRITE}, SKIP_FIGURE_EXPORT = {SKIP_FIGURE_EXPORT}"
    )
    log(f"[INFO] EXPORT_FORMATS = {','.join(sorted(EXPORT_FORMATS))}")
    log(f"[INFO] WRITE_COMBO_PAYLOADS = {int(WRITE_COMBO_PAYLOADS)}")

    combo_key_cols = [
        "Analysis_Tier",
        "Domain",
        "Subset",
        "Mode",
        "Diversity_Order",
        "Q_Label",
        "Estimate_Definition",
    ]
    if SKIP_SUMMARY_WRITE:
        combo_manifest = must_read_csv(DIR_CSV / "combo_manifest.csv")
        cells_df = must_read_csv(DIR_CSV / "comparison_cells_filtered.csv")
        if ONLY_FIGURE_06:
            export_heatmap_only_from_cells(cells_df)
            return
        if ONLY_FIGURE_06_2:
            export_mask_only_from_cells(cells_df)
            return
        elif LIGHTWEIGHT_FIGURE_FILTER:
            reference_bridge_summary = pd.DataFrame()
            sweep_sc_summary = pd.DataFrame()
            group_qc_summary = pd.DataFrame()
        else:
            reference_bridge_summary = must_read_csv(DIR_CSV / "reference_bridge_summary.csv")
            sweep_sc_summary = must_read_csv(DIR_CSV / "sweep_sc_summary.csv")
            group_qc_summary = must_read_csv(DIR_CSV / "group_qc_summary.csv")
    else:
        cells_df = pd.read_csv(
            CELLS_CSV,
            usecols=[
                "Domain",
                "Subset",
                "Mode",
                "Diversity_Order",
                "Q_Label",
                "Estimate_Definition",
                "Cutoff_Key",
                "Cutoff_Side",
                "Retain_Pct",
                "Retain_Ratio",
                "Cutoff_Signed_Pct",
                "Cutoff_Label",
                "Reference_Target_SC",
                "Sweep_SC",
                "Sweep_SC_Label",
                "Sweep_SC_Pct",
                "Sweep_Is_Asymptotic",
                "Max_Sweep_M_over_N",
                "Sweep_Status",
                "Primary_Eligible",
                "Metric",
                "Metric_Column",
                "Value",
            ],
        )
        cells_df["Diversity_Order"] = cells_df["Diversity_Order"].astype(int)
        cells_df["Retain_Pct"] = cells_df["Retain_Pct"].astype(int)
        cells_df["Retain_Ratio"] = cells_df["Retain_Ratio"].astype(float)
        cells_df["Cutoff_Signed_Pct"] = cells_df["Cutoff_Signed_Pct"].astype(float)
        cells_df["Reference_Target_SC"] = cells_df["Reference_Target_SC"].astype(float)
        cells_df["Sweep_SC"] = cells_df["Sweep_SC"].astype(float)
        cells_df["Sweep_SC_Pct"] = cells_df["Sweep_SC_Pct"].astype(int)
        cells_df["Sweep_Is_Asymptotic"] = cells_df["Sweep_Is_Asymptotic"].astype(bool)
        cells_df["Max_Sweep_M_over_N"] = cells_df["Max_Sweep_M_over_N"].astype(float)
        cells_df["Primary_Eligible"] = cells_df["Primary_Eligible"].astype("boolean")
        cells_df["Value"] = cells_df["Value"].astype(float)
        cells_df["Metric_Pretty"] = cells_df["Metric"].map(metric_pretty)
        cells_df["Estimate_Pretty"] = cells_df["Estimate_Definition"].map(estimate_pretty)
        cells_df["Analysis_Tier"] = cells_df["Q_Label"].map(analysis_tier)
        cells_df["Q0_Caution"] = [q0_caution_text(q, est) for q, est in zip(cells_df["Q_Label"], cells_df["Estimate_Definition"])]
        cells_df["Tail_Display"] = cells_df["Cutoff_Side"].map(tail_display)

        for col, vals in [
            ("Domain", DOMAIN_FILTER),
            ("Subset", SUBSET_FILTER),
            ("Mode", MODE_FILTER),
            ("Q_Label", Q_FILTER),
            ("Estimate_Definition", ESTIMATE_FILTER),
            ("Metric", METRIC_FILTER),
        ]:
            cells_df = maybe_filter_in(cells_df, col, vals)
        if cells_df.empty:
            raise RuntimeError("No rows remain in cell long data after filters.")

        sweep_levels = (
            cells_df[["Sweep_SC", "Sweep_SC_Label"]]
            .drop_duplicates()
            .sort_values("Sweep_SC")["Sweep_SC_Label"]
            .tolist()
        )
        cells_df["Sweep_SC_F"] = pd.Categorical(cells_df["Sweep_SC_Label"], categories=sweep_levels, ordered=True)

        ref_rows = (
            cells_df[
                [
                    "Analysis_Tier",
                    "Domain",
                    "Subset",
                    "Mode",
                    "Diversity_Order",
                    "Q_Label",
                    "Estimate_Definition",
                    "Reference_Target_SC",
                    "Sweep_SC",
                ]
            ]
            .drop_duplicates()
            .copy()
        )
        reference_lookup_rows: List[Dict[str, object]] = []
        ref_group_cols = [
            "Analysis_Tier",
            "Domain",
            "Subset",
            "Mode",
            "Diversity_Order",
            "Q_Label",
            "Estimate_Definition",
            "Reference_Target_SC",
        ]
        for keys, grp in ref_rows.groupby(ref_group_cols, dropna=False, sort=False):
            row = dict(zip(ref_group_cols, keys if isinstance(keys, tuple) else (keys,)))
            selected = pick_reference_sweep(grp["Sweep_SC"], float(grp["Reference_Target_SC"].iloc[0]), REFERENCE_SWEEP_OVERRIDE)
            row["Reference_Sweep_Selected"] = selected
            row["Reference_Sweep_Label"] = format_sweep_label(selected)
            reference_lookup_rows.append(row)
        reference_lookup = pd.DataFrame(reference_lookup_rows)

        group_reference_lookup = reference_lookup[
            ["Domain", "Subset", "Mode", "Reference_Target_SC", "Reference_Sweep_Selected", "Reference_Sweep_Label"]
        ].drop_duplicates()

        def join_reference_lookup(df: pd.DataFrame) -> pd.DataFrame:
            out = df.merge(
                reference_lookup,
                on=[
                    "Analysis_Tier",
                    "Domain",
                    "Subset",
                    "Mode",
                    "Diversity_Order",
                    "Q_Label",
                    "Estimate_Definition",
                    "Reference_Target_SC",
                ],
                how="left",
            )
            out["Is_Reference_Sweep_Selected"] = np.abs(out["Sweep_SC"] - out["Reference_Sweep_Selected"]) < PRIMARY_SC_TOL
            return out

        cells_df = join_reference_lookup(cells_df)
        if ONLY_FIGURE_06:
            export_heatmap_only_from_cells(cells_df)
            return
        if ONLY_FIGURE_06_2:
            export_mask_only_from_cells(cells_df)
            return
        selected_metric_keys = sorted(cells_df["Metric"].dropna().unique().tolist())
        selected_outer_value_cols = sorted({metric_value_col(metric_key) for metric_key in selected_metric_keys if metric_value_col(metric_key) is not None})

        outer_cols = [
            "OuterRep",
            "Cutoff_Key",
            "Cutoff_Side",
            "Retain_Pct",
            "Retain_Ratio",
            "Cutoff_Signed_Pct",
            "Cutoff_Label",
            "Reference_Target_SC",
            "Sweep_SC",
            "Sweep_SC_Label",
            "Sweep_SC_Pct",
            "Sweep_Is_Asymptotic",
            "Domain",
            "Subset",
            "Mode",
            "Diversity_Order",
            "Q_Label",
            "Estimate_Definition",
            "Max_Sweep_M_over_N",
            "Sweep_Status",
            "Primary_Eligible",
        ] + selected_outer_value_cols
        outer_wide = pd.read_csv(OUTER_CSV, usecols=outer_cols)
        outer_wide["OuterRep"] = outer_wide["OuterRep"].astype(int)
        outer_wide["Diversity_Order"] = outer_wide["Diversity_Order"].astype(int)
        outer_wide["Retain_Pct"] = outer_wide["Retain_Pct"].astype(int)
        outer_wide["Retain_Ratio"] = outer_wide["Retain_Ratio"].astype(float)
        outer_wide["Cutoff_Signed_Pct"] = outer_wide["Cutoff_Signed_Pct"].astype(float)
        outer_wide["Reference_Target_SC"] = outer_wide["Reference_Target_SC"].astype(float)
        outer_wide["Sweep_SC"] = outer_wide["Sweep_SC"].astype(float)
        outer_wide["Sweep_SC_Pct"] = outer_wide["Sweep_SC_Pct"].astype(int)
        outer_wide["Sweep_Is_Asymptotic"] = outer_wide["Sweep_Is_Asymptotic"].astype(bool)
        outer_wide["Max_Sweep_M_over_N"] = outer_wide["Max_Sweep_M_over_N"].astype(float)
        outer_wide["Primary_Eligible"] = outer_wide["Primary_Eligible"].astype("boolean")
        for col in selected_outer_value_cols:
            outer_wide[col] = pd.to_numeric(outer_wide[col], errors="coerce")
        outer_wide["Estimate_Pretty"] = outer_wide["Estimate_Definition"].map(estimate_pretty)
        outer_wide["Analysis_Tier"] = outer_wide["Q_Label"].map(analysis_tier)
        outer_wide["Q0_Caution"] = [q0_caution_text(q, est) for q, est in zip(outer_wide["Q_Label"], outer_wide["Estimate_Definition"])]
        outer_wide["Tail_Display"] = outer_wide["Cutoff_Side"].map(tail_display)
        outer_wide["Sweep_SC_F"] = pd.Categorical(outer_wide["Sweep_SC_Label"], categories=sweep_levels, ordered=True)
        outer_wide = join_reference_lookup(outer_wide)

        for col, vals in [
            ("Domain", DOMAIN_FILTER),
            ("Subset", SUBSET_FILTER),
            ("Mode", MODE_FILTER),
            ("Q_Label", Q_FILTER),
            ("Estimate_Definition", ESTIMATE_FILTER),
        ]:
            outer_wide = maybe_filter_in(outer_wide, col, vals)
        if outer_wide.empty:
            raise RuntimeError("No rows remain in outer result data after filters.")
        if LIGHTWEIGHT_FIGURE_FILTER:
            export_reference_and_heatmap_from_cells_outer(cells_df, outer_wide)
            return

        group_df = pd.read_csv(
            GROUP_CSV,
            usecols=[
                "OuterRep",
                "Cutoff_Key",
                "Cutoff_Side",
                "Retain_Pct",
                "Retain_Ratio",
                "Cutoff_Signed_Pct",
                "Cutoff_Label",
                "Reference_Target_SC",
                "Sweep_SC",
                "Sweep_SC_Label",
                "Sweep_SC_Pct",
                "Sweep_Is_Asymptotic",
                "Domain",
                "Subset",
                "Mode",
                "Group",
                "Reference_N",
                "Sweep_Planned_M",
                "Sweep_M_over_N",
                "Uses_Extrapolation",
                "Sweep_Status",
                "Primary_Eligible",
                "Realized_SC_Model",
                "Realized_SC_Empirical",
                "Mean_f1",
                "Mean_Realized_N",
            ],
        )
        group_df["OuterRep"] = group_df["OuterRep"].astype(int)
        group_df["Retain_Pct"] = group_df["Retain_Pct"].astype(int)
        group_df["Retain_Ratio"] = group_df["Retain_Ratio"].astype(float)
        group_df["Cutoff_Signed_Pct"] = group_df["Cutoff_Signed_Pct"].astype(float)
        group_df["Reference_Target_SC"] = group_df["Reference_Target_SC"].astype(float)
        group_df["Sweep_SC"] = group_df["Sweep_SC"].astype(float)
        group_df["Sweep_SC_Pct"] = group_df["Sweep_SC_Pct"].astype(int)
        group_df["Sweep_Is_Asymptotic"] = group_df["Sweep_Is_Asymptotic"].astype(bool)
        group_df["Primary_Eligible"] = group_df["Primary_Eligible"].astype("boolean")
        group_df["Reference_N"] = group_df["Reference_N"].astype(float)
        group_df["Sweep_Planned_M"] = group_df["Sweep_Planned_M"].astype(float)
        group_df["Sweep_M_over_N"] = group_df["Sweep_M_over_N"].astype(float)
        group_df["Uses_Extrapolation"] = group_df["Uses_Extrapolation"].astype("boolean")
        group_df["Realized_SC_Model"] = group_df["Realized_SC_Model"].astype(float)
        group_df["Realized_SC_Empirical"] = group_df["Realized_SC_Empirical"].astype(float)
        group_df["Mean_f1"] = group_df["Mean_f1"].astype(float)
        group_df["Mean_Realized_N"] = group_df["Mean_Realized_N"].astype(float)
        group_df["Tail_Display"] = group_df["Cutoff_Side"].map(tail_display)
        group_df["Sweep_SC_F"] = pd.Categorical(group_df["Sweep_SC_Label"], categories=sweep_levels, ordered=True)
        group_df["Analysis_Tier"] = "group_support"
        group_df["Diversity_Order"] = pd.Series([pd.NA] * len(group_df), dtype="Int64")
        group_df["Q_Label"] = "group_support"
        group_df["Estimate_Definition"] = "group_support"
        group_df = group_df.merge(
            group_reference_lookup,
            on=["Domain", "Subset", "Mode", "Reference_Target_SC"],
            how="left",
        )
        group_df["Is_Reference_Sweep_Selected"] = np.abs(group_df["Sweep_SC"] - group_df["Reference_Sweep_Selected"]) < PRIMARY_SC_TOL
        for col, vals in [("Domain", DOMAIN_FILTER), ("Subset", SUBSET_FILTER), ("Mode", MODE_FILTER)]:
            group_df = maybe_filter_in(group_df, col, vals)

        outer_export_long = pd.concat(
            [
                materialize_outer_metric(
                    outer_wide.merge(
                        cells_df[cells_df["Metric"] == metric_key][combo_key_cols].drop_duplicates(),
                        on=combo_key_cols,
                        how="inner",
                    ),
                    metric_key,
                )
                for metric_key in selected_metric_keys
            ],
            ignore_index=True,
        )

        reference_bridge_parts: List[pd.DataFrame] = []
        for metric_key in selected_metric_keys:
            key_df = cells_df[cells_df["Metric"] == metric_key][combo_key_cols].drop_duplicates()
            metric_wide = outer_wide.merge(key_df, on=combo_key_cols, how="inner")
            metric_ref = metric_wide[
                (np.abs(metric_wide["Reference_Target_SC"] - PRIMARY_REFERENCE_TARGET_SC) < PRIMARY_SC_TOL)
                & metric_wide["Is_Reference_Sweep_Selected"].astype(bool)
            ].copy()
            reference_bridge_parts.append(
                summarise_outer_metric(
                    metric_ref,
                    value_col=metric_value_col(metric_key),
                    metric_key=metric_key,
                )
            )
        reference_bridge_summary = pd.concat(reference_bridge_parts, ignore_index=True)

        sweep_sc_summary = pd.concat(
            [summarise_outer_metric(outer_wide, metric_value_col(metric_key), metric_key) for metric_key in selected_metric_keys],
            ignore_index=True,
        )
        group_qc_summary = summarise_group_qc(group_df)

        combo_manifest = (
            cells_df[
                [
                    "Analysis_Tier",
                    "Domain",
                    "Subset",
                    "Mode",
                    "Diversity_Order",
                    "Q_Label",
                    "Estimate_Definition",
                    "Estimate_Pretty",
                    "Q0_Caution",
                    "Metric",
                    "Metric_Pretty",
                    "Reference_Target_SC",
                    "Reference_Sweep_Selected",
                    "Reference_Sweep_Label",
                ]
            ]
            .drop_duplicates()
            .sort_values(["Analysis_Tier", "Domain", "Subset", "Mode", "Diversity_Order", "Metric"])
            .reset_index(drop=True)
        )
        combo_manifest["combo_id"] = np.arange(1, len(combo_manifest) + 1)
        combo_manifest["output_rel_dir"] = combo_manifest.apply(
            lambda row: str(
                Path("figures")
                / clean_id(row["Analysis_Tier"])
                / clean_id(row["Domain"])
                / clean_id(row["Subset"])
                / clean_id(row["Mode"])
                / clean_id(row["Q_Label"])
                / clean_id(row["Estimate_Definition"])
                / clean_id(row["Metric"])
            ),
            axis=1,
        )

    combo_manifest["combo_id"] = pd.to_numeric(combo_manifest["combo_id"], errors="raise").astype(int)
    if COMBO_START > 0:
        combo_end = COMBO_END if COMBO_END > 0 else int(combo_manifest["combo_id"].max())
        combo_manifest["export_selected"] = (
            (combo_manifest["combo_id"] >= COMBO_START)
            & (combo_manifest["combo_id"] <= combo_end)
        )
    else:
        combo_manifest["export_selected"] = (
            True if math.isinf(MAX_COMBOS_TO_EXPORT) else combo_manifest["combo_id"] <= MAX_COMBOS_TO_EXPORT
        )

    if not SKIP_SUMMARY_WRITE:
        combo_manifest.to_csv(DIR_CSV / "combo_manifest.csv", index=False)
        cells_df.to_csv(DIR_CSV / "comparison_cells_filtered.csv", index=False)
        reference_bridge_summary.to_csv(DIR_CSV / "reference_bridge_summary.csv", index=False)
        sweep_sc_summary.to_csv(DIR_CSV / "sweep_sc_summary.csv", index=False)
        group_qc_summary.to_csv(DIR_CSV / "group_qc_summary.csv", index=False)
        q0_note = [
            "q = 0 caution for v36.9 recommended outputs",
            "- q0 / Chao1_asymptotic is retained, but should be interpreted as a cautious lower-bound-oriented reference.",
            "- q0 / Observed_Standardized is empirical and is not an asymptotic endpoint.",
            "- Main recommended figures remain q1 and q2; q0 is exported alongside them with this note.",
        ]
        (DIR_CSV / "q0_caution_note.txt").write_text("\n".join(q0_note) + "\n", encoding="utf-8")

        if WRITE_COMBO_PAYLOADS:
            metric_combo_key_cols = combo_key_cols + ["Metric"]
            group_combo_key_cols = ["Domain", "Subset", "Mode"]
            for _, one in combo_manifest[combo_manifest["export_selected"]].iterrows():
                payload_dir = ensure_dir(combo_payload_dir(int(one["combo_id"])))
                filter_by_combo_keys(cells_df, one, metric_combo_key_cols).to_csv(payload_dir / "cells.csv", index=False)
                filter_by_combo_keys(outer_export_long, one, metric_combo_key_cols).to_csv(payload_dir / "outer.csv", index=False)
                filter_by_combo_keys(reference_bridge_summary, one, metric_combo_key_cols).to_csv(payload_dir / "reference.csv", index=False)
                filter_by_combo_keys(sweep_sc_summary, one, metric_combo_key_cols).to_csv(payload_dir / "sweep.csv", index=False)
                filter_by_combo_keys(group_qc_summary, one, group_combo_key_cols).to_csv(payload_dir / "group.csv", index=False)

    combos_to_export = combo_manifest[combo_manifest["export_selected"]].copy()
    if SKIP_FIGURE_EXPORT:
        combos_to_export = combos_to_export.iloc[0:0].copy()
    if sns is not None:
        sns.set_theme(style="whitegrid")
    if should_export_figure("01") and not combos_to_export.empty:
        render_fig01_legend_only(DIR_FIG / "figure01_boxplot_legend_style2_horizontal_single_row.png")

    metric_combo_key_cols = combo_key_cols + ["Metric"]
    group_combo_key_cols = ["Domain", "Subset", "Mode"]

    for _, one in combos_to_export.iterrows():
        if SKIP_SUMMARY_WRITE:
            payload_dir = combo_payload_dir(int(one["combo_id"]))
            cell_one = must_read_csv(payload_dir / "cells.csv")
            outer_one = must_read_csv(payload_dir / "outer.csv")
            reference_one = (
                must_read_csv(payload_dir / "reference.csv")
                if should_export_figure("02")
                else pd.DataFrame()
            )
            sweep_summary_one = (
                must_read_csv(payload_dir / "sweep.csv")
                if should_export_figure("05")
                else pd.DataFrame()
            )
            group_one = (
                must_read_csv(payload_dir / "group.csv")
                if should_export_figure("07") or should_export_figure("08")
                else pd.DataFrame()
            )
        else:
            cell_one = filter_by_combo_keys(cells_df, one, metric_combo_key_cols)
            outer_one = filter_by_combo_keys(outer_export_long, one, metric_combo_key_cols)
            reference_one = filter_by_combo_keys(reference_bridge_summary, one, metric_combo_key_cols)
            sweep_summary_one = filter_by_combo_keys(sweep_sc_summary, one, metric_combo_key_cols)
            group_one = filter_by_combo_keys(group_qc_summary, one, group_combo_key_cols)

        outdir_one = ensure_dir(
            DIR_FIG
            / clean_id(one["Analysis_Tier"])
            / clean_id(one["Domain"])
            / clean_id(one["Subset"])
            / clean_id(one["Mode"])
            / clean_id(one["Q_Label"])
            / clean_id(one["Estimate_Definition"])
            / clean_id(one["Metric"])
        )

        if should_export_figure("01"):
            p1 = plot_reference_bridge_boxplots(outer_one, include_legend=True)
            save_plot_pdf_png(
                p1,
                outdir_one / "01_reference_bridge_boxplots_paired.pdf",
                outdir_one / "01_reference_bridge_boxplots_paired.png",
                20.4,
                6.6,
            )
            p1_nolegend = plot_reference_bridge_boxplots(outer_one, include_legend=False)
            save_plot_pdf_png(
                p1_nolegend,
                outdir_one / "01_reference_bridge_boxplots_paired_nolegend.pdf",
                outdir_one / "01_reference_bridge_boxplots_paired_nolegend.png",
                17.6,
                6.6,
            )

        if should_export_figure("02"):
            p2 = plot_reference_bridge_quantile_lines(reference_one)
            save_plot_pdf_png(p2, outdir_one / "02_reference_bridge_quantile_lines.pdf", outdir_one / "02_reference_bridge_quantile_lines.png", 12, 6)

        if should_export_figure("03"):
            p3 = plot_reference_bridge_histograms(outer_one)
            save_plot_pdf_png(p3, outdir_one / "03_reference_bridge_histograms.pdf", outdir_one / "03_reference_bridge_histograms.png", 13, 10)

        if should_export_figure("04"):
            p4 = plot_representative_sweep_boxplots(outer_one)
            save_plot_pdf_png(p4, outdir_one / "04_representative_sweep_boxplots_paired.pdf", outdir_one / "04_representative_sweep_boxplots_paired.png", 18, 12)

        if should_export_figure("05"):
            p5 = plot_sweep_sc_support_lines(sweep_summary_one)
            save_plot_pdf_png(p5, outdir_one / "05_sweep_sc_support_lines.pdf", outdir_one / "05_sweep_sc_support_lines.png", 14, 10)

        if should_export_figure("06"):
            p6 = plot_sweep_sc_heatmap(cell_one)
            save_plot_pdf_png(p6, outdir_one / "06_sweep_sc_heatmap.pdf", outdir_one / "06_sweep_sc_heatmap.png", 15, 6.4)

        if should_export_figure("06_2"):
            p6_2 = plot_sweep_sc_large_effect_mask(cell_one)
            save_plot_pdf_png(p6_2, outdir_one / "06_2_sweep_sc_large_effect_mask.pdf", outdir_one / "06_2_sweep_sc_large_effect_mask.png", 15, 6.4)

        if should_export_figure("07"):
            p7 = plot_sweep_m_over_n_lines(group_one)
            save_plot_pdf_png(p7, outdir_one / "07_sweep_m_over_n_lines.pdf", outdir_one / "07_sweep_m_over_n_lines.png", 14, 9)

        if should_export_figure("08"):
            p8 = plot_realized_sc_qc(group_one)
            save_plot_pdf_png(p8, outdir_one / "08_realized_sc_qc.pdf", outdir_one / "08_realized_sc_qc.png", 14, 9)

    run_note = [
        f"RUN_DIR: {RUN_DIR}",
        f"RUN_DIR_SOURCE: {RUN_DIR_SOURCE}",
        f"CELLS_CSV: {CELLS_CSV}",
        f"OUTER_CSV: {OUTER_CSV}",
        f"GROUP_CSV: {GROUP_CSV}",
        f"OUT_DIR: {OUT_DIR}",
        f"ROOT_DIR: {ROOT_DIR}",
        f"DEVICE: {DEVICE}",
        f"FIGURE_FILTER: {','.join(sorted(FIGURE_FILTER)) if FIGURE_FILTER else 'ALL'}",
        f"EXPORT_FORMATS: {','.join(sorted(EXPORT_FORMATS))}",
        f"WRITE_COMBO_PAYLOADS: {int(WRITE_COMBO_PAYLOADS)}",
        "",
        "Interpretation notes:",
        "- v36.9 is explicitly two-stage: Stage 1 fixes the reference pseudo-sample at Reference_Target_SC; Stage 2 sweeps sample coverage from that reference sample.",
        f"- PRIMARY_REFERENCE_TARGET_SC = {PRIMARY_REFERENCE_TARGET_SC}",
        (
            f"- Reference-like bridge slice uses Sweep_SC nearest to the reference target (override requested: {REFERENCE_SWEEP_OVERRIDE})."
            if math.isfinite(REFERENCE_SWEEP_OVERRIDE)
            else "- Reference-like bridge slice uses Sweep_SC nearest to the reference target."
        ),
        "- Sample coverage s is the main interpretation axis. The Asymptotic row is a separate endpoint.",
        "- 06_sweep_sc_heatmap and 06_2_sweep_sc_large_effect_mask use a signed display cutoff coordinate: 0% is no cutoff, negative values cut high-intensity mass, and positive values cut low-intensity mass.",
        "- For Cliff's delta, 06_sweep_sc_heatmap splits positive and negative surfaces into separate panels.",
        f"- 06_2_sweep_sc_large_effect_mask marks medium/large Cliff's delta effects at |delta| >= {CLIFF_MEDIUM_THRESHOLD} / {CLIFF_LARGE_THRESHOLD}.",
        "- NegLogP and WelchNegLogP panels receive p=0.05 and p=0.01 reference lines.",
        "- A panels receive a no-difference guide at 0.5.",
        "- A = P(Syneco > Conv) on the probability scale; Cliff's delta = 2A - 1.",
        "- q = 0 outputs are retained, but q0_caution_note.txt should be read alongside them.",
        (
            f"- Exported combos: {len(combos_to_export)} / {len(combo_manifest)}"
            if math.isinf(MAX_COMBOS_TO_EXPORT)
            else f"- Exported combos: {len(combos_to_export)} / {len(combo_manifest)} (limited by TOMATO_MAX_COMBOS)"
        ),
    ]
    if not SKIP_SUMMARY_WRITE:
        (OUT_DIR / "run_note.txt").write_text("\n".join(run_note) + "\n", encoding="utf-8")
    log(f"[SUCCESS] Output written to: {OUT_DIR}")


def resolve_gpu_ids() -> List[int]:
    if os.getenv("TOMATO_FORCE_CPU", "0") == "1" or torch is None or not torch.cuda.is_available():
        return []
    visible = torch.cuda.device_count()
    parsed: List[int] = []
    if GPU_ID_RAW is not None:
        for item in GPU_ID_RAW:
            try:
                gpu_id = int(item)
            except (TypeError, ValueError):
                continue
            if 0 <= gpu_id < visible and gpu_id not in parsed:
                parsed.append(gpu_id)
    if not parsed:
        parsed = list(range(visible))
    return parsed[:REQUESTED_NUM_WORKERS]


def split_contiguous_ranges(total: int, num_parts: int) -> List[tuple[int, int]]:
    num_parts = max(1, min(num_parts, total))
    base = total // num_parts
    rem = total % num_parts
    out: List[tuple[int, int]] = []
    start = 1
    for idx in range(num_parts):
        size = base + (1 if idx < rem else 0)
        end = start + size - 1
        out.append((start, end))
        start = end + 1
    return out


def launch_worker(
    *,
    worker_name: str,
    gpu_id: int,
    skip_summary_write: bool,
    skip_figure_export: bool,
    combo_start: int = 0,
    combo_end: int = 0,
) -> tuple[subprocess.Popen, object, Path]:
    log_dir = OUT_DIR / "_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{worker_name}.log"
    env = os.environ.copy()
    env.update(
        {
            "TOMATO_4GPU_WORKER": "1",
            "TOMATO_4GPU_GPU_ID": str(gpu_id),
            "TOMATO_4GPU_SKIP_SUMMARY_WRITE": "1" if skip_summary_write else "0",
            "TOMATO_4GPU_SKIP_FIGURE_EXPORT": "1" if skip_figure_export else "0",
            "TOMATO_4GPU_COMBO_START": str(combo_start),
            "TOMATO_4GPU_COMBO_END": str(combo_end),
            "TOMATO_OUT_DIR": str(OUT_DIR),
        }
    )
    handle = log_path.open("w", encoding="utf-8")
    proc = subprocess.Popen(
        [sys.executable, str(Path(__file__).resolve())],
        cwd=str(ROOT_DIR),
        env=env,
        stdout=handle,
        stderr=subprocess.STDOUT,
    )
    log(
        f"[ORCH] worker={worker_name} gpu={gpu_id} "
        f"combo_range={combo_start}-{combo_end if combo_end > 0 else 'ALL'} log={log_path}"
    )
    return proc, handle, log_path


def overwrite_run_note(gpu_ids: Sequence[int], combo_ranges: Sequence[tuple[int, int]]) -> None:
    run_note = [
        f"RUN_DIR: {RUN_DIR}",
        f"RUN_DIR_SOURCE: {RUN_DIR_SOURCE}",
        f"CELLS_CSV: {CELLS_CSV}",
        f"OUTER_CSV: {OUTER_CSV}",
        f"GROUP_CSV: {GROUP_CSV}",
        f"OUT_DIR: {OUT_DIR}",
        f"ROOT_DIR: {ROOT_DIR}",
        "DEVICE: multi_gpu_orchestrated",
        f"GPU_IDS: {', '.join(str(v) for v in gpu_ids)}",
        f"COMBO_RANGES: {', '.join(f'{start}-{end}' for start, end in combo_ranges)}",
        "",
        "Interpretation notes:",
        "- v36.9 downstream 4gpu mode computes summaries once, then shards figure export across workers.",
        f"- PRIMARY_REFERENCE_TARGET_SC = {PRIMARY_REFERENCE_TARGET_SC}",
        (
            f"- Reference-like bridge slice uses Sweep_SC nearest to the reference target (override requested: {REFERENCE_SWEEP_OVERRIDE})."
            if math.isfinite(REFERENCE_SWEEP_OVERRIDE)
            else "- Reference-like bridge slice uses Sweep_SC nearest to the reference target."
        ),
        "- Sample coverage s is the main interpretation axis. The Asymptotic row is a separate endpoint.",
        "- 06_sweep_sc_heatmap and 06_2_sweep_sc_large_effect_mask use a signed display cutoff coordinate: 0% is no cutoff, negative values cut high-intensity mass, and positive values cut low-intensity mass.",
        "- NegLogP and WelchNegLogP panels receive p=0.05 and p=0.01 reference lines.",
        "- A panels receive a no-difference guide at 0.5.",
        "- A = P(Syneco > Conv) on the probability scale; Cliff's delta = 2A - 1.",
        "- q = 0 outputs are retained, but q0_caution_note.txt should be read alongside them.",
    ]
    (OUT_DIR / "run_note.txt").write_text("\n".join(run_note) + "\n", encoding="utf-8")


def orchestrate_4gpu() -> None:
    gpu_ids = resolve_gpu_ids()
    if len(gpu_ids) <= 1:
        log("[ORCH] Fewer than 2 usable GPUs detected. Falling back to single-worker execution.")
        main()
        return

    summary_proc, summary_handle, summary_log = launch_worker(
        worker_name="summary",
        gpu_id=gpu_ids[0],
        skip_summary_write=False,
        skip_figure_export=True,
    )
    summary_rc = summary_proc.wait()
    summary_handle.close()
    if summary_rc != 0:
        raise RuntimeError(f"Summary prepass failed. See log: {summary_log}")

    combo_manifest_path = DIR_CSV / "combo_manifest.csv"
    if not combo_manifest_path.exists():
        raise RuntimeError(f"Summary prepass did not produce combo manifest: {combo_manifest_path}")
    combo_manifest = pd.read_csv(combo_manifest_path)
    selected_mask = combo_manifest["export_selected"].astype(str).str.lower().isin({"true", "1"})
    selected = combo_manifest[selected_mask].copy()
    if selected.empty:
        raise RuntimeError("No combos selected for export after summary prepass.")

    combo_ranges = split_contiguous_ranges(len(selected), len(gpu_ids))
    launched: List[tuple[subprocess.Popen, object, Path]] = []
    try:
        for idx, (local_start, local_end) in enumerate(combo_ranges, start=1):
            combo_start = int(selected.iloc[local_start - 1]["combo_id"])
            combo_end = int(selected.iloc[local_end - 1]["combo_id"])
            proc, handle, log_path = launch_worker(
                worker_name=f"figures_{idx:02d}",
                gpu_id=gpu_ids[idx - 1],
                skip_summary_write=True,
                skip_figure_export=False,
                combo_start=combo_start,
                combo_end=combo_end,
            )
            launched.append((proc, handle, log_path))

        failed_logs: List[Path] = []
        for proc, handle, log_path in launched:
            rc = proc.wait()
            handle.close()
            if rc != 0:
                failed_logs.append(log_path)
        if failed_logs:
            joined = "\n".join(str(path) for path in failed_logs)
            raise RuntimeError(f"One or more figure workers failed. See logs:\n{joined}")

        overwrite_run_note(gpu_ids, [(int(selected.iloc[start - 1]['combo_id']), int(selected.iloc[end - 1]['combo_id'])) for start, end in combo_ranges])
        log(f"[SUCCESS] 4GPU downstream output written to: {OUT_DIR}")
    finally:
        for proc, handle, _ in launched:
            try:
                if proc.poll() is None:
                    proc.terminate()
            except Exception:
                pass
            try:
                handle.close()
            except Exception:
                pass


if __name__ == "__main__":
    if WORKER_MODE == 1:
        main()
    else:
        orchestrate_4gpu()
