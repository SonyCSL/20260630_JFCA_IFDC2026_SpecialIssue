#!/usr/bin/env python3
# ==============================================================================
# v36.9 Python / PyTorch port with 4GPU orchestration v4 Stage2 iNEXT q2-fast
# 日本語メモ:
#   2026-06-21 の v3 では、q=1/q=2 の Stage 2 diversity を
#   empirical Hill from bootstrap samples ではなく、Chao et al. Table 2 /
#   iNEXT TD.m.est 型の abundance rarefaction / extrapolation estimator で
#   計算する。bootstrap は Figure 01 / CI / stability diagnostics 用の
#   分布生成層として残し、各 bootstrap reference count vector にも
#   同じ estimator を適用する。
#
#   2026-04-08 の v2 では、Brunner-Munzel を主指標のまま維持しつつ、
#   Welch の p-value と -log10(p) を optional に追加した。
#   これは主結果を Welch へ切り替えるためではなく、BM と大きくは
#   変わらないことを確認する補助線として使うためである。
#   - Signed 199-cutoff design:
#       low 1..99%, none, high 1..99%
#   - Stage 1:
#       build one reference pseudo-sample at fixed target sample coverage = 0.965
#       using the model-based SC(m) on the intensity-derived probability vector
#   - Stage 2:
#       from that reference sample, evaluate coverage-based rarefaction /
#       extrapolation on a 1%..100% sample-coverage grid
#       using a manual Chao/iNEXT-type abundance coverage estimator and
#       Table 2 / TD.m.est-type q=1/q=2 diversity estimators
#       (100% is treated as a separate Chao/iNEXT asymptotic endpoint)
#   - Effect metrics:
#       BM p-value -> -log10(BM p)
#       optional Welch p-value -> -log10(Welch p)
#       A = P(Conv < Syneco) + 0.5 P(Conv = Syneco)
#       Cliff's delta = 2A - 1
#   - Performance notes:
#       remove SciPy from the hot path
#       precompute signed-cutoff probability tensors
#       batch stage-1 tag sampling and tag-to-domain aggregation
#       cache stage-2 coverage solves and repeated sweep-size sampling
#   - POC mode:
#       Formula x Combo6yr x All only
# ==============================================================================

from __future__ import annotations

import json
import math
import os
import random
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import torch


DTYPE = torch.float32
SOLVE_DTYPE = torch.float64
BASE_SEED = 12345
P_FLOOR = 1e-300
ALPHA_P = 0.05
STATUS_RANK = [
    "reference_empty_after_mapping",
    "no_valid_boot_replicates",
    "insufficient_group_replicates",
    "bm_undefined_zero_denominator",
    "bm_invalid_df",
    "welch_undefined_zero_denominator",
    "welch_invalid_df",
    "welch_disabled",
    "ok",
]


def parse_int_env(name: str, default: int) -> int:
    raw = os.getenv(name, "")
    try:
        return int(raw)
    except (TypeError, ValueError):
        return default


BASE_SEED = parse_int_env("TOMATO_BASE_SEED", BASE_SEED)


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


def parse_float_list_env(name: str, fallback: Optional[Sequence[float]] = None) -> Optional[List[float]]:
    raw = os.getenv(name, "").strip()
    if not raw:
        return list(fallback) if fallback is not None else None
    out: List[float] = []
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        try:
            val = float(item)
        except ValueError:
            continue
        if math.isfinite(val):
            out.append(val)
    return out if out else (list(fallback) if fallback is not None else None)


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
    "Early3": ["Period1_2015_2017"],
    "early3": ["Period1_2015_2017"],
    "Last3": ["Period2_2018_2020"],
    "last3": ["Period2_2018_2020"],
    "Late3": ["Period2_2018_2020"],
    "late3": ["Period2_2018_2020"],
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


def normalize_sc_grid(values: Sequence[float]) -> List[float]:
    out: List[float] = []
    for value in values:
        value = float(value)
        if not math.isfinite(value):
            continue
        if value > 1.0 + 1e-9:
            value = value / 100.0
        value = min(max(value, 0.01), 1.0)
        out.append(round(value, 6))
    if not out:
        return [1.0]
    out = sorted(set(out))
    if not any(abs(v - 1.0) < 1e-9 for v in out):
        out.append(1.0)
    return out


DEFAULT_ROOT_DIR = Path("/home/nonaka/work/nonaka/Chao1_Intensity")
LOCAL_FALLBACK_ROOT = Path(__file__).resolve().parent
ROOT_DIR = Path(os.getenv("TOMATO_ROOT_DIR", str(DEFAULT_ROOT_DIR)))
if not ROOT_DIR.exists():
    ROOT_DIR = LOCAL_FALLBACK_ROOT
DATA_DIR = Path(os.getenv("TOMATO_DATA_DIR", str(ROOT_DIR / "data")))
if not DATA_DIR.exists() and (ROOT_DIR / "data").exists():
    DATA_DIR = ROOT_DIR / "data"

POC_MODE = parse_int_env("TOMATO_POC_MODE", 0)
REFERENCE_TARGET_SC = parse_float_env("TOMATO_REFERENCE_TARGET_SC", 0.965)
DEFAULT_SWEEP_SC_LEVELS = [round(v, 2) for v in np.arange(0.01, 1.0, 0.01)] + [1.0]
SWEEP_SC_LEVELS = normalize_sc_grid(parse_float_list_env("TOMATO_SWEEP_SC_GRID", DEFAULT_SWEEP_SC_LEVELS) or DEFAULT_SWEEP_SC_LEVELS)
N_OUTER_REP = parse_int_env("TOMATO_OUTER_REP", 6 if POC_MODE == 1 else 20)
N_BOOT_MATRIX = max(1, parse_int_env("TOMATO_BOOT_MATRIX", 12 if POC_MODE == 1 else 30))
N_UPPER_CAP = parse_int_env("TOMATO_UPPER_CAP", int(1e7))
PRIMARY_MAX_M_OVER_N = parse_float_env("TOMATO_PRIMARY_MAX_M_OVER_N", 2.0)
DOMAIN_FILTER_RAW = parse_env_list("TOMATO_DOMAIN_FILTER")
DOMAIN_FILTER = expand_filter_aliases(DOMAIN_FILTER_RAW, DOMAIN_FILTER_ALIASES)
SUBSET_FILTER = parse_env_list("TOMATO_SUBSET_FILTER")
MODE_FILTER_RAW = parse_env_list("TOMATO_MODE_FILTER")
MODE_FILTER = expand_filter_aliases(MODE_FILTER_RAW, MODE_FILTER_ALIASES)
EXCLUDE_FORMULAS = parse_env_list("TOMATO_EXCLUDE_FORMULAS", [])
EXCLUDE_FORMULA_MODE = os.getenv("TOMATO_EXCLUDE_FORMULA_MODE", "zero").strip().lower()
INTENSITY_TRANSFORM = os.getenv("TOMATO_INTENSITY_TRANSFORM", "raw").strip().lower()
STAGE2_EXTRAPOLATION_METHOD = os.getenv("TOMATO_STAGE2_EXTRAPOLATION_METHOD", "inext_table2_bootstrap").strip().lower()
INEXT_UNSEEN_CAP = max(0, parse_int_env("TOMATO_INEXT_UNSEEN_CAP", 0))
INEXT_UNSEEN_MIN_CATEGORIES = max(1, parse_int_env("TOMATO_INEXT_UNSEEN_MIN_CATEGORIES", 1))
ENABLE_WELCH = parse_int_env("TOMATO_ENABLE_WELCH", 0) == 1
WRITE_CELLS_ONLY = parse_int_env("TOMATO_WRITE_CELLS_ONLY", 0) == 1
WORKER_MODE = parse_int_env("TOMATO_4GPU_WORKER", 0)
WORKER_GPU_ID = parse_int_env("TOMATO_4GPU_GPU_ID", 0)
REQUESTED_NUM_WORKERS = max(1, parse_int_env("TOMATO_4GPU_NUM_WORKERS", 4))
WORKER_CUTOFF_START = parse_int_env("TOMATO_4GPU_CUTOFF_START", 1)
WORKER_CUTOFF_END = parse_int_env("TOMATO_4GPU_CUTOFF_END", 0)
WORKER_CUTOFF_CHUNK_SIZE = max(0, parse_int_env("TOMATO_4GPU_CUTOFF_CHUNK_SIZE", 0))
OUT_DIR_OVERRIDE = os.getenv("TOMATO_OUT_DIR", "").strip()
GPU_ID_RAW = parse_env_list("TOMATO_4GPU_GPU_IDS")

if POC_MODE == 1:
    DOMAIN_FILTER = ["Formula"]
    DOMAIN_FILTER_RAW = ["Formula"]
    SUBSET_FILTER = ["All"]
    MODE_FILTER = ["Combo6yr"]
    MODE_FILTER_RAW = ["Combo6yr"]

PERIOD_MODE_YEARS = {
    "Period1_2015_2017": {"2015", "2016", "2017"},
    "Period2_2018_2020": {"2018", "2019", "2020"},
}

DEVICE = torch.device(
    f"cuda:{WORKER_GPU_ID}"
    if torch.cuda.is_available() and os.getenv("TOMATO_FORCE_CPU", "0") != "1"
    else "cpu"
)

METRIC_META = [
    {"Metric_Key": "q0_observed_standardized", "Diversity_Order": 0, "Q_Label": "q0", "Estimate_Definition": "Observed_Standardized"},
    {"Metric_Key": "q0_chao1_asymptotic", "Diversity_Order": 0, "Q_Label": "q0", "Estimate_Definition": "Chao1_asymptotic"},
    {"Metric_Key": "q1_inext_td_m_est", "Diversity_Order": 1, "Q_Label": "q1", "Estimate_Definition": "Stage2_iNEXT_TD_m_est"},
    {"Metric_Key": "q2_inext_td_m_est", "Diversity_Order": 2, "Q_Label": "q2", "Estimate_Definition": "Stage2_iNEXT_TD_m_est"},
]
Q_FILTER_RAW = parse_env_list("TOMATO_Q_FILTER")


def metric_matches_q_filter(meta: Dict[str, object], q_filter: Optional[List[str]]) -> bool:
    if q_filter is None:
        return True
    allowed = {str(item).strip().lower() for item in q_filter}
    q_label = str(meta["Q_Label"]).lower()
    q_order = str(meta["Diversity_Order"]).lower()
    return q_label in allowed or q_order in allowed


ACTIVE_METRIC_META = [meta for meta in METRIC_META if metric_matches_q_filter(meta, Q_FILTER_RAW)]
if not ACTIVE_METRIC_META:
    raise RuntimeError(f"No q metrics remain after TOMATO_Q_FILTER={Q_FILTER_RAW}")
ACTIVE_METRIC_KEYS = {str(meta["Metric_Key"]) for meta in ACTIVE_METRIC_META}
NEED_Q0_STAGE2 = bool({"q0_observed_standardized", "q0_chao1_asymptotic"} & ACTIVE_METRIC_KEYS)
NEED_Q1_INEXT = "q1_inext_td_m_est" in ACTIVE_METRIC_KEYS
NEED_Q2_INEXT = "q2_inext_td_m_est" in ACTIVE_METRIC_KEYS
Q2_FAST_PATH = (
    os.getenv("TOMATO_Q2_FAST_PATH", "1").strip() != "0"
    and NEED_Q2_INEXT
    and not NEED_Q1_INEXT
    and not NEED_Q0_STAGE2
)
Q1_FAST_PATH = (
    os.getenv("TOMATO_Q1_FAST_PATH", "1").strip() != "0"
    and NEED_Q1_INEXT
    and not NEED_Q2_INEXT
    and not NEED_Q0_STAGE2
)
Q0_FAST_PATH = (
    os.getenv("TOMATO_Q0_FAST_PATH", "1").strip() != "0"
    and NEED_Q0_STAGE2
    and not NEED_Q1_INEXT
    and not NEED_Q2_INEXT
)


def now_tag() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


TS_TAG = now_tag()
BASE_FOLDER_NAME = "PoC_v36_9_stage2_inext_4gpu" if POC_MODE == 1 else "Full_v36_9_stage2_inext_4gpu"
OUT_DIR = (
    Path(OUT_DIR_OVERRIDE)
    if OUT_DIR_OVERRIDE
    else ROOT_DIR / "out_sensitivity" / f"{BASE_FOLDER_NAME}_{TS_TAG}"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)


def log(msg: str) -> None:
    print(msg, flush=True)


def normalize_keys(series: Iterable[str]) -> List[str]:
    out: List[str] = []
    for item in series:
        x = str(item).strip()
        x = re.sub(r'^"+|"+$', "", x)
        x = re.sub(r"^'+|'+$", "", x)
        x = re.sub(r"\s+", " ", x)
        out.append(x)
    return out


def normalize_formula_list(values: Sequence[str]) -> List[str]:
    return normalize_keys(values)


def ensure_nonnegative_matrix(matrix: np.ndarray, transform: str) -> None:
    finite = matrix[np.isfinite(matrix)]
    if finite.size and float(np.nanmin(finite)) < 0:
        raise RuntimeError(
            f"TOMATO_INTENSITY_TRANSFORM={transform!r} requires nonnegative intensity values."
        )


def apply_intensity_transform(matrix: np.ndarray, transform: str) -> np.ndarray:
    transform = transform.lower()
    if transform in {"", "raw", "none"}:
        return matrix
    if transform == "sqrt":
        ensure_nonnegative_matrix(matrix, transform)
        return np.sqrt(matrix)
    if transform in {"fourthroot", "root4", "pow0.25"}:
        ensure_nonnegative_matrix(matrix, transform)
        return np.power(matrix, 0.25)
    if transform == "log1p":
        ensure_nonnegative_matrix(matrix, transform)
        return np.log1p(matrix)
    if transform == "winsor99":
        out = matrix.copy()
        for col_idx in range(out.shape[1]):
            col = out[:, col_idx]
            positive = col[np.isfinite(col) & (col > 0)]
            if positive.size == 0:
                continue
            cap = float(np.quantile(positive, 0.99))
            if math.isfinite(cap) and cap > 0:
                out[:, col_idx] = np.minimum(col, cap)
        return out
    raise RuntimeError(
        f"Unsupported TOMATO_INTENSITY_TRANSFORM={transform!r}. "
        "Expected raw, sqrt, fourthroot, root4, pow0.25, log1p, or winsor99."
    )


def mode_includes_year(mode: str, year: str) -> bool:
    if mode == "Combo6yr":
        return year in {str(value) for value in range(2015, 2021)}
    if mode in PERIOD_MODE_YEARS:
        return year in PERIOD_MODE_YEARS[mode]
    return mode == year


def quiet_read_csv(path: Path) -> Optional[pd.DataFrame]:
    if not path.exists():
        return None
    return pd.read_csv(path)


def read_form_matrix(path: Path) -> Optional[pd.DataFrame]:
    if not path.exists():
        return None
    df = pd.read_csv(path, index_col=0)
    df.index = normalize_keys(df.index.astype(str))
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return df


def median_or_na(values: Sequence[float]) -> float:
    arr = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.median(arr)) if arr.size else float("nan")


def quantile_or_na(values: Sequence[float], prob: float) -> float:
    arr = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.quantile(arr, prob)) if arr.size else float("nan")


def pick_first_non_na(values: Sequence[object]) -> object:
    for value in values:
        if pd.notna(value):
            return value
    return np.nan


def max_or_na(values: Sequence[float]) -> float:
    arr = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(arr.max()) if arr.size else float("nan")


def pick_status(values: Sequence[object]) -> object:
    candidates = [str(v) for v in values if pd.notna(v)]
    if not candidates:
        return np.nan
    ranked = [(STATUS_RANK.index(v) if v in STATUS_RANK else len(STATUS_RANK), v) for v in candidates]
    ranked.sort(key=lambda item: item[0])
    return ranked[0][1]


def clamp_p(p: float) -> float:
    if p is None or not math.isfinite(p):
        return float("nan")
    return min(max(float(p), P_FLOOR), 1.0)


def safe_neglog10p(p: float) -> float:
    p2 = clamp_p(p)
    if not math.isfinite(p2):
        return float("nan")
    return -math.log10(p2)


def format_sc_label(value: float) -> str:
    return "Asymptotic" if abs(value - 1.0) < 1e-9 else f"{value:.3f}"


def make_cutoff_plan() -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for pct in range(1, 100):
        rows.append(
            {
                "Cutoff_Key": f"low_{pct:02d}",
                "Cutoff_Side": "low",
                "Retain_Pct": pct,
                "Retain_Ratio": pct / 100.0,
                "Cutoff_Signed_Pct": -pct,
                "Cutoff_Label": f"L{pct:02d}",
            }
        )
    rows.append(
        {
            "Cutoff_Key": "none_100",
            "Cutoff_Side": "none",
            "Retain_Pct": 100,
            "Retain_Ratio": 1.0,
            "Cutoff_Signed_Pct": 0,
            "Cutoff_Label": "None",
        }
    )
    for pct in range(1, 100):
        rows.append(
            {
                "Cutoff_Key": f"high_{pct:02d}",
                "Cutoff_Side": "high",
                "Retain_Pct": pct,
                "Retain_Ratio": pct / 100.0,
                "Cutoff_Signed_Pct": pct,
                "Cutoff_Label": f"H{pct:02d}",
            }
        )
    return pd.DataFrame(rows)


FULL_CUTOFF_PLAN = make_cutoff_plan()
ACTIVE_CUTOFF_START = 1
ACTIVE_CUTOFF_END = len(FULL_CUTOFF_PLAN)
if WORKER_MODE == 1:
    ACTIVE_CUTOFF_START = min(max(WORKER_CUTOFF_START, 1), len(FULL_CUTOFF_PLAN))
    ACTIVE_CUTOFF_END = WORKER_CUTOFF_END if WORKER_CUTOFF_END > 0 else len(FULL_CUTOFF_PLAN)
    ACTIVE_CUTOFF_END = min(max(ACTIVE_CUTOFF_END, ACTIVE_CUTOFF_START), len(FULL_CUTOFF_PLAN))
    CUTOFF_PLAN = FULL_CUTOFF_PLAN.iloc[ACTIVE_CUTOFF_START - 1 : ACTIVE_CUTOFF_END].reset_index(drop=True).copy()
else:
    CUTOFF_PLAN = FULL_CUTOFF_PLAN.copy()
ACTIVE_CUTOFF_SLICE = slice(ACTIVE_CUTOFF_START - 1, ACTIVE_CUTOFF_END)


def build_signed_cutoff_prob_stack(vec: torch.Tensor) -> torch.Tensor:
    n = int(vec.numel())
    out_full = torch.zeros((len(FULL_CUTOFF_PLAN), n), dtype=DTYPE, device=DEVICE)
    if n == 0:
        return out_full[ACTIVE_CUTOFF_SLICE].clone()

    valid_mask = torch.isfinite(vec) & (vec > 0)
    if not bool(valid_mask.any().item()):
        return out_full[ACTIVE_CUTOFF_SLICE].clone()

    idx = torch.nonzero(valid_mask, as_tuple=False).squeeze(1)
    vals = vec[idx].to(dtype=DTYPE, device=DEVICE)
    total = vals.sum()
    if not bool((total > 0).item()):
        return out_full[ACTIVE_CUTOFF_SLICE].clone()

    retain_ratios = torch.arange(1, 100, device=DEVICE, dtype=DTYPE) / 100.0

    def build_one_side(descending: bool) -> torch.Tensor:
        order = torch.argsort(vals, descending=descending, stable=True)
        vals_ord = vals[order]
        idx_ord = idx[order]
        csum = torch.cumsum(vals_ord, dim=0)
        masses = total * retain_ratios
        keep_counts = torch.searchsorted(csum, masses, right=True).to(torch.int64)
        keep_counts = torch.where(keep_counts <= 0, torch.ones_like(keep_counts), keep_counts)
        col_ix = torch.arange(vals_ord.numel(), device=DEVICE, dtype=torch.int64).unsqueeze(0)
        keep_mask = col_ix < keep_counts.unsqueeze(1)
        kept_vals = keep_mask.to(DTYPE) * vals_ord.unsqueeze(0)
        side = torch.zeros((99, n), dtype=DTYPE, device=DEVICE)
        side[:, idx_ord] = kept_vals
        side_sum = side.sum(dim=1, keepdim=True)
        side = torch.where(side_sum > 0, side / side_sum, torch.zeros_like(side))
        return side

    low_side = build_one_side(descending=False)
    high_side = build_one_side(descending=True)
    none_vec = torch.zeros(n, dtype=DTYPE, device=DEVICE)
    none_vec[idx] = vals / total

    out_full[:99] = low_side
    out_full[99] = none_vec
    out_full[100:] = high_side
    return out_full[ACTIVE_CUTOFF_SLICE].clone()


def _betacf(a: float, b: float, x: float, max_iter: int = 200, eps: float = 3e-14) -> float:
    fpmin = 1e-300
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        h *= d * c

        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def _regularized_incomplete_beta(x: float, a: float, b: float) -> float:
    if not math.isfinite(x):
        return float("nan")
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    log_bt = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x)
    bt = math.exp(log_bt)
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def _student_t_two_sided_p(statistic: float, df: float) -> float:
    if not (math.isfinite(statistic) and math.isfinite(df)) or df <= 0:
        return float("nan")
    t_abs = abs(float(statistic))
    if t_abs <= 0:
        return 1.0
    x = df / (df + t_abs * t_abs)
    p = _regularized_incomplete_beta(x, df / 2.0, 0.5)
    return min(max(p, 0.0), 1.0)


def _rankdata_average_1d(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    n = arr.size
    if n == 0:
        return np.empty(0, dtype=float)
    sorter = np.argsort(arr, kind="mergesort")
    sorted_arr = arr[sorter]
    ranks_sorted = np.empty(n, dtype=float)
    i = 0
    while i < n:
        j = i + 1
        while j < n and sorted_arr[j] == sorted_arr[i]:
            j += 1
        rank = 0.5 * (i + j - 1) + 1.0
        ranks_sorted[i:j] = rank
        i = j
    out = np.empty(n, dtype=float)
    out[sorter] = ranks_sorted
    return out


def set_global_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def filter_by_cumulative_tail_tensor(vec: torch.Tensor, retain_ratio: float, tail: str) -> torch.Tensor:
    if vec.numel() == 0:
        return vec
    if tail == "none" or retain_ratio >= 1.0 - 1e-12:
        return vec.clone()
    if retain_ratio <= 1e-12 or torch.nansum(vec) <= 0:
        return torch.zeros_like(vec)

    valid_mask = torch.isfinite(vec) & (vec > 0)
    if not bool(valid_mask.any().item()):
        return torch.zeros_like(vec)

    idx = torch.nonzero(valid_mask, as_tuple=False).squeeze(1)
    vals = vec[idx]
    order = torch.argsort(vals, descending=(tail == "high"), stable=True)
    vals_ord = vals[order]
    idx_ord = idx[order]
    cutoff_mass = vals_ord.sum() * retain_ratio
    keep_mask = torch.cumsum(vals_ord, dim=0) <= cutoff_mass
    if not bool(keep_mask.any().item()) and retain_ratio > 0:
        keep_mask[0] = True
    out = torch.zeros_like(vec)
    out[idx_ord[keep_mask]] = vals_ord[keep_mask]
    return out


def make_model_based_cov_fn(p: torch.Tensor):
    p = p[p > 0].to(dtype=SOLVE_DTYPE, device=DEVICE)
    p = p / p.sum()

    def cov(m_values: Sequence[int] | int | torch.Tensor) -> torch.Tensor:
        m_tensor = torch.as_tensor(m_values, dtype=SOLVE_DTYPE, device=p.device)
        scalar_input = m_tensor.ndim == 0
        if scalar_input:
            m_tensor = m_tensor.reshape(1)
        out = 1.0 - torch.sum(p.unsqueeze(0) * torch.pow(1.0 - p.unsqueeze(0), m_tensor.unsqueeze(1)), dim=1)
        return out[0] if scalar_input else out

    return cov


def solve_m_expected(p: torch.Tensor, cstar: float, m_max: int = N_UPPER_CAP) -> int:
    if p.numel() == 0 or float(torch.nansum(p).item()) <= 0:
        return 0
    cov = make_model_based_cov_fn(p)
    lo = 0
    hi = 1
    while float(cov(hi).item()) < cstar and hi < m_max:
        hi *= 2
    if hi >= m_max:
        return int(m_max)
    for _ in range(30):
        mid = (lo + hi) // 2
        if float(cov(mid).item()) < cstar:
            lo = mid
        else:
            hi = mid
    return int(hi)


def make_chao_inext_cov_fn(counts: torch.Tensor):
    x = counts[counts > 0].to(dtype=SOLVE_DTYPE, device=DEVICE)
    if x.numel() == 0:
        def zero_cov(m_values: Sequence[int] | int | torch.Tensor) -> torch.Tensor:
            m_tensor = torch.as_tensor(m_values, dtype=SOLVE_DTYPE, device=DEVICE)
            return torch.zeros_like(m_tensor)
        return zero_cov

    n = x.sum()
    n_int = int(round(float(n.item())))
    f1 = (x == 1).sum().to(DTYPE)
    f2 = (x == 2).sum().to(DTYPE)
    if float(f2.item()) == 0:
        f0_hat = ((n - 1.0) / n) * f1 * torch.clamp(f1 - 1.0, min=0.0) / 2.0
    else:
        f0_hat = ((n - 1.0) / n) * (f1 ** 2) / (2.0 * f2)
    if float(f1.item()) > 0:
        a_hat = n * f0_hat / (n * f0_hat + f1)
    else:
        a_hat = torch.tensor(1.0, dtype=SOLVE_DTYPE, device=DEVICE)

    def cov_int(m_int_tensor: torch.Tensor) -> torch.Tensor:
        if m_int_tensor.numel() == 0:
            return torch.empty(0, dtype=DTYPE, device=DEVICE)
        mm = m_int_tensor.to(dtype=SOLVE_DTYPE, device=DEVICE)
        out = torch.empty(mm.shape[0], dtype=SOLVE_DTYPE, device=DEVICE)
        lt_mask = mm < n_int
        eq_mask = mm == n_int
        gt_mask = mm > n_int

        if bool(lt_mask.any().item()):
            mm_lt = mm[lt_mask]
            xx = x.unsqueeze(0)
            mm_col = mm_lt.unsqueeze(1)
            valid = (n - xx) >= mm_col
            term = torch.zeros((mm_lt.shape[0], x.shape[0]), dtype=SOLVE_DTYPE, device=DEVICE)
            if bool(valid.any().item()):
                log_term = (
                    torch.lgamma((n - xx) + 1.0)
                    - torch.lgamma((n - xx - mm_col) + 1.0)
                    - math.lgamma(float(n.item()))
                    + torch.lgamma(n - mm_col)
                )
                term = torch.where(valid, (xx / n) * torch.exp(log_term), torch.zeros_like(term))
            out[lt_mask] = 1.0 - term.sum(dim=1)

        if bool(eq_mask.any().item()):
            out[eq_mask] = 1.0 - (f1 / n) * a_hat

        if bool(gt_mask.any().item()):
            mm_gt = mm[gt_mask]
            out[gt_mask] = 1.0 - (f1 / n) * torch.pow(a_hat, mm_gt - n + 1.0)

        return torch.clamp(out, 0.0, 1.0)

    def cov(m_values: Sequence[int] | int | torch.Tensor) -> torch.Tensor:
        m_tensor = torch.as_tensor(m_values, dtype=SOLVE_DTYPE, device=DEVICE)
        scalar_input = m_tensor.ndim == 0
        if scalar_input:
            m_tensor = m_tensor.reshape(1)
        m_flat = m_tensor.reshape(-1)
        m_floor = torch.floor(m_flat).to(torch.int64)
        m_ceil = torch.ceil(m_flat).to(torch.int64)
        same_mask = m_floor == m_ceil
        out = torch.empty(m_flat.shape[0], dtype=SOLVE_DTYPE, device=DEVICE)
        if bool(same_mask.any().item()):
            out[same_mask] = cov_int(m_floor[same_mask])
        if bool((~same_mask).any().item()):
            floor_vals = cov_int(m_floor[~same_mask])
            ceil_vals = cov_int(m_ceil[~same_mask])
            frac = m_flat[~same_mask] - m_floor[~same_mask].to(DTYPE)
            out[~same_mask] = (1.0 - frac) * floor_vals + frac * ceil_vals
        out = out.reshape(m_tensor.shape)
        return out[0] if scalar_input else out

    return cov


def calc_sc_chao_inext(counts: torch.Tensor, m: Optional[Sequence[int] | int | float | torch.Tensor] = None):
    if counts.numel() == 0 or float(counts.sum().item()) <= 0:
        if m is None:
            return 0.0
        m_tensor = torch.as_tensor(m, dtype=SOLVE_DTYPE, device=DEVICE)
        return torch.zeros_like(m_tensor)
    cov = make_chao_inext_cov_fn(counts)
    if m is None:
        m = int(round(float(counts.sum().item())))
    out = cov(m)
    if isinstance(out, torch.Tensor) and out.ndim == 0:
        return float(out.item())
    return out


def solve_m_chao_inext_grid(counts: torch.Tensor, targets: Sequence[float], m_max: int = N_UPPER_CAP) -> List[int]:
    x = counts[counts > 0]
    if x.numel() == 0 or float(x.sum().item()) <= 0:
        return [0 for _ in targets]

    target_tensor = torch.as_tensor(targets, dtype=SOLVE_DTYPE, device=DEVICE)
    cov = make_chao_inext_cov_fn(x)
    lo = torch.zeros_like(target_tensor, dtype=torch.int64)
    hi = torch.ones_like(target_tensor, dtype=torch.int64)
    cov_cache: Dict[int, float] = {}

    def cached_cov(m_tensor: torch.Tensor) -> torch.Tensor:
        flat = m_tensor.detach().to("cpu", dtype=torch.int64).reshape(-1).tolist()
        missing = sorted({int(v) for v in flat if int(v) not in cov_cache})
        if missing:
            cov_vals = cov(torch.as_tensor(missing, dtype=SOLVE_DTYPE, device=DEVICE)).detach().to("cpu").tolist()
            for m_val, c_val in zip(missing, cov_vals):
                cov_cache[int(m_val)] = float(c_val)
        out = torch.as_tensor([cov_cache[int(v)] for v in flat], dtype=SOLVE_DTYPE, device=DEVICE)
        return out.reshape(m_tensor.shape)

    while True:
        cov_hi = cached_cov(hi)
        need = (cov_hi < target_tensor) & (hi < m_max)
        if not bool(need.any().item()):
            break
        hi = torch.where(need, hi * 2, hi)

    hi = torch.clamp(hi, max=m_max)
    cov_hi = cached_cov(hi)
    cap_mask = (hi >= m_max) & (cov_hi < target_tensor)
    lo = torch.where(cap_mask, hi, lo)

    for _ in range(30):
        mid = (lo + hi) // 2
        cov_mid = cached_cov(mid)
        lo = torch.where(cov_mid < target_tensor, mid, lo)
        hi = torch.where(cov_mid < target_tensor, hi, mid)

    out = torch.where(cap_mask, torch.full_like(hi, m_max), hi)
    return [int(v) for v in out.detach().cpu().tolist()]


def calc_sc_empirical_batch(count_batch: torch.Tensor) -> torch.Tensor:
    counts = count_batch.to(torch.int64)
    n = counts.sum(dim=1).to(SOLVE_DTYPE)
    f1 = (counts == 1).sum(dim=1).to(SOLVE_DTYPE)
    out = torch.zeros_like(n)
    valid = n > 0
    out[valid] = 1.0 - (f1[valid] / n[valid])
    no_singletons = valid & (f1 == 0)
    out[no_singletons] = 1.0
    return out


def calc_sc_chao_inext_batch(count_batch: torch.Tensor) -> torch.Tensor:
    counts = count_batch.to(torch.int64)
    n = counts.sum(dim=1).to(SOLVE_DTYPE)
    valid = n > 0
    out = torch.zeros(counts.shape[0], dtype=SOLVE_DTYPE, device=counts.device)
    if not bool(valid.any().item()):
        return out
    counts_v = counts[valid]
    n_v = n[valid]
    f1 = (counts_v == 1).sum(dim=1).to(SOLVE_DTYPE)
    f2 = (counts_v == 2).sum(dim=1).to(SOLVE_DTYPE)
    f0_hat = torch.where(
        f2 == 0,
        ((n_v - 1.0) / n_v) * f1 * torch.clamp(f1 - 1.0, min=0.0) / 2.0,
        ((n_v - 1.0) / n_v) * (f1 ** 2) / (2.0 * f2),
    )
    a_hat = torch.where(f1 > 0, n_v * f0_hat / (n_v * f0_hat + f1), torch.ones_like(f1))
    out_valid = 1.0 - (f1 / n_v) * a_hat
    out[valid] = torch.clamp(out_valid, 0.0, 1.0)
    return out


def calc_chao1_asymptotic_batch(count_batch: torch.Tensor) -> torch.Tensor:
    counts = count_batch.to(torch.int64)
    n = counts.sum(dim=1).to(SOLVE_DTYPE)
    s_obs = (counts > 0).sum(dim=1).to(SOLVE_DTYPE)
    f1 = (counts == 1).sum(dim=1).to(SOLVE_DTYPE)
    f2 = (counts == 2).sum(dim=1).to(SOLVE_DTYPE)
    with torch.no_grad():
        f0 = torch.where(
            f2 > 0,
            ((n - 1.0) / n) * (f1 ** 2) / (2.0 * f2),
            ((n - 1.0) / n) * (f1 * torch.clamp(f1 - 1.0, min=0.0) / 2.0),
        )
    val = s_obs + f0
    return torch.where(torch.isnan(val), s_obs, val)


def calc_hill_metrics_batch(count_batch: torch.Tensor) -> Dict[str, torch.Tensor]:
    counts = count_batch.to(DTYPE)
    positive = counts > 0
    n = counts.sum(dim=1, keepdim=True)
    q0 = positive.sum(dim=1).to(DTYPE)
    p = torch.where(positive, counts / torch.clamp(n, min=1.0), torch.zeros_like(counts))
    logp = torch.where(positive, torch.log(p), torch.zeros_like(p))
    entropy = -(p * logp).sum(dim=1)
    q1 = torch.exp(entropy)
    q2 = 1.0 / torch.clamp((p ** 2).sum(dim=1), min=1e-12)
    q0_chao = calc_chao1_asymptotic_batch(count_batch)
    return {
        "q0_observed_standardized": q0,
        "q0_chao1_asymptotic": q0_chao,
        "q1_empirical_hill": q1,
        "q2_empirical_hill": q2,
    }


def calc_q0_stage2_matrix(
    count_batch: torch.Tensor,
    sizes: Sequence[float],
    endpoint_mask: Sequence[bool],
) -> Dict[str, torch.Tensor]:
    counts = count_batch.to(SOLVE_DTYPE)
    device = counts.device
    n = counts.sum(dim=1)
    n_col = n.reshape(-1, 1)
    m, valid = make_stage2_m_matrix(n, sizes, endpoint_mask, device)
    endpoint_t = torch.as_tensor(endpoint_mask, dtype=torch.bool, device=device).reshape(1, -1)

    positive = counts > 0
    s_obs = positive.sum(dim=1).to(SOLVE_DTYPE)
    f1 = (counts == 1).sum(dim=1).to(SOLVE_DTYPE)
    f2 = (counts == 2).sum(dim=1).to(SOLVE_DTYPE)
    valid_rows = n > 0
    f0_hat = torch.zeros_like(n)
    f0_f2 = ((n - 1.0) / torch.clamp(n, min=1.0)) * (f1 ** 2) / (2.0 * torch.clamp(f2, min=1.0))
    f0_no_f2 = ((n - 1.0) / torch.clamp(n, min=1.0)) * f1 * torch.clamp(f1 - 1.0, min=0.0) / 2.0
    f0_hat = torch.where(valid_rows & (f2 > 0), f0_f2, f0_no_f2)
    f0_hat = torch.where((valid_rows & (f1 > 0)), torch.clamp(f0_hat, min=0.0), torch.zeros_like(f0_hat))
    q0_asym = s_obs + f0_hat

    q0_obs = s_obs.reshape(-1, 1).expand_as(m).clone()
    lt_mask = valid & (n_col > 1.0) & (m < n_col) & (~endpoint_t.expand_as(m))
    if bool(lt_mask.any().item()):
        x = counts.reshape(counts.shape[0], 1, counts.shape[1])
        n3 = n.reshape(-1, 1, 1)
        m3 = m.reshape(m.shape[0], m.shape[1], 1)
        possible = (x > 0) & (m3 <= (n3 - x))
        log_ratio = (
            torch.lgamma(n3 - x + 1.0)
            - torch.lgamma(n3 - x - m3 + 1.0)
            + torch.lgamma(n3 - m3 + 1.0)
            - torch.lgamma(n3 + 1.0)
        )
        ratio = torch.where(possible, torch.exp(log_ratio), torch.zeros_like(log_ratio))
        interp_terms = torch.where(x > 0, 1.0 - ratio, torch.zeros_like(ratio))
        interp_q0 = interp_terms.sum(dim=2)
        q0_obs = torch.where(lt_mask, interp_q0, q0_obs)

    gt_mask = valid & (n_col > 0.0) & (m > n_col) & (~endpoint_t.expand_as(m))
    if bool(gt_mask.any().item()):
        beta = torch.where(
            (f0_hat > 0.0) & (f1 > 0.0),
            f1 / torch.clamp(n * f0_hat + f1, min=1e-300),
            torch.zeros_like(f1),
        ).reshape(-1, 1)
        extra = torch.clamp(m - n_col, min=0.0)
        extra_q0 = s_obs.reshape(-1, 1) + f0_hat.reshape(-1, 1) * (1.0 - torch.pow(1.0 - beta, extra))
        q0_obs = torch.where(gt_mask, extra_q0, q0_obs)

    q0_obs = torch.where(valid, q0_obs, torch.full_like(q0_obs, float("nan")))
    q0_asym_matrix = q0_asym.reshape(-1, 1).expand_as(m)
    q0_asym_matrix = torch.where(valid, q0_asym_matrix, torch.full_like(q0_asym_matrix, float("nan")))
    return {
        "q0_observed_standardized": q0_obs,
        "q0_chao1_asymptotic": q0_asym_matrix,
    }


def entropy_hill_from_counts(counts: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
    x = counts[counts > 0].to(SOLVE_DTYPE)
    if x.numel() == 0:
        nan = torch.tensor(float("nan"), dtype=SOLVE_DTYPE, device=counts.device)
        return nan, nan
    n = torch.clamp(x.sum(), min=1.0)
    p = x / n
    q1 = torch.exp(-(p * torch.log(p)).sum())
    q2 = 1.0 / torch.clamp((p ** 2).sum(), min=1e-300)
    return q1, q2


def hill_q2_from_counts(counts: torch.Tensor) -> torch.Tensor:
    x = counts[counts > 0].to(SOLVE_DTYPE)
    if x.numel() == 0:
        return torch.tensor(float("nan"), dtype=SOLVE_DTYPE, device=counts.device)
    p = x / torch.clamp(x.sum(), min=1.0)
    return 1.0 / torch.clamp((p ** 2).sum(), min=1e-300)


def chao_f0_a_from_counts(counts: torch.Tensor) -> Tuple[float, float, float, float, float, float]:
    x = counts[counts > 0].to(SOLVE_DTYPE)
    if x.numel() == 0:
        return 0.0, 0.0, 0.0, 1.0, 0.0, 0.0
    n = float(x.sum().item())
    f1 = float((x == 1).sum().item())
    f2 = float((x == 2).sum().item())
    if n <= 0.0 or not math.isfinite(n):
        return 0.0, 0.0, 0.0, 1.0, 0.0, 0.0
    if f2 > 0:
        f0_hat = ((n - 1.0) / n) * (f1 ** 2) / (2.0 * f2)
    else:
        f0_hat = ((n - 1.0) / n) * f1 * max(f1 - 1.0, 0.0) / 2.0
    a_hat = (n * f0_hat / (n * f0_hat + f1)) if f1 > 0 else 1.0
    unseen_mass = (f1 / n) * a_hat if f1 > 0 else 0.0
    return n, f1, f2, a_hat, max(float(f0_hat), 0.0), min(max(float(unseen_mass), 0.0), 1.0)


def inext_q1_second_order(n: int, f1: float, f2: float, device: torch.device) -> torch.Tensor:
    if f1 <= 0 or n <= 1:
        return torch.tensor(0.0, dtype=SOLVE_DTYPE, device=device)
    if f2 > 0:
        a = 2.0 * f2 / ((n - 1.0) * f1 + 2.0 * f2)
    else:
        a = 2.0 / ((n - 1.0) * max(f1 - 1.0, 0.0) + 2.0)
    if a >= 1.0:
        return torch.tensor(0.0, dtype=SOLVE_DTYPE, device=device)
    r = torch.arange(1, n, dtype=SOLVE_DTYPE, device=device)
    partial = (torch.pow(1.0 - a, r) / r).sum()
    tail = max(-math.log(a) - float(partial.item()), 0.0)
    if tail <= 0.0:
        return torch.tensor(0.0, dtype=SOLVE_DTYPE, device=device)
    log_val = math.log(f1 / n) + (-n + 1) * math.log1p(-a) + math.log(tail)
    if log_val > 700:
        val = float("inf")
    else:
        val = math.exp(log_val)
    if not math.isfinite(val):
        return torch.tensor(0.0, dtype=SOLVE_DTYPE, device=device)
    return torch.tensor(max(val, 0.0), dtype=SOLVE_DTYPE, device=device)


def inext_asymptotic_q1_q2(counts: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
    x = counts[counts > 0].to(SOLVE_DTYPE)
    if x.numel() == 0:
        nan = torch.tensor(float("nan"), dtype=SOLVE_DTYPE, device=counts.device)
        return nan, nan
    n_float, f1, f2, _a_hat, _f0_hat, _unseen = chao_f0_a_from_counts(x)
    n = int(round(n_float))
    unique_x, freq = torch.unique(x, sorted=True, return_counts=True)
    tab = freq.to(SOLVE_DTYPE)
    term = (tab * unique_x / n_float * (torch.digamma(torch.tensor(n_float, dtype=SOLVE_DTYPE, device=counts.device)) - torch.digamma(unique_x))).sum()
    q1 = torch.exp(term + inext_q1_second_order(n, f1, f2, counts.device))
    denom = (x * torch.clamp(x - 1.0, min=0.0)).sum()
    if n <= 1 or float(denom.item()) <= 0.0:
        _q1_obs, q2_obs = entropy_hill_from_counts(x)
        q2 = q2_obs
    else:
        q2 = (n_float * (n_float - 1.0)) / denom
    return q1, q2


def inext_rtd_q1_interpolation(counts: torch.Tensor, m: int) -> torch.Tensor:
    x = counts[counts > 0].to(SOLVE_DTYPE)
    if x.numel() == 0 or m <= 0:
        return torch.tensor(float("nan"), dtype=SOLVE_DTYPE, device=counts.device)
    n = int(round(float(x.sum().item())))
    if m >= n:
        q1, _q2 = entropy_hill_from_counts(x)
        return q1
    unique_x, freq = torch.unique(x, sorted=True, return_counts=True)
    fhat = torch.zeros(m, dtype=SOLVE_DTYPE, device=counts.device)
    log_denom = torch.lgamma(torch.tensor(float(n + 1), dtype=SOLVE_DTYPE, device=counts.device)) - torch.lgamma(torch.tensor(float(m + 1), dtype=SOLVE_DTYPE, device=counts.device)) - torch.lgamma(torch.tensor(float(n - m + 1), dtype=SOLVE_DTYPE, device=counts.device))
    for z_raw, f_raw in zip(unique_x.detach().cpu().tolist(), freq.detach().cpu().tolist()):
        z = int(round(float(z_raw)))
        k_start = max(1, m - (n - z))
        k_end = min(m, z)
        if k_start > k_end:
            continue
        k = torch.arange(k_start, k_end + 1, dtype=SOLVE_DTYPE, device=counts.device)
        z_t = torch.full_like(k, float(z))
        log_num = (
            torch.lgamma(z_t + 1.0)
            - torch.lgamma(k + 1.0)
            - torch.lgamma(z_t - k + 1.0)
            + torch.lgamma(torch.tensor(float(n - z + 1), dtype=SOLVE_DTYPE, device=counts.device))
            - torch.lgamma(torch.tensor(float(n - z), dtype=SOLVE_DTYPE, device=counts.device) - (torch.tensor(float(m), dtype=SOLVE_DTYPE, device=counts.device) - k) + 1.0)
            - torch.lgamma(torch.tensor(float(m), dtype=SOLVE_DTYPE, device=counts.device) - k + 1.0)
        )
        fhat[k_start - 1 : k_end] += float(f_raw) * torch.exp(log_num - log_denom)
    k_all = torch.arange(1, m + 1, dtype=SOLVE_DTYPE, device=counts.device)
    prop = k_all / float(m)
    entropy = (-(prop * torch.log(prop)) * fhat).sum()
    return torch.exp(entropy)


def inext_rtd_nm1_q1(counts: torch.Tensor) -> torch.Tensor:
    x = counts[counts > 0].to(SOLVE_DTYPE)
    if x.numel() == 0:
        return torch.tensor(float("nan"), dtype=SOLVE_DTYPE, device=counts.device)
    n = float(x.sum().item())
    if n <= 1.0:
        q1, _q2 = entropy_hill_from_counts(x)
        return q1
    keep_same = (n - x) / n
    remove_one = x / n
    p_same = x / (n - 1.0)
    h_same = -p_same * torch.log(p_same)
    x_minus = torch.clamp(x - 1.0, min=0.0)
    p_minus = x_minus / (n - 1.0)
    h_minus = torch.where(x_minus > 0, -p_minus * torch.log(torch.clamp(p_minus, min=1e-300)), torch.zeros_like(p_minus))
    return torch.exp((keep_same * h_same + remove_one * h_minus).sum())


def inext_q2_interpolation(counts: torch.Tensor, m: int) -> torch.Tensor:
    x = counts[counts > 0].to(SOLVE_DTYPE)
    if x.numel() == 0 or m <= 0:
        return torch.tensor(float("nan"), dtype=SOLVE_DTYPE, device=counts.device)
    n = float(x.sum().item())
    if n <= 1.0:
        return hill_q2_from_counts(x)
    p = x / n
    s2 = (p ** 2).sum()
    denom = s2 + ((n - float(m)) / (float(m) * (n - 1.0))) * (1.0 - s2)
    return 1.0 / torch.clamp(denom, min=1e-300)


def inext_q1_extrapolation(counts: torch.Tensor, m_total: int) -> torch.Tensor:
    obs_q1, _obs_q2 = entropy_hill_from_counts(counts)
    asy_q1, _asy_q2 = inext_asymptotic_q1_q2(counts)
    rfd_nm1 = inext_rtd_nm1_q1(counts)
    denom = asy_q1 - rfd_nm1
    if not bool(torch.isfinite(denom).item()) or abs(float(denom.item())) < 1e-12:
        return obs_q1
    beta = float(((obs_q1 - rfd_nm1) / denom).item())
    n = int(round(float(counts[counts > 0].sum().item())))
    m_star = max(0, int(round(m_total)) - n)
    factor = 1.0 - ((1.0 - beta) ** m_star)
    val = float(obs_q1.item()) + (float(asy_q1.item()) - float(obs_q1.item())) * factor
    if not math.isfinite(val):
        return obs_q1
    return torch.tensor(max(val, 0.0), dtype=SOLVE_DTYPE, device=counts.device)


def inext_q2_extrapolation(counts: torch.Tensor, m_total: int) -> torch.Tensor:
    x = counts[counts > 0].to(SOLVE_DTYPE)
    if x.numel() == 0:
        return torch.tensor(float("nan"), dtype=SOLVE_DTYPE, device=counts.device)
    n = float(x.sum().item())
    if n <= 1.0:
        return hill_q2_from_counts(x)
    m_total = max(1, int(round(m_total)))
    common = (x / n * torch.clamp(x - 1.0, min=0.0) / (n - 1.0)).sum()
    denom = (1.0 / float(m_total)) + (1.0 - 1.0 / float(m_total)) * common
    return 1.0 / torch.clamp(denom, min=1e-300)


def calc_inext_td_metrics_single(
    counts: torch.Tensor,
    size: Optional[int],
    is_endpoint: bool,
    want_q1: bool = NEED_Q1_INEXT,
    want_q2: bool = NEED_Q2_INEXT,
) -> Tuple[torch.Tensor, torch.Tensor]:
    x = counts[counts > 0].to(SOLVE_DTYPE)
    if x.numel() == 0:
        nan = torch.tensor(float("nan"), dtype=SOLVE_DTYPE, device=counts.device)
        return nan, nan
    n = int(round(float(x.sum().item())))
    nan = torch.tensor(float("nan"), dtype=SOLVE_DTYPE, device=counts.device)
    q1 = nan
    q2 = nan
    if is_endpoint:
        asy_q1, asy_q2 = inext_asymptotic_q1_q2(x)
        if want_q1:
            q1 = asy_q1
        if want_q2:
            q2 = asy_q2
        return q1, q2
    m = n if size is None else max(1, int(round(size)))
    if m < n:
        if want_q1:
            q1 = inext_rtd_q1_interpolation(x, m)
        if want_q2:
            q2 = inext_q2_interpolation(x, m)
    elif m == n:
        if want_q1:
            q1, q2_emp = entropy_hill_from_counts(x)
            if want_q2:
                q2 = q2_emp
        elif want_q2:
            q2 = hill_q2_from_counts(x)
    else:
        if want_q1:
            q1 = inext_q1_extrapolation(x, m)
        if want_q2:
            q2 = inext_q2_extrapolation(x, m)
    return q1, q2


def calc_inext_td_metrics_batch(count_batch: torch.Tensor, size: Optional[int], is_endpoint: bool) -> Dict[str, torch.Tensor]:
    q1_vals: List[torch.Tensor] = []
    q2_vals: List[torch.Tensor] = []
    for row in count_batch:
        q1, q2 = calc_inext_td_metrics_single(row, size, is_endpoint)
        if NEED_Q1_INEXT:
            q1_vals.append(q1.to(DTYPE))
        if NEED_Q2_INEXT:
            q2_vals.append(q2.to(DTYPE))
    empty = torch.empty(0, dtype=DTYPE, device=count_batch.device)
    out: Dict[str, torch.Tensor] = {}
    if NEED_Q1_INEXT:
        out["q1_inext_td_m_est"] = torch.stack(q1_vals) if q1_vals else empty
    if NEED_Q2_INEXT:
        out["q2_inext_td_m_est"] = torch.stack(q2_vals) if q2_vals else empty
    return out


def make_stage2_m_matrix(
    n: torch.Tensor,
    sizes: Sequence[float],
    endpoint_mask: Sequence[bool],
    device: torch.device,
) -> Tuple[torch.Tensor, torch.Tensor]:
    sizes_t = torch.as_tensor(sizes, dtype=SOLVE_DTYPE, device=device)
    endpoint_t = torch.as_tensor(endpoint_mask, dtype=torch.bool, device=device)
    n_col = n.to(SOLVE_DTYPE).reshape(-1, 1)
    m = sizes_t.reshape(1, -1).expand(n_col.shape[0], -1).clone()
    m[:, endpoint_t] = n_col.expand(-1, int(endpoint_t.sum().item())) if bool(endpoint_t.any().item()) else m[:, endpoint_t]
    valid = torch.isfinite(m) & (m > 0) & (n_col > 0)
    m = torch.where(valid, torch.clamp(torch.round(m), min=1.0), torch.ones_like(m))
    return m, valid


def calc_inext_q2_matrix(
    count_batch: torch.Tensor,
    sizes: Sequence[float],
    endpoint_mask: Sequence[bool],
) -> torch.Tensor:
    counts = count_batch.to(SOLVE_DTYPE)
    device = counts.device
    n = counts.sum(dim=1)
    m, valid = make_stage2_m_matrix(n, sizes, endpoint_mask, device)
    n_col = n.reshape(-1, 1)
    p = counts / torch.clamp(n_col, min=1.0)
    s2 = (p ** 2).sum(dim=1).reshape(-1, 1)
    denom_emp = s2.expand_as(m)

    n_minus_1 = torch.clamp(n_col - 1.0, min=1.0)
    interp_denom = s2 + ((n_col - m) / (m * n_minus_1)) * (1.0 - s2)
    common = (counts * torch.clamp(counts - 1.0, min=0.0)).sum(dim=1).reshape(-1, 1) / torch.clamp(
        n_col * (n_col - 1.0),
        min=1.0,
    )
    extra_denom = (1.0 / m) + (1.0 - 1.0 / m) * common
    asym_denom = (
        (counts * torch.clamp(counts - 1.0, min=0.0)).sum(dim=1).reshape(-1, 1)
        / torch.clamp(n_col * (n_col - 1.0), min=1.0)
    )
    asym_q2 = torch.where(
        (n_col > 1.0) & (asym_denom > 0.0),
        1.0 / torch.clamp(asym_denom, min=1e-300),
        1.0 / torch.clamp(denom_emp, min=1e-300),
    )

    lt_mask = valid & (n_col > 1.0) & (m < n_col)
    gt_mask = valid & (n_col > 1.0) & (m > n_col)
    denom = torch.where(lt_mask, interp_denom, denom_emp)
    denom = torch.where(gt_mask, extra_denom, denom)
    q2 = 1.0 / torch.clamp(denom, min=1e-300)
    endpoint_t = torch.as_tensor(endpoint_mask, dtype=torch.bool, device=device).reshape(1, -1)
    q2 = torch.where(valid & endpoint_t.expand_as(q2), asym_q2.expand_as(q2), q2)
    return torch.where(valid, q2, torch.full_like(q2, float("nan")))


def calc_sc_chao_inext_matrix(
    count_batch: torch.Tensor,
    sizes: Sequence[float],
    endpoint_mask: Sequence[bool],
) -> torch.Tensor:
    counts = count_batch.to(SOLVE_DTYPE)
    device = counts.device
    n = counts.sum(dim=1)
    m, valid = make_stage2_m_matrix(n, sizes, endpoint_mask, device)
    n_col = n.reshape(-1, 1)
    f1 = (counts == 1).sum(dim=1).to(SOLVE_DTYPE).reshape(-1, 1)
    f2 = (counts == 2).sum(dim=1).to(SOLVE_DTYPE).reshape(-1, 1)
    f0_hat = torch.where(
        f2 > 0,
        ((n_col - 1.0) / torch.clamp(n_col, min=1.0)) * (f1 ** 2) / (2.0 * torch.clamp(f2, min=1.0)),
        ((n_col - 1.0) / torch.clamp(n_col, min=1.0)) * f1 * torch.clamp(f1 - 1.0, min=0.0) / 2.0,
    )
    a_hat = torch.where(f1 > 0, n_col * f0_hat / torch.clamp(n_col * f0_hat + f1, min=1e-300), torch.ones_like(f1))

    base = 1.0 - (f1 / torch.clamp(n_col, min=1.0)) * a_hat
    out = base.expand_as(m).clone()
    endpoint_t = torch.as_tensor(endpoint_mask, dtype=torch.bool, device=device).reshape(1, -1)
    out = torch.where(valid & endpoint_t.expand_as(out), torch.ones_like(out), out)

    lt_mask = valid & (n_col > 1.0) & (m < n_col)
    if bool(lt_mask.any().item()):
        x = counts.reshape(counts.shape[0], 1, counts.shape[1])
        n3 = n.reshape(-1, 1, 1)
        m3 = m.reshape(m.shape[0], m.shape[1], 1)
        positive = x > 0
        possible = positive & (m3 <= (n3 - x))
        log_ratio = (
            torch.lgamma(n3 - x + 1.0)
            - torch.lgamma(n3 - x - m3 + 1.0)
            + torch.lgamma(n3 - m3)
            - torch.lgamma(n3)
        )
        ratio = torch.where(possible, torch.exp(log_ratio), torch.zeros_like(log_ratio))
        terms = (x / torch.clamp(n3, min=1.0)) * ratio
        interp_cov = 1.0 - terms.sum(dim=2)
        out = torch.where(lt_mask, interp_cov, out)

    gt_mask = valid & (n_col > 1.0) & (m > n_col)
    if bool(gt_mask.any().item()):
        extra_cov = 1.0 - (f1 / torch.clamp(n_col, min=1.0)) * torch.pow(a_hat, m - n_col + 1.0)
        out = torch.where(gt_mask, extra_cov, out)

    out = torch.clamp(out, 0.0, 1.0)
    return torch.where(valid, out, torch.full_like(out, float("nan")))


def _get_scipy_gammaln():
    try:
        from scipy.special import gammaln  # type: ignore

        return gammaln
    except Exception:
        return None


def _q1_interp_np_for_row(
    positive_counts: np.ndarray,
    unique_x: np.ndarray,
    freq: np.ndarray,
    log_factorial: np.ndarray,
    n: int,
    m: int,
) -> float:
    if m <= 0 or n <= 0:
        return float("nan")
    if m >= n:
        p = positive_counts.astype(float) / float(n)
        return float(math.exp(float(-(p * np.log(np.clip(p, 1e-300, None))).sum())))

    log_denom = log_factorial[n] - log_factorial[m] - log_factorial[n - m]
    entropy = 0.0
    for z_raw, f_raw in zip(unique_x, freq):
        z = int(z_raw)
        k_start = max(1, m - (n - z))
        k_end = min(m, z)
        if k_start > k_end:
            continue
        k = np.arange(k_start, k_end + 1, dtype=np.int64)
        m_minus_k = m - k
        log_num = (
            log_factorial[z]
            - log_factorial[k]
            - log_factorial[z - k]
            + log_factorial[n - z]
            - log_factorial[m_minus_k]
            - log_factorial[n - z - m_minus_k]
        )
        prop = k.astype(float) / float(m)
        entropy += float(f_raw) * float(np.sum((-(prop * np.log(prop))) * np.exp(log_num - log_denom)))
    return float(math.exp(entropy)) if math.isfinite(entropy) else float("nan")


def calc_inext_q1_matrix(
    count_batch: torch.Tensor,
    sizes: Sequence[float],
    endpoint_mask: Sequence[bool],
) -> torch.Tensor:
    gammaln = _get_scipy_gammaln()
    device = count_batch.device
    rows = count_batch.detach().cpu().numpy()
    max_np_n = max(1, parse_int_env("TOMATO_Q1_FAST_MAX_NP_N", 2_000_000))
    out = np.full((rows.shape[0], len(sizes)), np.nan, dtype=float)

    for row_idx, row_np_raw in enumerate(rows):
        positive = row_np_raw[row_np_raw > 0].astype(np.int64, copy=False)
        if positive.size == 0:
            continue
        n = int(round(float(positive.sum())))
        if n <= 0:
            continue
        row_t = count_batch[row_idx]
        obs_q1, _obs_q2 = entropy_hill_from_counts(row_t)
        obs = float(obs_q1.detach().cpu().item())
        asy_q1, _asy_q2 = inext_asymptotic_q1_q2(row_t)
        asy = float(asy_q1.detach().cpu().item())
        unique_x, freq = np.unique(positive, return_counts=True)

        log_factorial: Optional[np.ndarray]
        if gammaln is not None and n <= max_np_n:
            log_factorial = np.asarray(gammaln(np.arange(n + 1, dtype=float) + 1.0), dtype=float)
        else:
            log_factorial = None

        for col_idx, (size, is_endpoint) in enumerate(zip(sizes, endpoint_mask)):
            if is_endpoint:
                out[row_idx, col_idx] = asy
                continue
            if not math.isfinite(float(size)):
                continue
            m = max(1, int(round(float(size))))
            if m < n:
                if log_factorial is None:
                    out[row_idx, col_idx] = float(inext_rtd_q1_interpolation(row_t, m).detach().cpu().item())
                else:
                    out[row_idx, col_idx] = _q1_interp_np_for_row(positive, unique_x, freq, log_factorial, n, m)
            elif m == n:
                out[row_idx, col_idx] = obs
            else:
                out[row_idx, col_idx] = float(inext_q1_extrapolation(row_t, m).detach().cpu().item())

    return torch.as_tensor(out, dtype=SOLVE_DTYPE, device=device)


def build_stage2_q0_fast_payload(
    ref_acc: Optional[torch.Tensor],
    planned_m_values: Sequence[float],
    endpoint_mask: Sequence[bool],
) -> Optional[Dict[str, object]]:
    if ref_acc is None:
        return None
    batch = sample_inext_bootstrap_reference_batch(ref_acc, N_BOOT_MATRIX)
    if batch is None:
        return None
    q0_matrices = calc_q0_stage2_matrix(batch, planned_m_values, endpoint_mask)
    point_matrices = calc_q0_stage2_matrix(ref_acc.unsqueeze(0), planned_m_values, endpoint_mask)
    sc_matrix = calc_sc_chao_inext_matrix(batch, planned_m_values, endpoint_mask)
    return {
        "q0_observed_standardized": q0_matrices["q0_observed_standardized"].detach().cpu().numpy().astype(float),
        "q0_chao1_asymptotic": q0_matrices["q0_chao1_asymptotic"].detach().cpu().numpy().astype(float),
        "point": {
            "q0_observed_standardized": point_matrices["q0_observed_standardized"].detach().cpu().numpy().astype(float).reshape(-1),
            "q0_chao1_asymptotic": point_matrices["q0_chao1_asymptotic"].detach().cpu().numpy().astype(float).reshape(-1),
        },
        "sc_model": sc_matrix.detach().cpu().numpy().astype(float),
        "sc_emp": sc_matrix.detach().cpu().numpy().astype(float),
        "f1": (batch == 1).sum(dim=1).detach().cpu().numpy().astype(float),
        "n": batch.sum(dim=1).detach().cpu().numpy().astype(float),
    }


def build_stage2_q1_fast_payload(
    ref_acc: Optional[torch.Tensor],
    planned_m_values: Sequence[float],
    endpoint_mask: Sequence[bool],
) -> Optional[Dict[str, object]]:
    if ref_acc is None:
        return None
    batch = sample_inext_bootstrap_reference_batch(ref_acc, N_BOOT_MATRIX)
    if batch is None:
        return None
    q1_matrix = calc_inext_q1_matrix(batch, planned_m_values, endpoint_mask)
    point_matrix = calc_inext_q1_matrix(ref_acc.unsqueeze(0), planned_m_values, endpoint_mask)
    sc_matrix = calc_sc_chao_inext_matrix(batch, planned_m_values, endpoint_mask)
    return {
        "q1_inext_td_m_est": q1_matrix.detach().cpu().numpy().astype(float),
        "point": {
            "q1_inext_td_m_est": point_matrix.detach().cpu().numpy().astype(float).reshape(-1),
        },
        "sc_model": sc_matrix.detach().cpu().numpy().astype(float),
        "sc_emp": sc_matrix.detach().cpu().numpy().astype(float),
        "f1": (batch == 1).sum(dim=1).detach().cpu().numpy().astype(float),
        "n": batch.sum(dim=1).detach().cpu().numpy().astype(float),
    }


def build_stage2_q2_fast_payload(
    ref_acc: Optional[torch.Tensor],
    planned_m_values: Sequence[float],
    endpoint_mask: Sequence[bool],
) -> Optional[Dict[str, object]]:
    if ref_acc is None:
        return None
    batch = sample_inext_bootstrap_reference_batch(ref_acc, N_BOOT_MATRIX)
    if batch is None:
        return None
    q2_matrix = calc_inext_q2_matrix(batch, planned_m_values, endpoint_mask)
    point_matrix = calc_inext_q2_matrix(ref_acc.unsqueeze(0), planned_m_values, endpoint_mask)
    sc_matrix = calc_sc_chao_inext_matrix(batch, planned_m_values, endpoint_mask)
    return {
        "q2_inext_td_m_est": q2_matrix.detach().cpu().numpy().astype(float),
        "point": {
            "q2_inext_td_m_est": point_matrix.detach().cpu().numpy().astype(float).reshape(-1),
        },
        "sc_model": sc_matrix.detach().cpu().numpy().astype(float),
        "sc_emp": sc_matrix.detach().cpu().numpy().astype(float),
        "f1": (batch == 1).sum(dim=1).detach().cpu().numpy().astype(float),
        "n": batch.sum(dim=1).detach().cpu().numpy().astype(float),
    }


def calc_sc_chao_inext_at_size_batch(count_batch: torch.Tensor, size: Optional[int], is_endpoint: bool) -> torch.Tensor:
    out: List[float] = []
    for row in count_batch:
        x = row[row > 0]
        if x.numel() == 0:
            out.append(float("nan"))
            continue
        if is_endpoint:
            out.append(1.0)
            continue
        n = int(round(float(x.sum().item())))
        m = n if size is None else max(1, int(round(size)))
        out.append(float(calc_sc_chao_inext(x, m)))
    return torch.as_tensor(out, dtype=SOLVE_DTYPE, device=count_batch.device)


def sample_multinomial_counts(prob: torch.Tensor, size: int) -> torch.Tensor:
    size = int(round(size))
    if size <= 0:
        return torch.zeros_like(prob, dtype=torch.int64)
    draws = torch.multinomial(prob, num_samples=size, replacement=True)
    return torch.bincount(draws, minlength=prob.numel()).to(torch.int64)


def sample_multinomial_counts_batch(prob_batch: torch.Tensor, size_vec: torch.Tensor) -> torch.Tensor:
    size_vec = size_vec.to(torch.int64)
    out = torch.zeros(prob_batch.shape, dtype=torch.int64, device=prob_batch.device)
    if prob_batch.numel() == 0:
        return out
    valid = (size_vec > 0) & (prob_batch.sum(dim=1) > 0)
    if not bool(valid.any().item()):
        return out
    probs_valid = prob_batch[valid]
    sizes_valid = size_vec[valid]
    max_size = int(sizes_valid.max().item())
    draws = torch.multinomial(probs_valid, num_samples=max_size, replacement=True)
    step_ix = torch.arange(max_size, device=prob_batch.device, dtype=torch.int64).unsqueeze(0)
    weights = (step_ix < sizes_valid.unsqueeze(1)).to(torch.int64)
    sampled = torch.zeros((probs_valid.shape[0], prob_batch.shape[1]), dtype=torch.int64, device=prob_batch.device)
    sampled.scatter_add_(1, draws, weights)
    out[valid] = sampled
    return out


def estimate_chao_unseen_tail(counts: torch.Tensor) -> Tuple[float, float, float]:
    x = counts[counts > 0].to(SOLVE_DTYPE)
    if x.numel() == 0:
        return 0.0, 1.0, 0.0
    n = float(x.sum().item())
    if n <= 0.0 or not math.isfinite(n):
        return 0.0, 1.0, 0.0
    f1 = float((x == 1).sum().item())
    f2 = float((x == 2).sum().item())
    if f1 <= 0.0:
        return 0.0, 1.0, 0.0
    if f2 > 0.0:
        f0_hat = ((n - 1.0) / n) * (f1 ** 2) / (2.0 * f2)
    else:
        f0_hat = ((n - 1.0) / n) * f1 * max(f1 - 1.0, 0.0) / 2.0
    if f0_hat <= 0.0 or not math.isfinite(f0_hat):
        return 0.0, 1.0, 0.0
    a_hat = (n * f0_hat) / (n * f0_hat + f1)
    unseen_mass = (f1 / n) * a_hat
    unseen_mass = min(max(unseen_mass, 0.0), 0.95)
    return float(f0_hat), float(a_hat), float(unseen_mass)


def make_inext_bootstrap_probs(counts: torch.Tensor) -> torch.Tensor:
    x = counts[counts > 0].to(DTYPE)
    if x.numel() == 0:
        return x
    n, _f1, _f2, _a_hat, f0_hat, unseen_mass = chao_f0_a_from_counts(x)
    total = torch.clamp(x.sum(), min=1.0)
    p = x / total
    b = (p * torch.pow(1.0 - p, torch.tensor(float(n), dtype=DTYPE, device=counts.device))).sum()
    if f0_hat <= 0.0 or float(b.item()) <= 0.0:
        w = 0.0
    else:
        w = unseen_mass / float(b.item())
    obs_probs = p * (1.0 - w * torch.pow(1.0 - p, torch.tensor(float(n), dtype=DTYPE, device=counts.device)))
    obs_probs = torch.clamp(obs_probs, min=0.0)
    if unseen_mass <= 0.0:
        probs = obs_probs
    else:
        unseen_k = max(math.ceil(f0_hat), INEXT_UNSEEN_MIN_CATEGORIES)
        if INEXT_UNSEEN_CAP > 0:
            unseen_k = min(unseen_k, INEXT_UNSEEN_CAP)
        unseen_probs = torch.full((unseen_k,), unseen_mass / unseen_k, dtype=DTYPE, device=counts.device)
        probs = torch.cat([obs_probs, unseen_probs], dim=0)
    probs = torch.clamp(probs, min=0.0)
    prob_sum = probs.sum()
    if float(prob_sum.item()) <= 0.0:
        return x / total
    return probs / prob_sum


def sample_inext_bootstrap_reference_batch(ref_counts: torch.Tensor, num_boot: int) -> Optional[torch.Tensor]:
    counts = ref_counts[ref_counts > 0].to(torch.int64)
    if counts.numel() == 0:
        return None
    total = int(counts.sum().item())
    if total <= 0:
        return None
    probs = make_inext_bootstrap_probs(counts)
    draws = torch.multinomial(probs.expand(num_boot, -1), num_samples=total, replacement=True)
    out = torch.zeros((num_boot, probs.numel()), device=DEVICE, dtype=torch.int64)
    out.scatter_add_(1, draws, torch.ones_like(draws, dtype=torch.int64))
    return out


def sample_stage2_batch(ref_counts: torch.Tensor, size: Optional[int], num_boot: int, is_asymptotic: bool) -> Optional[torch.Tensor]:
    counts = ref_counts[ref_counts > 0].to(torch.int64)
    if counts.numel() == 0:
        return None
    total = int(counts.sum().item())
    if is_asymptotic:
        return counts.unsqueeze(0).repeat(num_boot, 1)
    if size is None:
        return None
    size = int(round(size))
    if size <= 0 or total <= 0:
        return None
    if size <= total:
        population = torch.repeat_interleave(torch.arange(counts.numel(), device=DEVICE, dtype=torch.int64), counts)
        scores = torch.rand((num_boot, population.numel()), device=DEVICE, dtype=DTYPE)
        top_idx = torch.topk(scores, k=size, dim=1).indices
        picked = population[top_idx]
        out = torch.zeros((num_boot, counts.numel()), device=DEVICE, dtype=torch.int64)
        out.scatter_add_(1, picked, torch.ones_like(picked, dtype=torch.int64))
        return out

    if STAGE2_EXTRAPOLATION_METHOD in {"chao_unseen_bootstrap", "inext_bootstrap", "chao_inext_bootstrap", "inext_table2_bootstrap"}:
        probs = make_inext_bootstrap_probs(counts)
    else:
        probs = counts.to(DTYPE)
        probs = probs / probs.sum()
    draws = torch.multinomial(probs.expand(num_boot, -1), num_samples=size, replacement=True)
    out = torch.zeros((num_boot, probs.numel()), device=DEVICE, dtype=torch.int64)
    out.scatter_add_(1, draws, torch.ones_like(draws, dtype=torch.int64))
    return out


def infer_missing_reason(value: float, n_valid_conv: float, n_valid_syneco: float, cell_status: str, metric_name: str) -> Optional[str]:
    if value is not None and math.isfinite(value):
        return None
    if cell_status == "reference_empty_after_mapping":
        return cell_status
    if cell_status == "no_valid_boot_replicates":
        return cell_status
    if cell_status == "insufficient_group_replicates":
        return cell_status
    if metric_name == "NegLogP" and isinstance(cell_status, str) and cell_status.startswith("bm_"):
        return "bm_pvalue_undefined"
    if metric_name == "WelchNegLogP" and isinstance(cell_status, str) and cell_status.startswith("welch_"):
        return "welch_pvalue_undefined"
    if not math.isfinite(n_valid_conv) or not math.isfinite(n_valid_syneco):
        return "not_evaluated"
    return "metric_missing_unexpected"


def calc_bm_stats(x: np.ndarray, y: np.ndarray) -> Dict[str, float | str]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]

    if x.size < 2 or y.size < 2:
        return {
            "p": np.nan,
            "A": np.nan,
            "CliffsDelta": np.nan,
            "BM_Statistic": np.nan,
            "BM_DF": np.nan,
            "status": "insufficient_group_replicates",
        }

    n1 = x.size
    n2 = y.size
    combined = np.concatenate([x, y])
    r = _rankdata_average_1d(combined)
    r1 = _rankdata_average_1d(x)
    r2 = _rankdata_average_1d(y)
    idx1 = np.arange(n1)
    idx2 = n1 + np.arange(n2)
    m1 = float(np.mean(r[idx1]))
    m2 = float(np.mean(r[idx2]))
    a_val = (m2 - (n2 + 1.0) / 2.0) / n1
    a_val = min(max(a_val, 0.0), 1.0)
    cliff = 2.0 * a_val - 1.0

    v1 = np.sum((r[idx1] - r1 - m1 + (n1 + 1.0) / 2.0) ** 2) / (n1 - 1)
    v2 = np.sum((r[idx2] - r2 - m2 + (n2 + 1.0) / 2.0) ** 2) / (n2 - 1)
    denom = math.sqrt(n1 * v1 + n2 * v2) if (n1 * v1 + n2 * v2) > 0 else 0.0

    if not math.isfinite(denom) or denom <= 1e-12:
        p_val = 1.0 if abs(cliff) < 1e-12 else np.nan
        status = "ok" if abs(cliff) < 1e-12 else "bm_undefined_zero_denominator"
        return {
            "p": p_val,
            "A": a_val,
            "CliffsDelta": cliff,
            "BM_Statistic": 0.0 if abs(cliff) < 1e-12 else np.nan,
            "BM_DF": np.nan,
            "status": status,
        }

    statistic = n1 * n2 * (m2 - m1) / (n1 + n2) / denom
    df_num = (n1 * v1 + n2 * v2) ** 2
    df_den = ((n1 * v1) ** 2) / (n1 - 1) + ((n2 * v2) ** 2) / (n2 - 1)
    if not math.isfinite(df_den) or df_den <= 0:
        return {
            "p": np.nan,
            "A": a_val,
            "CliffsDelta": cliff,
            "BM_Statistic": statistic,
            "BM_DF": np.nan,
            "status": "bm_invalid_df",
        }

    df_bm = df_num / df_den
    p_val = _student_t_two_sided_p(statistic, df_bm)
    return {
        "p": p_val,
        "A": a_val,
        "CliffsDelta": cliff,
        "BM_Statistic": statistic,
        "BM_DF": df_bm,
        "status": "ok",
    }


def calc_welch_stats(x: np.ndarray, y: np.ndarray) -> Dict[str, float | str]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]

    if x.size < 2 or y.size < 2:
        return {
            "p": np.nan,
            "Welch_Statistic": np.nan,
            "Welch_DF": np.nan,
            "status": "insufficient_group_replicates",
        }

    n1 = x.size
    n2 = y.size
    mean1 = float(np.mean(x))
    mean2 = float(np.mean(y))
    var1 = float(np.var(x, ddof=1))
    var2 = float(np.var(y, ddof=1))
    denom_sq = var1 / n1 + var2 / n2
    denom = math.sqrt(denom_sq) if denom_sq > 0 else 0.0

    if not math.isfinite(denom) or denom <= 1e-12:
        if abs(mean2 - mean1) < 1e-12:
            return {
                "p": 1.0,
                "Welch_Statistic": 0.0,
                "Welch_DF": np.nan,
                "status": "ok",
            }
        return {
            "p": np.nan,
            "Welch_Statistic": np.nan,
            "Welch_DF": np.nan,
            "status": "welch_undefined_zero_denominator",
        }

    statistic = (mean2 - mean1) / denom
    df_num = denom_sq ** 2
    df_den = ((var1 / n1) ** 2) / (n1 - 1) + ((var2 / n2) ** 2) / (n2 - 1)
    if not math.isfinite(df_den) or df_den <= 0:
        return {
            "p": np.nan,
            "Welch_Statistic": statistic,
            "Welch_DF": np.nan,
            "status": "welch_invalid_df",
        }

    df_welch = df_num / df_den
    p_val = _student_t_two_sided_p(statistic, df_welch)
    return {
        "p": p_val,
        "Welch_Statistic": statistic,
        "Welch_DF": df_welch,
        "status": "ok",
    }


def classify_sweep_status(r_conv: float, r_syneco: float, is_asymptotic: bool = False) -> Dict[str, object]:
    if is_asymptotic:
        return {
            "Conv_Uses_Extrapolation": np.nan,
            "Syneco_Uses_Extrapolation": np.nan,
            "Sweep_Status": "both_asymptotic",
            "Primary_Eligible": False,
        }
    conv_ex = np.nan if not math.isfinite(r_conv) else (r_conv > 1.0)
    syn_ex = np.nan if not math.isfinite(r_syneco) else (r_syneco > 1.0)
    if pd.isna(conv_ex) or pd.isna(syn_ex):
        status = "unknown"
    elif conv_ex and syn_ex:
        status = "both_extrapolation"
    elif conv_ex or syn_ex:
        status = "mixed"
    else:
        status = "both_rarefaction"
    primary = (
        math.isfinite(r_conv)
        and math.isfinite(r_syneco)
        and max(r_conv, r_syneco) <= PRIMARY_MAX_M_OVER_N
    )
    return {
        "Conv_Uses_Extrapolation": conv_ex,
        "Syneco_Uses_Extrapolation": syn_ex,
        "Sweep_Status": status,
        "Primary_Eligible": bool(primary),
    }


def build_binary_matrix_from_pairs(
    pairs_df: pd.DataFrame,
    row_col: str,
    value_col: str,
    row_index: pd.Index,
) -> Optional[pd.DataFrame]:
    if pairs_df.empty:
        return None
    pairs_df = pairs_df[[row_col, value_col]].dropna().drop_duplicates()
    if pairs_df.empty:
        return None
    mat = (
        pairs_df.assign(v=1.0)
        .pivot_table(index=row_col, columns=value_col, values="v", fill_value=0.0, aggfunc="max")
        .reindex(row_index, fill_value=0.0)
    )
    mat.index = row_index
    return mat.astype(float)


KEGG_DRUG_PREFIX_RE = re.compile(r"^(?:D\d{5}|DG\d+)\b")
ATC_CODE_RE = re.compile(r"^([A-Z](?:\d{2}(?:[A-Z]{1,2}(?:\d{2})?)?)?)(?=\s|$)")


def extract_atc_code(name: object) -> Optional[str]:
    text = str(name).strip()
    if KEGG_DRUG_PREFIX_RE.match(text):
        return None
    match = ATC_CODE_RE.match(text)
    if not match:
        return None
    code = match.group(1)
    return code if code else None


def flatten_atc(node: Dict[str, object], current_code: Optional[str] = None) -> List[Dict[str, str]]:
    name = str(node.get("name", ""))
    my_code = extract_atc_code(name)
    children = node.get("children")
    if not children:
        did_match = re.search(r"D\d{5}", name)
        if did_match and current_code:
            return [{"DID": did_match.group(0), "ATC": current_code}]
        return []
    out: List[Dict[str, str]] = []
    next_code = my_code or current_code
    for child in children:
        out.extend(flatten_atc(child, next_code))
    return out


def main() -> None:
    set_global_seed(BASE_SEED)
    log(f"[CONFIG] ROOT_DIR={ROOT_DIR}")
    log(f"[CONFIG] DATA_DIR={DATA_DIR}")
    log(f"[CONFIG] DEVICE={DEVICE}")
    log(f"[CONFIG] Output Directory: {OUT_DIR}")
    log(f"[CONFIG] WORKER_MODE={WORKER_MODE}, WORKER_GPU_ID={WORKER_GPU_ID}")
    log(f"[CONFIG] ACTIVE_CUTOFF_RANGE={ACTIVE_CUTOFF_START}-{ACTIVE_CUTOFF_END} / {len(FULL_CUTOFF_PLAN)}")
    log(
        f"[CONFIG] POC_MODE={POC_MODE}, N_OUTER_REP={N_OUTER_REP}, "
        f"N_BOOT_MATRIX={N_BOOT_MATRIX}, ENABLE_WELCH={int(ENABLE_WELCH)}"
    )
    log(f"[CONFIG] REFERENCE_TARGET_SC={REFERENCE_TARGET_SC:.3f}, SWEEP_SC_LEVELS={len(SWEEP_SC_LEVELS)} values")
    log(
        f"[CONFIG] STAGE2_EXTRAPOLATION_METHOD={STAGE2_EXTRAPOLATION_METHOD}, "
        f"INEXT_UNSEEN_CAP={INEXT_UNSEEN_CAP}"
    )
    log(f"[CONFIG] DOMAIN_FILTER_RAW={DOMAIN_FILTER_RAW}, DOMAIN_FILTER_EXPANDED={DOMAIN_FILTER}")
    log(f"[CONFIG] MODE_FILTER_RAW={MODE_FILTER_RAW}, MODE_FILTER_EXPANDED={MODE_FILTER}")
    log(f"[CONFIG] Q_FILTER_RAW={Q_FILTER_RAW}, ACTIVE_Q={[meta['Q_Label'] for meta in ACTIVE_METRIC_META]}")
    log(f"[CONFIG] Q0_FAST_PATH={int(Q0_FAST_PATH)}, Q1_FAST_PATH={int(Q1_FAST_PATH)}, Q2_FAST_PATH={int(Q2_FAST_PATH)}")

    mtbl = quiet_read_csv(DATA_DIR / "SupplementaryMaterial 2 Intensity Data KEGG.csv")
    btbl = quiet_read_csv(DATA_DIR / "brite250819.csv")
    atc_json_path = DATA_DIR / "kegg_br08303.json"
    if mtbl is None:
        raise RuntimeError("[FATAL] mtbl is NULL. Check data path/files.")

    normalized_names = {col: re.sub(r"[^A-Za-z0-9]", "", col).lower() for col in mtbl.columns}
    high_targets = [
        "High-EI 2014",
        "High-EI 2015",
        "High-EI 2016",
        "High-EI 2017",
        "High-EI 2018",
        "High-EI 2019",
        "High-EI 2020",
    ]
    low_targets = [
        "Low-EI 2015",
        "Low-EI 2016",
        "Low-EI 2017",
        "Low-EI 2018",
        "Low-EI 2019",
        "Low-EI 2020",
    ]
    high_norm = {re.sub(r"[^A-Za-z0-9]", "", col).lower() for col in high_targets}
    low_norm = {re.sub(r"[^A-Za-z0-9]", "", col).lower() for col in low_targets}
    assemblages_resolved = {
        "syneco": [col for col, norm in normalized_names.items() if norm in high_norm],
        "conv": [col for col, norm in normalized_names.items() if norm in low_norm],
    }
    formula_keys = normalize_keys(mtbl["Formula"].astype(str))
    mtbl = mtbl.copy()
    mtbl["Formula_norm"] = formula_keys
    mtbl["idx"] = np.arange(1, len(mtbl) + 1)
    formula_index = pd.Index(formula_keys, name="Formula_norm")
    formula_to_pos = {formula: idx for idx, formula in enumerate(formula_index)}
    requested_exclude_formulas = normalize_formula_list(EXCLUDE_FORMULAS or [])
    exclude_formula_set = set(requested_exclude_formulas)
    matched_exclude_formulas = sorted(exclude_formula_set.intersection(set(formula_index)))
    missing_exclude_formulas = sorted(exclude_formula_set.difference(set(formula_index)))
    if requested_exclude_formulas and EXCLUDE_FORMULA_MODE not in {"zero", "none"}:
        raise RuntimeError(
            f"Unsupported TOMATO_EXCLUDE_FORMULA_MODE={EXCLUDE_FORMULA_MODE!r}. "
            "Expected zero or none."
        )
    if requested_exclude_formulas:
        log(
            "[CONFIG] LOFO exclude formulas requested="
            f"{len(requested_exclude_formulas)}, matched={len(matched_exclude_formulas)}, "
            f"missing={len(missing_exclude_formulas)}, mode={EXCLUDE_FORMULA_MODE}"
        )
        if missing_exclude_formulas:
            log(f"[WARN] Requested exclude formulas not found: {', '.join(missing_exclude_formulas)}")

    f2c_rows: List[Dict[str, object]] = []
    for idx, row in mtbl[["idx", "Formula_norm", "KEGG ID"]].iterrows():
        cids = re.findall(r"C\d{5}", str(row["KEGG ID"]))
        for cid in cids:
            f2c_rows.append({"idx": row["idx"], "Formula_norm": row["Formula_norm"], "CID": cid})
    f2c_all = pd.DataFrame(f2c_rows)

    M_brite = None
    if btbl is not None and "CID" in btbl.columns:
        br_map = btbl.copy()
        br_map["BRITE_ID"] = br_map["BRITE"].astype(str).str.extract(r"(?i)(br\d{5})", expand=False)
        br_map = br_map[["CID", "BRITE_ID"]].dropna().drop_duplicates()
        br_mapped = (
            f2c_all.merge(br_map, on="CID", how="inner")[["Formula_norm", "BRITE_ID"]]
            .drop_duplicates()
        )
        M_brite = build_binary_matrix_from_pairs(br_mapped, "Formula_norm", "BRITE_ID", formula_index)

    M_atc_L1 = M_atc_L2 = M_atc_L3 = None
    if atc_json_path.exists() and btbl is not None:
        cid_did_map = btbl.copy()
        cid_did_map["DID"] = cid_did_map["BRITE"].astype(str).str.extract(r"(D\d{5})", expand=False)
        cid_did_map = cid_did_map[["CID", "DID"]].dropna().drop_duplicates()
        with atc_json_path.open("r", encoding="utf-8") as handle:
            json_data = json.load(handle)
        atc_table = pd.DataFrame(flatten_atc(json_data))
        if not atc_table.empty and {"DID", "ATC"}.issubset(atc_table.columns):
            atc_mapped = f2c_all.merge(cid_did_map, on="CID", how="inner").merge(atc_table, on="DID", how="inner")
            atc_mapped["L1"] = atc_mapped["ATC"].astype(str).str.slice(0, 1)
            atc_mapped["L2"] = atc_mapped["ATC"].astype(str).str.slice(0, 3)
            atc_mapped["L3"] = atc_mapped["ATC"].astype(str).str.slice(0, 4)
            M_atc_L1 = build_binary_matrix_from_pairs(atc_mapped, "Formula_norm", "L1", formula_index)
            M_atc_L2 = build_binary_matrix_from_pairs(atc_mapped, "Formula_norm", "L2", formula_index)
            M_atc_L3 = build_binary_matrix_from_pairs(atc_mapped, "Formula_norm", "L3", formula_index)
        else:
            log("[WARN] ATC table did not expose expected DID/ATC columns; skipping ATC domains.")

    M_path = read_form_matrix(DATA_DIR / "form_path250821.csv")
    M_nt = read_form_matrix(DATA_DIR / "form_ne250821.csv")
    M_ds_ne = read_form_matrix(DATA_DIR / "form_ds_by_ne250821.csv")
    M_ds_icd11 = read_form_matrix(DATA_DIR / "form_ds_icd11cl250821.csv")
    M_ds_pathcl = read_form_matrix(DATA_DIR / "form_ds_pathcl250821.csv")
    M_pm_def = read_form_matrix(DATA_DIR / "form_prmmtbl251022.csv")
    M_sm_def = read_form_matrix(DATA_DIR / "form_scdmtbl251022.csv")

    all_f_norm = list(formula_index)
    pm_flag = pd.Series(0, index=formula_index, dtype=int)
    sm_flag = pd.Series(0, index=formula_index, dtype=int)
    if M_pm_def is not None:
        pm_names = [idx for idx, value in M_pm_def.sum(axis=1).items() if value > 0]
        pm_flag.loc[pm_flag.index.intersection(pm_names)] = 1
    if M_sm_def is not None:
        sm_names = [idx for idx, value in M_sm_def.sum(axis=1).items() if value > 0]
        sm_flag.loc[sm_flag.index.intersection(sm_names)] = 1

    def align_domain_matrix(df: Optional[pd.DataFrame]) -> Optional[torch.Tensor]:
        if df is None:
            return None
        aligned = df.groupby(level=0).sum().reindex(formula_index, fill_value=0.0)
        return torch.as_tensor(aligned.to_numpy(dtype=np.float64), dtype=DTYPE, device=DEVICE)

    domain_entries = [
        {"matrix": None, "slug": "Formula"},
        {"matrix": align_domain_matrix(M_brite), "slug": "Brite"},
        {"matrix": align_domain_matrix(M_path), "slug": "Pathway"},
        {"matrix": align_domain_matrix(M_atc_L1), "slug": "ATC_L1"},
        {"matrix": align_domain_matrix(M_atc_L2), "slug": "ATC_L2"},
        {"matrix": align_domain_matrix(M_atc_L3), "slug": "ATC_L3"},
        {"matrix": align_domain_matrix(M_nt), "slug": "Network"},
        {"matrix": align_domain_matrix(M_ds_ne), "slug": "Disease_NE"},
        {"matrix": align_domain_matrix(M_ds_icd11), "slug": "Disease_ICD11"},
        {"matrix": align_domain_matrix(M_ds_pathcl), "slug": "Disease_PathCL"},
    ]
    domain_entries = [entry for entry in domain_entries if entry["matrix"] is not None or entry["slug"] == "Formula"]
    if DOMAIN_FILTER is not None:
        domain_entries = [entry for entry in domain_entries if entry["slug"] in DOMAIN_FILTER]
    if not domain_entries:
        raise RuntimeError("No domains remain after filtering.")

    tag_specs: List[Dict[str, str]] = []
    for grp_name, cols in assemblages_resolved.items():
        for col in cols:
            year_match = re.search(r"20\d{2}", col)
            if year_match is None:
                continue
            tag_specs.append(
                {
                    "tag": f"{grp_name} {year_match.group(0)}",
                    "group": grp_name,
                    "year": year_match.group(0),
                    "col": col,
                }
            )
    if not tag_specs:
        raise RuntimeError("No assemblage tags could be resolved.")
    all_tags = [spec["tag"] for spec in tag_specs]
    all_cols_raw = [spec["col"] for spec in tag_specs]

    x_raw_np = (
        mtbl[all_cols_raw]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0.0)
        .clip(lower=0.0)
        .to_numpy(dtype=np.float64)
    )
    if matched_exclude_formulas and EXCLUDE_FORMULA_MODE == "zero":
        exclude_mask_np = np.array([formula in exclude_formula_set for formula in formula_index], dtype=bool)
        x_raw_np[exclude_mask_np, :] = 0.0
    x_raw_np = apply_intensity_transform(x_raw_np, INTENSITY_TRANSFORM)
    log(f"[CONFIG] INTENSITY_TRANSFORM={INTENSITY_TRANSFORM}")
    x_raw = torch.as_tensor(x_raw_np, dtype=DTYPE, device=DEVICE)
    tag_raw_matrix = x_raw.transpose(0, 1).contiguous()

    subset_masks = {
        "All": torch.ones(len(formula_index), dtype=torch.int64, device=DEVICE),
        "PM": torch.as_tensor(pm_flag.to_numpy(dtype=np.int64), dtype=torch.int64, device=DEVICE),
        "SM": torch.as_tensor(sm_flag.to_numpy(dtype=np.int64), dtype=torch.int64, device=DEVICE),
    }

    subsets = ["All", "PM", "SM"]
    modes = [str(year) for year in range(2015, 2021)] + list(PERIOD_MODE_YEARS.keys()) + ["Combo6yr"]
    if SUBSET_FILTER is not None:
        subsets = [sb for sb in subsets if sb in SUBSET_FILTER]
    if MODE_FILTER is not None:
        modes = [md for md in modes if md in MODE_FILTER]
    if not subsets:
        raise RuntimeError("No subsets remain after filtering.")
    if not modes:
        raise RuntimeError("No modes remain after filtering.")

    mode_group_weights: Dict[str, torch.Tensor] = {}
    mode_tag_manifest_rows: List[Dict[str, object]] = []
    for md in modes:
        weights = torch.zeros((2, len(all_tags)), dtype=DTYPE, device=DEVICE)
        for idx, spec in enumerate(tag_specs):
            if not mode_includes_year(md, spec["year"]):
                continue
            row = 0 if spec["group"] == "conv" else 1
            weights[row, idx] = 1.0
            mode_tag_manifest_rows.append(
                {
                    "Mode": md,
                    "Tag": spec["tag"],
                    "Group": spec["group"],
                    "Year": spec["year"],
                    "Column": spec["col"],
                }
            )
        mode_group_weights[md] = weights
    pd.DataFrame(mode_tag_manifest_rows).to_csv(OUT_DIR / "mode_tag_manifest_v36_9.csv", index=False)

    log(f"[CONFIG] Precomputing signed-cutoff probability cube for {len(all_tags)} tags ...")
    cutoff_prob_cube = torch.stack(
        [build_signed_cutoff_prob_stack(tag_raw_matrix[idx]) for idx in range(len(all_tags))],
        dim=1,
    ).contiguous()

    log("[CONFIG] Solving stage-1 reference sizes for all cutoff x tag cells ...")
    reference_planned_m_mat = torch.zeros((len(CUTOFF_PLAN), len(all_tags)), dtype=torch.int64, device=DEVICE)
    for cut_idx in range(len(CUTOFF_PLAN)):
        for tag_idx in range(len(all_tags)):
            p = cutoff_prob_cube[cut_idx, tag_idx]
            if float(p.sum().item()) <= 0:
                reference_planned_m_mat[cut_idx, tag_idx] = 0
            else:
                reference_planned_m_mat[cut_idx, tag_idx] = max(
                    1,
                    solve_m_expected(p, REFERENCE_TARGET_SC, m_max=N_UPPER_CAP),
                )

    all_results_outer: List[Dict[str, object]] = []
    all_group_outer: List[Dict[str, object]] = []

    finite_sweep_levels = [value for value in SWEEP_SC_LEVELS if abs(value - 1.0) >= 1e-9]

    for cut_i, cut_def in enumerate(CUTOFF_PLAN.to_dict(orient="records"), start=1):
        global_cut_i = ACTIVE_CUTOFF_START + cut_i - 1
        log(
            f"[Cutoff {global_cut_i}/{len(FULL_CUTOFF_PLAN)} | shard {cut_i}/{len(CUTOFF_PLAN)}] {cut_def['Cutoff_Key']} "
            f"({cut_def['Cutoff_Side']} {cut_def['Retain_Pct']}%)"
        )

        cut_idx0 = cut_i - 1
        prob_matrix = cutoff_prob_cube[cut_idx0]
        planned_m_vec = reference_planned_m_mat[cut_idx0]
        if float(prob_matrix.sum().item()) <= 0:
            continue

        for outer in range(1, N_OUTER_REP + 1):
            log(f"  [Outer {outer}/{N_OUTER_REP}]")
            seed = BASE_SEED + outer + cut_i * 1000
            set_global_seed(seed)

            reference_counts_matrix = sample_multinomial_counts_batch(prob_matrix, planned_m_vec)

            for dom in domain_entries:
                dom_slug = str(dom["slug"])
                dom_matrix = dom["matrix"]
                for sb in subsets:
                    subset_mask = subset_masks[sb]
                    filtered_matrix = reference_counts_matrix * subset_mask.unsqueeze(0)
                    if dom_slug == "Formula":
                        mapped_matrix = filtered_matrix
                    else:
                        mapped_matrix = torch.matmul(filtered_matrix.to(DTYPE), dom_matrix).round().to(torch.int64)

                    for md in modes:
                        weights = mode_group_weights[md]
                        aggregated = torch.matmul(weights, mapped_matrix.to(DTYPE)).round().to(torch.int64)
                        ref_acc_conv = aggregated[0]
                        ref_acc_syneco = aggregated[1]
                        if int(ref_acc_conv.sum().item()) <= 0:
                            ref_acc_conv = None
                        if int(ref_acc_syneco.sum().item()) <= 0:
                            ref_acc_syneco = None

                        reference_empty = ref_acc_conv is None or ref_acc_syneco is None
                        ref_sc_conv_model = calc_sc_chao_inext(ref_acc_conv) if ref_acc_conv is not None else float("nan")
                        ref_sc_syneco_model = calc_sc_chao_inext(ref_acc_syneco) if ref_acc_syneco is not None else float("nan")
                        ref_sc_conv_emp = (
                            float(calc_sc_empirical_batch(ref_acc_conv.unsqueeze(0)).item())
                            if ref_acc_conv is not None
                            else float("nan")
                        )
                        ref_sc_syneco_emp = (
                            float(calc_sc_empirical_batch(ref_acc_syneco.unsqueeze(0)).item())
                            if ref_acc_syneco is not None
                            else float("nan")
                        )
                        ref_n_conv = float(ref_acc_conv.sum().item()) if ref_acc_conv is not None else float("nan")
                        ref_n_syneco = float(ref_acc_syneco.sum().item()) if ref_acc_syneco is not None else float("nan")

                        finite_planned_conv = (
                            solve_m_chao_inext_grid(ref_acc_conv, finite_sweep_levels, m_max=N_UPPER_CAP)
                            if (not reference_empty and ref_acc_conv is not None and math.isfinite(ref_n_conv) and ref_n_conv > 0)
                            else [None] * len(finite_sweep_levels)
                        )
                        finite_planned_syneco = (
                            solve_m_chao_inext_grid(ref_acc_syneco, finite_sweep_levels, m_max=N_UPPER_CAP)
                            if (not reference_empty and ref_acc_syneco is not None and math.isfinite(ref_n_syneco) and ref_n_syneco > 0)
                            else [None] * len(finite_sweep_levels)
                        )
                        finite_plan_map_conv = {level: finite_planned_conv[idx] for idx, level in enumerate(finite_sweep_levels)}
                        finite_plan_map_syneco = {level: finite_planned_syneco[idx] for idx, level in enumerate(finite_sweep_levels)}
                        conv_boot_cache: Dict[tuple, Optional[Dict[str, np.ndarray]]] = {}
                        syneco_boot_cache: Dict[tuple, Optional[Dict[str, np.ndarray]]] = {}

                        def get_stage2_payload(
                            ref_acc: Optional[torch.Tensor],
                            planned_m: float,
                            sweep_is_asymptotic: bool,
                            cache: Dict[tuple, Optional[Dict[str, np.ndarray]]],
                        ) -> Optional[Dict[str, np.ndarray]]:
                            if ref_acc is None:
                                return None
                            cache_key = ("asym",) if sweep_is_asymptotic else ("finite", int(round(planned_m)))
                            if cache_key in cache:
                                return cache[cache_key]
                            planned_m_int = None if not math.isfinite(planned_m) else int(round(planned_m))
                            batch = sample_inext_bootstrap_reference_batch(ref_acc, N_BOOT_MATRIX)
                            if batch is None:
                                cache[cache_key] = None
                                return None
                            if NEED_Q1_INEXT or NEED_Q2_INEXT:
                                metrics = calc_inext_td_metrics_batch(batch, planned_m_int, sweep_is_asymptotic)
                                point_metrics = calc_inext_td_metrics_batch(ref_acc.unsqueeze(0), planned_m_int, sweep_is_asymptotic)
                            else:
                                metrics = {}
                                point_metrics = {}
                            sc_model = calc_sc_chao_inext_at_size_batch(batch, planned_m_int, sweep_is_asymptotic)
                            payload: Dict[str, object] = {
                                "point": {},
                                "sc_model": sc_model.detach().cpu().numpy().astype(float),
                                "sc_emp": sc_model.detach().cpu().numpy().astype(float),
                                "f1": (batch == 1).sum(dim=1).detach().cpu().numpy().astype(float),
                                "n": batch.sum(dim=1).detach().cpu().numpy().astype(float),
                            }
                            point_payload = payload["point"]
                            assert isinstance(point_payload, dict)
                            if NEED_Q0_STAGE2:
                                q0_metrics = calc_q0_stage2_matrix(batch, [planned_m], [sweep_is_asymptotic])
                                q0_point_metrics = calc_q0_stage2_matrix(ref_acc.unsqueeze(0), [planned_m], [sweep_is_asymptotic])
                                for key in ("q0_observed_standardized", "q0_chao1_asymptotic"):
                                    payload[key] = q0_metrics[key][:, 0].detach().cpu().numpy().astype(float)
                                    point_payload[key] = float(q0_point_metrics[key][0, 0].detach().cpu().item())
                            for key, values in metrics.items():
                                payload[key] = values.detach().cpu().numpy().astype(float)
                            for key, values in point_metrics.items():
                                point_payload[key] = float(values.detach().cpu().numpy()[0])
                            cache[cache_key] = payload
                            return payload

                        sweep_endpoint_mask = [abs(value - 1.0) < 1e-9 for value in SWEEP_SC_LEVELS]
                        fast_planned_m_conv: List[float] = []
                        fast_planned_m_syneco: List[float] = []
                        for fast_sweep_sc, fast_is_endpoint in zip(SWEEP_SC_LEVELS, sweep_endpoint_mask):
                            if reference_empty or fast_is_endpoint or not math.isfinite(ref_n_conv) or ref_n_conv <= 0:
                                fast_planned_m_conv.append(np.nan)
                            else:
                                fast_planned_m_conv.append(float(max(1, int(finite_plan_map_conv[fast_sweep_sc]))))
                            if reference_empty or fast_is_endpoint or not math.isfinite(ref_n_syneco) or ref_n_syneco <= 0:
                                fast_planned_m_syneco.append(np.nan)
                            else:
                                fast_planned_m_syneco.append(float(max(1, int(finite_plan_map_syneco[fast_sweep_sc]))))

                        fast_conv_payload = (
                            build_stage2_q2_fast_payload(ref_acc_conv, fast_planned_m_conv, sweep_endpoint_mask)
                            if Q2_FAST_PATH and not reference_empty
                            else None
                        )
                        fast_syneco_payload = (
                            build_stage2_q2_fast_payload(ref_acc_syneco, fast_planned_m_syneco, sweep_endpoint_mask)
                            if Q2_FAST_PATH and not reference_empty
                            else None
                        )
                        if Q1_FAST_PATH and not reference_empty:
                            fast_conv_payload = build_stage2_q1_fast_payload(ref_acc_conv, fast_planned_m_conv, sweep_endpoint_mask)
                            fast_syneco_payload = build_stage2_q1_fast_payload(ref_acc_syneco, fast_planned_m_syneco, sweep_endpoint_mask)
                        if Q0_FAST_PATH and not reference_empty:
                            fast_conv_payload = build_stage2_q0_fast_payload(ref_acc_conv, fast_planned_m_conv, sweep_endpoint_mask)
                            fast_syneco_payload = build_stage2_q0_fast_payload(ref_acc_syneco, fast_planned_m_syneco, sweep_endpoint_mask)

                        for sweep_idx, sweep_sc in enumerate(SWEEP_SC_LEVELS):
                            sweep_is_asymptotic = abs(sweep_sc - 1.0) < 1e-9
                            sweep_label = format_sc_label(sweep_sc)
                            sweep_pct = int(round(sweep_sc * 100))

                            if reference_empty or sweep_is_asymptotic or not math.isfinite(ref_n_conv) or ref_n_conv <= 0:
                                planned_m_conv = np.nan
                            else:
                                planned_m_conv = float(max(1, int(finite_plan_map_conv[sweep_sc])))
                            if reference_empty or sweep_is_asymptotic or not math.isfinite(ref_n_syneco) or ref_n_syneco <= 0:
                                planned_m_syneco = np.nan
                            else:
                                planned_m_syneco = float(max(1, int(finite_plan_map_syneco[sweep_sc])))

                            ratio_conv = (
                                np.nan
                                if sweep_is_asymptotic or not math.isfinite(ref_n_conv) or ref_n_conv <= 0 or not math.isfinite(planned_m_conv)
                                else float(planned_m_conv / ref_n_conv)
                            )
                            ratio_syneco = (
                                np.nan
                                if sweep_is_asymptotic or not math.isfinite(ref_n_syneco) or ref_n_syneco <= 0 or not math.isfinite(planned_m_syneco)
                                else float(planned_m_syneco / ref_n_syneco)
                            )
                            sweep_info = classify_sweep_status(ratio_conv, ratio_syneco, is_asymptotic=sweep_is_asymptotic)

                            boot_metric_c = {meta["Metric_Key"]: np.array([], dtype=float) for meta in ACTIVE_METRIC_META}
                            boot_metric_s = {meta["Metric_Key"]: np.array([], dtype=float) for meta in ACTIVE_METRIC_META}
                            boot_sc_model_c = np.array([], dtype=float)
                            boot_sc_model_s = np.array([], dtype=float)
                            boot_sc_emp_c = np.array([], dtype=float)
                            boot_sc_emp_s = np.array([], dtype=float)
                            boot_f1_c = np.array([], dtype=float)
                            boot_f1_s = np.array([], dtype=float)
                            boot_n_c = np.array([], dtype=float)
                            boot_n_s = np.array([], dtype=float)
                            point_metric_c = {meta["Metric_Key"]: np.nan for meta in ACTIVE_METRIC_META}
                            point_metric_s = {meta["Metric_Key"]: np.nan for meta in ACTIVE_METRIC_META}

                            if not reference_empty:
                                if Q1_FAST_PATH:
                                    if fast_conv_payload is not None:
                                        q1_matrix_c = np.asarray(fast_conv_payload["q1_inext_td_m_est"], dtype=float)
                                        boot_metric_c["q1_inext_td_m_est"] = q1_matrix_c[:, sweep_idx]
                                        point_metric_c["q1_inext_td_m_est"] = float(
                                            np.asarray(fast_conv_payload["point"]["q1_inext_td_m_est"], dtype=float)[sweep_idx]
                                        )
                                        boot_sc_model_c = np.asarray(fast_conv_payload["sc_model"], dtype=float)[:, sweep_idx]
                                        boot_sc_emp_c = np.asarray(fast_conv_payload["sc_emp"], dtype=float)[:, sweep_idx]
                                        boot_f1_c = np.asarray(fast_conv_payload["f1"], dtype=float)
                                        boot_n_c = np.asarray(fast_conv_payload["n"], dtype=float)

                                    if fast_syneco_payload is not None:
                                        q1_matrix_s = np.asarray(fast_syneco_payload["q1_inext_td_m_est"], dtype=float)
                                        boot_metric_s["q1_inext_td_m_est"] = q1_matrix_s[:, sweep_idx]
                                        point_metric_s["q1_inext_td_m_est"] = float(
                                            np.asarray(fast_syneco_payload["point"]["q1_inext_td_m_est"], dtype=float)[sweep_idx]
                                        )
                                        boot_sc_model_s = np.asarray(fast_syneco_payload["sc_model"], dtype=float)[:, sweep_idx]
                                        boot_sc_emp_s = np.asarray(fast_syneco_payload["sc_emp"], dtype=float)[:, sweep_idx]
                                        boot_f1_s = np.asarray(fast_syneco_payload["f1"], dtype=float)
                                        boot_n_s = np.asarray(fast_syneco_payload["n"], dtype=float)
                                elif Q2_FAST_PATH:
                                    if fast_conv_payload is not None:
                                        q2_matrix_c = np.asarray(fast_conv_payload["q2_inext_td_m_est"], dtype=float)
                                        boot_metric_c["q2_inext_td_m_est"] = q2_matrix_c[:, sweep_idx]
                                        point_metric_c["q2_inext_td_m_est"] = float(
                                            np.asarray(fast_conv_payload["point"]["q2_inext_td_m_est"], dtype=float)[sweep_idx]
                                        )
                                        boot_sc_model_c = np.asarray(fast_conv_payload["sc_model"], dtype=float)[:, sweep_idx]
                                        boot_sc_emp_c = np.asarray(fast_conv_payload["sc_emp"], dtype=float)[:, sweep_idx]
                                        boot_f1_c = np.asarray(fast_conv_payload["f1"], dtype=float)
                                        boot_n_c = np.asarray(fast_conv_payload["n"], dtype=float)

                                    if fast_syneco_payload is not None:
                                        q2_matrix_s = np.asarray(fast_syneco_payload["q2_inext_td_m_est"], dtype=float)
                                        boot_metric_s["q2_inext_td_m_est"] = q2_matrix_s[:, sweep_idx]
                                        point_metric_s["q2_inext_td_m_est"] = float(
                                            np.asarray(fast_syneco_payload["point"]["q2_inext_td_m_est"], dtype=float)[sweep_idx]
                                        )
                                        boot_sc_model_s = np.asarray(fast_syneco_payload["sc_model"], dtype=float)[:, sweep_idx]
                                        boot_sc_emp_s = np.asarray(fast_syneco_payload["sc_emp"], dtype=float)[:, sweep_idx]
                                        boot_f1_s = np.asarray(fast_syneco_payload["f1"], dtype=float)
                                        boot_n_s = np.asarray(fast_syneco_payload["n"], dtype=float)
                                elif Q0_FAST_PATH:
                                    if fast_conv_payload is not None:
                                        for key in ("q0_observed_standardized", "q0_chao1_asymptotic"):
                                            q0_matrix_c = np.asarray(fast_conv_payload[key], dtype=float)
                                            boot_metric_c[key] = q0_matrix_c[:, sweep_idx]
                                            point_metric_c[key] = float(
                                                np.asarray(fast_conv_payload["point"][key], dtype=float)[sweep_idx]
                                            )
                                        boot_sc_model_c = np.asarray(fast_conv_payload["sc_model"], dtype=float)[:, sweep_idx]
                                        boot_sc_emp_c = np.asarray(fast_conv_payload["sc_emp"], dtype=float)[:, sweep_idx]
                                        boot_f1_c = np.asarray(fast_conv_payload["f1"], dtype=float)
                                        boot_n_c = np.asarray(fast_conv_payload["n"], dtype=float)

                                    if fast_syneco_payload is not None:
                                        for key in ("q0_observed_standardized", "q0_chao1_asymptotic"):
                                            q0_matrix_s = np.asarray(fast_syneco_payload[key], dtype=float)
                                            boot_metric_s[key] = q0_matrix_s[:, sweep_idx]
                                            point_metric_s[key] = float(
                                                np.asarray(fast_syneco_payload["point"][key], dtype=float)[sweep_idx]
                                            )
                                        boot_sc_model_s = np.asarray(fast_syneco_payload["sc_model"], dtype=float)[:, sweep_idx]
                                        boot_sc_emp_s = np.asarray(fast_syneco_payload["sc_emp"], dtype=float)[:, sweep_idx]
                                        boot_f1_s = np.asarray(fast_syneco_payload["f1"], dtype=float)
                                        boot_n_s = np.asarray(fast_syneco_payload["n"], dtype=float)
                                else:
                                    conv_payload = get_stage2_payload(ref_acc_conv, planned_m_conv, sweep_is_asymptotic, conv_boot_cache)
                                    syneco_payload = get_stage2_payload(ref_acc_syneco, planned_m_syneco, sweep_is_asymptotic, syneco_boot_cache)

                                    if conv_payload is not None:
                                        for key in boot_metric_c:
                                            boot_metric_c[key] = conv_payload[key]
                                            point_metric_c[key] = conv_payload["point"].get(key, np.nan)
                                        boot_sc_model_c = conv_payload["sc_model"]
                                        boot_sc_emp_c = conv_payload["sc_emp"]
                                        boot_f1_c = conv_payload["f1"]
                                        boot_n_c = conv_payload["n"]

                                    if syneco_payload is not None:
                                        for key in boot_metric_s:
                                            boot_metric_s[key] = syneco_payload[key]
                                            point_metric_s[key] = syneco_payload["point"].get(key, np.nan)
                                        boot_sc_model_s = syneco_payload["sc_model"]
                                        boot_sc_emp_s = syneco_payload["sc_emp"]
                                        boot_f1_s = syneco_payload["f1"]
                                        boot_n_s = syneco_payload["n"]

                            group_common = {
                                "OuterRep": outer,
                                "Cutoff_Key": cut_def["Cutoff_Key"],
                                "Cutoff_Side": cut_def["Cutoff_Side"],
                                "Retain_Pct": int(cut_def["Retain_Pct"]),
                                "Retain_Ratio": float(cut_def["Retain_Ratio"]),
                                "Cutoff_Signed_Pct": float(cut_def["Cutoff_Signed_Pct"]),
                                "Cutoff_Label": cut_def["Cutoff_Label"],
                                "Reference_Target_SC": REFERENCE_TARGET_SC,
                                "Sweep_SC": float(sweep_sc),
                                "Sweep_SC_Label": sweep_label,
                                "Sweep_SC_Pct": sweep_pct,
                                "Sweep_Is_Asymptotic": bool(sweep_is_asymptotic),
                                "Domain": dom_slug,
                                "Subset": sb,
                                "Mode": md,
                                "Sweep_Status": sweep_info["Sweep_Status"],
                                "Primary_Eligible": sweep_info["Primary_Eligible"],
                            }
                            all_group_outer.append(
                                {
                                    **group_common,
                                    "Group": "Conv",
                                    "Reference_SC_Model": ref_sc_conv_model,
                                    "Reference_SC_Empirical": ref_sc_conv_emp,
                                    "Reference_N": ref_n_conv,
                                    "Sweep_Planned_M": planned_m_conv,
                                    "Sweep_M_over_N": ratio_conv,
                                    "Uses_Extrapolation": sweep_info["Conv_Uses_Extrapolation"],
                                    "Realized_SC_Model": float(np.nanmean(boot_sc_model_c)) if boot_sc_model_c.size else np.nan,
                                    "Realized_SC_Empirical": float(np.nanmean(boot_sc_emp_c)) if boot_sc_emp_c.size else np.nan,
                                    "Mean_f1": float(np.nanmean(boot_f1_c)) if boot_f1_c.size else np.nan,
                                    "Mean_Realized_N": float(np.nanmean(boot_n_c)) if boot_n_c.size else np.nan,
                                }
                            )
                            all_group_outer.append(
                                {
                                    **group_common,
                                    "Group": "Syneco",
                                    "Reference_SC_Model": ref_sc_syneco_model,
                                    "Reference_SC_Empirical": ref_sc_syneco_emp,
                                    "Reference_N": ref_n_syneco,
                                    "Sweep_Planned_M": planned_m_syneco,
                                    "Sweep_M_over_N": ratio_syneco,
                                    "Uses_Extrapolation": sweep_info["Syneco_Uses_Extrapolation"],
                                    "Realized_SC_Model": float(np.nanmean(boot_sc_model_s)) if boot_sc_model_s.size else np.nan,
                                    "Realized_SC_Empirical": float(np.nanmean(boot_sc_emp_s)) if boot_sc_emp_s.size else np.nan,
                                    "Mean_f1": float(np.nanmean(boot_f1_s)) if boot_f1_s.size else np.nan,
                                    "Mean_Realized_N": float(np.nanmean(boot_n_s)) if boot_n_s.size else np.nan,
                                }
                            )

                            for meta in ACTIVE_METRIC_META:
                                metric_key = meta["Metric_Key"]
                                x = boot_metric_c[metric_key]
                                y = boot_metric_s[metric_key]
                                n_valid_conv = int(np.sum(np.isfinite(x)))
                                n_valid_syneco = int(np.sum(np.isfinite(y)))

                                if reference_empty:
                                    st = {
                                        "p": np.nan,
                                        "A": np.nan,
                                        "CliffsDelta": np.nan,
                                        "BM_Statistic": np.nan,
                                        "BM_DF": np.nan,
                                        "status": "reference_empty_after_mapping",
                                    }
                                elif n_valid_conv == 0 and n_valid_syneco == 0:
                                    st = {
                                        "p": np.nan,
                                        "A": np.nan,
                                        "CliffsDelta": np.nan,
                                        "BM_Statistic": np.nan,
                                        "BM_DF": np.nan,
                                        "status": "no_valid_boot_replicates",
                                    }
                                else:
                                    st = calc_bm_stats(x, y)

                                if reference_empty:
                                    welch_st = {
                                        "p": np.nan,
                                        "Welch_Statistic": np.nan,
                                        "Welch_DF": np.nan,
                                        "status": "reference_empty_after_mapping",
                                    }
                                elif n_valid_conv == 0 and n_valid_syneco == 0:
                                    welch_st = {
                                        "p": np.nan,
                                        "Welch_Statistic": np.nan,
                                        "Welch_DF": np.nan,
                                        "status": "no_valid_boot_replicates",
                                    }
                                elif ENABLE_WELCH:
                                    welch_st = calc_welch_stats(x, y)
                                else:
                                    welch_st = {
                                        "p": np.nan,
                                        "Welch_Statistic": np.nan,
                                        "Welch_DF": np.nan,
                                        "status": "welch_disabled",
                                    }

                                p_raw = st["p"]
                                p_val = clamp_p(float(p_raw)) if pd.notna(p_raw) else np.nan
                                neglogp = safe_neglog10p(p_val) if math.isfinite(p_val) else np.nan
                                welch_p_raw = welch_st["p"]
                                welch_p = clamp_p(float(welch_p_raw)) if pd.notna(welch_p_raw) else np.nan
                                welch_neglogp = safe_neglog10p(welch_p) if math.isfinite(welch_p) else np.nan
                                all_results_outer.append(
                                    {
                                        "OuterRep": outer,
                                        "Cutoff_Key": cut_def["Cutoff_Key"],
                                        "Cutoff_Side": cut_def["Cutoff_Side"],
                                        "Retain_Pct": int(cut_def["Retain_Pct"]),
                                        "Retain_Ratio": float(cut_def["Retain_Ratio"]),
                                        "Cutoff_Signed_Pct": float(cut_def["Cutoff_Signed_Pct"]),
                                        "Cutoff_Label": cut_def["Cutoff_Label"],
                                        "Reference_Target_SC": REFERENCE_TARGET_SC,
                                        "Sweep_SC": float(sweep_sc),
                                        "Sweep_SC_Label": sweep_label,
                                        "Sweep_SC_Pct": sweep_pct,
                                        "Sweep_Is_Asymptotic": bool(sweep_is_asymptotic),
                                        "Domain": dom_slug,
                                        "Subset": sb,
                                        "Mode": md,
                                        "Diversity_Order": int(meta["Diversity_Order"]),
                                        "Q_Label": meta["Q_Label"],
                                        "Estimate_Definition": meta["Estimate_Definition"],
                                        "Reference_SC_Conv_Model": ref_sc_conv_model,
                                        "Reference_SC_Syneco_Model": ref_sc_syneco_model,
                                        "Reference_SC_Conv_Empirical": ref_sc_conv_emp,
                                        "Reference_SC_Syneco_Empirical": ref_sc_syneco_emp,
                                        "Reference_N_Conv": ref_n_conv,
                                        "Reference_N_Syneco": ref_n_syneco,
                                        "Sweep_Planned_M_Conv": planned_m_conv,
                                        "Sweep_Planned_M_Syneco": planned_m_syneco,
                                        "Sweep_M_over_N_Conv": ratio_conv,
                                        "Sweep_M_over_N_Syneco": ratio_syneco,
                                        "Max_Sweep_M_over_N": max_or_na([ratio_conv, ratio_syneco]),
                                        "Conv_Uses_Extrapolation": sweep_info["Conv_Uses_Extrapolation"],
                                        "Syneco_Uses_Extrapolation": sweep_info["Syneco_Uses_Extrapolation"],
                                        "Sweep_Status": sweep_info["Sweep_Status"],
                                        "Primary_Eligible": sweep_info["Primary_Eligible"],
                                        "p_raw": p_raw,
                                        "p": p_val,
                                        "NegLogP": neglogp,
                                        "Welch_p_raw": welch_p_raw,
                                        "Welch_p": welch_p,
                                        "WelchNegLogP": welch_neglogp,
                                        "A": st["A"],
                                        "CliffsDelta": st["CliffsDelta"],
                                        "BM_Statistic": st["BM_Statistic"],
                                        "BM_DF": st["BM_DF"],
                                        "Cell_Status": st["status"],
                                        "Welch_Statistic": welch_st["Welch_Statistic"],
                                        "Welch_DF": welch_st["Welch_DF"],
                                        "Welch_Status": welch_st["status"],
                                        "N_Valid_Conv": n_valid_conv,
                                        "N_Valid_Syneco": n_valid_syneco,
                                        "Point_Conv": point_metric_c.get(metric_key, np.nan),
                                        "Point_Syneco": point_metric_s.get(metric_key, np.nan),
                                        "Point_Diff_SynecoMinusConv": (
                                            point_metric_s.get(metric_key, np.nan) - point_metric_c.get(metric_key, np.nan)
                                            if math.isfinite(point_metric_c.get(metric_key, np.nan)) and math.isfinite(point_metric_s.get(metric_key, np.nan))
                                            else np.nan
                                        ),
                                        "SD_Conv": float(np.std(x[np.isfinite(x)], ddof=1)) if n_valid_conv >= 2 else np.nan,
                                        "SD_Syneco": float(np.std(y[np.isfinite(y)], ddof=1)) if n_valid_syneco >= 2 else np.nan,
                                        "Mean_Conv": float(np.mean(x[np.isfinite(x)])) if n_valid_conv >= 1 else np.nan,
                                        "Mean_Syneco": float(np.mean(y[np.isfinite(y)])) if n_valid_syneco >= 1 else np.nan,
                                    }
                                )

    results_outer_df = pd.DataFrame(all_results_outer)
    group_outer_df = pd.DataFrame(all_group_outer)
    if results_outer_df.empty:
        raise RuntimeError("No comparison rows were generated.")

    comparison_group_cols = [
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
    ]

    comparison_agg_rows: List[Dict[str, object]] = []
    for keys, grp in results_outer_df.groupby(comparison_group_cols, dropna=False, sort=False):
        row = dict(zip(comparison_group_cols, keys if isinstance(keys, tuple) else (keys,)))
        row.update(
            {
                "Reference_SC_Conv_Model": median_or_na(grp["Reference_SC_Conv_Model"]),
                "Reference_SC_Syneco_Model": median_or_na(grp["Reference_SC_Syneco_Model"]),
                "Reference_SC_Conv_Empirical": median_or_na(grp["Reference_SC_Conv_Empirical"]),
                "Reference_SC_Syneco_Empirical": median_or_na(grp["Reference_SC_Syneco_Empirical"]),
                "Reference_N_Conv": median_or_na(grp["Reference_N_Conv"]),
                "Reference_N_Syneco": median_or_na(grp["Reference_N_Syneco"]),
                "Sweep_Planned_M_Conv": median_or_na(grp["Sweep_Planned_M_Conv"]),
                "Sweep_Planned_M_Syneco": median_or_na(grp["Sweep_Planned_M_Syneco"]),
                "Sweep_M_over_N_Conv": median_or_na(grp["Sweep_M_over_N_Conv"]),
                "Sweep_M_over_N_Syneco": median_or_na(grp["Sweep_M_over_N_Syneco"]),
                "Max_Sweep_M_over_N": median_or_na(grp["Max_Sweep_M_over_N"]),
                "Conv_Uses_Extrapolation": pick_first_non_na(grp["Conv_Uses_Extrapolation"]),
                "Syneco_Uses_Extrapolation": pick_first_non_na(grp["Syneco_Uses_Extrapolation"]),
                "Sweep_Status": pick_first_non_na(grp["Sweep_Status"]),
                "Primary_Eligible": pick_first_non_na(grp["Primary_Eligible"]),
                "Prob_BM_Sig": float(((grp["p"].astype(float) < ALPHA_P) & grp["p"].notna()).mean()),
                "NegLogP_Median": median_or_na(grp["NegLogP"]),
                "NegLogP_Q025": quantile_or_na(grp["NegLogP"], 0.025),
                "NegLogP_Q975": quantile_or_na(grp["NegLogP"], 0.975),
                "Prob_Welch_Sig": float(((grp["Welch_p"].astype(float) < ALPHA_P) & grp["Welch_p"].notna()).mean()),
                "WelchNegLogP_Median": median_or_na(grp["WelchNegLogP"]),
                "WelchNegLogP_Q025": quantile_or_na(grp["WelchNegLogP"], 0.025),
                "WelchNegLogP_Q975": quantile_or_na(grp["WelchNegLogP"], 0.975),
                "A_Median": median_or_na(grp["A"]),
                "A_Q025": quantile_or_na(grp["A"], 0.025),
                "A_Q975": quantile_or_na(grp["A"], 0.975),
                "Cliff_Median": median_or_na(grp["CliffsDelta"]),
                "Cliff_Q025": quantile_or_na(grp["CliffsDelta"], 0.025),
                "Cliff_Q975": quantile_or_na(grp["CliffsDelta"], 0.975),
                "Cell_Status": pick_status(grp["Cell_Status"]),
                "Welch_Status": pick_status(grp["Welch_Status"]),
                "N_Valid_Conv": median_or_na(grp["N_Valid_Conv"]),
                "N_Valid_Syneco": median_or_na(grp["N_Valid_Syneco"]),
                "Point_Conv": median_or_na(grp["Point_Conv"]),
                "Point_Syneco": median_or_na(grp["Point_Syneco"]),
                "Point_Diff_SynecoMinusConv": median_or_na(grp["Point_Diff_SynecoMinusConv"]),
                "SD_Conv": median_or_na(grp["SD_Conv"]),
                "SD_Syneco": median_or_na(grp["SD_Syneco"]),
                "Mean_Conv": median_or_na(grp["Mean_Conv"]),
                "Mean_Syneco": median_or_na(grp["Mean_Syneco"]),
            }
        )
        comparison_agg_rows.append(row)
    comparison_agg_df = pd.DataFrame(comparison_agg_rows)

    def build_cells_long(agg_df: pd.DataFrame, metric_col: str, metric_name: str, status_col: str = "Cell_Status") -> pd.DataFrame:
        df = agg_df.copy()
        df["Missing_Reason"] = [
            infer_missing_reason(v, c, s, st, metric_name)
            for v, c, s, st in zip(df[metric_col], df["N_Valid_Conv"], df["N_Valid_Syneco"], df[status_col])
        ]
        out_cols = [
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
            "Reference_SC_Conv_Model",
            "Reference_SC_Syneco_Model",
            "Reference_SC_Conv_Empirical",
            "Reference_SC_Syneco_Empirical",
            "Reference_N_Conv",
            "Reference_N_Syneco",
            "Sweep_Planned_M_Conv",
            "Sweep_Planned_M_Syneco",
            "Sweep_M_over_N_Conv",
            "Sweep_M_over_N_Syneco",
            "Max_Sweep_M_over_N",
            "Conv_Uses_Extrapolation",
            "Syneco_Uses_Extrapolation",
            "Sweep_Status",
            "Primary_Eligible",
            "Cell_Status",
            "Missing_Reason",
            "N_Valid_Conv",
            "N_Valid_Syneco",
            "Point_Conv",
            "Point_Syneco",
            "Point_Diff_SynecoMinusConv",
            "SD_Conv",
            "SD_Syneco",
            "Mean_Conv",
            "Mean_Syneco",
        ]
        out = df[out_cols].copy()
        out["Cell_Status"] = df[status_col]
        out["Metric"] = metric_name
        out["Metric_Column"] = metric_col
        out["Value"] = df[metric_col]
        return out

    cells_long_parts = [
        build_cells_long(comparison_agg_df, "NegLogP_Median", "NegLogP"),
        build_cells_long(comparison_agg_df, "A_Median", "A"),
        build_cells_long(comparison_agg_df, "Cliff_Median", "CliffsDelta"),
    ]
    if ENABLE_WELCH:
        cells_long_parts.append(
            build_cells_long(comparison_agg_df, "WelchNegLogP_Median", "WelchNegLogP", status_col="Welch_Status")
        )
    comparison_cells_long_df = pd.concat(cells_long_parts, ignore_index=True)

    if WRITE_CELLS_ONLY:
        cutoff_manifest_path = OUT_DIR / "cutoff_manifest_v36_9.csv"
        CUTOFF_PLAN.to_csv(cutoff_manifest_path, index=False)

        sweep_manifest_df = pd.DataFrame(
            {
                "Reference_Target_SC": REFERENCE_TARGET_SC,
                "Sweep_SC": SWEEP_SC_LEVELS,
                "Sweep_SC_Label": [format_sc_label(v) for v in SWEEP_SC_LEVELS],
                "Sweep_SC_Pct": [int(round(v * 100)) for v in SWEEP_SC_LEVELS],
                "Sweep_Is_Asymptotic": [abs(v - 1.0) < 1e-9 for v in SWEEP_SC_LEVELS],
            }
        )
        sweep_manifest_df.to_csv(OUT_DIR / "sweep_sc_manifest_v36_9.csv", index=False)
        comparison_agg_df.to_csv(OUT_DIR / "comparison_agg_v36_9_q012.csv", index=False)
        comparison_cells_long_df.to_csv(OUT_DIR / "comparison_cells_long_all_q012_v36_9.csv", index=False)

        run_meta = {
            "version": "v36.9-python-4gpu-v4-stage2-inext-q2fast",
            "output_directory": str(OUT_DIR),
            "root_dir": str(ROOT_DIR),
            "data_dir": str(DATA_DIR),
            "device": str(DEVICE),
            "write_cells_only": True,
            "omitted_outputs": ["comparison_outer_v36_9_q012.csv", "group_outer_v36_9_q012.csv", "group_summary_v36_9_q012.csv"],
            "poc_mode": POC_MODE,
            "n_outer_rep": N_OUTER_REP,
            "n_boot_matrix": N_BOOT_MATRIX,
            "reference_target_sc": REFERENCE_TARGET_SC,
            "q_filter_raw": Q_FILTER_RAW,
            "active_metric_meta": ACTIVE_METRIC_META,
            "q0_fast_path": Q0_FAST_PATH,
            "q0_estimator_lineage": "2026-06-30 planned-m corrected q0 richness: finite q0_observed_standardized uses the Stage2 planned sample size; q0_chao1_asymptotic is reference-count Chao1 richness",
            "sweep_sc_levels": SWEEP_SC_LEVELS,
            "primary_max_m_over_n": PRIMARY_MAX_M_OVER_N,
            "domain_filter_raw": DOMAIN_FILTER_RAW,
            "domain_filter": DOMAIN_FILTER,
            "subset_filter": SUBSET_FILTER,
            "mode_filter_raw": MODE_FILTER_RAW,
            "mode_filter": MODE_FILTER,
            "intensity_transform": INTENSITY_TRANSFORM,
            "enable_welch": ENABLE_WELCH,
            "seed": BASE_SEED,
        }
        with (OUT_DIR / "run_metadata_v36_9.json").open("w", encoding="utf-8") as handle:
            json.dump(run_meta, handle, ensure_ascii=False, indent=2)
        run_note = [
            f"OUT_DIR: {OUT_DIR}",
            f"ROOT_DIR: {ROOT_DIR}",
            f"DEVICE: {DEVICE}",
            f"WRITE_CELLS_ONLY: {int(WRITE_CELLS_ONLY)}",
            f"N_OUTER_REP: {N_OUTER_REP}",
            f"N_BOOT_MATRIX: {N_BOOT_MATRIX}",
            f"REFERENCE_TARGET_SC: {REFERENCE_TARGET_SC}",
            f"DOMAIN_FILTER_RAW: {DOMAIN_FILTER_RAW}",
            f"DOMAIN_FILTER_EXPANDED: {DOMAIN_FILTER}",
            f"MODE_FILTER_RAW: {MODE_FILTER_RAW}",
            f"MODE_FILTER_EXPANDED: {MODE_FILTER}",
            f"Q_FILTER_RAW: {Q_FILTER_RAW}",
            f"Q0_FAST_PATH: {int(Q0_FAST_PATH)}",
            f"INTENSITY_TRANSFORM: {INTENSITY_TRANSFORM}",
            "",
            "Notes:",
            "- Lightweight cells-only output omits outer-replicate and group-level CSVs.",
            "- This mode is intended for downstream mask/overlay diagnostics that only require comparison_cells_long_all_q012_v36_9.csv.",
        ]
        (OUT_DIR / "run_note.txt").write_text("\n".join(run_note) + "\n", encoding="utf-8")
        log(f"[SUCCESS] v36.9 cells-only output written to: {OUT_DIR}")
        return

    group_group_cols = [
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
    ]
    group_summary_rows: List[Dict[str, object]] = []
    for keys, grp in group_outer_df.groupby(group_group_cols, dropna=False, sort=False):
        row = dict(zip(group_group_cols, keys if isinstance(keys, tuple) else (keys,)))
        row.update(
            {
                "Reference_SC_Model": median_or_na(grp["Reference_SC_Model"]),
                "Reference_SC_Empirical": median_or_na(grp["Reference_SC_Empirical"]),
                "Reference_N": median_or_na(grp["Reference_N"]),
                "Sweep_Planned_M": median_or_na(grp["Sweep_Planned_M"]),
                "Sweep_M_over_N": median_or_na(grp["Sweep_M_over_N"]),
                "Uses_Extrapolation": pick_first_non_na(grp["Uses_Extrapolation"]),
                "Sweep_Status": pick_first_non_na(grp["Sweep_Status"]),
                "Primary_Eligible": pick_first_non_na(grp["Primary_Eligible"]),
                "Realized_SC_Model_Median": median_or_na(grp["Realized_SC_Model"]),
                "Realized_SC_Model_Q025": quantile_or_na(grp["Realized_SC_Model"], 0.025),
                "Realized_SC_Model_Q975": quantile_or_na(grp["Realized_SC_Model"], 0.975),
                "Realized_SC_Empirical_Median": median_or_na(grp["Realized_SC_Empirical"]),
                "Mean_f1_Median": median_or_na(grp["Mean_f1"]),
                "Mean_Realized_N_Median": median_or_na(grp["Mean_Realized_N"]),
            }
        )
        group_summary_rows.append(row)
    group_summary_df = pd.DataFrame(group_summary_rows)

    cutoff_manifest_path = OUT_DIR / "cutoff_manifest_v36_9.csv"
    CUTOFF_PLAN.to_csv(cutoff_manifest_path, index=False)

    sweep_manifest_df = pd.DataFrame(
        {
            "Reference_Target_SC": REFERENCE_TARGET_SC,
            "Sweep_SC": SWEEP_SC_LEVELS,
            "Sweep_SC_Label": [format_sc_label(v) for v in SWEEP_SC_LEVELS],
            "Sweep_SC_Pct": [int(round(v * 100)) for v in SWEEP_SC_LEVELS],
            "Sweep_Is_Asymptotic": [abs(v - 1.0) < 1e-9 for v in SWEEP_SC_LEVELS],
        }
    )
    sweep_manifest_df.to_csv(OUT_DIR / "sweep_sc_manifest_v36_9.csv", index=False)

    results_outer_df.to_csv(OUT_DIR / "comparison_outer_v36_9_q012.csv", index=False)
    group_outer_df.to_csv(OUT_DIR / "group_outer_v36_9_q012.csv", index=False)
    group_summary_df.to_csv(OUT_DIR / "group_summary_v36_9_q012.csv", index=False)
    comparison_agg_df.to_csv(OUT_DIR / "comparison_agg_v36_9_q012.csv", index=False)
    comparison_cells_long_df.to_csv(OUT_DIR / "comparison_cells_long_all_q012_v36_9.csv", index=False)

    run_meta = {
        "version": "v36.9-python-4gpu-v4-stage2-inext-q2fast",
        "output_directory": str(OUT_DIR),
        "root_dir": str(ROOT_DIR),
        "data_dir": str(DATA_DIR),
        "device": str(DEVICE),
        "poc_mode": POC_MODE,
        "n_outer_rep": N_OUTER_REP,
        "n_boot_matrix": N_BOOT_MATRIX,
        "reference_target_sc": REFERENCE_TARGET_SC,
        "reference_stage1_sc_method": "model_based_SCm_on_intensity_pdf",
        "stage2_sweep_sc_method": "manual_Chao_iNEXT_abundance_coverage_no_scipy_gpu_optimized",
        "stage2_sweep_sc_grid_policy": "use TOMATO_SWEEP_SC_GRID or the default 0.01..1.00 grid as given; do not auto-add Reference_Target_SC",
            "stage2_diversity_estimator": "Table2/iNEXT TD.m.est abundance estimator for q1/q2; q0 finite rows use planned-m richness rarefaction/extrapolation; bootstrap vectors are evaluated with the same estimator family",
            "stage2_extrapolation_method": STAGE2_EXTRAPOLATION_METHOD,
            "stage2_rarefaction_method": "Table2/iNEXT interpolation estimator for q1/q2",
            "stage2_endpoint_policy": "Sweep_SC=1.000 is labeled Asymptotic and evaluated with Chao/iNEXT asymptotic q1/q2 estimators, not as an ordinary finite coverage cell",
            "stage2_boot_meaning": "iNEXT-style EstiBootComm.Ind bootstrap reference count vectors conditional on each Stage1 reference count vector; not experimental replicates",
            "q0_fast_path": Q0_FAST_PATH,
            "q0_estimator_lineage": "2026-06-30 planned-m corrected q0 richness: finite q0_observed_standardized uses the Stage2 planned sample size; q0_chao1_asymptotic is reference-count Chao1 richness",
            "q1_fast_path": Q1_FAST_PATH,
            "q1_fast_path_meaning": "q1-only runs reuse one bootstrap count batch across all Sweep_SC values and evaluate q1 Table2 interpolation with a NumPy/SciPy vectorized hypergeometric sum",
            "q2_fast_path": Q2_FAST_PATH,
        "q2_fast_path_meaning": "q2-only runs reuse one bootstrap count batch across all Sweep_SC values and evaluate q2/coverage as Boot x Sweep_SC matrices",
        "inext_unseen_cap": INEXT_UNSEEN_CAP,
        "inext_unseen_min_categories": INEXT_UNSEEN_MIN_CATEGORIES,
        "q_filter_raw": Q_FILTER_RAW,
        "active_metric_meta": ACTIVE_METRIC_META,
        "dtype": str(DTYPE),
        "solve_dtype": str(SOLVE_DTYPE),
        "sweep_sc_levels": SWEEP_SC_LEVELS,
        "primary_max_m_over_n": PRIMARY_MAX_M_OVER_N,
        "domain_filter_raw": DOMAIN_FILTER_RAW,
        "domain_filter": DOMAIN_FILTER,
        "subset_filter": SUBSET_FILTER,
        "mode_filter_raw": MODE_FILTER_RAW,
        "mode_filter": MODE_FILTER,
        "period_mode_years": {key: sorted(value) for key, value in PERIOD_MODE_YEARS.items()},
        "combo6yr_years": [str(year) for year in range(2015, 2021)],
        "combo6yr_excludes_unpaired_2014": True,
        "intensity_transform": INTENSITY_TRANSFORM,
        "exclude_formulas_requested": requested_exclude_formulas,
        "exclude_formula_mode": EXCLUDE_FORMULA_MODE,
        "exclude_formulas_matched": matched_exclude_formulas,
        "exclude_formulas_missing": missing_exclude_formulas,
        "enable_welch": ENABLE_WELCH,
        "a_definition": "P(Conv < Syneco) + 0.5 P(Conv = Syneco)",
        "cliffs_definition": "2A - 1",
        "seed": BASE_SEED,
    }
    with (OUT_DIR / "run_metadata_v36_9.json").open("w", encoding="utf-8") as handle:
        json.dump(run_meta, handle, ensure_ascii=False, indent=2)

    run_note = [
        f"OUT_DIR: {OUT_DIR}",
        f"ROOT_DIR: {ROOT_DIR}",
        f"DATA_DIR: {DATA_DIR}",
        f"DEVICE: {DEVICE}",
        f"POC_MODE: {POC_MODE}",
        f"N_OUTER_REP: {N_OUTER_REP}",
        f"N_BOOT_MATRIX: {N_BOOT_MATRIX}",
        f"ENABLE_WELCH: {int(ENABLE_WELCH)}",
        f"REFERENCE_TARGET_SC: {REFERENCE_TARGET_SC}",
        "REFERENCE_STAGE1_SC_METHOD: model_based_SCm_on_intensity_pdf",
        "STAGE2_SWEEP_SC_METHOD: manual_Chao_iNEXT_abundance_coverage_no_scipy_gpu_optimized",
        f"STAGE2_EXTRAPOLATION_METHOD: {STAGE2_EXTRAPOLATION_METHOD}",
        f"INEXT_UNSEEN_CAP: {INEXT_UNSEEN_CAP}",
        f"DOMAIN_FILTER_RAW: {DOMAIN_FILTER_RAW}",
        f"DOMAIN_FILTER_EXPANDED: {DOMAIN_FILTER}",
        f"MODE_FILTER_RAW: {MODE_FILTER_RAW}",
        f"MODE_FILTER_EXPANDED: {MODE_FILTER}",
        f"Q_FILTER_RAW: {Q_FILTER_RAW}",
        f"ACTIVE_Q: {', '.join(str(meta['Q_Label']) for meta in ACTIVE_METRIC_META)}",
            f"Q0_FAST_PATH: {int(Q0_FAST_PATH)}",
            f"Q1_FAST_PATH: {int(Q1_FAST_PATH)}",
            f"Q2_FAST_PATH: {int(Q2_FAST_PATH)}",
        f"DTYPE: {DTYPE}",
        f"SOLVE_DTYPE: {SOLVE_DTYPE}",
        f"SWEEP_SC_LEVELS: {', '.join(format_sc_label(v) for v in SWEEP_SC_LEVELS)}",
        f"PRIMARY_MAX_M_OVER_N: {PRIMARY_MAX_M_OVER_N}",
        f"INTENSITY_TRANSFORM: {INTENSITY_TRANSFORM}",
        f"EXCLUDE_FORMULAS_REQUESTED: {', '.join(requested_exclude_formulas) if requested_exclude_formulas else '(none)'}",
        f"EXCLUDE_FORMULAS_MATCHED: {', '.join(matched_exclude_formulas) if matched_exclude_formulas else '(none)'}",
        f"EXCLUDE_FORMULA_MODE: {EXCLUDE_FORMULA_MODE}",
        "COMBO6YR_YEARS: 2015, 2016, 2017, 2018, 2019, 2020",
        "PERIOD_MODES: Period1_2015_2017=2015,2016,2017; Period2_2018_2020=2018,2019,2020",
        "",
        "Notes:",
        "- Signed cutoff design = low 1..99, none, high 1..99.",
        "- Combo6yr excludes unpaired High-EI 2014 and uses only paired 2015-2020 tags.",
        "- Stage 1 uses fixed reference target SC to generate one pseudo-sample per tag.",
        "- Stage 1 still solves m_ref with model-based SC(m) on the intensity-derived probability vector.",
        "- Stage 2 sweeps sample coverage from that reference sample.",
        "- Stage 2 Sweep_SC grid is independent of Reference_Target_SC and is not auto-augmented with the Stage 1 target.",
            "- Stage 2 planned_m and coverage QC use a manual Chao/iNEXT-type abundance coverage estimator.",
            "- Finite q0_observed_standardized uses the planned Stage2 sample size m; q0_chao1_asymptotic remains a reference-count Chao1 richness estimate.",
            "- Finite rarefaction q1/q2 sweeps use the Table 2 / iNEXT interpolation estimator.",
        "- Finite extrapolation q1/q2 sweeps use the Table 2 / iNEXT extrapolation estimator.",
        "- Bootstrap draws follow an EstiBootComm.Ind-style unseen-tail assemblage and are evaluated with the same estimator.",
        "- Sweep_SC = 1.000 is labeled Asymptotic and evaluated with Chao/iNEXT asymptotic q1/q2 estimators.",
        "- Stage 2 Boot values are Monte Carlo bootstrap count vectors, not experimental replicates.",
        "- A is defined so larger values mean Syneco tends to exceed Conv.",
        "- Cliff's delta is derived as 2A - 1.",
        "- NegLogP is -log10(BM p).",
        (
            "- WelchNegLogP is -log10(Welch p) and is exported as a diagnostic-only companion metric."
            if ENABLE_WELCH
            else "- Welch diagnostic output is disabled for this run."
        ),
        "- Primary_Eligible means max(finite sweep m / reference n) <= 2.",
    ]
    (OUT_DIR / "run_note.txt").write_text("\n".join(run_note) + "\n", encoding="utf-8")
    log(f"[SUCCESS] v36.9 output written to: {OUT_DIR}")


def resolve_gpu_ids() -> List[int]:
    if os.getenv("TOMATO_FORCE_CPU", "0") == "1" or not torch.cuda.is_available():
        return []
    visible = torch.cuda.device_count()
    parsed: List[int] = []
    if GPU_ID_RAW is not None:
        for item in GPU_ID_RAW:
            try:
                gpu_id = int(item)
            except (TypeError, ValueError):
                continue
            if 0 <= gpu_id < visible:
                parsed.append(gpu_id)
    if not parsed:
        parsed = list(range(visible))
    if len(parsed) < REQUESTED_NUM_WORKERS:
        base = list(parsed)
        while len(parsed) < REQUESTED_NUM_WORKERS:
            parsed.extend(base)
    return parsed[:REQUESTED_NUM_WORKERS]


def split_contiguous_ranges(total: int, num_parts: int) -> List[Tuple[int, int]]:
    num_parts = max(1, min(num_parts, total))
    base = total // num_parts
    rem = total % num_parts
    out: List[Tuple[int, int]] = []
    start = 1
    for idx in range(num_parts):
        size = base + (1 if idx < rem else 0)
        end = start + size - 1
        out.append((start, end))
        start = end + 1
    return out


def split_fixed_size_ranges(total: int, chunk_size: int) -> List[Tuple[int, int]]:
    chunk_size = max(1, chunk_size)
    out: List[Tuple[int, int]] = []
    start = 1
    while start <= total:
        end = min(total, start + chunk_size - 1)
        out.append((start, end))
        start = end + 1
    return out


def launch_worker(
    worker_idx: int,
    gpu_id: int,
    cutoff_start: int,
    cutoff_end: int,
    shard_dir: Path,
) -> Tuple[subprocess.Popen, object, Path]:
    shard_dir.mkdir(parents=True, exist_ok=True)
    log_path = shard_dir / "worker.log"
    env = os.environ.copy()
    env.update(
        {
            "TOMATO_4GPU_WORKER": "1",
            "TOMATO_4GPU_GPU_ID": str(gpu_id),
            "TOMATO_4GPU_CUTOFF_START": str(cutoff_start),
            "TOMATO_4GPU_CUTOFF_END": str(cutoff_end),
            "TOMATO_4GPU_NUM_WORKERS": str(REQUESTED_NUM_WORKERS),
            "TOMATO_OUT_DIR": str(shard_dir),
            # Keep cutoff-local seeds identical to the single-GPU schedule.
            "TOMATO_BASE_SEED": str(BASE_SEED + (cutoff_start - 1) * 1000),
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
        f"[ORCH] worker={worker_idx} gpu={gpu_id} "
        f"cutoffs={cutoff_start}-{cutoff_end} log={log_path}"
    )
    return proc, handle, log_path


def stream_concat_shard_csvs(shard_dirs: Sequence[Path], filename: str, out_path: Path) -> int:
    expected_header: Optional[str] = None
    row_count = 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    found = False
    with out_path.open("w", encoding="utf-8", newline="") as out_handle:
        for shard_dir in shard_dirs:
            path = shard_dir / filename
            if not path.exists():
                raise RuntimeError(f"Missing shard output: {path}")
            found = True
            with path.open("r", encoding="utf-8", newline="") as in_handle:
                header = in_handle.readline()
                if not header:
                    raise RuntimeError(f"Empty shard output: {path}")
                if expected_header is None:
                    expected_header = header
                    out_handle.write(header)
                elif header.rstrip("\r\n") != expected_header.rstrip("\r\n"):
                    raise RuntimeError(f"Header mismatch while merging {filename}: {path}")
                while True:
                    chunk = in_handle.read(1024 * 1024)
                    if not chunk:
                        break
                    row_count += chunk.count("\n")
                    out_handle.write(chunk)
    if not found:
        raise RuntimeError(f"No shard outputs found for {filename}")
    return row_count

def write_final_outputs(shard_dirs: Sequence[Path], gpu_ids: Sequence[int], cutoff_ranges: Sequence[Tuple[int, int]]) -> None:
    merge_row_counts = {
        "comparison_outer_v36_9_q012.csv": stream_concat_shard_csvs(
            shard_dirs,
            "comparison_outer_v36_9_q012.csv",
            OUT_DIR / "comparison_outer_v36_9_q012.csv",
        ),
        "group_outer_v36_9_q012.csv": stream_concat_shard_csvs(
            shard_dirs,
            "group_outer_v36_9_q012.csv",
            OUT_DIR / "group_outer_v36_9_q012.csv",
        ),
        "group_summary_v36_9_q012.csv": stream_concat_shard_csvs(
            shard_dirs,
            "group_summary_v36_9_q012.csv",
            OUT_DIR / "group_summary_v36_9_q012.csv",
        ),
        "comparison_agg_v36_9_q012.csv": stream_concat_shard_csvs(
            shard_dirs,
            "comparison_agg_v36_9_q012.csv",
            OUT_DIR / "comparison_agg_v36_9_q012.csv",
        ),
        "comparison_cells_long_all_q012_v36_9.csv": stream_concat_shard_csvs(
            shard_dirs,
            "comparison_cells_long_all_q012_v36_9.csv",
            OUT_DIR / "comparison_cells_long_all_q012_v36_9.csv",
        ),
    }

    FULL_CUTOFF_PLAN.to_csv(OUT_DIR / "cutoff_manifest_v36_9.csv", index=False)
    sweep_manifest_df = pd.DataFrame(
        {
            "Reference_Target_SC": REFERENCE_TARGET_SC,
            "Sweep_SC": SWEEP_SC_LEVELS,
            "Sweep_SC_Label": [format_sc_label(v) for v in SWEEP_SC_LEVELS],
            "Sweep_SC_Pct": [int(round(v * 100)) for v in SWEEP_SC_LEVELS],
            "Sweep_Is_Asymptotic": [abs(v - 1.0) < 1e-9 for v in SWEEP_SC_LEVELS],
        }
    )
    sweep_manifest_df.to_csv(OUT_DIR / "sweep_sc_manifest_v36_9.csv", index=False)
    mode_tag_manifest_parts = [
        pd.read_csv(shard_dir / "mode_tag_manifest_v36_9.csv")
        for shard_dir in shard_dirs
        if (shard_dir / "mode_tag_manifest_v36_9.csv").exists()
    ]
    if mode_tag_manifest_parts:
        (
            pd.concat(mode_tag_manifest_parts, ignore_index=True)
            .drop_duplicates()
            .sort_values(["Mode", "Group", "Year", "Tag"])
            .to_csv(OUT_DIR / "mode_tag_manifest_v36_9.csv", index=False)
        )

    run_meta = {
        "version": "v36.9-python-4gpu-v4-stage2-inext-q2fast",
        "output_directory": str(OUT_DIR),
        "root_dir": str(ROOT_DIR),
        "data_dir": str(DATA_DIR),
        "device": "multi_gpu_orchestrated",
        "gpu_ids": list(gpu_ids),
        "worker_count": len(gpu_ids),
        "cutoff_chunk_size": WORKER_CUTOFF_CHUNK_SIZE,
        "orchestration_scheduler": "dynamic_queue",
        "final_merge_method": "streaming_header_preserving_csv_concat",
        "final_merge_row_counts": merge_row_counts,
        "execution_unit_count": len(cutoff_ranges),
        "cutoff_ranges": [{"start": start, "end": end} for start, end in cutoff_ranges],
        "poc_mode": POC_MODE,
        "n_outer_rep": N_OUTER_REP,
        "n_boot_matrix": N_BOOT_MATRIX,
        "reference_target_sc": REFERENCE_TARGET_SC,
        "reference_stage1_sc_method": "model_based_SCm_on_intensity_pdf",
        "stage2_sweep_sc_method": "manual_Chao_iNEXT_abundance_coverage_no_scipy_gpu_optimized",
        "stage2_sweep_sc_grid_policy": "use TOMATO_SWEEP_SC_GRID or the default 0.01..1.00 grid as given; do not auto-add Reference_Target_SC",
        "stage2_diversity_estimator": "Table2/iNEXT TD.m.est abundance estimator for q1/q2; q0 finite rows use planned-m richness rarefaction/extrapolation; bootstrap vectors are evaluated with the same estimator family",
        "stage2_extrapolation_method": STAGE2_EXTRAPOLATION_METHOD,
        "stage2_rarefaction_method": "Table2/iNEXT interpolation estimator for q1/q2",
        "stage2_endpoint_policy": "Sweep_SC=1.000 is labeled Asymptotic and evaluated with Chao/iNEXT asymptotic q1/q2 estimators, not as an ordinary finite coverage cell",
        "stage2_boot_meaning": "iNEXT-style EstiBootComm.Ind bootstrap reference count vectors conditional on each Stage1 reference count vector; not experimental replicates",
        "q0_fast_path": Q0_FAST_PATH,
        "q0_estimator_lineage": "2026-06-30 planned-m corrected q0 richness: finite q0_observed_standardized uses the Stage2 planned sample size; q0_chao1_asymptotic is reference-count Chao1 richness",
        "q1_fast_path": Q1_FAST_PATH,
        "q1_fast_path_meaning": "q1-only runs reuse one bootstrap count batch across all Sweep_SC values and evaluate q1 Table2 interpolation with a NumPy/SciPy vectorized hypergeometric sum",
        "q2_fast_path": Q2_FAST_PATH,
        "q2_fast_path_meaning": "q2-only runs reuse one bootstrap count batch across all Sweep_SC values and evaluate q2/coverage as Boot x Sweep_SC matrices",
        "inext_unseen_cap": INEXT_UNSEEN_CAP,
        "inext_unseen_min_categories": INEXT_UNSEEN_MIN_CATEGORIES,
        "q_filter_raw": Q_FILTER_RAW,
        "active_metric_meta": ACTIVE_METRIC_META,
        "dtype": str(DTYPE),
        "solve_dtype": str(SOLVE_DTYPE),
        "sweep_sc_levels": SWEEP_SC_LEVELS,
        "primary_max_m_over_n": PRIMARY_MAX_M_OVER_N,
        "domain_filter_raw": DOMAIN_FILTER_RAW,
        "domain_filter": DOMAIN_FILTER,
        "subset_filter": SUBSET_FILTER,
        "mode_filter_raw": MODE_FILTER_RAW,
        "mode_filter": MODE_FILTER,
        "period_mode_years": {key: sorted(value) for key, value in PERIOD_MODE_YEARS.items()},
        "combo6yr_years": [str(year) for year in range(2015, 2021)],
        "combo6yr_excludes_unpaired_2014": True,
        "intensity_transform": INTENSITY_TRANSFORM,
        "exclude_formulas_requested": normalize_formula_list(EXCLUDE_FORMULAS or []),
        "exclude_formula_mode": EXCLUDE_FORMULA_MODE,
        "enable_welch": ENABLE_WELCH,
        "a_definition": "P(Conv < Syneco) + 0.5 P(Conv = Syneco)",
        "cliffs_definition": "2A - 1",
        "seed": BASE_SEED,
    }
    with (OUT_DIR / "run_metadata_v36_9.json").open("w", encoding="utf-8") as handle:
        json.dump(run_meta, handle, ensure_ascii=False, indent=2)
    with (OUT_DIR / "run_metadata_v36_9_4gpu.json").open("w", encoding="utf-8") as handle:
        json.dump(run_meta, handle, ensure_ascii=False, indent=2)

    run_note = [
        f"OUT_DIR: {OUT_DIR}",
        f"ROOT_DIR: {ROOT_DIR}",
        f"DATA_DIR: {DATA_DIR}",
        "DEVICE: multi_gpu_orchestrated",
        f"GPU_IDS: {', '.join(str(v) for v in gpu_ids)}",
        f"WORKER_COUNT: {len(gpu_ids)}",
        f"CUTOFF_CHUNK_SIZE: {WORKER_CUTOFF_CHUNK_SIZE}",
        "ORCHESTRATION_SCHEDULER: dynamic_queue",
        "FINAL_MERGE_METHOD: streaming_header_preserving_csv_concat",
        f"CUTOFF_RANGES: {', '.join(f'{start}-{end}' for start, end in cutoff_ranges)}",
        f"POC_MODE: {POC_MODE}",
        f"N_OUTER_REP: {N_OUTER_REP}",
        f"N_BOOT_MATRIX: {N_BOOT_MATRIX}",
        f"ENABLE_WELCH: {int(ENABLE_WELCH)}",
        f"REFERENCE_TARGET_SC: {REFERENCE_TARGET_SC}",
        "REFERENCE_STAGE1_SC_METHOD: model_based_SCm_on_intensity_pdf",
        "STAGE2_SWEEP_SC_METHOD: manual_Chao_iNEXT_abundance_coverage_no_scipy_gpu_optimized",
        f"STAGE2_EXTRAPOLATION_METHOD: {STAGE2_EXTRAPOLATION_METHOD}",
        f"INEXT_UNSEEN_CAP: {INEXT_UNSEEN_CAP}",
        f"DOMAIN_FILTER_RAW: {DOMAIN_FILTER_RAW}",
        f"DOMAIN_FILTER_EXPANDED: {DOMAIN_FILTER}",
        f"MODE_FILTER_RAW: {MODE_FILTER_RAW}",
        f"MODE_FILTER_EXPANDED: {MODE_FILTER}",
        f"Q_FILTER_RAW: {Q_FILTER_RAW}",
        f"ACTIVE_Q: {', '.join(str(meta['Q_Label']) for meta in ACTIVE_METRIC_META)}",
        f"Q0_FAST_PATH: {int(Q0_FAST_PATH)}",
        f"Q1_FAST_PATH: {int(Q1_FAST_PATH)}",
        f"Q2_FAST_PATH: {int(Q2_FAST_PATH)}",
        f"DTYPE: {DTYPE}",
        f"SOLVE_DTYPE: {SOLVE_DTYPE}",
        f"SWEEP_SC_LEVELS: {', '.join(format_sc_label(v) for v in SWEEP_SC_LEVELS)}",
        f"PRIMARY_MAX_M_OVER_N: {PRIMARY_MAX_M_OVER_N}",
        f"INTENSITY_TRANSFORM: {INTENSITY_TRANSFORM}",
        f"EXCLUDE_FORMULAS_REQUESTED: {', '.join(normalize_formula_list(EXCLUDE_FORMULAS or [])) if EXCLUDE_FORMULAS else '(none)'}",
        f"EXCLUDE_FORMULA_MODE: {EXCLUDE_FORMULA_MODE}",
        "COMBO6YR_YEARS: 2015, 2016, 2017, 2018, 2019, 2020",
        "PERIOD_MODES: Period1_2015_2017=2015,2016,2017; Period2_2018_2020=2018,2019,2020",
        "",
        "Notes:",
        "- 4GPU orchestration shards cutoff ranges across workers, then stream-concatenates shard outputs without pandas all-in-memory merge.",
        "- Combo6yr excludes unpaired High-EI 2014 and uses only paired 2015-2020 tags.",
        "- Stage 1 definition is unchanged from v36.9 v2.",
        "- Stage 2 Sweep_SC grid is independent of Reference_Target_SC and is not auto-augmented with the Stage 1 target.",
        "- Finite q0_observed_standardized uses the planned Stage2 sample size m; q0_chao1_asymptotic remains a reference-count Chao1 richness estimate.",
        "- Stage 2 q1/q2 finite rarefaction/extrapolation uses Table 2 / iNEXT TD.m.est-style estimators.",
        "- Bootstrap draws follow an EstiBootComm.Ind-style unseen-tail assemblage and are evaluated with the same estimator.",
        "- Sweep_SC = 1.000 is labeled Asymptotic and evaluated with Chao/iNEXT asymptotic q1/q2 estimators.",
        "- OuterRep and Boot remain Monte Carlo repetitions, not experimental replicates.",
        "- NegLogP is -log10(BM p), A = P(Conv < Syneco) + 0.5 P(Conv = Syneco), Cliff's delta = 2A - 1.",
        (
            "- WelchNegLogP is -log10(Welch p) and is exported as a diagnostic-only companion metric."
            if ENABLE_WELCH
            else "- Welch diagnostic output is disabled for this run."
        ),
    ]
    (OUT_DIR / "run_note.txt").write_text("\n".join(run_note) + "\n", encoding="utf-8")


def orchestrate_4gpu() -> None:
    gpu_ids = resolve_gpu_ids()
    if len(gpu_ids) <= 1:
        log("[ORCH] Fewer than 2 usable GPUs detected. Falling back to single-worker execution.")
        main()
        return

    if WORKER_CUTOFF_CHUNK_SIZE > 0:
        cutoff_ranges = split_fixed_size_ranges(len(FULL_CUTOFF_PLAN), WORKER_CUTOFF_CHUNK_SIZE)
        log(
            f"[ORCH] Chunked orchestration enabled: chunk_size={WORKER_CUTOFF_CHUNK_SIZE}, "
            f"execution_units={len(cutoff_ranges)}, concurrent_workers={len(gpu_ids)}"
        )
    else:
        cutoff_ranges = split_contiguous_ranges(len(FULL_CUTOFF_PLAN), len(gpu_ids))
    shard_root = OUT_DIR / "_shards"
    shard_root.mkdir(parents=True, exist_ok=True)

    launched: List[Tuple[subprocess.Popen, object, Path, Path]] = []
    try:
        next_range_idx = 0
        active: List[Tuple[subprocess.Popen, object, Path, Path, int]] = []
        log(
            f"[ORCH] Dynamic queue enabled: execution_units={len(cutoff_ranges)}, "
            f"concurrent_workers={len(gpu_ids)}"
        )

        def launch_next_on_gpu(gpu_id: int) -> Optional[Tuple[subprocess.Popen, object, Path, Path, int]]:
            nonlocal next_range_idx
            if next_range_idx >= len(cutoff_ranges):
                return None
            cutoff_start, cutoff_end = cutoff_ranges[next_range_idx]
            worker_idx = next_range_idx + 1
            shard_dir = shard_root / f"worker_{worker_idx:03d}_gpu{gpu_id}_cutoffs_{cutoff_start:03d}_{cutoff_end:03d}"
            proc, handle, log_path = launch_worker(worker_idx, gpu_id, cutoff_start, cutoff_end, shard_dir)
            launched.append((proc, handle, log_path, shard_dir))
            next_range_idx += 1
            return (proc, handle, log_path, shard_dir, gpu_id)

        for gpu_id in gpu_ids:
            record = launch_next_on_gpu(gpu_id)
            if record is not None:
                active.append(record)

        failed_logs: List[Path] = []
        finished_count = 0
        while active:
            next_active: List[Tuple[subprocess.Popen, object, Path, Path, int]] = []
            completed_this_poll = 0
            for proc, handle, log_path, shard_dir, gpu_id in active:
                return_code = proc.poll()
                if return_code is None:
                    next_active.append((proc, handle, log_path, shard_dir, gpu_id))
                    continue
                completed_this_poll += 1
                finished_count += 1
                handle.close()
                if return_code != 0:
                    failed_logs.append(log_path)
                    continue
                log(f"[ORCH] Finished worker {finished_count}/{len(cutoff_ranges)} on gpu={gpu_id}: {shard_dir.name}")
                record = launch_next_on_gpu(gpu_id)
                if record is not None:
                    next_active.append(record)
            active = next_active
            if failed_logs:
                joined = "\n".join(str(path) for path in failed_logs)
                raise RuntimeError(f"One or more 4GPU workers failed. See logs:\n{joined}")
            if active and completed_this_poll == 0:
                time.sleep(5)

        write_final_outputs([item[3] for item in launched], gpu_ids, cutoff_ranges)
        log(f"[SUCCESS] v36.9 4GPU output written to: {OUT_DIR}")
    finally:
        for proc, handle, _, _ in launched:
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
