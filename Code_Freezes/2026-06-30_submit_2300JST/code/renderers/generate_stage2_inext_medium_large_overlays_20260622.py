#!/usr/bin/env python3
"""Generate Stage2 iNEXT q2 medium/large Cliff's delta consensus overlays."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


DOMAINS = [
    "Formula",
    "Brite",
    "Pathway",
    "ATC_L1",
    "ATC_L2",
    "ATC_L3",
    "Network",
    "Disease_NE",
    "Disease_ICD11",
    "Disease_PathCL",
]
SUBSETS = ["All", "PM", "SM"]
Q_LABEL = "q2"
ESTIMATE = "Stage2_iNEXT_TD_m_est"
METRIC = "CliffsDelta"
MEDIUM_THRESHOLD = 0.33
LARGE_THRESHOLD = 0.474
PERIOD_SCENARIOS = ["Period1_2015_2017", "Period2_2018_2020"]
ANNUAL_SCENARIOS = ["2015", "2016", "2017", "2018", "2019", "2020"]
USECOLS = [
    "Domain",
    "Subset",
    "Mode",
    "Q_Label",
    "Estimate_Definition",
    "Cutoff_Side",
    "Retain_Pct",
    "Cutoff_Signed_Pct",
    "Sweep_SC",
    "Sweep_Is_Asymptotic",
    "Sweep_Status",
    "Primary_Eligible",
    "Metric",
    "Value",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--domains", default=",".join(DOMAINS))
    parser.add_argument("--subsets", default=",".join(SUBSETS))
    parser.add_argument("--q-label", default=Q_LABEL)
    parser.add_argument("--medium-threshold", type=float, default=MEDIUM_THRESHOLD)
    parser.add_argument("--large-threshold", type=float, default=LARGE_THRESHOLD)
    parser.add_argument("--threshold-kinds", default="medium,large")
    parser.add_argument("--chunksize", type=int, default=250_000)
    parser.add_argument("--skip-annual", action="store_true")
    return parser.parse_args()


def parse_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def font(size: int, bold: bool = False, italic: bool = False) -> ImageFont.ImageFont:
    candidates: List[str] = []
    if bold and italic:
        candidates.extend(
            [
                "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf",
                "/Library/Fonts/Arial Bold Italic.ttf",
            ]
        )
    if bold:
        candidates.extend(
            [
                "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
                "/Library/Fonts/Arial Bold.ttf",
            ]
        )
    if italic:
        candidates.extend(
            [
                "/System/Library/Fonts/Supplemental/Arial Italic.ttf",
                "/Library/Fonts/Arial Italic.ttf",
            ]
        )
    candidates.extend(
        [
            "/System/Library/Fonts/Supplemental/Arial.ttf",
            "/System/Library/Fonts/Supplemental/Helvetica.ttf",
            "/Library/Fonts/Arial.ttf",
        ]
    )
    for path in candidates:
        if Path(path).exists():
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


F_TITLE = font(24, True)
F_SUBTITLE = font(14)
F_AXIS = font(13)
F_AXIS_ITALIC = font(13, italic=True)
F_PANEL_TITLE = font(15)
F_TICK = font(11)
F_LEGEND = font(11)
F_FOOT = font(10)
SIGNED_CUTOFF_DISPLAY_COL = "Signed_Cutoff_Display_Pct"
SIGNED_CUTOFF_AXIS_LABEL_PREFIX = "Intensity cutoff "
SAMPLE_COVERAGE_AXIS_LABEL_PREFIX = "Sample coverage "
SIGNED_CUTOFF_DISPLAY_VALUES = list(range(-99, 100))
SIGNED_CUTOFF_TICK_TARGETS = [-95, -75, -50, -25, 0, 25, 50, 75, 95]


def split_csv(text: str) -> List[str]:
    return [item.strip() for item in text.split(",") if item.strip()]


def present_usecols(path: Path) -> List[str]:
    header = pd.read_csv(path, nrows=0).columns.tolist()
    return [col for col in USECOLS if col in header]


def cells_path(root: Path, domain: str) -> Path:
    domain_path = root / domain / "csv" / "comparison_cells_filtered.csv"
    if domain_path.exists():
        return domain_path
    combined_path = root / "csv" / "comparison_cells_filtered.csv"
    if combined_path.exists():
        return combined_path
    return domain_path


def load_cells(
    root: Path,
    domain: str,
    subsets: Sequence[str],
    modes: Sequence[str],
    chunksize: int,
    q_label: str,
) -> pd.DataFrame:
    path = cells_path(root, domain)
    if not path.exists():
        raise FileNotFoundError(path)
    frames: List[pd.DataFrame] = []
    for chunk in pd.read_csv(path, usecols=present_usecols(path), chunksize=chunksize, low_memory=False):
        mask = (
            chunk["Domain"].astype(str).eq(domain)
            & chunk["Subset"].astype(str).isin(subsets)
            & chunk["Mode"].astype(str).isin(modes)
            & chunk["Q_Label"].astype(str).eq(q_label)
            & chunk["Estimate_Definition"].astype(str).eq(ESTIMATE)
            & chunk["Metric"].astype(str).eq(METRIC)
        )
        sub = chunk.loc[mask].copy()
        if not sub.empty:
            frames.append(sub)
    if not frames:
        raise RuntimeError(f"No rows extracted for {domain}")
    out = pd.concat(frames, ignore_index=True)
    out["Domain"] = out["Domain"].astype(str)
    out["Subset"] = out["Subset"].astype(str)
    out["Mode"] = out["Mode"].astype(str)
    out["Cutoff_Side"] = out["Cutoff_Side"].astype(str)
    out["Retain_Pct"] = pd.to_numeric(out["Retain_Pct"], errors="coerce")
    out["Cutoff_Signed_Pct"] = pd.to_numeric(out["Cutoff_Signed_Pct"], errors="coerce")
    out["Sweep_SC"] = pd.to_numeric(out["Sweep_SC"], errors="coerce")
    if "Sweep_Is_Asymptotic" in out.columns:
        out["Sweep_Is_Asymptotic"] = out["Sweep_Is_Asymptotic"].map(parse_bool)
    else:
        out["Sweep_Is_Asymptotic"] = np.isclose(out["Sweep_SC"].astype(float), 1.0)
    out["Value"] = pd.to_numeric(out["Value"], errors="coerce")
    out["Primary_Eligible"] = out["Primary_Eligible"].map(parse_bool)
    out["Sweep_Status"] = out["Sweep_Status"].astype(str)
    return out


def support_ok(df: pd.DataFrame) -> pd.Series:
    finite_value = df["Value"].notna() & np.isfinite(df["Value"].astype(float))
    asymptotic = df.get("Sweep_Is_Asymptotic", pd.Series(False, index=df.index)).astype(bool)
    finite_primary = df["Primary_Eligible"].astype(bool)
    return finite_value & (finite_primary | asymptotic)


def draw_phrase_with_italic_symbol(
    draw: ImageDraw.ImageDraw,
    xy: Tuple[float, float],
    prefix: str,
    symbol: str,
    *,
    fill: Tuple[int, int, int],
    anchor: str = "left",
    suffix: str = "",
) -> None:
    prefix_w = draw.textlength(prefix, font=F_AXIS)
    symbol_w = draw.textlength(symbol, font=F_AXIS_ITALIC)
    suffix_w = draw.textlength(suffix, font=F_AXIS)
    x, y = xy
    if anchor == "center":
        x -= (prefix_w + symbol_w + suffix_w) / 2
    draw.text((x, y), prefix, fill=fill, font=F_AXIS)
    draw.text((x + prefix_w, y), symbol, fill=fill, font=F_AXIS_ITALIC)
    if suffix:
        draw.text((x + prefix_w + symbol_w, y), suffix, fill=fill, font=F_AXIS)


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


def grid_values(df: pd.DataFrame, scenarios: Sequence[str]) -> Tuple[List[float], List[float]]:
    sub = add_signed_cutoff_display_pct(df[df["Mode"].isin(scenarios)].copy())
    x_vals = [float(v) for v in SIGNED_CUTOFF_DISPLAY_VALUES]
    y_vals = sorted(sub["Sweep_SC"].dropna().astype(float).unique().tolist())
    return x_vals, y_vals


def scenario_arrays(
    df: pd.DataFrame,
    *,
    subset: str,
    scenario: str,
    x_vals: Sequence[float],
    y_vals: Sequence[float],
    medium: float,
    large: float,
) -> Dict[str, np.ndarray]:
    sub = df[df["Subset"].eq(subset) & df["Mode"].eq(scenario)].copy()
    if sub.empty:
        raise RuntimeError(f"No rows for subset={subset}, scenario={scenario}")
    sub = add_signed_cutoff_display_pct(sub)
    ok = support_ok(sub)
    state = pd.DataFrame(
        {
            SIGNED_CUTOFF_DISPLAY_COL: sub[SIGNED_CUTOFF_DISPLAY_COL].astype(float),
            "Sweep_SC": sub["Sweep_SC"].astype(float),
            "ok": ok.astype(bool),
            "pos_medium": (ok & (sub["Value"] >= medium)).astype(bool),
            "pos_large": (ok & (sub["Value"] >= large)).astype(bool),
            "neg_medium": (ok & (sub["Value"] <= -medium)).astype(bool),
            "neg_large": (ok & (sub["Value"] <= -large)).astype(bool),
        }
    )
    grouped = state.groupby(["Sweep_SC", SIGNED_CUTOFF_DISPLAY_COL], dropna=False).agg(
        ok=("ok", "any"),
        pos_medium=("pos_medium", "any"),
        pos_large=("pos_large", "any"),
        neg_medium=("neg_medium", "any"),
        neg_large=("neg_large", "any"),
    )
    idx = pd.MultiIndex.from_product([y_vals, x_vals], names=["Sweep_SC", SIGNED_CUTOFF_DISPLAY_COL])
    grouped = grouped.reindex(idx).fillna(False)
    shape = (len(y_vals), len(x_vals))
    return {col: grouped[col].to_numpy(dtype=bool).reshape(shape) for col in grouped.columns}


def overlay_arrays(grids: Sequence[Dict[str, np.ndarray]], direction: str) -> Dict[str, np.ndarray]:
    medium_key = "pos_medium" if direction == "positive" else "neg_medium"
    large_key = "pos_large" if direction == "positive" else "neg_large"
    ok_stack = np.stack([grid["ok"] for grid in grids], axis=0)
    med_stack = np.stack([grid[medium_key] for grid in grids], axis=0)
    large_stack = np.stack([grid[large_key] for grid in grids], axis=0)
    ok_count = ok_stack.sum(axis=0).astype(int)
    all_ok = ok_count == len(grids)
    all_bad = ok_count == 0
    partial_bad = (ok_count > 0) & (ok_count < len(grids))
    medium_count = med_stack.sum(axis=0).astype(int)
    large_count = large_stack.sum(axis=0).astype(int)
    n = len(grids)
    return {
        "ok_count": ok_count,
        "all_ok": all_ok,
        "all_bad": all_bad,
        "partial_bad": partial_bad,
        "any_bad": all_bad | partial_bad,
        "medium_count": medium_count,
        "large_count": large_count,
        "any_medium": medium_count > 0,
        "any_large": large_count > 0,
        "all_medium": all_ok & (medium_count == n),
        "all_large": all_ok & (large_count == n),
    }


def threshold_label(kind: str) -> str:
    if kind == "large":
        return "large+"
    if kind == "medium":
        return "medium+"
    raise ValueError(f"Unsupported threshold kind: {kind}")


def threshold_value(kind: str, medium: float, large: float) -> float:
    return large if kind == "large" else medium


def threshold_count_key(kind: str) -> str:
    return "large_count" if kind == "large" else "medium_count"


def threshold_any_key(kind: str) -> str:
    return "any_large" if kind == "large" else "any_medium"


def threshold_all_key(kind: str) -> str:
    return "all_large" if kind == "large" else "all_medium"


def lerp_color(a: Tuple[int, int, int], b: Tuple[int, int, int], t: float) -> Tuple[int, int, int]:
    return tuple(int(round(a[i] + (b[i] - a[i]) * t)) for i in range(3))


def count_color(direction: str, count: int, n: int) -> Tuple[int, int, int]:
    if count <= 0:
        return (255, 255, 255)
    t = (count - 1) / max(1, n - 1)
    if direction == "positive":
        return lerp_color((178, 226, 226), (5, 113, 176), t)
    return lerp_color((241, 182, 218), (197, 27, 125), t)


def draw_hatch(
    draw: ImageDraw.ImageDraw,
    box: Tuple[int, int, int, int],
    color: Tuple[int, int, int] = (205, 205, 205),
) -> None:
    x0, y0, x1, y1 = box
    step = 13
    for offset in range(-int(y1 - y0), int(x1 - x0) + step, step):
        sx = x0 + max(offset, 0)
        sy = y1 - max(-offset, 0)
        ex = x0 + min(offset + int(y1 - y0), int(x1 - x0))
        ey = y1 - min(max(int(y1 - y0) - offset, 0), int(y1 - y0))
        draw.line((sx, sy, ex, ey), fill=color, width=1)


def tick_positions(values: Sequence[float], targets: Iterable[float]) -> List[Tuple[int, str]]:
    out: List[Tuple[int, str]] = []
    used: set[int] = set()
    if not values:
        return out
    for target in targets:
        idx = min(range(len(values)), key=lambda i: abs(values[i] - target))
        if idx in used:
            continue
        used.add(idx)
        value = values[idx]
        if abs(value) >= 1:
            label = f"{value:.0f}"
        else:
            label = f"{value:.2f}".rstrip("0").rstrip(".")
        out.append((idx, label))
    return out


def format_sweep_tick_label(value: float) -> str:
    if math.isclose(float(value), 1.0):
        return "Asym."
    return f"{int(round(float(value) * 100))}"


def draw_y_ticks(
    draw: ImageDraw.ImageDraw,
    *,
    left: int,
    top: int,
    rows: int,
    cell_h: float,
    y_vals: Sequence[float],
) -> None:
    ticks: List[Tuple[int, str, float]] = []
    for idx, _ in tick_positions(y_vals, [0.01, 0.25, 0.50, 0.75, 0.99, 1.00]):
        y = int(round(top + (rows - idx - 0.5) * cell_h))
        ticks.append((y, format_sweep_tick_label(y_vals[idx]), float(y_vals[idx])))

    last_text_bottom = -10_000
    for y, y_label, _ in sorted(ticks, key=lambda item: item[0]):
        draw.line((left - 5, y, left, y), fill=(70, 70, 70), width=1)
        tw = draw.textlength(y_label, font=F_TICK)
        text_y = y - 7
        if text_y < last_text_bottom + 2:
            text_y = last_text_bottom + 2
        draw.text((left - 10 - tw, text_y), y_label, fill=(50, 50, 50), font=F_TICK)
        last_text_bottom = text_y + 14


def draw_vertical_axis_label_text(image: Image.Image, *, center: Tuple[int, int]) -> None:
    prefix = "Sample coverage "
    symbol = "s"
    suffix = " [%]"
    probe = ImageDraw.Draw(Image.new("RGBA", (1, 1)))
    prefix_w = probe.textlength(prefix, font=F_AXIS)
    symbol_w = probe.textlength(symbol, font=F_AXIS_ITALIC)
    suffix_w = probe.textlength(suffix, font=F_AXIS)
    text_w = int(prefix_w + symbol_w + suffix_w)
    label = Image.new("RGBA", (text_w + 8, 26), (255, 255, 255, 0))
    label_draw = ImageDraw.Draw(label)
    x = 4
    label_draw.text((x, 4), prefix, fill=(50, 50, 50, 255), font=F_AXIS)
    x += prefix_w
    label_draw.text((x, 4), symbol, fill=(50, 50, 50, 255), font=F_AXIS_ITALIC)
    x += symbol_w
    label_draw.text((x, 4), suffix, fill=(50, 50, 50, 255), font=F_AXIS)
    rotated = label.rotate(90, expand=True)
    image.paste(rotated, (center[0] - rotated.width // 2, center[1] - rotated.height // 2), rotated)


def draw_y_axis_label(image: Image.Image, *, left: int, top: int, bottom: int) -> None:
    draw_vertical_axis_label_text(
        image,
        center=(left - 72, int(round((top + bottom) / 2))),
    )


def draw_direction_title(draw: ImageDraw.ImageDraw, *, box: Tuple[int, int, int, int], direction: str) -> None:
    left, top, right, _ = box
    title = "Effect size δ > 0" if direction == "positive" else "Effect size δ < 0"
    tw = draw.textlength(title, font=F_PANEL_TITLE)
    draw.text(((left + right) / 2 - tw / 2, top - 30), title, fill=(30, 30, 30), font=F_PANEL_TITLE)


def axis_targets(values: Sequence[float]) -> List[float]:
    if not values:
        return []
    lo, hi = min(values), max(values)
    preferred = SIGNED_CUTOFF_TICK_TARGETS
    in_range = [v for v in preferred if lo <= v <= hi]
    if in_range:
        return in_range
    return [lo, hi]


def draw_panel(
    image: Image.Image,
    draw: ImageDraw.ImageDraw,
    overlay: Dict[str, np.ndarray],
    *,
    box: Tuple[int, int, int, int],
    x_vals: Sequence[float],
    y_vals: Sequence[float],
    direction: str,
    scenario_n: int,
    threshold_kind: str,
    show_y_label: bool = False,
) -> None:
    left, top, right, bottom = box
    count_key = threshold_count_key(threshold_kind)
    rows, cols = overlay[count_key].shape
    cell_w = (right - left) / cols
    cell_h = (bottom - top) / rows
    draw_direction_title(draw, box=box, direction=direction)
    draw.rectangle(box, fill=(255, 255, 255))
    for r in range(rows):
        for c in range(cols):
            x0 = int(round(left + c * cell_w))
            x1 = int(round(left + (c + 1) * cell_w))
            y0 = int(round(top + (rows - r - 1) * cell_h))
            y1 = int(round(top + (rows - r) * cell_h))
            if overlay["all_bad"][r, c]:
                color = (217, 217, 217)
            else:
                color = count_color(direction, int(overlay[count_key][r, c]), scenario_n)
            draw.rectangle((x0, y0, x1, y1), fill=color)
            if overlay["partial_bad"][r, c]:
                draw_hatch(draw, (x0, y0, x1, y1))
    for idx, label in tick_positions(x_vals, axis_targets(x_vals)):
        x = int(round(left + (idx + 0.5) * cell_w))
        draw.line((x, bottom, x, bottom + 5), fill=(70, 70, 70), width=1)
        tw = draw.textlength(label, font=F_TICK)
        draw.text((x - tw / 2, bottom + 8), label, fill=(50, 50, 50), font=F_TICK)
    draw_y_ticks(draw, left=left, top=top, rows=rows, cell_h=cell_h, y_vals=y_vals)
    if show_y_label:
        draw_y_axis_label(image, left=left, top=top, bottom=bottom)
    if x_vals:
        zero_idx = min(range(len(x_vals)), key=lambda i: abs(x_vals[i]))
        zx = int(round(left + (zero_idx + 0.5) * cell_w))
        draw.line((zx, top, zx, bottom), fill=(130, 130, 130), width=1)
    draw.rectangle(box, outline=(50, 50, 50), width=2)
    draw_phrase_with_italic_symbol(
        draw,
        ((left + right) / 2, bottom + 34),
        SIGNED_CUTOFF_AXIS_LABEL_PREFIX,
        "r",
        fill=(50, 50, 50),
        anchor="center",
        suffix=" [%]",
    )


def draw_legend(
    draw: ImageDraw.ImageDraw,
    *,
    x: int,
    y: int,
    scenario_n: int,
    threshold_kind: str,
) -> None:
    draw.text((x, y), "Consensus classes", fill=(40, 40, 40), font=F_LEGEND)
    y += 22
    items: List[Tuple[Tuple[int, int, int], str]] = []
    label = threshold_label(threshold_kind)
    for count in range(1, scenario_n + 1):
        items.append((count_color("positive", count, scenario_n), f"positive {label} in {count}/{scenario_n}"))
    for count in range(1, scenario_n + 1):
        items.append((count_color("negative", count, scenario_n), f"negative {label} in {count}/{scenario_n}"))
    items.append(((185, 185, 185), "some scenarios unavailable hatch"))
    items.append(((217, 217, 217), "all scenarios unavailable"))
    for color, label in items:
        if "hatch" in label:
            draw.rectangle((x, y + 4, x + 20, y + 14), fill=(255, 255, 255), outline=(90, 90, 90))
            draw_hatch(draw, (x, y + 4, x + 20, y + 14), color=color)
        else:
            draw.rectangle((x, y + 4, x + 20, y + 14), fill=color, outline=color)
        draw.text((x + 28, y), label, fill=(40, 40, 40), font=F_LEGEND)
        y += 20


def draw_overlay(
    *,
    domain: str,
    subset: str,
    scenarios: Sequence[str],
    overlay_id: str,
    title_label: str,
    df: pd.DataFrame,
    out_path: Path,
    medium: float,
    large: float,
    threshold_kind: str,
) -> List[Dict[str, object]]:
    x_vals, y_vals = grid_values(df[df["Subset"].eq(subset)], scenarios)
    if not x_vals or not y_vals:
        raise RuntimeError(f"No grid for {domain}/{subset}/{overlay_id}")
    grids = [
        scenario_arrays(
            df,
            subset=subset,
            scenario=scenario,
            x_vals=x_vals,
            y_vals=y_vals,
            medium=medium,
            large=large,
        )
        for scenario in scenarios
    ]
    overlays = {
        "positive": overlay_arrays(grids, "positive"),
        "negative": overlay_arrays(grids, "negative"),
    }
    width, height = 1600, 820
    image = Image.new("RGB", (width, height), (255, 255, 255))
    draw = ImageDraw.Draw(image)
    label = threshold_label(threshold_kind)
    cutoff = threshold_value(threshold_kind, medium, large)
    draw_panel(
        image,
        draw,
        overlays["positive"],
        box=(120, 135, 660, 680),
        x_vals=x_vals,
        y_vals=y_vals,
        direction="positive",
        scenario_n=len(scenarios),
        threshold_kind=threshold_kind,
        show_y_label=True,
    )
    draw_panel(
        image,
        draw,
        overlays["negative"],
        box=(720, 135, 1260, 680),
        x_vals=x_vals,
        y_vals=y_vals,
        direction="negative",
        scenario_n=len(scenarios),
        threshold_kind=threshold_kind,
    )
    draw_legend(draw, x=1310, y=275, scenario_n=len(scenarios), threshold_kind=threshold_kind)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    image.save(out_path, optimize=True)
    rows: List[Dict[str, object]] = []
    for direction, overlay in overlays.items():
        any_key = threshold_any_key(threshold_kind)
        all_key = threshold_all_key(threshold_kind)
        rows.append(
            {
                "Overlay_ID": overlay_id,
                "Threshold_Kind": threshold_kind,
                "Threshold_Label": label,
                "Threshold_Value": cutoff,
                "Domain": domain,
                "Subset": subset,
                "Q_Label": Q_LABEL,
                "Estimate_Definition": ESTIMATE,
                "Metric": METRIC,
                "Scenario_List": ",".join(scenarios),
                "Direction": direction,
                "Medium_Threshold": medium,
                "Large_Threshold": large,
                "Scenario_Count": len(scenarios),
                "Grid_Cells": int(overlay[threshold_count_key(threshold_kind)].size),
                "Any_MediumOrLarger_Cells": int(overlay["any_medium"].sum()),
                "AllScenario_MediumOrLarger_StatusOk_Cells": int(overlay["all_medium"].sum()),
                "AllScenario_Large_StatusOk_Cells": int(overlay["all_large"].sum()),
                "Any_ThresholdOrLarger_Cells": int(overlay[any_key].sum()),
                "AllScenario_ThresholdOrLarger_StatusOk_Cells": int(overlay[all_key].sum()),
                "Any_Missing_NotOk_Cells": int(overlay["any_bad"].sum()),
                "Partial_Missing_NotOk_Cells": int(overlay["partial_bad"].sum()),
                "All_Missing_NotOk_Cells": int(overlay["all_bad"].sum()),
                "PNG": str(out_path),
            }
        )
    return rows


def make_contact_sheet(paths: Sequence[Path], out_path: Path, title: str) -> None:
    if not paths:
        return
    cols = 3
    thumb_w, thumb_h = 480, 246
    pad = 18
    title_h = 52
    rows = (len(paths) + cols - 1) // cols
    image = Image.new("RGB", (cols * thumb_w + (cols + 1) * pad, title_h + rows * (thumb_h + 34) + (rows + 1) * pad), "white")
    draw = ImageDraw.Draw(image)
    draw.text((pad, 14), title, fill=(20, 20, 20), font=F_TITLE)
    for idx, path in enumerate(paths):
        r, c = divmod(idx, cols)
        x = pad + c * (thumb_w + pad)
        y = title_h + pad + r * (thumb_h + 34 + pad)
        label = path.stem.replace("_", " ")
        draw.text((x, y), label[:72], fill=(40, 40, 40), font=F_FOOT)
        thumb = Image.open(path).convert("RGB")
        thumb.thumbnail((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        image.paste(thumb, (x + (thumb_w - thumb.width) // 2, y + 26))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    image.save(out_path, optimize=True)


def write_outputs(
    out_dir: Path,
    *,
    root: Path,
    manifest: Sequence[Dict[str, object]],
    summary: Sequence[Dict[str, object]],
    medium: float,
    large: float,
    q_label: str,
) -> None:
    manifest_path = out_dir / "stage2_inext_medium_large_overlay_manifest.csv"
    summary_path = out_dir / "stage2_inext_medium_large_overlay_summary.csv"
    if manifest:
        with manifest_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(manifest[0].keys()))
            writer.writeheader()
            writer.writerows(manifest)
    if summary:
        with summary_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
            writer.writeheader()
            writer.writerows(summary)
    lines = [
        "# Stage2 iNEXT q2 threshold-variant consensus overlays",
        "",
        f"Source root: `{root}`",
        "",
        "Inputs:",
        "",
        "- Per-domain `csv/comparison_cells_filtered.csv` files from the Stage2 iNEXT q2 lightweight downstream run.",
        "- Rendered PNGs were not used as numeric input.",
        "",
        "Scope:",
        "",
        f"- Domains: `{', '.join(DOMAINS)}`",
        f"- Subsets: `{', '.join(SUBSETS)}`",
        f"- Hill order: `{q_label}`",
        f"- Estimate definition: `{ESTIMATE}`",
        f"- Metric: `{METRIC}`",
        f"- Medium-or-larger threshold variant: `|delta| >= {medium}`",
        f"- Large-or-larger threshold variant: `|delta| >= {large}`",
        "",
        "Interpretation boundary:",
        "",
        "- These are descriptive coverage-cutoff mask-overlap diagnostics.",
        "- The x coordinate is intensity cutoff r: 0% is no cutoff, negative values cut high-intensity mass, and positive values cut low-intensity mass.",
        "- The y coordinate is sample coverage s; the Asymptotic row is a separate endpoint and is not excluded by the finite m/n primary-range criterion.",
        "- They are not independent cell-wise hypothesis tests.",
        "- `OuterRep` and bootstrap layers remain Monte Carlo layers, not biological replicates.",
        "- No contour outlines are drawn. Fill color alone encodes how many scenarios pass the selected threshold.",
        "- Solid gray is reserved for cells where all scenarios are unavailable; partial unavailable cells retain their threshold-count fill and are marked with diagonal hatching.",
        "",
    ]
    (out_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    if args.large_threshold < args.medium_threshold:
        raise ValueError("large threshold must be >= medium threshold")
    root = args.root
    out_dir = args.out_dir or root / "summaries" / "stage2_inext_medium_large_overlays_20260622"
    domains = split_csv(args.domains)
    subsets = split_csv(args.subsets)
    q_label = str(args.q_label)
    threshold_kinds = split_csv(args.threshold_kinds)
    for kind in threshold_kinds:
        if kind not in {"medium", "large"}:
            raise ValueError(f"Unsupported threshold kind: {kind}")
    modes = list(dict.fromkeys(PERIOD_SCENARIOS + ([] if args.skip_annual else ANNUAL_SCENARIOS)))
    manifest: List[Dict[str, object]] = []
    summary: List[Dict[str, object]] = []
    period_paths: List[Path] = []
    annual_paths: List[Path] = []
    for domain in domains:
        print(f"[LOAD] {domain}", flush=True)
        cells = load_cells(root, domain, subsets, modes, args.chunksize, q_label)
        print(f"[LOAD DONE] {domain} rows={len(cells)}", flush=True)
        for subset in subsets:
            for threshold_kind in threshold_kinds:
                label = threshold_label(threshold_kind).replace("+", "plus")
                out_path = out_dir / f"period_overlay_{label}" / f"{domain.lower()}_{subset.lower()}_{q_label}_period1_period2_{label}_overlay_posneg.png"
                rows = draw_overlay(
                    domain=domain,
                    subset=subset,
                    scenarios=PERIOD_SCENARIOS,
                    overlay_id="period1_period2",
                    title_label="Period1 vs Period2",
                    df=cells,
                    out_path=out_path,
                    medium=args.medium_threshold,
                    large=args.large_threshold,
                    threshold_kind=threshold_kind,
                )
                summary.extend(rows)
                period_paths.append(out_path)
                manifest.append(
                    {
                        "Figure_Type": f"period_overlay_{label}",
                        "Threshold_Kind": threshold_kind,
                        "Domain": domain,
                        "Subset": subset,
                        "Q_Label": q_label,
                        "Estimate_Definition": ESTIMATE,
                        "Metric": METRIC,
                        "PNG": str(out_path),
                        "Status": "generated",
                    }
                )
                if not args.skip_annual:
                    annual_path = out_dir / f"annual_overlay_{label}" / f"{domain.lower()}_{subset.lower()}_{q_label}_annual6_{label}_overlay_posneg.png"
                    rows = draw_overlay(
                        domain=domain,
                        subset=subset,
                        scenarios=ANNUAL_SCENARIOS,
                        overlay_id="annual6",
                        title_label="six annual years",
                        df=cells,
                        out_path=annual_path,
                        medium=args.medium_threshold,
                        large=args.large_threshold,
                        threshold_kind=threshold_kind,
                    )
                    summary.extend(rows)
                    annual_paths.append(annual_path)
                    manifest.append(
                        {
                            "Figure_Type": f"annual_overlay_{label}",
                            "Threshold_Kind": threshold_kind,
                            "Domain": domain,
                            "Subset": subset,
                            "Q_Label": q_label,
                            "Estimate_Definition": ESTIMATE,
                            "Metric": METRIC,
                            "PNG": str(annual_path),
                            "Status": "generated",
                        }
                    )
    for threshold_kind in threshold_kinds:
        label = threshold_label(threshold_kind).replace("+", "plus")
        make_contact_sheet(
            [path for path in period_paths if f"_{label}_" in path.name],
            out_dir / "contact_sheets" / f"period_overlay_{label}_all_domains_subsets_contact_sheet.png",
            f"Period {label} overlays",
        )
        if annual_paths:
            make_contact_sheet(
                [path for path in annual_paths if f"_{label}_" in path.name],
                out_dir / "contact_sheets" / f"annual_overlay_{label}_all_domains_subsets_contact_sheet.png",
                f"Annual {label} overlays",
            )
    write_outputs(out_dir, root=root, manifest=manifest, summary=summary, medium=args.medium_threshold, large=args.large_threshold, q_label=q_label)
    print(f"[OK] wrote {len(manifest)} overlays under {out_dir}", flush=True)


if __name__ == "__main__":
    main()
