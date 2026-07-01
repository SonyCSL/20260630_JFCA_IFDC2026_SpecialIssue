#!/usr/bin/env python3
"""Materialize small per-domain Stage2 iNEXT downstream inputs.

The v6 downstream plotter intentionally preserves the established figure logic,
but its full summary prepass reads large upstream CSVs into memory.  This helper
streams the completed upstream CSVs and writes filtered run directories that keep
only the q2 Cliff's delta scope needed for the urgent figure export.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from pathlib import Path
from typing import Iterable

import pandas as pd


CELLS_NAME = "comparison_cells_long_all_q012_v36_9.csv"
OUTER_NAME = "comparison_outer_v36_9_q012.csv"
GROUP_NAME = "group_outer_v36_9_q012.csv"

DEFAULT_DOMAINS = ["Formula", "Brite", "Network", "Disease_NE"]
DEFAULT_SUBSETS = ["All", "PM", "SM"]
DEFAULT_MODES = [
    "2015",
    "2016",
    "2017",
    "2018",
    "2019",
    "2020",
    "Period1_2015_2017",
    "Period2_2018_2020",
    "Combo6yr",
]

CELL_USECOLS = [
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
]

OUTER_USECOLS = [
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
    "CliffsDelta",
]

GROUP_USECOLS = [
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
]


def parse_csv_list(raw: str) -> list[str]:
    return [item.strip() for item in raw.split(",") if item.strip()]


def now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(message: str) -> None:
    print(f"[{now()}] {message}", flush=True)


def ensure_columns(path: Path, usecols: Iterable[str]) -> None:
    header = pd.read_csv(path, nrows=0)
    missing = [col for col in usecols if col not in header.columns]
    if missing:
        raise RuntimeError(f"{path} is missing required columns: {missing}")


def write_empty_outputs(out_root: Path, domains: list[str], filename: str, columns: list[str]) -> None:
    for domain in domains:
        out_dir = out_root / domain
        out_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=columns).to_csv(out_dir / filename, index=False)


def append_domain_parts(
    chunk: pd.DataFrame,
    *,
    filename: str,
    out_root: Path,
    domains: list[str],
    initialized: dict[str, bool],
    counts: dict[str, int],
) -> None:
    if chunk.empty:
        return
    for domain, part in chunk.groupby("Domain", sort=False):
        if domain not in initialized:
            continue
        out_path = out_root / domain / filename
        mode = "a" if initialized[domain] else "w"
        part.to_csv(out_path, mode=mode, header=not initialized[domain], index=False)
        initialized[domain] = True
        counts[domain] += int(len(part))


def filter_common(
    chunk: pd.DataFrame,
    *,
    domains: list[str],
    subsets: list[str],
    modes: list[str],
) -> pd.DataFrame:
    mask = chunk["Domain"].isin(domains)
    if "Subset" in chunk.columns:
        mask &= chunk["Subset"].isin(subsets)
    if "Mode" in chunk.columns:
        mask &= chunk["Mode"].isin(modes)
    return chunk.loc[mask].copy()


def split_one_csv(
    *,
    source_path: Path,
    filename: str,
    usecols: list[str],
    out_root: Path,
    domains: list[str],
    subsets: list[str],
    modes: list[str],
    q_label: str | None,
    metric: str | None,
    chunksize: int,
) -> dict[str, int]:
    ensure_columns(source_path, usecols)
    initialized = {domain: False for domain in domains}
    counts = {domain: 0 for domain in domains}
    chunk_count = 0
    log(f"split start: {source_path}")
    for chunk in pd.read_csv(source_path, usecols=usecols, chunksize=chunksize, low_memory=False):
        chunk_count += 1
        filtered = filter_common(chunk, domains=domains, subsets=subsets, modes=modes)
        if q_label is not None and "Q_Label" in filtered.columns:
            filtered = filtered[filtered["Q_Label"] == q_label].copy()
        if metric is not None and "Metric" in filtered.columns:
            filtered = filtered[filtered["Metric"] == metric].copy()
        append_domain_parts(
            filtered,
            filename=filename,
            out_root=out_root,
            domains=domains,
            initialized=initialized,
            counts=counts,
        )
        if chunk_count % 10 == 0:
            log(f"  {filename}: chunks={chunk_count}, rows={counts}")
    for domain in domains:
        if not initialized[domain]:
            # Keep downstream failure mode explicit but schema-compatible.
            pd.DataFrame(columns=usecols).to_csv(out_root / domain / filename, index=False)
    log(f"split complete: {filename}, chunks={chunk_count}, rows={counts}")
    return counts


def file_sizes(paths: Iterable[Path]) -> dict[str, int]:
    return {str(path): path.stat().st_size for path in paths if path.exists()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-run", required=True, type=Path)
    parser.add_argument("--out-root", required=True, type=Path)
    parser.add_argument("--domains", default=",".join(DEFAULT_DOMAINS))
    parser.add_argument("--subsets", default=",".join(DEFAULT_SUBSETS))
    parser.add_argument("--modes", default=",".join(DEFAULT_MODES))
    parser.add_argument("--q-label", default="q2")
    parser.add_argument("--metric", default="CliffsDelta")
    parser.add_argument("--chunksize", default=500_000, type=int)
    parser.add_argument("--replace", action="store_true")
    args = parser.parse_args()

    source_run = args.source_run.resolve()
    out_root = args.out_root.resolve()
    domains = parse_csv_list(args.domains)
    subsets = parse_csv_list(args.subsets)
    modes = parse_csv_list(args.modes)
    if not domains:
        raise RuntimeError("At least one domain is required.")
    if args.chunksize <= 0:
        raise RuntimeError("--chunksize must be positive.")
    if out_root.exists() and args.replace:
        shutil.rmtree(out_root)
    if out_root.exists() and any(out_root.iterdir()):
        raise RuntimeError(f"Output root already exists and is not empty: {out_root}")
    out_root.mkdir(parents=True, exist_ok=True)
    for domain in domains:
        (out_root / domain).mkdir(parents=True, exist_ok=True)

    source_paths = {
        CELLS_NAME: source_run / CELLS_NAME,
        OUTER_NAME: source_run / OUTER_NAME,
        GROUP_NAME: source_run / GROUP_NAME,
    }
    for path in source_paths.values():
        if not path.exists():
            raise RuntimeError(f"Missing source CSV: {path}")

    started = now()
    manifest: dict[str, object] = {
        "started_at": started,
        "source_run": str(source_run),
        "out_root": str(out_root),
        "domains": domains,
        "subsets": subsets,
        "modes": modes,
        "q_label": args.q_label,
        "metric": args.metric,
        "chunksize": args.chunksize,
        "inputs": file_sizes(source_paths.values()),
        "row_counts": {},
    }

    manifest["row_counts"][CELLS_NAME] = split_one_csv(
        source_path=source_paths[CELLS_NAME],
        filename=CELLS_NAME,
        usecols=CELL_USECOLS,
        out_root=out_root,
        domains=domains,
        subsets=subsets,
        modes=modes,
        q_label=args.q_label,
        metric=args.metric,
        chunksize=args.chunksize,
    )
    manifest["row_counts"][OUTER_NAME] = split_one_csv(
        source_path=source_paths[OUTER_NAME],
        filename=OUTER_NAME,
        usecols=OUTER_USECOLS,
        out_root=out_root,
        domains=domains,
        subsets=subsets,
        modes=modes,
        q_label=args.q_label,
        metric=None,
        chunksize=args.chunksize,
    )
    manifest["row_counts"][GROUP_NAME] = split_one_csv(
        source_path=source_paths[GROUP_NAME],
        filename=GROUP_NAME,
        usecols=GROUP_USECOLS,
        out_root=out_root,
        domains=domains,
        subsets=subsets,
        modes=modes,
        q_label=None,
        metric=None,
        chunksize=args.chunksize,
    )

    output_files: list[Path] = []
    for domain in domains:
        domain_dir = out_root / domain
        (domain_dir / "RUN_COMPLETE.txt").write_text("filtered downstream split complete\n", encoding="utf-8")
        output_files.extend([domain_dir / CELLS_NAME, domain_dir / OUTER_NAME, domain_dir / GROUP_NAME])
    manifest["completed_at"] = now()
    manifest["outputs"] = file_sizes(output_files)
    (out_root / "split_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    log(f"split all complete: {out_root}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr, flush=True)
        raise
