#!/usr/bin/env python3
"""Safely merge q2fast shard CSVs without loading all shards into memory."""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path
from typing import Iterable, List

import pandas as pd


REQUIRED_CSVS = [
    "comparison_outer_v36_9_q012.csv",
    "group_outer_v36_9_q012.csv",
    "group_summary_v36_9_q012.csv",
    "comparison_agg_v36_9_q012.csv",
    "comparison_cells_long_all_q012_v36_9.csv",
]
STATIC_FILES = [
    "cutoff_manifest_v36_9.csv",
    "sweep_sc_manifest_v36_9.csv",
    "run_metadata_v36_9.json",
    "run_metadata_v36_9_4gpu.json",
    "run_notes_v36_9.txt",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--expected-shards", type=int, default=0)
    return parser.parse_args()


def shard_key(path: Path) -> int:
    match = re.search(r"worker_(\d+)_", path.name)
    if not match:
        return 10**9
    return int(match.group(1))


def list_shards(run_dir: Path) -> List[Path]:
    shard_root = run_dir / "_shards"
    if not shard_root.is_dir():
        raise RuntimeError(f"Missing shard root: {shard_root}")
    shards = sorted((p for p in shard_root.iterdir() if p.is_dir()), key=shard_key)
    if not shards:
        raise RuntimeError(f"No shard directories under {shard_root}")
    return shards


def require_files(shards: Iterable[Path], filenames: Iterable[str]) -> None:
    missing: List[str] = []
    for shard in shards:
        for filename in filenames:
            if not (shard / filename).is_file():
                missing.append(str(shard / filename))
    if missing:
        preview = "\n".join(missing[:20])
        more = "" if len(missing) <= 20 else f"\n... and {len(missing) - 20} more"
        raise RuntimeError(f"Missing shard CSV(s):\n{preview}{more}")


def stream_concat_csv(shards: List[Path], filename: str, out_path: Path) -> None:
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    expected_header: str | None = None
    with tmp_path.open("w", encoding="utf-8", newline="") as out:
        for idx, shard in enumerate(shards):
            path = shard / filename
            with path.open("r", encoding="utf-8", newline="") as src:
                header = src.readline()
                if not header:
                    raise RuntimeError(f"Empty shard CSV: {path}")
                if expected_header is None:
                    expected_header = header
                    out.write(header)
                elif header != expected_header:
                    raise RuntimeError(f"Header mismatch in {path}")
                shutil.copyfileobj(src, out, length=1024 * 1024)
            print(f"[MERGE] {filename}: {idx + 1}/{len(shards)}", flush=True)
    tmp_path.replace(out_path)


def copy_static_files(shards: List[Path], run_dir: Path) -> None:
    first = shards[0]
    for filename in STATIC_FILES:
        src = first / filename
        if src.is_file():
            shutil.copy2(src, run_dir / filename)


def merge_mode_tag_manifest(shards: List[Path], run_dir: Path) -> None:
    parts = [
        pd.read_csv(path)
        for shard in shards
        if (path := shard / "mode_tag_manifest_v36_9.csv").is_file()
    ]
    if not parts:
        return
    out = pd.concat(parts, ignore_index=True).drop_duplicates()
    sort_cols = [col for col in ["Mode", "Group", "Year", "Tag"] if col in out.columns]
    if sort_cols:
        out = out.sort_values(sort_cols)
    out.to_csv(run_dir / "mode_tag_manifest_v36_9.csv", index=False)


def main() -> int:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    shards = list_shards(run_dir)
    if args.expected_shards and len(shards) != args.expected_shards:
        raise RuntimeError(f"Expected {args.expected_shards} shards, found {len(shards)}")
    require_files(shards, REQUIRED_CSVS)
    copy_static_files(shards, run_dir)
    merge_mode_tag_manifest(shards, run_dir)
    for filename in REQUIRED_CSVS:
        stream_concat_csv(shards, filename, run_dir / filename)
    (run_dir / "RUN_COMPLETE.txt").write_text("safe shard merge complete\n", encoding="utf-8")
    print(f"[OK] safe-merged {len(shards)} shards under {run_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
