# Code freeze: 2026-06-30 submission state

This directory preserves a lightweight code snapshot for the analysis state used
around the 2026-06-30 23:00 JST manuscript submission.

The purpose is reproducibility of the submitted analysis lineage, not a new
analysis run. Large output tables, figure bundles, rendered manuscript pages,
and S3 payload archives are intentionally excluded.

## Remote target

This freeze is intended for the personal publication repository inside the
nested `GitHub/` repository, not the SonyCSL organization repository:

```text
origin = git@github.com:Ariakei/20260630_JFCA_IFDC2026_SpecialIssue.git
company = git@github.com:SonyCSL/20260630_JFCA_IFDC2026_SpecialIssue.git
```

## Contents

- `code/upstream/`: primary Stage 2 q-fast upstream source.
- `code/downstream/`: primary Stage 2 iNEXT downstream source.
- `code/wrappers/`: June 29-30 execution wrappers and watcher scripts.
- `code/renderers/`: submission-adjacent figure and deck rendering helpers.
- `code/support/`: safe merge, split, and Retention/Jaccard support utilities.
- `manifest.tsv`: checksum, source path, mtime, tracked state, and role for each
  frozen file.
- `source_root_git_status_at_freeze.txt`: source-root working-tree state at the
  time this freeze was assembled.
- `source_root_git_diff_at_freeze.patch`: source-root tracked-file diff summary
  and text diff at freeze time.
- `github_repo_status_at_freeze.txt`: nested `GitHub/` repository status before
  staging this code freeze.
- `github_repo_diff_at_freeze.patch`: nested `GitHub/` repository tracked-file
  diff before staging this code freeze.

## Time boundary

The target boundary is `2026-06-30 23:00 JST`.

Most frozen files have mtimes at or before that boundary. One wrapper is kept
because it belongs to the same submission-adjacent code path and records the
Formula q2 raw Bootstrap=50 urgent check:

```text
scripts/run_stage2_q2_formula_raw_boot50_gpupc2_20260630.sh
mtime: 2026-06-30 23:13:37 JST
```

This post-23:00 file is explicitly marked as `post_2300_jst = yes` in
`manifest.tsv`. It should not be interpreted as evidence that all code here was
strictly frozen before 23:00.

## Scientific interpretation

This freeze preserves code paths that estimate or summarize:

- Stage 1 finite reference count vectors at `Reference_Target_SC = 0.965`;
- Stage 2 Chao/iNEXT-style coverage sweeps over finite reference count vectors;
- q-specific upstream outputs for the submitted workflow;
- downstream Cliff's delta summaries and figure generation;
- SSC sensitivity wrappers relevant to the submission-facing conclusions.

It does not create a new analysis, validate every submitted figure, or convert
Monte Carlo repetitions into independent biological replicates. `OuterRep`,
`Boot`, and shard-level repetitions remain workflow-level Monte Carlo
components unless separately justified.

## Reproduction notes

Use the files under `code/` as the preserved source of the submission-adjacent
execution path. For actual reruns, copy or reference them in an environment with
the expected project layout, data files, Python environment, and GPU resources.

Downstream wrappers should set `TOMATO_RUN_DIR` explicitly. Do not rely on
auto-detection when multiple upstream runs exist.

The manifest checksums verify the frozen copies, not the original source files
after later edits.
