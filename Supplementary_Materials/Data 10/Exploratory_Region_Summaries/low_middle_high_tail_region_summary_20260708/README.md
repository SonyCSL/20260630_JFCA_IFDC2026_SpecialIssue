# Low / middle / high signed-tail region summary

This directory contains a downstream-only diagnostic summary of existing q=2
Cliff's delta coverage-cutoff surfaces.

## Source

Source lineage:

```text
/Users/nonaka/Documents/work/Rarefaction/Chao1_Intensity/out_heatmap_cutoff_quantification_from_csv/stage2_q2_mask_labelrev4_20260627/raw_all10
```

Primary filter:

```text
Q_Label = q2
Estimate_Definition = Stage2_iNEXT_TD_m_est
Metric = CliffsDelta
Domains = Formula, Brite, Pathway
Subsets = All, PM, SM
Modes = Combo6yr, Period1_2015_2017, Period2_2018_2020
Large-effect threshold = |Cliff's delta| >= 0.474
```

## Region bins

The summaries split the signed `Cutoff_Signed_Pct` axis into contiguous
signed-rank bins for `n = 1, 2, 3, 4, 5, 6`.
The exact bin membership is in `cutoff_region_membership_n1_to_n6.tsv`.

`Middle-tail` is an operational label for this diagnostic, not an established
statistical or ecological term.  The more precise phrase is `middle
signed-cutoff region`.

## Denominator

The primary summaries use the finite-grid denominator:

```text
19,701 = 199 cutoff levels x 99 finite Sweep_SC levels
```

The `Sweep_Set=all100` rows include `Sweep_SC=1.0` / `Asymptotic`; these rows
are provided for traceability but should not be treated as ordinary finite
target-coverage cells.

## Interpretation boundary

These files summarize evaluation-grid area in existing downstream surfaces.
They do not rerun Stage 2 estimation, create new p-values, create confidence
intervals, estimate biological replicate-level uncertainty, or treat
`cutoff x Sweep_SC` cells as independent experimental replicates.

The `n=1..6` split is an exploratory partition-sensitivity check.  It is not a
set of independent confirmatory tests.

## Main files

- `region_summary_all_n.tsv`: complete region summaries for all `n=1..6`.
- `dominant_region_summary_all_n.tsv`: dominant signed large-effect region per
  domain/subset/mode/cultivation.
- `table1_style_region_labels_n3.tsv`: Combo6yr broad labels for `n=3`.
- `table2_style_period_common_region_labels_n3.tsv`: Period1/Period2 common
  signed large-effect labels for `n=3`.
- `region_support_heatmap_n3.png` and `region_support_heatmap_n6.png`: first
  review figures.
