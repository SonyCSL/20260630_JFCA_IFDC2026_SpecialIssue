# Supplementary Data 10

Figure-output manifests and figure-support tables.

## Files

- `Data_10_manifest.tsv`: manifest of large Data 10 figure bundles stored outside GitHub.
- `Figures/PPT/`: PowerPoint review decks tracked directly in GitHub.
- `manifests/`: small figure, overlay, and combo manifests included directly in GitHub.
- `small_tables/`: small compressed support tables for selected Figure 01 boxplot inputs.
- `SHA256SUMS.txt`: SHA-256 checksums for files stored in this directory.

## Description

Data 10 records visualization and result-reporting outputs. The current deposited scope includes q=2 figure-support material for the analyzed domain-subset combinations. A q=1 raw Formula/Brite/Pathway figure bundle has also been prepared locally after main-downstream style synchronization.

The `Figures/PPT/` directory contains GitHub-tracked review PowerPoint decks for q=0, q=1, and q=2. These decks are below the per-file GitHub size limit and are listed in `manifests/Data_10_figures_ppt_manifest_20260705.tsv`.

The q=2 Coverage-cutoff contrast surface plot deck was rendered downstream from the existing labelrev4/raw_all10 cell-level CSV outputs. No upstream diversity estimates were recomputed for that deck.

Large figure bundles are not stored directly in GitHub because they exceed the 100 MB repository threshold. Their external download locations are listed in `Data_10_manifest.tsv`.

## External Data

Large Data 10 objects are available through the links below. Existing deposited object URLs retain their original external object names to avoid link breakage.

| File | SHA-256 |
| --- | --- |
| [Data_10_q2_local_figure_bundles_20260629.tar.zst](https://d38f5mdcvtp0z3.cloudfront.net/JFCA_IFDC2026_SpecialIssue/D10/D10_q2_local_figure_bundles_20260629.tar.zst) | [sha256](https://d38f5mdcvtp0z3.cloudfront.net/JFCA_IFDC2026_SpecialIssue/D10/D10_q2_local_figure_bundles_20260629.tar.zst.sha256) |

## Local Payloads Pending External Upload

The following q=1 three-domain figure bundle has been compressed under
`_S3_payloads/Data 10/` and is not tracked by Git:

| File | Scope | Status |
| --- | --- | --- |
| `Data_10_q1_formula_brite_pathway_style_sync_figures_20260630.tar.zst` | q=1, raw, Formula/Brite/Pathway, All/PM/SM, annual/period/Combo, figure families `01`, `06`, `06_2` | `local_payload_pending_external_upload` |
