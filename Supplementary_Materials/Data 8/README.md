# Supplementary Data 8

Synecoculture versus conventional comparison summaries.

## Files

- `Data_8_manifest.tsv`: manifest of large Data 8 payloads stored outside GitHub.
- `SHA256SUMS.txt`: SHA-256 checksums for files stored in this directory.

## Description

Data 8 refers to statistical comparison outputs between Synecoculture and conventional tea profiles. The deposited payloads include aggregated comparison summaries (`comparison_agg`) and, where needed for replicate-level review or figure regeneration, Monte Carlo repetition-level comparison outputs (`comparison_outer`).

The q=2 all-domain comparison summary files and the q=1 Formula/Brite/Pathway comparison payloads are larger than 100 MB after compression and are therefore not stored directly in GitHub. Their external download locations or pending external-upload targets are listed in `Data_8_manifest.tsv`.

The aggregated `comparison_agg` files provide median and quantile summaries over the 24 Monte Carlo repetitions. The `comparison_outer` files preserve the 24 repetition-level contrast values and are required for regenerating box-and-whisker distribution plots from underlying contrast values.

## External Data

Large Data 8 objects are available through the links below. Existing deposited object URLs retain their original external object names to avoid link breakage.

| File | SHA-256 |
| --- | --- |
| [Data_8_fourthroot_all10_comparison_agg_v36_9_q2.csv.gz](https://d38f5mdcvtp0z3.cloudfront.net/JFCA_IFDC2026_SpecialIssue/D8/D8_fourthroot_all10_comparison_agg_v36_9_q2.csv.gz) | [sha256](https://d38f5mdcvtp0z3.cloudfront.net/JFCA_IFDC2026_SpecialIssue/D8/D8_fourthroot_all10_comparison_agg_v36_9_q2.csv.gz.sha256) |
| [Data_8_fourthroot_all10_comparison_outer_v36_9_q2.csv.zst](https://d38f5mdcvtp0z3.cloudfront.net/JFCA_IFDC2026_SpecialIssue/D8/D8_fourthroot_all10_comparison_outer_v36_9_q2.csv.zst) | [sha256](https://d38f5mdcvtp0z3.cloudfront.net/JFCA_IFDC2026_SpecialIssue/D8/D8_fourthroot_all10_comparison_outer_v36_9_q2.csv.zst.sha256) |
| [Data_8_log1p_all10_comparison_agg_v36_9_q2.csv.gz](https://d38f5mdcvtp0z3.cloudfront.net/JFCA_IFDC2026_SpecialIssue/D8/D8_log1p_all10_comparison_agg_v36_9_q2.csv.gz) | [sha256](https://d38f5mdcvtp0z3.cloudfront.net/JFCA_IFDC2026_SpecialIssue/D8/D8_log1p_all10_comparison_agg_v36_9_q2.csv.gz.sha256) |
| [Data_8_log1p_all10_comparison_outer_v36_9_q2.csv.zst](https://d38f5mdcvtp0z3.cloudfront.net/JFCA_IFDC2026_SpecialIssue/D8/D8_log1p_all10_comparison_outer_v36_9_q2.csv.zst) | [sha256](https://d38f5mdcvtp0z3.cloudfront.net/JFCA_IFDC2026_SpecialIssue/D8/D8_log1p_all10_comparison_outer_v36_9_q2.csv.zst.sha256) |
| [Data_8_sqrt_all10_comparison_agg_v36_9_q2.csv.gz](https://d38f5mdcvtp0z3.cloudfront.net/JFCA_IFDC2026_SpecialIssue/D8/D8_sqrt_all10_comparison_agg_v36_9_q2.csv.gz) | [sha256](https://d38f5mdcvtp0z3.cloudfront.net/JFCA_IFDC2026_SpecialIssue/D8/D8_sqrt_all10_comparison_agg_v36_9_q2.csv.gz.sha256) |
| [Data_8_sqrt_all10_comparison_outer_v36_9_q2.csv.zst](https://d38f5mdcvtp0z3.cloudfront.net/JFCA_IFDC2026_SpecialIssue/D8/D8_sqrt_all10_comparison_outer_v36_9_q2.csv.zst) | [sha256](https://d38f5mdcvtp0z3.cloudfront.net/JFCA_IFDC2026_SpecialIssue/D8/D8_sqrt_all10_comparison_outer_v36_9_q2.csv.zst.sha256) |

## Local Payloads Pending External Upload

The following q=1 three-domain payloads have been compressed under
`_S3_payloads/Data 8/` and are not tracked by Git:

| File | Scope | Status |
| --- | --- | --- |
| `Data_8_q1_formula_brite_pathway_comparison_agg_v36_9.csv.gz` | q=1 aggregated comparison summaries | `local_payload_pending_external_upload` |
| `Data_8_q1_formula_brite_pathway_comparison_outer_v36_9.csv.zst` | q=1 `OuterRep=24` repetition-level comparison values for boxplot regeneration | `local_payload_pending_external_upload` |
