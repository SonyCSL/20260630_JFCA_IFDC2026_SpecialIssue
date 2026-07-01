# Supplementary Data 7

Stage 2 Hill-number estimate summaries.

## Files

- `Data_7_manifest.tsv`: manifest of large Data 7 payloads stored outside GitHub.
- `SHA256SUMS.txt`: SHA-256 checksums for files stored in this directory.

## Description

Data 7 refers to Stage 2 diversity-estimation outputs and coverage diagnostics over the coverage grid for each domain, metabolite subset, cutoff condition, and Monte Carlo repetition.

The q=2 all-domain group summary files and the q=1 Formula/Brite/Pathway group summary file are larger than 100 MB after compression and are therefore not stored directly in GitHub. Their external download locations or pending external-upload targets are listed in `Data_7_manifest.tsv`.

The current q=2 pipeline does not retain bootstrap replicate-level rows or formula-level Stage 1 count vectors as primary deposited data objects. Group-level coverage diagnostics are provided in Data 7. Hill-number point and bootstrap-summary columns used for statistical comparison are included in the Data 8 comparison payloads.

## External Data

Large Data 7 objects are available through the links below. Existing deposited object URLs retain their original external object names to avoid link breakage.

| File | SHA-256 |
| --- | --- |
| [Data_7_fourthroot_all10_group_summary_v36_9_q2.csv.gz](https://d38f5mdcvtp0z3.cloudfront.net/JFCA_IFDC2026_SpecialIssue/D7/D7_fourthroot_all10_group_summary_v36_9_q2.csv.gz) | [sha256](https://d38f5mdcvtp0z3.cloudfront.net/JFCA_IFDC2026_SpecialIssue/D7/D7_fourthroot_all10_group_summary_v36_9_q2.csv.gz.sha256) |
| [Data_7_log1p_all10_group_summary_v36_9_q2.csv.gz](https://d38f5mdcvtp0z3.cloudfront.net/JFCA_IFDC2026_SpecialIssue/D7/D7_log1p_all10_group_summary_v36_9_q2.csv.gz) | [sha256](https://d38f5mdcvtp0z3.cloudfront.net/JFCA_IFDC2026_SpecialIssue/D7/D7_log1p_all10_group_summary_v36_9_q2.csv.gz.sha256) |
| [Data_7_sqrt_all10_group_summary_v36_9_q2.csv.gz](https://d38f5mdcvtp0z3.cloudfront.net/JFCA_IFDC2026_SpecialIssue/D7/D7_sqrt_all10_group_summary_v36_9_q2.csv.gz) | [sha256](https://d38f5mdcvtp0z3.cloudfront.net/JFCA_IFDC2026_SpecialIssue/D7/D7_sqrt_all10_group_summary_v36_9_q2.csv.gz.sha256) |

## Local Payloads Pending External Upload

The following q=1 three-domain payload has been compressed under
`_S3_payloads/Data 7/` and is not tracked by Git:

| File | Scope | Status |
| --- | --- | --- |
| `Data_7_q1_formula_brite_pathway_group_summary_v36_9.csv.gz` | q=1, raw, Formula/Brite/Pathway, All/PM/SM, 199 cutoffs | `local_payload_pending_external_upload` |
