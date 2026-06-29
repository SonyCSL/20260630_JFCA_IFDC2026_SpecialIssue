# Supplementary Data D7

Stage 2 Hill-number estimate summaries.

## Files

- `D7_manifest.tsv`: manifest of large D7 payloads stored outside GitHub.
- `SHA256SUMS.txt`: SHA-256 checksums for files stored in this directory.

## Description

D7 refers to Stage 2 diversity-estimation outputs and coverage diagnostics over the coverage grid for each domain, metabolite subset, cutoff condition, and Monte Carlo repetition.

The q=2 all-domain group summary files are larger than 100 MB after compression and are therefore not stored directly in GitHub. Their S3 object locations are listed in `D7_manifest.tsv`.

The current q=2 pipeline does not retain bootstrap replicate-level rows or formula-level Stage 1 count vectors as primary deposited data objects. Group-level coverage diagnostics are provided in D7. Hill-number point and bootstrap-summary columns used for statistical comparison are included in the D8 comparison payloads.

## External Data

Large D7 objects are listed under the S3 prefix:

`s3://20260630-jfca-ifdc2026-specialissue/Chao1_Intensity/Supplementary_Materials/D7/`

Manifest rows marked `pending_external_upload` identify objects that still need to be deposited at the listed URI.
