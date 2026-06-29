# Supplementary Data D8

Synecoculture versus conventional comparison summaries.

## Files

- `D8_manifest.tsv`: manifest of large D8 payloads stored outside GitHub.
- `SHA256SUMS.txt`: SHA-256 checksums for files stored in this directory.

## Description

D8 refers to statistical comparison outputs between Synecoculture and conventional tea profiles. The deposited payloads include aggregated comparison summaries (`comparison_agg`) and, where needed for replicate-level review or figure regeneration, Monte Carlo repetition-level comparison outputs (`comparison_outer`).

The q=2 all-domain comparison summary files are larger than 100 MB after compression and are therefore not stored directly in GitHub. Their S3 object locations are listed in `D8_manifest.tsv`.

The aggregated `comparison_agg` files provide median and quantile summaries over the 24 Monte Carlo repetitions. The `comparison_outer` files preserve the 24 repetition-level contrast values and are required for regenerating box-and-whisker distribution plots from underlying contrast values.

## External Data

Large D8 objects are listed under the S3 prefix:

`s3://20260630-jfca-ifdc2026-specialissue/Chao1_Intensity/Supplementary_Materials/D8/`

Manifest rows marked `pending_external_upload` identify objects that still need to be deposited at the listed URI.
