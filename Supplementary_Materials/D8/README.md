# Supplementary Data D8

Synecoculture versus conventional comparison summaries.

## Files

- `D8_manifest.tsv`: manifest of large D8 payloads stored outside GitHub.
- `SHA256SUMS.txt`: SHA-256 checksums for files stored in this directory.

## Description

D8 refers to statistical comparison summaries between Synecoculture and conventional tea profiles. The primary comparison output includes effect-size summaries such as Cliff's delta and related diagnostic statistics over the coverage-cutoff grid.

The q=2 all-domain comparison summary files are larger than 100 MB after compression and are therefore not stored directly in GitHub. Their S3 object locations are listed in `D8_manifest.tsv`.

Full `comparison_outer` tables are not treated as the default D8 payload because they are substantially larger than the aggregated comparison summaries.

## External Data

Large D8 objects are listed under the S3 prefix:

`s3://20260630-jfca-ifdc2026-specialissue/Chao1_Intensity/Supplementary_Materials/D8/`

Manifest rows marked `pending_external_upload` identify objects that still need to be deposited at the listed URI.
