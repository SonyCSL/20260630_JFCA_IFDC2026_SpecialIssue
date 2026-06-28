# Supplementary Data D10

Figure-output manifests and figure-support tables.

## Files

- `D10_manifest.tsv`: manifest of large D10 figure bundles stored outside GitHub.
- `manifests/`: small figure, overlay, and combo manifests included directly in GitHub.
- `small_tables/`: small compressed support tables for selected Figure 01 boxplot inputs.
- `SHA256SUMS.txt`: SHA-256 checksums for files stored in this directory.

## Description

D10 records visualization and result-reporting outputs. The current deposited scope is q=2 figure-support material for the analyzed domain-subset combinations. q=0 and q=1 figure outputs may be added later if generated and retained in the Methods.

Large figure bundles are not stored directly in GitHub because they exceed the 100 MB repository threshold. Their S3 object locations are listed in `D10_manifest.tsv`.

## External Data

Large D10 objects are listed under the S3 prefix:

`s3://20260630-jfca-ifdc2026-specialissue/Chao1_Intensity/Supplementary_Materials/D10/`

Manifest rows marked `pending_external_upload` identify objects that still need to be deposited at the listed URI.
