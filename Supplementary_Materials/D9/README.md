# Supplementary Data D9

Summary modes and sensitivity-analysis records.

## Files

- `D9_manifest.tsv`: manifest of large D9 payloads stored outside GitHub.
- `manifests/`: small input and figure manifests included directly in GitHub.
- `SHA256SUMS.txt`: SHA-256 checksums for files stored in this directory.

## Description

D9 records the summary-mode and sensitivity-analysis materials associated with the q=2 workflow, including input manifests, figure manifests, run-scope descriptions, and links to large summary payloads.

The repository stores small manifests directly. Larger summary archives are referenced in `D9_manifest.tsv` and should be retrieved from S3 once deposited.

Intermediate run directories, shard directories, and full `comparison_outer` payloads are not included as the default D9 data package because they are much larger than the summary objects needed to identify and review the reported analyses.

## External Data

Large D9 objects are listed under the S3 prefix:

`s3://20260630-jfca-ifdc2026-specialissue/Chao1_Intensity/Supplementary_Materials/D9/`

Manifest rows marked `pending_external_upload` identify objects that still need to be deposited at the listed URI.
