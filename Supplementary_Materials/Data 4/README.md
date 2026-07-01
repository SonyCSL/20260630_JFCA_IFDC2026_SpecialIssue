# Supplementary Data 4

Relative-abundance distributions for signed intensity-cutoff conditions.

## Files

- `Data_4_relative_abundance_distributions.csv.gz`: dense relative-abundance table for 199 cutoff conditions, 13 LC-MS profiles, and 541 formula-level features.
- `Data_4_cutoff_manifest.csv`: cutoff identifiers and signed cutoff definitions.
- `Data_4_profile_manifest.csv`: LC-MS profile labels and paired-comparison inclusion fields.
- `Data_4_validation_summary.json`: row-count and probability-sum validation summary.
- `SHA256SUMS.txt`: SHA-256 checksums.

## Description

Data 4 provides the fixed relative-abundance distributions used in downstream analyses. For each LC-MS profile and cutoff condition, retained intensities are normalized to sum to one. Formulae not retained under a given profile-cutoff condition are retained in the dense table with `retained=false` and `relative_abundance_p=0`.

The cutoff axis contains 199 conditions: `low_01` to `low_99`, `none_100`, and `high_01` to `high_99`. Low cutoffs retain the low-intensity tail, high cutoffs retain the high-intensity tail, and `none_100` retains all positive-intensity formulae.

## Validation

The generated table contains 1,399,567 data rows, corresponding to 199 cutoff conditions x 13 profiles x 541 formula-level features. The maximum absolute error in per-distribution probability sums is recorded in `Data_4_validation_summary.json`.

## Integrity

File checksums are provided in `SHA256SUMS.txt`.
