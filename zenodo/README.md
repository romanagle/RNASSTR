# RNASSTR v2 Zenodo archive

The Zenodo record should be titled **RNASSTR v2: RNA sequence–secondary-structure pairs and model benchmark outputs** and released under CC BY 4.0.

## Files to archive

### Dataset partitions

- `train.csv` — RNASSTR training partition
- `val.csv` — RNASSTR validation partition
- `test.csv` — RNASSTR held-out test partition
- `csv_split_summary.tsv`
- `family_sequence_counts.tsv`

### Model predictions and scores

- published and RNASSTR-retrained SincFold raw test predictions
- combined per-sequence SincFold scores
- Lyra-TransPred raw test predictions
- normalized per-sequence Lyra-TransPred scores
- RNAfold subsample, predictions, and per-sequence scores
- matched RNAfold/SincFold predicted-pair-density analysis

### Metadata and reproducibility

- family and sequence metadata tables
- sampled sequences, Infernal `cmscan` tabular output, and structural-hit records used to construct the final partitions
- validation reports
- release manifest with SHA-256 checksums
- source archive of the matching GitHub release

Do not upload manuscript drafts, reviewer correspondence, `.DS_Store` files, or Python caches.

The DOI and checksums will be added after the final files are frozen.
