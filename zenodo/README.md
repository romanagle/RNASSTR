# RNASSTR v2 Zenodo archive

The RNASSTR dataset and model-output files are archived at [https://doi.org/10.5281/zenodo.15319167](https://doi.org/10.5281/zenodo.15319167) under CC BY 4.0.

## Archived files

### Dataset partitions

- `train.csv` — RNASSTR training partition
- `val.csv` — RNASSTR validation partition
- `test.csv` — RNASSTR held-out test partition

### Model predictions and scores

- published and RNASSTR-retrained SincFold raw test predictions
- combined per-sequence SincFold scores
- Lyra-TransPred raw test predictions
- normalized per-sequence Lyra-TransPred scores
- RNAfold subsample, predictions, and per-sequence scores

### Structural-split provenance

- sampled sequences, Infernal `cmscan` tabular output, and structural-hit records used to construct the final partitions
