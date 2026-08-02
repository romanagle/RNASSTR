# RNASSTR

RNASSTR is a large, structure-aware collection of RNA sequence–secondary-structure pairs derived from Rfam covariance models, reference genomes, and Rfam full alignments. RNASSTR is described in *Improving RNA Secondary Structure Prediction Through Expanded Training Data*.

The complete dataset and large per-sequence model outputs are distributed through Zenodo: [https://doi.org/10.5281/zenodo.15319167](https://doi.org/10.5281/zenodo.15319167).

## Contents

- `scripts/shared/rfam_utils.py`: parse Rfam covariance-model metadata and family-to-clan assignments.
- `scripts/dataset/make_struct_clan_splits_v3.py`: construct structural- and clan-aware family partitions.
- `scripts/dataset/apply_family_splits_to_sto_with_family.py`: project Stockholm consensus structures onto individual sequences and write the final partition CSVs.
- `scripts/dataset/filter_rnasstr_candidates.py`: apply the sequence and structure quality-control criteria used to generate RNASSTR.
- `results/global/`: compact model-wide summary tables.
- `results/per_family/`: family-level performance tables used in the manuscript and supplementary figures.
- `release/`: dataset and benchmark manifests.
- `zenodo/`: metadata and inventory for the archived data release.

## Models evaluated

The revised analysis reports:

- SincFold using its published parameters;
- SincFold retrained using RNASSTR;
- Lyra-TransPred trained using RNASSTR; and
- RNAfold on a deterministic, family-stratified subset of the test partition.

## Installation

Create a Python 3.10 or later environment and install the metadata-parsing dependency:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Build a combined Rfam covariance-model and clan metadata table with:

```bash
python scripts/shared/rfam_utils.py --cm Rfam.cm --clans Rfam.clanin --out rfam_metadata.tsv
```

The large dataset and model-output files used in the paper are archived on [Zenodo](https://doi.org/10.5281/zenodo.15319167).

## Dataset-generation workflow

Dataset generation used Rfam v14.10, GTDB release 214, NCBI RefSeq release 229, and Infernal v1.1.5. See the manuscript Methods for the complete workflow and filtering criteria. The initial searches were run manually for each Rfam family. A representative command was:

```bash
cmsearch --cpu 32 --cut_ga --tblout RFxxxxx.tblout -A RFxxxxx.sto RFxxxxx.cm reference_sequences.fna > RFxxxxx.cmsearch.log
```

`filter_rnasstr_candidates.py` implements the quality-control criteria described in the manuscript Methods. No phylogenetic filter is applied.

```bash
python scripts/dataset/filter_rnasstr_candidates.py --candidates merged_candidates.csv --rfam-reference rfam_reference_sequences.csv --output filtered_candidates.csv --rejections qc_rejections.csv --reference-summary rfam_qc_thresholds.csv
```

The final structural partitions were generated with:

```bash
python scripts/dataset/make_struct_clan_splits_v3.py --sto-dir STO_dedupe --rfam-cm Rfam.cm --clanin Rfam.clanin --outdir splits_struct_clan_v3_fast10 --max-seqs-per-family 10 --evalue 0.01 --cpu 64 --min-val-families 300 --min-test-families 300 --min-val-components 150 --min-test-components 150
```

The family assignments were applied to the Stockholm files with:

```bash
python scripts/dataset/apply_family_splits_to_sto_with_family.py --sto-dir STO_dedupe --train-families splits_struct_clan_v3_fast10/train_families.txt --val-families splits_struct_clan_v3_fast10/val_families.txt --test-families splits_struct_clan_v3_fast10/test_families.txt --outdir csv_splits_v3_fast10
```

## Citation

Please cite the manuscript and the [RNASSTR Zenodo record](https://doi.org/10.5281/zenodo.15319167).

## Licenses

Code in this repository is released under the [MIT License](LICENSE). The RNASSTR data archive is intended for release under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); see [`DATA_LICENSE.md`](DATA_LICENSE.md).
