# RNASSTR

RNASSTR is a large, structure-aware collection of RNA sequence–secondary-structure pairs derived from Rfam covariance models, reference genomes, and Rfam full alignments. This repository accompanies the revised manuscript, *Improving RNA Secondary Structure Prediction Through Expanded Training Data*.

The complete dataset and large per-sequence model outputs are distributed through Zenodo. The DOI will be added here when the RNASSTR v2 record is published.

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
python scripts/shared/rfam_utils.py \
  --cm Rfam.cm \
  --clans Rfam.clanin \
  --out rfam_metadata.tsv
```

The exact large inputs and outputs used for the paper will be provided in the associated Zenodo record. See [`zenodo/README.md`](zenodo/README.md) for the planned archive contents.

## Dataset-generation workflow

Rfam v14.10 covariance models were searched against GTDB release 214 and NCBI RefSeq release 229 using Infernal v1.1.5 and family-specific gathering thresholds (`--cut_ga`). Rfam full-alignment sequences were incorporated to retain representation of families not recovered from the reference-genome searches. Hits were realigned to their family covariance models, converted to sequence–structure pairs, filtered using the criteria described in the manuscript, and assigned to structure-aware partitions.

The initial searches were run manually for each Rfam family rather than through a single pipeline script. For each family, the corresponding covariance model downloaded from the Rfam FTP site was searched against the appropriate reference-sequence database. A representative command was:

```bash
cmsearch \
  --cpu 32 \
  --cut_ga \
  --tblout results/RFxxxxx.tblout \
  -A results/RFxxxxx.sto \
  models/RFxxxxx.cm \
  databases/reference_sequences.fna \
  > results/RFxxxxx.cmsearch.log
```

`RFxxxxx.cm` and `reference_sequences.fna` were replaced with the family covariance model and corresponding sequence database for each search. Reported hits were subsequently filtered to retain E-values of 0.01 or less, together with the sequence and structural criteria described in the manuscript.

`filter_rnasstr_candidates.py` implements the quality-control procedure used to generate the released dataset. It removes exact sequence duplicates, sequences outside the family reference-length mean by more than two standard deviations, sequences more than two standard deviations below the family reference mean for annotated or canonical base-pair counts, and lower-ranked overlapping genomic hits. No phylogenetic filter is applied.

The filter accepts normalized candidate and Rfam-reference CSV files containing `id`, `family` or `rfam_id`, `sequence`, and either `structure` or `base_pairs`. Genomic overlap resolution additionally uses `accession`, `start`, `end`, and `evalue` when available:

```bash
python scripts/dataset/filter_rnasstr_candidates.py \
  --candidates merged_candidates.csv \
  --rfam-reference rfam_reference_sequences.csv \
  --output filtered_candidates.csv \
  --rejections qc_rejections.csv \
  --reference-summary rfam_qc_thresholds.csv
```

The final structural partitions were generated with:

```bash
python scripts/dataset/make_struct_clan_splits_v3.py \
  --sto-dir STO_dedupe \
  --rfam-cm Rfam.cm \
  --clanin Rfam.clanin \
  --outdir splits_struct_clan_v3_fast10 \
  --max-seqs-per-family 10 \
  --evalue 0.01 \
  --cpu 64 \
  --min-val-families 300 \
  --min-test-families 300 \
  --min-val-components 150 \
  --min-test-components 150
```

The family assignments were applied to the Stockholm files with:

```bash
python scripts/dataset/apply_family_splits_to_sto_with_family.py \
  --sto-dir STO_dedupe \
  --train-families splits_struct_clan_v3_fast10/train_families.txt \
  --val-families splits_struct_clan_v3_fast10/val_families.txt \
  --test-families splits_struct_clan_v3_fast10/test_families.txt \
  --outdir csv_splits_v3_fast10
```

## Citation

Citation metadata is provided in [`CITATION.cff`](CITATION.cff). Please cite both the manuscript and the versioned Zenodo dataset record once the DOI is available.

## Licenses

Code in this repository is released under the [MIT License](LICENSE). The RNASSTR data archive is intended for release under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); see [`DATA_LICENSE.md`](DATA_LICENSE.md).
