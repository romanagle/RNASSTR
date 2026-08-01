# RNASSTR

RNASSTR is a large, structure-aware collection of RNA sequence–secondary-structure pairs derived from Rfam covariance models, reference genomes, and Rfam full alignments. This repository contains the reproducible scoring and figure-generation code used for the revised manuscript, *Improving RNA Secondary Structure Prediction Through Expanded Training Data*.

The complete dataset and large per-sequence model outputs are distributed through Zenodo. The DOI will be added here when the RNASSTR v2 record is published.

## Contents

- `scripts/scoring/`: normalize and score published/retrained SincFold and Lyra-TransPred outputs.
- `scripts/rnafold/`: create the deterministic family-stratified RNAfold subset, score RNAfold output, and analyze predicted-pair density.
- `scripts/figures/`: generate the revised dataset and model-comparison figures.
- `scripts/shared/`: common Rfam, feature, and plotting helpers.
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

MXFold2 and Lyra-UFold are not part of the revised model comparison.

## Installation

Create a Python 3.10 or later environment and install the analysis dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

RNAfold must be installed separately through ViennaRNA for the minimum-free-energy benchmark. SincFold and Lyra model inference require their respective model implementations; this repository scores their exported predictions rather than redistributing third-party source code.

## Reproducing the analyses

All scoring scripts provide command-line documentation:

```bash
python scripts/scoring/score_sincfold_models.py --help
python scripts/scoring/score_lyra_models.py --help
python scripts/rnafold/make_rnafold_subsample.py --help
python scripts/rnafold/score_rnafold_output.py --help
python scripts/rnafold/analyze_rnafold_sincfold_pair_bias.py --help
```

The scoring convention treats each possible unordered nucleotide pair as a binary candidate. F1 uses exact matching base pairs, and MCC includes true negatives among all possible unordered pairs.

The exact large inputs and outputs used for the paper will be provided in the associated Zenodo record. See [`zenodo/README.md`](zenodo/README.md) for the planned archive contents.

## Dataset-generation workflow

Rfam v14.10 covariance models were searched against GTDB release 214 and NCBI RefSeq release 229 using Infernal v1.1.5 and family-specific gathering thresholds (`--cut_ga`). Rfam full-alignment sequences were incorporated to retain representation of families not recovered from the reference-genome searches. Hits were realigned to their family covariance models, converted to sequence–structure pairs, filtered using the criteria described in the manuscript, and assigned to structure-aware partitions. The exact search and filtering scripts are not yet present in this branch; they must be added before the repository can be described as reproducing dataset generation end to end.

## Citation

Citation metadata is provided in [`CITATION.cff`](CITATION.cff). Please cite both the manuscript and the versioned Zenodo dataset record once the DOI is available.

## Licenses

Code in this repository is released under the [MIT License](LICENSE). The RNASSTR data archive is intended for release under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); see [`DATA_LICENSE.md`](DATA_LICENSE.md).
