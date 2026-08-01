# RNASSTR

RNASSTR is a large, structure-aware collection of RNA sequence–secondary-structure pairs derived from Rfam covariance models, reference genomes, and Rfam full alignments. This repository accompanies the revised manuscript, *Improving RNA Secondary Structure Prediction Through Expanded Training Data*.

The complete dataset and large per-sequence model outputs are distributed through Zenodo. The DOI will be added here when the RNASSTR v2 record is published.

## Contents

- `scripts/shared/rfam_utils.py`: parse Rfam covariance-model metadata and family-to-clan assignments.
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

## Citation

Citation metadata is provided in [`CITATION.cff`](CITATION.cff). Please cite both the manuscript and the versioned Zenodo dataset record once the DOI is available.

## Licenses

Code in this repository is released under the [MIT License](LICENSE). The RNASSTR data archive is intended for release under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); see [`DATA_LICENSE.md`](DATA_LICENSE.md).
