# RNASSTR: RNA Secondary Structure Repository

**RNASSTR** is a curated benchmark suite and reproducible pipeline for evaluating RNA secondary structure prediction models. It provides standardized datasets, conversion utilities, and evaluation scripts to facilitate consistent benchmarking across various models, including [sincFold](https://github.com/sinc-lab/sincFold) and [MXfold2](https://github.com/mxfold/mxfold2).

If you use RNASSTR in your research, please cite the following paper:

> @article{langeberg2025improving,  
> title={Improving RNA Secondary Structure Prediction Through Expanded Training Data},  
> author={Langeberg, Conner J and Kim, Taehan and Nagle, Roma and Meredith, Charlotte and Garuadapuri, Dimple Amitha and Doudna, Jennifer A and Cate, Jamie HD},  
> journal={bioRxiv},  
> pages={2025--05},  
> year={2025},  
> publisher={Cold Spring Harbor Laboratory}  
> }

---

## Table of Contents

- [About](#about)
- [Installation](#installation)
- [Format Setup for Prediction and Evaluation](#format-setup-for-prediction-and-evaluation)
- [MXfold2 Setup](#mxfold2-setup)
- [sincFold Evaluation](#sincfold-evaluation)
- [Contact](#contact)

---

## About

In recent years, deep learning has revolutionized protein structure prediction, achieving remarkable speed and accuracy. RNA structure prediction, however, has lagged behind. Although several methods have shown moderate success in predicting RNA secondary and tertiary structures, none have reached the accuracy observed with contemporary protein models. The lack of success of these RNA structure prediction models has been proposed to be due to limited high-quality structural information that can be used as training data. To probe
this proposed limitation, we developed a large and diverse dataset comprising paired RNA sequences and their corresponding secondary structures.

### Dataset Format

Each sample in our dataset is stored in CSV files with the following columns:
- `id`: Unique identifier of the RNA sequence.
- `sequence`: The RNA nucleotide sequence.
- `structure`: The dot-bracket notation of the RNA structure.
- `base_pairs`: A list of paired indices in square brackets (e.g., `[[0, 9], [1, 8]]`).
- `len`: The length of the sequence.

Datasets are split into:
- `train/`
- `validate/`
- `test/`

---

## Installation

Each model may require a different environment setup. See below for details.

### For sincFold

sincFold's environment and dependencies are self-contained. Please refer directly to the [sincFold repository](https://github.com/sinc-lab/sincFold) for installation. If using RNASSTR for evaluation with sincFold:
- You must modify the weights in `model.py` of sincFold to switch between training checkpoints.
- The script `pred_eval_sincfold.py` will use the CLI version to predict structures.

### For MXfold2

MXfold2 has a more complex setup and requires specific formats and dependencies. This repo includes the following files for reproducibility:
- `mxfold2_requirements.txt` listing key Python dependencies
- `mxfold2_install_dependencies.sh` to install MXfold2 from a prebuilt wheel (for Linux or Mac)

Use this setup:

```bash
# Clone the repository
git clone https://github.com/romanagle/RNASSTR.git
cd RNASSTR

# Create and activate a conda environment
conda create --name mxfold2_env python=3.10
conda activate mxfold2_env

# Install dependencies
pip install -r setting/mxfold2_requirements.txt

# Run custom installation script for MXfold2 wheel
bash setting/mxfold2_install_dependencies.sh
```

---

## Format Setup for Prediction and Evaluation

Different models require different input formats. The following steps show how to convert and prepare input data appropriately:

### For MXfold2 and Ufold (.lst with .bpseq)

These models require the input in `.bpseq` format, and `.lst` files that list paths to these `.bpseq` files.

#### Step 1: Convert CSV to BPSEQ

Use `curate_bpseq.py` with one of the dataset splits (e.g., `rna_train.csv`, `rna_test.csv`, or `rna_validate.csv`):

```bash
python curate_bpseq.py --input_csv data/rna_train.csv --output_dir data/bpseq/train
```

This script is designed with argparse and looks like:

```python
def main():
    parser = argparse.ArgumentParser(description="Generate BPSEQ files from a Dataset CSV.")
    parser.add_argument('-i', '--input_csv', required=True, help="Path to the input CSV file (e.g., combined_output.csv).")
    parser.add_argument('-o', '--output_dir', required=True, help="Path to the output directory where BPSEQ files will be saved.")
    args = parser.parse_args()
    process_csv(args.input_csv, args.output_dir)
```

#### Step 2: Create .lst File

Once BPSEQ files are created, generate `.lst` files required by MXfold2:

```bash
python curate_lst.py --bpseq_dir data/bpseq/train --output_lst data/train.lst
```

> 🔧 Note: For `curate_lst.py` and `pred_eval_sincfold.py`, you may need to modify the path variables directly in the script unless you add argparse handling.

---

### Script Summary

| Script                  | Use Case                                                                 |
|------------------------|--------------------------------------------------------------------------|
| `curate_bpseq.py`       | Converts RNA CSV to `.bpseq` format for models like MXfold2 or Ufold.     |
| `curate_lst.py`         | Creates `.lst` file that lists `.bpseq` files for MXfold2 input.         |
| `pred_eval_sincfold.py` | Runs sincFold predictions and evaluates them using F1 and MCC metrics.   |

---

## MXfold2 Setup

MXfold2 requires `.lst` files containing paths to `.bpseq` files for both training and evaluation.

We include wheel files for ease of installation:
- `setting/mxfold2-0.1.2-cp310-cp310-manylinux_2_17_x86_64.whl` (Linux)
- `setting/mxfold2-0.1.2-cp310-cp310-macosx_13_0_arm64.whl` (MacOS)

Refer to the official repo for details: [https://github.com/mxfold/mxfold2](https://github.com/mxfold/mxfold2)

---

## sincFold Evaluation

To evaluate sincFold using RNASSTR benchmark datasets:

```bash
python pred_eval_sincfold.py --input_csv data/rna_validate.csv --output_csv results/sincfold_predictions.csv
```

This script performs the following:
- Runs sincFold structure prediction using the CLI.
- Parses the predicted dot-bracket strings.
- Computes F1 score (with tolerance) and Matthews Correlation Coefficient (MCC) for each sequence.
- Taking input of the original csv, writes results to a new csv with a column containing the Predicted structure

> 📌 **Note**: Be sure to modify `sincfold/model.py` to point to your desired checkpoint weights before running prediction. This ensures you're using the correct trained model during evaluation, and run pip install again to ensure update.

---

## Contact

For inquiries, open a GitHub issue or reach out directly. For script-specific help, check the docstrings or headers in each script.