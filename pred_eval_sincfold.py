# Below Script gives a function for RNA secondary structure using sincFold, ouputting a csv with the 
# 
import pandas as pd
import subprocess
from subprocess import run, PIPE
from math import sqrt


# NOTE: CHANGE THIS TO YOUR OWN VALIDATION CSV LOCATION 
df = pd.read_csv("/data/RUDD/RUDD/2d-prediction/terry/curated_data_oct/rna_validate.csv")

def predict_structure_sincfold(sequence):
    try:
        # Run sincfold prediction command and capture all output
        result = run(f'sincFold pred {sequence}', shell=True, stdout=PIPE, text=True)
        prediction_output = result.stdout.strip()
        
        # Print the raw output for debugging purposes
        print(f"Full sincfold output for sequence '{sequence}':\n{prediction_output}")
        
        # Extract lines from the output
        lines = prediction_output.split('\n')
        
        # Find the sequence in the output and capture the next line (dot-bracket structure)
        for i, line in enumerate(lines):
            if line.strip() == sequence.strip():
                if i + 1 < len(lines):  # Ensure there's another line after the sequence
                    predicted_structure = lines[i + 1].strip()
                    print(f"Extracted output: {predicted_structure}")
                    return predicted_structure
        
        # If the structure isn't found, return a placeholder
        print("Extracted output: ...")
        return "..."
    except Exception as e:
        # Handle any errors and return a placeholder
        print(f"Error running sincfold: {e}")
        print("Extracted output: ...")
        return "..."

# Convert dot-bracket notation to base pairs
def dot_bracket_to_pairs(dot_bracket):
    stack = []
    pairs = []
    for i, char in enumerate(dot_bracket):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                pairs.append((stack.pop(), i))
    return pairs

# Calculate F1 score with tolerance on base pairs
def f1_shift(ref_bp, pred_bp):
    tp1 = sum(1 for rbp in ref_bp if rbp in pred_bp or
              (rbp[0], rbp[1] - 1) in pred_bp or
              (rbp[0], rbp[1] + 1) in pred_bp or
              (rbp[0] + 1, rbp[1]) in pred_bp or
              (rbp[0] - 1, rbp[1]) in pred_bp)
    tp2 = sum(1 for pbp in pred_bp if pbp in ref_bp or
              (pbp[0], pbp[1] - 1) in ref_bp or
              (pbp[0], pbp[1] + 1) in ref_bp or
              (pbp[0] + 1, pbp[1]) in ref_bp or
              (pbp[0] - 1, pbp[1]) in ref_bp)
    
    fn = len(ref_bp) - tp1
    fp = len(pred_bp) - tp1

    precision = tp2 / float(tp1 + fp) if (tp1 + fp) > 0 else 0
    recall = tp1 / float(tp1 + fn) if (tp1 + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    return f1

# Calculate MCC based on base pairs
def calculate_mcc_from_dot_bracket(pred_dot_bracket, true_dot_bracket):
    pred_pairs = dot_bracket_to_pairs(pred_dot_bracket)
    true_pairs = dot_bracket_to_pairs(true_dot_bracket)
    total_nucleotides = len(pred_dot_bracket)

    pred_pairs_set = set(pred_pairs)
    true_pairs_set = set(true_pairs)

    TP = len(pred_pairs_set.intersection(true_pairs_set))
    FP = len(pred_pairs_set) - TP
    FN = len(true_pairs_set) - TP

    total_possible_pairs = total_nucleotides * (total_nucleotides - 1) / 2
    TN = total_possible_pairs - len(true_pairs) - FP - FN

    mcc_numerator = TP * TN - FP * FN
    mcc_denominator = sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

    return mcc_numerator / mcc_denominator if mcc_denominator != 0 else 0

# Loop over the validation set, make predictions, and calculate F1 and MCC scores

# < Prediction Loop >
predictions = []
f1_scores = []
mcc_scores = []

for index, row in df.iterrows():
    sequence = row['sequence']
    true_dot_bracket = row['structure']

    # Predict structure using sincfold
    predicted_dot_bracket = predict_structure_sincfold(sequence)
    predictions.append(predicted_dot_bracket)
    
    # Calculate F1 score
    ref_pairs = dot_bracket_to_pairs(true_dot_bracket)
    pred_pairs = dot_bracket_to_pairs(predicted_dot_bracket)
    f1 = f1_shift(ref_pairs, pred_pairs)
    f1_scores.append(f1)
    
    # Calculate MCC score
    mcc = calculate_mcc_from_dot_bracket(predicted_dot_bracket, true_dot_bracket)
    mcc_scores.append(mcc)
    
    print(f"Sequence {index+1}: F1 Score = {f1:.3f}, MCC = {mcc:.3f}")

# Below --- Save results with filename including F1 and MCC ---
# --- Summary and Output ---
mean_f1 = sum(f1_scores) / len(f1_scores)
mean_mcc = sum(mcc_scores) / len(mcc_scores)

formatted_f1 = f"{mean_f1:.3f}"
formatted_mcc = f"{mean_mcc:.3f}"

# NOTE: CHANGE THIS TO YOUR OWN experiment model_label and output_path
# This is where the final result CSV will be saved with model and scores embedded

# Define file label or name
model_label = "epoch10_Trial2_sincfold"

# Construct the full output path
output_path = f"/data/RUDD/RUDD/2d-prediction/terry/curated_data_oct/sincfold/{model_label}_f1_{formatted_f1}_mcc_{formatted_mcc}.csv"

# Save to CSV
df.to_csv(output_path, index=False)
print(f"\n Results saved to: {output_path}")