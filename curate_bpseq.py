import csv
import ast
import re
import os
import argparse

# If facing trouble with this script or have inquiries:
# feel free to make the github issue or contact through github issue or: conner.langeberg@berkeley.edu, terry.kim@berkeley.edu

def sanitize_filename(name):
    """Sanitize the sequence name to create a valid filename."""
    # Remove any characters that are not alphanumeric, underscore, hyphen, dot, or space
    sanitized_name = re.sub(r'[^\w\-_\. ]', '_', name)
    return sanitized_name

def parse_numerical_pairing(pairing_str):
    """Parse the numerical pairing string into a dictionary of pairings."""
    try:
        # Remove surrounding quotes if present
        pairing_str = pairing_str.strip('\"')
        # Use ast.literal_eval to safely parse the string into a list
        pairing_list = ast.literal_eval(pairing_str)
        pairing_dict = {}
        for i, j in pairing_list:
            pairing_dict[i] = j
            pairing_dict[j] = i
        return pairing_dict
    except (SyntaxError, ValueError) as e:
        print(f"Error parsing pairing string: {pairing_str}\n{e}")
        return {}

def process_csv(input_csv, output_dir):
    """Read the Dataset CSV and generate BPSEQ files."""
    if not os.path.isfile(input_csv):
        print(f"Error: The input file '{input_csv}' does not exist.")
        return

    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
            print(f"Created output directory: {output_dir}")
        except OSError as e:
            print(f"Error creating output directory '{output_dir}': {e}")
            return

    with open(input_csv, 'r', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader, None)  # Skip header if present
        for row_num, row in enumerate(reader, start=2):  # Start at 2 considering header
            # Unpack the columns
            if len(row) < 5:
                print(f"Skipping incomplete row {row_num}: {row}")
                continue
            # Note: our dataset has columns id,sequence,structure,base_pairs,len
            name, sequence, dot_bracket, numerical_pairing_str, length = row[:5]
            # Sanitize the name to create a valid filename
            filename = sanitize_filename(name) + '.bpseq'
            filepath = os.path.join(output_dir, filename)
            # Parse the numerical pairing
            pairing_dict = parse_numerical_pairing(numerical_pairing_str)
            # Initialize pairing indices
            sequence_length = len(sequence)
            pairing_indices = [0] * sequence_length
            # Update pairing indices based on the pairing dictionary
            for i in range(1, sequence_length + 1):
                if i in pairing_dict:
                    pairing_indices[i - 1] = pairing_dict[i]
            # Write the BPSEQ file
            try:
                with open(filepath, 'w', encoding='utf-8') as bpseq_file:
                    for index, (nucleotide, pairing_index) in enumerate(zip(sequence, pairing_indices), start=1):
                        bpseq_file.write(f"{index} {nucleotide} {pairing_index}\n")
                print(f"Generated BPSEQ file: {filepath}")
            except IOError as e:
                print(f"Error writing to file '{filepath}': {e}")

def main():
    parser = argparse.ArgumentParser(description="Generate BPSEQ files from a Dataset CSV.")
    parser.add_argument(
        '-i', '--input_csv',
        required=True,
        help="Path to the input CSV file (e.g., combined_output.csv)."
    )
    parser.add_argument(
        '-o', '--output_dir',
        required=True,
        help="Path to the output directory where BPSEQ files will be saved."
    )
    args = parser.parse_args()

    input_csv = args.input_csv
    output_dir = args.output_dir

    print(f"Input CSV File: {input_csv}")
    print(f"Output Directory: {output_dir}")

    process_csv(input_csv, output_dir)
    print("BPSEQ file generation completed.")

if __name__ == "__main__":
    main()

