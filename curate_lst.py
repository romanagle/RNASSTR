import os

def compile_bpseq_paths(directory, output_file, base_path=None, limit=None):
    """
    Compile paths to .bpseq files into a .lst file, prepending base_path and the directory type (train/test) to each file.
    
    Parameters:
    - directory: Directory containing .bpseq files (train or test folder).
    - output_file: Path to the .lst file to save the file paths.
    - base_path: The base directory to prepend to each .bpseq file path. If None, use the current working directory.
    - limit: Optional, limits the number of .bpseq files to include.
    
    Returns:
    - A list of file paths that were written to the .lst file.
    """
    # If base_path is None, use the current working directory
    if base_path is None:
        base_path = os.getcwd()

    # Ensure base_path doesn't end with a slash
    base_path = base_path.rstrip('/')

    # Determine if we are processing train or test data
    if 'train' in directory:
        sub_folder = 'train'
    elif 'test' in directory:
        sub_folder = 'test'
    elif 'validate' in directory:
        sub_folder = 'validate'
    else:
        sub_folder = ''  # Default if no train/test distinction

    # Initialize an empty list to store the file paths
    bpseq_paths = []

    # Walk through the directory to find all .bpseq files
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".bpseq"):
                bpseq_relative_path = os.path.relpath(os.path.join(root, file), start=directory)
                # Prepend the base path and the sub-folder to each file
                full_bpseq_path = os.path.join(base_path, sub_folder, bpseq_relative_path)
                bpseq_paths.append(full_bpseq_path)
    
    # If a limit is specified, slice the list
    if limit is not None:
        bpseq_paths = bpseq_paths[:limit]
    
    # Write the paths to the output .lst file
    with open(output_file, 'w') as f:
        for path in bpseq_paths:
            f.write(f"{path}\n")
    
    print(f"Compiled {len(bpseq_paths)} .bpseq file paths into {output_file}")
    
    # Return the list of paths that were written to the file
    return bpseq_paths

# Usage:
if __name__ == "__main__":
    # NOTE: Define your data directories (ensure they contain only .bpseq files or subfolders of them)
    # REPLACE: the initial path defined below as you need
    test_dir = '/home/exx/terry/curated_bpseq/test'
    train_dir = '/home/exx/terry/curated_bpseq/train'
    validate_dir = '/home/exx/terry/curated_bpseq/validate'
    train_output_file = '/home/exx/terry/curated_bpseq/train.lst'
    test_output_file = '/home/exx/terry/curated_bpseq/test.lst'
    validate_output_file  = '/home/exx/terry/curated_bpseq/validate.lst'
    # Base path that will be prepended to each .bpseq file path inside the .lst (must match how MXfold2 expects paths)
    base_path = '/home/exx/terry/curated_bpseq'
    
    # OPTIONAL: Compile paths with a limit of 100 files for testing purposes, prepending the current directory
    #changed_train_files = compile_bpseq_paths(train_dir, train_output_file, base_path=base_path, limit=100)
    #changed_test_files = compile_bpseq_paths(test_dir, test_output_file, base_path=base_path, limit=100)
    
    changed_train_files = compile_bpseq_paths(train_dir, train_output_file, base_path=base_path)
    changed_test_files = compile_bpseq_paths(test_dir, test_output_file, base_path=base_path)
    changed_validate_files = compile_bpseq_paths(validate_dir, validate_output_file, base_path=base_path)
    # Print out the changed files
    print("Train files that were written to the .lst file:")
    for file in changed_train_files:
        print(file)
    
    print("\nTest files that were written to the .lst file:")
    for file in changed_test_files:
        print(file)
    print("\nValidation files that were written to the .lst file:")
    for file in changed_validate_files:
        print(file)
