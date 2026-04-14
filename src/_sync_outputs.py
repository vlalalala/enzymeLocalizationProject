import os
import glob
import shutil
import filecmp
from pathlib import Path
from tqdm import tqdm

"""
This file was mostly copy-pasted from Claude.
"""

def get_user_input_files(folder):
    """Return only the input files (not output files) within a combined_* folder.
    """
    input_filenames = [
        "parameters_solver_input.yaml",
        "parameters_solver_output.yaml",
        "parameters_geometry.yaml",
        "parameters_discretization.yaml",
        "parameters_value_conditions.yaml",
        "parameters_optimization.yaml",
        "species.csv",
        "spontaneous_reactions.csv",
        "enzymes.csv",
        "enzymatic_reactions.csv"
        # add all your input filenames here
    ]
    return {f: os.path.join(folder, f) for f in input_filenames
            if os.path.exists(os.path.join(folder, f))}

def folders_have_same_inputs(folder_a, folder_b):
    """Check if all input files in folder_a are identical to those in folder_b."""
    inputs_a = get_user_input_files(folder_a)
    inputs_b = get_user_input_files(folder_b)

    # Must have the same input files present
    if set(inputs_a.keys()) != set(inputs_b.keys()):
        return False

    # Compare contents of each input file
    for filename in inputs_a:
        if not filecmp.cmp(inputs_a[filename], inputs_b[filename], shallow=False):
            return False

    return True

def get_output_files(folder, input_filenames):
    """Return all files in folder that are not input files."""
    all_files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]
    return [f for f in all_files if f not in input_filenames]

def copy_outputs_from_b_to_a(folder_a, folder_b, input_filenames):
    """Copy output files and folders from folder_b to folder_a."""
    all_items = os.listdir(folder_b)
    
    for item in all_items:
        # Skip input files
        if item in input_filenames:
            continue
        
        src = os.path.join(folder_b, item)
        dst = os.path.join(folder_a, item)
        
        if os.path.isfile(src):
            if not os.path.exists(dst):
                shutil.copy2(src, dst)
                print(f"  Copied file: {item}")
            else:
                print(f"  Skipped file (already exists): {item}")
        
        elif os.path.isdir(src):
            if not os.path.exists(dst):
                shutil.copytree(src, dst)
                print(f"  Copied folder: {item}")
            else:
                print(f"  Skipped folder (already exists): {item}")

def sync_outputs(folder_a, folder_b):
    """
    For each combined_* folder in A, find a matching combined_* folder in B
    (same input files) and copy its outputs to A.
    """
    combined_a = sorted(glob.glob(os.path.join(folder_a, "combined_*")))
    combined_b = sorted(glob.glob(os.path.join(folder_b, "combined_*")))

    input_filenames = set(get_user_input_files(combined_a[0]).keys()) if combined_a else set()

    for folder_a_sub in tqdm(combined_a, desc="Syncing folders"):
        print(f"\nChecking {folder_a_sub}...")
        matched = False
        for folder_b_sub in combined_b:
            if folders_have_same_inputs(folder_a_sub, folder_b_sub):
                print(f"  Match found: {folder_b_sub}")
                copy_outputs_from_b_to_a(folder_a_sub, folder_b_sub, input_filenames)
                matched = True
                break
        if not matched:
            print(f"  No match found in {folder_b}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Sync outputs from folder B to folder A based on matching input files.")
    parser.add_argument("--folder_a", help="Folder A (destination)")
    parser.add_argument("--folder_b", help="Folder B (source of outputs)")
    args = parser.parse_args()

    sync_outputs(args.folder_a, args.folder_b)