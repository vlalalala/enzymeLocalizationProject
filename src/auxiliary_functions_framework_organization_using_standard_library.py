import os
import re
from pathlib import Path

def get_concentrations_files_within_folder(folder):
    """Returns a list with the complete paths of files that correspond to
    concentrations (or an empty list, if none of these files exist).

    Compatible with Snakemake wildcards.
    """
    # Find all concentration files matching pattern
    files = [f for f in os.listdir(folder) if re.match(r".iteration_nr_\d+_concentrations", f)]
    return files

def rename_iteration_files(folder, min_digits) -> int:
    """
    Renames iteration files (e.g. .iteration_nr_3_concentrations.json → .iteration_nr_00003_concentrations.json)
    so that *all* files have a consistent number of digits in their iteration number.
    Returns the final number of digits.

    Rules:
    - Only adds zeros (never removes).
    - The final width is max(min_digits, max_digits_found_in_folder).
    - Ensures all files use consistent zero-padding.

    Parameters
    ----------
    folder : str
        Folder containing the files.
    min_digits : int
        Minimum desired number of digits.
    """
    files = get_concentrations_files_within_folder(folder)
    folder = Path(folder)

    if not files:
        print(f"No files found in {folder} matching iterations pattern.")
        return min_digits

    # Ensure full paths
    files = [Path(f) if os.path.isabs(f) else folder / f for f in files]

    # Detect maximum existing digit count in filenames
    max_found_digits = 0
    for f in files:
        match = re.search(r"(\d+)", os.path.basename(f))
        if match:
            max_found_digits = max(max_found_digits, len(match.group(1)))

    # Determine final padding width
    final_digits = max(min_digits, max_found_digits)

    for f in files:
        basename = os.path.basename(f)
        match = re.search(r"(\d+)", basename)
        if not match:
            continue
        old_num_str = match.group(1)
        num = int(old_num_str)
        old_digits = len(old_num_str)
        # Skip if already correct
        if old_digits == final_digits:
            continue
        # Create new padded string (adding zeros if needed)
        new_num_str = f"{num:0{final_digits}d}"
        new_basename = re.sub(r"(\d+)", new_num_str, basename, count=1)
        new_path = os.path.join(folder, new_basename)
        if f != new_path:
            old_path = os.path.join(folder, basename)   # <-- add this
            os.rename(old_path, new_path)
            print(f"Renamed: {basename} → {new_basename}")
    return final_digits


def find_latest_solution(folder):
    """
    Return the path to the latest {folder}/iteration_nr_XX_concentrations.json file.
    If there are no concentrations files, returns None. 
    (which is created when the file structure is first created.)
    Handles zero-padded iteration numbers (e.g. 00, 001, etc.).
    """
    files = get_concentrations_files_within_folder(folder)
    if not files:
        return None
    # Extract iteration numbers as integers
    iter_nums = [int(re.search(r"(\d+)", f).group(1)) for f in files]
    latest_iter = max(iter_nums)
    # Keep original zero-padding width
    width = len(re.search(r"(\d+)", files[0]).group(1))
    return os.path.join(folder, f".iteration_nr_{latest_iter:0{width}d}_concentrations.json")
