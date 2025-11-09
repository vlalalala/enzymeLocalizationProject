import os
import re
import glob
from auxiliary_functions_using_standard_library import load_json

config_info = load_json("config/config.json")
solver_data_folder = config_info["solver_data_folder"]

def get_solver_data_folder(wildcards):
    """Returns the path of the folder in which the iteration data is saved."""
    results_folder = f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/{solver_data_folder}"
    return results_folder

def get_concentrations_files_within_folder(folder):
    """Returns a list with the complete paths of files that correspond to
    concentrations (or an empty list, if none of these files exist).

    Compatible with Snakemake wildcards.
    """
    # Find all concentration files matching pattern
    files = [f for f in os.listdir(folder) if re.match(r"\.iteration_nr_\d+_concentrations", f)]
    return files

def rename_iteration_files(folder, min_digits) -> None:
    """
    Renames iteration files (e.g. .iteration_nr_3_concentrations.json → .iteration_nr_00003_concentrations.json)
    so that *all* files have a consistent number of digits in their iteration number.

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

    if not files:
        print(f"No files found in {folder} matching iterations pattern.")

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
            os.rename(f, new_path)
            print(f"Renamed: {basename} → {new_basename}")


def find_latest_solution(wildcards):
    """
    Return the path to the latest {results_folder}/iteration_nr_XX_concentrations.json file.
    If there are no concentrations files, returns the path to the .initialized file 
    (which is created when the file structure is first created.)
    Handles zero-padded iteration numbers (e.g. 00, 001, etc.).

    Compatible with Snakemake wildcards.
    """
    results_folder = get_solver_data_folder(wildcards)
    files = get_concentrations_files_within_folder(results_folder)
    if not files:
        return "none"
    # Extract iteration numbers as integers
    iter_nums = [int(re.search(r"(\d+)", f).group(1)) for f in files]
    latest_iter = max(iter_nums)
    # Keep original zero-padding width
    width = len(re.search(r"(\d+)", files[0]).group(1))
    return os.path.join(results_folder, f".iteration_nr_{latest_iter:0{width}d}_concentrations.json")

def find_iteration_to_plot(wildcards, iteration=None):
    """
    Return the path to the requested or latest iteration file:
        If config["iteration"] is given, return that iteration’s file.
        Otherwise, return the latest iteration file.
    Important: all iteration files should have the same number of digits!
    (if unsure, run rename_iteration_files)
    Raises an error if there are no iterations to plot.
    
    Compatible with Snakemake wildcards.
    """
    results_folder = get_solver_data_folder(wildcards)
    files = get_concentrations_files_within_folder(results_folder)
    # Just in case, rename files so that all files have the same number of digits
    rename_iteration_files(results_folder, min_digits=0)
    width = len(re.search(r"(\d+)", files[0]).group(1))
    latest_solution = find_latest_solution(wildcards)
    if ".initialized" in os.path.basename(latest_solution):
        raise FileNotFoundError("There are no iterations from which to plot anything.")

    # Extract numeric iteration identifiers
    iter_nums = [int(re.search(r"(\d+)", f).group(1)) for f in files]
    if iteration is not None:
        if int(iteration) not in iter_nums:
            raise FileNotFoundError(f"Requested iteration {iter_num} not found in {results_folder}")
        # Specific iteration requested
        iter_num = int(iteration)
        fname = f".iteration_nr_{iter_num:0{width}d}_concentrations.json"
        return os.path.join(results_folder, fname)      
    return latest_solution

def define_plot_output(base_folder, target_iteration=None):
    """
    Defines output filename for the plot.
    If target_iteration is given, use that number.
    If not, detect the latest iteration file and use its number.
    """
    iteration_folder = os.path.join(base_folder, "solver_iteration_data")
    pattern = os.path.join(iteration_folder, ".iteration_nr_*_concentration.json")
    files = glob.glob(pattern)

    if not files:
        raise FileNotFoundError(f"No iteration files found in {iteration_folder}")

    def get_iter_num(f):
        m = re.search(r"(\d+)", os.path.basename(f))
        return int(m.group(1)) if m else -1

    if target_iteration is not None:
        iter_num = int(float(target_iteration))
    else:
        latest = max(files, key=get_iter_num)
        iter_num = get_iter_num(latest)

    # Preserve zero-padding width from any file
    example_file = files[0]
    width = len(re.search(r"(\d+)", os.path.basename(example_file)).group(1))

    return os.path.join(base_folder, f".plot_iteration_{iter_num:0{width}d}.png")