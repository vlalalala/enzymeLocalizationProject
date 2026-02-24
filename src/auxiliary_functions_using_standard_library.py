"""
Functions (to use e.g. for cleaning config inputs) that only require Python's standard library
"""
import os
import re
import ast
import json
import pickle
import glob
import csv
import math

def is_int_value(value):
    return float(int(value)) == float(value)

def dump_json_base(dump_directory: str, file_basename: str, dict_to_dump: dict):
    """
    Dump a dictionary (possibly nested) as JSON. 
    No type conversions.
    Set file_basename without the .json ending.
    """
    os.makedirs(dump_directory, exist_ok=True)
    path = os.path.join(dump_directory, f"{file_basename}.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(dict_to_dump, f, indent=4)

def int_from_sci(x):
    """Returns the integer out of scientific writing
    (e.g. 1e2 is interpreted as a float)
    """
    return int(float(x))

def all_non_negative(dictionary):
    """Checks in a nested dictionary whether all values are non-negative.
    """
    for v in dictionary.values():
        if isinstance(v, dict):
            if not all_non_negative(v):
                return False
        elif v < 0:
            return False
    return True

def format_sci(x: float) -> str:
    """Returns a string, which is a positive float shown through
    scientific notation with 1 decimal and 2-digit exponent.
    """
    if not math.isfinite(x):
        return str(x)
    s = f"{x:.1e}"           # e.g. "3.6e-16"
    base, exp = s.split("e")
    exp_num = int(exp)
    return f"{base}e{'+' if exp_num >= 0 else '-'}{abs(exp_num):02d}"

def load_json(path):
    """path should have the .json extension.
    in json files, the keys are converted to strings always.
    When loading, the keys are converted to integers, if possible.
    """
    _, file_extension = os.path.splitext(os.path.basename(path))
    if not os.path.isfile(path) or file_extension != ".json":
        raise ValueError(f"The file {path} does not exist or is not a .json file.")
    with open(path, "r", encoding="utf-8") as f:
        contents = json.load(f)

    def convert_keys(obj):
        if isinstance(obj, dict):
            new_dict = {}
            for k, v in obj.items():
                # Try converting the key to int
                try:
                    new_key = int(k)
                except (ValueError, TypeError):
                    new_key = k
                # Recursively convert nested dictionaries
                new_dict[new_key] = convert_keys(v)
            return new_dict
        if isinstance(obj, list):
            return [convert_keys(i) for i in obj]
        return obj

    return convert_keys(contents)

def pickle_dump_binary(path, variable) -> None:
    """Saves a variable in a path.
    """
    with open(path, 'wb') as f:
        pickle.dump(variable, f)

def pickle_load_binary(path):
    """Returns a pickled file.
    """
    with open(path, 'rb') as f:
        loaded_variable = pickle.load(f)
    return loaded_variable

def closest_value(my_list, target):
    """Returns the closest value within my_list to the target.
    """
    return min(my_list, key=lambda x: abs(x - target))

def as_list(value, type_cast=str):
    """
    Converts various types of input into a list of elements of type `type_cast`.

    Handles:
    - Single values (e.g., 5 or "abc")
    - Space-separated strings (e.g., "1 2 3")
    - Real lists (e.g., [1, 2])
    - Strings that represent Python lists (e.g., "[1, 2]")

    Raises:
        ValueError: If a string looks like a list (contains '[' or ']') but is invalid
    """
    if isinstance(value, list):
        return [type_cast(v) for v in value]

    elif isinstance(value, str):
        val = value.strip()

        # Case: string that looks like a list (e.g. "[1, 2]")
        if "[" in val or "]" in val:
            try:
                # Try JSON parsing first
                parsed = json.loads(val)
                if not isinstance(parsed, list):
                    raise ValueError(f"Expected list in string, got: {val}")
                return [type_cast(v) for v in parsed]
            except json.JSONDecodeError: # Try JSON parsing first
                # Fallback: Try Python's literal_eval
                try:
                    parsed = ast.literal_eval(val)
                    if not isinstance(parsed, list):
                        raise ValueError(f"Expected list-like string, got: {val}")
                    return [type_cast(v) for v in parsed]
                except Exception as e:
                    raise ValueError(f"Invalid list syntax in string: {val!r}") from e
        
        # Case: space-separated string
        return [type_cast(v) for v in val.split()]

    else:
        # Single value (e.g. int or float)
        return [type_cast(value)]



def find_sorted_unique_files_with_max_digits_and_max_value(folder, pattern_to_find, max_iteration_value):
    """
    Returns a tuple:
        (unique_files_sorted, max_digits)

    - unique_files_sorted: list of files sorted numerically by iteration number,
      keeping only one file per iteration number (the one with the most digits / leading zeros),
      and only including iterations <= max_iteration_value.
    - max_digits: the highest number of digits in the iteration numbers
    """
    pattern = os.path.join(folder, pattern_to_find)
    files = glob.glob(pattern)
    
    iter_map = {}  # iteration number -> list of files
    max_digits = 0

    for f in files:
        basename = os.path.basename(f)
        m = re.search(r"(\d+)", basename)
        if not m:
            continue

        iter_num_str = m.group(1)
        iter_num = int(iter_num_str)

        # skip anything above the requested max iteration value
        if iter_num > max_iteration_value:
            continue

        max_digits = max(max_digits, len(iter_num_str))
        iter_map.setdefault(iter_num, []).append(f)

    unique_files = []
    for iter_num, flist in iter_map.items():
        # Sort by number of digits (descending), keep the most zero-padded one
        flist.sort(
            key=lambda x: len(re.search(r"(\d+)", os.path.basename(x)).group(1)),
            reverse=True
        )
        best_file = flist[0]
        unique_files.append(best_file)

        # Delete duplicates (keep only the most zero-padded one)
        for duplicate in flist[1:]:
            try:
                os.remove(duplicate)
                print(f"Deleted duplicate: {duplicate}")
            except OSError:
                print(f"Warning: could not delete {duplicate}")

    # Sort the unique files numerically
    unique_files.sort(
        key=lambda x: int(re.search(r"(\d+)", os.path.basename(x)).group(1))
    )

    return unique_files, max_digits

def find_max_in_nested_dict(dictionary):
    """Finds the maximum value within a nested dictionary.
    """
    max_val = float('-inf')
    for v in dictionary.values():
        if isinstance(v, dict):
            # Recurse into nested dictionary
            max_val = max(max_val, find_max_in_nested_dict(v))
        elif isinstance(v, (int, float)):
            max_val = max(max_val, v)
    return max_val



class CSVLogger:
    """Class written by ChatGPT."""
    def __init__(self, path):
        self.path = path
        self._header_written = os.path.exists(path)

    def log(self, iteration: int, values: dict):
        """
        Append one row to the CSV log.

        iteration: int
        values: dict[str, scalar]
        """
        fieldnames = ["iteration", *values.keys()]

        with open(self.path, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)

            if not self._header_written:
                writer.writeheader()
                self._header_written = True

            row = {"iteration": iteration, **values}
            writer.writerow(row)