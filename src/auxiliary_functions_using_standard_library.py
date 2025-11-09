"""
Functions (to use e.g. for cleaning config inputs) that only require Python's standard library
"""
import os
import re
import ast
import json
import pickle
import glob

def int_from_sci(x):
    return int(float(x))

def nested_max(dictionary):
    max_val = float("-inf")
    for v in dictionary.values():
        if isinstance(v, dict):
            max_val = max(max_val, nested_max(v))
        else:
            max_val = max(max_val, v)
    return max_val

def all_non_negative(dictionary):
    """Checks in a nested dictionary whether all values are non-negative"""
    for v in dictionary.values():
        if isinstance(v, dict):
            if not all_non_negative(v):
                return False
        elif v < 0:
            return False
    return True

def format_sci(x: float) -> str:
    """Format a positive float as scientific notation with 1 decimal and 2-digit exponent."""
    s = f"{x:.1e}"           # e.g. "3.6e-16"
    base, exp = s.split("e")
    exp_num = int(exp)
    return f"{base}e{'+' if exp_num >= 0 else '-'}{abs(exp_num):02d}"

def load_json(path):
    _, file_extension = os.path.splitext(os.path.basename(path))
    if not os.path.isfile(path) or file_extension != ".json":
        raise ValueError(f"The file {path} does not exist or is not a .json file.")
    with open(path, "r") as f:
        contents = json.load(f)
    return contents

def dump_json(dump_directory: str, file_basename: str, dict_to_dump: dict):
    """
    Dump a dictionary (possibly nested) as JSON, converting Species objects to their .name.
    Handles Species in keys, values, lists, sets, tuples, etc.
    """
    
    def convert_species(obj):
        # Convert Species objects (keys or values)
        if hasattr(obj, "__class__") and obj.__class__.__name__ == "Species":
            return obj.name

        # Handle dictionaries (convert both keys and values)
        elif isinstance(obj, dict):
            return {str(convert_species(k)): convert_species(v) for k, v in obj.items()}

        # Handle lists, tuples, sets
        elif isinstance(obj, (list, tuple, set)):
            return [convert_species(i) for i in obj]

        # Base case: leave unchanged
        else:
            return obj

    converted = convert_species(dict_to_dump)
    os.makedirs(dump_directory, exist_ok=True)
    path = os.path.join(dump_directory, f"{file_basename}.json")
    with open(path, "w") as f:
        json.dump(converted, f, indent=4)

def pickle_dump_binary(path, variable):
    with open(path, 'wb') as f:
        pickle.dump(variable, f)

def pickle_load_binary(path):
    with open(path, 'rb') as f:
        loaded_variable = pickle.load(f)
    return loaded_variable

def closest_value(my_list, target):
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



def find_max_in_nested_dict(d):
    max_val = float('-inf')
    for v in d.values():
        if isinstance(v, dict):
            # Recurse into nested dictionary
            max_val = max(max_val, find_max_in_nested_dict(v))
        elif isinstance(v, (int, float)):
            max_val = max(max_val, v)
    return max_val