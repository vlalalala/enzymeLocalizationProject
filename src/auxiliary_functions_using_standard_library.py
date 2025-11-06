"""
Functions (to use e.g. for cleaning config inputs) that only require Python's standard library
"""
import os
import ast
import json
import pickle
import glob
import re

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
    
def find_sorted_file_names(folder, pattern_to_find):
    """ pattern_to_find must look somewhat like this: ".iteration_nr_*_concentration.png"
    To have them sorted correctly, the file names should have the correct zero-padding.
    """
    pattern = os.path.join(folder, pattern_to_find)
    files = glob.glob(pattern)
    files.sort()
    return files

def find_max_in_nested_dict(d):
    max_val = float('-inf')
    for v in d.values():
        if isinstance(v, dict):
            # Recurse into nested dictionary
            max_val = max(max_val, find_max_in_nested_dict(v))
        elif isinstance(v, (int, float)):
            max_val = max(max_val, v)
    return max_val