"""
Functions (to use e.g. for cleaning config inputs) that only require Python's standard library
"""
import os
import ast
import json
import pickle
import glob
import re


def load_json(path):
    _, file_extension = os.path.splitext(os.path.basename(path))
    if not os.path.isfile(path) or file_extension != ".json":
        raise ValueError(f"The file {path} does not exist or is not a .json file.")
    with open(path, "r") as f:
        contents = json.load(f)
    return contents

def dump_json(dump_directory: str, file_basename: str, dict_to_dump: dict):
    with open(os.path.join(dump_directory, f"{file_basename}.json"), "w") as f:
        json.dump(dict_to_dump, f, indent=4)

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