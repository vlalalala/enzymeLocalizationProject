"""
Functions (to use e.g. for cleaning config inputs) that only require Python's standard library
"""
import os
import ast
import json

def load_json(path):
    print(path, "exists", os.path.isfile(path))
    with open(path, "r") as f:
        contents = json.load(f)
    return contents

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