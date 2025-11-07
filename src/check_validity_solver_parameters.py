import sys
import os
import numpy as np
from auxiliary_functions_using_standard_library import load_json, pickle_dump_binary, closest_value

def check_validity_solver_parameters(case_directory):
    """Checks that in the solver info, everything has the correct type
    """
    solver_info_nested_dict = load_json(os.path.join(case_directory, "solver_info.json"))
    
    for section, section_dict in solver_info_nested_dict.items():
        for key, value in section_dict.items():
            # Rule 1: Keys containing "num" or "every" must be int.
            if any(integer_substring in key.lower() for integer_substring in ["num", "every"]):
                # If a float is given e.g. in scientific notation, then the value is converted to int
                if isinstance(value, float) and value.is_integer():
                    value = int(value)
                    section_dict[key] = value
                # After conversion attempts, check type
                if not isinstance(value, int):
                    raise TypeError(
                        f"Invalid type for '{key}' in '{section}': expected int, got {type(value).__name__}"
                    )
            # Rule 2: (only if num or every are not present) check if the key has a verb and saves True or False depending on "True" or "False" written
            elif any(integer_substring in key.lower() for integer_substring in ["save", "print", "override", "create", "delete"]):
                if isinstance(value, str):
                    if value.lower() == "true":
                        value = True
                        section_dict[key] = value
                    elif value.lower() == "false":
                        value = False
                        section_dict[key] = value
                    else:
                        raise ValueError(
                            f"Invalid string for boolean '{key}' in '{section}': expected 'True' or 'False', got '{value}'"
                        )

                if not isinstance(value, bool):
                    raise TypeError(
                        f"Invalid type for '{key}' in '{section}': expected bool, got {type(value).__name__}"
                    )                

            # Rule 3: Keys containing "alpha", "gamma", "tol" must be int or float and equal or larger than zero
            elif any(float_substring in key.lower() for float_substring in ["alpha", "gamma", "tol"]):
                if not isinstance(value, (int, float)) or not value>0:
                    raise TypeError(
                        f"Invalid type for '{key}' in '{section}': expected positive int or float, got {type(value).__name__}"
                    )
                # Rule 4: if it has inc, only values above 1 are permitted
                if any(relative_substring in key.lower() for relative_substring in ["inc"]):
                    if value<=1:
                        raise ValueError(
                            f"Invalid value for '{key}' in '{section}': expected value higher than 1, got {value}"
                        )
                elif any(relative_substring in key.lower() for relative_substring in ["dec"]):
                    if value>=1:
                        raise ValueError(
                            f"Invalid value for '{key}' in '{section}': expected value lower than 1, got {value}"
                        )
        # --- Enforce min/max consistency ---
        for key in list(section_dict.keys()):
            if key.endswith("_min"):
                base_key = key[:-4]  # remove "_min"
                max_key = base_key + "_max"
                if max_key in section_dict:
                    min_val = section_dict[key]
                    max_val = section_dict[max_key]

                    # Try to convert both to float for comparison (handles int, float)
                    try:
                        min_val_f = float(min_val)
                        max_val_f = float(max_val)
                    except (TypeError, ValueError):
                        raise TypeError(
                            f"Values for '{key}' and '{max_key}' in '{section}' must be numeric for comparison."
                        )

                    if min_val_f > max_val_f:
                        raise ValueError(
                            f"In '{section}', '{key}' ({min_val_f}) cannot be greater than '{max_key}' ({max_val_f})."
                        )
   
    pickle_dump_binary(
        os.path.join(case_directory, ".solver_info_pickle"),
        solver_info_nested_dict
    )

if __name__ == "__main__":
    folder_to_check_validity = sys.argv[1]
    check_validity_solver_parameters(folder_to_check_validity)