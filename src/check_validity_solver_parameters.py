import sys
import os
from auxiliary_functions_using_standard_library import load_json, pickle_dump_binary, closest_value

def check_validity_solver_parameters(case_directory):
    """Checks that in the solver info, everything has the correct type
    """
    solver_info_nested_dict = load_json(os.path.join(case_directory, "solver_info.json"))

    for section, section_dict in solver_info_nested_dict.items():
        for key, value in section_dict.items():
            key_lower = key.lower()

            # --- Rule 1: Keys containing "num" or "every" must be int ---
            if any(sub in key_lower for sub in ["num", "every"]):
                # Convert numeric strings or floats (like 1e2) to int if possible
                if isinstance(value, (float, int)):
                    if isinstance(value, float) and value.is_integer():
                        value = int(value)
                elif isinstance(value, str):
                    # Skip if it's clearly a boolean-looking string
                    if value.lower() in ["true", "false"]:
                        raise TypeError(
                            f"Invalid type for '{key}' in '{section}': expected int, got str ('{value}')"
                        )
                    try:
                        float_val = float(value)
                        if float_val.is_integer():
                            value = int(float_val)
                        else:
                            raise TypeError(
                                f"Invalid numeric string for '{key}' in '{section}': expected int-like value, got '{value}'"
                            )
                    except ValueError:
                        raise TypeError(
                            f"Invalid type for '{key}' in '{section}': expected int, got str ('{value}')"
                        )
                else:
                    raise TypeError(
                        f"Invalid type for '{key}' in '{section}': expected int, got {type(value).__name__}"
                    )

                section_dict[key] = value  # store possibly converted int

            # --- Rule 2: Keys with action verbs must be bool ---
            elif any(sub in key_lower for sub in ["plot", "save", "print", "override", "create", "delete"]):
                if isinstance(value, str):
                    if value.lower() == "true":
                        value = True
                    elif value.lower() == "false":
                        value = False
                    else:
                        raise ValueError(
                            f"Invalid string for boolean '{key}' in '{section}': expected 'True' or 'False', got '{value}'"
                        )
                elif not isinstance(value, bool):
                    raise TypeError(
                        f"Invalid type for '{key}' in '{section}': expected bool, got {type(value).__name__}"
                    )

                section_dict[key] = value

            # --- Rule 3: "alpha", "gamma", "tol" must be positive numbers ---
            elif any(sub in key_lower for sub in ["alpha", "gamma", "tol"]):
                if not isinstance(value, (int, float)):
                    raise TypeError(
                        f"Invalid type for '{key}' in '{section}': expected positive int or float, got {type(value).__name__}"
                    )
                if value <= 0:
                    raise ValueError(
                        f"Invalid value for '{key}' in '{section}': expected positive number, got {value}"
                    )

                # --- Rule 4: Relative "inc"/"dec" checks ---
                if "inc" in key_lower and value <= 1:
                    raise ValueError(
                        f"Invalid value for '{key}' in '{section}': expected value higher than 1, got {value}"
                    )
                elif "dec" in key_lower and value >= 1:
                    raise ValueError(
                        f"Invalid value for '{key}' in '{section}': expected value lower than 1, got {value}"
                    )

        # --- Rule 5: Enforce _min/_max consistency within each section ---
        for key in list(section_dict.keys()):
            if key.endswith("_min"):
                base = key[:-4]  # remove "_min"
                max_key = base + "_max"

                if max_key in section_dict:
                    min_val = section_dict[key]
                    max_val = section_dict[max_key]

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