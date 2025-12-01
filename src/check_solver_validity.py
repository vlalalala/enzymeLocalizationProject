import sys
import os
from auxiliary_functions_using_standard_library import load_json

def check_validity_file(case_directory, file_name):
    """For a particular file (file_name, with ending) in the directory case_directory,
    checks whether the input data is valid.
    """
    solver_info_nested_dict = load_json(os.path.join(case_directory, file_name))

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
                    try:
                        float_val = float(value)
                        if float_val.is_integer():
                            value = int(float_val)
                        else:
                            raise TypeError(
                                f"Invalid numeric string for '{key}' in '{section}': expected int-like value, got '{value}'"
                            )
                    except ValueError as exc:
                        raise TypeError(
                            f"Invalid type for '{key}' in '{section}': expected int, got str ('{value}')"
                        ) from exc

                section_dict[key] = value  # store possibly converted int

            # --- Rule 2: Keys with action verbs must be bool ---
            elif any(sub in key_lower for sub in ["plot", "save", "print", "override", "create", "delete"]):
                if not isinstance(value, bool):
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
                    except (TypeError, ValueError) as exc:
                        raise TypeError(
                            f"Values for '{key}' and '{max_key}' in '{section}' must be numeric for comparison."
                        ) from exc

                    if min_val_f > max_val_f:
                        raise ValueError(
                            f"In '{section}', '{key}' ({min_val_f}) cannot be greater than '{max_key}' ({max_val_f})."
                        )
    return solver_info_nested_dict

if __name__ == "__main__":
    folder_to_check_validity = sys.argv[1]
    for file_name in ["solver_input", "solver_params"]:
        solver_info_nested_dict = check_validity_file(folder_to_check_validity, f"{file_name}.json")