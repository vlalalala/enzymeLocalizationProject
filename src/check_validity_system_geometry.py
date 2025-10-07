import sys
import os
from auxiliary_functions_using_standard_library import load_json

def check_validity_system_geometry_info(case_directory):
    
    system_geometry_nested_dict = load_json(os.path.join(case_directory, "SYSTEM_GEOMETRY.json"))
    
    for section, section_dict in system_geometry_nested_dict.items():
        for key, value in section_dict.items():
            # Rule 1: Keys containing "num" must be int
            if any(integer_substring in key.lower() for integer_substring in ["num"]):
                if not isinstance(value, int):
                    raise TypeError(
                        f"Invalid type for '{key}' in '{section}': expected int, got {type(value).__name__}"
                    )
            # Rule 2: Keys containing "radius" must be int or float
            elif any(float_substring in key.lower() for float_substring in ["radius"]):
                if not isinstance(value, (int, float)):
                    raise TypeError(
                        f"Invalid type for '{key}' in '{section}': expected int or float, got {type(value).__name__}"
                    )
                # Rule 3: Keys containing "radius" must be strictly larger than 0
                if any(positive_value_key in key.lower() for positive_value_key in ["radius"]):
                    if not value>0:
                        raise ValueError(
                            f"Invalid type for '{key}' in '{section}': expected positive value, got {value}"
                        )
            # Rule 4: if anything is not yet programmed, raise NotImplementedError
            if section=="GEOMETRY_CONFIG" and key=="num_dimensions" and value!=3:
                raise NotImplementedError("Only 3-dimensional systems are implemented to date.")

if __name__ == "__main__":
    folder_to_check_validity = sys.argv[1]
    check_validity_system_geometry_info(folder_to_check_validity)