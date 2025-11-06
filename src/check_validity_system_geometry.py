import sys
import os
import numpy as np
from auxiliary_functions_using_standard_library import load_json, pickle_dump_binary, closest_value

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
                # Rule 4: if it has radius and relative in it, values above 1 are also not permitted
                if any(relative_substring in key.lower() for relative_substring in ["relative"]):
                    if value>=1:
                        raise ValueError(
                            f"Invalid value for '{key}' in '{section}': expected value under 1, got {value}"
                        )                        
            # Rule 4: check type of values in lists
            elif any(float_list_substring in key.lower() for float_list_substring in ["radii"]):
                if not isinstance(value, (list, np.ndarray)):
                    raise TypeError(
                        f"Invalid type for '{key}' in '{section}': expected list or np.array, got {type(value).__name__}"
                    )
                if any(relative_substring in key.lower() for relative_substring in ["relative"]):
                    if any(element>1 or element < 0 for element in value):
                        raise ValueError(
                        f"Invalid value for '{key}' in '{section}': expected all values between 0 and 1 (inclusive)"
                        )
                if any(not isinstance(element, (int, float)) for element in value):
                    raise TypeError(
                        f"Invalid type for '{key}' in '{section}': expected all elements with type int or float"
                    )
            # Rule 5: if anything is not yet programmed, raise NotImplementedError
            if section=="GEOMETRY_CONFIG" and key=="num_dimensions" and value!=3:
                raise NotImplementedError("Only 3-dimensional systems are implemented to date.")

    #if (system_geometry_nested_dict["GEOMETRY_CONFIG"]["inner_mesh_radius_relative_to_smallest_membrane"]
    #    >= system_geometry_nested_dict["GEOMETRY_CONFIG"]["inner_mesh_radius_relative_to_smallest_membrane"]):

    R = system_geometry_nested_dict["GEOMETRY_CONFIG"]["outer_membrane_radius"]
    MEMBRANE_INPUT_RADII = [relative * R for relative in system_geometry_nested_dict["GEOMETRY_CONFIG"]["internal_membrane_relative_radii"] + [1]]
    system_geometry_nested_dict["GEOMETRY_CONFIG"]["MEMBRANE_INPUT_RADII"] = MEMBRANE_INPUT_RADII
    MESH_POINTS = np.linspace(0, R, num = system_geometry_nested_dict["GEOMETRY_CONFIG"]["num_mesh_points"])
    system_geometry_nested_dict["GEOMETRY_CONFIG"]["MESH_POINTS"] = MESH_POINTS
    MEMBRANE_RADII = [closest_value(MESH_POINTS, membrane_input_radius) for membrane_input_radius in MEMBRANE_INPUT_RADII]
    system_geometry_nested_dict["GEOMETRY_CONFIG"]["MEMBRANE_RADII"] = MEMBRANE_RADII
    NUM_REGIONS = len(MEMBRANE_RADII)
    system_geometry_nested_dict["GEOMETRY_CONFIG"]["NUM_REGIONS"]  = NUM_REGIONS
    BOUNDARY_RADII = [0] + MEMBRANE_RADII
    system_geometry_nested_dict["GEOMETRY_CONFIG"]["BOUNDARY_RADII"]  = BOUNDARY_RADII
    MESH_POINTS_IN_REGIONS = {
        region_idx : [mesh_point for mesh_point in MESH_POINTS if BOUNDARY_RADII[region_idx]<=mesh_point<=BOUNDARY_RADII[region_idx+1]]
        for region_idx in range(NUM_REGIONS)
    }
    system_geometry_nested_dict["GEOMETRY_CONFIG"]["MESH_POINTS_IN_REGIONS"]  = MESH_POINTS_IN_REGIONS
    NUM_MESH_POINTS_IN_REGIONS = {
        region_idx: len(mesh_points_in_region)
        for region_idx, mesh_points_in_region in MESH_POINTS_IN_REGIONS.items()
    }
    system_geometry_nested_dict["GEOMETRY_CONFIG"]["NUM_MESH_POINTS_IN_REGIONS"]  = NUM_MESH_POINTS_IN_REGIONS

    pickle_dump_binary(
        os.path.join(case_directory, ".SYSTEM_GEOMETRY_pickle"),
        system_geometry_nested_dict
    )
    


if __name__ == "__main__":
    folder_to_check_validity = sys.argv[1]
    check_validity_system_geometry_info(folder_to_check_validity)