import sys
import os
from typing import Dict, Any, Union
import numpy as np
from auxiliary_functions_using_standard_library import load_json, closest_value
from auxiliary_functions import dump_json

def check_validity_system_geometry_info(case_directory):
    """
    """
    
    system_geometry_nested_dict : Dict[Union[str, int], Dict[str, Any]] = load_json(os.path.join(case_directory, "system_geometry.json"))
    solver_input_nested_dict : Dict[Union[str, int], Dict[str, Any]] = load_json(os.path.join(case_directory, "solver_input.json"))
    
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
            if section=="geometry_config" and key=="num_dimensions" and value!=3:
                raise NotImplementedError("Only 3-dimensional systems are implemented to date.")

    #if (system_geometry_nested_dict["GEOMETRY_CONFIG"]["inner_mesh_radius_relative_to_smallest_membrane"]
    #    >= system_geometry_nested_dict["GEOMETRY_CONFIG"]["inner_mesh_radius_relative_to_smallest_membrane"]):

    external_radius = system_geometry_nested_dict["geometry_config"]["outer_membrane_radius"]
    membrane_input_radii = [relative * external_radius for relative in system_geometry_nested_dict["geometry_config"]["internal_membrane_relative_radii"] + [1]]
    system_geometry_nested_dict["geometry_config"]["membrane_input_radii"] = membrane_input_radii
    mesh_points = np.linspace(0, external_radius, num = solver_input_nested_dict["geometry_parameters"]["num_mesh_points"])
    system_geometry_nested_dict["geometry_config"]["mesh_points"] = mesh_points
    membrane_radii = [closest_value(mesh_points, membrane_input_radius) for membrane_input_radius in membrane_input_radii]
    system_geometry_nested_dict["geometry_config"]["membrane_radii"] = membrane_radii
    num_regions = len(membrane_radii)
    system_geometry_nested_dict["geometry_config"]["num_regions"]  = num_regions
    boundary_radii = [0] + membrane_radii
    system_geometry_nested_dict["geometry_config"]["boundary_radii"]  = boundary_radii
    mesh_points_in_regions = {
        region_idx : [mesh_point for mesh_point in mesh_points if boundary_radii[region_idx]<=mesh_point<=boundary_radii[region_idx+1]]
        for region_idx in range(num_regions)
    }
    system_geometry_nested_dict["geometry_config"]["mesh_points_in_regions"]  = mesh_points_in_regions
    num_mesh_points_in_regions = {
        region_idx: len(mesh_points_in_region)
        for region_idx, mesh_points_in_region in mesh_points_in_regions.items()
    }
    system_geometry_nested_dict["geometry_config"]["num_mesh_points_in_regions"]  = num_mesh_points_in_regions

    dump_json(
        case_directory,
        ".expanded_system_geometry",
        system_geometry_nested_dict
    )

if __name__ == "__main__":
    folder_to_check_validity = sys.argv[1]
    check_validity_system_geometry_info(folder_to_check_validity)

