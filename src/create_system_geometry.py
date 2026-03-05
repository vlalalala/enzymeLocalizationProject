import sys
import os
from typing import Dict, Any, Union
import numpy as np
from auxiliary_functions_using_standard_library import closest_value
from auxiliary_functions import dump_json, read_yaml_file

def check_validity_system_geometry_info(nested_dict_to_check):
    """Raises errors if there 
    """
    for section, section_dict in nested_dict_to_check.items():
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

def construct_baseline_mesh_points(
    case_directory,
    system_geometry_nested_dict,
    discretization_params_dict
):  
    """
    Expands and saves file with the system geometry expanded to accomodate
    for the discretization (the membrane positions are set to the closest 
    mesh position given the (baseline) discretization chosen.)
    """
    external_radius = system_geometry_nested_dict["geometry_config"]["outer_membrane_radius"]
    membrane_input_radii = [relative * external_radius for relative in system_geometry_nested_dict["geometry_config"]["internal_membrane_relative_radii"] + [1]]
    system_geometry_nested_dict["geometry_config"]["membrane_input_radii"] = membrane_input_radii
    mesh_points = np.linspace(0, external_radius, num = discretization_params_dict["discretization_parameters"]["baseline_num_mesh_points"])
    system_geometry_nested_dict["geometry_config"]["baseline_mesh_points"] = mesh_points
    membrane_radii = [closest_value(mesh_points, membrane_input_radius) for membrane_input_radius in membrane_input_radii]
    system_geometry_nested_dict["geometry_config"]["membrane_radii"] = membrane_radii
    num_regions = len(membrane_radii)
    system_geometry_nested_dict["geometry_config"]["num_regions"]  = num_regions
    boundary_radii = [0] + membrane_radii
    system_geometry_nested_dict["geometry_config"]["boundary_radii"]  = boundary_radii
    volume_regions = {region: 4*np.pi/3 * (boundary_radii[region+1]**3 - boundary_radii[region]**3)
                      for region in range(num_regions)}
    system_geometry_nested_dict["geometry_config"]["volume_regions"]  = volume_regions

    dump_json(
        case_directory,
        ".system_geometry",
        system_geometry_nested_dict
    )

if __name__ == "__main__":
    case_directory = sys.argv[1]
    # Import files with user input
    system_geometry_nested_dict = read_yaml_file(os.path.join(case_directory, "parameters_geometry.yaml"))
    discretization_params_dict = read_yaml_file(os.path.join(case_directory, "parameters_discretization.yaml"))
    # Check that user input is valid
    check_validity_system_geometry_info(system_geometry_nested_dict)
    check_validity_system_geometry_info(discretization_params_dict)
    # Save baseline system geometry
    construct_baseline_mesh_points(
        case_directory,
        system_geometry_nested_dict,
        discretization_params_dict
    )
    

