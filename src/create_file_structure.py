"""
To create .csv templates, run from snakemake root directory
src/create_file_structure.py

"""
import sys
import os
from itertools import product
from pathlib import Path
import pandas as pd
from auxiliary_functions_using_standard_library import load_json, as_list, dump_json

def create_csv_header_file(file_path, header_list):
    """Creates a .csv file in the location file_path (which does not require the .csv ending)
    with each element in header_list in a cell within the first row of the .csv file.
    """
    Path(file_path).parent.mkdir(parents=True, exist_ok=True) # Creates all the parent/intermediate files if not already existent
    dataframe = pd.DataFrame(columns=header_list)
    dataframe.to_csv(file_path, encoding='utf-8-sig', index=False)

def get_nested_dict_with_None_values(path_to_json_file_describing_dict_structure):
    """Takes a dictionary where the values are lists and creates a .json file
    that stores a "nested" dictionary.
    """
    dict_from_original_json_file = load_json(path_to_json_file_describing_dict_structure)
    # Create the nested dictionary
    nested = {
        parent_key: {child_key: None for child_key in child_keys}
        for parent_key, child_keys in dict_from_original_json_file.items()
    }
    return nested

def create_system_geometry_json_file(path_to_json_file_describing_dict_structure, membrane_type, dump_directory, dump_basename):
    """Creates the json file with the geometry and in case the membrane type is enzymatic, adds the necessary parameters.
    """
    nested = get_nested_dict_with_None_values(path_to_json_file_describing_dict_structure)
    if membrane_type == "enzymatic":
        nested.update({"MEMBRANE_PROPERTIES": {"pore_density": None}})
    # Save the new nested dictionary to a new JSON file
    dump_json(dump_directory, dump_basename, nested)

def create_solver_parameters_json_file(path_to_json_file_describing_dict_structure, dump_directory, dump_basename):
    """Creates the json file with the solver parameters.
    """
    nested = get_nested_dict_with_None_values(path_to_json_file_describing_dict_structure)
    dump_json(dump_directory, dump_basename, nested)

if __name__ == "__main__":
    # Take as an input on the command line which membrane model we are using
    membrane_type = sys.argv[1]
    membrane_type_header_map = {
        "permeability": ["permeability_constant"],
        "enzymatic" : ["k_on", "k_off"]
    }
    if membrane_type not in membrane_type_header_map:
        raise ValueError(
            f"The membrane types available are {membrane_type_header_map.keys()}, and {membrane_type} is not one of them.")

    # Define which output paths are to be created (if they do not exist yet)
    config_info = load_json("config/config.json")
    data_folder = as_list(config_info["data_folder"])
    base_name = as_list(config_info["base_name"])
    case_numbers = as_list(config_info["case_numbers"], int)
    digits_case_numbers = int(config_info["digits_case_numbers"])
    padded_case_numbers = [str(n).zfill(digits_case_numbers) for n in case_numbers]
    # Get information about .csv files to create
    reaction_network_info_dict = load_json("src/reaction_network_info.json")

    # Modify headers for information of species about interaction with semipermeable membranes 
    reaction_network_info_dict["SPECIES"].remove("*membrane_parameters")
    reaction_network_info_dict["SPECIES"] += membrane_type_header_map[membrane_type]

    # Create all necessary files in which to fill in specific system data
    for df, bn, cn in product(data_folder, base_name, padded_case_numbers):
        for file_name in reaction_network_info_dict.keys():
            file_path = f"{df}/{bn}_{cn}/{file_name}.csv"
            if os.path.isfile(file_path):
                continue
            create_csv_header_file(file_path, reaction_network_info_dict[file_name])
        solver_iteration_data_folder = f"{df}/{bn}_{cn}/solver_iteration_data"
        os.makedirs(solver_iteration_data_folder, exist_ok=True)
        with open(os.path.join(solver_iteration_data_folder, ".iteration_data"), 'w') as f:
            pass  # Ensure no content is written
        create_system_geometry_json_file(
            "src/geometry_info.json", membrane_type, f"{df}/{bn}_{cn}", "SYSTEM_GEOMETRY")
        create_solver_parameters_json_file("src/solver_info.json", f"{df}/{bn}_{cn}", "solver_info")

