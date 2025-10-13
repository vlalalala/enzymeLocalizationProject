"""
To create .csv templates, run from snakemake root directory
src/create_file_structure.py

"""
import os
from itertools import product
from pathlib import Path
import pandas as pd
from auxiliary_functions_using_standard_library import load_json, as_list, dump_json

def create_csv_header_file(file_path, header_list):
    """Takes the files_to_create variable and creates a .csv file for each key,
    with each column having headers as specified in the dictionary values.
    """
    Path(file_path).parent.mkdir(parents=True, exist_ok=True) # Creates all the parent/intermediate files if not already existent
    df = pd.DataFrame(columns=header_list)
    df.to_csv(file_path, encoding='utf-8-sig', index=False)

def create_data_json_file(geometry_info_json_path, dump_directory, dump_basename):
    """
    Takes a dictionary where the values are lists and creates a .json file
    that stores a "nested" dictionary.
    """
    dict_from_original_json_file = load_json(geometry_info_json_path)
    # Create the nested dictionary
    nested = {
        parent_key: {child_key: None for child_key in child_keys}
        for parent_key, child_keys in dict_from_original_json_file.items()
    }
    # Save the new nested dictionary to a new JSON file
    dump_json(dump_directory, dump_basename, nested)

if __name__ == "__main__":
    # Define which output paths are to be created (if they do not exist yet)
    config_info = load_json("config/config.json")
    data_folder = as_list(config_info["data_folder"])
    base_name = as_list(config_info["base_name"])
    case_numbers = as_list(config_info["case_numbers"], int)
    digits_case_numbers = int(config_info["digits_case_numbers"])
    padded_case_numbers = [str(n).zfill(digits_case_numbers) for n in case_numbers]

    # Get information about .csv files to create
    reaction_network_info_dict = load_json("src/reaction_network_info.json")
    # Create all necessary files in which to fill in specific system data
    for df, bn, cn in product(data_folder, base_name, padded_case_numbers):
        for file_name in reaction_network_info_dict.keys():
            file_path = f"{df}/{bn}_{cn}/{file_name}.csv"
            if os.path.isfile(file_path):
                continue
            create_csv_header_file(file_path, reaction_network_info_dict[file_name])
        create_data_json_file(
            "src/geometry_info.json", f"{df}/{bn}_{cn}", "SYSTEM_GEOMETRY")

