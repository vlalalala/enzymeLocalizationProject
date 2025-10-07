"""
To create .csv templates, run from snakemake root directory
src/create_file_structure.py

"""
import sys
import os
import pandas as pd
from itertools import product
from pathlib import Path
from auxiliary_functions_using_standard_libraries import load_json, as_list

def create_csv_header_file(file_path, header_list):
    """Takes the files_to_create variable and creates a csv file for each key,
    with each column having headers as specified in the dictionary values.
    """
    Path(file_path).parent.mkdir(parents=True, exist_ok=True) # Creates all the parent/intermediate files if not already existent
    df = pd.DataFrame(columns=header_list)
    df.to_csv(file_path, encoding='utf-8-sig', index=False)

if __name__ == "__main__":
    # Get information about .csv files to create
    chemical_network_info_dict = load_json("src/chemical_network_info.json")
    geometry_info_dict = load_json("src/geometry_info.json")
    # Merge the two dictionaries for files with headers
    complete_info_dict = chemical_network_info_dict | geometry_info_dict
    # Define which output paths are to be created (if they do not exist yet)
    config_info = load_json("config/config.json")
    data_folder = as_list(config_info["data_folder"])
    base_name = as_list(config_info["base_name"])
    case_numbers = as_list(config_info["case_numbers"], int)
    digits_case_numbers = int(config_info["digits_case_numbers"])
    padded_case_numbers = [str(n).zfill(digits_case_numbers) for n in case_numbers]

    for df, bn, cn in product(data_folder, base_name, padded_case_numbers):
        for file_name in complete_info_dict.keys():
            file_path = f"{df}/{bn}_{cn}/{file_name}.csv"
            if os.path.isfile(file_path):
                continue
            create_csv_header_file(file_path, complete_info_dict[file_name])

    """
    # Create the csv files one by one
    for file_path in output_files:
        # If the file already exists, do nothing
        if os.path.isfile(file_path):
            print(f"{file_path} already exists")
            continue
        # Get which csv file is being created
        basename = os.path.basename(file_path).removesuffix(".csv")
        # Pass the path and headers required for that file
        create_csv_header_file(file_path, complete_info_dict[basename])
    """