import sys
import os
import pandas as pd
from pathlib import Path
from standardLibraryAuxFcts import as_list, load_json
#from src.auxFcts import load_yaml_as_dict

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
    # Get the full path of the files to create
    output_files = sys.argv[1:]
    # Create the csv files one by one
    for file_path in output_files:
        # Get which csv file is being created
        basename = os.path.basename(file_path).removesuffix(".csv")
        # Pass the path and headers required for that file
        create_csv_header_file(file_path, complete_info_dict[basename])