import sys
import shutil
import os
from pathlib import Path
import pandas as pd
from auxiliary_functions_using_standard_library import load_json

def copy_and_rename_files(target_dir, source_dir, old_substr, new_substr):
    """
    Copies files from source_dir to target_dir, changing a substring in the filenames.

    Parameters:
    - target_dir: str, path to the folder where copies will go
    - source_dir: str, path to folder where original files are
    - old_substr: str, substring to find in filenames
    - new_substr: str, substring to replace it with
    """
    os.makedirs(target_dir, exist_ok=True)  # create target folder if it doesn't exist
    for filename in os.listdir(source_dir):
        if old_substr in filename and filename.endswith('.yaml'):
            new_filename = filename.replace(old_substr, new_substr)
            src_path = os.path.join(source_dir, filename)
            dst_path = os.path.join(target_dir, new_filename)
            shutil.copy2(src_path, dst_path)  # preserves metadata
            print(f"Copied: {filename} → {new_filename}")

def create_csv_header_file(file_path, header_list):
    """Creates a .csv file in the location file_path (which does not require the .csv ending)
    with each element in header_list in a cell within the first row of the .csv file.
    """
    Path(file_path).parent.mkdir(parents=True, exist_ok=True) # Creates all the parent/intermediate files if not already existent
    dataframe = pd.DataFrame(columns=header_list)
    dataframe.to_csv(file_path, encoding='utf-8-sig', index=False)

if __name__ == "__main__":
    new_folder_with_parameter_ranges = sys.argv[1]
    src_path = Path(__file__).parent
    copy_and_rename_files(
        new_folder_with_parameter_ranges,
        src_path,
        "_template",
        "parameters"
    )
    for filename in os.listdir(src_path):
        if not ("_template" in filename and filename.endswith('.json')):
            continue
        reaction_network_header_dict = load_json(src_path / filename)
        for file_name in reaction_network_header_dict.keys():
            file_path = os.path.join(new_folder_with_parameter_ranges, file_name)
            create_csv_header_file(f"{file_path}.csv", reaction_network_header_dict[file_name])
            print(f"Created: csv file for {os.path.basename(file_path)}")
