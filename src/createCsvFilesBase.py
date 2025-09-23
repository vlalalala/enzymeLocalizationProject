import os
import pandas as pd

def find_new_case(case_parent_directory: str, digits_case_numbers:int, case_number) -> str:
    """Returns the absolute path of the new case folder to be created.
    If a case number is an int, it tries to create that one (fails if it doesn't fit or has 
    already been used). If case_number = "none", it takes the next one 1 higher up than the existing max.
    """
    # Step 1: List all case integers in the data directory that match 'caseNNN'
    subfolders = [subfolder for subfolder in os.listdir(case_parent_directory)
                    if "case" in subfolder]
    existing_cases = [int(case_string.replace("case", ""))
                        for case_string in subfolders
                        if "case" in case_string]
    # Step 2: Get next case number
    if case_number == "none": # If "?" passed make the new case number the next available (from the end)
        if existing_cases:
            new_case_number = max(existing_cases) + 1
        else:
            new_case_number = 0
    elif not (0 <= int(case_number) < 10**digits_case_numbers): # must fit within the digits used
        raise ValueError("The case_number given cannot be converted to a valid number.")
    elif int(case_number) in existing_cases: # must not be already used
        raise ValueError("The case_number given already exists.")
    else:
        new_case_number = int(case_number)
    
    # Step 3: Format the new folder name
    new_case = f"case{str(new_case_number).zfill(digits_case_numbers)}"
    new_case_dir = os.path.join(case_parent_directory, new_case)
    
    return new_case_dir

def create_csv_files(case_directory, files_to_create_dict) -> None:
    """Takes the files_to_create variable and creates a csv file for each key,
    with each column having headers as specified in the dictionary values.
    """
    for csv_file_name, header_list in files_to_create_dict.items():
        csv_file_path = os.path.join(case_directory, f'{csv_file_name}.csv')
        df = pd.DataFrame(columns=header_list)
        df.to_csv(csv_file_path, encoding='utf-8-sig', index=False)
        print(f"{csv_file_name}.csv", "created")

#%%
