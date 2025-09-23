import os

from src.createCsvFilesBase import find_new_case, create_csv_files
from src.checkValidityCsvFiles import check_validity_of_csv_files

data_folder = "data"
digits_case_numbers = 3

rule create_next_case: # python -m snakemake -s Snakefile.smk create_next_case --config case_number="none" --cores 1
    run:
        new_case_path = find_new_case(data_folder, digits_case_numbers, config['case_number'])
        os.mkdir(new_case_path)
        create_csv_files(new_case_path, files_to_create_dict)

case_folder = os.path.join(
    data_folder,
    f"case{str(config['case_number']).zfill(digits_case_numbers)}"
)

files_to_create_dict = {
    "ENZYMATIC_REACTIONS": ["start_species", "end_species", "ratio_endtostart_species", "enzyme", "k_cat", "k_M", "hill"],
    "SPONTANEOUS_REACTIONS": ["start_species", "end_species", "ratio_endtostart_species", "k"],
    "SPECIES": ["name", "diffusion_constant", "permeability_constant"],
    "ENZYMES": ["name", "quantity", "localization"],
    #"SIMULATION_CONFIG": ["radius"]
}

csv_file_names = files_to_create_dict.keys()

rule check_case_input_validity: # python -m snakemake -s Snakefile.smk check_case_input_validity --config case_number=0 --cores 1
    input:
        [os.path.join(case_folder, f"{name}.csv") for name in csv_file_names]
    output:
        [os.path.join(case_folder, f"{name}_dataframe_pickle") for name in csv_file_names]
    run:
        check_validity_of_csv_files(case_folder, csv_file_names)


"""
rule generate_reaction_network:
    input: 
        data/case{case_number}/enzymatic_reactions.csv data/case{case_number}/spontaneous_reactions.csv
    output:
        data/case{case_number}/reaction_network
    run:
        generate_reaction_network(input)

# Get the absolute path to this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Go up to the project root (one or more levels as needed)
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))

# Access other folders relative to the root
DATA_PATH = os.path.join(PROJECT_ROOT, "data", "case000")
DATA_PATH
"""
        
