import os
from itertools import product

# ast is a module in the python standard library
# change of plan: outside the rules, only use and import stuff that builds on the standard library.

from src.auxFcts import as_list, create_csv_header_file
from src.checkValidityCsvFiles import check_validity_of_csv_files
from src.defineReactionNetwork import generate_reaction_network

configfile: os.path.join("config", "config.yaml")

data_folders = as_list(config["data_folders"], str)
base_names = as_list(config["base_names"], str)
case_numbers = as_list(config["case_numbers"], int)
digits_case_numbers = int(config["digits_case_numbers"])

padded_case_numbers = [str(n).zfill(digits_case_numbers) for n in case_numbers]

chemical_network_dict = {
    "ENZYMATIC_REACTIONS": ["start_species", "end_species", "ratio_endtostart_species", "enzyme", "k_cat", "k_M", "hill"],
    "SPONTANEOUS_REACTIONS": ["start_species", "end_species", "ratio_endtostart_species", "k"],
    "ENZYMES": ["name", "quantity", "localization"],
    "SPECIES": ["name", "diffusion_constant", "permeability_constant"],
}

system_parameters_dict = {
    "SIMULATION_CONFIG": ["radius"],
}

chemical_network_csv_file_names = chemical_network_dict.keys()
system_parameters_csv_file_names = system_parameters_dict.keys()

files_to_create_dict = chemical_network_csv_file_names | system_parameters_csv_file_names # merge the two dictionaries
# changes to the mutable values in files_to_create_dict will reflect in the original dictionaries, and vice versa

rule create_cases: # snakemake -s Snakefile.smk create_cases --config case_numbers = --cores 1
    input:
        data_folders
    output:
        [expand("{df}/{bn}_{cn}/{cf}.csv",
        df = data_folders, bn = base_names, cn = padded_case_numbers, cf = files_to_create_dict)]
    run:
        for df, bn, cn, cf in product(data_folders, base_names, padded_case_numbers, files_to_create_dict):
            folder = os.path.join(df, f"{bn}_{cn}")
            os.makedirs(folder, exist_ok=True)
            create_csv_header_file(os.path.join(folder, f"{cf}.csv"), files_to_create_dict[cf])

rule check_chemical_network_data_validity:
    # snakemake -s Snakefile check_chemical_network_data_validity --config case_numbers=0 base_names=violacein --cores 1
    input:
        [expand("{df}/{bn}_{cn}/{cf}.csv",
        df = data_folders, bn = base_names, cn = padded_case_numbers, cf = chemical_network_csv_file_names)]
    output:
        [expand("{df}/{bn}_{cn}/.{cf}_dataframe_pickle",
        df = data_folders, bn = base_names, cn = padded_case_numbers, cf = chemical_network_csv_file_names)]
    run:
        for df, bn, cn in product(data_folders, base_names, padded_case_numbers):
            check_validity_of_csv_files(f"{df}/{bn}_{cn}/", chemical_network_csv_file_names)

rule create_system_network:
    # snakemake -s Snakefile.smk create_system_network --config base_names = violacein case_numbers=0 --cores 1
    input:
        [expand("{df}/{bn}_{cn}/.{cf}_dataframe_pickle",
        df = data_folders, bn = base_names, cn = padded_case_numbers, cf = chemical_network_csv_file_names)]
    output:
        outputs = [file
                for suffix in [".NETWORK_system_pickle", "network_graph.png"]
                for file in expand(
                    "{df}/{bn}_{cn}/" + suffix,
                    df=data_folders, bn=base_names, cn=padded_case_numbers)]
    run:
        for df, bn, cn in product(data_folders, base_names, padded_case_numbers):
            generate_reaction_network(f"{df}/{bn}_{cn}/", chemical_network_csv_file_names)


rule test_conda:
    # conda config --set channel_priority strict
    conda:
        "config/snakerunner_environment.yaml"
    run:
        import sys
        import subprocess
        print("Python executable:", sys.executable)
        subprocess.run(["python", "-m", "pip", "list"])
        subprocess.run(["conda", "list"])
        import networkx
        print("✅ networkx is available!")


"""
# Get the absolute path to this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Go up to the project root (one or more levels as needed)
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))

# Access other folders relative to the root
DATA_PATH = os.path.join(PROJECT_ROOT, "data", "case000")
DATA_PATH

case_folder = os.path.join(
    data_folder,
    f"case{str(config['case_number']).zfill(digits_case_numbers)}"
)


"""
        
