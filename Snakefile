import os
from itertools import product

from src.auxFcts import as_list, create_csv_header_file
from src.checkValidityCsvFiles import check_validity_of_csv_files
#from src.defineReactionNetwork import generate_reaction_network

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
    # snakemake -s Snakefile check_chemical_network_data_validity --config case_numbers=0 base_names=violacein --cores 1 --use-conda
    input:
        [expand("{df}/{bn}_{cn}/{cf}.csv",
        df = data_folders, bn = base_names, cn = padded_case_numbers, cf = chemical_network_csv_file_names)]
    output:
        [expand("{df}/{bn}_{cn}/.{cf}_dataframe_pickle",
        df = data_folders, bn = base_names, cn = padded_case_numbers, cf = chemical_network_csv_file_names)]
    conda:
        "config/environment.yaml"
    run:
        for df, bn, cn in product(data_folders, base_names, padded_case_numbers):
            check_validity_of_csv_files(f"{df}/{bn}_{cn}/", chemical_network_csv_file_names)

#rule test_conda: # snakemake -s Snakefile test_conda --config case_number=0 --cores 1 --use-conda
#    conda: # conda env create -f config/environment.yaml
#        "config/environment.yaml"
#    run:
#        import networkx
rule test_conda:
    conda: "config\environment.yaml"
    run:
        import sys
        import subprocess
        print("Python executable:", sys.executable)
        subprocess.run(["python", "-m", "pip", "list"])
        subprocess.run(["conda", "list"])
        import networkx
        print("✅ networkx is available!")



rule create_system_network:
    # snakemake -s Snakefile.smk create_system_network --config case_number=0 --cores 1 --use-conda
    input:
        [expand("{df}/{bn}_{cn}/.{cf}_dataframe_pickle",
        df = data_folders, bn = base_names, cn = padded_case_numbers, cf = chemical_network_csv_file_names)]
    output:
        outputs = [file
                for suffix in [".NETWORK_system_pickle", "network_graph.png"]
                for file in expand(
                    "{df}/{bn}_{cn}/" + suffix,
                    df=data_folders, bn=base_names, cn=padded_case_numbers
                )
            ]

    conda:
        "config/environment.yaml"
    run:
        for df, bn, cn in product(data_folders, base_names, padded_case_numbers):
            generate_reaction_network(f"{df}/{bn}_{cn}/", chemical_network_csv_file_names)

#rule all:






"""
rule create_specific_cases:
    output:
        output_paths
    run:
        for df, bn, cn, csv in zip(data_folders, base_names, case_numbers, csv_files):
            folder = os.path.join(df, f"{bn}_{cn}")
            os.makedirs(folder, exist_ok=True)
            create_csv_files(folder, {f"{csv}.csv": files_to_create_dict[f"{csv}.csv"]})



rule create_specific_cases:
    output:
        output_paths
    run:
        for df, bn, cn, csv in zip(data_folders, base_names, case_numbers, csv_files):
            folder = os.path.join(df, f"{bn}_{cn}")
            os.makedirs(folder, exist_ok=True)
            create_csv_files(folder, {f"{csv}.csv": files_to_create_dict[f"{csv}.csv"]})
"""


"""
rule check_case_input_validity: # python -m snakemake -s Snakefile.smk check_case_input_validity --config case_number=0 --cores 1 --use-conda
    input:
        [os.path.join(case_folder, f"{name}.csv") for name in csv_file_names]
    output:
        [os.path.join(case_folder, *["intermediate_files", f"{name}_dataframe_pickle"]) for name in csv_file_names]
    conda:
        "config/defaultEnvironment.yaml"
    run:
        check_validity_of_csv_files(case_folder, csv_file_names)


rule create_system_network: # python -m snakemake -s Snakefile.smk create_system_network --config case_number=0 --cores 1 --use-conda
    input:
        [os.path.join(case_folder, *["intermediate_files", f"{name}_dataframe_pickle"]) for name in csv_file_names]
    output:
        os.path.join(case_folder, *["intermediate_files", "NETWORK_system_pickle"]), os.path.join(case_folder,"network_graph.png")
    conda:
        "config/defaultEnvironment.yaml"
    run:
        generate_reaction_network(case_folder, csv_file_names)



rule test_input_and_output:
    input:
        [expand("{data_folder}/case{caseNr}/ENZYMATIC_REACTIONS.csv", 
        data_folder=data_folder, caseNr = padded_case_numbers)]
    output:
        [expand("{data_folder}/case{caseNr}/test.txt",
        data_folder=data_folder, caseNr = padded_case_numbers)]
    run:
        new_files = expand("{data_folder}/case{caseNr}/test.txt",
        data_folder=data_folder, caseNr = padded_case_numbers)
        for file in new_files:
            with open(file, 'w') as fp:
                pass
        print("case folder", case_folder)
"""
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

case_folder = os.path.join(
    data_folder,
    f"case{str(config['case_number']).zfill(digits_case_numbers)}"
)


"""
        
