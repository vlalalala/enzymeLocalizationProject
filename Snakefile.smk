### Import libraries from Python standard libraries ###

import os
from itertools import product
from src.standardLibraryAuxFcts import as_list, load_json

### Setup import of inputs ###

configfile: "config/config.yaml"

data_folder = config["data_folder"]
base_name = config["base_name"]

case_numbers = as_list(config["case_numbers"], int)
digits_case_numbers = int(config["digits_case_numbers"])
padded_case_numbers = [str(n).zfill(digits_case_numbers) for n in case_numbers]

### Get template for information from input files ###

chemical_network_info_dict = load_json("src/chemical_network_info.json")
geometry_info_dict = load_json("src/geometry_info.json")

complete_info_dict = chemical_network_info_dict | geometry_info_dict # merge the two dictionaries

"""
Important: Snakemake reruns the entire rule if any output file is missing (aka will delete output files and remake them if any one of them is missing!)
"""


rule create_cases:
    # snakemake -s Snakefile.smk create_cases --cores 1
    output:
        [expand("{df}/{bn}_{cn}/{cf}.csv",
        df = data_folder, bn = base_name, cn = padded_case_numbers, cf = complete_info_dict.keys())]
    conda:
        "config/environment.yaml"
    shell:
        "python src/create_file_structure.py {output}"

rule check_chemical_network_data_validity:
    # snakemake -s Snakefile.smk check_chemical_network_data_validity  --cores 1 --use-conda
    input:
        lambda wildcards: [f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/{cf}.csv"
                           for cf in chemical_network_info_dict.keys()]
        #lambda wildcards: [f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/{cf}.csv" 
        #                   for cf in chemical_network_info_dict.keys()]
        #["{df}/{bn}_{cn}/{cf}.csv" for cf in chemical_network_info_dict.keys()]
        #[expand("{df}/{bn}_{cn}/{cf}.csv",
        #df = data_folder, bn = base_name, cn = padded_case_numbers, cf = chemical_network_info_dict.keys())]
    output:
        touch("{df}/{bn}_{cn}/.chemical_network_validated")
        #lambda wildcards: [f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.{cf}_dataframe_pickle"
        #                   for cf in chemical_network_info_dict.keys()]
        #["{df}/{bn}_{cn}/.{cf}_dataframe_pickle" for cf in chemical_network_info_dict.keys()]
        #[expand("{df}/{bn}_{cn}/.{cf}_dataframe_pickle",
        #df = data_folder, bn = base_name, cn = padded_case_numbers, cf = chemical_network_info_dict.keys())]
    conda:
        "config/environment.yaml"
    shell:
        "python src/check_validity_chemical_network.py {wildcards.df}/{wildcards.bn}_{wildcards.cn}"
        #    check_validity_of_csv_files(f"{df}/{bn}_{cn}/", chemical_network_info_dict.keys())

rule all:
    # snakemake -s Snakefile.smk data/case_00/.chemical_network_validated --cores 1 --use-conda
    input:
        expand("data/{bn}_{cn}/.chemical_network_validated",
               bn=["case"],
               cn=["00"])

rule create_system_network:
    # snakemake -s Snakefile create_system_network --config case_number=0 --cores 1 --use-conda
    input:
        [expand("{df}/{bn}_{cn}/.{cf}_dataframe_pickle",
        df = data_folder, bn = base_name, cn = padded_case_numbers, cf = chemical_network_info_dict.keys())]
    output:
        outputs = [file for suffix in [".NETWORK_system_pickle", "network_graph.png"]
                for file in expand("{df}/{bn}_{cn}/" + suffix,
                    df=data_folder, bn=base_name, cn=padded_case_numbers)]
    conda:
        "config/environment.yaml"
    run:
        from src.defineReactionNetwork import generate_reaction_network
        for df, bn, cn in product(data_folder, base_name, padded_case_numbers):
            generate_reaction_network(f"{df}/{bn}_{cn}/", chemical_network_csv_file_names.keys())

rule test_conda_with_output:
    # To run shell giving the output: e.g.
    # snakemake -s Snakefile.smk Alice.txt --cores 1 --use-conda
    conda:
        "config/environment.yaml"
    params:
        exclamation_mark = "True"
    output:
        "{name}.txt" # name is a wildcard
    shell:
        "python helloworld.py {wildcards.name} {params.exclamation_mark}"

rule all_tests:
    # snakemake -s Snakefile.smk --cores 1 --use-conda
    input:
        "Bob.txt", "Alice.txt"


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
        


#chemical_network_csv_file_names = chemical_network_dict.keys()
#system_parameters_csv_file_names = system_parameters_dict.keys()


#complete_file_names = complete_file_headers_dict.keys()
