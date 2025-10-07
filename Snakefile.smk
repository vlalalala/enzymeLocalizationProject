# conda activate snakemake-runner

### Import libraries from Python standard libraries ###

import os
from itertools import product
from src.auxiliary_functions_using_standard_library import as_list, load_json

### Setup import of inputs ###

configfile: "config/config.json"

data_folder = config["data_folder"]
base_name = config["base_name"]

case_numbers = as_list(config["case_numbers"], int)
digits_case_numbers = int(config["digits_case_numbers"])
padded_case_numbers = [str(n).zfill(digits_case_numbers) for n in case_numbers]

### Get template for input files about chemical network ###

chemical_network_info_dict = load_json("src/chemical_network_info.json")

### Define rules to run. The templates must already have been created and all data filled in. ###
# Create the template files through
# python src/create_file_structure.py

rule check_system_geometry_data_validity:
    # snakemake -s Snakefile.smk data/violacein_0/.chemical_network_validated --cores 1 --use-conda
    input:
        lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/SYSTEM_GEOMETRY.json"
    output:
        touch("{df}/{bn}_{cn}/.system_network_validated")
    conda:
        "config/environment.yaml"
    shell:
        "python src/check_validity_system_geometry.py {wildcards.df}/{wildcards.bn}_{wildcards.cn}"

rule check_chemical_network_data_validity:
    # snakemake -s Snakefile.smk data/violacein_0/.chemical_network_validated --cores 1 --use-conda
    input:
        lambda wildcards: [f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/{cf}.csv"
                           for cf in chemical_network_info_dict.keys()]
    output:
        touch("{df}/{bn}_{cn}/.chemical_network_validated")
    conda:
        "config/environment.yaml"
    shell:
        "python src/check_validity_chemical_network.py {wildcards.df}/{wildcards.bn}_{wildcards.cn}"

rule create_reaction_network:
    # snakemake -s Snakefile.smk data/violacein_0/network_graph.png --cores 1 --use-conda
    input:
        lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.chemical_network_validated"
    output:
        ["{df}/{bn}_{cn}/.NETWORK_system_pickle", "{df}/{bn}_{cn}/network_graph.png"]
    conda:
        "config/environment.yaml"
    shell:
        "python src/create_reaction_network.py {wildcards.df}/{wildcards.bn}_{wildcards.cn}"

rule define_kinetics_in_reaction_network:
    # snakemake -s Snakefile.smk data/violacein_0/.kinetics_defined --cores 1 --use-conda
    input:
        lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.NETWORK_system_pickle"
    output:
        touch("{df}/{bn}_{cn}/.chemical_network_completed")
    conda:
        "config/environment.yaml"
    shell:
        "python src/define_kinetics_in_reaction_network {wildcards.df}/{wildcards.bn}_{wildcards.cn}"



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

"""
rule all:
    # snakemake -s Snakefile.smk data/case_00/.chemical_network_validated --cores 1 --use-conda
    input:
        expand("data/{bn}_{cn}/.chemical_network_validated",
               bn=["case"],
               cn=["00"])
"""    


#chemical_network_csv_file_names = chemical_network_dict.keys()
#system_parameters_csv_file_names = system_parameters_dict.keys()


#complete_file_names = complete_file_headers_dict.keys()
