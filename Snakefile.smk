# conda activate snakemake-runner

### Import libraries from Python standard libraries ###

import os
from itertools import product
from src.auxiliary_functions_using_standard_library import as_list, load_json

### Get template for input files about reaction network ###

reaction_network_info_dict = load_json("src/reaction_network_info.json")

### Define rules to run. The templates must already have been created and all data filled in. ###
# Create the template files through
# python src/create_file_structure.py permeability 
# python src/create_file_structure.py enzymatic

rule check_solver_data_validity:
    # snakemake -s Snakefile.smk data/violacein_0/.solver_info_pickle --cores 1 --use-conda
    input:
        lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/solver_info.json"
    output:
        touch("{df}/{bn}_{cn}/.solver_info_pickle")
    conda:
        "config/environment.yaml"
    shell:
        "python src/check_validity_solver_parameters.py {wildcards.df}/{wildcards.bn}_{wildcards.cn}"


rule check_system_geometry_data_validity:
    """ To define the system geometry, the number of mesh points for the solver already
    has to be read, in order to shift the membrane positions to the closest mesh positions
    """
    # snakemake -s Snakefile.smk data/violacein_0/.SYSTEM_GEOMETRY_pickle --cores 1 --use-conda
    input:
        lambda wildcards: [f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/SYSTEM_GEOMETRY.json",
                           f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/solver_info.json"]
    output:
        touch("{df}/{bn}_{cn}/.SYSTEM_GEOMETRY_pickle")
    conda:
        "config/environment.yaml"
    shell:
        "python src/check_validity_system_geometry.py {wildcards.df}/{wildcards.bn}_{wildcards.cn}"

rule check_reaction_network_data_validity:
    # snakemake -s Snakefile.smk data/violacein_0/.reaction_network_validated --cores 1 --use-conda
    input:
        lambda wildcards: [f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/{rn}.csv"
                           for rn in reaction_network_info_dict.keys()] + [f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.SYSTEM_GEOMETRY_pickle"]
    output:
        touch("{df}/{bn}_{cn}/.reaction_network_validated")
    conda:
        "config/environment.yaml"
    shell:
        "python src/check_validity_reaction_network.py {wildcards.df}/{wildcards.bn}_{wildcards.cn}"

rule create_reaction_network:
    # snakemake -s Snakefile.smk data/violacein_0/reaction_network_graph.png --cores 1 --use-conda
    input:
        lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.reaction_network_validated"
    output:
        ["{df}/{bn}_{cn}/.REACTION_NETWORK_pickle", "{df}/{bn}_{cn}/reaction_network_graph.png"]
    conda:
        "config/environment.yaml"
    shell:
        "python src/create_reaction_network.py {wildcards.df}/{wildcards.bn}_{wildcards.cn}"

rule solve_boundary_value_problem:
    # snakemake -s Snakefile.smk data/violacein_0/.species_steady_state_concentrations.json --cores 1 --use-conda
    input:
        lambda wildcards: [
            f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.REACTION_NETWORK_pickle",
            f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.SYSTEM_GEOMETRY_pickle",
            f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.solver_info_pickle"]
    output:
        "{df}/{bn}_{cn}/.species_steady_state_concentrations.json"
    conda:
        "config/environment.yaml"
    shell:
        "python src/solve_boundary_value_problem.py {wildcards.df}/{wildcards.bn}_{wildcards.cn}"




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

