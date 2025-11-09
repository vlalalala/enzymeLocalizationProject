# conda activate snakemake-runner

### Import libraries from Python standard libraries ###

import os
import re
from itertools import product
from src.auxiliary_functions_using_standard_library import as_list, load_json
from src.auxiliary_functions_framework_organization_using_standard_library import(
    find_latest_solution)

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

rule cleanup_old_iterations:
    """In case the input files for a simulation have been changed, all of the
    files with .iteration_nr_* have to be deleted, as well as the log file created
    previously.
    """
    input:
        network = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.REACTION_NETWORK_pickle",
        geometry = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.SYSTEM_GEOMETRY_pickle",
        solver_info = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.solver_info_pickle",
    output:
        touch("{df}/{bn}_{cn}/.clean_iterations")
    run:
        import os, glob
        folder = f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}"
        print(f"Cleaning up {folder}")
        # delete all old iteration files
        for f in glob.glob(os.path.join(folder, "solver_iteration_data/.iteration_nr_*")):
            os.remove(f)
        log_file = os.path.join(folder, ".newton_solver.log")
        if os.path.exists(log_file):
            os.remove(log_file)
        # mark cleanup as done
        with open(output[0], "w") as f:
            f.write("done\n")

rule solve_boundary_value_problem:
    """The max-iterations condition can be changed as required without deleting anything.
    Automatically finds the latest iteration saved.
    """
    # snakemake -s Snakefile.smk data/exampleToManuallyCheck_0/.species_steady_state_concentrations.json --config max_iterations=1e6 --cores 1 --use-conda
    input:
        network = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.REACTION_NETWORK_pickle",
        geometry = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.SYSTEM_GEOMETRY_pickle",
        solver_info = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.solver_info_pickle",
        cleanup = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.clean_iterations",
    output:
        "{df}/{bn}_{cn}/.species_steady_state_concentrations.json"
    params:
        max_iterations = lambda wildcards: int(config.get("max_iterations", 1e6)),
        previous_solution = find_latest_solution
    conda:
        "config/environment.yaml"
    shell:
        """
        python src/solve_boundary_value_problem.py \
            {wildcards.df}/{wildcards.bn}_{wildcards.cn} \
            --max-iterations {params.max_iterations} \
            --previous-solution {params.previous_solution}"
        """

rule plot_boundary_value_problem_iteration:
    """Plots the concentration of the species. If the number of the iteration
    --config iteration=<n> is given, that iteration is plotted.
    Otherwise, the latest iteration file found is used.
    """
    # snakemake -s Snakefile.smk data/exampleToManuallyCheck_0/.plot_iteration --config iteration=42 --cores 1 --use-conda
    input:
        network = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.REACTION_NETWORK_pickle",
        geometry = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.SYSTEM_GEOMETRY_pickle",
        solver_info = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.solver_info_pickle",
    params:
        iteration_file = lambda wildcards: find_iteration_to_plot(wildcards, config.get("iteration", None))
    output:
        lambda w: os.path.splitext(params.iteration_file(w))[0] + ".png"
    conda:
        "config/environment.yaml"
    shell:
        """
        python src/plot_boundary_value_problem.py \
            {wildcards.df}/{wildcards.bn}_{wildcards.cn} \

        
        """


        




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

