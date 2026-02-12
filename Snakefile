# conda activate snakemake_env

### Import libraries from Python standard libraries ###

import os
import re
import glob
from itertools import product
from src.auxiliary_functions_using_standard_library import as_list, load_json
from src.auxiliary_functions_framework_organization_using_standard_library import(
    find_latest_solution)

### Get template for input files about reaction network ###

reaction_network_info_dict = load_json("src/_template_reaction_network.json")

############################################################
# Step 1: Create a new folder
# python src/_create_parameters_template.py path_to_new_folder

# Step 2: Write in all of the parameters that should be tested out

# Step 3: Create combination files
# python src/_create_templates_expanded.py path_to_new_folder
# python src/_create_phase_space.py path_to_new_folder

# Step 4: Run rule all
# snakemake --use-conda --cores 2
############################################################

df = "data/test_phase_space"
df = "examples/simple_decay_without_inner_boundaries"
df = "examples/simple_decay_with_one_inner_boundary"
sim_folders = sorted(glob.glob(os.path.join(df, "combined_*")))
all_outputs = [os.path.join(f, ".species_steady_state_concentrations.json") for f in sim_folders]

rule all:
    input:
        all_outputs

rule check_solver_info_validity:
    # snakemake -s Snakefile.smk data/test_0/.solver_input_pickle --cores 1 --use-conda
    input:
        lambda wildcards: [
            f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_solver_input.yaml",
            f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_solver_params.yaml"
        ]
    output:
        touch("{df}/{bn}_{cn}/.validated_solver")
    conda:
        "config/environment.yaml"
    shell:
        "python src/check_solver_validity.py {wildcards.df}/{wildcards.bn}_{wildcards.cn}"


rule check_system_geometry_info_validity:
    """ To define the system geometry, the number of mesh points for the solver already
    has to be read, in order to shift the membrane positions to the closest mesh positions
    """
    # snakemake -s Snakefile.smk data/test_0/.system_geometry_expanded.json --cores 1 --use-conda
    input:
        lambda wildcards: [f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_geometry.yaml",
                           f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.validated_solver",
                           ]
    output:
        "{df}/{bn}_{cn}/.expanded_system_geometry.json"
    conda:
        "config/environment.yaml"
    shell:
        "python src/check_system_geometry_validity.py {wildcards.df}/{wildcards.bn}_{wildcards.cn}"

rule check_reaction_network_info_validity:
    # snakemake -s Snakefile.smk data/test_0/.validated_reaction_network --cores 1 --use-conda
    input:
        lambda wildcards: [
            f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/{rn}.csv"
            for rn in reaction_network_info_dict.keys()
            ] + [f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.expanded_system_geometry.json"]
    output:
        touch("{df}/{bn}_{cn}/.validated_reaction_network")
    conda:
        "config/environment.yaml"
    shell:
        "python src/check_reaction_network_validity.py {wildcards.df}/{wildcards.bn}_{wildcards.cn}"

rule create_reaction_network:
    # snakemake -s Snakefile.smk data/violacein_0/reaction_network_graph.png --cores 1 --use-conda
    input:
        lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.validated_reaction_network"
    output:
        ["{df}/{bn}_{cn}/.pickled_reaction_network", "{df}/{bn}_{cn}/reaction_network_graph.png"]
    conda:
        "config/environment.yaml"
    shell:
        "python src/create_reaction_network.py {wildcards.df}/{wildcards.bn}_{wildcards.cn}"

rule create_system_mesh:
    # snakemake -s Snakefile.smk data/test_0/.expanded_system_mesh.json --cores 1 --use-conda
    input:
        lambda wildcards: [
            f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.expanded_system_geometry.json",
            f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.pickled_reaction_network",
            f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_solver_input.yaml"]
    output:
        "{df}/{bn}_{cn}/.expanded_system_mesh.json"
    conda:
        "config/environment.yaml"
    shell:
        "python src/create_system_mesh.py {wildcards.df}/{wildcards.bn}_{wildcards.cn}"

rule cleanup_old_iterations:
    """In case the input files for a simulation have been changed, all of the
    files with .iteration_nr_* have to be deleted, as well as the log file created
    previously.
    """
    input:
        network = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.pickled_reaction_network",
        geometry = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.expanded_system_geometry.json",
        solver_input = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_solver_input.yaml",
    output:
        touch("{df}/{bn}_{cn}/.validated_iterations")
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
    # snakemake -s Snakefile.smk data/test_0/.species_steady_state_concentrations.json --config max_iterations=1e4 --cores 1 --use-conda
    input:
        network = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.pickled_reaction_network",
        geometry = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.expanded_system_geometry.json",
        solver_input = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_solver_input.yaml",
        system_mesh = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.expanded_system_mesh.json",
        cleanup = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.validated_iterations",
    output:
        "{df}/{bn}_{cn}/.species_steady_state_concentrations.json"
    params:
        max_iterations = lambda wildcards: int(config.get("max_iterations", 1e5)),
        solver_params = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_solver_params.yaml"
    conda:
        "config/environment.yaml"
    shell:
        """
        python src/run_bvp_solver.py \
            {wildcards.df}/{wildcards.bn}_{wildcards.cn} \
            --max_iterations {params.max_iterations} \
            --solver_input_file {input.solver_input} \
            --solver_params_file {params.solver_params} \
        """

rule plot_boundary_value_problem:
    """Rule is not meant to be chained to other rules.
    """
    # snakemake -s Snakefile data/test_phase_space/combined_000001/.completed_visualization --cores 1 --use-conda
    output:
        touch("{df}/{bn}_{cn}/.completed_visualization")
    conda:
        "config/environment.yaml"
    shell:
        """
        python src/plot_bvp_solution.py \
            {wildcards.df}/{wildcards.bn}_{wildcards.cn} \
        """



### To have a specific iteration of the solver plotted, run on terminal ###
# python src/plot_bvp_solution.py data/test_0 --plot_iteration 40
### To have a gif of the iterations of the solver (already before the solver has converged), run on terminal ###
# python src/plot_bvp_solution.py data/test_0
