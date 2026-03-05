# conda activate snakemake_env

### Import libraries from Python standard libraries ###

import os
import re
import glob
from itertools import product
from src.auxiliary_functions_using_standard_library import as_list, load_json

#########
# TODO: new snakemake rule. inputs factor for mesh and resolves everything
# with that factor mesh. (int). takes .species_steady_state_concentrations.json
# and interpolates
#########

############################################################
# Step 1: Create a new folder
# python src/_create_parameters_template.py path_to_new_folder

# Step 2: Write in all of the parameters that should be tested out

# Step 3: Create combination files
# python src/_create_templates_expanded.py path_to_new_folder
# python src/_create_phase_space.py path_to_new_folder

# Step 4: Run rule all
# snakemake --use-conda --cores 2

# or
"""
snakemake \
  --profile config/slurm \
  --jobs 216 \
  --rerun-incomplete \
  --keep-going \
  --use-conda \
  --forcerun solve_boundary_value_problem
"""
# Get the number of jobs running through squeue --me -h | wc -l
# --keep-going stops snakemake from submitting jobs once one has not worked

############################################################

df = "data/test_phase_space"
df = "examples/simple_decay_without_inner_boundaries"
df = "examples/simple_decay_with_one_inner_boundary"
df = "examples/simple_decay_with_two_inner_boundaries"
#df = "data/simple_cycle_system"
df = "data_private/slurm_test2"
#df = "data_private/test1"

sim_folders = sorted(glob.glob(os.path.join(df, "combined_*")))
all_outputs = [os.path.join(f, ".validated_iterations") for f in sim_folders]
all_outputs = [os.path.join(f, ".species_steady_state_concentrations.json") for f in sim_folders]
all_outputs = [os.path.join(f, ".completed_visualization") for f in sim_folders]

rule all:
    input:
        all_outputs

### Get template for input files about reaction network ###
reaction_network_info_dict = load_json("src/_template_reaction_network.json")

# Maximum resident set size (kbytes) computed with 2500 mesh points, 4 species, 3 enzymes and 6 reactions.

############################################
# RULES FOR CHECKING VALIDITY OF USER INPUT
############################################

rule check_solver_input_validity:
    # snakemake -s Snakefile data/test_0/.validated_solver_input --cores 1 --use-conda
    input:
        lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_solver_input.yaml"
    output:
        touch("{df}/{bn}_{cn}/.validated_solver_input")
    threads: 1
    resources:
        mem_mb=300, # using /usr/bin/time -v gave me Maximum resident set size (kbytes): 71480
        runtime=6 # 3 minutes
    conda:
        "config/environment.yaml"
    shell:
        "python src/check_solver_validity.py {wildcards.df}/{wildcards.bn}_{wildcards.cn} parameters_solver_input"

rule check_solver_output_validity:
    # snakemake -s Snakefile data/test_0/.validated_solver_output --cores 1 --use-conda
    input:
        lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_solver_input.yaml"
    output:
        touch("{df}/{bn}_{cn}/.validated_solver_output")
    threads: 1
    resources:
        mem_mb=300, # using /usr/bin/time -v gave me Maximum resident set size (kbytes): 71480
        runtime=6 # 3 minutes
    conda:
        "config/environment.yaml"
    shell:
        "python src/check_solver_validity.py {wildcards.df}/{wildcards.bn}_{wildcards.cn} parameters_solver_output"

rule check_reaction_network_info_validity:
    # snakemake -s Snakefile data/test_0/.validated_reaction_network_input --cores 1 --use-conda
    input:
        lambda wildcards: [
            f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/{rn}.csv"
            for rn in reaction_network_info_dict.keys()
            ] + [f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_geometry.yaml"]
    output:
        touch("{df}/{bn}_{cn}/.validated_reaction_network_input")
    threads: 1
    resources:
        mem_mb=300, # using /usr/bin/time -v gave me Maximum resident set size (kbytes): 73560
        runtime=7 # 1 minute
    conda:
        "config/environment.yaml"
    shell:
        "python src/check_reaction_network_validity.py {wildcards.df}/{wildcards.bn}_{wildcards.cn}"

#################################################
# RULES FOR DEFINING SYSTEM
#################################################

rule create_system_geometry:
    """ To define the system geometry, the baseline number of mesh points for the solver already
    has to be read, in order to shift the membrane positions to the closest mesh positions
    """
    # snakemake -s Snakefile data/test_0/.system_geometry.json --cores 1 --use-conda
    input:
        lambda wildcards: [f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_geometry.yaml",
                           f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_discretization.yaml",
                           ]
    output:
        "{df}/{bn}_{cn}/.system_geometry.json"
    threads: 1
    resources:
        mem_mb=300, # using /usr/bin/time -v gave me Maximum resident set size (kbytes): 71900
        runtime=12 # 1 minute
    conda:
        "config/environment.yaml"
    shell:
        "python src/create_system_geometry.py {wildcards.df}/{wildcards.bn}_{wildcards.cn}"

rule create_reaction_network:
    # snakemake -s Snakefile data/violacein_0/.pickled_reaction_network --cores 1 --use-conda
    input:
        lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.validated_reaction_network_input"
    output:
        ["{df}/{bn}_{cn}/.pickled_reaction_network_without_enzyme_concentration"]
    threads: 1
    resources:
        mem_mb=500, # using /usr/bin/time -v gave me Maximum resident set size (kbytes): 136612
        runtime=10 # 1 minute
    conda:
        "config/environment.yaml"
    shell:
        "python src/create_reaction_network.py {wildcards.df}/{wildcards.bn}_{wildcards.cn}"

rule define_enzyme_concentrations:
    # snakemake -s Snakefile data/violacein_0/.pickled_reaction_network --cores 1 --use-conda
    input:
        lambda wildcards: [
            f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.validated_reaction_network_input",
            f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.pickled_reaction_network_without_enzyme_concentration",
            f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.system_geometry.json",
            f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_value_conditions.yaml",
        ]
    output:
        ["{df}/{bn}_{cn}/.pickled_reaction_network"]
    threads: 1
    resources:
        mem_mb=500, # using /usr/bin/time -v gave me Maximum resident set size (kbytes): 136612
        runtime=10 # 1 minute
    conda:
        "config/environment.yaml"
    shell:
        "python src/define_enzyme_concentrations.py {wildcards.df}/{wildcards.bn}_{wildcards.cn}"
    
####################################################
# RULES TO FIND AND PLOT SOLUTION
####################################################

rule cleanup_old_iterations:
    """In case any of the input files for a simulation have been changed, all of the
    files with .*iteration_nr_* have to be deleted, as well as the log file created
    previously.
    """
    input:
        discretization_yaml = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_discretization.yaml",
        geometry_yaml = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_geometry.yaml",
        solver_input_yaml = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_solver_input.yaml",
        solver_output_yaml = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_solver_output.yaml",
        value_conditions_yaml = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_value_conditions.yaml",
        geometry = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.system_geometry.json",
        network = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.pickled_reaction_network",
    output:
        touch("{df}/{bn}_{cn}/.validated_iterations")
    threads: 1
    resources:
        mem_mb=1000, # estimating based on the others.
        runtime=10 # 1 minute
    run:
        import os, glob
        folder = f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}"
        print(f"Cleaning up {folder}")
        # delete all old iteration files
        for f in glob.glob(os.path.join(folder, "solver_iteration_data/*iteration_nr_*")):
            os.remove(f)
        log_patterns = patterns = ["*.log", "*_log_*", ".*.log", ".*_log_*", ".progress_log_*"]
        log_files = []
        for pattern in patterns:
            log_files.extend(glob.glob(os.path.join(folder, pattern)))

        if log_files:
            for log_file in log_files:
                if os.path.exists(log_file):
                    os.remove(log_file)
                    print(f"Removing {log_file} file.")
        else:
            print("No log files found.")
        # mark cleanup as done
        with open(output[0], "w") as f:
            f.write("done\n")

rule solve_boundary_value_problem_with_mesh_adaptation:
    """The max-iterations condition can be changed as required without deleting anything.
    Automatically finds the latest iteration saved.
    """
    # snakemake -s Snakefile data/test_0/.species_steady_state_concentrations.json --cores 1 --use-conda
    input:
        discretization_yaml = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_discretization.yaml",
        geometry_yaml = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_geometry.yaml",
        solver_input_yaml = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_solver_input.yaml",
        solver_output_yaml = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_solver_output.yaml",
        value_conditions_yaml = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_value_conditions.yaml",
        geometry = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.system_geometry.json",
        network = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.pickled_reaction_network",
        cleanup = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.validated_iterations"
    output:
        "{df}/{bn}_{cn}/.species_steady_state_concentrations.json"
    params:
        max_num_Newton_iterations = lambda wildcards: int(config.get("max_num_Newton_iterations", 10000)),
        max_num_interpolation_times = lambda wildcards: int(config.get("max_num_interpolation_times", 8)),
        max_relative_species_concentrations_difference = lambda wildcards: config.get("max_relative_species_concentrations_difference", 1.0e-2)
    conda:
        "config/environment.yaml"
    threads: 1
    resources:
        mem_mb=1000,
        runtime= 20
    shell:
        """
        python src/run_bvp_solver_mesh_adaptation.py \
            {wildcards.df}/{wildcards.bn}_{wildcards.cn} \
            --max_num_Newton_iterations {params.max_num_Newton_iterations} \
            --max_num_interpolation_times {params.max_num_interpolation_times} \
            --max_relative_species_concentrations_difference {params.max_relative_species_concentrations_difference} \
        """


rule plot_boundary_value_problem:
    """Rule is not meant to be chained to other rules.
    """
    # snakemake -s Snakefile data/test_phase_space/combined_000001/.completed_visualization --cores 1 --use-conda
    input:
        lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/.species_steady_state_concentrations.json"
    output:
        touch("{df}/{bn}_{cn}/.completed_visualization")
    conda:
        "config/environment.yaml"
    threads: 1
    resources:
        mem_mb=1000,
        runtime= 20
    shell:
        """
        python src/plot_bvp_solution.py \
            {wildcards.df}/{wildcards.bn}_{wildcards.cn} \
        """


### To have a specific iteration of the solver plotted, run on terminal ###
# python src/plot_bvp_solution.py data/test_0 --plot_iteration 40
### To have a gif of the iterations of the solver (already before the solver has converged), run on terminal ###
# python src/plot_bvp_solution.py data/test_0


#https://collab.dvb.bayern/spaces/TUMnat/pages/431097554/SLURM+Queuing+system
#https://collab.dvb.bayern/pages/viewpage.action?spaceKey=TUMnat&title=PH+Theory+Cluster