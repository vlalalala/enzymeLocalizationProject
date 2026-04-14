# conda activate snakemake_env

### Import libraries from Python standard libraries ###

import os
import re
import glob
import hashlib
import shutil
from itertools import product
from src.auxiliary_functions_using_standard_library import as_list, load_json

############################################################

# Step 1: Create a new folder
# python src/_create_parameters_template.py path_to_new_folder

# Step 2: Write in all of the parameters that should be tested out

# Step 3: Create combination files
# python src/_create_templates_expanded.py path_to_new_folder
# python src/_create_phase_space.py path_to_new_folder

# Step 4: Run rule all
# snakemake --use-conda --cores 10 --scheduler greedy --latency-wait 120 --retries 2 --config fast_mode=true
# scheduler greedy to schedule short jobs first
# or
"""
snakemake \
  --profile config/slurm \
  --jobs 50 \
  --keep-going \
  --use-conda \
  --latency-wait 120 \
  --retries 2 \
  --config fast_mode=true
"""
# Get the number of jobs running through squeue --me -h | wc -l
# --keep-going stops snakemake from submitting jobs once one has not worked

# Note: The maximum resident set size (kbytes) was computed with 2500 mesh points, 4 species, 3 enzymes and 6 reactions.

# Reminder: du -sh . for getting current directory size

############################################################
# Create environment
#rule create_environment:
#    # snakemake -s Snakefile config/.environment_with_snakemake_created --cores 1 --use-conda
#    output:
#        touch("config/.environment_with_snakemake_created")
#    conda:
#        "config/environment_with_snakemake.yaml"

############################################################

FAST_MODE = config.get("fast_mode", False)

print(f"Running with FAST_MODE: {FAST_MODE} (type {type(FAST_MODE)}).")

df = "data_private/optimizationSpeedTest_spontaneousXdecaysToY_1InnerBoundary"
df = "data_private/optimization_spontaneousXYZchain_1InnerBoundary"
df = "data_private/optimization_enzymaticXtoY_spontaneousXtoZ_1InnerBoundary"

if not os.path.isdir(df):
    raise ValueError(df, "does not exist.")

sim_folders = sorted(glob.glob(os.path.join(df, "combined_*")))


# Different outputs dependent on mode
# Without optimization
all_outputs = [os.path.join(f, ".validated_iterations") for f in sim_folders]
all_outputs = [os.path.join(f, ".species_steady_state_concentrations.json") for f in sim_folders]
all_outputs = [os.path.join(f, ".completed_visualization") for f in sim_folders]

# With optimization
# all_outputs = [os.path.join(f, "best_result.json") for f in sim_folders]
all_outputs = [os.path.join(f, ".modifications_of_best_result_created") for f in sim_folders]
all_outputs = [os.path.join(f, "optimization_modifications_analysis.json") for f in sim_folders]

rule all:
    input:
        all_outputs

########################################
# MODE CONTROL
########################################

def trial_path(wildcards, filename=""):
    """Returns the correct folder path depending on whether the
    optimal solution should be found or not 
    (if in optimization mode there are round/trial wildcards,
     or modification mode there are check/modification wildcards)
    """
    round_ = getattr(wildcards, "round", None)
    trial  = getattr(wildcards, "trial", None)
    mod    = getattr(wildcards, "modification", None)

    if round_ is not None and trial is not None:
        base = (f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/"
                f"optimization_round_{round_}/trial_{trial}")
    elif mod is not None:
        base = (f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/"
                f"optimization_check/modification_{mod}")
    else:
        base = f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}"

    return f"{base}/{filename}" if filename else base

N_TRIALS = 20    # parallel solver calls per round
N_ROUNDS = 20   # optimization rounds
ROUND_WIDTH = len(str(N_ROUNDS - 1))
TRIAL_WIDTH = len(str(N_TRIALS - 1))
#TRIALS = [str(i).zfill(TRIAL_WIDTH) for i in range(N_TRIALS)]
#ROUNDS = [str(i).zfill(ROUND_WIDTH) for i in range(N_ROUNDS)]
TRIALS = list(range(N_TRIALS))
ROUNDS = list(range(N_ROUNDS))

FILE_NAMES = [
    "parameters_discretization.yaml",
    "parameters_geometry.yaml",
    "parameters_solver_input.yaml",
    "parameters_solver_output.yaml",
    "parameters_value_conditions.yaml",
    "enzymes.csv",
    "enzymatic_reactions.csv",
    "species.csv",
    "spontaneous_reactions.csv"
]

wildcard_constraints:
#    #round = r"\d{" + str(ROUND_WIDTH) + r"}",  # exactly 2 digits: 00-99,
#    #trial = r"\d{" + str(TRIAL_WIDTH) + r"}",
#    bn    = r"[a-zA-Z]+",
    cn    = r"\d{6}", # cn has 6 digits
#    df    = r"(?!.*optimization).*"  # negative lookahead: no "optimization" anywhere in df


############################################
# RULES FOR CHECKING VALIDITY OF USER INPUT
############################################

### Get template for input files about reaction network ###
reaction_network_info_dict = load_json("src/_template_reaction_network.json")

if FAST_MODE:
    rule get_result_from_complete_input:
        # snakemake -s Snakefile data_private/errorTest_enzymatic_allocation_wrong/combined_000001/.result_computed_fast --cores 1 --use-conda --config fast_mode=true
        input:
            lambda wildcards: trial_path(wildcards, "parameters_solver_input.yaml"),
            lambda wildcards: trial_path(wildcards, "parameters_solver_output.yaml"),
            lambda wildcards: [
                    trial_path(wildcards, f"{rn}.csv")
                    for rn in reaction_network_info_dict.keys()
                ],
            lambda wildcards: trial_path(wildcards, "parameters_geometry.yaml"),
            lambda wildcards: trial_path(wildcards, "parameters_discretization.yaml"),
            lambda wildcards: trial_path(wildcards, "parameters_value_conditions.yaml")
        output:
            touch("{df}/{bn}_{cn}/.result_computed_fast")
        params:
            folder = lambda wildcards: trial_path(wildcards),
            max_num_Newton_iterations_creeping_reaction    = lambda wildcards: int(config.get("max_num_Newton_iterations", 2500)),
            max_num_creeping_reaction_simulations          = lambda wildcards: int(config.get("max_num_creeping_reaction_simulations", 3)),
            override_creeping_reaction_simulations         = True,
            max_num_Newton_iterations                      = lambda wildcards: int(config.get("max_num_Newton_iterations", 1000)),
            max_num_interpolation_times                    = lambda wildcards: int(config.get("max_num_interpolation_times", 3)),
            max_relative_species_concentrations_difference = lambda wildcards: config.get("max_relative_species_concentrations_difference", 1.0e-2),
            max_relative_flux_difference                   = lambda wildcards: config.get("max_relative_flux_difference", 1.0e-2),
            min_relative_concentration_difference_considered_relevant = lambda wildcards: config.get("min_relative_concentration_difference_considered_relevant", 1.0e-2)
        threads: 1
        resources:
            mem_mb=1000,
            runtime=10 # runtime resource for whole group
        group: lambda wildcards: f"sp_{hashlib.md5(trial_path(wildcards).encode()).hexdigest()[:8]}"
        priority:
            100  # run first
        conda:
            "config/environment.yaml"
        shell:
            """
            python src/check_solver_validity.py {params.folder} parameters_solver_input
            python src/check_solver_validity.py {params.folder} parameters_solver_output
            python src/check_reaction_network_validity.py {params.folder}
            python src/create_system_geometry.py {params.folder}
            python src/create_reaction_network.py {params.folder}
            python src/define_enzyme_concentrations.py {params.folder}
            python src/run_bvp_solver_initial_guess_calculator.py \
                --folder {params.folder} \
                --max_num_Newton_iterations {params.max_num_Newton_iterations_creeping_reaction} \
                --max_num_creeping_reaction_simulations {params.max_num_creeping_reaction_simulations} \
                --override {params.override_creeping_reaction_simulations}
            python src/run_bvp_solver_mesh_adaptation.py \
                --folder {params.folder} \
                --max_num_Newton_iterations {params.max_num_Newton_iterations} \
                --max_num_interpolation_times {params.max_num_interpolation_times} \
                --max_relative_species_concentrations_difference {params.max_relative_species_concentrations_difference} \
                --max_relative_flux_difference {params.max_relative_flux_difference} \
                --min_relative_concentration_difference_considered_relevant {params.min_relative_concentration_difference_considered_relevant}
            python src/study_bvp_solution.py {params.folder}/
            """
   
    use rule get_result_from_complete_input as get_result_from_complete_input_within_optimization with:
        output:
            touch("{df}/{bn}_{cn}/optimization_round_{round}/trial_{trial}/.result_computed_fast")

else:
    rule check_solver_input_validity:
        # snakemake -s Snakefile data/test_0/.validated_solver_input --cores 1 --use-conda
        input:
            lambda wildcards: trial_path(wildcards, "parameters_solver_input.yaml")
        output:
            touch("{df}/{bn}_{cn}/.validated_solver_input")
        params:
            folder = lambda wildcards: trial_path(wildcards)
        threads: 1
        resources:
            mem_mb=300,
            runtime=20 # runtime resource for whole group
        group: lambda wildcards: f"sp_{hashlib.md5(trial_path(wildcards).encode()).hexdigest()[:8]}"
        priority:
            100  # run first
        conda:
            "config/environment.yaml"
        shell:
            "python src/check_solver_validity.py {params.folder} parameters_solver_input"

    use rule check_solver_input_validity as check_solver_input_validity_within_optimization with:
        group: lambda wildcards: f"sp_{hashlib.md5(trial_path(wildcards).encode()).hexdigest()[:8]}"
        output:
            touch("{df}/{bn}_{cn}/optimization_round_{round}/trial_{trial}/.validated_solver_input")

    rule check_solver_output_validity:
        # snakemake -s Snakefile data/test_0/.validated_solver_output --cores 1 --use-conda
        input:
            lambda wildcards: trial_path(wildcards, "parameters_solver_output.yaml")
        output:
            touch("{df}/{bn}_{cn}/.validated_solver_output")
        params:
            folder = lambda wildcards: trial_path(wildcards)
        threads: 1
        resources:
            mem_mb=300,
            runtime=5
        group: lambda wildcards: f"sp_{hashlib.md5(trial_path(wildcards).encode()).hexdigest()[:8]}"
        priority:
            100  # run first
        conda:
            "config/environment.yaml"
        shell:
            "python src/check_solver_validity.py {params.folder} parameters_solver_output"

    use rule check_solver_output_validity as check_solver_output_validity_within_optimization with:
        group: lambda wildcards: f"sp_{hashlib.md5(trial_path(wildcards).encode()).hexdigest()[:8]}"
        output:
            touch("{df}/{bn}_{cn}/optimization_round_{round}/trial_{trial}/.validated_solver_output")

    rule check_reaction_network_info_validity:
        # snakemake -s Snakefile data/test_0/.validated_reaction_network_input --cores 1 --use-conda
        input:
            lambda wildcards: [
                trial_path(wildcards, f"{rn}.csv")
                for rn in reaction_network_info_dict.keys()
            ] + [trial_path(wildcards, "parameters_geometry.yaml")]
        output:
            touch("{df}/{bn}_{cn}/.validated_reaction_network_input")
        params:
            folder = lambda wildcards: trial_path(wildcards)
        threads: 1
        resources:
            mem_mb=300,
            runtime=5
        group: lambda wildcards: f"sp_{hashlib.md5(trial_path(wildcards).encode()).hexdigest()[:8]}"
        priority:
            100  # run first
        conda:
            "config/environment.yaml"
        shell:
            "python src/check_reaction_network_validity.py {params.folder}"

    use rule check_reaction_network_info_validity as check_reaction_network_info_validity_within_optimization with:
        group: lambda wildcards: f"sp_{hashlib.md5(trial_path(wildcards).encode()).hexdigest()[:8]}"
        output:
            touch("{df}/{bn}_{cn}/optimization_round_{round}/trial_{trial}/.validated_reaction_network_input")

    #################################################
    # RULES FOR DEFINING SYSTEM
    #################################################

    rule create_system_geometry:
        """ To define the system geometry, the baseline number of mesh points for the solver already
        has to be read, in order to shift the membrane positions to the closest mesh positions
        """
        # snakemake -s Snakefile data/test_0/.system_geometry.json --cores 1 --use-conda
        input:
            lambda wildcards: [
                trial_path(wildcards, "parameters_geometry.yaml"),
                trial_path(wildcards, "parameters_discretization.yaml"),
            ]
        output:
            "{df}/{bn}_{cn}/.system_geometry.json"
        params:
            folder = lambda wildcards: trial_path(wildcards)
        threads: 1
        resources:
            mem_mb=300,
            runtime=12
        group: lambda wildcards: f"sp_{hashlib.md5(trial_path(wildcards).encode()).hexdigest()[:8]}"
        priority:
            100  # run first
        conda:
            "config/environment.yaml"
        shell:
            "python src/create_system_geometry.py {params.folder}"

    use rule create_system_geometry as create_system_geometry_within_optimization with:
        group: lambda wildcards: f"sp_{hashlib.md5(trial_path(wildcards).encode()).hexdigest()[:8]}"
        output:
            "{df}/{bn}_{cn}/optimization_round_{round}/trial_{trial}/.system_geometry.json"

    rule create_reaction_network:
        # snakemake -s Snakefile data/violacein_0/.pickled_reaction_network --cores 1 --use-conda
        input:
            lambda wildcards: trial_path(wildcards, ".validated_reaction_network_input")
        output:
            "{df}/{bn}_{cn}/.pickled_reaction_network_without_enzyme_concentration"
        params:
            folder = lambda wildcards: trial_path(wildcards)
        threads: 1
        resources:
            mem_mb=500,
            runtime=10
        priority:
            100  # run first
        group: lambda wildcards: f"sp_{hashlib.md5(trial_path(wildcards).encode()).hexdigest()[:8]}"
        conda:
            "config/environment.yaml"
        shell:
            "python src/create_reaction_network.py {params.folder}"

    use rule create_reaction_network as create_reaction_network_within_optimization with:
        group: lambda wildcards: f"sp_{hashlib.md5(trial_path(wildcards).encode()).hexdigest()[:8]}"
        output:
            "{df}/{bn}_{cn}/optimization_round_{round}/trial_{trial}/.pickled_reaction_network_without_enzyme_concentration"

    rule define_enzyme_concentrations:
        # snakemake -s Snakefile data/violacein_0/.pickled_reaction_network --cores 1 --use-conda
        input:
            lambda wildcards: [
                trial_path(wildcards, ".validated_reaction_network_input"),
                trial_path(wildcards, ".pickled_reaction_network_without_enzyme_concentration"),
                trial_path(wildcards, ".system_geometry.json"),
                trial_path(wildcards, "parameters_value_conditions.yaml"),
            ]
        output:
            [
                "{df}/{bn}_{cn}/.pickled_reaction_network",
                "{df}/{bn}_{cn}/enzyme_concentrations.json",
            ]
        params:
            folder = lambda wildcards: trial_path(wildcards)
        threads: 1
        resources:
            mem_mb=500,
            runtime=10
        priority:
            100  # run first
        group: lambda wildcards: f"sp_{hashlib.md5(trial_path(wildcards).encode()).hexdigest()[:8]}"
        conda:
            "config/environment.yaml"
        shell:
            "python src/define_enzyme_concentrations.py {params.folder}"

    use rule define_enzyme_concentrations as define_enzyme_concentrations_within_optimization with:
        # hashlib to keep it shorter: group: lambda wildcards: f"solver_preparation_{trial_path(wildcards).replace("/", "_")}"
        group: lambda wildcards: f"sp_{hashlib.md5(trial_path(wildcards).encode()).hexdigest()[:8]}"
        output:
            [
                "{df}/{bn}_{cn}/optimization_round_{round}/trial_{trial}/.pickled_reaction_network",
                "{df}/{bn}_{cn}/optimization_round_{round}/trial_{trial}/enzyme_concentrations.json"
            ]
    
    ####################################################
    # RULES TO FIND AND PLOT SOLUTION
    ####################################################

    rule create_initial_guess:
        # snakemake -s Snakefile data_private/reaction_scaling_test/combined_000001/species_initial_guess.json --cores 1 --use-conda
        input:
            discretization_yaml   = lambda wildcards: trial_path(wildcards, "parameters_discretization.yaml"),
            geometry_yaml         = lambda wildcards: trial_path(wildcards, "parameters_geometry.yaml"),
            solver_input_yaml     = lambda wildcards: trial_path(wildcards, "parameters_solver_input.yaml"),
            solver_output_yaml    = lambda wildcards: trial_path(wildcards, "parameters_solver_output.yaml"),
            value_conditions_yaml = lambda wildcards: trial_path(wildcards, "parameters_value_conditions.yaml"),
            geometry              = lambda wildcards: trial_path(wildcards, ".system_geometry.json"),
            network               = lambda wildcards: trial_path(wildcards, ".pickled_reaction_network"),
        output:
            "{df}/{bn}_{cn}/species_initial_guess.json"
        params:
            folder                                         = lambda wildcards: trial_path(wildcards),
            max_num_Newton_iterations                      = lambda wildcards: int(config.get("max_num_Newton_iterations", 2500)),
            max_num_creeping_reaction_simulations          = lambda wildcards: int(config.get("max_num_creeping_reaction_simulations", 3)),
            override                                       = True
        conda:
            "config/environment.yaml"
        threads: 1
        resources:
            mem_mb=5000,
            runtime=300
        priority:
            0  # run LAST
        shell:
            """
            python src/run_bvp_solver_initial_guess_calculator.py \
                --folder {params.folder} \
                --max_num_Newton_iterations {params.max_num_Newton_iterations} \
                --max_num_creeping_reaction_simulations {params.max_num_creeping_reaction_simulations} \
                --override {params.override}
            """

    use rule create_initial_guess as create_initial_guess_within_optimization with:
        # snakemake -s Snakefile data_private/case_01/combined_000015/optimization_round_0/trial_0/species_initial_guess.json --cores 1 --use-conda
        output:
            "{df}/{bn}_{cn}/optimization_round_{round}/trial_{trial}/species_initial_guess.json"

    rule cleanup_old_iterations:
        """In case any of the input files for a simulation have been changed, all of the
        files with .*iteration_nr_* have to be deleted, as well as the log file created
        previously.
        (This is needed because those files are not outputs of any previous snakemake rule,
        so snakemake does not automatically delete then when some input file is modified.)
        "concentration_convergence", "flux_convergence_without_concentration_convergence",
        "no_flux_nor_concentration_convergence"
        """
        input:
            discretization_yaml   = lambda wildcards: trial_path(wildcards, "parameters_discretization.yaml"),
            geometry_yaml         = lambda wildcards: trial_path(wildcards, "parameters_geometry.yaml"),
            solver_input_yaml     = lambda wildcards: trial_path(wildcards, "parameters_solver_input.yaml"),
            solver_output_yaml    = lambda wildcards: trial_path(wildcards, "parameters_solver_output.yaml"),
            value_conditions_yaml = lambda wildcards: trial_path(wildcards, "parameters_value_conditions.yaml"),
            geometry              = lambda wildcards: trial_path(wildcards, ".system_geometry.json"),
            network               = lambda wildcards: trial_path(wildcards, ".pickled_reaction_network"),
        output:
            touch("{df}/{bn}_{cn}/.validated_iterations")
        params:
            folder = lambda wildcards: trial_path(wildcards)
        threads: 1
        resources:
            mem_mb=1000,
            runtime=5
        priority:
            100  # run first
        run:
            import os, glob
            folder = params.folder
            print(f"Cleaning up {folder}")
            # Remove files saved during previous simulations with the non-modified rates
            solver_iteration_data_dir = os.path.join(folder, "solver_iteration_data")
            if os.path.isdir(solver_iteration_data_dir):
                print(f"Removing {solver_iteration_data_dir} directory.")
                shutil.rmtree(solver_iteration_data_dir)
            # Remove files saved during previous simulations with the modified rates
            creeping_reaction_simulations_dir = os.path.join(folder, "creeping_reaction_simulations")
            if os.path.isdir(creeping_reaction_simulations_dir):
                print(f"Removing {creeping_reaction_simulations_dir} directory.")
                shutil.rmtree(creeping_reaction_simulations_dir)
            # Remove the final result from the simulations with modified rates
            initial_guess_file = os.path.join(folder, "species_initial_guess.json")
            if os.path.isfile(initial_guess_file):
                os.remove(initial_guess_file)
            # Remove all log files from previous (now obsolete) simulations
            log_patterns = ["*.log", "*_log_*", ".*.log", ".*_log_*", ".progress_log_*"]
            log_files = []
            for pattern in log_patterns:
                log_files.extend(glob.glob(os.path.join(folder, pattern)))
            if log_files:
                for log_file in log_files:
                    if os.path.exists(log_file):
                        os.remove(log_file)
                        print(f"Removing {log_file} file.")
            else:
                print("No log files found.")
            # Remove all plots and files related to them
            for file in os.listdir(folder):
                if any(filename in file for filename in ["newton_iterations.gif", "max_y", "convergence.png"]):
                    os.remove(os.path.join(folder, file))
            # Remove all files that inform about type of convergence
            for file in os.listdir(folder):
                if any(filename in file for filename in [
                    "concentration_convergence", "flux_convergence_without_concentration_convergence",
                    "no_flux_nor_concentration_convergence"
                ]):
                    os.remove(os.path.join(folder, file))
            # Write completed file with current metadata
            with open(output[0], "w") as f:
                f.write("done\n")

    use rule cleanup_old_iterations as cleanup_old_iterations_within_optimization with:
        output:
            touch("{df}/{bn}_{cn}/optimization_round_{round}/trial_{trial}/.validated_iterations")


    rule solve_boundary_value_problem_with_mesh_adaptation:
        # The max-iterations condition can be changed as required without deleting anything.
        # Automatically finds the latest iteration saved.
        # snakemake -s Snakefile data/test_0/.species_steady_state_concentrations.json --cores 1 --use-conda
        input:
            discretization_yaml   = lambda wildcards: trial_path(wildcards, "parameters_discretization.yaml"),
            geometry_yaml         = lambda wildcards: trial_path(wildcards, "parameters_geometry.yaml"),
            solver_input_yaml     = lambda wildcards: trial_path(wildcards, "parameters_solver_input.yaml"),
            solver_output_yaml    = lambda wildcards: trial_path(wildcards, "parameters_solver_output.yaml"),
            value_conditions_yaml = lambda wildcards: trial_path(wildcards, "parameters_value_conditions.yaml"),
            geometry              = lambda wildcards: trial_path(wildcards, ".system_geometry.json"),
            network               = lambda wildcards: trial_path(wildcards, ".pickled_reaction_network"),
            cleanup               = lambda wildcards: trial_path(wildcards, ".validated_iterations"),
            initial_guess_concentrations = lambda wildcards: trial_path(wildcards, "species_initial_guess.json")
        output:
            "{df}/{bn}_{cn}/.species_steady_state_concentrations.json"
        params:
            folder                                         = lambda wildcards: trial_path(wildcards),
            max_num_Newton_iterations                      = lambda wildcards: int(config.get("max_num_Newton_iterations", 1000)),
            max_num_interpolation_times                    = lambda wildcards: int(config.get("max_num_interpolation_times", 3)),
            max_relative_species_concentrations_difference = lambda wildcards: config.get("max_relative_species_concentrations_difference", 1.0e-2),
            max_relative_flux_difference = lambda wildcards: config.get("max_relative_flux_difference", 1.0e-2),
            min_relative_concentration_difference_considered_relevant = lambda wildcards: config.get("min_relative_concentration_difference_considered_relevant", 1.0e-2)
        conda:
            "config/environment.yaml"
        threads: 1
        resources:
            mem_mb=5000,
            runtime=15
        priority:
            0  # run LAST
        shell:
            """
            python src/run_bvp_solver_mesh_adaptation.py \
                --folder {params.folder} \
                --max_num_Newton_iterations {params.max_num_Newton_iterations} \
                --max_num_interpolation_times {params.max_num_interpolation_times} \
                --max_relative_species_concentrations_difference {params.max_relative_species_concentrations_difference} \
                --max_relative_flux_difference {params.max_relative_flux_difference} \
                --min_relative_concentration_difference_considered_relevant {params.min_relative_concentration_difference_considered_relevant}
            """

    use rule solve_boundary_value_problem_with_mesh_adaptation as solve_boundary_value_problem_with_mesh_adaptation_within_optimization with:
        # snakemake -s Snakefile data_private/case_01/combined_000015/optimization_round_0/trial_0/.species_steady_state_concentrations.json --cores 1 --use-conda
        output:
            [
            "{df}/{bn}_{cn}/optimization_round_{round}/trial_{trial}/.species_steady_state_concentrations.json",
            "{df}/{bn}_{cn}/optimization_round_{round}/trial_{trial}/system_geometry_for_convergence.json",
            "{df}/{bn}_{cn}/optimization_round_{round}/trial_{trial}/.expanded_system_mesh_for_convergence.json",
            ]


    rule study_bvp_solution:
        input:
            lambda wildcards: trial_path(wildcards, ".species_steady_state_concentrations.json")
        output:
            "{df}/{bn}_{cn}/fluxes.json"
        params:
            folder = lambda wildcards: trial_path(wildcards)
        threads: 1
        resources:
            mem_mb=1000,
            runtime=5
        priority:
            100  # run first
        conda:
            "config/environment.yaml"
        shell:
            "python src/study_bvp_solution.py {params.folder}/"

    use rule study_bvp_solution as study_bvp_solution_within_optimization with:
        output:
            "{df}/{bn}_{cn}/optimization_round_{round}/trial_{trial}/fluxes.json"

    rule plot_boundary_value_problem:
        """Rule is not meant to be chained to other rules.
        """
        # snakemake -s Snakefile data/test_phase_space/combined_000001/.completed_visualization --cores 1 --use-conda
        input:
            lambda wildcards: trial_path(wildcards, ".species_steady_state_concentrations.json")
        output:
            touch("{df}/{bn}_{cn}/.completed_visualization")
        params:
            folder = lambda wildcards: trial_path(wildcards)
        conda:
            "config/environment.yaml"
        threads: 1
        resources:
            mem_mb=1000,
            runtime= 20
        priority:
            100  # run first
        shell:
            """
            python src/plot_bvp_solver_mesh_adaptation_progress.py {params.folder}/
            """

    use rule plot_boundary_value_problem as plot_boundary_value_problem_within_optimization with:
        # snakemake -s Snakefile data_private/case_01/combined_000015/optimization_round_0/trial_0/.completed_visualization --cores 1 --use-conda
        output:
            touch("{df}/{bn}_{cn}/optimization_round_{round}/trial_{trial}/.completed_visualization")



#######################################################
# ORGANIZE OPTIMIZATION
#######################################################

# --- Round 0: suggest from base folder (no previous .done file) ---
rule suggest_optimization_params_round_0:
    # snakemake -s Snakefile data_private/optimization_spontaneousXdecaysToY_1InnerBoundary/combined_000001/optimization_round_0/.trial_files_created --cores 1 --use-conda
    wildcard_constraints:
        round = r"0+"  # only zeros: "0" or "00" etc.
    input:
        discretization_yaml       = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_discretization.yaml",
        geometry_yaml             = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_geometry.yaml",
        solver_input_yaml         = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_solver_input.yaml",
        solver_output_yaml        = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_solver_output.yaml",
        value_conditions_yaml     = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/parameters_value_conditions.yaml",
        enzymes_csv               = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/enzymes.csv",
        enzymatic_reactions_csv   = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/enzymatic_reactions.csv",
        species_csv               = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/species.csv",
        spontaneous_reactions_csv = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/spontaneous_reactions.csv",
    output:
        files = expand("{{df}}/{{bn}}_{{cn}}/optimization_round_{round}/trial_{trial}/{file_name}",
                    round=ROUNDS[0],
                    trial=TRIALS,
                    file_name=FILE_NAMES),
        flag  = touch("{df}/{bn}_{cn}/optimization_round_" + str(ROUNDS[0]) + "/.trial_files_created")
    conda:
        "config/environment.yaml"
    threads: 1
    resources:
        mem_mb=1000,
        runtime= 5
    priority:
        101  # run first
    shell:
        """
        python src/optimizer_suggest_trials.py \
            --folder_to_solve {wildcards.df}/{wildcards.bn}_{wildcards.cn} \
            --round 0 \
            --n_trials {N_TRIALS} \
            --n_rounds {N_ROUNDS}
        """

# --- Rounds 1..N-1: suggest from previous round's .done file ---
rule suggest_optimization_params:
    wildcard_constraints:
        round = r"0*[1-9]\d*"  # any number except 0 (with or without leading zeros)
    input:
        previous_round_best = lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/optimization_round_{int(wildcards.round)-1}_best.json"
    output:
        files = expand("{{df}}/{{bn}}_{{cn}}/optimization_round_{{round}}/trial_{trial}/{file_name}",
                    trial=TRIALS,
                    file_name=FILE_NAMES),
        flag  = touch("{df}/{bn}_{cn}/optimization_round_{round}/.trial_files_created")
    conda:
        "config/environment.yaml"
    priority:
        101  # run first
    threads: 1
    resources:
        mem_mb=1000,
        runtime= 5
    shell:
        """
        python src/optimizer_suggest_trials.py \
            --folder_to_solve {wildcards.df}/{wildcards.bn}_{wildcards.cn} \
            --round {wildcards.round} \
            --n_trials {N_TRIALS} \
            --n_rounds {N_ROUNDS}
        """

def run_result_input(wildcards):
    if FAST_MODE:
        return expand(
            f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/optimization_round_{wildcards.round}/trial_{{trial}}/.result_computed_fast",
            trial=TRIALS
        )
    else:
        return expand(
            f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/optimization_round_{wildcards.round}/trial_{{trial}}/fluxes.json",
            trial=TRIALS
        )

rule collect_and_update:
    # snakemake -s Snakefile data_private/optuna_test/combined_000001/optimization_round_4.done --cores 1 --use-conda
    input:
        results = run_result_input,
        params_opt = lambda wildcards: f"{wildcards.df}/parameters_optimization.yaml"
    output:
        "{df}/{bn}_{cn}/optimization_round_{round}_best.json"
    conda:
        "config/environment.yaml"
    threads: 1
    resources:
        mem_mb=1000,
        runtime= 5
    priority:
        100  # run first
    shell:
        """
        python src/optimizer_collect_trial_results.py \
            --folder_to_solve {wildcards.df}/{wildcards.bn}_{wildcards.cn} \
            --round {wildcards.round} \

        python src/optimizer_get_best_result.py \
            --folder_to_solve {wildcards.df}/{wildcards.bn}_{wildcards.cn} \
            --round {wildcards.round}
        """

rule best_result:
    # snakemake -s Snakefile data_private/optuna_test/combined_000001/best_result.json --cores 1 --use-conda
    input:
        lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/optimization_round_{ROUNDS[-1]}_best.json"
    output:
        ["{df}/{bn}_{cn}/best_result.json",
        "{df}/{bn}_{cn}/optimization_progress.png"
        ]
    conda:
        "config/environment.yaml"
    params:
        last_round = N_ROUNDS-1
    threads: 1
    resources:
        mem_mb=1000,
        runtime= 3
    priority:
        100  # run first
    shell:
        """
        python src/optimizer_get_best_result.py \
            --folder_to_solve {wildcards.df}/{wildcards.bn}_{wildcards.cn} \
            --round {params.last_round}

        python src/optimizer_plot_progress.py {wildcards.df}/{wildcards.bn}_{wildcards.cn} {params.last_round}
        """

############################################################
# CHECK OPTIMIZATION
############################################################
# Number of total files to create: not combinatorics!
# (changing each parameter and adding the options all together)

# Factor of n_inner_membranes * 2: move each membrane to one mesh point to the left and one to the right

# Factor of n_enzymes * 2: multiply the allocation of each enzyme by 0.99 and 1.01
# (making sure to then normalize total quantity)

# Allocation of enzymes to different regions: n_enzymes * n_regions * 2:
# multiply the allocation of each enzyme to each region by 0.99 and 1.01 
# (again, making sure to then normalize the total quantity)

NUMBER_ENZYMES = 0
NUMBER_INNER_MEMBRANES = 1
TOTAL_NUMBER_MODIFICATIONS = 2 * (
    NUMBER_INNER_MEMBRANES
    + NUMBER_ENZYMES
    + NUMBER_ENZYMES * (NUMBER_INNER_MEMBRANES + 1)
)

MODIFICATIONS =  list(range(TOTAL_NUMBER_MODIFICATIONS))

rule create_modifications_of_best_result:
    # snakemake -s Snakefile data_private/optuna_test/combined_000001/.modifications_of_best_result_created --cores 1 --use-conda
    input:
        lambda wildcards: f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/best_result.json"
    output:
        touch("{df}/{bn}_{cn}/.modifications_of_best_result_created"),
        files = expand("{{df}}/{{bn}}_{{cn}}/optimization_check/modification_{modification}/{file_name}",
                    modification=MODIFICATIONS,
                    file_name=FILE_NAMES),
    params:
        expected_number_modifications = TOTAL_NUMBER_MODIFICATIONS
    conda:
        "config/environment.yaml"
    threads: 1
    resources:
        mem_mb=1000,
        runtime= 3
    priority:
        100  # run first
    shell:
        """
        python src/optimizer_create_modifications_to_extremum.py \
            --folder {wildcards.df}/{wildcards.bn}_{wildcards.cn} \
            --expected_number_modifications {params.expected_number_modifications}
        """    

def run_result_input_modifications(wildcards):
    if FAST_MODE:
        return expand(
            f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/optimization_check/modification_{{modification}}/.result_computed_fast",
            modification=MODIFICATIONS
        )
    else:
        return expand(
            f"{wildcards.df}/{wildcards.bn}_{wildcards.cn}/optimization_check/modification_{{modification}}/fluxes.json",
            modification=MODIFICATIONS
        )

rule analyze_results_of_modifications_of_best_result:
    # snakemake -s Snakefile data_private/optimization_spontaneousXdecaysToY_1InnerBoundary/combined_000001/optimization_modifications_analysis.json --cores 1 --use-conda --config fast_mode=true
    input:
        run_result_input_modifications
    output:
        "{df}/{bn}_{cn}/optimization_modifications_analysis.json",
        "{df}/{bn}_{cn}/optimization_modifications_summary.json",
    conda:
        "config/environment.yaml"
    threads: 1
    resources:
        mem_mb=1000,
        runtime= 3
    priority:
        100  # run first
    shell:
        """
        python src/optimizer_analyze_modifications_to_extremum.py \
            --folder {wildcards.df}/{wildcards.bn}_{wildcards.cn} \
        """   

rule summarize_modifications:
    # snakemake -s Snakefile data_private/optimization_spontaneousXdecaysToY_1InnerBoundary/modifications_overall_summary.json --cores 1 --use-conda --config fast_mode=true
    input:
        lambda wildcards: expand(
            "{folder}/optimization_modifications_summary.json",
            folder=sorted(glob.glob(os.path.join(wildcards.df, "combined_*")))
        )
    output:
        "{df}/modifications_overall_summary.json"
    run:
        import json
        flux_increased_folders = []
        for folder, filepath in zip(sim_folders, input):
            with open(filepath) as f:
                data = json.load(f)
            flux_increase = data["flux increase under modification"]
            if flux_increase[0] != 0:
                flux_increased_folders.append(folder)
        summary = {
            "number of combinations where the flux increased": len(flux_increased_folders),
            "indices where the flux increased": flux_increased_folders,
        }
        with open(output[0], "w") as f:
            json.dump(summary, f, indent=4)


if FAST_MODE:
    use rule get_result_from_complete_input as get_result_from_complete_input_within_optimization_modification with:
        output:
            touch("{df}/{bn}_{cn}/optimization_check/modification_{modification}/.result_computed_fast")
else:
    use rule check_reaction_network_info_validity as check_reaction_network_info_validity_within_optimization_modification with:
        group: lambda wildcards: f"sp_{hashlib.md5(trial_path(wildcards).encode()).hexdigest()[:8]}"
        output:
            touch("{df}/{bn}_{cn}/optimization_check/modification_{modification}/.validated_reaction_network_input")

    use rule check_solver_input_validity as check_solver_input_validity_within_optimization_modification with:
        group: lambda wildcards: f"sp_{hashlib.md5(trial_path(wildcards).encode()).hexdigest()[:8]}"
        output:
            touch("{df}/{bn}_{cn}/optimization_check/modification_{modification}/.validated_solver_input")

    use rule check_solver_output_validity as check_solver_output_validity_within_optimization_modification with:
        group: lambda wildcards: f"sp_{hashlib.md5(trial_path(wildcards).encode()).hexdigest()[:8]}"
        output:
            touch("{df}/{bn}_{cn}/optimization_check/modification_{modification}/.validated_solver_output")

    use rule create_system_geometry as create_system_geometry_within_optimization_modification with:
        group: lambda wildcards: f"sp_{hashlib.md5(trial_path(wildcards).encode()).hexdigest()[:8]}"
        output:
            "{df}/{bn}_{cn}/optimization_check/modification_{modification}/.system_geometry.json"

    use rule create_reaction_network as create_reaction_network_within_optimization_modification with:
        group: lambda wildcards: f"sp_{hashlib.md5(trial_path(wildcards).encode()).hexdigest()[:8]}"
        output:
            "{df}/{bn}_{cn}/optimization_check/modification_{modification}/.pickled_reaction_network_without_enzyme_concentration"

    use rule create_initial_guess as create_initial_guess_within_optimization_modification with:
        # snakemake -s Snakefile data_private/case_01/combined_000015/optimization_round_0/trial_0/species_initial_guess.json --cores 1 --use-conda
        output:
            "{df}/{bn}_{cn}/optimization_check/modification_{modification}/species_initial_guess.json"

    use rule cleanup_old_iterations as cleanup_old_iterations_within_optimization_modification with:
        output:
            touch("{df}/{bn}_{cn}/optimization_check/modification_{modification}/.validated_iterations")

    use rule define_enzyme_concentrations as define_enzyme_concentrations_within_optimization_modification with:
        # hashlib to keep it shorter: group: lambda wildcards: f"solver_preparation_{trial_path(wildcards).replace("/", "_")}"
        group: lambda wildcards: f"sp_{hashlib.md5(trial_path(wildcards).encode()).hexdigest()[:8]}"
        output:
            [
                "{df}/{bn}_{cn}/optimization_check/modification_{modification}/.pickled_reaction_network",
                "{df}/{bn}_{cn}/optimization_check/modification_{modification}/enzyme_concentrations.json"
            ]

    use rule solve_boundary_value_problem_with_mesh_adaptation as solve_boundary_value_problem_with_mesh_adaptation_within_optimization_modification with:
        # snakemake -s Snakefile data_private/case_01/combined_000015/optimization_round_0/trial_0/.species_steady_state_concentrations.json --cores 1 --use-conda
        output:
            [
            "{df}/{bn}_{cn}/optimization_check/modification_{modification}/.species_steady_state_concentrations.json",
            "{df}/{bn}_{cn}/optimization_check/modification_{modification}/system_geometry_for_convergence.json",
            "{df}/{bn}_{cn}/optimization_check/modification_{modification}/.expanded_system_mesh_for_convergence.json",
            ]

    use rule study_bvp_solution as study_bvp_solution_within_optimization_modification with:
        output:
            "{df}/{bn}_{cn}/optimization_check/modification_{modification}/fluxes.json"

    use rule plot_boundary_value_problem as plot_boundary_value_problem_within_optimization_modification with:
        # snakemake -s Snakefile data_private/case_01/combined_000015/optimization_round_0/trial_0/.completed_visualization --cores 1 --use-conda
        output:
            touch("{df}/{bn}_{cn}/optimization_check/modification_{modification}/.completed_visualization")
#"""

### To have a specific iteration of the solver plotted, run on terminal ###
# python src/plot_bvp_solution.py data/test_0 --plot_iteration 40
### To have a gif of the iterations of the solver (already before the solver has converged), run on terminal ###
# python src/plot_bvp_solution.py data/test_0


#https://collab.dvb.bayern/spaces/TUMnat/pages/431097554/SLURM+Queuing+system
#https://collab.dvb.bayern/pages/viewpage.action?spaceKey=TUMnat&title=PH+Theory+Cluster

