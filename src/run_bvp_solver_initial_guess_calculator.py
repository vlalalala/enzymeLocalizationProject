import os
import sys
import argparse
import re
import math
import numpy as np
from contextlib import redirect_stdout
import time
from define_initial_guess_calculator_steps import create_creeping_reaction_folder, calculate_timescales
from auxiliary_functions import read_yaml_file, dump_json
from auxiliary_functions_using_standard_library import load_json, pickle_load_binary
from run_bvp_solver_mesh_adaptation import (
    get_initial_values, make_iteration_filename
    )
from create_system_mesh import build_system_mesh
from run_bvp_solver import solve_newton
import shutil
from pathlib import Path
from auxiliary_functions_framework_organization import get_dict_with_correct_key_types_from_json_file
from plot_bvp_solution import plot_steady_state_concentrations
from plot_bvp_solver_initial_guess_calculator import make_newton_iterations_gif

ITERATIONS_SAVING_PATTERN = re.compile(
    r"interpolation_iteration_nr_(\d+)_Newton_iteration_nr_(\d+)_concentrations\.json"
)
def str2bool(v):
    return v.lower() in ("yes", "true", "1")

def get_non_reaction_concentrations(reaction_network, system_geometry):
    max_external_concentration = max([species.external_concentration for species in reaction_network.species])
    initial_species_concentrations = {
        region_idx : {
            mesh_point_idx : {
                species : species.external_concentration + max_external_concentration*0.01 # in order for the concentrations not to be exactly 0
                for species in reaction_network.species}
            for mesh_point_idx in range(system_geometry["geometry_config"]["num_mesh_points_in_regions"][region_idx])}
        for region_idx in range(system_geometry["geometry_config"]["num_regions"])
    }
    return initial_species_concentrations

if __name__ == "__main__":
    # Parse arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--folder", type=str)
    # int(float()) used to be able to pass scientific notation (these are originally read as a string)
    parser.add_argument("--max_num_Newton_iterations", type=lambda x: int(float(x)), help="Maximum Newton iterations") 
    parser.add_argument("--max_num_creeping_reaction_simulations", type=int)
    parser.add_argument("--override", type=str2bool)

    args = parser.parse_args()
    FOLDER_TO_SOLVE = args.folder
    MAX_NUM_NEWTON_ITERATIONS = args.max_num_Newton_iterations
    min_num_iterations_digits = int(math.log10(MAX_NUM_NEWTON_ITERATIONS)+1)
    MAX_NUM_CREEPING_REACTION_SIMULATIONS = args.max_num_creeping_reaction_simulations
    OVERRIDE = args.override

    SOLVER_INPUT_INFO = read_yaml_file(os.path.join(FOLDER_TO_SOLVE, "parameters_solver_input.yaml"))
    WANTED_MAXIMUM_INITIAL_RATIO = SOLVER_INPUT_INFO["initial_concentration_guess_parameters"]["wanted_maximum_initial_ratio_diffusion_to_fastest_reaction_timescales"]
    WANTED_RATIO_GROWTH_FACTOR = SOLVER_INPUT_INFO["initial_concentration_guess_parameters"]["gamma_inc"]

    SYSTEM_GEOMETRY_DICT = load_json(os.path.join(FOLDER_TO_SOLVE, ".system_geometry.json"))
    SOLVER_OUTPUT_INFO = read_yaml_file(os.path.join(FOLDER_TO_SOLVE, "parameters_solver_output.yaml"))
    PARAMETER_VALUE_CONDITIONS = read_yaml_file(os.path.join(FOLDER_TO_SOLVE, "parameters_value_conditions.yaml")) 
    
    
    if os.path.isfile(os.path.join(FOLDER_TO_SOLVE, "pruned.json")):
        dump_json(FOLDER_TO_SOLVE, "species_initial_guess", {"pruned": True})
        sys.exit(0)

    ORIGINAL_REACTION_NETWORK = pickle_load_binary(os.path.join(FOLDER_TO_SOLVE, ".pickled_reaction_network"))
    SPECIES_LOOKUP = {sp.name: sp for sp in ORIGINAL_REACTION_NETWORK.species}
    CREEPING_REACTION_SIMULATIONS_FOLDER = os.path.join(FOLDER_TO_SOLVE, "creeping_reaction_simulations")

    # This try statement is here for the case in which a modification has been forced
    # and this is the first point at which it is checked whether the geometry is actually
    # solvable. Creates the prune file automatically
    try:
        expanded_system_geometry_dict, expanded_system_mesh_dict = build_system_mesh(
                    SYSTEM_GEOMETRY_DICT, 
                    ORIGINAL_REACTION_NETWORK,
                    0
                )
    except ValueError as e:
        if ("less than 3 points" in str(e)
            and os.path.isfile(os.path.join(FOLDER_TO_SOLVE, "info_on_modification.txt"))):
            pruned_file_dict = {"prune": True, "reason": f"The modification means that a region has less than 3 points ({str(e)})"}
            dump_json(FOLDER_TO_SOLVE, "pruned", {"pruned": pruned_file_dict})
            dump_json(FOLDER_TO_SOLVE, "species_initial_guess", {"pruned": True})
            sys.exit(0)
        else:
            raise  # re-raise if it's a different ValueError
       
    # Inform about the original timescales
    timescales_log_path = os.path.join(FOLDER_TO_SOLVE, f"reaction_timescales.log")
    with open(timescales_log_path, "a") as f, redirect_stdout(f): 
        for species in ORIGINAL_REACTION_NETWORK.species:
            calculate_timescales(species, expanded_system_geometry_dict["geometry_config"]["outer_membrane_radius"])
    
    if OVERRIDE is True:
        initial_species_concentrations = get_non_reaction_concentrations(ORIGINAL_REACTION_NETWORK, expanded_system_geometry_dict)
        dump_json(FOLDER_TO_SOLVE, "species_initial_guess", initial_species_concentrations)
        # Inform about the timescales
        
        
        sys.exit()

    # Check that it is possible for the maximum ratio to reach 1
    if WANTED_MAXIMUM_INITIAL_RATIO * WANTED_RATIO_GROWTH_FACTOR ** (MAX_NUM_CREEPING_REACTION_SIMULATIONS-1) < 1:
        raise ValueError("The max number of creeping reaction simulations must be larger.")

    # Step 1: Create or open the log of maximum ratio of diffusion to reaction timescales and success log
    log_file_maximum_ratio = os.path.join(
        CREEPING_REACTION_SIMULATIONS_FOLDER,
        "wanted_maximum_ratio_diffusion_to_fastest_reaction_timescales.json"                              
    )
    if not os.path.isfile(log_file_maximum_ratio):
        # If the file does not exist, create it, with the initial value as the one passed
        wanted_maximum_ratio_log = {0: [WANTED_MAXIMUM_INITIAL_RATIO, "simulation_unfinished"]}
        dump_json(
            CREEPING_REACTION_SIMULATIONS_FOLDER,
            "wanted_maximum_ratio_diffusion_to_fastest_reaction_timescales",
            wanted_maximum_ratio_log
        )
    else:
        wanted_maximum_ratio_log = load_json(log_file_maximum_ratio)

    # Step 2: start checking/running simulations starting from the largest key
    last_creeping_simulation_started = max(wanted_maximum_ratio_log) # gets max key
    creeping_reaction_simulation_idx_larger_than_1 = None
    for creeping_reaction_simulation_idx in range(
        last_creeping_simulation_started, MAX_NUM_CREEPING_REACTION_SIMULATIONS):
        creeping_reaction_simulation_folder = os.path.join(
            CREEPING_REACTION_SIMULATIONS_FOLDER,
            f"creeping_reaction_simulation_{creeping_reaction_simulation_idx}")
        #final_creeping_simulation_file = os.path.join(
        #    creeping_reaction_simulation_folder, "species_steady_state_concentrations.json")
        # Step 2.1: Create necessary data (the folder is created if it does not already exist internally)
        # modification saves whether the timescale for diffusion is shorter than whan the maximum timescale expected in this simulation index
        modification = create_creeping_reaction_folder(
            FOLDER_TO_SOLVE,
            creeping_reaction_simulation_folder,
            wanted_maximum_ratio_diffusion_to_fastest_reaction_timescales = wanted_maximum_ratio_log[creeping_reaction_simulation_idx][0]
        )

        # Step 2.2. Create iteration data path if it doesn't exist (I'm writing this as a separate step in case some KeyboardInterrupt
        # occurs)
        iteration_data_path = os.path.join(creeping_reaction_simulation_folder, "solver_iteration_data")
        if not os.path.isdir(iteration_data_path):
            os.makedirs(iteration_data_path)

        progress_log_path = os.path.join(
            creeping_reaction_simulation_folder,
            f".progress_log.csv"
        )

        # Step 2.3: Figure out whether we can continue an existing simulation
        initial_iteration_number, initial_species_concentrations, initial_t_n, initial_residual_norm, initial_runtime, num_iterations_digits  = get_initial_values(
            iteration_data_path,
            0,
            SPECIES_LOOKUP,
            ITERATIONS_SAVING_PATTERN,
            progress_log_path,
            min_num_iterations_digits,
            SYSTEM_GEOMETRY_DICT
        )
        # Find the latest simulation that was successful in converging
        latest_successful_simulation_index = max(
                (k for k, v in wanted_maximum_ratio_log.items() if v[1] == "simulation_early_exit"),
                default=None
            )
        # Step 2.4: If we cannot start from an existing simulation with the same parameters,
        # define the initial concentration for the simulation
        if initial_iteration_number is None:
            initial_iteration_number = 0
            # Step 2.4.1: Define the initial concentrations
            # If it is the first simulation with modified kinetic rates or no other simulation has converged,
            # define the initial concentration as the external concentration for each species
            if creeping_reaction_simulation_idx == 0 or latest_successful_simulation_index == None:
                initial_species_concentrations = get_non_reaction_concentrations(ORIGINAL_REACTION_NETWORK, expanded_system_geometry_dict)
            # Else: there exists some simulation that has converged. Use the latest one as a starting point.
            else:   
                latest_successful_creeping_reaction_simulation_folder = os.path.join(
                    CREEPING_REACTION_SIMULATIONS_FOLDER,
                    f"creeping_reaction_simulation_{latest_successful_simulation_index}")
                latest_successful_steady_state_concentrations_file = os.path.join(latest_successful_creeping_reaction_simulation_folder,
                    ".species_steady_state_concentrations.json")
                initial_species_concentrations = load_json(latest_successful_steady_state_concentrations_file)
                # make the concentrations file have the keys in the correct format
                initial_species_concentrations = get_dict_with_correct_key_types_from_json_file(
                    initial_species_concentrations, SPECIES_LOOKUP)

        if initial_t_n is None:
            initial_t_n = 1
        if initial_residual_norm is None:
            initial_residual_norm = np.inf
        if initial_runtime is None:
            initial_runtime = 0
        
        # Step 2.5: Run the simulation
        filename_for_newton_function = lambda newton_iteration: make_iteration_filename(
            0, newton_iteration, num_iterations_digits)
        if modification:
            # Modified reaction network
            REACTION_NETWORK = pickle_load_binary(
                os.path.join(creeping_reaction_simulation_folder, ".pickled_reaction_network"))
        log_path = os.path.join(CREEPING_REACTION_SIMULATIONS_FOLDER, f".newton_solver_creeping_simulation_{creeping_reaction_simulation_idx}.log")
        with open(log_path, "a") as f, redirect_stdout(f):
            if modification:
                print(f"Starting solver from iteration number {initial_iteration_number} \n")
                # Run the solver
                start_time = time.time()
                species_concentrations_final = solve_newton(
                    # timing
                    simulation_start_time=start_time,
                    # system information
                    reaction_network = REACTION_NETWORK,
                    system_geometry_dict=SYSTEM_GEOMETRY_DICT,
                    expanded_system_mesh_dict=expanded_system_mesh_dict,
                    adaptive_step_parameters=SOLVER_INPUT_INFO["adaptive_step_parameters"],
                    max_num_newton_iterations=MAX_NUM_NEWTON_ITERATIONS,
                    # simulation initial values
                    initial_iteration_number=initial_iteration_number,
                    initial_t_n = initial_t_n,
                    initial_species_concentrations=initial_species_concentrations,
                    initial_residual_norm=initial_residual_norm,
                    initial_runtime = initial_runtime,
                    # saving
                    iteration_data_path = iteration_data_path,
                    progress_log_path=progress_log_path,
                    filename_for_newton_function=filename_for_newton_function,
                    variables_to_save_dictionary = SOLVER_OUTPUT_INFO["variables_to_save"],
                    save_data_every=SOLVER_OUTPUT_INFO["output_options"]["save_data_every"],
                    log_progress_every = SOLVER_OUTPUT_INFO["output_options"]["log_progress_every"],
                    plot_iteration_data_during_simulation = SOLVER_OUTPUT_INFO["output_options"]["plot_iteration_data_during_simulation"]
                )
                end_time = time.time()
            else:
                print(f"Simulation did not have to run, since diffusion is fast enough.")
                species_concentrations_final = initial_species_concentrations

        
        
        # Step 2.6: update state of simulation in log dictionary
        with open(log_path, "r") as f:
            lines = f.readlines()
            last_line = lines[-1].strip() if lines else ""
        print(last_line)
        # Check whether an early convergence was logged
        if "numerical limit" in last_line or "negative" in last_line: ########################### not sure how to separate cases... 
            # Early convergence
            wanted_maximum_ratio_log[creeping_reaction_simulation_idx][1] = "simulation_early_exit"
        elif "diffusion is fast enough" in last_line:
            wanted_maximum_ratio_log[creeping_reaction_simulation_idx][1] = "simulation_not_required"
        else:
            wanted_maximum_ratio_log[creeping_reaction_simulation_idx][1] = "simulation_non_early_exit"
        
        # Step 2.7: save concentrations
        dump_json(creeping_reaction_simulation_folder,
                ".species_steady_state_concentrations",
                species_concentrations_final)
        plot_steady_state_concentrations(
            ORIGINAL_REACTION_NETWORK,
            num_regions=SYSTEM_GEOMETRY_DICT["geometry_config"]["num_regions"],
            num_mesh_points_in_regions=SYSTEM_GEOMETRY_DICT["geometry_config"]["num_mesh_points_in_regions"],
            radii=expanded_system_mesh_dict["radii"],
            membrane_radii=SYSTEM_GEOMETRY_DICT["geometry_config"]["membrane_radii"],
            output_file_name=os.path.join(CREEPING_REACTION_SIMULATIONS_FOLDER, f"final_concentration_creeping_reaction_idx_{creeping_reaction_simulation_idx}.png"),
            species_concentrations_to_plot=species_concentrations_final,
            system_geometry_dict=SYSTEM_GEOMETRY_DICT["geometry_config"],
        )
        make_newton_iterations_gif(FOLDER_TO_SOLVE, creeping_reaction_simulation_folder)
        
        # Step 2.7: Define the maximum ratio wanted for the next simulation
        # If the simulation run in this iteration of the for loop exited early, increase the growth factor
        if (wanted_maximum_ratio_log[creeping_reaction_simulation_idx][1] == "simulation_early_exit"
            or wanted_maximum_ratio_log[creeping_reaction_simulation_idx][1] == "simulation_not_required"):
            print("Increasing the reaction rates!")
            new_factor = wanted_maximum_ratio_log[creeping_reaction_simulation_idx][0] * WANTED_RATIO_GROWTH_FACTOR
        else:
            print("Decreasing the reaction rates!")
            if latest_successful_simulation_index is not None:
                latest_successful_ratio = wanted_maximum_ratio_log[latest_successful_simulation_index][0]
                new_factor = abs((wanted_maximum_ratio_log[creeping_reaction_simulation_idx][1] - latest_successful_ratio)/2)
            else:
                # If there hasn't been any successful simulation, take the initial one and half it
                # as many times as needed
                new_factor = WANTED_MAXIMUM_INITIAL_RATIO / 2**(creeping_reaction_simulation_idx+1)

        if new_factor > 1:
            creeping_reaction_simulation_idx_larger_than_1 = creeping_reaction_simulation_idx
            break
        
        # Step 2.8: update the dictionary in temporary and permanent storage
        wanted_maximum_ratio_log.update({
                creeping_reaction_simulation_idx+1: [new_factor, "simulation_unfinished"]})
        dump_json(
            CREEPING_REACTION_SIMULATIONS_FOLDER,
            "wanted_maximum_ratio_diffusion_to_fastest_reaction_timescales",
            wanted_maximum_ratio_log
        )
    
    # Step 3: If new_factor larger than 1: save initial guess as a concentration file
    # If for loop broke: raise ValueError: maybe not enough iterations 
    # Build in a check: if all simulations work, the number gamma should be enough to have the factor reach over 1. 
    if creeping_reaction_simulation_idx_larger_than_1 is not None:
        solution_path = Path(CREEPING_REACTION_SIMULATIONS_FOLDER) / f"creeping_reaction_simulation_{creeping_reaction_simulation_idx}" / ".species_steady_state_concentrations.json"
        destination_path = Path(FOLDER_TO_SOLVE) / "species_initial_guess.json"
        if not os.path.isfile(solution_path):
            raise ValueError(f"Path {solution_path} does not exist.")
        shutil.copy(src = solution_path, dst=destination_path)
        if not os.path.isfile(destination_path):
            raise ValueError(f"Path {destination_path} does not exist. Shutil failed silently.")
    else:
        raise ValueError("Could not reach a factor larger than 1. Maybe the max number of Newton iterations is too low or the max number of creeping simulation steps is too low.")
        
        
        

"""
# If the folder does not exist, create it
if not os.path.isdir(creeping_reaction_simulation_folder):
    create_creeping_reaction_folder(
        FOLDER_TO_SOLVE,
        creeping_reaction_simulation_folder,
        wanted_maximum_ratio_diffusion_to_fastest_reaction_timescales=?
    )
    if creeping_reaction_simulation_idx == 0:
        wanted_maximum_ratio_diffusion_to_fastest_reaction_timescales = WANTED_MAXIMUM_INITIAL_RATIO
    else:
        previous_creeping_reaction_simulation_folder = os.path.join(
            CREEPING_REACTION_SIMULATIONS_FOLDER,
            f"creeping_reaction_simulation_{creeping_reaction_simulation_idx-1}")
        previous_steady_state_concentrations_file = os.path.join(previous_creeping_reaction_simulation_folder,
            ".species_steady_state_concentrations.json")
        previous_wanted_maximum_ratio = wanted_maximum_ratio[creeping_reaction_simulation_idx-1]
        if not os.path.isfile(previous_steady_state_concentrations_file):
            # If .species_steady_state_concentrations.json for the previous simulation does not exist,
            # then that means that the previous simulation did not converge and so 
            # the wanted_maximum_ratio has to be reduced
            if creeping_reaction_simulation_idx == 1:
                wanted_maximum_ratio_diffusion_to_fastest_reaction_timescales = wanted_maximum_ratio[creeping_reaction_simulation_idx-1] * 

# Step 1: Check how many creeping reaction simulation folders already exist
folders = [os.path.join(FOLDER_TO_SOLVE, path)
        for path in os.listdir(FOLDER_TO_SOLVE)
        if path.startswith("creeping_reaction")
]

# If this number is already the max num creeping reaction simulations
# Create initial guess concentrations file with concentrations equal
# to external concentration
if n_previous_creeping_reaction_folders >= MAX_NUM_CREEPING_REACTION_SIMULATIONS:
    raise ValueError("The maximum number of creeping reaction simulations has already been met.")

# Step 3:

# Step 2: Create the creeping reaction folder

initial_species_concentrations = {
    region_idx : {
        mesh_point_idx : {
            species : species.external_concentration
            for species in REACTION_NETWORK.species}
        for mesh_point_idx in range(system_geometry_dict["geometry_config"]["num_mesh_points_in_regions"][region_idx])}
    for region_idx in range(system_geometry_dict["geometry_config"]["num_regions"])
}
"""