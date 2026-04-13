import argparse
import sys
import csv
import os
import re
import shutil
import gc
import matplotlib.pyplot as plt
import numpy as np
from contextlib import redirect_stdout
import time
from auxiliary_functions import read_yaml_file
import math
from create_system_mesh import build_system_mesh
from auxiliary_functions_using_standard_library import (
    pickle_load_binary, load_json)
from auxiliary_functions_framework_organization_using_standard_library import (
    get_concentrations_files_within_folder, rename_iteration_files)
from auxiliary_functions import dump_json
from run_bvp_solver import solve_newton
from plot_bvp_solution import plot_steady_state_concentrations
from auxiliary_functions_framework_organization import get_dict_with_correct_key_types_from_json_file
from study_bvp_solution import get_outward_fluxes
from plot_bvp_solver_mesh_adaptation_progress import plot_convergence_progress, make_newton_iterations_gif

def find_latest_solution_of_given_interpolation_iteration(
        concentration_files, interpolation_iteration,
        pattern):
    """
    Returns YY and the full path to the latest file of the form:
    something_XX_somethingother_YY
    (e.g. interpolation_iteration_nr_XX_Newton_iteration_nr_YY_concentrations.json)
    where XX matches the given interpolation_iteration (ignoring leading zeros),
    and YY is maximal.
    """
    matching_files = []
    for path in concentration_files:
        filename = os.path.basename(path)
        match = pattern.fullmatch(filename)
        if not match:
            continue
        interp_str, newton_str = match.groups()
        if int(interp_str) == interpolation_iteration:
            matching_files.append((int(newton_str), path))
    if not matching_files:
        return None, None
    # pick file with largest Newton iteration which is not empty
    # (Snakemake may create the latest file but it might be empty, so we guard against that)
    for newton_iteration, path in sorted(matching_files, key=lambda x: x[0], reverse=True):
        if os.path.getsize(path) > 0:
            return newton_iteration, path

    return None, None
    #latest_newton_iteration, latest_path = max(matching_files, key=lambda x: x[0])
    #return latest_newton_iteration, latest_path

def get_t_n_and_residual_and_runtime_from_progress_log(progress_log_path, iteration_number):
    """
    Reads F_vector_norm (residual norm) from the CSV progress log.
    Uses the row with the largest iteration number that is strictly less than
    iteration_number, as a fallback if the exact iteration is not found.
    """
    if not os.path.exists(progress_log_path):
        raise FileNotFoundError(f"Progress log not found: {progress_log_path}")

    with open(progress_log_path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    # Find the row with the largest iteration strictly less than iteration_number
    candidate_rows = [row for row in rows if int(row["iteration"]) < iteration_number]
    if not candidate_rows:
        raise ValueError(
            f"No row with iteration smaller than {iteration_number} found in {progress_log_path}"
        )

    best_row = max(candidate_rows, key=lambda row: int(row["iteration"]))
    return float(best_row["t_n"]), float(best_row["F_vector_norm"]), float(best_row["runtime"].replace(" seconds", ""))

def get_initial_values(
        iterations_folder,
        interpolation_iteration,
        species_lookup,
        iterations_saving_pattern,
        progress_log_path,
        min_num_iterations_digits,
        system_geometry_dict
    ):
    """
    Returns 
    initial_iteration_number, initial_species_concentrations, initial_tau, initial_residual_norm, initial_runtime
    """
    # Figure out whether another simulation has already been started
    # for the same iterpolation iteration
    concentration_files = get_concentrations_files_within_folder(iterations_folder)
    latest_newton_iteration, path_to_latest_newton_iteration = find_latest_solution_of_given_interpolation_iteration(
        concentration_files, interpolation_iteration, iterations_saving_pattern)
    if latest_newton_iteration is not None:
        previous_solution_species_concentrations_dict_with_strings = load_json(
            path_to_latest_newton_iteration)
        initial_species_concentrations = get_dict_with_correct_key_types_from_json_file(
            previous_solution_species_concentrations_dict_with_strings, species_lookup)
        num_iterations_digits = rename_iteration_files(iterations_folder, min_digits=min_num_iterations_digits)
        initial_t_n, initial_residual_norm, initial_runtime = get_t_n_and_residual_and_runtime_from_progress_log(progress_log_path, latest_newton_iteration)
        print(f"Continuing simulation for from iteration {latest_newton_iteration} within {iterations_folder}.")
        return latest_newton_iteration, initial_species_concentrations, initial_t_n, initial_residual_norm, initial_runtime, num_iterations_digits
    
    # If there are no simulations with the current interpolation
    # iteration, construct the initial concentration out of the
    # previous, if it exists
    if interpolation_iteration > 0:
        previous_solution_species_concentrations_dict_with_strings = load_json(
            os.path.join(iterations_folder, f"interpolation_iteration_nr_{interpolation_iteration-1}_final_concentrations.json"))
        initial_species_concentrations = {
            region_idx : {
                mesh_point_idx : {
                    species_object : 0.5 * (
                        previous_solution_species_concentrations_dict_with_strings[region_idx][mesh_point_idx//2][species_name]
                        + previous_solution_species_concentrations_dict_with_strings[region_idx][(mesh_point_idx+1)//2][species_name]
                    )
                    for species_name, species_object in species_lookup.items()}
                for mesh_point_idx in range(system_geometry_dict["geometry_config"]["num_mesh_points_in_regions"][region_idx])}
            for region_idx in range(system_geometry_dict["geometry_config"]["num_regions"])
        }
        print(f"Continuing from final concentrations from iteration {interpolation_iteration-1}")
        return 0, initial_species_concentrations, None, None, None, min_num_iterations_digits
    
    return None, None, None, None, None, min_num_iterations_digits

def map_radius_to_broader_and_finer_keys_for_region(broader_dict, finer_dict):
    """
    For each shared value between broader_dict and finer_dict, returns a dict
    mapping value -> (broader_key, finer_key) tuple.
    """
    finer_values_to_keys = {v: k for k, v in finer_dict.items()}
    result = {}
    for broader_key, broader_value in broader_dict.items():
        finer_key = finer_values_to_keys.get(broader_value)
        if finer_key is not None:
            result[broader_value] = (broader_key, finer_key)
    return result

def get_max_nested(d):
    """Returns the maximum value within a nested dict"""
    values = []
    for v in d.values():
        if isinstance(v, dict):
            values.append(get_max_nested(v))
        elif v is not None:
            values.append(v)
    return max(values)

def get_solutions_convergence(
        broad_radii, fine_radii,
        broad_species_concentrations,
        fine_species_concentrations,
        species_lookup,
        max_relative_species_concentrations_difference,
        min_relative_concentration_difference_considered_relevant
    ):
    """ Returns True if the solutions computed with the
    different number of mesh points are close enough,
    according to max_relative_species_concentrations_change.

    The reason for introducing an absolute tolerance is that,
    if the values are very close to 0, increasing the refinement
    will change the values a lot, which means that the relative
    deviation will be pretty much 1.0
    """
    absolute_tolerance = min_relative_concentration_difference_considered_relevant * max(
        get_max_nested(fine_species_concentrations),
        get_max_nested(broad_species_concentrations)
    )
    convergence = True
    convergence_values = {}
    for region in range(len(broad_radii)):
        node_mapping = map_radius_to_broader_and_finer_keys_for_region(
            broad_radii[region], fine_radii[region]
        )
        for node_tuple in node_mapping.values():
            for species in species_lookup.keys():
                concentration_broad_system = broad_species_concentrations[region][node_tuple[0]][species]
                concentration_fine_system = fine_species_concentrations[region][node_tuple[1]][species]
                abs_diff = abs(concentration_broad_system - concentration_fine_system)
                magnitude = max(abs(concentration_broad_system), abs(concentration_fine_system))
                if magnitude < absolute_tolerance:
                    # Values are negligibly small — treat as converged
                    rel_deviation = None
                else:
                    rel_deviation = 1 - (min(concentration_broad_system, concentration_fine_system) / magnitude)
                    if rel_deviation > max_relative_species_concentrations_difference:
                        convergence = False

                # Build nested dict without overwriting existing keys
                if region not in convergence_values:
                    convergence_values[region] = {}
                if node_tuple[0] not in convergence_values[region]:
                    convergence_values[region][node_tuple[0]] = {}
                convergence_values[region][node_tuple[0]][species] = rel_deviation
    return convergence, convergence_values






def make_iteration_filename(interpolation_iteration, newton_iteration, min_digits):
    """
    Returns a filename like:
    interpolation_iteration_nr_3_Newton_iteration_nr_007_concentrations.json
    """
    return (
        f"interpolation_iteration_nr_{interpolation_iteration}"
        f"_Newton_iteration_nr_{newton_iteration:0{min_digits}d}"
    )

def dicts_are_similar(dict1, dict2,
    threshold):
    for key in dict1:
        v1, v2 = dict1[key], dict2[key]
        if min(v1, v2) / max(v1, v2) < threshold:
            return False
    
    return True

def save_last_interpolation_iteration_files(
    folder_to_solve, fine_system_mesh_dict, fine_species_concentrations,
    iteration_data_path, interpolation_iteration
):
        dump_json(
            folder_to_solve,
            f".expanded_system_mesh_for_convergence",
            fine_system_mesh_dict
        )
        dump_json(
            folder_to_solve,
            f".species_steady_state_concentrations",
            fine_species_concentrations
        )
        shutil.copy(
            os.path.join(iteration_data_path,
                            f".system_geometry_interpolating_{interpolation_iteration}_times.json"),
            os.path.join(folder_to_solve, "system_geometry_for_convergence.json")
        )
        shutil.copy(
            os.path.join(iteration_data_path,
                            f".system_geometry_interpolating_{interpolation_iteration}_times.json"),
            os.path.join(folder_to_solve, "system_geometry_for_convergence.json")
        )

if __name__ == "__main__":
    # Parse arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--folder", type=str)
    # int(float()) used to be able to pass scientific notation (these are originally read as a string)
    parser.add_argument("--max_num_Newton_iterations", type=lambda x: int(float(x)), help="Maximum Newton iterations") 
    parser.add_argument("--max_num_interpolation_times", type=int)
    parser.add_argument("--max_relative_species_concentrations_difference", type=lambda x: float(x)) 
    parser.add_argument("--max_relative_flux_difference", type=lambda x: float(x))
    parser.add_argument("--min_relative_concentration_difference_considered_relevant", type=lambda x: float(x)) 

    args = parser.parse_args()

    # Load all the passed information, create folder for background info
    FOLDER_TO_SOLVE = args.folder

    # Do not run if the simulation has been pruned. Directly create dummy output file and exit.
    if os.path.isfile(os.path.join(FOLDER_TO_SOLVE, "pruned.json")):
        dump_json(FOLDER_TO_SOLVE, ".species_steady_state_concentrations", {"pruned": True})
        dump_json(FOLDER_TO_SOLVE, "system_geometry_for_convergence", {"pruned": True})
        dump_json(FOLDER_TO_SOLVE, ".expanded_system_mesh_for_convergence", {"pruned": True})
        sys.exit(0)

    ITERATION_DATA_PATH = os.path.join(FOLDER_TO_SOLVE, "solver_iteration_data")
    os.makedirs(ITERATION_DATA_PATH, exist_ok=True)
    max_relative_flux_difference = args.max_relative_flux_difference
    min_relative_concentration_difference_considered_relevant = args.min_relative_concentration_difference_considered_relevant
    print(f"Solving the problem defined in {FOLDER_TO_SOLVE}.")

    SOLVER_INPUT_INFO = read_yaml_file(os.path.join(FOLDER_TO_SOLVE, "parameters_solver_input.yaml"))
    SOLVER_OUTPUT_INFO = read_yaml_file(os.path.join(FOLDER_TO_SOLVE, "parameters_solver_output.yaml"))
    PARAMETER_VALUE_CONDITIONS = read_yaml_file(os.path.join(FOLDER_TO_SOLVE, "parameters_value_conditions.yaml")) 
    MAX_NUM_NEWTON_ITERATIONS = args.max_num_Newton_iterations
    MAX_NUM_INTERPOLATION_TIMES = args.max_num_interpolation_times
    MAX_RELATIVE_SPECIES_CONCENTRATIONS_DIFFERENCE = args.max_relative_species_concentrations_difference
    min_num_iterations_digits = int(math.log10(MAX_NUM_NEWTON_ITERATIONS)+1)
    
    # Load common reaction information
    REACTION_NETWORK = pickle_load_binary(os.path.join(FOLDER_TO_SOLVE, ".pickled_reaction_network"))
    SPECIES_LOOKUP = {sp.name: sp for sp in REACTION_NETWORK.species}
    SYSTEM_GEOMETRY_DICT = load_json(os.path.join(FOLDER_TO_SOLVE, ".system_geometry.json"))
    
    if (SOLVER_OUTPUT_INFO["output_options"]["create_gif_with_saved_data"] is True
        and SOLVER_OUTPUT_INFO["variables_to_save"]["save_concentrations"] is False):
        raise ValueError("Cannot make the gif if the concentrations are not saved.")

    ITERATIONS_SAVING_PATTERN = re.compile(
        r"interpolation_iteration_nr_(\d+)_Newton_iteration_nr_(\d+)_concentrations\.json"
    )

    for interpolation_iteration in range(MAX_NUM_INTERPOLATION_TIMES):
        iteration_final_file = os.path.join(
            ITERATION_DATA_PATH,
            f"interpolation_iteration_nr_{interpolation_iteration}_final_concentrations.json")
        
        if not os.path.isfile(iteration_final_file):
            print(f"Running simulation for interpolation iteration number {interpolation_iteration}.")
            progress_log_path = os.path.join(
                FOLDER_TO_SOLVE,
                f".progress_log_interpolating_{interpolation_iteration}_times.csv"
            )
            # Create system and mesh geometry
            system_geometry_dict, expanded_system_mesh_dict = build_system_mesh(
                SYSTEM_GEOMETRY_DICT, 
                REACTION_NETWORK,
                interpolation_iteration
            )
            # Save system and mesh geometry
            dump_json(
                ITERATION_DATA_PATH,
                f".system_geometry_interpolating_{interpolation_iteration}_times",
                system_geometry_dict
            )
            dump_json(
                ITERATION_DATA_PATH,
                f".expanded_system_mesh_interpolating_{interpolation_iteration}_times",
                expanded_system_mesh_dict
            )

            # Get initial conditions, if these exist with these same parameters from another simulation
            # with the same interpolation iteration or earlier
            initial_iteration_number, initial_species_concentrations, initial_t_n, initial_residual_norm, initial_runtime, num_iterations_digits  = get_initial_values(
                ITERATION_DATA_PATH,
                interpolation_iteration,
                SPECIES_LOOKUP,
                ITERATIONS_SAVING_PATTERN,
                progress_log_path,
                min_num_iterations_digits,
                system_geometry_dict
            )
            # In case no other simulation is found, create initial values
            if initial_iteration_number is None:
                #previous_solution = None
                initial_iteration_number = 0
                #max_external_concentration = max([species.external_concentration for species in REACTION_NETWORK.species])
                #initial_species_concentrations = {
                #    region_idx : {
                #        mesh_point_idx : {
                #            species : species.external_concentration + max_external_concentration*0.1 # max_external_concentration #* RADII[region_idx][mesh_point_idx] / RADII[NUM_REGIONS-1][NUM_MESH_POINTS_IN_REGIONS[NUM_REGIONS-1]-1]
                #            for species in REACTION_NETWORK.species}
                #        for mesh_point_idx in range(system_geometry_dict["geometry_config"]["num_mesh_points_in_regions"][region_idx])}
                #    for region_idx in range(system_geometry_dict["geometry_config"]["num_regions"])
                #}
                initial_species_concentrations = load_json(os.path.join(FOLDER_TO_SOLVE, "species_initial_guess.json"))
                # make the concentrations file have the keys in the correct format
                initial_species_concentrations = get_dict_with_correct_key_types_from_json_file(
                    initial_species_concentrations, SPECIES_LOOKUP)
            if initial_t_n is None:
                initial_t_n = 1
            if initial_residual_norm is None:
                initial_residual_norm = np.inf
            if initial_runtime is None:
                initial_runtime = 0

            filename_for_newton_function = lambda newton_iteration: make_iteration_filename(
                interpolation_iteration, newton_iteration, num_iterations_digits)

            with open(
                os.path.join(FOLDER_TO_SOLVE, f".newton_solver_interpolating_{interpolation_iteration}_times.log"), "a") as f, redirect_stdout(f):
                print(f"Starting solver from iteration number {initial_iteration_number} \n")
                # Run the solver
                start_time = time.time()
                species_concentrations_final = solve_newton(
                    # timing
                    simulation_start_time=start_time,
                    # system information
                    reaction_network = REACTION_NETWORK,
                    system_geometry_dict=system_geometry_dict,
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
                    iteration_data_path = ITERATION_DATA_PATH,
                    progress_log_path=progress_log_path,
                    filename_for_newton_function=filename_for_newton_function,
                    variables_to_save_dictionary = SOLVER_OUTPUT_INFO["variables_to_save"],
                    save_data_every=SOLVER_OUTPUT_INFO["output_options"]["save_data_every"],
                    log_progress_every = SOLVER_OUTPUT_INFO["output_options"]["log_progress_every"],
                    plot_iteration_data_during_simulation = SOLVER_OUTPUT_INFO["output_options"]["plot_iteration_data_during_simulation"]
                )
                end_time = time.time()
            
            dump_json(ITERATION_DATA_PATH,
                    f"interpolation_iteration_nr_{interpolation_iteration}_final_concentrations",
                    species_concentrations_final)
            # draw final concentration of this round
            fig, _ = plot_steady_state_concentrations(
                reaction_network=REACTION_NETWORK,
                num_regions=system_geometry_dict["geometry_config"]["num_regions"],
                num_mesh_points_in_regions=system_geometry_dict["geometry_config"]["num_mesh_points_in_regions"],
                radii=expanded_system_mesh_dict["radii"],
                membrane_radii=system_geometry_dict["geometry_config"]["membrane_radii"],
                output_file_name=os.path.join(ITERATION_DATA_PATH, f"interpolation_iteration_nr_{interpolation_iteration}_final_concentrations.png"),
                species_concentrations_to_plot=species_concentrations_final,
                system_geometry_dict=system_geometry_dict["geometry_config"],
            )
            plt.close(fig) 
            del species_concentrations_final
            gc.collect()

        if interpolation_iteration>0:
            print(f"Comparing results from interpolation iteration numbers {interpolation_iteration-1} and {interpolation_iteration}.")
            # Get maps from radii to node
            broad_system_mesh_dict = load_json(
                os.path.join(ITERATION_DATA_PATH, f".expanded_system_mesh_interpolating_{interpolation_iteration-1}_times.json"))
            fine_system_mesh_dict = load_json(
                os.path.join(ITERATION_DATA_PATH, f".expanded_system_mesh_interpolating_{interpolation_iteration}_times.json"))
            broad_radii = broad_system_mesh_dict["radii"]
            fine_radii = fine_system_mesh_dict["radii"]
            # Load concentrations map from node to concentration
            broad_species_concentrations = load_json(
                os.path.join(
                    ITERATION_DATA_PATH,
                    f"interpolation_iteration_nr_{interpolation_iteration-1}_final_concentrations.json"
                )
            )
            fine_species_concentrations = load_json(
                os.path.join(
                    ITERATION_DATA_PATH,
                    f"interpolation_iteration_nr_{interpolation_iteration}_final_concentrations.json"
                )
            )
            convergence, convergence_dict = get_solutions_convergence(
                broad_radii, fine_radii,
                broad_species_concentrations,
                fine_species_concentrations,
                SPECIES_LOOKUP,
                MAX_RELATIVE_SPECIES_CONCENTRATIONS_DIFFERENCE,
                min_relative_concentration_difference_considered_relevant
            )
            dump_json(
                ITERATION_DATA_PATH,
                f".convergence_data_between_{interpolation_iteration-1}_and_{interpolation_iteration}_interpolations",
                convergence_dict
            )
            if convergence:
                print("Converged! Saving.")
                open(os.path.join(FOLDER_TO_SOLVE, "concentration_convergence"), "w").close()
                save_last_interpolation_iteration_files(
                    FOLDER_TO_SOLVE,
                    fine_system_mesh_dict,
                    fine_species_concentrations,
                    ITERATION_DATA_PATH,
                    interpolation_iteration
                )
                print("Converged! Saved. Plotting.")
                plot_convergence_progress(FOLDER_TO_SOLVE, REACTION_NETWORK)
                make_newton_iterations_gif(FOLDER_TO_SOLVE, REACTION_NETWORK, SPECIES_LOOKUP)
                print("Plotted.")
                sys.exit()
            else:
                print("Did not converge. Refine the mesh.")
    # In case there has been no convergence:
    # (this can happen if reactions are extremely fast compared to diffusion,
    # such that there are large gradients close to the membranes)
    
    # calculate the fluxes within the last 2 interpolation iterations
    fluxes_list = []
    max_concentration = 0
    for interpolation_iteration in [MAX_NUM_INTERPOLATION_TIMES-2, MAX_NUM_INTERPOLATION_TIMES-1]:
        species_concentrations_dict = load_json(
            os.path.join(ITERATION_DATA_PATH, f"interpolation_iteration_nr_{interpolation_iteration}_final_concentrations.json")
        )
        system_geometry = load_json(
            os.path.join(ITERATION_DATA_PATH, f".system_geometry_interpolating_{interpolation_iteration}_times.json")
        )
        num_regions = system_geometry["geometry_config"]["num_regions"]
        num_mesh_points_in_regions = system_geometry["geometry_config"]["num_mesh_points_in_regions"]
        fluxes = get_outward_fluxes(
            species_concentrations_dict,
            REACTION_NETWORK,
            num_regions,
            num_mesh_points_in_regions
        )
        fluxes_list.append(fluxes)
        max_concentration = max(max_concentration, get_max_nested(species_concentrations_dict))
    
    if dicts_are_similar(
        fluxes_list[0], fluxes_list[1],
        threshold=max_relative_flux_difference,
    ):
        open(os.path.join(FOLDER_TO_SOLVE, "flux_convergence_without_concentration_convergence"), "w").close()
    else:
        open(os.path.join(FOLDER_TO_SOLVE, "no_flux_nor_concentration_convergence"), "w").close()

    # save anyways
    save_last_interpolation_iteration_files(
        FOLDER_TO_SOLVE,
        fine_system_mesh_dict,
        fine_species_concentrations,
        ITERATION_DATA_PATH,
        interpolation_iteration
    )
    
    plot_convergence_progress(FOLDER_TO_SOLVE, REACTION_NETWORK)
    make_newton_iterations_gif(FOLDER_TO_SOLVE, REACTION_NETWORK, SPECIES_LOOKUP)
    print("Plotted.")
