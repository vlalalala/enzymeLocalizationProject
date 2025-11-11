#%%
import sys
import os
import copy
from typing import Dict, Any
import time
import glob
import re
from contextlib import redirect_stdout
import argparse
from tqdm import tqdm
import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_matrix
import matplotlib.pyplot as plt
from auxiliary_functions_using_standard_library import (rename_iteration_files,
    nested_max, all_non_negative, format_sci, pickle_load_binary,
    load_json)
from auxiliary_functions import dump_json
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme
from plot_boundary_value_problem import plot_steady_state_concentrations, make_newton_iterations_gif
from auxiliary_functions_framework_organization import (
    get_species_concentrations_from_json_file, save_newton_iteration_data)
from create_system_mesh import (build_point_ids_dict, build_reverse_point_ids_dict, 
    build_radii_dict, build_point_infos_dict, build_point_neighbor_dict)

def calculate_reaction_term(current_species_concentrations, region, n, species):
    """ Gives the reaction term for F_i.
    """
    reaction_term = 0
    for reaction in species.as_reactant_in + species.as_product_in:
        if isinstance(reaction, SpontaneousReaction):
            term = reaction.k * current_species_concentrations[region][n][reaction.start_species]
        else:
            term = reaction.k_cat * ENZYMES_CONCENTRATIONS[region][reaction.enzyme] * current_species_concentrations[region][n][reaction.start_species] / (reaction.k_M + current_species_concentrations[region][n][reaction.start_species])
        if reaction in species.as_reactant_in: # if acts as reactant, diminishes
            term *= -1
        reaction_term += term
    return reaction_term

def calculate_reaction_partial_derivative(current_species_concentrations, reaction_to_check, partial_derivative_species, region, n):
    """ Gives the partial derivative of a reaction to a concentration of a species
    that is involved in the reaction.
    """
    if isinstance(reaction_to_check, SpontaneousReaction):
        derivative = reaction_to_check.k
    elif isinstance(reaction_to_check, EnzymaticReaction):
        derivative = reaction_to_check.k_cat * ENZYMES_CONCENTRATIONS[region][reaction_to_check.enzyme] * reaction_to_check.k_M / ( reaction_to_check.k_M + current_species_concentrations[region][n][partial_derivative_species])
    if partial_derivative_species == reaction_to_check.start_species:
        derivative *= -1
    return derivative

def get_pore_density_occupation_information(current_species_concentrations, info_minus_side, info_plus_side):
    """ Gives necessary information about the density of occupied pores by each species, when given the current
    concentrations of species dictionary. It also requires 
    info_minus_side and info_plus_side, which are each tuples (region, n).
    In case the info_plus_side is such that the region would correspond to the exterior (and thus the 
    current_species_concentrations would not be able to read a concentration), it gives the exterior concentration
    """
    # concentration_rate_ratio_factor = kon/2koff * ( M- + M+ )
    region_minus, region_minus_last_n = info_minus_side
    region_plus, region_plus_first_n = info_plus_side
    concentration_rate_ratio_factor = {specific_species: None for specific_species in REACTION_NETWORK.species}
    for specific_species in REACTION_NETWORK.species:
        concentration_left = current_species_concentrations[region_minus][region_minus_last_n][specific_species]
        if region_plus == NUM_REGIONS:
            concentration_right = specific_species.external_concentration
        else:
            concentration_right = current_species_concentrations[region_plus][region_plus_first_n][specific_species]
        concentration_rate_ratio_factor[specific_species] = specific_species.k_on/specific_species.k_off * (concentration_left + concentration_right)
    sum_concentration_rate_ratio_factor = sum(concentration_rate_ratio_factor.values())
    occupied_pore_density = {
        specific_species: PORE_DENSITY * concentration_rate_ratio_factor[specific_species] / (1 + sum_concentration_rate_ratio_factor)
        for specific_species in REACTION_NETWORK.species
    }
    total_occupied_pore_density = sum(occupied_pore_density.values())
    return concentration_rate_ratio_factor, sum_concentration_rate_ratio_factor, occupied_pore_density, total_occupied_pore_density

def define_newton_residual_and_optionally_jacobian(current_species_concentrations, fill_jacobian = True):
    """Defines the residual vector F and the jacobian matrix J (not sparse) 
    (the latter only if fill_jacobian is set to True (default)).
    Returns either F, _ or F, J.    
    """
    F = np.zeros(NUM_POINTS)
    if fill_jacobian:
        J = lil_matrix((NUM_POINTS, NUM_POINTS))# np.zeros((NUM_POINTS, NUM_POINTS)) 
    for i in range(NUM_POINTS):
        (region, n, species) = REVERSE_POINT_IDS[i]
        r = RADII[region][n]
        diff = species.diffusion_constant
        point_type = POINT_INFOS[region][n]
        if MEMBRANE_TYPE == "enzymatic":
            k_on = species.k_on
            k_off = species.k_off
        # CONSTRUCT F_i
        # FOR EACH POINT WITHIN THE BULK
        if point_type == "i":
            (_, left_n), (_, center_n), (_, right_n) = NEIGHBORS[(region, n)]
            c_left = current_species_concentrations[region][left_n][species]
            c_center = current_species_concentrations[region][center_n][species]
            c_right = current_species_concentrations[region][right_n][species]
            diffusion_term = diff * (1/ DELTA_R**2 * (c_right - 2* c_center + c_left) + 1 /(DELTA_R*r) * (c_right - c_left))
            reaction_term = calculate_reaction_term(current_species_concentrations, region, center_n, species)
            F[i] = diffusion_term + reaction_term
            # FILL IN J_ij
            if not fill_jacobian:
                continue
            for j in range(NUM_POINTS):
                (j_region, j_n, j_species) = REVERSE_POINT_IDS[j]
                # Contributions from diffusion
                if j_region == region and j_n == n and j_species == species: # j == i, basically
                    J[i,j] += diff * (1/DELTA_R**2 * (-2))
                elif j_region==region and j==right_n and j_species == species: # same species, right or left
                    J[i,j] += diff * (1/DELTA_R**2 + 1/(DELTA_R*r))
                elif j_region==region and j==left_n and j_species == species: # same species, right or left
                    J[i,j] += diff * (1/DELTA_R**2 - 1/(DELTA_R*r))
                # Contributions from reactions
                if j_region == region and j_n == center_n: # if on the same place but not necessarily the same species
                    for reaction in species.as_reactant_in + species.as_product_in:
                        if j_species in [reaction.start_species, reaction.end_species]:
                            J[i,j] += calculate_reaction_partial_derivative(current_species_concentrations, reaction, j_species, region, center_n)
        elif point_type == "l":
            if region==0: # deal with r=0 point, no membrane
                (_, r0_n), (_, r0_neighbor_n) = NEIGHBORS[(region, n)]
                c_r0 = current_species_concentrations[region][r0_n][species]
                c_r0_neighbor = current_species_concentrations[region][r0_neighbor_n][species]
                diffusion_term = 3 * diff / DELTA_R**2 * 2 * (c_r0_neighbor - c_r0)
                reaction_term = calculate_reaction_term(current_species_concentrations, region, r0_n, species)
                F[i] = diffusion_term + reaction_term
                if not fill_jacobian:
                    continue
                for j in range(NUM_POINTS):
                    (j_region, j_n, j_species) = REVERSE_POINT_IDS[j]
                    # Contributions from diffusion
                    if j_region == region and j_n == 0 and j_species == species: # j == i, basically
                        J[i,j] += -3 * diff / DELTA_R**2 * 2
                    elif j_region == region and j_n == 1 and j_species == species: # partial derivative to the one on the right
                        J[i,j] += 3 * diff / DELTA_R**2 * 2
                    # Contributions from reactions
                    if j_region == region and j_n == n: # if on the same place but not necessarily the same species
                        for reaction in species.as_reactant_in + species.as_product_in:
                            if j_species in [reaction.start_species, reaction.end_species]:
                                J[i,j] += calculate_reaction_partial_derivative(current_species_concentrations, reaction, j_species, region, n)
            else: # deal with left-most point within region (except r=0)
                (prev_region, prev_region_last_n), (_, _), (_, _) = NEIGHBORS[(region, n)]
                c_prev_region_last = current_species_concentrations[prev_region][prev_region_last_n][species]
                c_region_first = current_species_concentrations[region][0][species]
                c_region_second = current_species_concentrations[region][1][species]
                if MEMBRANE_TYPE == "permeability":
                    F[i] = diff  * (c_region_second - c_region_first) / DELTA_R - species.permeability_constant * (c_region_first - c_prev_region_last)
                    if not fill_jacobian:
                        continue
                    for j in range(NUM_POINTS):
                        (j_region, j_n, j_species) = REVERSE_POINT_IDS[j]
                        # Contributions from diffusion
                        if j_region == region and j_n == n and j_species == species:
                            J[i,j] += -diff/DELTA_R - species.permeability_constant
                        elif j_region == region and j_species == species and j_n == 1:
                            J[i,j] += diff/DELTA_R
                        elif j_region == prev_region and j_species == species and j_n == prev_region_last_n:
                            J[i,j] += -species.permeability_constant
                        # No contributions from reactions (flux considered)
                elif MEMBRANE_TYPE == "enzymatic":
                    # membrane is at the left of the segment, dM_+/dt
                    (concentration_rate_ratio_factor, sum_concentration_rate_ratio_factor,
                        occupied_pore_density, total_occupied_pore_density) = get_pore_density_occupation_information(current_species_concentrations, (prev_region, prev_region_last_n), (region, 0))
                    flux_term = -k_on * (PORE_DENSITY - total_occupied_pore_density) * current_species_concentrations[region][0][species] + k_off * occupied_pore_density[species]
                    F[i] = diff * (c_region_second - c_region_first) / DELTA_R - flux_term
                    if not fill_jacobian:
                        continue
                    for j in range(NUM_POINTS):
                        (j_region, j_n, j_species) = REVERSE_POINT_IDS[j]
                        # Diffusion contribution
                        if j_region == region and j_n == n and j_species == species: # derivative of c_region_first
                            J[i,j] += -diff/DELTA_R
                        elif j_region == region and j_species == species and j_n == 1: # derivative of c_region_second
                            J[i,j] += diff/DELTA_R
                        # Flux contribution
                        if j_region == region and j_n == n and j_species == species: # derivative to concentration on right of flux term
                            derivative_occupied_pore_density = 1 ##########
                            complete_derivative = (
                                -k_on * (PORE_DENSITY - total_occupied_pore_density) # product rule!
                                -k_on * current_species_concentrations[region][0][species] * (-1) #################################
                            )
                            J[i,j] += complete_derivative

                    
        else: # point_type == "r"
            if region == NUM_REGIONS-1: # deal with r=R point
                (_, rR_neighbor_n), (_, rR_n) = NEIGHBORS[(region, n)]
                c_rR_neighbor = current_species_concentrations[region][rR_neighbor_n][species]
                c_rR = current_species_concentrations[region][rR_n][species]
                if MEMBRANE_TYPE == "permeability":
                    F[i] = diff * (c_rR - c_rR_neighbor) / DELTA_R - species.permeability_constant * (species.external_concentration - c_rR)
                    if not fill_jacobian:
                        continue
                    # CONSTRUCT J_ij
                    for j in range(NUM_POINTS):
                        (j_region, j_n, j_species) = REVERSE_POINT_IDS[j]
                        if j_region == region and j_n == n and j_species == species: # basically i=j
                            J[i,j] += diff/DELTA_R + species.permeability_constant
                        elif j_region == region and j_species == species and j_n == rR_neighbor_n:
                            J[i,j] += -diff/DELTA_R
                elif MEMBRANE_TYPE == "enzymatic":
                    raise NotImplementedError("Enzymatic")
            else: # deal with right-most point within region (except r=R)
                (_, _), (_, _), (next_region, _) = NEIGHBORS[(region, n)]
                c_second_to_last = current_species_concentrations[region][n-1][species]
                c_last = current_species_concentrations[region][n][species]
                c_next_region_first = current_species_concentrations[next_region][0][species]
                if MEMBRANE_TYPE == "permeability":
                    F[i] = diff * (c_last - c_second_to_last) / DELTA_R - species.permeability_constant * (c_next_region_first - c_last)
                    if not fill_jacobian:
                        continue
                    for j in range(NUM_POINTS):
                        (j_region, j_n, j_species) = REVERSE_POINT_IDS[j]
                        if j_region == region and j_n == n and j_species == species: # basically i=j
                            J[i,j] += diff/DELTA_R + species.permeability_constant
                        elif j_region == region and j_species == species and j_n == n-1:
                            J[i,j] += -diff/DELTA_R
                        elif j_region == next_region and j_species == species and j_n == 0:
                            J[i,j] += species.permeability_constant
                elif MEMBRANE_TYPE == "enzymatic":
                    raise NotImplementedError("Enzymatic")
    if fill_jacobian:
        return F, J
    else:
        return F, _

def compute_newton_step(species_concentrations):
    """ Returns the vector of the residual, the jacobian as a sparse matrix, and delta u vector
    found by solving for the change in the concentrations.
    """
    # Step 1: Compute residual F and jacobian J
    F_vector, J_matrix = define_newton_residual_and_optionally_jacobian(species_concentrations)
    # Step 2: Assemble J (sparse), convert to solver-friendly format
    J_sparse = J_matrix.tocsc()
    # Step 3: Solve J du = -F
    du = spsolve(J_sparse, -F_vector)
    return F_vector, J_sparse, du

def adaptive_newton_step(
    species_concentrations,
    alpha_current,
    successive_unsuccessful_steps,
    adaptive_step_parameters,
    ):
    """
    Perform one Newton step with adaptive step length (alpha may be >1).
    Returns (next_species_concentrations, alpha_current) where info is dict with diagnostics.
    """
    # Step 0: Unpack parameters
    initial_alpha = adaptive_step_parameters.get("initial_alpha")
    alpha_min = adaptive_step_parameters.get("alpha_min")
    alpha_max = adaptive_step_parameters.get("alpha_max")
    gamma_inc = adaptive_step_parameters.get("gamma_inc")
    gamma_dec = adaptive_step_parameters.get("gamma_dec")
    max_backtrack = adaptive_step_parameters.get("max_num_backtrack")
    max_accepted_successive_unsuccessful_steps = adaptive_step_parameters.get("max_num_accepted_successive_unsuccessful_steps")
    # Steps 1-3: Compute du (see compute_newton_step)
    F_vector, _, du = compute_newton_step(species_concentrations)
    norm_F_vector = np.linalg.norm(F_vector)
    if norm_F_vector == 0:
        return species_concentrations, alpha_current, 0, 0
    # Step 4: Attempt alpha > 1 first (grow from previous alpha_current)
    alpha_try = min(alpha_current * gamma_inc, alpha_max)
    success = False
    for _ in range(max_backtrack):
        species_concentrations_try = copy.deepcopy(species_concentrations)
        for i, du_value in enumerate(du):
            (i_region, i_n, i_species) = REVERSE_POINT_IDS[i]
            species_concentrations_try[i_region][i_n][i_species] += alpha_try * du_value
        # Step 5: Check that all new concentrations are still positive (=0 included)
        # If some negative concentrations, decrease alpha. If alpha is already really small,
        # will later exit without success
        if not all_non_negative(species_concentrations_try):
            alpha_try *= gamma_dec
            if alpha_try < alpha_min:
                break
            continue
        # Step 6: Compute trial residual
        F_vector_try, _ = define_newton_residual_and_optionally_jacobian(
            species_concentrations_try, fill_jacobian=False)
        norm_F_vector_try = np.linalg.norm(F_vector_try)
        # Step 7: Accept if residual decreased
        if norm_F_vector_try < norm_F_vector:
            species_concentrations = species_concentrations_try
            alpha_current = min(alpha_current * gamma_inc, alpha_max)
            success = True
            successive_unsuccessful_steps = 0
            norm_F_to_return = norm_F_vector_try
            break
        # otherwise shrink alpha and retry
        alpha_try *= gamma_dec
        if alpha_try < alpha_min:
            break

    if success is False:
        successive_unsuccessful_steps += 1
        if successive_unsuccessful_steps > max_accepted_successive_unsuccessful_steps:
            raise ValueError("Newton failed")
        # in case that the backtracking did not work, set alpha_current to initial value
        alpha_current = initial_alpha
        for i, du_value in enumerate(du):
            (region, n, species) = REVERSE_POINT_IDS[i]
            species_concentrations[region][n][species] += alpha_current * du_value
        norm_F_to_return = norm_F_vector
        
    return species_concentrations, alpha_current, successive_unsuccessful_steps, norm_F_to_return

def check_convergence(current_species_concentrations, convergence_parameters, print_info):
    """
    Returns true if convergence fulfilled (see below); false if not
    """
    info = {"max_relative_change": 0, "max_Delta_u":0.0, "F_vector_norm":0.0}
    convergence = True
    # Step 0: Unpack parameters
    tol_rel = convergence_parameters.get("tol_relative", 1)
    tol_abs = convergence_parameters.get("tol_absolute", 1)
    tol_res = convergence_parameters.get("tol_residual", 1)
    # Step 1: Get du from concentrations
    F_vector, _, du = compute_newton_step(current_species_concentrations)
    # Step 2: Check that the norm of the residual is small
    F_vector_norm = np.linalg.norm(F_vector)
    if print_info:
        info["F_vector_norm"] = F_vector_norm
        info["max_Delta_u"] = max(du)
    if F_vector_norm > tol_res:
        convergence = False
    # Step 3: Check that each node has had a very small relative change
    # (In case the node has a very small value, have the change be smaller than some absolute value)
    for i, du_value in enumerate(du):
        (region, n, species) = REVERSE_POINT_IDS[i]
        node_u = current_species_concentrations[region][n][species]
        max_tolerated_relative_change = tol_rel * node_u
        if print_info:
            info["max_relative_change"] = max(info["max_relative_change"], max_tolerated_relative_change)
        if du_value > max(tol_abs, max_tolerated_relative_change):
            convergence = False
    if print_info:
        print({k: f"{v:.2e}" for k, v in info.items()})
    return convergence

def solve_newton(
        simulation_start_time,
        initial_iteration_number,
        max_num_newton_iterations,
        initial_species_concentrations,
        adaptive_step_parameters,
        convergence_parameters,
        variables_to_save_dictionary,
        save_data_every=1000,
        check_convergence_every=1000,
        adaptive=True,
        print_convergence_info=False,
        print_iteration_info_every=0,
        plot_iteration_data_during_simulation=False
    ):
    """
    save_data_every and check_convergence_every N iterations. If not to be done, set each to 0.
    """
    current_species_concentrations = initial_species_concentrations
    current_alpha = adaptive_step_parameters["initial_alpha"]
    current_successive_unsuccessful_steps = 0
    early_convergence = False
    for iter in tqdm(range(initial_iteration_number, int(max_num_newton_iterations)),
                     file=sys.stderr,
                     total=int(max_num_newton_iterations),
                     initial=initial_iteration_number):
        # Improve species concentration estimate
        if adaptive == False:
            current_F, _, du = compute_newton_step(current_species_concentrations)
            for i, du_value in enumerate(du):
                (region, n, species) = REVERSE_POINT_IDS[i]
                current_species_concentrations[region][n][species] +=  du_value
            if print_iteration_info_every != 0 and iter%print_iteration_info_every==0 :
                print(f"No step adaptation:\n"
                      f"iteration: {iter}, norm of F: {format_sci(np.linalg.norm(current_F))}\n"
                      f"after {time.time() - simulation_start_time:.3f} seconds of runtime.\n", flush=True
                )
        else:
            current_species_concentrations, current_alpha, current_successive_unsuccessful_steps, current_F_norm = adaptive_newton_step(
                current_species_concentrations, current_alpha, current_successive_unsuccessful_steps, adaptive_step_parameters)
            if print_iteration_info_every != 0 and iter%print_iteration_info_every==0 :
                print(f"Step adaptation:\n"
                      f"iteration: {iter}, norm of F: {current_F_norm},\n"
                      f"alpha: {current_alpha}, current successive unsuccessful steps: {current_successive_unsuccessful_steps}\n",
                      f"after {time.time() - simulation_start_time:.3f} seconds of runtime.\n", flush=True
                )
        # Save result if needed
        if save_data_every !=0 and (iter+1)%save_data_every==0:
            F_vector, J_matrix, du = compute_newton_step(current_species_concentrations)
            iter_string = str(iter).zfill(NUM_NEWTON_ITERATIONS_DIGITS)
            save_newton_iteration_data(iter_string, J_matrix, F_vector, current_species_concentrations, du, variables_to_save_dictionary)
            if plot_iteration_data_during_simulation:
                plot_steady_state_concentrations(
                    os.path.join(ITERATION_DATA_PATH, f".iteration_nr_{iter_string}_plot.png"),
                    current_species_concentrations,
                    title = (
                            f"iteration #{iter}\n"
                            f"residual norm: {format_sci(np.linalg.norm(F_vector))}\n"
                            f"max absolute step: {format_sci(max(du))}"
                        )
                )
        # Stop iterating if criterion for convergence fulfilled
        if check_convergence_every !=0 and iter%check_convergence_every==0:
            convergence = check_convergence(
                current_species_concentrations,
                convergence_parameters,
                print_info=print_convergence_info)
            if convergence:
                print(f"Convergence after {iter} iterations.")
                early_convergence = True
                break
    return current_species_concentrations, early_convergence





if __name__ == "__main__":
    # Parse arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("folder_to_solve", type=str, help="Path to folder with system info")
    parser.add_argument("solver_input_file", type=str)
    parser.add_argument("solver_params_file", type=str)
    # use int(float()) to be able to pass scientific notation
    parser.add_argument("--max-iterations", type=lambda x: int(float(x)), help="Maximum Newton iterations") 
    args = parser.parse_args()

    # Load all the passed information, create folder for background info
    FOLDER_TO_SOLVE = args.folder_to_solve
    ITERATION_DATA_PATH = os.path.join(FOLDER_TO_SOLVE, "solver_iteration_data")
    SOLVER_PARAMS = load_json(args.solver_params_file)
    SOLVER_INPUT = load_json(args.solver_input_file)
    MAX_NUM_NEWTON_ITERATIONS = args.max_iterations

    # Load inputs and define global parameters
    REACTION_NETWORK = pickle_load_binary(os.path.join(FOLDER_TO_SOLVE, ".pickled_REACTION_NETWORK"))
    SYSTEM_GEOMETRY_DICT = load_json(os.path.join(FOLDER_TO_SOLVE, ".expanded_system_geometry"))
    SYSTEM_MESH_DICT= load_json(os.path.join(FOLDER_TO_SOLVE, ".expanded_system_mesh"))
    
    if SOLVER_INPUT["OUTPUT_OPTIONS"]["create_gif_with_saved_data"] is True and SOLVER_INPUT["VARIABLES_TO_SAVE"]["save_concentrations"] is False:
        raise ValueError("Cannot make the gif if the concentrations are not saved.")


    # Lookup to be able to match the species from the species name saved in .json files  
    SPECIES_LOOKUP = {sp.name: sp for sp in REACTION_NETWORK.species}

    # Deal with case permeability vs enzymatic
    # Read out type of membrane
    if hasattr(REACTION_NETWORK.species[0], "permeability_constant"):
        MEMBRANE_TYPE = "permeability"
    elif hasattr(REACTION_NETWORK.species[0], "k_on"):
        MEMBRANE_TYPE = "enzymatic"
    else:
        raise ValueError("Membrane type not correctly specified.")
    if MEMBRANE_TYPE == "enzymatic":
        PORE_DENSITY = SYSTEM_GEOMETRY_DICT["MEMBRANE_PROPERTIES"]["pore_density"]
    
    # Step 1: Define all geometry variables
    R = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["outer_membrane_radius"]
    MESH_POINTS_IN_REGIONS = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["MESH_POINTS_IN_REGIONS"]
    NUM_MESH_POINTS_IN_REGIONS = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["NUM_MESH_POINTS_IN_REGIONS"]
    NUM_REGIONS = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["NUM_REGIONS"]
    MEMBRANE_RADII = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["MEMBRANE_RADII"]
   
    # Step 2: Define structures to access geometry information
    POINT_IDS = SYSTEM_MESH_DICT["POINT_IDS"]
    REVERSE_POINT_IDS = SYSTEM_MESH_DICT["REVERSE_POINT_IDS"]
    RADII = SYSTEM_MESH_DICT["RADII"]
    DELTA_R = SYSTEM_MESH_DICT["DELTA_R"]
    NUM_POINTS = SYSTEM_MESH_DICT["NUM_POINTS"]
    POINT_INFOS = SYSTEM_MESH_DICT["POINT_INFOS"]
    NEIGHBORS = SYSTEM_MESH_DICT["NEIGHBORS"]

    # Step 3: Put enzyme location information
    ENZYMES_CONCENTRATIONS = {
    region_idx : {
        enzyme : enzyme.concentration if region_idx in enzyme.regions else 0
        for enzyme in REACTION_NETWORK.enzymes
    }
    for region_idx in range(NUM_REGIONS)
    }

    # Step 4: Define structure that saves concentrations at each point and which
    # is updated with every iteration of Newton
    species_concentrations_guess = {
        region_idx : {
            mesh_point_idx : {
                species : species.external_concentration * RADII[region_idx][mesh_point_idx] / RADII[NUM_REGIONS-1][NUM_MESH_POINTS_IN_REGIONS[NUM_REGIONS-1]-1]
                for species in REACTION_NETWORK.species}
            for mesh_point_idx in range(NUM_MESH_POINTS_IN_REGIONS[region_idx])}
        for region_idx in range(NUM_REGIONS)
    }
    
    # Guess is used to find the order of magnitude of convergence conditions
    max_guess_concentration = nested_max(species_concentrations_guess)
    F_vector_guess, _, _ = compute_newton_step(species_concentrations_guess)
    # Load any previous solution
    BASENAME_PREVIOUS_SOLUTION = os.path.basename(PREVIOUS_SOLUTION)
    if "none" not in BASENAME_PREVIOUS_SOLUTION:
        previous_solution_species_concentrations_dict_with_strings = load_json(PREVIOUS_SOLUTION)
        species_concentrations_initial = get_species_concentrations_from_json_file(
            previous_solution_species_concentrations_dict_with_strings)
        match = re.search(r"(\d+)", BASENAME_PREVIOUS_SOLUTION)
        if match:
            INITIAL_ITERATION_NUMBER = int(match.group(1))
        else:
            raise ValueError("Could not get the iter from the previous solution")
        rename_iteration_files(ITERATION_DATA_PATH, max_digits=NUM_NEWTON_ITERATIONS_DIGITS, dry_run=False)
        
    else:
        species_concentrations_initial = species_concentrations_guess
        INITIAL_ITERATION_NUMBER = 0
    
    # Step 5: Define convergence criterion
    convergence_parameters = {
        "tol_relative":SOLVER_PARAMS["CONVERGENCE_PARAMETERS"]["tol_relative_value"],
        # tol_absolute is the tolerance of the maximum value of Delta u;
        # max_guess_concentration gives the order of magnitude in which solutions are expected to be
        "tol_absolute":max_guess_concentration*SOLVER_PARAMS["CONVERGENCE_PARAMETERS"]["tol_absolute_factor"],
        "tol_residual":np.linalg.norm(F_vector_guess)*SOLVER_PARAMS["CONVERGENCE_PARAMETERS"]["tol_residual_factor"],
    }

    # Important: save output for checking. tqdm is excluded
    with open(os.path.join(FOLDER_TO_SOLVE, ".newton_solver.log"), "a") as f, redirect_stdout(f):
        print(f"Starting solver from iteration number {INITIAL_ITERATION_NUMBER} \n")
        print("convergence parameters",
            {k: float(f"{v:.2e}") for k, v in convergence_parameters.items()}, "\n")
        # Step 6: Run solver (timed)
        start_time = time.time()
        species_concentrations_final, early_convergence = solve_newton(
            simulation_start_time=start_time,
            initial_iteration_number=INITIAL_ITERATION_NUMBER,
            max_num_newton_iterations=MAX_NUM_NEWTON_ITERATIONS,
            initial_species_concentrations=species_concentrations_initial,
            adaptive_step_parameters=SOLVER_PARAMS["ADAPTIVE_STEP_PARAMETERS"],
            convergence_parameters=convergence_parameters,
            variables_to_save_dictionary = SOLVER_PARAMS["VARIABLES_TO_SAVE"],
            save_data_every=SOLVER_PARAMS["OUTPUT_OPTIONS"]["save_data_every"],
            check_convergence_every=SOLVER_PARAMS["NEWTON_PARAMETERS"]["check_convergence_every"],
            adaptive = not SOLVER_PARAMS["NEWTON_PARAMETERS"]["override_adaptive_method"],
            print_iteration_info_every = SOLVER_PARAMS["OUTPUT_OPTIONS"]["print_iteration_info_every"],
            print_convergence_info = SOLVER_PARAMS["OUTPUT_OPTIONS"]["print_convergence_progress"],
            plot_iteration_data_during_simulation = SOLVER_PARAMS["OUTPUT_OPTIONS"]["plot_iteration_data_during_simulation"]
        )
        end_time = time.time()

        # Log relevant data (if given)
        F_vector_final, _, du_final = compute_newton_step(species_concentrations_final)

        print(f"Runtime was {end_time - start_time:.3f} s\n"
            f"for a residual norm of {format_sci(np.linalg.norm(F_vector_final))},\n"
            f"and a maximum absolute step in concentration of {format_sci(max(du_final))}\n"
            f"with early convergence: {early_convergence}")

    dump_json(FOLDER_TO_SOLVE, ".species_steady_state_concentrations", species_concentrations_final)
    plot_steady_state_concentrations(
        os.path.join(FOLDER_TO_SOLVE, "species_steady_state_concentrations.png"),
        species_concentrations_final)

    # Make gif
    if SOLVER_PARAMS["OUTPUT_OPTIONS"]["create_gif_with_saved_data"]:
        make_newton_iterations_gif(ITERATION_DATA_PATH, FOLDER_TO_SOLVE)
    if SOLVER_PARAMS["OUTPUT_OPTIONS"]["delete_data_at_the_end"]:
        files = glob.glob(os.path.join(ITERATION_DATA_PATH, "*"))
        for f in files:
            os.remove(f)
