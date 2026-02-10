#%%
import sys
import os
import copy
import time
import math
import re
from contextlib import redirect_stdout
import argparse
from tqdm import tqdm
import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_matrix
import matplotlib.pyplot as plt
from auxiliary_functions_using_standard_library import (
    find_max_in_nested_dict, all_non_negative, format_sci, pickle_load_binary,
    load_json, CSVLogger)
from auxiliary_functions import dump_json
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme
from plot_bvp_solution import plot_steady_state_concentrations, make_newton_iterations_gif
from auxiliary_functions_framework_organization import (
    get_species_concentrations_from_json_file,
    get_correct_point_ids_dict, get_correct_reverse_point_ids_dict, get_correct_neighbors_dict)
from auxiliary_functions_framework_organization_using_standard_library import (
    find_latest_solution, rename_iteration_files)
from auxiliary_functions_using_scipy import save_newton_iteration_data
from auxiliary_functions import read_yaml_file

def calculate_reaction_term(current_species_concentrations, region, n, species):
    """ Gives the reaction term for F_i (for a specific species at a specific point in the mesh).
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
    Used for defining the jacobian.
    """
    if isinstance(reaction_to_check, SpontaneousReaction):
        derivative = reaction_to_check.k
    elif isinstance(reaction_to_check, EnzymaticReaction):
        derivative = reaction_to_check.k_cat * ENZYMES_CONCENTRATIONS[region][reaction_to_check.enzyme] * reaction_to_check.k_M / ( reaction_to_check.k_M + current_species_concentrations[region][n][partial_derivative_species])
    if partial_derivative_species == reaction_to_check.start_species:
        derivative *= -1
    return derivative

def get_pore_density_occupation_information(current_species_concentrations, info_minus_side, info_plus_side):
    """ Used for the enzymatic membrane model.
    Gives necessary information about the density of occupied pores by each species, when given the current
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

def check_convergence_via_jacobian_and_residual(
        current_species_concentrations, convergence_parameters, get_full_info):
    """
    Returns true if convergence fulfilled (see below); false if not
    """
    if get_full_info:
        info = {"max_relative_change": 0, "max_Delta_u":0.0, "F_vector_norm":0.0}
    else:
        info = {}
    convergence = True
    # Step 0: Unpack parameters
    tol_rel = convergence_parameters.get("tol_relative", 1)
    tol_abs = convergence_parameters.get("tol_absolute", 1)
    tol_res = convergence_parameters.get("tol_residual", 1)
    # Step 1: Get du from concentrations
    F_vector, _, du = compute_newton_step(current_species_concentrations)
    # Step 2: Check that the norm of the residual is small
    F_vector_norm = np.linalg.norm(F_vector)
    if get_full_info:
        info["F_vector_norm"] = F_vector_norm
        info["max_Delta_u"] = max(du)
    if F_vector_norm > tol_res:
        convergence = False
        if not get_full_info: # early return in case no convergence and no need for full info
            return convergence, {}
    # Step 3: Check that each node has had a very small relative change
    # (In case the node has a very small value, have the change be smaller than some absolute value)
    for i, du_value in enumerate(du):
        (region, n, species) = REVERSE_POINT_IDS[i]
        node_u = current_species_concentrations[region][n][species]
        max_tolerated_relative_change = tol_rel * node_u
        if get_full_info:
            info["max_relative_change"] = max(info["max_relative_change"], max_tolerated_relative_change)
        if du_value > max(tol_abs, max_tolerated_relative_change):
            convergence = False
            if not get_full_info: # early return in case no convergence and no need for full info
                return convergence, {}
    return convergence, info


def check_convergence_via_flux_equilibrium(
        current_species_concentrations, convergence_parameters, get_full_info):
    """
    Returns True if the relative excess of flux in any direction is smaller than some
    deviation, given as the value to the key tol_relative_flux_deviation in
    the dictionary convergence_parameters
    """
    if get_full_info:
        info = {species: 0 for species in REACTION_NETWORK.species}
    else:
        info = {}
    convergence = True
    tol_relative_flux_deviation = convergence_parameters.get("tol_relative_flux_deviation", 1)
    # First, calculate net reaction fluxes within the sphere
    reaction_fluxes = {species: 0
        for species in REACTION_NETWORK.species}
    for i in range(NUM_POINTS):
        (region, n, species) = REVERSE_POINT_IDS[i]
        r = RADII[region][n]
        reaction_flux = calculate_reaction_term(current_species_concentrations, region, n, species)
        reaction_fluxes[species] += 4 * np.pi * reaction_flux * r**2
    # Second, calculate flux from boundary with exterior
    # the flux is positive if the concentration on the exterior is larger than on the interior at r=R
    boundary_fluxes = {species: 0
        for species in REACTION_NETWORK.species}
    for species in REACTION_NETWORK.species:
        boundary_fluxes[species] = species.permeability_constant * (
            species.external_concentration
            - current_species_concentrations[NUM_REGIONS-1][NUM_MESH_POINTS_IN_REGIONS[NUM_REGIONS-1]-1][species])
    # Since we are simulating the steady state, we want the total net flux to be 0 for each species
    # Because of numerics, we need some tolerance
    for species in REACTION_NETWORK.species:
        relative_deviation = abs(reaction_fluxes[species] + boundary_fluxes[species]) / max(abs(boundary_fluxes[species]), abs(reaction_fluxes[species]))
        if get_full_info:
            info[species] = relative_deviation
        if relative_deviation > tol_relative_flux_deviation:
            convergence = False
            if not get_full_info: # if we do not need the full info, early return
                return convergence, {}
    return convergence, info


def solve_newton(
        simulation_start_time,
        initial_iteration_number,
        max_num_newton_iterations,
        initial_species_concentrations,
        adaptive_step_parameters,
        convergence_parameters,
        variables_to_save_dictionary,
        save_data_every,
        check_convergence_every,
        adaptive,
        log_convergence_info,
        convergence_info_logger_path,
        log_iteration_info_every,
        plot_iteration_data_during_simulation
    ):
    """
    save_data_every and check_convergence_every N iterations. If not to be done, set each to 0.
    """
    current_species_concentrations = initial_species_concentrations
    current_alpha = adaptive_step_parameters["initial_alpha"]
    current_successive_unsuccessful_steps = 0
    early_convergence = False # tracks whether the system fin
    convergence_info_logger = CSVLogger(convergence_info_logger_path)
    for iter in tqdm(range(initial_iteration_number, int(max_num_newton_iterations)),
                     file=sys.stderr,
                     total=int(max_num_newton_iterations),
                     initial=initial_iteration_number):
        # Improve species concentration estimate
        if adaptive == False:
            _, _, du = compute_newton_step(current_species_concentrations)
            for i, du_value in enumerate(du):
                (region, n, species) = REVERSE_POINT_IDS[i]
                current_species_concentrations[region][n][species] +=  du_value
            if log_iteration_info_every != 0 and iter%log_iteration_info_every==0 :
                print(f"No step adaptation:\n"
                      f"iteration: {iter}\n"
                      f"after {time.time() - simulation_start_time:.3f} seconds of runtime.\n", flush=True
                )
        else:
            current_species_concentrations, current_alpha, current_successive_unsuccessful_steps, current_F_norm = adaptive_newton_step(
                current_species_concentrations, current_alpha, current_successive_unsuccessful_steps, adaptive_step_parameters)
            if log_iteration_info_every != 0 and iter%log_iteration_info_every==0 :
                print(f"Step adaptation:\n"
                      f"iteration: {iter}\n"
                      f"alpha: {current_alpha}, current successive unsuccessful steps: {current_successive_unsuccessful_steps}\n",
                      f"after {time.time() - simulation_start_time:.3f} seconds of runtime.\n", flush=True
                )
        # Save result if needed
        if save_data_every !=0 and (iter+1)%save_data_every==0:
            F_vector, J_matrix, du = compute_newton_step(current_species_concentrations)
            iter_string = str(iter).zfill(num_iterations_digits)
            save_newton_iteration_data(ITERATION_DATA_PATH, iter_string,
                J_matrix, F_vector, current_species_concentrations, du, variables_to_save_dictionary)
            if plot_iteration_data_during_simulation:
                plot_steady_state_concentrations(
                    reaction_network=REACTION_NETWORK,
                    num_regions=NUM_REGIONS,
                    num_mesh_points_in_regions=NUM_MESH_POINTS_IN_REGIONS,
                    radii=RADII,
                    membrane_radii=MEMBRANE_RADII,
                    output_file_name=os.path.join(ITERATION_DATA_PATH, f".iteration_nr_{iter_string}_iteration.png"),
                    species_concentrations_to_plot=current_species_concentrations,
                    title = (
                            f"iteration #{iter}\n"
                            f"residual norm: {format_sci(np.linalg.norm(F_vector))}\n"
                            f"max absolute step: {format_sci(max(du))}"
                        )
                )
        # Stop iterating if criterion for convergence fulfilled
        if check_convergence_every !=0 and iter%check_convergence_every==0:
            convergence_flux_equilibration, info_flux_equilibration = check_convergence_via_flux_equilibrium(
                current_species_concentrations,
                convergence_parameters,
                get_full_info=log_convergence_info
            )
            convergence_jacobian_residual, info_jacobian_residual = check_convergence_via_jacobian_and_residual(
                current_species_concentrations,
                convergence_parameters,
                get_full_info=log_convergence_info
            )
            if log_convergence_info:
                info = info_flux_equilibration | info_jacobian_residual
                convergence_info_logger.log(iter, info)
            if convergence_flux_equilibration and convergence_jacobian_residual:
                early_convergence = True
                break
    return current_species_concentrations, early_convergence


if __name__ == "__main__":
    # Parse arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("folder_to_solve", type=str, help="Path to folder with system info")
    parser.add_argument("--solver_input_file", type=str)
    parser.add_argument("--solver_params_file", type=str)
    # use int(float()) to be able to pass scientific notation
    parser.add_argument("--max_iterations", type=lambda x: int(float(x)), help="Maximum Newton iterations") 
    args = parser.parse_args()

    # Load all the passed information, create folder for background info
    FOLDER_TO_SOLVE = args.folder_to_solve
    ITERATION_DATA_PATH = os.path.join(FOLDER_TO_SOLVE, "solver_iteration_data")
    SOLVER_PARAMS = read_yaml_file(args.solver_params_file)
    SOLVER_INPUT = read_yaml_file(args.solver_input_file)
    MAX_NUM_NEWTON_ITERATIONS = args.max_iterations
    min_num_iterations_digits = int(math.log10(MAX_NUM_NEWTON_ITERATIONS)+1)
    os.makedirs(ITERATION_DATA_PATH, exist_ok=True)

    # Load inputs and define global parameters
    REACTION_NETWORK = pickle_load_binary(os.path.join(FOLDER_TO_SOLVE, ".pickled_reaction_network"))
    SYSTEM_GEOMETRY_DICT = load_json(os.path.join(FOLDER_TO_SOLVE, ".expanded_system_geometry.json"))
    SYSTEM_MESH_DICT= load_json(os.path.join(FOLDER_TO_SOLVE, ".expanded_system_mesh.json"))
    
    if SOLVER_INPUT["output_options"]["create_gif_with_saved_data"] is True and SOLVER_INPUT["variables_to_save"]["save_concentrations"] is False:
        raise ValueError("Cannot make the gif if the concentrations are not saved.")

    # Deal with case permeability vs enzymatic
    # Read out type of membrane
    if hasattr(REACTION_NETWORK.species[0], "permeability_constant"):
        MEMBRANE_TYPE = "permeability"
    elif hasattr(REACTION_NETWORK.species[0], "k_on"):
        MEMBRANE_TYPE = "enzymatic"
    else:
        raise ValueError("Membrane type not correctly specified.")
    if MEMBRANE_TYPE == "enzymatic":
        PORE_DENSITY = SYSTEM_GEOMETRY_DICT["membrane_properties"]["pore_density"]
    
    # Step 1: Define all geometry variables
    R = SYSTEM_GEOMETRY_DICT["geometry_config"]["outer_membrane_radius"]
    MESH_POINTS_IN_REGIONS = SYSTEM_GEOMETRY_DICT["geometry_config"]["mesh_points_in_regions"]
    NUM_MESH_POINTS_IN_REGIONS = SYSTEM_GEOMETRY_DICT["geometry_config"]["num_mesh_points_in_regions"]
    NUM_REGIONS = SYSTEM_GEOMETRY_DICT["geometry_config"]["num_regions"]
    MEMBRANE_RADII = SYSTEM_GEOMETRY_DICT["geometry_config"]["membrane_radii"]

    SPECIES_LOOKUP = {sp.name: sp for sp in REACTION_NETWORK.species}
  
    # Step 2: Define structures to access geometry information
    POINT_IDS = get_correct_point_ids_dict(SYSTEM_MESH_DICT["point_ids"], SPECIES_LOOKUP)
    REVERSE_POINT_IDS = get_correct_reverse_point_ids_dict(SYSTEM_MESH_DICT["reverse_point_ids"], SPECIES_LOOKUP)
    RADII = SYSTEM_MESH_DICT["radii"]
    DELTA_R = SYSTEM_MESH_DICT["delta_r"]
    NUM_POINTS = SYSTEM_MESH_DICT["num_points"]
    POINT_INFOS = SYSTEM_MESH_DICT["point_infos"] 
    NEIGHBORS = get_correct_neighbors_dict(SYSTEM_MESH_DICT["neighbors"])
    
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
    # currently, the concentration is initially set to increase linearly with the 
    # radius up until its external concentration
    species_concentrations_guess = {
        region_idx : {
            mesh_point_idx : {
                species : species.external_concentration * RADII[region_idx][mesh_point_idx] / RADII[NUM_REGIONS-1][NUM_MESH_POINTS_IN_REGIONS[NUM_REGIONS-1]-1]
                for species in REACTION_NETWORK.species}
            for mesh_point_idx in range(NUM_MESH_POINTS_IN_REGIONS[region_idx])}
        for region_idx in range(NUM_REGIONS)
    }
    
    print("Computing order of magnitude of residual.")
    # Guess is used to find the order of magnitude of convergence conditions
    max_guess_concentration = find_max_in_nested_dict(species_concentrations_guess)
    F_vector_guess, _, _ = compute_newton_step(species_concentrations_guess)

    # Find out whether there is any previous solution from which to continue
    previous_solution = find_latest_solution(ITERATION_DATA_PATH)
    if previous_solution is None:
        print(f"Starting simulation for folder {FOLDER_TO_SOLVE} from scratch.")
        species_concentrations_initial = species_concentrations_guess
        initial_iteration_number = 0
        num_iterations_digits = min_num_iterations_digits
    else:
        previous_solution_species_concentrations_dict_with_strings = load_json(previous_solution)
        species_concentrations_initial = get_species_concentrations_from_json_file(
            previous_solution_species_concentrations_dict_with_strings, SPECIES_LOOKUP)
        basename_previous_solution = os.path.basename(previous_solution)
        match = re.search(r"(\d+)", basename_previous_solution)
        if match:
            initial_iteration_number = int(match.group(1))
        else:
            raise ValueError("Could not get the iter from the previous solution")
        num_iterations_digits = rename_iteration_files(ITERATION_DATA_PATH, min_digits=min_num_iterations_digits)
        print(f"Continuing simulation for folder {FOLDER_TO_SOLVE} from iteration {initial_iteration_number}.")

    
    # Step 5: Define convergence criterion
    convergence_parameters = {
        "tol_relative":SOLVER_PARAMS["convergence_parameters"]["tol_relative_value"],
        # tol_absolute is the tolerance of the maximum value of Delta u;
        # max_guess_concentration gives the order of magnitude in which solutions are expected to be
        "tol_absolute":max_guess_concentration*SOLVER_PARAMS["convergence_parameters"]["tol_absolute_factor"],
        "tol_residual":np.linalg.norm(F_vector_guess)*SOLVER_PARAMS["convergence_parameters"]["tol_residual_factor"],
        "tol_relative_flux_deviation": SOLVER_PARAMS["convergence_parameters"]["tol_relative_flux_deviation"],
    }

    # Important: save output for checking. tqdm is excluded
    print(f"Opening log file for {FOLDER_TO_SOLVE}.")
    with open(os.path.join(FOLDER_TO_SOLVE, ".newton_solver.log"), "a") as f, redirect_stdout(f):
        print(f"Starting solver from iteration number {initial_iteration_number} \n")
        print("convergence parameters",
            {k: float(f"{v:.2e}") for k, v in convergence_parameters.items()}, "\n")
        # Step 6: Run solver (timed)
        start_time = time.time()
        species_concentrations_final, early_convergence = solve_newton(
            simulation_start_time=start_time,
            initial_iteration_number=initial_iteration_number,
            max_num_newton_iterations=MAX_NUM_NEWTON_ITERATIONS,
            initial_species_concentrations=species_concentrations_initial,
            adaptive_step_parameters=SOLVER_INPUT["adaptive_step_parameters"],
            convergence_parameters=convergence_parameters,
            variables_to_save_dictionary = SOLVER_INPUT["variables_to_save"],
            save_data_every=SOLVER_INPUT["output_options"]["save_data_every"],
            check_convergence_every=SOLVER_PARAMS["newton_parameters"]["check_convergence_every"],
            adaptive = not SOLVER_INPUT["newton_parameters"]["override_adaptive_method"],
            log_iteration_info_every = SOLVER_INPUT["output_options"]["log_iteration_info_every"],
            log_convergence_info = SOLVER_INPUT["output_options"]["log_convergence_progress"],
            convergence_info_logger_path=os.path.join(FOLDER_TO_SOLVE, ".convergence_logger.csv"),
            plot_iteration_data_during_simulation = SOLVER_INPUT["output_options"]["plot_iteration_data_during_simulation"]
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
        reaction_network=REACTION_NETWORK,
        num_regions=NUM_REGIONS,
        num_mesh_points_in_regions=NUM_MESH_POINTS_IN_REGIONS,
        radii=RADII,
        membrane_radii=MEMBRANE_RADII,
        output_file_name=os.path.join(FOLDER_TO_SOLVE, "species_steady_state_concentrations.png"),
        species_concentrations_to_plot = species_concentrations_final)

    # Make gif
    if SOLVER_INPUT["output_options"]["create_gif_with_saved_data"]:
        make_newton_iterations_gif(
            reaction_network=REACTION_NETWORK,
            num_regions=NUM_REGIONS,
            num_mesh_points_in_regions=NUM_MESH_POINTS_IN_REGIONS,
            radii=RADII,
            membrane_radii=MEMBRANE_RADII,
            iteration_data_folder=ITERATION_DATA_PATH,
            gif_output_folder=FOLDER_TO_SOLVE,
            species_lookup_dict=SPECIES_LOOKUP)

    #if SOLVER_INPUT["output_options"]["delete_data_at_the_end"]:
    #    files = glob.glob(os.path.join(ITERATION_DATA_PATH, "*"))
    #    for f in files:
    #        os.remove(f)
