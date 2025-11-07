#%%
import sys
import os
import copy
import math
import time
import shutil
from tqdm import tqdm
import numpy as np
from itertools import count
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_matrix, csr_matrix
import matplotlib.pyplot as plt
from auxiliary_functions_using_standard_library import nested_max, all_non_negative, format_sci, pickle_load_binary, closest_value, dump_json, find_sorted_file_names, load_json, find_max_in_nested_dict
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme
from auxiliary_functions import save_matrix_as_sparse_txt
import imageio.v3 as iio
import re

def build_point_ids_dict() -> dict:
    """ Build a nested dict mapping (region, mesh_point, species) to unique IDs.
    """
    counter = count()  # local counter — resets each time you call the function
    point_ids_dict = {
        region_idx: {
            mesh_point_idx: {
                species: next(counter)
                for species in REACTION_NETWORK.species
            }
            for mesh_point_idx in range(NUM_MESH_POINTS_IN_REGIONS[region_idx])
        }
        for region_idx in range(NUM_REGIONS)
    }
    return point_ids_dict

def build_reverse_point_ids_dict(point_ids_dict) -> dict:
    """ Takes the return dict from build_point_ids_dict and constructs a
    inverse dictionary.
    The key is the index of the node and the value is (region, n, species).
    """
    reverse_point_ids_dict = {
        value: (region_idx, mesh_point_idx, species)
        for region_idx, mesh_points in point_ids_dict.items()
        for mesh_point_idx, species_map in mesh_points.items()
        for species, value in species_map.items()
    }
    return reverse_point_ids_dict

def build_radii_dict() -> dict:
    """ Returns a dictionary dict_name[region][n] : radius_to_origin
    """
    radii_dict = {
        region_idx : {
            mesh_point_idx : MESH_POINTS_IN_REGIONS[region_idx][mesh_point_idx]
            for mesh_point_idx in range(NUM_MESH_POINTS_IN_REGIONS[region_idx])
        }
        for region_idx in range(NUM_REGIONS)
    }
    return radii_dict

def build_point_infos_dict() -> dict:
    """ Gives information on whether the node is within the bulk of the region,
    the left-most node within the region or the right-most node within the region.
    dict_name[region][n] : "i" or "l" or "r", respectively
    """
    point_infos_dict = {
        region_idx : {
            mesh_point_idx : "l" if mesh_point_idx==0 else ("r" if mesh_point_idx==NUM_MESH_POINTS_IN_REGIONS[region_idx]-1 else "i")
            for mesh_point_idx in range(NUM_MESH_POINTS_IN_REGIONS[region_idx])
        }
        for region_idx in range(NUM_REGIONS)
    }
    return point_infos_dict

def build_point_neighbor_dict() -> dict:
    """ Gives a list of (region, n) tuples for each [region][n] which specifies
    the node information of spatial neighbors (and itself).
    dict_name[region][n] : [(region, n-1), (region, n), (region, n+1)] e.g. for
    "i" nodes. for "l" and "r" nodes it depends on whether the node is at an inner
    boundary or at r=0 or r=R. Always goes from left to right.
    """
    neighbors_dict = {}

    for region_idx, mesh_points in POINT_INFOS.items():
        for n, kind in mesh_points.items():
            if kind == "i":
                # Internal: previous, self, next
                neighbors_dict[(region_idx, n)] = [
                    (region_idx, n - 1),
                    (region_idx, n),
                    (region_idx, n + 1),
                ]
            elif kind == "l":
                # Left boundary: connect to previous region's rightmost point (if it exists)
                if region_idx > 0:
                    prev_region_last = NUM_MESH_POINTS_IN_REGIONS[region_idx - 1] - 1
                    neighbors_dict[(region_idx, n)] = [
                        (region_idx - 1, prev_region_last),
                    ]
                else:
                    neighbors_dict[(region_idx, n)] = []
                neighbors_dict[(region_idx, n)].append( (region_idx, 0) )
                neighbors_dict[(region_idx, n)].append( (region_idx, 1) )
            elif kind == "r":
                # Right boundary: self’s region last, then 0 and 1 in same region
                last_in_region = NUM_MESH_POINTS_IN_REGIONS[region_idx] - 1
                neighbors_dict[(region_idx, n)] = [
                    (region_idx, last_in_region-1),
                    (region_idx, last_in_region),
                ]
                if region_idx < NUM_REGIONS-1:
                    neighbors_dict[(region_idx, n)].append( (region_idx+1, 0) )
    return neighbors_dict

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

def save_newton_iteration_data(
    iter_string, J, F, species_concentrations_to_save, max_Du_to_save):
    dump_json(ITERATION_DATA_PATH, f".iteration_nr_{iter_string}_concentration",
            species_concentrations_to_save)
    dump_json(ITERATION_DATA_PATH, f".iteration_nr_{iter_string}_max_Du",
            max_Du_to_save)
    np.savetxt(os.path.join(ITERATION_DATA_PATH, f".iteration_nr_{iter_string}_F.txt"), F, fmt="%.15e", delimiter="\n")
    save_matrix_as_sparse_txt(J, os.path.join(ITERATION_DATA_PATH, f".iteration_nr_{iter_string}_J"))

def compute_newton_step(species_concentrations):
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
        return species_concentrations, alpha_current
    # Step 4: Attempt alpha > 1 first (grow from previous alpha_current)
    alpha_try = min(alpha_current * gamma_inc, alpha_max)
    success = False
    for _ in range(max_backtrack):
        species_concentrations_try = copy.deepcopy(species_concentrations)
        for i, du_value in enumerate(du):
            (region, n, species) = REVERSE_POINT_IDS[i]
            species_concentrations_try[region][n][species] += alpha_try * du_value   
        # Step 5: Check that all new concentrations are still positive (=0 included)
        # If some negative concentrations, decrease alpha. If alpha is already really small,
        # will later exit without success
        if not all_non_negative(species_concentrations_try):
            print("negative_values!?")
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
            alpha_current = min(alpha_try, alpha_max)
            # optionally enlarge alpha for next iter
            alpha_current = min(alpha_current * gamma_inc, alpha_max)
            success = True
            successive_unsuccessful_steps = 0
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
        
        
    return species_concentrations, alpha_current, successive_unsuccessful_steps

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
        max_num_newton_iterations,
        initial_species_concentrations_guess,
        adaptive_step_parameters,
        convergence_parameters,
        save_data_every=1000,
        check_convergence_every=1000,
        adaptive=True,
        print_convergence_info=False
    ):
    """
    save_data_every and check_convergence_every N iterations. If not to be done, set each to 0.
    """
    current_species_concentrations = initial_species_concentrations_guess
    current_alpha = adaptive_step_parameters["initial_alpha"]
    current_successive_unsuccessful_steps = 0
    early_convergence = False
    for iter in tqdm(range(int(max_num_newton_iterations))):
        # Improve species concentration estimate
        if adaptive == False:
            _, _, du = compute_newton_step(current_species_concentrations)
            for i, du_value in enumerate(du):
                (region, n, species) = REVERSE_POINT_IDS[i]
                current_species_concentrations[region][n][species] +=  du_value
        else:
            current_species_concentrations, current_alpha, current_successive_unsuccessful_steps = adaptive_newton_step(
                current_species_concentrations, current_alpha, current_successive_unsuccessful_steps, adaptive_step_parameters)   
        # Save result if needed
        if save_data_every !=0 and iter%save_data_every==0:
            F_vector, J_matrix, du = compute_newton_step(current_species_concentrations)
            iter_string = str(iter).zfill(int(math.log10(max_num_newton_iterations)+1))
            save_newton_iteration_data(iter_string, J_matrix, F_vector, current_species_concentrations, max(du))
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

def plot_steady_state_concentrations(output_file_name, species_concentrations_to_plot, title = None, ymax = None):
    x_values = []
    y_values = {}
    for species_idx, species in enumerate(REACTION_NETWORK.species):
        species_y_values = []
        for region in range(NUM_REGIONS):
            for n in range(NUM_MESH_POINTS_IN_REGIONS[region]):
                if species_idx == 0:
                    x_values.append(RADII[region][n])
                species_y_values.append(species_concentrations_to_plot[region][n][species])
        y_values[species] = species_y_values

    fig, ax = plt.subplots(1,1, figsize = (5,3))
    for x_value in x_values:
        ax.axvline(x_value/max(x_values), ymin = 0.95, ymax = 1, color="k")
    for species in REACTION_NETWORK.species:
        curve, = ax.plot(x_values/max(x_values), y_values[species], label=species.name)
        color = curve.get_color()
        ax.hlines(species.external_concentration, xmin=1, xmax = 1.1, color = color)
    ax.set_ylabel("concentration / M")
    ax.set_xlabel("relative distance to origin / r/R")
    ax.legend(
        loc='upper center',      # anchor point of legend
        bbox_to_anchor=(0.5, -0.25),  # (x, y) position in figure coordinates
        ncol=3,                  # number of columns
        frameon=False
    )
    for membrane_radius in MEMBRANE_RADII:
        ax.axvline(membrane_radius/max(MEMBRANE_RADII), linestyle = "--", alpha = 0.5, c = "k")

    max_value = max(max(y_values[species]) for species in REACTION_NETWORK.species)
    if ymax == None:
        ymax = max_value * 1.05
    if title != None:
        ax.set_title(title, loc="left")
    ax.set_ylim(ymin=0, ymax = ymax)
    ax.set_xlim(xmin=0, xmax = 1.1)
    fig.savefig(output_file_name, dpi = 300, bbox_inches='tight')
    plt.close(fig)

def make_newton_iterations_gif(iteration_data_folder, gif_output_folder):
    """
    Important: delete any .json files from previous simulations that may not be overwritten.
    """
    sorted_files = find_sorted_file_names(iteration_data_folder, ".iteration_nr_*_concentration.json")
    max_concentration_value = 0
    for file in sorted_files:
        concentration_dict = load_json(file)
        max_value = find_max_in_nested_dict(concentration_dict)
        if max_value > max_concentration_value:
            max_concentration_value = max_value
    for species in REACTION_NETWORK.species:
        if species.external_concentration > max_concentration_value:
            max_concentration_value = species.external_concentration
    max_concentration_value *= 1.1
    png_files_created = []
    # to put the species object back in the dictionary  
    species_lookup = {sp.name: sp for sp in REACTION_NETWORK.species}
    for file in sorted_files:
        species_concentrations_to_plot_dict_with_strings = load_json(file)
        max_Du_file = file.replace("_concentration.json", "_max_Du.json")
        max_Du = load_json(max_Du_file)
        # Make keys that got converted to strings instead of integers into integers again
        species_concentrations_to_plot_dict = {
            int(region_idx): {
                int(mesh_point_idx): {
                    species_lookup[species_name]: data
                    for species_name, data in mesh_point_info.items()}
                for mesh_point_idx, mesh_point_info in region_info.items()}
            for region_idx, region_info in species_concentrations_to_plot_dict_with_strings.items()
        }
        png_file = os.path.splitext(file)[0] + ".png" # remove .json
        number = re.findall(r"\d+", os.path.basename(png_file))[0]
        max_Du_plottable = format_sci(max_Du)
        plot_steady_state_concentrations(png_file, species_concentrations_to_plot_dict, title = f"iteration # {number}, max(Delta concentration) = {max_Du_plottable}", ymax = max_concentration_value)
        png_files_created.append(png_file)
    file_to_create = os.path.join(gif_output_folder, "newton_iterations.gif")
    print("creating gif", file_to_create)
    with iio.imopen(file_to_create, "w") as writer:
        for filename in png_files_created:
            image = iio.imread(filename)
            writer.write(image)

if __name__ == "__main__":
    # Load all the information
    FOLDER_TO_SOLVE = sys.argv[1]
    # Create a folder in which to save iteration data
    ITERATION_DATA_PATH = os.path.join(FOLDER_TO_SOLVE, ".solver_iteration_data")
    os.makedirs(ITERATION_DATA_PATH, exist_ok=True)

    REACTION_NETWORK = pickle_load_binary(os.path.join(FOLDER_TO_SOLVE, ".REACTION_NETWORK_pickle"))
    SYSTEM_GEOMETRY_DICT = pickle_load_binary(os.path.join(FOLDER_TO_SOLVE, ".SYSTEM_GEOMETRY_pickle"))
    # Read out type of membrane
    if hasattr(REACTION_NETWORK.species[0], "permeability_constant"):
        MEMBRANE_TYPE = "permeability"
    elif hasattr(REACTION_NETWORK.species[0], "k_on"):
        MEMBRANE_TYPE = "enzymatic"
    else:
        raise ValueError("Membrane type not correctly specified.")

    # Step 0: Get all solver parameters
    SOLVER_PARAMS = pickle_load_binary(os.path.join(FOLDER_TO_SOLVE, ".solver_info_pickle"))

    # Step 1: Define all geometry variables
    R = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["outer_membrane_radius"]
    MESH_POINTS_IN_REGIONS = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["MESH_POINTS_IN_REGIONS"]
    NUM_MESH_POINTS_IN_REGIONS = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["NUM_MESH_POINTS_IN_REGIONS"]
    NUM_REGIONS = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["NUM_REGIONS"]
    MEMBRANE_RADII = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["MEMBRANE_RADII"]
    if MEMBRANE_TYPE == "enzymatic":
        PORE_DENSITY = SYSTEM_GEOMETRY_DICT["MEMBRANE_PROPERTIES"]["pore_density"]
    

    # Step 2: Define structures to access geometry information
    POINT_IDS = build_point_ids_dict()
    REVERSE_POINT_IDS = build_reverse_point_ids_dict(POINT_IDS)
    RADII = build_radii_dict()
    DELTA_R = RADII[0][1]-RADII[0][0] # the different points within a region are equally spaced
    NUM_POINTS = len(REVERSE_POINT_IDS) # each point saves the concentration for one species at one node
    POINT_INFOS = build_point_infos_dict()
    NEIGHBORS = build_point_neighbor_dict()
    # Save dictionaries in .json files for readability
    dump_json(FOLDER_TO_SOLVE, ".solver_POINT_IDS", POINT_IDS)
    dump_json(FOLDER_TO_SOLVE, ".solver_REVERSE_POINT_IDS", REVERSE_POINT_IDS)
    dump_json(FOLDER_TO_SOLVE, ".solver_RADII", RADII)
    dump_json(FOLDER_TO_SOLVE, ".solver_POINT_INFOS", POINT_INFOS)
    dump_json(FOLDER_TO_SOLVE, ".solver_NEIGHBORS", NEIGHBORS)

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
    max_guess_concentration = nested_max(species_concentrations_guess)
    F_vector_guess, _, _ = compute_newton_step(species_concentrations_guess)
    
    # Step 5: Define convergence criterion
    convergence_parameters = {
        "tol_relative":SOLVER_PARAMS["CONVERGENCE_PARAMETERS"]["tol_relative_value"],
        # tol_absolute is the tolerance of the maximum value of Delta u;
        # max_guess_concentration gives the order of magnitude in which solutions are expected to be
        "tol_absolute":max_guess_concentration*SOLVER_PARAMS["CONVERGENCE_PARAMETERS"]["tol_absolute_factor"],
        "tol_residual":np.linalg.norm(F_vector_guess)*SOLVER_PARAMS["CONVERGENCE_PARAMETERS"]["tol_residual_factor"],
    }

    print("convergence parameters",
          {k: float(f"{v:.2e}") for k, v in convergence_parameters.items()})

    # Step 6: Run solver (timed)
    start_time = time.time()
    species_concentrations_final, early_convergence = solve_newton(
        max_num_newton_iterations=SOLVER_PARAMS["NEWTON_PARAMS"]["max_num_newton_iterations"],
        initial_species_concentrations_guess=species_concentrations_guess,
        adaptive_step_parameters=SOLVER_PARAMS["ADAPTIVE_STEP_PARAMETERS"],
        convergence_parameters=convergence_parameters,
        save_data_every=SOLVER_PARAMS["OUTPUT_OPTIONS"]["save_data_every"],
        check_convergence_every=SOLVER_PARAMS["NEWTON_PARAMS"]["check_convergence_every"],
        adaptive = not SOLVER_PARAMS["NEWTON_PARAMS"]["override_adaptive_method"],
        print_convergence_info = SOLVER_PARAMS["OUTPUT_OPTIONS"]["print_convergence_progress"]
    )
    end_time = time.time()
    
    # Print run time (if given)
    F_vector_final, _, _ = compute_newton_step(species_concentrations_final)
    print(f"Runtime was {end_time - start_time:.3f} s for a residual norm of {np.linalg.norm(F_vector_final)}, with early convergence: {early_convergence}")
    
    # Save final concentration
    dump_json(FOLDER_TO_SOLVE, ".species_steady_state_concentrations", species_concentrations_final)
    # Make gif
    if SOLVER_PARAMS["OUTPUT_OPTIONS"]["create_gif_with_saved_data"]:
        make_newton_iterations_gif(ITERATION_DATA_PATH, FOLDER_TO_SOLVE)
    if SOLVER_PARAMS["OUTPUT_OPTIONS"]["delete_data_at_the_end"]:
        shutil.rmtree(ITERATION_DATA_PATH)



