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
import scipy
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_matrix
import matplotlib.pyplot as plt
from auxiliary_functions_using_standard_library import (
    find_max_in_nested_dict, all_non_negative, format_sci, pickle_load_binary,
    load_json, CSVLogger, pickle_dump_binary)
from auxiliary_functions import dump_json
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme
from plot_bvp_solution import plot_steady_state_concentrations, make_newton_iterations_gif, plot_convergence_progress
from auxiliary_functions_framework_organization import (
    get_species_concentrations_from_json_file,
    get_correct_point_ids_dict, get_correct_reverse_point_ids_dict, get_correct_neighbors_dict)
from auxiliary_functions_framework_organization_using_standard_library import (
    find_latest_solution, rename_iteration_files)
from auxiliary_functions import read_yaml_file
from enforce_parameter_value_conditions import (
    return_reaction_network_with_total_fixed_quantity_asserted,
    return_enzyme_concentrations,
    assert_no_conflicts_in_enzyme_positioning
)
from auxiliary_functions_using_scipy import save_newton_iteration_data
from study_bvp_solution import get_outward_fluxes



def calculate_reaction_partial_derivative(current_species_concentrations, reaction_to_check, partial_derivative_species, region, n):
    """ Gives the partial derivative of a reaction to a concentration of a species
    that is involved in the reaction.
    Used for defining the jacobian.
    """
    if isinstance(reaction_to_check, SpontaneousReaction):
        derivative = reaction_to_check.k
    elif isinstance(reaction_to_check, EnzymaticReaction):
        S = current_species_concentrations[region][n][reaction_to_check.start_species]
        if partial_derivative_species == reaction_to_check.start_species:
            derivative = (
                reaction_to_check.k_cat
                * ENZYME_CONCENTRATIONS[region][reaction_to_check.enzyme]
                * reaction_to_check.k_M
                / (reaction_to_check.k_M + S)**2
            )
        else:
            derivative = 0.0
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

def calculate_reaction_term(current_species_concentrations, region, n, species):
    """ Gives the reaction term for F_i (for a specific species at a specific point in the mesh).
    """
    reaction_term = 0
    for reaction in species.as_reactant_in + species.as_product_in:
        if isinstance(reaction, SpontaneousReaction):
            term = reaction.k * current_species_concentrations[region][n][reaction.start_species]
        else:
            term = reaction.k_cat * ENZYME_CONCENTRATIONS[region][reaction.enzyme] * current_species_concentrations[region][n][reaction.start_species] / (reaction.k_M + current_species_concentrations[region][n][reaction.start_species])
        if reaction in species.as_reactant_in: # if acts as reactant, diminishes
            term *= -1
        reaction_term += term
    return reaction_term

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
            # j is that by which we are applying the derivative to F_i
            for j in range(NUM_POINTS):
                (j_region, j_n, j_species) = REVERSE_POINT_IDS[j]
                # Contributions from diffusion
                if j_region == region and j_n == n and j_species == species: # j == i, basically
                    J[i,j] += diff * (1/DELTA_R**2 * (-2))
                elif j_region==region and j_n==right_n and j_species == species: # same species, right or left ###################
                    J[i,j] += diff * (1/DELTA_R**2 + 1/(DELTA_R*r))
                elif j_region==region and j_n==left_n and j_species == species: # same species, right or left #####################
                    J[i,j] += diff * (1/DELTA_R**2 - 1/(DELTA_R*r))
                # Contributions from reactions
                # if we are looking at the right position in the mesh (position in the mesh of j is the same as of i)
                if j_region == region and j_n == center_n:
                    # we look at all reactions in which the species from i is involved
                    for reaction in species.as_reactant_in + species.as_product_in:
                        # but only the derivative with respect to the substrate matters
                        if j_species != reaction.start_species:
                            continue
                        # get concentration of substrate
                        S = current_species_concentrations[region][center_n][reaction.start_species]
                        # if the reaction is spontaneous
                        if isinstance(reaction, SpontaneousReaction):
                            base = reaction.k # sign is placed afterwards
                        else: #EnzymaticReaction
                            base = (
                                reaction.k_cat
                                * ENZYME_CONCENTRATIONS[region][reaction.enzyme]
                                * reaction.k_M
                                / (reaction.k_M + S)**2
                            )
                        # Apply sign
                        if species == reaction.start_species:
                            J[i,j] += -base
                        elif species == reaction.end_species:
                            J[i,j] += +base
 
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
                    if j_region == region and j_n == n:# if on the same place but not necessarily the same species
                        # we look at all reactions in which the species from i is involved
                        for reaction in species.as_reactant_in + species.as_product_in:
                            # but only the derivative with respect to the substrate matters
                            if j_species != reaction.start_species:
                                continue
                            # get concentration of substrate
                            S = current_species_concentrations[region][n][reaction.start_species]
                            # if the reaction is spontaneous
                            if isinstance(reaction, SpontaneousReaction):
                                base = reaction.k # sign is placed afterwards
                            else: #EnzymaticReaction
                                base = (
                                    reaction.k_cat
                                    * ENZYME_CONCENTRATIONS[region][reaction.enzyme]
                                    * reaction.k_M
                                    / (reaction.k_M + S)**2
                                )
                            # Apply sign
                            if species == reaction.start_species:
                                J[i,j] += -base
                            elif species == reaction.end_species:
                                J[i,j] += +base



            else: # deal with left-most point within region (except r=0)
                (prev_region, prev_region_last_n), (_, _), (_, _) = NEIGHBORS[(region, n)]
                c_prev_region_last = current_species_concentrations[prev_region][prev_region_last_n][species]
                c_region_first = current_species_concentrations[region][0][species]
                c_region_second = current_species_concentrations[region][1][species]
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

        else: # point_type == "r"
            if region == NUM_REGIONS-1: # deal with r=R point
                (_, rR_neighbor_n), (_, rR_n) = NEIGHBORS[(region, n)]
                c_rR_neighbor = current_species_concentrations[region][rR_neighbor_n][species]
                c_rR = current_species_concentrations[region][rR_n][species]
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

            else: # deal with right-most point within region (except r=R)
                (_, _), (_, _), (next_region, _) = NEIGHBORS[(region, n)]
                c_second_to_last = current_species_concentrations[region][n-1][species]
                c_last = current_species_concentrations[region][n][species]
                c_next_region_first = current_species_concentrations[next_region][0][species]
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

    if fill_jacobian:
        return F, J
    else:
        return F, _

def define_newton_residual_and_optionally_jacobian_efficient(current_species_concentrations, fill_jacobian = True):
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
            # j is that by which we are applying the derivative to F_i
            # diffusion center
            j = POINT_IDS[(region, n, species)]
            J[i, j] += diff * (-2 / DELTA_R**2)

            # diffusion right
            j = POINT_IDS[(region, right_n, species)]
            J[i, j] += diff * (1/DELTA_R**2 + 1/(DELTA_R*r))

            # diffusion left
            j = POINT_IDS[(region, left_n, species)]
            J[i, j] += diff * (1/DELTA_R**2 - 1/(DELTA_R*r))

            # reactions
            for reaction in species.as_reactant_in + species.as_product_in:
                j = POINT_IDS[(region, n, reaction.start_species)]
                # get concentration of substrate
                S = current_species_concentrations[region][center_n][reaction.start_species]
                # if the reaction is spontaneous
                if isinstance(reaction, SpontaneousReaction):
                    base = reaction.k # sign is placed afterwards
                else: #EnzymaticReaction
                    base = (
                        reaction.k_cat
                        * ENZYME_CONCENTRATIONS[region][reaction.enzyme]
                        * reaction.k_M
                        / (reaction.k_M + S)**2
                    )
                # Apply sign
                if species == reaction.start_species:
                    J[i,j] += -base
                elif species == reaction.end_species:
                    J[i,j] += +base
        
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
                    if j_region == region and j_n == n:# if on the same place but not necessarily the same species
                        # we look at all reactions in which the species from i is involved
                        for reaction in species.as_reactant_in + species.as_product_in:
                            # but only the derivative with respect to the substrate matters
                            if j_species != reaction.start_species:
                                continue
                            # get concentration of substrate
                            S = current_species_concentrations[region][n][reaction.start_species]
                            # if the reaction is spontaneous
                            if isinstance(reaction, SpontaneousReaction):
                                base = reaction.k # sign is placed afterwards
                            else: #EnzymaticReaction
                                base = (
                                    reaction.k_cat
                                    * ENZYME_CONCENTRATIONS[region][reaction.enzyme]
                                    * reaction.k_M
                                    / (reaction.k_M + S)**2
                                )
                            # Apply sign
                            if species == reaction.start_species:
                                J[i,j] += -base
                            elif species == reaction.end_species:
                                J[i,j] += +base



            else: # deal with left-most point within region (except r=0)
                (prev_region, prev_region_last_n), (_, _), (_, _) = NEIGHBORS[(region, n)]
                c_prev_region_last = current_species_concentrations[prev_region][prev_region_last_n][species]
                c_region_first = current_species_concentrations[region][0][species]
                c_region_second = current_species_concentrations[region][1][species]
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

        else: # point_type == "r"
            if region == NUM_REGIONS-1: # deal with r=R point
                (_, rR_neighbor_n), (_, rR_n) = NEIGHBORS[(region, n)]
                c_rR_neighbor = current_species_concentrations[region][rR_neighbor_n][species]
                c_rR = current_species_concentrations[region][rR_n][species]
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

            else: # deal with right-most point within region (except r=R)
                (_, _), (_, _), (next_region, _) = NEIGHBORS[(region, n)]
                c_second_to_last = current_species_concentrations[region][n-1][species]
                c_last = current_species_concentrations[region][n][species]
                c_next_region_first = current_species_concentrations[next_region][0][species]
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
    tau_current,
    adaptive_step_parameters,
    last_F_norm
    ):
    # Unpack parameters
    tau_min = adaptive_step_parameters.get("tau_min")
    tau_max = adaptive_step_parameters.get("tau_max")
    tau_max = 1.0e-1 ####################################################################remove
    gamma_inc = adaptive_step_parameters.get("gamma_inc")
    gamma_dec = adaptive_step_parameters.get("gamma_dec")
    
    # From "An adaptive Newton-method based on a dynamical systems approach" paper
    F_vector, J_matrix = define_newton_residual_and_optionally_jacobian(species_concentrations)
    J_sparse = J_matrix.tocsc()
    NF = -spsolve(J_sparse, F_vector)
    norm_NF = np.linalg.norm(NF)
    t_n = min(np.sqrt(2 * tau_current / norm_NF), 1) # step size
    du = t_n * NF
    species_concentrations_try = copy.deepcopy(species_concentrations)
    for i, du_value in enumerate(du):
        (region, n, species) = REVERSE_POINT_IDS[i]
        species_concentrations_try[region][n][species] +=  du_value
        if species_concentrations_try[region][n][species] < 0:
            tau_new = tau_current * gamma_dec
            if tau_new < tau_min:
                raise ValueError("negative values!")
            return species_concentrations, tau_new, last_F_norm
    F_vector_new, _ = define_newton_residual_and_optionally_jacobian(species_concentrations_try, fill_jacobian=False)
    F_norm_new = np.linalg.norm(F_vector_new)

    if F_norm_new < last_F_norm:
        tau_new = min(tau_max, tau_current * gamma_inc)
        return species_concentrations_try, tau_new, F_norm_new

    else:
        tau_new = tau_current * gamma_dec
        if tau_new < tau_min:
           raise ValueError("Newton could not decrease the norm of the residual any more.")
        return species_concentrations, tau_new, F_norm_new

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
        point_type = POINT_INFOS[region][n]
        if point_type != "i":
            continue
        r = RADII[region][n]
        reaction_flux = calculate_reaction_term(current_species_concentrations, region, n, species)
        reaction_fluxes[species] += 4 * np.pi * reaction_flux * r**2 * DELTA_R
    # Second, calculate flux from boundary with exterior
    # the flux is positive if the concentration on the exterior is larger than on the interior at r=R
    boundary_fluxes = {species: 0
        for species in REACTION_NETWORK.species}
    for species in REACTION_NETWORK.species:
        #print(species.name, species.permeability_constant, species.external_concentration, current_species_concentrations[NUM_REGIONS-1][NUM_MESH_POINTS_IN_REGIONS[NUM_REGIONS-1]-1][species])
        boundary_fluxes[species] = 4 * np.pi * R**2 * species.permeability_constant * (
            species.external_concentration
            - current_species_concentrations[NUM_REGIONS-1][NUM_MESH_POINTS_IN_REGIONS[NUM_REGIONS-1]-1][species])
    # Since we are simulating the steady state, we want the total net flux to be 0 for each species
    # Because of numerics, we need some tolerance
    species_convergence = {species: None for species in REACTION_NETWORK.species}
    for species in REACTION_NETWORK.species:
        #print(species.name, reaction_fluxes[species], boundary_fluxes[species], flush=True)
        relative_deviation = abs(reaction_fluxes[species] + boundary_fluxes[species]) / max(abs(boundary_fluxes[species]), abs(reaction_fluxes[species]))
        if get_full_info:
            info[species] = relative_deviation
        if relative_deviation <= tol_relative_flux_deviation:
            species_convergence[species] = 0
        elif (1 - tol_relative_flux_deviation <= relative_deviation <= 1 + tol_relative_flux_deviation):
            species_convergence[species] = 1
        else:
            convergence = False
            if not get_full_info: # if we do not need the full info, early return
                return convergence, {}
    # In case all species have a relative deviation of around 1, this should not count
    # towards convergence, since this may occur at the beginning
    if all(species_convergence[species]==1 for species in REACTION_NETWORK.species):
        convergence = False
            
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
        log_convergence_info,
        convergence_info_logger_path,
        log_iteration_info_every,
        plot_iteration_data_during_simulation
    ):
    """
    save_data_every and check_convergence_every N iterations. If not to be done, set each to 0.
    """
    current_species_concentrations = initial_species_concentrations
    early_convergence = False # tracks whether the system fin
    convergence_info_logger = CSVLogger(convergence_info_logger_path)
    tau_current = convergence_parameters["initial_tau"]
    last_F_norm = 1e10
    for iter in tqdm(range(initial_iteration_number, int(max_num_newton_iterations)),
                     file=sys.stderr,
                     total=int(max_num_newton_iterations),
                     initial=initial_iteration_number):
        try:
            current_species_concentrations, tau_current, last_F_norm = adaptive_newton_step(
                current_species_concentrations,
                tau_current,
                adaptive_step_parameters,
                last_F_norm
                )
        except ValueError as e: # once the adaptive method cannot further decrease the norm of the residual, break
            if "Newton" in str(e):
                early_convergence = "Newton failed to decrease the norm any further"
            elif "negative" in str(e):
                early_convergence = "failure due to negative concentration"
            break

        # Save result if needed
        if save_data_every !=0 and (iter+1)%save_data_every==0:
            F_vector, J_matrix, du = compute_newton_step(current_species_concentrations)
            iter_string = str(iter).zfill(num_iterations_digits)
            save_newton_iteration_data(ITERATION_DATA_PATH, iter_string,
                J_matrix, F_vector, current_species_concentrations, du, variables_to_save_dictionary)
            if plot_iteration_data_during_simulation:
                fig, ax = plot_steady_state_concentrations(
                    reaction_network=REACTION_NETWORK,
                    num_regions=NUM_REGIONS,
                    num_mesh_points_in_regions=NUM_MESH_POINTS_IN_REGIONS,
                    radii=RADII,
                    membrane_radii=MEMBRANE_RADII,
                    output_file_name=os.path.join(ITERATION_DATA_PATH, f".iteration_nr_{iter_string}_concentrations.png"),
                    species_concentrations_to_plot=current_species_concentrations,
                    system_geometry_dict=SYSTEM_GEOMETRY_DICT["geometry_config"],
                    title = f"iteration #{iter}\n"
                )
                plt.close(fig)
        # Stop iterating if criterion for convergence fulfilled
        if check_convergence_every !=0 and iter%check_convergence_every==0:
            convergence_flux_equilibration, info_flux_equilibration = check_convergence_via_flux_equilibrium(
                current_species_concentrations,
                convergence_parameters,
                get_full_info=log_convergence_info
            )
            info_flux_equilibration.update({
                "runtime": f"{time.time() - simulation_start_time:.3f} seconds",
                "tau": float(tau_current),
                "F_vector_norm": float(last_F_norm)
            })
            if log_convergence_info:
                convergence_info_logger.log(iter, info_flux_equilibration)
            if convergence_flux_equilibration:
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
    PARAMETER_VALUE_CONDITIONS = read_yaml_file(os.path.join(FOLDER_TO_SOLVE, "parameters_value_conditions.yaml")) 
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
    print(NUM_POINTS)
    
    # Step 3: Put enzyme location information and assert that conditions are met
    try:
        REACTION_NETWORK = return_reaction_network_with_total_fixed_quantity_asserted(
            PARAMETER_VALUE_CONDITIONS["enzyme_total_fixed_quantity"],
            REACTION_NETWORK,
            PARAMETER_VALUE_CONDITIONS["enzyme_whose_quantity_to_modify_when_total_fixed_quantity"]
        )
        assert assert_no_conflicts_in_enzyme_positioning(
            REACTION_NETWORK,
            NUM_REGIONS,
            PARAMETER_VALUE_CONDITIONS["enzyme_impossible_combinations"]
        )
        ENZYME_CONCENTRATIONS, REACTION_NETWORK = return_enzyme_concentrations(
            REACTION_NETWORK,
            MEMBRANE_RADII,
            PARAMETER_VALUE_CONDITIONS["enzyme_maximum_concentration"]
        )
        pickle_dump_binary(
            os.path.join(FOLDER_TO_SOLVE, ".pickled_reaction_network_final"),
            REACTION_NETWORK
        )

    except Exception as e:
        # EARLY TERMINATION
        dump_json(FOLDER_TO_SOLVE,
                  ".species_steady_state_concentrations",
                  {"error": f"{e}"})
        sys.exit()


    # Step 4: Define structure that saves concentrations at each point and which
    # is updated with every iteration of Newton
    # currently, the concentration is initially set to increase linearly with the 
    # radius up until its external concentration
    max_external_concentration = max([species.external_concentration for species in REACTION_NETWORK.species])
    species_concentrations_guess = {
        region_idx : {
            mesh_point_idx : {
                species : max_external_concentration #* RADII[region_idx][mesh_point_idx] / RADII[NUM_REGIONS-1][NUM_MESH_POINTS_IN_REGIONS[NUM_REGIONS-1]-1]
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
        "tol_relative_flux_deviation": SOLVER_PARAMS["convergence_parameters"]["tol_relative_flux_deviation"],
        "initial_tau": 1#F_vector_guess*1.0e-5, # settting the initial value of tau close to F_vector_guess
    }

    # Important: save output for checking. tqdm is excluded
    print(f"Opening log file for {FOLDER_TO_SOLVE}.")
    with open(os.path.join(FOLDER_TO_SOLVE, ".newton_solver.log"), "a") as f, redirect_stdout(f):
        print(f"Starting solver from iteration number {initial_iteration_number} \n")
        #print("convergence parameters",
        #    {k: float(f"{v:.2e}") for k, v in convergence_parameters.items()}, "\n")
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
            log_iteration_info_every = SOLVER_INPUT["output_options"]["log_iteration_info_every"],
            log_convergence_info = SOLVER_INPUT["output_options"]["log_convergence_progress"],
            convergence_info_logger_path=os.path.join(FOLDER_TO_SOLVE, ".convergence_logger.csv"),
            plot_iteration_data_during_simulation = SOLVER_INPUT["output_options"]["plot_iteration_data_during_simulation"]
        )
        end_time = time.time()

        # Log data of final step
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
        species_concentrations_to_plot = species_concentrations_final,
        system_geometry_dict=SYSTEM_GEOMETRY_DICT["geometry_config"]
    )