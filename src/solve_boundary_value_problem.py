#%%
import sys
import os
import math
import numpy as np
from itertools import count
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_matrix, csr_matrix
import matplotlib.pyplot as plt
from auxiliary_functions_using_standard_library import pickle_load_binary, closest_value, dump_json
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme
from auxiliary_functions import save_matrix_as_sparse_txt

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

def return_saveable_species_concentrations_dict(original_species_concentrations_dict):
    """ Returns a dictionary equal to species_concentrations, but with species.name as keys
    to be able to save them as a .json file
    """
    species_concentrations_saveable_keys = {
        region_idx : {
            mesh_point_idx : {
                species.name : original_species_concentrations_dict[region_idx][mesh_point_idx][species]
                for species in REACTION_NETWORK.species}
            for mesh_point_idx in range(NUM_MESH_POINTS_IN_REGIONS[region_idx])}
        for region_idx in range(NUM_REGIONS)
    }
    return species_concentrations_saveable_keys

def calculate_reaction_term(region, n, species):
    """ Gives the reaction term for F_i.
    """
    reaction_term = 0
    for reaction in species.as_reactant_in + species.as_product_in:
        if isinstance(reaction, SpontaneousReaction):
            term = reaction.k * species_concentrations[region][n][reaction.start_species]
        else:
            term = reaction.k_cat * ENZYMES_CONCENTRATIONS[region][reaction.enzyme] * species_concentrations[region][n][reaction.start_species] / (reaction.k_M + species_concentrations[region][n][reaction.start_species])
        if reaction in species.as_reactant_in: # if acts as reactant, diminishes
            term *= -1
        reaction_term += term
    return reaction_term

def calculate_reaction_partial_derivative(reaction_to_check, partial_derivative_species, region, n):
    """ Gives the partial derivative of a reaction to a concentration of a species
    that is involved in the reaction.
    """
    if isinstance(reaction_to_check, SpontaneousReaction):
        derivative = reaction_to_check.k
    elif isinstance(reaction_to_check, EnzymaticReaction):
        derivative = reaction_to_check.k_cat * ENZYMES_CONCENTRATIONS[region][reaction_to_check.enzyme] * reaction_to_check.k_M / ( reaction_to_check.k_M + species_concentrations[region][n][partial_derivative_species])
    if partial_derivative_species == reaction_to_check.start_species:
        derivative *= -1
    return derivative

def solve_newton(max_newton_iterations, print_info=False):
    du_norm = np.inf
    for iter in range(max_newton_iterations):
        F = np.zeros(NUM_POINTS)
        J = lil_matrix((NUM_POINTS, NUM_POINTS))# np.zeros((NUM_POINTS, NUM_POINTS)) 
        for i in range(NUM_POINTS):
            (region, n, species) = REVERSE_POINT_IDS[i]
            r = RADII[region][n]
            D = species.diffusion_constant
            point_type = POINT_INFOS[region][n]
            # CONSTRUCT F_i
            # FOR EACH POINT WITHIN THE BULK
            if point_type == "i":
                (_, left_n), (_, center_n), (_, right_n) = NEIGHBORS[(region, n)]
                c_left = species_concentrations[region][left_n][species]
                c_center = species_concentrations[region][center_n][species]
                c_right = species_concentrations[region][right_n][species]
                diffusion_term = D * (1/ DELTA_R**2 * (c_right - 2* c_center + c_left) + 1 /(DELTA_R*r) * (c_right - c_left))
                reaction_term = calculate_reaction_term(region, center_n, species)
                F[i] = diffusion_term + reaction_term
                # FILL IN J_ij
                for j in range(NUM_POINTS):
                    (j_region, j_n, j_species) = REVERSE_POINT_IDS[j]
                    if j_region == region and j_n == n and j_species == species: # j == i, basically
                        J[i,j] += D * (1/DELTA_R**2 * (-2))
                    elif j_region==region and j==right_n and j_species == species: # same species, right or left 
                        J[i,j] += D * (1/DELTA_R**2 + 1/(DELTA_R*r))
                    elif j_region==region and j==left_n and j_species == species: # same species, right or left 
                        J[i,j] += D * (1/DELTA_R**2 - 1/(DELTA_R*r))
                    if j_region == region and j_n == center_n: # if on the same place but not necessarily the same species
                        for reaction in species.as_reactant_in + species.as_product_in:
                            if j_species in [reaction.start_species, reaction.end_species]:
                                J[i,j] += calculate_reaction_partial_derivative(reaction, j_species, region, center_n)
            elif point_type == "l":
                if region==0: # deal with r=0 point
                    (_, r0_n), (_, r0_neighbor_n) = NEIGHBORS[(region, n)]
                    c_r0 = species_concentrations[region][r0_n][species]
                    c_r0_neighbor = species_concentrations[region][r0_neighbor_n][species]
                    diffusion_term = 3 * D / DELTA_R**2 * 2 * (c_r0_neighbor - c_r0)
                    reaction_term = calculate_reaction_term(region, r0_n, species)
                    F[i] = diffusion_term + reaction_term
                    for j in range(NUM_POINTS):
                        (j_region, j_n, j_species) = REVERSE_POINT_IDS[j]
                        if j_region == region and j_n == 0 and j_species == species: # j == i, basically
                            J[i,j] += -3 * D / DELTA_R**2 * 2
                        elif j_region == region and j_n == 1 and j_species == species: # partial derivative to the one on the right
                            J[i,j] += 3 * D / DELTA_R**2 * 2
                        if j_region == region and j_n == n: # if on the same place but not necessarily the same species
                            for reaction in species.as_reactant_in + species.as_product_in:
                                if j_species in [reaction.start_species, reaction.end_species]:
                                    J[i,j] += calculate_reaction_partial_derivative(reaction, j_species, region, n)
                else: # deal with left-most point within region (except r=0)
                    (prev_region, prev_region_last_n), (_, _), (_, _) = NEIGHBORS[(region, n)]
                    c_prev_region_last = species_concentrations[prev_region][prev_region_last_n][species]
                    c_region_first = species_concentrations[region][0][species]
                    c_region_second = species_concentrations[region][1][species]
                    F[i] = D  * (c_region_second - c_region_first) / DELTA_R - species.permeability_constant * (c_region_first - c_prev_region_last)
                    for j in range(NUM_POINTS):
                        (j_region, j_n, j_species) = REVERSE_POINT_IDS[j]
                        if j_region == region and j_n == n and j_species == species:
                            J[i,j] = -D/DELTA_R - species.permeability_constant
                        elif j_region == region and j_species == species and j_n == 1:
                            J[i,j] = D/DELTA_R
                        elif j_region == prev_region and j_species == species and j_n == prev_region_last_n:
                            J[i,j] = -species.permeability_constant
                        
            else: # point_type == "r"
                if region == NUM_REGIONS-1: # deal with r=R point
                    (_, rR_neighbor_n), (_, rR_n) = NEIGHBORS[(region, n)]
                    c_rR_neighbor = species_concentrations[region][rR_neighbor_n][species]
                    c_rR = species_concentrations[region][rR_n][species]
                    F[i] = D * (c_rR - c_rR_neighbor) / DELTA_R - species.permeability_constant * (species.external_concentration - c_rR)
                    # CONSTRUCT J_ij
                    for j in range(NUM_POINTS):
                        (j_region, j_n, j_species) = REVERSE_POINT_IDS[j]
                        if j_region == region and j_n == n and j_species == species: # basically i=j
                            J[i,j] = D/DELTA_R + species.permeability_constant
                        elif j_region == region and j_species == species and j_n == rR_neighbor_n:
                            J[i,j] = -D/DELTA_R
                else: # deal with right-most point within region (except r=R)
                    (_, _), (_, _), (next_region, _) = NEIGHBORS[(region, n)]
                    c_second_to_last = species_concentrations[region][n-1][species]
                    c_last = species_concentrations[region][n][species]
                    c_next_region_first = species_concentrations[next_region][0][species]
                    F[i] = D  * (c_last - c_second_to_last) / DELTA_R - species.permeability_constant * (c_next_region_first - c_last)
                    for j in range(NUM_POINTS):
                        (j_region, j_n, j_species) = REVERSE_POINT_IDS[j]
                        if j_region == region and j_n == n and j_species == species: # basically i=j
                            J[i,j] = D/DELTA_R + species.permeability_constant
                        elif j_region == region and j_species == species and j_n == n-1:
                            J[i,j] = -D/DELTA_R
                        elif j_region == next_region and j_species == species and j_n == 0:
                            J[i,j] = species.permeability_constant
        # Newton update
        J_csr = J.tocsr()
        du = spsolve(J_csr, -F) # tocsr converts to CSR or CSC
        for i in range(len(du)):
            (region, n, species) = REVERSE_POINT_IDS[i]
            species_concentrations[region][n][species] += ALPHA * du[i]
        new_du_norm =  np.linalg.norm(du, np.inf)
        if print_info:
            if iter%1000==0:
                print(iter, new_du_norm)
                iter_string = str(iter).zfill(int(math.log10(max_newton_iterations)+1))
                plot_steady_state_concentrations(ITERATION_DATA_PATH, iter_string)
                #species_concentrations_saveable_keys = return_saveable_species_concentrations_dict(species_concentrations)
                #dump_json(ITERATION_DATA_PATH,
                #        f".iteration_nr_{iter_string}_concentration",
                #        species_concentrations_saveable_keys
                #)
                #np.savetxt(os.path.join(ITERATION_DATA_PATH, f".iteration_nr_{iter_string}_F.txt"), F, fmt="%.15e", delimiter="\n")
                #save_matrix_as_sparse_txt(J, os.path.join(ITERATION_DATA_PATH, f".iteration_nr_{iter_string}_J"))

        #if new_du_norm>du_norm:
        #    raise ValueError("The du norm increased!")
        du_norm = new_du_norm
        if np.linalg.norm(du, np.inf) < 1e-20:
            print(f"Converged in {iter+1} iterations.")
            break

def plot_steady_state_concentrations(folder, iter_string):
    x_values = []
    y_values = {}
    for species_idx, species in enumerate(REACTION_NETWORK.species):
        species_y_values = []
        for region in range(NUM_REGIONS):
            for n in range(NUM_MESH_POINTS_IN_REGIONS[region]):
                if species_idx == 0:
                    x_values.append(RADII[region][n])
                species_y_values.append(species_concentrations[region][n][species])
        y_values[species] = species_y_values

    fig, ax = plt.subplots(1,1, figsize = (5,3))
    for species in REACTION_NETWORK.species:
        ax.plot(x_values/max(x_values), y_values[species], label=species.name)
    ax.set_ylabel("concentration / M")
    ax.set_xlabel("relative distance to origin / r/R")
    ax.legend(
        loc='upper center',      # anchor point of legend
        bbox_to_anchor=(0.5, -0.25),  # (x, y) position in figure coordinates
        ncol=3,                  # number of columns
        frameon=False
    )
    for x_value in SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["internal_membrane_relative_radii"]:
        ax.axvline(x_value, linestyle = "--", alpha = 0.5, c = "k")
    max_value = max(max(y_values[species]) for species in REACTION_NETWORK.species)
    ax.set_ylim(ymin=0, ymax = max_value * 1.05)
    fig.savefig(os.path.join(folder, f".iteration_nr_{iter_string}_concentration.png"), dpi = 300, bbox_inches='tight')
    plt.close(fig)

if __name__ == "__main__":
    # Load all the information
    FOLDER_TO_SOLVE = sys.argv[1]
    REACTION_NETWORK = pickle_load_binary(os.path.join(FOLDER_TO_SOLVE, ".REACTION_NETWORK_pickle"))
    SYSTEM_GEOMETRY_DICT = pickle_load_binary(os.path.join(FOLDER_TO_SOLVE, ".SYSTEM_GEOMETRY_pickle"))

    # Step 1: Define all geometry variables
    R = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["outer_membrane_radius"]
    MESH_POINTS_IN_REGIONS = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["MESH_POINTS_IN_REGIONS"]
    NUM_MESH_POINTS_IN_REGIONS = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["NUM_MESH_POINTS_IN_REGIONS"]
    NUM_REGIONS = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["NUM_REGIONS"]

    # Step 2: Define structures to access geometry information
    POINT_IDS = build_point_ids_dict()
    REVERSE_POINT_IDS = build_reverse_point_ids_dict(POINT_IDS)
    RADII = build_radii_dict()
    DELTA_R = RADII[0][1]-RADII[0][0] # the different points within a region are equally spaced
    NUM_POINTS = len(REVERSE_POINT_IDS) # each point saves the concentration for one species at one node
    POINT_INFOS = build_point_infos_dict()
    NEIGHBORS = build_point_neighbor_dict()

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
    species_concentrations = {
        region_idx : {
            mesh_point_idx : {
                species : species.external_concentration * RADII[region_idx][mesh_point_idx] / RADII[NUM_REGIONS-1][NUM_MESH_POINTS_IN_REGIONS[NUM_REGIONS-1]-1]
                for species in REACTION_NETWORK.species}
            for mesh_point_idx in range(NUM_MESH_POINTS_IN_REGIONS[region_idx])}
        for region_idx in range(NUM_REGIONS)
    }
    # Step 5: Run solver and save result; plot
    ALPHA = 1

    ITERATION_DATA_PATH = os.path.join(FOLDER_TO_SOLVE, "solver_iteration_data")
    if not os.path.exists(ITERATION_DATA_PATH):
        os.makedirs(ITERATION_DATA_PATH)
    solve_newton(100000, True)
    # Modify species_concentrations such that we save the species through species.name
    species_concentrations_saveable_keys = {
        region_idx : {
            mesh_point_idx : {
                species.name : species_concentrations[region_idx][mesh_point_idx][species]
                for species in REACTION_NETWORK.species}
            for mesh_point_idx in range(NUM_MESH_POINTS_IN_REGIONS[region_idx])}
        for region_idx in range(NUM_REGIONS)
    }
    dump_json(FOLDER_TO_SOLVE, ".species_steady_state_concentrations", species_concentrations_saveable_keys)
    #plot_steady_state_concentrations()


