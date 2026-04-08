import sys
import os
import copy
import time
from tqdm import tqdm
import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_matrix
import matplotlib.pyplot as plt
from auxiliary_functions_using_standard_library import CSVLogger, find_max_in_nested_dict
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme
from plot_bvp_solution import plot_steady_state_concentrations
from auxiliary_functions_using_scipy import save_newton_iteration_data


def calculate_reaction_term(current_species_concentrations, region, n, species):
    """ Gives the reaction term for F_i (for a specific species at a specific point in the mesh).
    """
    reaction_term = 0
    for reaction in species.as_reactant_in + species.as_product_in:
        if isinstance(reaction, SpontaneousReaction):
            term = reaction.k * current_species_concentrations[region][n][reaction.start_species]
        else:
            term = reaction.k_cat * reaction.enzyme.regional_concentrations[region] * current_species_concentrations[region][n][reaction.start_species] / (reaction.k_M + current_species_concentrations[region][n][reaction.start_species])
        if reaction in species.as_reactant_in: # if acts as reactant, diminishes
            term *= -1
        reaction_term += term
    return reaction_term

def define_newton_residual_and_optionally_jacobian(
        current_species_concentrations,
        num_points,
        num_regions,
        radii,
        reverse_point_ids,
        point_infos,
        neighbors,
        Delta_r,
        fill_jacobian = True):
    """Defines the residual vector F and the jacobian matrix J (not sparse) 
    (the latter only if fill_jacobian is set to True (default)).
    Returns either F, _ or F, J.    
    """
    F = np.zeros(num_points)
    if fill_jacobian:
        J = lil_matrix((num_points, num_points))# np.zeros((num_points, num_points)) 
    for i in range(num_points):
        (region, n, species) = reverse_point_ids[i]
        r = radii[region][n]
        diff = species.diffusion_constant
        point_type = point_infos[region][n]
        # CONSTRUCT F_i
        # FOR EACH POINT WITHIN THE BULK
        if point_type == "i":
            (_, left_n), (_, center_n), (_, right_n) = neighbors[(region, n)]
            c_left = current_species_concentrations[region][left_n][species]
            c_center = current_species_concentrations[region][center_n][species]
            c_right = current_species_concentrations[region][right_n][species]
            diffusion_term = diff * (1/ Delta_r**2 * (c_right - 2* c_center + c_left) + 1 /(Delta_r*r) * (c_right - c_left))
            reaction_term = calculate_reaction_term(current_species_concentrations, region, center_n, species)
            F[i] = diffusion_term + reaction_term
            # FILL IN J_ij
            if not fill_jacobian:
                continue
            # j is that by which we are applying the derivative to F_i
            for j in range(num_points):
                (j_region, j_n, j_species) = reverse_point_ids[j]
                # Contributions from diffusion
                if j_region == region and j_n == n and j_species == species: # j == i, basically
                    J[i,j] += diff * (1/Delta_r**2 * (-2))
                elif j_region==region and j_n==right_n and j_species == species: # same species, right or left ###################
                    J[i,j] += diff * (1/Delta_r**2 + 1/(Delta_r*r))
                elif j_region==region and j_n==left_n and j_species == species: # same species, right or left #####################
                    J[i,j] += diff * (1/Delta_r**2 - 1/(Delta_r*r))
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
                                * reaction.enzyme.regional_concentrations[region]#ENZYME_CONCENTRATIONS[region][reaction.enzyme]
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
                (_, r0_n), (_, r0_neighbor_n) = neighbors[(region, n)]
                c_r0 = current_species_concentrations[region][r0_n][species]
                c_r0_neighbor = current_species_concentrations[region][r0_neighbor_n][species]
                diffusion_term = 3 * diff / Delta_r**2 * 2 * (c_r0_neighbor - c_r0)
                reaction_term = calculate_reaction_term(current_species_concentrations, region, r0_n, species)
                F[i] = diffusion_term + reaction_term
                if not fill_jacobian:
                    continue
                for j in range(num_points):
                    (j_region, j_n, j_species) = reverse_point_ids[j]
                    # Contributions from diffusion
                    if j_region == region and j_n == 0 and j_species == species: # j == i, basically
                        J[i,j] += -3 * diff / Delta_r**2 * 2
                    elif j_region == region and j_n == 1 and j_species == species: # partial derivative to the one on the right
                        J[i,j] += 3 * diff / Delta_r**2 * 2
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
                                    * reaction.enzyme.regional_concentrations[region]#[region][reaction.enzyme]
                                    * reaction.k_M
                                    / (reaction.k_M + S)**2
                                )
                            # Apply sign
                            if species == reaction.start_species:
                                J[i,j] += -base
                            elif species == reaction.end_species:
                                J[i,j] += +base



            else: # deal with left-most point within region (except r=0)
                (prev_region, prev_region_last_n), (_, _), (_, _) = neighbors[(region, n)]
                c_prev_region_last = current_species_concentrations[prev_region][prev_region_last_n][species]
                c_region_first = current_species_concentrations[region][0][species]
                c_region_second = current_species_concentrations[region][1][species]
                F[i] = diff  * (c_region_second - c_region_first) / Delta_r - species.permeability_constant * (c_region_first - c_prev_region_last)
                if not fill_jacobian:
                    continue
                for j in range(num_points):
                    (j_region, j_n, j_species) = reverse_point_ids[j]
                    # Contributions from diffusion
                    if j_region == region and j_n == n and j_species == species:
                        J[i,j] += -diff/Delta_r - species.permeability_constant
                    elif j_region == region and j_species == species and j_n == 1:
                        J[i,j] += diff/Delta_r
                    elif j_region == prev_region and j_species == species and j_n == prev_region_last_n:
                        J[i,j] += +species.permeability_constant ###########################SHOULD BE +###########################################I had a minus sign here.
                    # No contributions from reactions (flux considered)

        else: # point_type == "r"
            if region == num_regions-1: # deal with r=R point
                (_, rR_neighbor_n), (_, rR_n) = neighbors[(region, n)]
                c_rR_neighbor = current_species_concentrations[region][rR_neighbor_n][species]
                c_rR = current_species_concentrations[region][rR_n][species]
                F[i] = diff * (c_rR - c_rR_neighbor) / Delta_r - species.permeability_constant * (species.external_concentration - c_rR)
                if not fill_jacobian:
                    continue
                # CONSTRUCT J_ij
                for j in range(num_points):
                    (j_region, j_n, j_species) = reverse_point_ids[j]
                    if j_region == region and j_n == n and j_species == species: # basically i=j
                        J[i,j] += diff/Delta_r + species.permeability_constant
                    elif j_region == region and j_species == species and j_n == rR_neighbor_n:
                        J[i,j] += -diff/Delta_r

            else: # deal with right-most point within region (except r=R)
                (_, _), (_, _), (next_region, _) = neighbors[(region, n)]
                c_second_to_last = current_species_concentrations[region][n-1][species]
                c_last = current_species_concentrations[region][n][species]
                c_next_region_first = current_species_concentrations[next_region][0][species]
                F[i] = diff * (c_last - c_second_to_last) / Delta_r - species.permeability_constant * (c_next_region_first - c_last)
                if not fill_jacobian:
                    continue
                for j in range(num_points):
                    (j_region, j_n, j_species) = reverse_point_ids[j]
                    if j_region == region and j_n == n and j_species == species: # basically i=j
                        J[i,j] += diff/Delta_r + species.permeability_constant
                    elif j_region == region and j_species == species and j_n == n-1:
                        J[i,j] += -diff/Delta_r
                    elif j_region == next_region and j_species == species and j_n == 0:
                        J[i,j] += -species.permeability_constant ############################################## added the zero here

    if fill_jacobian:
        return F, J
    else:
        return F, _

def compute_newton_step(
        species_concentrations,
        num_points,
        num_regions,
        radii,
        reverse_point_ids,
        point_infos,
        neighbors,
        Delta_r
    ):
    """ Returns the vector of the residual, the jacobian as a sparse matrix, and delta u vector
    found by solving for the change in the concentrations.
    """
    # Step 1: Compute residual F and jacobian J
    F_vector, J_matrix = define_newton_residual_and_optionally_jacobian(
        species_concentrations,
        num_points=num_points,
        num_regions=num_regions,
        radii=radii,
        reverse_point_ids=reverse_point_ids,
        point_infos=point_infos,
        neighbors=neighbors,
        Delta_r=Delta_r
    )
    # Step 2: Assemble J (sparse), convert to solver-friendly format
    J_sparse = J_matrix.tocsc()
    # Step 3: Solve J du = -F
    du = spsolve(J_sparse, -F_vector)
    return F_vector, J_sparse, du

def adaptive_newton_step(
    current_iteration,
    species_concentrations,
    t_n,
    adaptive_step_parameters,
    last_F_norm,
    num_points,
    num_regions,
    radii,
    reverse_point_ids,
    point_infos,
    neighbors,
    Delta_r
    ):
    # Unpack parameters
    t_n_min = adaptive_step_parameters.get("t_n_min")
    gamma_inc = adaptive_step_parameters.get("gamma_inc")
    gamma_dec = adaptive_step_parameters.get("gamma_dec")
    
    # From "An adaptive Newton-method based on a dynamical systems approach" paper
    F_vector, J_sparse = define_newton_residual_and_optionally_jacobian(
        species_concentrations,
        num_points=num_points,
        num_regions=num_regions,
        radii=radii,
        reverse_point_ids=reverse_point_ids,
        point_infos=point_infos,
        neighbors=neighbors,
        Delta_r=Delta_r
    )
    J_matrix = J_sparse.toarray()
    #print(f"J[0,0] = {J_matrix[0,0]:.6e},  J[1,1] = {J_matrix[1,1]:.6e}", flush=True)
    kappa = np.linalg.cond(J_matrix) # calculate the condition number
    J_sparse = J_sparse.tocsc()
    NF = -spsolve(J_sparse, F_vector)
    du = t_n * NF
    species_concentrations_try = copy.deepcopy(species_concentrations)
    tolerance = kappa * np.finfo(F_vector.dtype).eps # absolute floor (this will never be crossed!)
    ###################################################################################################
    # The absolute floor cannot be crossed. We set 2 orders of magnitude higher
    #print(f"The tolerance is {tolerance}")
    ##########################################################
    # FOR CHECKING WHETHER IT IS ALWAYS THE SAME SPECIES THAT CAUSES TROUBLE
    ##########################################################
    #problematic_species = []
    #for i, du_value in enumerate(du):
    #    (region, n, species) = reverse_point_ids[i]
    #    if species_concentrations_try[region][n][species] +  du_value < 0:
    #        if species not in problematic_species:
    #            problematic_species.append(species)
    #if len(problematic_species) != 0:
    #    print("problematic_species: ", problematic_species, flush=True)
    ########################################
    ##########################################
    #max_concentration = find_max_in_nested_dict(species_concentrations)
    tol_abs = 0 #1e-4 * max_concentration #######################################################################
    relative_change_in_u = {i: None for i in range(len(du))}
    for i, du_value in enumerate(du):
        (region, n, species) = reverse_point_ids[i]
        u_old = species_concentrations[region][n][species] ############## instead of species_concentrations_try
        #if abs(u_old) > 0:
        #    relative_change_in_u[i] = abs(du_value / u_old)
        #else:
        #    relative_change_in_u[i] = abs(du_value)
        relative_change_in_u[i] = abs(du_value) / (tol_abs + tolerance * abs(u_old))
        species_concentrations_try[region][n][species] +=  du_value
        if species_concentrations_try[region][n][species] < 0:
            #print(species, flush=True)
            t_n_new = t_n * gamma_dec
            if t_n_new < t_n_min:
                raise ValueError("t_n is under the minimum.") # is only called
            #print("negative concentrations", flush=True)
            return species_concentrations, t_n_new, last_F_norm, t_n, kappa, np.nan
        #relative_change_in_u[i] = du_value / species_concentrations_try[region][n][species]
    #print("The max value is", max(relative_change_in_u.values()), flush=True)
    #if max(relative_change_in_u.values()) < 10 and current_iteration>10:
        #tolerance: #####################################################
        # AND: always run at least 100 iterations. (so that it does not immediately say it has converged.)
    #    raise ValueError("The numerical limit was found.")

    F_vector_new, _ = define_newton_residual_and_optionally_jacobian(
        species_concentrations_try,
        num_points=num_points,
        num_regions=num_regions,
        radii=radii,
        reverse_point_ids=reverse_point_ids,
        point_infos=point_infos,
        neighbors=neighbors,
        Delta_r=Delta_r,
        fill_jacobian=False)
    F_norm_new = np.linalg.norm(F_vector_new)

    if F_norm_new < last_F_norm:
        t_n_new = min(1, t_n * gamma_inc)
        #print(f"decreased norm, good! with step size t_n = {t_n}", flush=True)
        return species_concentrations_try, t_n_new, F_norm_new, kappa, max(relative_change_in_u.values())

    else:
        t_n_new = t_n * gamma_dec
        #print(f"did not decrease norm! with step size t_n = {t_n}", flush=True)
        if t_n_new < t_n_min:
           raise ValueError("t_n is under the minimum.")
        return species_concentrations, t_n_new, F_norm_new, kappa, max(relative_change_in_u.values()) #########################################

def get_info_flux_equilibrium(
        current_species_concentrations,
        reaction_network,
        num_points,
        reverse_point_ids,
        point_infos,
        radii,
        Delta_r,
        num_regions,
        R,
        num_mesh_points_in_regions
    ):
    """Returns a dictionary with information"""
    info_absolute = {species: 0 for species in reaction_network.species}
    info_relative = {species: 0 for species in reaction_network.species}

    # First, calculate net reaction fluxes within the sphere
    reaction_fluxes = {species: 0
        for species in reaction_network.species}
    for i in range(num_points):
        (region, n, species) = reverse_point_ids[i]
        point_type = point_infos[region][n]
        if point_type != "i":
            continue
        r = radii[region][n]
        reaction_flux = calculate_reaction_term(current_species_concentrations, region, n, species)
        reaction_fluxes[species] += 4 * np.pi * reaction_flux * r**2 * Delta_r
    
    # Second, calculate flux from boundary with exterior
    # the flux is positive if the concentration on the exterior is larger than on the interior at r=R
    boundary_fluxes = {species: 0
        for species in reaction_network.species}
    for species in reaction_network.species:
        #print(species.name, species.permeability_constant, species.external_concentration, current_species_concentrations[NUM_REGIONS-1][NUM_MESH_POINTS_IN_REGIONS[NUM_REGIONS-1]-1][species])
        boundary_fluxes[species] = 4 * np.pi * R**2 * species.permeability_constant * (
            species.external_concentration
            - current_species_concentrations[num_regions-1][num_mesh_points_in_regions[num_regions-1]-1][species])

    #  Since we are simulating the steady state, we want the total net flux to be 0 for each species
    # Because of numerics, we need some tolerance
    for species in reaction_network.species:
        info_absolute[species] = abs(reaction_fluxes[species] + boundary_fluxes[species])
        relative_deviation = info_absolute[species] / max(abs(boundary_fluxes[species]), abs(reaction_fluxes[species]))
        info_relative[species] = relative_deviation

    info_absolute = {f"{species}_absolute": value for species, value in info_absolute.items()}
    info_relative = {f"{species}_relative": value for species, value in info_relative.items()}     
    info = info_absolute | info_relative # merge the two dictionaries

    return info


def solve_newton(
        # timing
        simulation_start_time,
        # system information
        reaction_network,
        system_geometry_dict,
        expanded_system_mesh_dict,
        adaptive_step_parameters,
        max_num_newton_iterations,
        # simulation initial values
        initial_iteration_number,
        initial_species_concentrations,
        initial_residual_norm,
        initial_runtime,
        # saving
        iteration_data_path,
        progress_log_path,
        filename_for_newton_function,
        variables_to_save_dictionary,
        save_data_every,
        log_progress_every,
        plot_iteration_data_during_simulation
    ):
    """
    save_data_every and check_convergence_every N iterations. If not to be done, set each to 0.
    """
    # Alias for information in system_geometry_dict
    R = system_geometry_dict["geometry_config"]["outer_membrane_radius"]
    mesh_points_in_regions = system_geometry_dict["geometry_config"]["mesh_points_in_regions"]
    num_mesh_points_in_regions = system_geometry_dict["geometry_config"]["num_mesh_points_in_regions"]
    num_regions = system_geometry_dict["geometry_config"]["num_regions"]
    membrane_radii = system_geometry_dict["geometry_config"]["membrane_radii"]

    # Alias for information in expanded_system_mesh_dict
    radii = expanded_system_mesh_dict["radii"]
    Delta_r = expanded_system_mesh_dict["delta_r"]
    num_points = expanded_system_mesh_dict["num_points"]
    point_infos = expanded_system_mesh_dict["point_infos"]
    reverse_point_ids = expanded_system_mesh_dict["reverse_point_ids"]
    neighbors = expanded_system_mesh_dict["neighbors"]

    # Start 
    current_species_concentrations = initial_species_concentrations
    progress_logger = CSVLogger(progress_log_path)
    t_n = 1
    last_F_norm = initial_residual_norm
    for iter in tqdm(range(initial_iteration_number, int(max_num_newton_iterations)),
                     file=sys.stderr,
                     total=int(max_num_newton_iterations),
                     initial=initial_iteration_number):
        try:
            current_species_concentrations, t_n, last_F_norm, kappa, max_rel_concentration_change = adaptive_newton_step(
                iter,
                current_species_concentrations,
                t_n,
                adaptive_step_parameters,
                last_F_norm,
                num_points=num_points,
                num_regions=num_regions,
                radii=radii,
                reverse_point_ids=reverse_point_ids,
                point_infos=point_infos,
                neighbors=neighbors,
                Delta_r=Delta_r
                )
        except ValueError as e: # once the adaptive method cannot further decrease the norm of the residual, break
            if "The numerical limit was found." in str(e):
                print(f"The numerical limit was found in iteration {iter}", flush=True)
            elif "t_n is under the minimum" in str(e):
                print(f"t_n reached the minimum in iteration {iter}", flush=True)
            break

        # Save result if needed
        if save_data_every !=0 and (iter+1)%save_data_every==0:
            F_vector, J_matrix, du = compute_newton_step(
                current_species_concentrations,
                num_points=num_points,
                num_regions=num_regions,
                radii=radii,
                reverse_point_ids=reverse_point_ids,
                point_infos=point_infos,
                neighbors=neighbors,
                Delta_r=Delta_r
            )
            save_newton_iteration_data(iteration_data_path, filename_for_newton_function(iter),
                J_matrix, F_vector, current_species_concentrations, du, variables_to_save_dictionary)
            if plot_iteration_data_during_simulation:
                fig, _ = plot_steady_state_concentrations(
                    reaction_network=reaction_network,
                    num_regions=num_regions,
                    num_mesh_points_in_regions=num_mesh_points_in_regions,
                    radii=radii,
                    membrane_radii=membrane_radii,
                    output_file_name=os.path.join(iteration_data_path, f".iteration_nr_{str(iter).zfill(6)}_concentrations.png"),
                    species_concentrations_to_plot=current_species_concentrations,
                    system_geometry_dict=system_geometry_dict["geometry_config"],
                    title = f"iteration #{iter}\n"
                )
                plt.close(fig)
        # Stop iterating if criterion for convergence fulfilled
        if log_progress_every !=0 and iter%log_progress_every==0:
            info_flux_equilibration = get_info_flux_equilibrium(
                current_species_concentrations,
                reaction_network=reaction_network,
                num_points=num_points,
                reverse_point_ids=reverse_point_ids,
                point_infos=point_infos,
                radii=radii,
                Delta_r=Delta_r,
                num_regions=num_regions,
                R=R,
                num_mesh_points_in_regions=num_mesh_points_in_regions
                )
            info_flux_equilibration.update({
                "runtime": f"{initial_runtime + time.time() - simulation_start_time:.3f} seconds",
                "F_vector_norm": float(last_F_norm),
                "t_n": float(t_n),
                "condition number": float(kappa),
                "max_rel_concentration_change": float(max_rel_concentration_change)
            })
            progress_logger.log(iter, info_flux_equilibration)

    return current_species_concentrations



"""
def define_newton_residual_and_optionally_jacobian_efficient(current_species_concentrations, fill_jacobian = True):
    #Defines the residual vector F and the jacobian matrix J (not sparse) 
    #(the latter only if fill_jacobian is set to True (default)).
    #Returns either F, _ or F, J.    
    #
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

"""