#%%
import sys
import os
import numpy as np
import pprint
from scipy import integrate
from auxiliary_functions_using_standard_library import pickle_load_binary, closest_value, dump_json
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme

# Get the absolute path to this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Go up to the project root (one or more levels as needed)
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
# Access other folders relative to the root
DATA_PATH = os.path.join(PROJECT_ROOT, "data")

case_folder = os.path.join(
    DATA_PATH, "minimal_test_0"
)

folder_to_solve = case_folder
reaction_network = pickle_load_binary(os.path.join(folder_to_solve, ".REACTION_NETWORK_pickle"))
system_geometry_dict = pickle_load_binary(os.path.join(folder_to_solve, ".SYSTEM_GEOMETRY_pickle"))

# Step 1: Define all geometry variables
R = system_geometry_dict["GEOMETRY_CONFIG"]["outer_membrane_radius"]
MESH_POINTS_IN_REGIONS = system_geometry_dict["GEOMETRY_CONFIG"]["MESH_POINTS_IN_REGIONS"]
NUM_MESH_POINTS_IN_REGIONS = system_geometry_dict["GEOMETRY_CONFIG"]["NUM_MESH_POINTS_IN_REGIONS"]
NUM_REGIONS = system_geometry_dict["GEOMETRY_CONFIG"]["NUM_REGIONS"]

species_concentrations = {
    region_idx : {
        mesh_point_idx : {
            species : 0
            for species in reaction_network.species}
        for mesh_point_idx in range(NUM_MESH_POINTS_IN_REGIONS[region_idx])}
    for region_idx in range(NUM_REGIONS)
}

enzymes_concentrations = {
    region_idx : {
        enzyme : enzyme.concentration if region_idx in enzyme.regions else 0
        for enzyme in reaction_network.enzymes
    }
    for region_idx in range(NUM_REGIONS)
}

from itertools import count
counter = count()
point_ids = {
    region_idx : {
        mesh_point_idx : {
            species : next(counter)
            for species in reaction_network.species}
        for mesh_point_idx in range(NUM_MESH_POINTS_IN_REGIONS[region_idx])}
    for region_idx in range(NUM_REGIONS)
}

reverse_point_ids = {
    value: (region_idx, mesh_point_idx, species)
    for region_idx, mesh_points in point_ids.items()
    for mesh_point_idx, species_map in mesh_points.items()
    for species, value in species_map.items()
}

radii = {
    region_idx : {
        mesh_point_idx : MESH_POINTS_IN_REGIONS[region_idx][mesh_point_idx]
        for mesh_point_idx in range(NUM_MESH_POINTS_IN_REGIONS[region_idx])
    }
    for region_idx in range(NUM_REGIONS)
}
DELTA_R = radii[0][1]-radii[0][0]

NUM_POINTS = len(reverse_point_ids)

point_infos = {
    region_idx : {
        mesh_point_idx : "l" if mesh_point_idx==0 else ("r" if mesh_point_idx==NUM_MESH_POINTS_IN_REGIONS[region_idx]-1 else "i")
        for mesh_point_idx in range(NUM_MESH_POINTS_IN_REGIONS[region_idx])
    }
    for region_idx in range(NUM_REGIONS)
}

neighbors = {}

for region_idx, mesh_points in point_infos.items():
    for n, kind in mesh_points.items():
        if kind == "i":
            # Internal: previous, self, next
            neighbors[(region_idx, n)] = [
                (region_idx, n - 1),
                (region_idx, n),
                (region_idx, n + 1),
            ]
        elif kind == "l":
            # Left boundary: connect to previous region's rightmost point (if it exists)
            if region_idx > 0:
                prev_region_last = NUM_MESH_POINTS_IN_REGIONS[region_idx - 1] - 1
                neighbors[(region_idx, n)] = [
                    (region_idx - 1, prev_region_last),
                ]
            else:
                 neighbors[(region_idx, n)] = []
            neighbors[(region_idx, n)].append( (region_idx, 0) )
            neighbors[(region_idx, n)].append( (region_idx, 1) )
        elif kind == "r":
            # Right boundary: self’s region last, then 0 and 1 in same region
            last_in_region = NUM_MESH_POINTS_IN_REGIONS[region_idx] - 1
            neighbors[(region_idx, n)] = [
                (region_idx, last_in_region-1),
                (region_idx, last_in_region),
            ]
            if region_idx < NUM_REGIONS-1:
                neighbors[(region_idx, n)].append( (region_idx+1, 0) )

def calculate_reaction_term(region, n, species):
    reaction_term = 0
    for reaction in species.as_reactant_in + species.as_product_in:
        if isinstance(reaction, SpontaneousReaction):
            term = reaction.k * species_concentrations[region][n][reaction.start_species]
        elif isinstance(reaction, EnzymaticReaction):
            term = reaction.k_cat * enzymes_concentrations[region][reaction.enzyme] * species_concentrations[region][n][reaction.start_species] / (reaction.k_M + species_concentrations[region][n][reaction.start_species])
        else:
            raise ValueError("somehow", reaction, "is neither Enzymatic nor Spontaneous...")
        if reaction in species.as_reactant_in: # if acts as reactant, diminishes
            term *= -1
        reaction_term += term
    return reaction_term

def calculate_reaction_partial_derivative(reaction_to_check, partial_derivative_species, region, n):
    if isinstance(reaction_to_check, SpontaneousReaction):
        derivative = reaction_to_check.k
    elif isinstance(reaction_to_check, EnzymaticReaction):
        derivative = reaction_to_check.k_cat * enzymes_concentrations[region][reaction_to_check.enzyme] * reaction_to_check.k_M / ( reaction_to_check.k_M + species_concentrations[region][n][partial_derivative_species])
    if partial_derivative_species == reaction_to_check.start_species:
        derivative *= -1
    return derivative


from scipy.sparse.linalg import spsolve
max_newton_iterations = 30
for iter in range(max_newton_iterations):
    F = np.zeros(NUM_POINTS)
    J = np.zeros((NUM_POINTS, NUM_POINTS))
    for i in range(NUM_POINTS):
        (region, n, species) = reverse_point_ids[i]
        r = radii[region][n]
        D = species.diffusion_constant
        point_type = point_infos[region][n]
        # CONSTRUCT F_i
        # FOR EACH POINT WITHIN THE BULK
        if point_type == "i": 
            (_, left_n), (_, center_n), (_, right_n) = neighbors[(region, n)]
            c_left = species_concentrations[region][left_n][species]
            c_center = species_concentrations[region][center_n][species]
            c_right = species_concentrations[region][right_n][species]
            diffusion_term = r**2 / DELTA_R**2 * (c_right - 2* c_center + c_left) + r /DELTA_R * (c_right - c_left)
            reaction_term = calculate_reaction_term(region, center_n, species)
            F[i] = diffusion_term + r**2/D * reaction_term
            # FILL IN J_ij
            for j in range(NUM_POINTS):
                (j_region, j_n, j_species) = reverse_point_ids[j]
                if j_region == region and j_n == n and j_species == species: # j == i, basically
                    J[i][j] += r**2 / DELTA_R**2 * (-2)
                elif j_region==region and (j_n==left_n or j==right_n) and j_species == species: # same species, right or left 
                    J[i][j] += r**2 / DELTA_R**2 + r/DELTA_R
                if j_region == region and j_n == center_n: # if on the same place but not necessarily the same species
                    for reaction in species.as_reactant_in + species.as_product_in:
                        if j_species in [reaction.start_species, reaction.end_species]:
                            J[i][j] += calculate_reaction_partial_derivative(reaction, j_species, region, center_n)
        elif point_type == "l":
            if region==0: # deal with r=0 point
                (_, r0_n), (_, r0_neighbor_n) = neighbors[(region, n)]
                c_r0 = species_concentrations[region][r0_n][species]
                c_r0_neighbor = species_concentrations[region][r0_neighbor_n][species]
                diffusion_term = 3 * D / DELTA_R**2 * (c_r0_neighbor - c_r0)
                reaction_term = calculate_reaction_term(region, r0_n, species)
                F[i] = diffusion_term + reaction_term
                for j in range(NUM_POINTS):
                    (j_region, j_n, j_species) = reverse_point_ids[j]
                    if j_region == region and j_n == n and j_species == species: # j == i, basically
                        J[i][j] += -3 * D / DELTA_R**2 * 2
                    elif j_region == region and j_n == n+1 and j_species == species: # partial derivative to the one on the right
                        J[i][j] += 3 * D / DELTA_R**2 * 2
                    if j_region == region and j_n == n: # if on the same place but not necessarily the same species
                        for reaction in species.as_reactant_in + species.as_product_in:
                            if j_species in [reaction.start_species, reaction.end_species]:
                                J[i][j] += calculate_reaction_partial_derivative(reaction, j_species, region, n)
            else: # deal with left-most point within region (except r=0)
                (prev_region, prev_region_last_n), (_, _), (_, _) = neighbors[(region, n)]
                c_prev_region_last = species_concentrations[prev_region][prev_region_last_n][species]
                c_region_first = species_concentrations[region][0][species]
                c_region_second = species_concentrations[region][1][species]
                F[i] = D  * (c_region_second - c_region_first) / DELTA_R - species.permeability_constant * (c_region_first - c_prev_region_last)
                for j in range(NUM_POINTS):
                    (j_region, j_n, j_species) = reverse_point_ids[j]
                    if j_region == region and j_n == n and j_species == species:
                        J[i][j] = -D/DELTA_R - species.permeability_constant
                    elif j_region == region and j_species == species and j_n == 1:
                        J[i][j] = D/DELTA_R
                    elif j_region == prev_region and j_species == species and j_n == prev_region_last_n:
                        J[i][j] = -species.permeability_constant
                    
        else: # point_type == "r"
            if region == NUM_REGIONS-1: # deal with r=R point
                (_, rR_neighbor_n), (_, rR_n) = neighbors[(region, n)]
                c_rR_neighbor = species_concentrations[region][rR_neighbor_n][species]
                c_rR = species_concentrations[region][rR_n][species]
                F[i] = D * (c_rR - c_rR_neighbor) / DELTA_R - species.permeability_constant * (species.external_concentration - c_rR)
                # CONSTRUCT J_ij
                for j in range(NUM_POINTS):
                    (j_region, j_n, j_species) = reverse_point_ids[j]
                    if j_region == region and j_n == n and j_species == species: # basically i=j
                        J[i][j] = D/DELTA_R + species.permeability_constant
                    elif j_region == region and j_species == species and j_n == rR_neighbor_n:
                        J[i][j] = -D/DELTA_R
            else: # deal with right-most point within region (except r=R)
                (_, _), (_, _), (next_region, _) = neighbors[(region, n)]
                c_second_to_last = species_concentrations[region][n-1][species]
                c_last = species_concentrations[region][n][species]
                c_next_region_first = species_concentrations[next_region][0][species]
                F[i] = D  * (c_last - c_second_to_last) / DELTA_R - species.permeability_constant * (c_next_region_first - c_last)
                for j in range(NUM_POINTS):
                    (j_region, j_n, j_species) = reverse_point_ids[j]
                    if j_region == region and j_n == n and j_species == species: # basically i=j
                        J[i][j] = D/DELTA_R + species.permeability_constant
                    elif j_region == region and j_species == species and j_n == n-1:
                        J[i][j] = -D/DELTA_R
                    elif j_region == next_region and j_species == species and j_n == 0:
                        J[i][j] = species.permeability_constant
    # Newton update
    du = spsolve(J, -F) # tocsr converts to CSR or CSC
    for i in range(len(du)):
        (region, n, species) = reverse_point_ids[i]
        species_concentrations[region][n][species] += du[i]

    # Convergence check
    if np.linalg.norm(du, np.inf) < 1e-8:
       print(f"Converged in {iter+1} iterations.")
       break

dump_json(folder_to_solve, "species_steady_state_concentrations", species_concentrations)
dump_json(folder_to_solve, "enzymes_concentrations", enzymes_concentrations)
dump_json(folder_to_solve, "mesh_radii", radii)

#%%
system_geometry_dict
#%%
import matplotlib.pyplot as plt
x_values = []
y_values = {}
for species_idx, species in enumerate(reaction_network.species):
    species_y_values = []
    for region in range(NUM_REGIONS):
        for n in range(NUM_MESH_POINTS_IN_REGIONS[region]):
            if species_idx == 0:
                x_values.append(radii[region][n])
            species_y_values.append(species_concentrations[region][n][species])
    y_values[species] = species_y_values

fig, ax = plt.subplots(1,1, figsize = (5,3))
for species in reaction_network.species:    
    ax.plot(x_values/max(x_values), y_values[species], label=species.name)
ax.set_ylabel("concentration / M")
ax.set_xlabel("relative distance to origin / r/R")
ax.legend(
    loc='upper center',      # anchor point of legend
    bbox_to_anchor=(0.5, -0.25),  # (x, y) position in figure coordinates
    ncol=3,                  # number of columns
    frameon=False
)
for x_value in system_geometry_dict["GEOMETRY_CONFIG"]["internal_membrane_relative_radii"]:
    ax.axvline(x_value, linestyle = "--", alpha = 0.5, c = "k")
fig.savefig(os.path.join(folder_to_solve, "solution.png"), dpi = 300, bbox_inches='tight')






#%%
def create_bvp_result_plots(bvp_result, bvp_variable_indices, reaction_network):
    
    x_values = np.linspace(0, max(bvp_result.x), num = 100)
    integrals = {}
    num_enzymes = len(reaction_network.enzymes)

    enzyme_max_concentration = max(enzyme.concentration for enzyme in reaction_network.enzymes)
    norm = mcolors.Normalize(vmin=0, vmax=enzyme_max_concentration)
    cmap = matplotlib.colormaps['Oranges']
    scalar_map = cm.ScalarMappable(norm=norm, cmap=cmap)

    # MAIN PLOT
    fig, ax = plt.subplots(num_enzymes + 1,1, figsize = (5,3), sharex = True,
                        gridspec_kw={'height_ratios': [1]*num_enzymes + [num_enzymes*2]})
    fig.subplots_adjust(hspace=0)
    ax[0].set_title(f"steady state distribution, computed with {len(bvp_result.x)} nodes")
    # Main plot enzyme concentration
    for enzyme_idx, enzyme in enumerate(reaction_network.enzymes):
        ax[enzyme_idx].set_ylabel(enzyme.name, rotation=0, labelpad=20)
        # Get the label Text object and adjust its position
        label = ax[enzyme_idx].yaxis.get_label()
        # Move it vertically to center (around 0.5 in axes coords) and keep horizontal offset from labelpad
        label.set_verticalalignment('center')  # ensure vertical alignment
        # Manually set position: x controls horizontal offset, y controls vertical position (0 bottom, 0.5 center, 1 top)
        label.set_position((label.get_position()[0], 0.5))
        # Colorbar color
        color = scalar_map.to_rgba(enzyme.concentration)
        # Fill x-ranges where enzymes are at
        for localizationTuple in enzyme.localization:
            min_range = localizationTuple.minMaxLoc[0]
            max_range = localizationTuple.minMaxLoc[1]
            ax[enzyme_idx].fill_between(
                x_values/max(x_values), 0, 1, where=(
                    x_values/max(x_values) >= min_range)
                    & (x_values/max(x_values) <= max_range),
                color=color)
        ax[enzyme_idx].set_yticks([])
    
    # Main plot species concentration
    for species in reaction_network.species:
        # bvp_result.sol is the found solution for y as scipy.interpolate.PPoly instance, a C1 continuous
        # cubic spline.
        integrand = lambda r: bvp_result.sol(r)[bvp_variable_indices[species]["prim"]]* 10**3 * 4 * np.pi * r**2
        integral, _ = integrate.quad(integrand, 0, max(bvp_result.x))
        integrals[species.name] = integral
    
    # Find smallest exponent (in base 10) across all values
    if integrals:
        min_exp = min(int(np.floor(np.log10(abs(val)))) if val != 0 else 0 for val in integrals.values())
    else:
        min_exp = 0  # fallback

    # Format legend labels
    labels_with_integrals = {}
    for species in reaction_network.species:
        val = integrals[species.name]
        scaled_val = val / (10 ** min_exp)  # Adjust to common exponent
        label = f"{species.name}: {scaled_val:.0f}e{min_exp} mol"
        labels_with_integrals[species.name] = label
    
    # Plot and assign labels
    for species in reaction_network.species:    
        solution_y_values_to_plot = bvp_result.sol(x_values)[bvp_variable_indices[species]["prim"]]
        ax[-1].plot(x_values/max(x_values), solution_y_values_to_plot, label=labels_with_integrals[species.name])
    
    # Main plot mesh points visualization
    for node_x in bvp_result.x:
        ax[-1].axvline(node_x/max(x_values), ymin = 0.95, ymax = 1, c = "k", linewidth = 1)
    # Main plot labels
    ax[-1].set_xlabel("r/R")
    ax[-1].set_xlim([0, 1])
    ax[-1].set_ylabel("concentration / M")

    # ENZYME CONCENTRATION COLORBAR
    # Assume scalar_map is your ScalarMappable (from cmap and norm)
    fig_colorbar = plt.figure(figsize=(2, 4))  # tall figure for vertical bar
    # Add axes for colorbar (full figure area)
    cbar_ax = fig_colorbar.add_axes([0.2, 0.05, 0.3, 0.9])  # [left, bottom, width, height]
    # Create the colorbar
    cbar = fig_colorbar.colorbar(scalar_map, cax=cbar_ax, orientation='vertical')
    cbar.set_label("enzyme concentration")
    # SPECIES CONCENTRATION LEGEND
    fig_legend = plt.figure(figsize=(3, 2))
    ax_legend = fig_legend.add_subplot(111)
    ax_legend.axis('off')
    handles, labels = ax[-1].get_legend_handles_labels()
    ax_legend.legend(handles, labels, loc='center')

    return fig, fig_colorbar, fig_legend     



#%%
if __name__ == "__main__":
    # Load all the information
    folder_to_solve = sys.argv[1]
    reaction_network = pickle_load_binary(os.path.join(folder_to_solve, ".REACTION_NETWORK_pickle"))
    system_geometry_dict = pickle_load_binary(os.path.join(folder_to_solve, ".SYSTEM_GEOMETRY_pickle"))

    # Step 1: Define all geometry variables
    membrane_radii = system_geometry_dict["GEOMETRY_CONFIG"]["membrane_radii"]
    num_mesh_points = system_geometry_dict["GEOMETRY_CONFIG"]["num_mesh_points"]
    initial_r_mesh = np.linspace(min(membrane_radii)*1e-3, max(membrane_radii), num = num_mesh_points) #nm # doesn't work within singularity

    # Step 2: Do preliminary work to make reaction diffusion system solvable
    """We need to define 2*n variables, where n is the number of species, since we
    have to convert the system with second derivatives to one with first derivatives
    with double the variables
    """
    species_variable_indices = { # each variable is mapped to an integer
        species: {"prim": species_idx*2,
                "1stDer": species_idx*2 + 1}
        for species_idx, species in enumerate(reaction_network.species)
    }
    # Set some initial values within the nodes for each variable 
    species_variables_inital_values_dict = {
        species: {key: np.zeros(num_mesh_points) for key in inner_dict}
        for species, inner_dict in species_variable_indices.items()
    }
    # Initial guesses can be modified like this, if needed
    species_variables_inital_values_dict[reaction_network.species["Trp"]]["prim"] = np.zeros(num_mesh_points)
    # Create the 2d array map
    initial_values_2d_array = np.zeros((2*len(reaction_network.species), num_mesh_points))
    for species in reaction_network.species:
        for variable_type in species_variable_indices[species].keys():
            # writes the species_variable_indices[species][variable_type]th row of the 2d array
            initial_values_2d_array[species_variable_indices[species][variable_type]] = species_variables_inital_values_dict[species][variable_type]

    # Step 3: run solver
    print(f"Solving bvp in {folder_to_solve}")
    bvp_result = solve_bvp(reaction_diffusion_system, boundary_conditions,
                    initial_r_mesh, initial_values_2d_array,
                    max_nodes=1000)
    print(f"Solved bvp in {folder_to_solve}")
    
    # Step 4: pickle
    pickle_dump_binary(
        os.path.join(folder_to_solve, ".BOUNDARY_VALUE_PROBLEM_RESULT_pickle"),
        bvp_result
    )
    pickle_dump_binary(
        os.path.join(folder_to_solve, ".BOUNDARY_VALUE_PROBLEM_VARIABLE_INDICES_pickle"),
        species_variable_indices
    )

#%%
reaction_network.enzymes["VioA"]
# %%
