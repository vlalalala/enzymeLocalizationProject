#%%
import os
import itertools
import numpy as np
from scipy.integrate import solve_bvp
from auxiliary_functions_using_standard_library import load_json, pickle_load_binary
from auxiliary_functions import LocationTuple
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme

#%%
# Get the absolute path to this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Go up to the project root (one or more levels as needed)
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))

violacein_folder = os.path.join(PROJECT_ROOT, *["data", "violacein_0"])

#%%
reaction_network = pickle_load_binary(os.path.join(violacein_folder, ".REACTION_NETWORK_pickle"))
system_geometry_dict = pickle_load_binary(os.path.join(violacein_folder, ".SYSTEM_GEOMETRY_pickle"))

#%% Step 1: calculate total volume within which enzyme is found
radius = system_geometry_dict["GEOMETRY_CONFIG"]["radius"]
num_mesh_points = system_geometry_dict["GEOMETRY_CONFIG"]["num_mesh_points"]

r_mesh = np.linspace(1e-6, radius, num = num_mesh_points) #nm # doesn't work within singularity

#%%
# We need to define 2*n variables, where n is the number of species, since we
# have to convert the system with second derivatives to one with first derivatives
# with double the variables
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

# modify initial guesses like this
species_variables_inital_values_dict[reaction_network.species["Trp"]]["prim"] = np.zeros(num_mesh_points)

# Create the 2d array map
initial_values_2d_array = np.zeros((2*len(reaction_network.species), num_mesh_points))
for species in reaction_network.species:
    for variable_type in species_variable_indices[species].keys():
        # writes the species_variable_indices[species][variable_type]th row of the 2d array
        initial_values_2d_array[species_variable_indices[species][variable_type]] = species_variables_inital_values_dict[species][variable_type]

#%% Define function to get enzyme concentration at each node
def enzyme_concentration(r_mesh):
    enzyme_conc = {
        enzyme: []
        for enzyme in reaction_network.enzymes
    }
    # Get the enzyme concentration for each mesh node PUT THIS IN A FUNCTION THAT TAKES R_MESH
    for enzyme in reaction_network.enzymes:
        enzyme_mesh_occupied_bool_dict = {r: False for r in r_mesh} # Initialize whether mesh points have enzyme
        total_occupied_volume = 0
        for locTuple in enzyme.localization:
            volume_part = 4/3 * np.pi * ((locTuple.minMaxLoc[1] * radius)**3 - (locTuple.minMaxLoc[0] * radius)**3)
            total_occupied_volume += volume_part
            for r_key in enzyme_mesh_occupied_bool_dict.keys():
                if locTuple.return_within_tuple(r_key):
                    enzyme_mesh_occupied_bool_dict[r_key] = True
        enzyme_concentration = [
            enzyme.quantity/total_occupied_volume 
            if enzyme_mesh_occupied_bool_dict[r] == True else 0
            for r in r_mesh
        ]################################ not really correct... but will work for now
        enzyme_conc[enzyme] = np.array(enzyme_concentration)
        # save enzyme.concentration
        enzyme.concentration = enzyme.quantity/total_occupied_volume
    return enzyme_conc

#%%
def michaelis_menten_term(k_cat, k_M, c_enzyme, c_substrate, hill=1):
    return k_cat * c_enzyme * c_substrate**hill / (k_M**hill + c_substrate**hill)

def spontaneous_reaction_term(k, c_reactant):
    return k*c_reactant

#%%
def reaction_diffusion_system(r_mesh, variable_values_2d_array):
    """ Returns right-hand side.
    Return array with shape (n,m), n = number of variables, m = number of nodes
    r_mesh has shape (m,)
    value has shape (n, m) # for each substance n at each node m a specific value (float)

    THE VALUES IN R_MESH CAN CHANGE, so the enzymes concentration per mesh
    point have to be reevaluated each time according to the r_mesh being used...
    """
    print(len(r_mesh))
    # Get the information from the 2d array for each species
    species_conc = { # concentration
        species.name : variable_values_2d_array[species_variable_indices[species]["prim"]]
        for species in reaction_network.species
    }
    species_conc_der = { # 1st derivative of concentration
        species.name : variable_values_2d_array[species_variable_indices[species]["1stDer"]]
        for species in reaction_network.species
    }

    # Get the concentration of enzymes at the r_mesh nodes
    enzyme_conc = enzyme_concentration(r_mesh)

    system_right_hand_side = { # includes no reaction terms here
        species: {"prim": species_conc_der[species.name],
                       "1stDer": -2/r_mesh * species_conc_der[species.name]}
        for species in reaction_network.species
    }

    # Calculate the reaction terms
    for species in reaction_network.species:
        reaction_term_sum = 0
        for reaction in species.as_reactant_in + species.as_product_in:
            if isinstance(reaction, SpontaneousReaction):
                term = reaction.k * species_conc[reaction.start_species.name] #################################
            elif isinstance(reaction, EnzymaticReaction):
                term = reaction.k_cat * enzyme_conc[reaction.enzyme] * species_conc[reaction.start_species.name] / (reaction.k_M + species_conc[reaction.start_species.name])
            if reaction in species.as_reactant_in: # if acts as reactant, diminishes
                term *= -1
            reaction_term_sum += term
        # Add reaction term to diffusion term
        system_right_hand_side[species]["1stDer"] += -1/species.diffusion_constant * reaction_term_sum
    
    # Create the right-hand side 2d array
    values_2d_array = np.zeros((2*len(reaction_network.species), len(r_mesh)))
    for species in reaction_network.species:
        for variable_type in species_variable_indices[species].keys():
            # writes the species_variable_indices[species][variable_type]th row of the 2d array
            values_2d_array[species_variable_indices[species][variable_type]] = system_right_hand_side[species][variable_type]

    return values_2d_array

def boundary_conditions(y_a, y_b):
    """ At origin a we have reflexion, so the derivative variable is 0.
    At r = R a flux is given.
    y_a and y_b are the vectors with shape (n, ) (n=2*num_species)
    """
    # origin
    origin_conditions = [
        y_a[species_variable_indices[species]["1stDer"]] - 0
        for species in reaction_network.species
    ]
    # edge
    border_conditions = [(
        y_b[species_variable_indices[species]["1stDer"]]
        - species.permeability_constant / species.diffusion_constant
        * (species.external_concentration - y_b[species_variable_indices[species]["prim"]]))
        for species in reaction_network.species
    ]
    conditions = np.array(origin_conditions + border_conditions)
    return conditions
#%%
#def solve_boundary_value_problem()

res_a = solve_bvp(reaction_diffusion_system, boundary_conditions,
                  r_mesh, initial_values_2d_array)
#%%
for enzyme_idx, enzyme in enumerate(reaction_network.enzymes):
    print(enzyme.name)
    print(enzyme.localization)
    print(enzyme.concentration)
#%%
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from scipy import integrate

final_mesh = res_a.x
solution_values_at_mesh = res_a.y
solution_derivatives_at_mesh = res_a.yp

success = res_a.success

integrals = {species.name: 0 for species in reaction_network.species}
x_values = np.linspace(0, radius, num = 100)
num_enzymes = len(reaction_network.enzymes)

enzyme_max_concentration = max([enzyme.concentration for enzyme in reaction_network.enzymes])

norm = mcolors.Normalize(vmin=0, vmax=enzyme_max_concentration)
cmap = cm.get_cmap('Oranges')
scalar_map = cm.ScalarMappable(norm=norm, cmap=cmap)

fig, ax = plt.subplots(num_enzymes + 1,1, figsize = (5,3), sharex = True, constrained_layout=True,
                       gridspec_kw={'height_ratios': [1]*num_enzymes + [num_enzymes]})
ax[0].set_title(f"steady state distribution, computed with {len(res_a.x)} nodes")

for enzyme_idx, enzyme in enumerate(reaction_network.enzymes):
    ax[enzyme_idx].set_ylabel(enzyme.name, rotation=0, labelpad=20)
    color = scalar_map.to_rgba(enzyme.concentration)
    for localizationTuple in enzyme.localization:
        min_range = localizationTuple.minMaxLoc[0]
        max_range = localizationTuple.minMaxLoc[1]
        ax[enzyme_idx].fill_between(
            x_values, 0, 1, where=(x_values >= min_range) & (x_values <= max_range),
            color=color)
    ax[enzyme_idx].set_yticks([])
for species in reaction_network.species:
    sol = res_a.sol(x_values)[species_variable_indices[species]["prim"]] 
        # Found solution for y as scipy.interpolate.PPoly instance, a C1 continuous
        # cubic spline.
    integral = integrate.quad(lambda r: res_a.sol(r)[species_variable_indices[species]["prim"]] * 4 * np.pi * r**2, 0, radius)
    integrals[species.name] = integral[0] # [0] is value [1] is error
    ax[-1].plot(x_values, sol, label=f"{species.name}: {int(np.round(integral[0]))} nM")
for node_x in res_a.x:
    ax[-1].axvline(node_x, ymin = 0.95, ymax = 1, c = "k", linewidth = 1)

ax[-1].set_xlabel("radius")
ax[-1].set_xlim([0, radius])
cbar = fig.colorbar(scalar_map, ax=ax[:num_enzymes], orientation='vertical', fraction=0.02, pad=0.04)
cbar.set_label("concentration")

species_integral_sum = np.sum([
    integrals[species.name]
    for species in reaction_network.species])

#annotation_list = [
#    f"{species.name}: {int(np.round(integrals[species.name]/species_integral_sum * 100))}%"
#    for species in reaction_network.species]
#annotation = ", ".join(annotation_list)

#ax[-1].annotate(annotation, (0.05, 0.9), xycoords = "axes fraction", va = "top", ha = "left", fontsize = 9)
ax[-1].set_xlabel("r / nm")
ax[-1].set_ylabel("nM / nm")

fig.savefig(os.path.join(violacein_folder, "test.png"), dpi = 300,
    bbox_inches = "tight")

# Assume scalar_map is your ScalarMappable (from cmap and norm)
fig_colorbar = plt.figure(figsize=(2, 4))  # tall figure for vertical bar

# Add axes for colorbar (full figure area)
cbar_ax = fig_colorbar.add_axes([0.2, 0.05, 0.3, 0.9])  # [left, bottom, width, height]

# Create the colorbar
cbar = fig_colorbar.colorbar(scalar_map, cax=cbar_ax, orientation='vertical')
cbar.set_label("Enzyme concentration")

# Save the colorbar figure
fig_colorbar.savefig("colorbar_only.png", bbox_inches='tight')
plt.close(fig_colorbar)

fig_legend = plt.figure(figsize=(3, 2))
ax_legend = fig_legend.add_subplot(111)
ax_legend.axis('off')

# Example dummy patches
red_patch = mpatches.Patch(color='red', label='Class A')
blue_patch = mpatches.Patch(color='blue', label='Class B')

legend = ax_legend.legend(handles=[red_patch, blue_patch], loc='center')
fig_legend.savefig('legend_only.png', bbox_inches='tight')
plt.close(fig_legend)



#ax[-1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.35), ncol = 3)

# Create legend fig separate
handles, labels = ax[-1].get_legend_handles_labels()
# Create a new figure just for the legend
fig_legend, ax_legend = plt.subplots(1,1, figsize = (3,2))

# Hide the axes frame and ticks
ax_legend.axis('off')
# Create the legend on this new axes
legend = ax_legend.legend(handles, labels, loc='center')

# Save the legend figure
fig_legend.savefig(os.path.join(violacein_folder, "test_legend.png"), dpi = 300,
    bbox_inches = "tight")

#%%
if __name__ == "__main__":
    folder_to_check_validity = sys.argv[1]
    reaction_network_info_file_names = load_json("src/reaction_network_info.json").keys()
    create_reaction_network(folder_to_check_validity, reaction_network_info_file_names)


