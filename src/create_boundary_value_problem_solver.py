#%%
import sys
import os
import numpy as np
from scipy import integrate
from auxiliary_functions_using_standard_library import pickle_load_binary, closest_value
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme
#%%
# Get the absolute path to this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Go up to the project root (one or more levels as needed)
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
# Access other folders relative to the root
DATA_PATH = os.path.join(PROJECT_ROOT, "data")

case_folder = os.path.join(
    DATA_PATH, "minimal_test_0"
)
#%%
folder_to_solve = case_folder
reaction_network = pickle_load_binary(os.path.join(folder_to_solve, ".REACTION_NETWORK_pickle"))
system_geometry_dict = pickle_load_binary(os.path.join(folder_to_solve, ".SYSTEM_GEOMETRY_pickle"))
#%%
# Step 1: Define all geometry variables
R = system_geometry_dict["GEOMETRY_CONFIG"]["outer_membrane_radius"]
MESH_POINTS_IN_REGIONS = system_geometry_dict["GEOMETRY_CONFIG"]["MESH_POINTS_IN_REGIONS"]
NUM_MESH_POINTS_IN_REGIONS = system_geometry_dict["GEOMETRY_CONFIG"]["NUM_MESH_POINTS_IN_REGIONS"]
NUM_REGIONS = system_geometry_dict["GEOMETRY_CONFIG"]["NUM_REGIONS"]

#%%
species_concentrations = {
    region_idx : {
        mesh_point_idx : {
            species : 0
            for species in reaction_network.species}
        for mesh_point_idx in range(NUM_MESH_POINTS_IN_REGIONS[region_idx])}
    for region_idx in range(NUM_REGIONS)
}
#%%
enzymes_concentrations = {
    region_idx : {
        enzyme : enzyme.concentration if region_idx in enzyme.regions else 0
        for enzyme in reaction_network.enzymes
    }
    for region_idx in range(NUM_REGIONS)
}
#%%
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
#%%
reverse_point_ids = {
    value: (region_idx, mesh_point_idx, species)
    for region_idx, mesh_points in point_ids.items()
    for mesh_point_idx, species_map in mesh_points.items()
    for species, value in species_map.items()
}
#%%
reverse_point_ids
#%%
point_info = {
    region_idx : {
        mesh_point_idx : "l" if mesh_point_idx==0 else ("r" if mesh_point_idx==NUM_MESH_POINTS_IN_REGIONS[region_idx]-1 else "i")
        for mesh_point_idx in range(NUM_MESH_POINTS_IN_REGIONS[region_idx])
    }
    for region_idx in range(NUM_REGIONS)
}
#%%
point_spatial_dependencies = 
#%%
point_info
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

    