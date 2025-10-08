#%%
import os
from auxiliary_functions_using_standard_library import pickle_load_binary
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme

# Get the absolute path to this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Go up to the project root (one or more levels as needed)
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))


violacein_folder = os.path.join(
    PROJECT_ROOT, *["data", "violacein_0"]
)

#%%
reaction_network = pickle_load_binary(os.path.join(violacein_folder, ".REACTION_NETWORK_pickle"))
system_geometry_dict = pickle_load_binary(os.path.join(violacein_folder, ".SYSTEM_GEOMETRY_pickle"))
#%% Step 1: calculate total volume within which enzyme is found
radius = system_geometry_dict["GEOMETRY_CONFIG"]["radius"]
num_mesh_points = system_geometry_dict["GEOMETRY_CONFIG"]["num_mesh_points"]

#%%
reaction_network.enzymes["VioA"].localization
#%%
enzyme_localized_concentrations = { # initialize with uniform enzyme concentration within range
    key: (lambda r, key = key, value=value: 3 * value / (enzyme_ranges[key][1]**3 - enzyme_ranges[key][0]**3) * r**2
    if (r>=enzyme_ranges[key][0] and r<=enzyme_ranges[key][1]) else 0)
    for key, value in enzyme_total_concentrations.items()
}