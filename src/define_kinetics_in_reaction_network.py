
#%%
import pickle
import os
#%%
import json
def load_json(path):
    _, file_extension = os.path.splitext(os.path.basename(path))
    if not os.path.isfile(path) or file_extension != ".json":
        raise ValueError(f"The file {path} does not exist or is not a .json file.")
    with open(path, "r") as f:
        contents = json.load(f)
    return contents

def pickle_load_binary(path):
    with open(path, 'rb') as f:
        loaded_variable = pickle.load(f)
    return loaded_variable
#%%
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme
#%%
def pickle_load_binary(path):
    with open(path, 'rb') as f:
        loaded_variable = pickle.load(f)
    return loaded_variable
#%%
# Get the absolute path to this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Go up to the project root (one or more levels as needed)
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))


violacein_folder = os.path.join(
    PROJECT_ROOT, *["data", "violacein_0"]
)

#%% Open .Network system
reaction_network = pickle_load_binary(os.path.join(violacein_folder, *[".REACTION_NETWORK_pickle"]))
#%%
system_geometry_dict = pickle_load_binary(os.path.join(violacein_folder, ".SYSTEM_GEOMETRY_pickle"))
#%%
system_geometry_dict


#%%
for species in reaction_network:
    species.concentration = np.zeros(system_geometry_dict["GEOMETRY_CONFIG"]["num_mesh_points"])
#%%
reaction_network.species["Trp"].as_reactant_in

#%%
for species in reaction_network.species:
    for reaction in species.as_reactant_in + species.as_product_in:
        if isinstance(reaction, SpontaneousReaction):
            term = reaction.k * reaction.start_species #################################
        elif isinstance(reaction, EnzymaticReaction):
            term = reaction.k_cat * reaction.enzyme * reaction.start_species / (reaction.k_M + reaction.start_species)
        if reaction in species.as_reactant_in: # if acts as reactant, concentration diminishes
            term *= -1
        species.first_time_derivative_terms.append(term)



#%%
def define_differential_equations(case_folder):
    system = pickle_load_binary(os.path.join(case_folder, ".NETWORK_system_pickle"))

#%%
define_differential_equations(violacein_folder)

# %%
