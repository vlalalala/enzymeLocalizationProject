import argparse
from auxiliary_functions_framework_organization import (
    get_dict_with_correct_key_types_from_json_file
)
import os
from auxiliary_functions import dump_json
from auxiliary_functions_using_standard_library import pickle_load_binary, load_json
import sys

def get_outward_fluxes(species_concentrations_dict, reaction_network, num_regions, num_mesh_points_in_regions):
    """
    """
    fluxes = {}
    for species in reaction_network.species:
        try:
            concentration_at_R = species_concentrations_dict[num_regions-1][num_mesh_points_in_regions[num_regions-1]-1][species]
        except:
            concentration_at_R = species_concentrations_dict[num_regions-1][num_mesh_points_in_regions[num_regions-1]-1][species.name]
        fluxes.update({species.name : species.permeability_constant * (concentration_at_R - species.external_concentration)})
    return fluxes

if __name__ == "__main__":
    # Parse arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("folder_to_solve", type=str, help="Path to folder with system info")
    args = parser.parse_args()
    FOLDER_TO_SOLVE = args.folder_to_solve

    # Do not run if the simulation has been pruned. Directly create dummy output file and exit.
    if os.path.isfile(os.path.join(FOLDER_TO_SOLVE, "pruned.json")):
        dump_json(FOLDER_TO_SOLVE, "fluxes", {"pruned": True})
        sys.exit(0)

    # Load inputs and define global parameters
    REACTION_NETWORK = pickle_load_binary(os.path.join(FOLDER_TO_SOLVE, ".pickled_reaction_network"))
    SYSTEM_GEOMETRY_DICT = load_json(os.path.join(FOLDER_TO_SOLVE, "system_geometry_for_convergence.json"))
    SYSTEM_MESH_DICT= load_json(os.path.join(FOLDER_TO_SOLVE, ".expanded_system_mesh_for_convergence.json"))

    SPECIES_LOOKUP = {sp.name: sp for sp in REACTION_NETWORK.species}

    R = SYSTEM_GEOMETRY_DICT["geometry_config"]["outer_membrane_radius"]
    NUM_MESH_POINTS_IN_REGIONS = SYSTEM_GEOMETRY_DICT["geometry_config"]["num_mesh_points_in_regions"]
    NUM_REGIONS = SYSTEM_GEOMETRY_DICT["geometry_config"]["num_regions"]
    MEMBRANE_RADII = SYSTEM_GEOMETRY_DICT["geometry_config"]["membrane_radii"]
    
    
    species_concentrations_with_strings = load_json(
        os.path.join(FOLDER_TO_SOLVE, ".species_steady_state_concentrations.json")
    )
    species_concentrations_dict = get_dict_with_correct_key_types_from_json_file(
        species_concentrations_with_strings, SPECIES_LOOKUP)

    fluxes = get_outward_fluxes(species_concentrations_dict, REACTION_NETWORK, NUM_REGIONS, NUM_MESH_POINTS_IN_REGIONS)
    dump_json(FOLDER_TO_SOLVE, "fluxes", fluxes)