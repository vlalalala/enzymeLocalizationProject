import argparse
from auxiliary_functions_framework_organization import (
    get_species_concentrations_from_json_file
)
from auxiliary_functions import dump_json


def get_outward_fluxes(species_concentrations_dict, reaction_network, num_regions, num_mesh_points_in_regions):
    """
    """
    fluxes = {}
    for species in reaction_network.species:
        concentration_at_R = species_concentrations_dict[num_regions-1][num_mesh_points_in_regions[num_regions-1]-1][species]
        fluxes.update({species.name : species.permeability_constant * (concentration_at_R - species.external_concentration)})
    return fluxes

if __name__ == "__main__":
    # Parse arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("folder_to_solve", type=str, help="Path to folder with system info")
    args = parser.parse_args()

    # Load all the passed information
    FOLDER_TO_SOLVE = args.folder_to_solve

    # Load inputs and define global parameters
    REACTION_NETWORK = pickle_load_binary(os.path.join(FOLDER_TO_SOLVE, ".pickled_reaction_network_final"))
    SYSTEM_GEOMETRY_DICT = load_json(os.path.join(FOLDER_TO_SOLVE, ".expanded_system_geometry.json"))
    SYSTEM_MESH_DICT= load_json(os.path.join(FOLDER_TO_SOLVE, ".expanded_system_mesh.json"))

    SPECIES_LOOKUP = {sp.name: sp for sp in REACTION_NETWORK.species}

    R = SYSTEM_GEOMETRY_DICT["geometry_config"]["outer_membrane_radius"]
    NUM_MESH_POINTS_IN_REGIONS = SYSTEM_GEOMETRY_DICT["geometry_config"]["num_mesh_points_in_regions"]
    NUM_REGIONS = SYSTEM_GEOMETRY_DICT["geometry_config"]["num_regions"]
    MEMBRANE_RADII = SYSTEM_GEOMETRY_DICT["geometry_config"]["membrane_radii"]
    
    species_concentrations_with_strings = load_json(
        os.path.join(FOLDER_TO_SOLVE, ".species_steady_state_concentrations.json")
    )
    species_concentrations_dict = get_species_concentrations_from_json_file(
        species_concentrations_with_strings, SPECIES_LOOKUP)

    fluxes = get_outward_fluxes(species_concentrations_dict, REACTION_NETWORK, NUM_REGIONS, NUM_MESH_POINTS_IN_REGIONS)
    dump_json(FOLDER_TO_SOLVE, "fluxes", fluxes)