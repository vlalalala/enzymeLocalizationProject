"""
This script is for estimating reaction and diffusion timescales.
"""
import os
import argparse
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme
from auxiliary_functions_using_standard_library import (
    pickle_load_binary, load_json, pickle_dump_binary
)
import shutil
import pandas as pd
from contextlib import redirect_stdout

def calculate_timescales(species, external_radius):
    """
    Returns a dictionary that gives the diffusion timescale and the reaction
    timescales for the different reactions (given as a dictionary, with the key being
    a Reaction instance)
    Does not account for sign in reactions (whether it is a reactant or product)
    species: instance of Species.
    external_radius: in m.
    """
    timescales = {"diffusion": None, "reaction": {}}
    timescales["diffusion"] = external_radius**2 / species.diffusion_constant
    for reaction in species.as_reactant_in + species.as_product_in:
        if isinstance(reaction, SpontaneousReaction):
            timescale = 1 / reaction.k # sign is placed afterwards
        else: #EnzymaticReaction
            max_enzyme_concentration = max(reaction.enzyme.regional_concentrations.values())
            timescale = reaction.k_M / (max_enzyme_concentration * reaction.k_cat)
        timescales["reaction"].update({reaction: timescale})
    print(f"species {species.name} timescales: \n", timescales, "\n",  flush = True)
    return timescales

def find_ratio_diffusion_to_fastest_reaction_timescales(timescales):
    """ Returns diffusion timescale divided by smallest reaction timescale (from fastest reaction)
    Basically, this function returns the maximum ratio.
    """
    smallest_reaction_timescale = min(timescales["reaction"].values())
    print(f"The smallest reaction timescale is {smallest_reaction_timescale}.")
    print("ratio diffusion to fastest reaction timescale: ", timescales["diffusion"] / smallest_reaction_timescale,
          flush=True)
    return timescales["diffusion"] / smallest_reaction_timescale

def find_factor_species_kinetics(
        wanted_maximum_initial_ratio_diffusion_to_fastest_reaction_timescales,
        maximum_timescale_ratio
    ):
    """ The initial simulation that should be close to the external concentrations
    should have a small ratio of timescales of diffusion to fastest reaction.
    If diffusion already dominates for that species (aka maximum_timescale_ratio is already small),
    nothing should be changed (aka maximum_timescale_ratio is already smaller than wanted_maximum_initial_ratio_diffusion_to_fastest_reaction_timescales)
    """
    ratio_timescale_to_goal_timescale = max(
        1,
        maximum_timescale_ratio/wanted_maximum_initial_ratio_diffusion_to_fastest_reaction_timescales
    )
    print("ratio timescale to goal timescale: ", ratio_timescale_to_goal_timescale, "\n", flush=True)
    return ratio_timescale_to_goal_timescale


def get_modified_reaction_network(
        reaction_network,
        external_radius,
        wanted_maximum_initial_ratio_diffusion_to_fastest_reaction_timescales = 0.1
    ):
    """    
    """
    # First figure out which factor each species requires at the least
    species_factor = {}
    for species in reaction_network.species:
        timescales = calculate_timescales(species, external_radius)
        maximum_timescale_ratio = find_ratio_diffusion_to_fastest_reaction_timescales(
            timescales)
        initial_factor = find_factor_species_kinetics(
            wanted_maximum_initial_ratio_diffusion_to_fastest_reaction_timescales,
            maximum_timescale_ratio
        )
        species_factor.update({species: initial_factor})
    
    print("species factors: \n ", species_factor, "\n", flush=True)
    max_factor = max(species_factor.values())
    print(f"The factor by which all of the reactions will be scaled will be 1 over {max_factor}.")
    if max_factor == 1:
        return reaction_network, False # modification is False
    # In case any modifications are happening
    print("\n enzymatic reactions: ", flush=True)
    for reaction in reaction_network.enzymatic_reactions:
        print("reaction: ", reaction.name, "previous k_cat: ", reaction.k_cat, flush=True)
        reaction.k_cat /= max_factor
        print("current k_cat: ", reaction.k_cat)
    print("\n spontaneous reactions: ", flush=True)
    for reaction in reaction_network.spontaneous_reactions:
        print("reaction: ", reaction.name, "previous k: ", reaction.k, flush=True)
        reaction.k /= max_factor   
        print("current k: ", reaction.k)
    return reaction_network, True
    # Then, find the largest of the factors. All of the kinetic rates will be changed
    # by that factor (so that the ratios in concentrations are kept relatively similar)
    # NOT THIS: THIS DOES NOT KEEP RATIOS Then, go through the different reactions and redefine the kinetic parameter
    # with the smallest factor of the species involved
    #print("\n enzymatic reactions: ", flush=True)
    #modification = False
    #for reaction in reaction_network.enzymatic_reactions:
    #    min_factor = min(species_factor[reaction.start_species], species_factor[reaction.end_species])
    #    print("reaction: ", reaction.name, "factor: ", min_factor, "previous k_cat: ", reaction.k_cat, flush=True)
    #    reaction.k_cat /= min_factor
    #    if min_factor != 1:
    #        modification = True
    #    print("current k_cat: ", reaction.k_cat)
    #print("\n spontaneous reactions: ", flush=True)
    #for reaction in reaction_network.spontaneous_reactions:
    #    min_factor = min(species_factor[reaction.start_species], species_factor[reaction.end_species])
    #    print("reaction: ", reaction.name, "factor: ",  min_factor, "previous k: ", reaction.k, flush=True)   
    #    reaction.k /= min_factor
    #    print("current k: ", reaction.k)
    #    if min_factor != 1:
    #        modification = True
    #return reaction_network, modification

def create_folder_with_scaled_data(folder_to_solve, new_reaction_network, path_to_new_folder):
    # Copy unmodified files
    # (species.csv, enzymes.csv, parameters_discretization.yaml, parameters_solver_input.yaml,
    # parameters_solver_output.yaml, parameters_value_conditions.yaml,
    # parameters_geometry.yaml)
    # onto new folder
    for file in [
        "species.csv", "enzymes.csv",
        "parameters_discretization.yaml", "parameters_solver_input.yaml",
        "parameters_solver_output.yaml", "parameters_value_conditions.yaml",
        "parameters_geometry.yaml"
    ]:  
        src = os.path.join(folder_to_solve, file)
        dst = os.path.join(path_to_new_folder, file)
        if not os.path.isfile(src):
            raise FileNotFoundError(f"Source file not found: {src}")
        shutil.copy(src, dst)
        if not os.path.isfile(dst): # In case shutil failed silently
            raise FileNotFoundError(f"Destination file not found: {dst}")
    
    # Change enzymatic reactions
    enzymatic_reactions_df = pd.read_csv(os.path.join(folder_to_solve, "enzymatic_reactions.csv"))
    for reaction in new_reaction_network.enzymatic_reactions:
        enzymatic_reactions_df.loc[
                (enzymatic_reactions_df['start_species'] == reaction.start_species.name)
                & (enzymatic_reactions_df['end_species'] == reaction.end_species.name)
                & (enzymatic_reactions_df['enzyme'] == reaction.enzyme),
            'k_cat'
        ] = reaction.k_cat # rewrite k_cat within the dataframe
    # Save the modified enzymatic_reactions dataframe
    enzymatic_reactions_df.to_csv(
        os.path.join(path_to_new_folder, "enzymatic_reactions.csv"),
        index=False
    )

    # Change spontaneous reactions
    spontaneous_reactions_df = pd.read_csv(os.path.join(folder_to_solve, "spontaneous_reactions.csv"))
    for reaction in new_reaction_network.spontaneous_reactions:
        spontaneous_reactions_df.loc[
                (spontaneous_reactions_df['start_species'] == reaction.start_species.name)
                & (spontaneous_reactions_df['end_species'] == reaction.end_species.name),
            'k'
        ] = reaction.k # rewrite k within the dataframe
    # Save the modified spontaneous_reactions dataframe
    spontaneous_reactions_df.to_csv(
        os.path.join(path_to_new_folder, "spontaneous_reactions.csv"),
        index=False
    )
    # Save the pickled_reaction_network
    pickle_dump_binary(os.path.join(path_to_new_folder, ".pickled_reaction_network"), new_reaction_network)

def create_creeping_reaction_folder(
        folder_to_solve,
        path_to_creeping_reaction_folder,
        wanted_maximum_ratio_diffusion_to_fastest_reaction_timescales
    ):
    os.makedirs(path_to_creeping_reaction_folder, exist_ok=True)
    # In order to print all of the information to a file
    with open(os.path.join(
            path_to_creeping_reaction_folder,
            f"initial_guess_calculator.log"), "a"
        ) as f, redirect_stdout(f):
        print(f"The wanted maximum ratio of diffusion to reaction timescales is {wanted_maximum_ratio_diffusion_to_fastest_reaction_timescales}. \n")
        # import the reaction network from the base parameters and the rest of the info
        reaction_network = pickle_load_binary(os.path.join(folder_to_solve, ".pickled_reaction_network"))
        system_geometry_dict = load_json(os.path.join(folder_to_solve, ".system_geometry.json"))
        external_radius = system_geometry_dict["geometry_config"]["outer_membrane_radius"]
        # Define the new reaction network and save it
        new_reaction_network, modification = get_modified_reaction_network(
            reaction_network, external_radius, wanted_maximum_ratio_diffusion_to_fastest_reaction_timescales
        )
        create_folder_with_scaled_data(
            folder_to_solve, new_reaction_network,
            path_to_creeping_reaction_folder
        )
        # Calculate the timescales again in order to check whether the current timescales adhere to what is needed
        print("\n After modifying the rates, the new timescales are: ")
        for species in new_reaction_network.species:
            new_timescales = calculate_timescales(species, external_radius)
            find_ratio_diffusion_to_fastest_reaction_timescales(new_timescales)
    return modification
        

#if __name__ == "__main__":
#    # Parse arguments from command line
#    parser = argparse.ArgumentParser()
#    parser.add_argument("--folder", type=str)
#    args = parser.parse_args()
#    FOLDER_TO_SOLVE = args.folder
#    reaction_network = pickle_load_binary(os.path.join(FOLDER_TO_SOLVE, ".pickled_reaction_network"))
#    SPECIES_LOOKUP = {sp.name: sp for sp in reaction_network.species}
#    SYSTEM_GEOMETRY_DICT = load_json(os.path.join(FOLDER_TO_SOLVE, ".system_geometry.json"))
#    external_radius = SYSTEM_GEOMETRY_DICT["geometry_config"]["outer_membrane_radius"]
#    new_reaction_network = get_modified_reaction_network(reaction_network, external_radius)
#    create_folder_with_scaled_data(FOLDER_TO_SOLVE, new_reaction_network)
