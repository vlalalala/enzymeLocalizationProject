import sys
import os
import numpy as np
from auxiliary_functions_using_standard_library import (
    pickle_load_binary, load_json, pickle_dump_binary)
from auxiliary_functions import read_yaml_file, dump_json
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme
import pandas as pd

def return_reaction_network_with_total_fixed_quantity_asserted(enzyme_total_fixed_quantity, reaction_network, enzyme_whose_quantity_to_modify):
    """ If enzyme_total_fixed_quantity is None, returns the reaction network inputted
    Else, enzyme_total_fixed_quantity is a float.
    """
    if enzyme_total_fixed_quantity is None:
        return reaction_network
    if enzyme_whose_quantity_to_modify not in [enzyme.name for enzyme in reaction_network.enzymes]:
        raise ValueError("could not find enzyme whose quantity to modify within the network")
    sum_of_enzyme_quantity_for_fixed_enzymes = sum([
        enzyme.quantity for enzyme in reaction_network.enzymes
        if enzyme.name != enzyme_whose_quantity_to_modify
    ])
    if sum_of_enzyme_quantity_for_fixed_enzymes > enzyme_total_fixed_quantity:
        raise ValueError("the sum of the quantity of fixed enzymes is larger than the specified total sum")
    for enzyme in reaction_network.enzymes:
        if enzyme.name == enzyme_whose_quantity_to_modify:
            enzyme.quantity = enzyme_total_fixed_quantity - sum_of_enzyme_quantity_for_fixed_enzymes
            enzyme.quantity_allocation = enzyme.allocate_quantity(enzyme.quantity, enzyme.allocation)
    return reaction_network

def define_regional_enzyme_concentrations(
        reaction_network, system_geometry_dict, maximum_concentration_allowed):

    # define the concentrations at each enzyme position
    for enzyme in reaction_network.enzymes:
        enzyme.regional_concentrations = {
            region: enzyme.quantity_allocation[region]/region_volume
            for region, region_volume in system_geometry_dict["geometry_config"]["volume_regions"].items()
        }
        if (maximum_concentration_allowed is not None
            and any(v > maximum_concentration_allowed for v in enzyme.regional_concentrations.values())):
            raise ValueError(f"A regional concentration of enzyme {enzyme.name}: {enzyme.regional_concentrations} is above {maximum_concentration_allowed}.")     

    return reaction_network


if __name__ == "__main__":
    folder = sys.argv[1]
    reaction_network = pickle_load_binary(
        os.path.join(folder, ".pickled_reaction_network_without_enzyme_concentration"))
    system_geometry_dict = load_json(os.path.join(folder, ".system_geometry.json"))
    parameter_value_conditions = read_yaml_file(os.path.join(folder, "parameters_value_conditions.yaml")) 

    if os.path.isfile(os.path.join(folder, "pruned.json")):
        # create dummy files for snakemake
        pickle_dump_binary(
            os.path.join(folder, ".pickled_reaction_network"),
            pd.DataFrame()
        )
        dump_json(folder, "enzyme_concentrations", {"pruned": "Enzyme allocation not valid. Pruned."})
        sys.exit()

    # Create final .species_steady_state_concentrations in case the parameter values are not viable
    # with the constraints given
    try:
        if "trial" not in folder:
            reaction_network = return_reaction_network_with_total_fixed_quantity_asserted(
                parameter_value_conditions["enzyme_total_fixed_quantity"],
                reaction_network,
                parameter_value_conditions["enzyme_whose_quantity_to_modify_when_total_fixed_quantity"]
            )
            print(f"modifying enzyme quantity of {parameter_value_conditions["enzyme_whose_quantity_to_modify_when_total_fixed_quantity"]}")
        else:
            print("defining enzyme concentrations without modifying quantities.")
        reaction_network = define_regional_enzyme_concentrations(
            reaction_network,
            system_geometry_dict,
            parameter_value_conditions["enzyme_maximum_concentration"]
        )
    except Exception as e:
        # EARLY TERMINATION
        dump_json(folder,
                ".species_steady_state_concentrations",
                {"error": f"{e}"})

    pickle_dump_binary(
            os.path.join(folder, ".pickled_reaction_network"),
            reaction_network
        )
    enzyme_concentrations = {
        enzyme: enzyme.regional_concentrations
        for enzyme in reaction_network.enzymes}
    
    # Save the enzyme concentrations in a human-readable format
    # (Important: save outside of enzymes.csv (yes, it would be tempting) but Snakemake will not like it)
    dump_json(folder, "enzyme_concentrations", enzyme_concentrations)