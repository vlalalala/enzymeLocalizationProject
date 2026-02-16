import numpy as np
from auxiliary_functions import read_yaml_file

"""
incompatible_enzyme_groups_list: [[A,B], [A,B,C]] # for example
"""

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
    return reaction_network

def assert_no_conflicts_in_enzyme_positioning(reaction_network, num_regions, enzyme_impossible_combinations):
    enzymes_in_regions = {region_idx: [] for region_idx in range(num_regions)}
    for enzyme in reaction_network.enzymes:
        for region in enzyme.regions:
            enzymes_in_regions[region].append(enzyme.name)
    sorted_enzyme_impossible_combinations = [
        sorted(combination) for combination in enzyme_impossible_combinations
    ]
    for region_idx, enzymes_in_region in enzymes_in_regions.items():
        if sorted(enzymes_in_region) in sorted_enzyme_impossible_combinations:
            raise ValueError(f"Combination of enzymes {enzymes_in_region} within the region {region_idx} incompatible.")
    return True

def return_enzyme_concentrations(reaction_network, membrane_radii, maximum_enzyme_concentration):
    """Returns a dictionary with keys [region_idx][species]
    and the modified reaction_network with the concentrations.
    """
    enzyme_concentrations = {
        region_idx: {
            enzyme: None
            for enzyme in reaction_network.enzymes
        }
        for region_idx in range(len(membrane_radii))
    }
    for enzyme in reaction_network.enzymes:
        regions_volume = 0
        for region in enzyme.regions:
            if region == 0:
                regions_volume += 4/3 * np.pi * membrane_radii[0]**3
            else:
                max_radius = membrane_radii[region]
                min_radius = membrane_radii[region-1]
                regions_volume += 4/3 * np.pi * (max_radius**3 - min_radius**3)
        concentration = enzyme.quantity / regions_volume
        reaction_network.enzymes[enzyme.name].concentration = concentration
        if maximum_enzyme_concentration != None and concentration > maximum_enzyme_concentration:
            raise ValueError(f"Enzyme {enzyme.name} has a larger concentration than allowed.")
        for region_idx in enzyme_concentrations.keys():
            if region_idx in enzyme.regions:
                enzyme_concentrations[region_idx][enzyme] = concentration
            else:
                enzyme_concentrations[region_idx][enzyme] = 0
    
    return enzyme_concentrations, reaction_network
