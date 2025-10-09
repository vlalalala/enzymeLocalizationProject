#%%
import os
import itertools
import numpy as np
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
a = [1,2,3]
b = [4,5,6]
c = [7,8,9]
np.vstack((a,b,c))
#%%







#%%
for enzyme in reaction_network.enzymes:
    enzyme_mesh_occupied_bool_dict = {r: False for r in r_mesh} # Initialize whether mesh points have enzyme
    total_occupied_volume = 0
    for locTuple in enzyme.localization:
        volume_part = 4/3 * np.pi * ((locTuple.minMaxLoc[1] * radius)**3 - (locTuple.minMaxLoc[0] * radius)**3)
        total_occupied_volume += volume_part
        for r_key in enzyme_mesh_occupied_bool_dict.keys():
            if locTuple.return_within_tuple(r_key):
                enzyme_mesh_occupied_bool_dict[r_key] = True
    enzyme.total_volume = total_occupied_volume    
    enzyme.concentration = [
        enzyme.quantity/enzyme.total_volume for r in r_mesh
        if enzyme_mesh_occupied_bool_dict[r] == True else 0
    ]################################ not really correct... but will work for now

#%% SPECIES VARIABLES
reaction_network.species
#%%
species_variables = [[f"{s}_prim", f"{s}_1stDer"] for s in reaction_network.species]
species_variables = list(itertools.chain.from_iterable(species_variables))
species_variables
#%%
reaction_network.species[0].diffusion_constant
#%%
type(species_variables[0])
#%%
initial_values_species_guess = np.zeros((len(species_variables), r_mesh.size))
print(initial_values_species_guess.shape, "is (n,m); see documentation of solve_bvp" )
#%%
#%%
def reaction_diffusion_system_1(r_mesh, values):
    """ Right hand side.
    Return array with shape (n,m), n = number of variables, m = number of nodes
    r_mesh has shape (m,)
    value has shape (n, m) # for each substance n at each node m a specific value (float)
    """
    c = {substance: values[sub_idx[substance]] for substance in substances}
    e = {key: None for key in enzyme_total_concentrations.keys()}
    for key, total_concentration in enzyme_total_concentrations.items():
        nonnormalized_values = np.array([enzyme_localized_concentrations[key](r)
                        for r in r_mesh])
        e[key] = nonnormalized_values * total_concentration / np.sum(nonnormalized_values)
    #for key in e.keys():
    #    print(key, np.sum(e[key]), e[key])
    #print(e["VioA"])
    # TRP
    der_Trp = values[sub_idx["l-Trp_1"]]
    #print(len(e["VioA"]), len(c["l-Trp"]))
    F_Trp = -k_cat["VioA"] * e["VioA"] * c["l-Trp"] / (K_m["VioA"] + c["l-Trp"] )
    der_Trp_1 = (-F_Trp) / D - 2 * values[sub_idx["l-Trp_1"]] / r_mesh
    # Imine
    der_imine = values[sub_idx["IPA_imine_1"]]
    F_imine = (k_cat["VioA"] * e["VioA"] * c["l-Trp"] / (K_m["VioA"] + c["l-Trp"]) # process from reaction by VioA 
            - 2 * k_cat["VioB"] *  e["VioB"] * c["IPA_imine"]**hill["VioB"] / (K_m["VioB"]**hill["VioB"] + c["IPA_imine"]**hill["VioB"]))
    der_imine_1 = (-F_imine) / D - 2 * values[sub_idx["IPA_imine_1"]] / r_mesh
    # Dimer
    der_dimer = values[sub_idx["IPA_imine_dimer_1"]]
    F_dimer = (k_cat["VioB"] *  e["VioB"] * c["IPA_imine"]**hill["VioB"] / (K_m["VioB"]**hill["VioB"] + c["IPA_imine"]**hill["VioB"]) # process from reaction by VioB
        - k_n["IPA_imine_dimer->CPA"] * c["IPA_imine_dimer"] # spontaneous reaction to CPA
        - k_cat["VioE"] * e["VioE"] * c["IPA_imine_dimer"] / (K_m["VioE"] + c["IPA_imine_dimer"])) # from reaction from VioE
    der_dimer_1 = (-F_dimer) / D - 2 * values[sub_idx["IPA_imine_dimer_1"]] / r_mesh

    return np.vstack(
        (der_Trp, der_Trp_1,
         der_imine, der_imine_1,
         der_dimer, der_dimer_1,
         ))





#%%
substances = ["A", "B", "C"]

substances_variables = [[substance, substance + "_1"] for substance in substances]
substances_variables = list(itertools.chain.from_iterable(substances_variables))
sub_idx = {variable: index for index, variable in enumerate(substances_variables)}
derivative_variables = [key for key in sub_idx.keys() if key not in substances]

print("sub_idx", sub_idx)
print("derivative_variables", derivative_variables)

#%%
def reaction_diffusion_system_1(r_mesh, values):
    """ Right hand side.
    Return array with shape (n,m), n = number of variables, m = number of nodes
    r_mesh has shape (m,)
    value has shape (n, m) # for each substance n at each node m a specific value (float)
    """
    concentration = {substance: values[sub_idx[substance]] for substance in substances}
    e = {key: None for key in enzyme_total_concentrations.keys()}
    for key, total_concentration in enzyme_total_concentrations.items():
        nonnormalized_values = np.array([enzyme_localized_concentrations[key](r)
                        for r in r_mesh])
        e[key] = nonnormalized_values * total_concentration / np.sum(nonnormalized_values)
    #for key in e.keys():
    #    print(key, np.sum(e[key]), e[key])
    #print(e["VioA"])
    # TRP
    der_Trp = values[sub_idx["l-Trp_1"]]
    #print(len(e["VioA"]), len(c["l-Trp"]))
    F_Trp = -k_cat["VioA"] * e["VioA"] * c["l-Trp"] / (K_m["VioA"] + c["l-Trp"] )
    der_Trp_1 = (-F_Trp) / D - 2 * values[sub_idx["l-Trp_1"]] / r_mesh
    # Imine
    der_imine = values[sub_idx["IPA_imine_1"]]
    F_imine = (k_cat["VioA"] * e["VioA"] * c["l-Trp"] / (K_m["VioA"] + c["l-Trp"]) # process from reaction by VioA 
            - 2 * k_cat["VioB"] *  e["VioB"] * c["IPA_imine"]**hill["VioB"] / (K_m["VioB"]**hill["VioB"] + c["IPA_imine"]**hill["VioB"]))
    der_imine_1 = (-F_imine) / D - 2 * values[sub_idx["IPA_imine_1"]] / r_mesh
    # Dimer
    der_dimer = values[sub_idx["IPA_imine_dimer_1"]]
    F_dimer = (k_cat["VioB"] *  e["VioB"] * c["IPA_imine"]**hill["VioB"] / (K_m["VioB"]**hill["VioB"] + c["IPA_imine"]**hill["VioB"]) # process from reaction by VioB
        - k_n["IPA_imine_dimer->CPA"] * c["IPA_imine_dimer"] # spontaneous reaction to CPA
        - k_cat["VioE"] * e["VioE"] * c["IPA_imine_dimer"] / (K_m["VioE"] + c["IPA_imine_dimer"])) # from reaction from VioE
    der_dimer_1 = (-F_dimer) / D - 2 * values[sub_idx["IPA_imine_dimer_1"]] / r_mesh
    # CPA
    der_CPA = values[sub_idx["CPA_1"]]
    F_CPA = k_n["IPA_imine_dimer->CPA"] * c["IPA_imine_dimer"]
    der_CPA_1 = (-F_CPA) / D - 2 * values[sub_idx["CPA_1"]] / r_mesh
    # PDVA
    der_PDVA = values[sub_idx["PDVA_1"]]
    F_PDVA = (k_cat["VioE"] * e["VioE"] * c["IPA_imine_dimer"] / (K_m["VioE"] + c["IPA_imine_dimer"]) # process from reaction by VioE
        - k_cat["VioC_PDVA"] * e["VioC"] * c["PDVA"] * K_m["VioC_PVA"]/ (K_m["VioC_PDVA"] * (c["PVA"] + K_m["VioC_PVA"]) + c["PDVA"] * K_m["VioC_PVA"])
        - k_cat["VioD"] * e["VioD"] * c["PDVA"] / (K_m["VioD"] + c["PDVA"])
        - k_n["PDVA->PDV"] * c["PDVA"])
    der_PDVA_1 = (-F_PDVA) / D - 2 * values[sub_idx["PDVA_1"]] / r_mesh
    # PDV
    der_PDV = values[sub_idx["PDV_1"]]
    F_PDV = k_n["PDVA->PDV"] * c["PDVA"]
    der_PDV_1 = (-F_PDV) / D - 2 * values[sub_idx["PDV_1"]] / r_mesh
    # DVA
    der_DVA = values[sub_idx["DVA_1"]]
    F_DVA = (k_cat["VioC_PDVA"] * e["VioC"] * c["PDVA"] * K_m["VioC_PVA"]/ (K_m["VioC_PDVA"] * (c["PVA"] + K_m["VioC_PVA"]) + c["PDVA"] * K_m["VioC_PVA"])
        - k_n["DVA->DV"] * c["DVA"])
    der_DVA_1 = (-F_DVA) / D - 2 * values[sub_idx["DVA_1"]] / r_mesh
    # DV
    der_DV = values[sub_idx["DV_1"]]
    F_DV = k_n["DVA->DV"] * c["DVA"]
    der_DV_1 = (-F_DV) / D - 2 * values[sub_idx["DV_1"]] / r_mesh
    # PVA
    der_PVA = values[sub_idx["PVA_1"]]
    F_PVA = (k_cat["VioD"] * e["VioD"] * c["PDVA"] / (K_m["VioD"] + c["PDVA"])
        - k_cat["VioC_PVA"] * e["VioC"] * c["PVA"] * K_m["VioC_PDVA"]/ (K_m["VioC_PVA"] * (c["PDVA"] + K_m["VioC_PDVA"]) + c["PVA"] * K_m["VioC_PDVA"])
        - k_n["PVA->PV"] * c["PVA"])
    der_PVA_1 = (-F_PVA) / D - 2 * values[sub_idx["PVA_1"]] / r_mesh
    # PV
    der_PV = values[sub_idx["PV_1"]]
    F_PV = k_n["PVA->PV"] * c["PVA"]
    der_PV_1 = (-F_PV) / D - 2 * values[sub_idx["PV_1"]] / r_mesh
    # VA
    der_VA = values[sub_idx["VA_1"]]
    F_VA = (k_cat["VioC_PVA"] * e["VioC"] * c["PVA"] * K_m["VioC_PDVA"]/ (K_m["VioC_PVA"] * (c["PDVA"] + K_m["VioC_PDVA"]) + c["PVA"] * K_m["VioC_PDVA"])
        - k_n["VA->violacein"] * c["VA"])
    der_VA_1 = (-F_VA) / D - 2 * values[sub_idx["VA_1"]] / r_mesh
    # violacein
    der_violacein = values[sub_idx["violacein_1"]]
    F_violacein = k_n["VA->violacein"] * c["VA"]
    der_violacein_1 = (-F_violacein) / D - 2 * values[sub_idx["violacein_1"]] / r_mesh

    return np.vstack(
        (der_Trp, der_Trp_1,
         der_imine, der_imine_1,
         der_dimer, der_dimer_1,
         der_CPA, der_CPA_1,
         der_PDVA, der_PDVA_1,
         der_PDV, der_PDV_1,
         der_DVA, der_DVA_1,
         der_DV, der_DV_1,
         der_PVA, der_PVA_1,
         der_PV, der_PV_1,
         der_VA, der_VA_1,
         der_violacein, der_violacein_1
         ))






#%%
y_a = ["A_y", "A_z", "B_y", "B_z", "C_y", "C_z"]
print([y_a[sub_idx[derivative_variable]]
                         for derivative_variable in derivative_variables])

#%%
for species in reaction_network:
    species.concentration = np.zeros(system_geometry_dict["GEOMETRY_CONFIG"]["num_mesh_points"])

def boundary_conditions_1(y_a, y_b):
    """ At origin a we have reflexion, so the derivative variable is 0.
    At r = R a flux is given.
    """
    # origin
    origin_conditions = [y_a[sub_idx[derivative_variable]]
                         for derivative_variable in derivative_variables]
    # edge
    border_conditions = [(
        y_b[sub_idx[substance + "_1"]]
        - substance_permeability[substance] / D
        * (substance_concentration_in_medium[substance] - y_b[sub_idx[substance]]))
        for substance in substances]
    conditions = np.array(origin_conditions + border_conditions)
    return conditions


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
