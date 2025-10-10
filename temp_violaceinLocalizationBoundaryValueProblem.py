# pylint: disable = invalid-name
# pylint: disable = line-too-long
# pylint: disable = too-many-lines
# pylint: disable = too-many-arguments
# pylint: disable = too-many-instance-attributes
# pylint: disable = unused-argument
# pylint: disable = too-many-locals
# pylint: disable = too-many-statements

"""
This is the source code for the theory part on characterization of the pathway
that leads to violacein and other side products.
"""

# %% IMPORTS
import os
import copy
import dill as pickle
import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy import integrate
import sys
#from tqdm import tqdm
import itertools
from scipy.integrate import solve_bvp
import heapq

#%%####################################################
FOLDER_WITH_DATA = os.path.join(os.getcwd(), "data")
print(FOLDER_WITH_DATA)
#####################################################

# %% ENZYME CHARACTERISTICS
MICHAELIS_CONSTANT_VALUE_DICTIONARY = { # given in nM
    "VioA": 125e3, 
    "VioB": 150e3, 
    "VioC_PDVA": 100e3,
    "VioC_PVA": 100e3, 
    "VioD":  75e3, 
    "VioE":  175e3,
}
CATALYTIC_RATES_VALUE_DICTIONARY = {
    "VioA": 0.75 * 60, #1/min,
    "VioB": 0.5 * 60,
    "VioC_PDVA": 1.15 * 60,
    "VioC_PVA": 1.15 * 60,
    "VioD": 1.25 * 60,
    "VioE": 1 * 60,
}

ENZYMATIC_LINEAR_REACTION_RATES_VALUE_DICTIONARY = {
    key: CATALYTIC_RATES_VALUE_DICTIONARY[key]/MICHAELIS_CONSTANT_VALUE_DICTIONARY[key]
    for key in CATALYTIC_RATES_VALUE_DICTIONARY # pylint: disable = consider-using-dict-items
}

HILL_COEFFICIENTS_VALUE_DICTIONARY = {"VioB": 1.75}
ENZYMATIC_LINEAR_REACTION_RATES_VALUE_DICTIONARY["VioB"] /= MICHAELIS_CONSTANT_VALUE_DICTIONARY["VioB"]**(HILL_COEFFICIENTS_VALUE_DICTIONARY["VioB"]-1)

# %% NONENZYME CHARACTERISTICS
NONENZYMATIC_REACTION_RATES_VALUE_DICTIONARY = {
    "IPA_imine_dimer->CPA": 0.001,
    "PDVA->PDV": 0.001,
    "DVA->DV": 0.2 * 60,
    "PVA->PV": 0.0003 * 60,
    "VA->violacein": 0.001 * 60
}

#%% BOUNDARY VALUE PROBLEM. CONDITIONS.
D = 1000e-9 * 60 #nm**2 / min
p_s = 1e-9 * 60 #nm / min
p_i = 18e3 * 60 #nm / min
s_0 = 25e-3 #nM
R = 100 # nm
r_mesh = np.linspace(0.01, R, num = 100) #nm

K_m = MICHAELIS_CONSTANT_VALUE_DICTIONARY
k_cat = CATALYTIC_RATES_VALUE_DICTIONARY
k_n = NONENZYMATIC_REACTION_RATES_VALUE_DICTIONARY
hill = HILL_COEFFICIENTS_VALUE_DICTIONARY

# %% ####################################################################################################
# Make all enzymes pretty much uniform
for key in K_m.keys():
    K_m[key] = 100e3
    k_cat[key] = 10
hill["VioB"] = 1
for key in k_n.keys():
    k_n[key] = 0.0000001

VIOA_CONCENTRATION = 10#100 # nM
VIOB_CONCENTRATION = 100
VIOC_CONCENTRATION = 10
VIOD_CONCENTRATION = 20
VIOE_CONCENTRATION = 10

enzyme_total_concentrations = {
    "VioA": VIOA_CONCENTRATION,
    "VioB": VIOB_CONCENTRATION,
    "VioC": VIOC_CONCENTRATION,
    "VioD": VIOD_CONCENTRATION,
    "VioE": VIOE_CONCENTRATION,
}
#%% TESTING CONFIGURATIONS OF ENZYMES
INDEX = 0
enzyme_ranges = {
    "VioA": (R*0.9,R*1),
    "VioB": (R*0.7, R*0.8),
    "VioC": (R*0.0, R*0.1),
    "VioD": (R*0.0, R*0.1),
    "VioE": (R*0.4, R*0.5),
}

INDEX = 1
enzyme_ranges = {
    "VioA": (R*0.9,R*1),
    "VioB": (R*0.7, R*0.8),
    "VioC": (R*0.0, R*0.1),
    "VioD": (R*0.0, R*0.1),
    "VioE": (R*0.6, R*0.7),
}

INDEX = 2
enzyme_ranges = {
    "VioA": (R*0.9,R*1),
    "VioB": (R*0.7, R*0.8),
    "VioC": (R*0.0, R*0.1),
    "VioD": (R*0.5, R*0.6),
    "VioE": (R*0.6, R*0.7),
}

INDEX = 3
enzyme_ranges = {
    "VioA": (R*0.9,R*1),
    "VioB": (R*0.7, R*0.8),
    "VioC": (R*0.4, R*0.5),
    "VioD": (R*0.5, R*0.6),
    "VioE": (R*0.6, R*0.7),
}


# %% DEFINE ENZYME DENSITY
#enzyme_localized_concentrations = { # initialize with uniform enzyme concentration
#    key: (lambda r, value=value: value / (4*np.pi * r**2))
#    for key, value in enzyme_total_concentrations.items()
#}
# DENSITY ON 1D
enzyme_localized_concentrations = { # initialize with uniform enzyme concentration within range
    key: (lambda r, key = key, value=value: 3 * value / (enzyme_ranges[key][1]**3 - enzyme_ranges[key][0]**3) #* r**2 #########################################
    if (r>=enzyme_ranges[key][0] and r<=enzyme_ranges[key][1]) else 0)
    for key, value in enzyme_total_concentrations.items()
}

substances = ["l-Trp", "IPA_imine", "IPA_imine_dimer", "CPA", "PDVA", "PDV",
    "DVA", "DV", "PVA", "PV", "VA", "violacein"]

final_products = ["CPA", "PDV", "DV", "PV", "violacein"]

substance_concentration_in_medium = {key: s_0 if key == "l-Trp" else 0
    for key in substances}

substance_permeability = {key: p_s if key == "l-Trp" else p_i
    for key in substances}

initial_values_guess = np.zeros((len(substances)*2, r_mesh.size))
#initial_values_guess[sub_idx["l-Trp"]][-1] = 0

# %% SYSTEM WITH S1
substances_variables = [[substance, substance + "_1"] for substance in substances]
substances_variables = list(itertools.chain.from_iterable(substances_variables))
sub_idx = {variable: index for index, variable in enumerate(substances_variables)}
derivative_variables = [key for key in sub_idx.keys() if key not in substances]

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
    print("Trp", c["l-Trp"])
    print("VioA", e["VioA"])
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

res_a = solve_bvp(reaction_diffusion_system_1, boundary_conditions_1,
                  r_mesh, initial_values_guess)
final_mesh = res_a.x
solution_values_at_mesh = res_a.y
solution_derivatives_at_mesh = res_a.yp

success = res_a.success

#print("final mesh", len(final_mesh))
#print("solution values at mesh", solution_values_at_mesh[0][::10])
#print("solution_derivatives_at_mesh", solution_derivatives_at_mesh[0][::10])
#print("success", success)

integrals = {substance: 0 for substance in substances}
x_plot = np.linspace(0, R, num = 100)

fig, ax = plt.subplots(6,1, figsize = (5,3), sharex = True, gridspec_kw={'height_ratios': [1,1,1,1,1, 5]})
ax[0].set_title(f"steady state distribution, computed with {len(res_a.x)} nodes; (s, s1)")
for idx, (enzyme, enzyme_range) in enumerate(enzyme_ranges.items()):
    ax[idx].set_ylabel(enzyme, rotation=0, labelpad=20)
    min_range = enzyme_range[0]
    max_range = enzyme_range[1]
    ax[idx].fill_between(
        x_plot, 0, 1, where=(x_plot >= min_range) & (x_plot <= max_range),
        color='orange', alpha=0.5)
    ax[idx].set_yticks([])
    #ax[idx].annotate()
for substance in substances:
    idx = sub_idx[substance]
    sol = res_a.sol(x_plot)[idx] 
        # Found solution for y as scipy.interpolate.PPoly instance, a C1 continuous
        # cubic spline.
    integral = integrate.quad(lambda r: res_a.sol(r)[idx] * 4 * np.pi * r**2, 0, R)
    integrals[substance] = integral[0] # [0] is value [1] is error
    ax[-1].plot(x_plot, sol, label=f"{substance}: {int(np.round(integral[0]))} nM")
for node_x in res_a.x:
    ax[-1].axvline(node_x, ymin = 0.95, ymax = 1, c = "k", linewidth = 1)

final_products_integral = np.sum([integrals[product] for product in final_products])
CPA_percentage = int(np.round(integrals["CPA"]/final_products_integral * 100))
PDV_percentage = int(np.round(integrals["PDV"]/final_products_integral * 100))
DV_percentage = int(np.round(integrals["DV"]/final_products_integral * 100))
PV_percentage = int(np.round(integrals["PV"]/final_products_integral * 100))
violacein_percentage = int(np.round(integrals["violacein"]/final_products_integral * 100))
annotation = f"CPA: {CPA_percentage}%, PDV: {PDV_percentage}%, DV: {DV_percentage}%\nPV: {PV_percentage}%, violacein: {violacein_percentage}%"
ax[-1].annotate(annotation, (0.05, 0.9), xycoords = "axes fraction", va = "top", ha = "left", fontsize = 9)
ax[-1].set_xlabel("r / nm")
ax[-1].set_ylabel("concentration / nM ")
ax[-1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.35), ncol = 3)
plt.subplots_adjust(hspace=0)
fig.savefig(os.path.join(FOLDER_WITH_DATA,
    f"{INDEX}_steady_state_distribution_with_{len(res_a.x)}_nodes_s1.png"),
    dpi = 300,
    bbox_inches = "tight")

#integrals
# %% SYSTEM WITH S2
substances_variables = [[substance, substance + "_2"] for substance in substances]
substances_variables = list(itertools.chain.from_iterable(substances_variables))
sub_idx = {variable: index for index, variable in enumerate(substances_variables)}
derivative_variables = [key for key in sub_idx.keys() if key not in substances]

def boundary_conditions_2(y_a, y_b):
    """ At origin a we have reflexion, so the derivative variable is 0.
    At r = R a flux is given.
    """
    # origin
    origin_conditions = [y_a[sub_idx[derivative_variable]]
                         for derivative_variable in derivative_variables]
    # edge
    border_conditions = [(
        y_b[sub_idx[substance + "_2"]]
        - substance_permeability[substance]
        * (substance_concentration_in_medium[substance] - y_b[sub_idx[substance]])
        * R**2 / D)
        for substance in substances]
    conditions = np.array(origin_conditions + border_conditions)
    return conditions

def reaction_diffusion_system_2(r_mesh, values):
    """ Right hand side.
    Return array with shape (n,m), n = number of variables, m = number of nodes
    r_mesh has shape (m,)
    value has shape (n, m) # for each substance n at each node m a specific value (float)
    """
    #c = {substance: values[sub_idx[substance]] for substance in substances}
    #e = {key: np.array([value(r) for r in r_mesh])
    #     for key, value in enzyme_localized_concentrations.items()}
    c = {substance: values[sub_idx[substance]] for substance in substances}
    e = {key: None for key in enzyme_total_concentrations.keys()}
    for key, total_concentration in enzyme_total_concentrations.items():
        nonnormalized_values = np.array([enzyme_localized_concentrations[key](r)
                        for r in r_mesh])
        e[key] = nonnormalized_values * total_concentration / np.sum(nonnormalized_values)
    # TRP
    der_Trp = values[sub_idx["l-Trp_2"]] / r_mesh**2
    #print(len(e["VioA"]), len(c["l-Trp"]))
    F_Trp = -k_cat["VioA"] * e["VioA"] * c["l-Trp"] / (K_m["VioA"] + c["l-Trp"] )
    der_Trp_2 = r_mesh**2 / D * (-F_Trp)
    # Imine
    der_imine = values[sub_idx["IPA_imine_2"]] / r_mesh**2
    F_imine = (k_cat["VioA"] * e["VioA"] * c["l-Trp"] / (K_m["VioA"] + c["l-Trp"]) # process from reaction by VioA 
            - 2 * k_cat["VioB"] *  e["VioB"] * c["IPA_imine"]**hill["VioB"] / (K_m["VioB"]**hill["VioB"] + c["IPA_imine"]**hill["VioB"]))
    der_imine_2 = r_mesh**2 / D * (-F_imine)
    # Dimer
    der_dimer = values[sub_idx["IPA_imine_dimer_2"]] / r_mesh**2
    F_dimer = (k_cat["VioB"] *  e["VioB"] * c["IPA_imine"]**hill["VioB"] / (K_m["VioB"]**hill["VioB"] + c["IPA_imine"]**hill["VioB"]) # process from reaction by VioB
        - k_n["IPA_imine_dimer->CPA"] * c["IPA_imine_dimer"] # spontaneous reaction to CPA
        - k_cat["VioE"] * e["VioE"] * c["IPA_imine_dimer"] / (K_m["VioE"] + c["IPA_imine_dimer"])) # from reaction from VioE
    der_dimer_2 = r_mesh**2 / D * (-F_dimer)
    # CPA
    der_CPA = values[sub_idx["CPA_2"]] / r_mesh**2
    F_CPA = k_n["IPA_imine_dimer->CPA"] * c["IPA_imine_dimer"]
    der_CPA_2 = r_mesh**2 / D * (-F_CPA)
    # PDVA
    der_PDVA = values[sub_idx["PDVA_2"]] / r_mesh**2
    F_PDVA = (k_cat["VioE"] * e["VioE"] * c["IPA_imine_dimer"] / (K_m["VioE"] + c["IPA_imine_dimer"]) # process from reaction by VioE
        - k_cat["VioC_PDVA"] * e["VioC"] * c["PDVA"] * K_m["VioC_PVA"]/ (K_m["VioC_PDVA"] * (c["PVA"] + K_m["VioC_PVA"]) + c["PDVA"] * K_m["VioC_PVA"])
        - k_cat["VioD"] * e["VioD"] * c["PDVA"] / (K_m["VioD"] + c["PDVA"])
        - k_n["PDVA->PDV"] * c["PDVA"])
    der_PDVA_2 = r_mesh**2 / D * (-F_PDVA)
    # PDV
    der_PDV = values[sub_idx["PDV_2"]] / r_mesh**2
    F_PDV = k_n["PDVA->PDV"] * c["PDVA"]
    der_PDV_2 = r_mesh**2 / D * (-F_PDV)
    # DVA
    der_DVA = values[sub_idx["DVA_2"]] / r_mesh**2
    F_DVA = (k_cat["VioC_PDVA"] * e["VioC"] * c["PDVA"] * K_m["VioC_PVA"]/ (K_m["VioC_PDVA"] * (c["PVA"] + K_m["VioC_PVA"]) + c["PDVA"] * K_m["VioC_PVA"])
        - k_n["DVA->DV"] * c["DVA"])
    der_DVA_2 = r_mesh**2 / D * (-F_DVA)
    # DV
    der_DV = values[sub_idx["DV_2"]] / r_mesh**2
    F_DV = k_n["DVA->DV"] * c["DVA"]
    der_DV_2 = r_mesh**2 / D * (-F_DV)
    # PVA
    der_PVA = values[sub_idx["PVA_2"]] / r_mesh**2
    F_PVA = (k_cat["VioD"] * e["VioD"] * c["PDVA"] / (K_m["VioD"] + c["PDVA"])
        - k_cat["VioC_PVA"] * e["VioC"] * c["PVA"] * K_m["VioC_PDVA"]/ (K_m["VioC_PVA"] * (c["PDVA"] + K_m["VioC_PDVA"]) + c["PVA"] * K_m["VioC_PDVA"])
        - k_n["PVA->PV"] * c["PVA"])
    der_PVA_2 = r_mesh**2 / D * (-F_PVA)
    # PV
    der_PV = values[sub_idx["PV_2"]] / r_mesh**2
    F_PV = k_n["PVA->PV"] * c["PVA"]
    der_PV_2 = r_mesh**2 / D * (-F_PV)
    # VA
    der_VA = values[sub_idx["VA_2"]] / r_mesh**2
    F_VA = (k_cat["VioC_PVA"] * e["VioC"] * c["PVA"] * K_m["VioC_PDVA"]/ (K_m["VioC_PVA"] * (c["PDVA"] + K_m["VioC_PDVA"]) + c["PVA"] * K_m["VioC_PDVA"])
        - k_n["VA->violacein"] * c["VA"])
    der_VA_2 = r_mesh**2 / D * (-F_VA)
    # violacein
    der_violacein = values[sub_idx["violacein_2"]] / r_mesh**2
    F_violacein = k_n["VA->violacein"] * c["VA"]
    der_violacein_2 = r_mesh**2 / D * (-F_violacein)

    return np.vstack(
        (der_Trp, der_Trp_2,
         der_imine, der_imine_2,
         der_dimer, der_dimer_2,
         der_CPA, der_CPA_2,
         der_PDVA, der_PDVA_2,
         der_PDV, der_PDV_2,
         der_DVA, der_DVA_2,
         der_DV, der_DV_2,
         der_PVA, der_PVA_2,
         der_PV, der_PV_2,
         der_VA, der_VA_2,
         der_violacein, der_violacein_2
         ))

res_a = solve_bvp(reaction_diffusion_system_2, boundary_conditions_2,
                  r_mesh, initial_values_guess)
final_mesh = res_a.x
solution_values_at_mesh = res_a.y
solution_derivatives_at_mesh = res_a.yp

success = res_a.success

#print("final mesh", len(final_mesh))
#print("solution values at mesh", solution_values_at_mesh[0][::10])
#print("solution_derivatives_at_mesh", solution_derivatives_at_mesh[0][::10])
#print("success", success)

integrals = {substance: 0 for substance in substances}
x_plot = np.linspace(0, R, num = 100)

fig, ax = plt.subplots(6,1, figsize = (5,3), sharex = True, gridspec_kw={'height_ratios': [1,1,1,1,1, 5]})
ax[0].set_title(f"steady state distribution, computed with {len(res_a.x)} nodes; (s, s2)")
for idx, (enzyme, enzyme_range) in enumerate(enzyme_ranges.items()):
    ax[idx].set_ylabel(enzyme, rotation=0, labelpad=20)
    min_range = enzyme_range[0]
    max_range = enzyme_range[1]
    ax[idx].fill_between(
        x_plot, 0, 1, where=(x_plot >= min_range) & (x_plot <= max_range),
        color='orange', alpha=0.5)
    ax[idx].set_yticks([])
    #ax[idx].annotate()
for substance in substances:
    idx = sub_idx[substance]
    sol = res_a.sol(x_plot)[idx] 
        # Found solution for y as scipy.interpolate.PPoly instance, a C1 continuous
        # cubic spline.
    integral = integrate.quad(lambda r: res_a.sol(r)[idx] * 4 * np.pi * r**2, 0, R)
    integrals[substance] = integral[0] # [0] is value [1] is error
    ax[-1].plot(x_plot, sol, label=f"{substance}: {int(np.round(integral[0]))} nM")
for node_x in res_a.x:
    ax[-1].axvline(node_x, ymin = 0.95, ymax = 1, c = "k", linewidth = 1)
final_products_integral = np.sum([integrals[product] for product in final_products])
CPA_percentage = int(np.round(integrals["CPA"]/final_products_integral * 100))
PDV_percentage = int(np.round(integrals["PDV"]/final_products_integral * 100))
DV_percentage = int(np.round(integrals["DV"]/final_products_integral * 100))
PV_percentage = int(np.round(integrals["PV"]/final_products_integral * 100))
violacein_percentage = int(np.round(integrals["violacein"]/final_products_integral * 100))
annotation = f"CPA: {CPA_percentage}%, PDV: {PDV_percentage}%, DV: {DV_percentage}%\nPV: {PV_percentage}%, violacein: {violacein_percentage}%"
ax[-1].annotate(annotation, (0.05, 0.9), xycoords = "axes fraction", va = "top", ha = "left", fontsize = 9)

ax[-1].set_xlabel("r / nm")
ax[-1].set_ylabel("nM / nm")
ax[-1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.35), ncol = 3)
plt.subplots_adjust(hspace=0)
fig.savefig(os.path.join(FOLDER_WITH_DATA,
    f"{INDEX}_steady_state_distribution_with_{len(res_a.x)}_nodes_s2.png"),
    dpi = 300,
    bbox_inches = "tight")
#integrals


def closest_2_nodes(nodes, approx_min, approx_max):
    """
    find the element in the list with the smallest absolute difference with approx values
    """
    closest_nodes_to_min = heapq.nsmallest(2, nodes, key=lambda x: abs(x-approx_min))
    closest_nodes_to_max = heapq.nsmallest(2, nodes, key=lambda x: abs(x-approx_max))
    closest_node_to_min = closest_nodes_to_min[0]
    closest_node_to_max = closest_nodes_to_max[0]
    if not closest_node_to_min == closest_node_to_max:
        return closest_node_to_min, closest_node_to_max
    print("not yet implemented what happens if the node distance too small")
    raise NotImplementedError

def compute_enzyme_concentrations_at_nodes(
        nodes_r, total_concentration, r_min, r_max, change_limits = False):
    """ Distribute the enzyme uniformly across a ring on the sphere of
    radii between approx_r_min and approx_r_max. The r_min and r_max are made to 
    match the closest nodes.
    """
    if change_limits:
        r_min, r_max = closest_2_nodes(nodes_r, r_min, r_max)
    nodes_with_enzymes = [node for node in nodes_r
        if (node <= r_max and node >= r_min)]
    density_at_nodes = np.array(
        [total_concentration * 3 * r**2 / (r_max**3 - r_min**3) if r in nodes_with_enzymes else 0
         for r in nodes_r] 
    )
    density_at_nodes_normalized = density_at_nodes * total_concentration / np.sum(density_at_nodes)
    return density_at_nodes_normalized
# %%
