import sys
import matplotlib.pyplot as plt
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src'))
from auxiliary_functions_using_standard_library import load_json
import re
import pandas as pd
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
from auxiliary_functions import read_yaml_file
import ast

def calculate_average_concentration_of_X_in_inner_compartment_with_enzyme_A(combined_folder):
    system_geometry_file = os.path.join(combined_folder, "system_geometry_for_convergence.json")
    concentrations_file = os.path.join(combined_folder, ".species_steady_state_concentrations.json")
    if os.path.isfile(system_geometry_file) and os.path.isfile(concentrations_file):
        system_geometry = load_json(system_geometry_file)
        concentrations = load_json(concentrations_file)
    else:
        return None

    radii_in_region_1 = system_geometry["geometry_config"]["mesh_points_in_regions"][1]
    X_concentrations_in_region_1 = []
    for radius_idx in range(len(radii_in_region_1)):
        X_concentrations_in_region_1.append(concentrations[1][radius_idx]["X"])
    
    concentration = 0
    for i in range(len(radii_in_region_1) - 1):
        r0 = radii_in_region_1[i]
        r1 = radii_in_region_1[i+1]
        c0 = X_concentrations_in_region_1[i]
        c1 = X_concentrations_in_region_1[i+1]
        
        # Volume-weighted average concentration in shell (assuming linear c(r))
        c_avg = (c0 * (3*r0**2 + 2*r0*r1 + r1**2) + c1 * (r0**2 + 2*r0*r1 + 3*r1**2)) \
                / (4 * (r0**2 + r0*r1 + r1**2))
        
        shell_volume = r1**3 - r0**3  # 4π/3 cancels with denominator
        concentration += c_avg * shell_volume

    concentration /= (radii_in_region_1[-1]**3 - radii_in_region_1[0]**3)
    return concentration

def get_data(folder):
    data = {}
    for folder_name in os.listdir(folder):
        match = re.match(r'^combined_(\d{6})$', folder_name)
        combined_folder = os.path.join(folder, folder_name)
        if match:
            index = match.group(1)  # keeps it as a string, preserving leading zeros
            geometry = read_yaml_file(os.path.join(combined_folder, "parameters_geometry.yaml"))
            internal_radius = geometry["geometry_config"]["internal_membrane_relative_radii"][0]
            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            if os.path.isfile(fluxes_file):
                flux = load_json(fluxes_file)["Y"]
            else:
                flux = None
            enzymatic_reactions_df = pd.read_csv(os.path.join(combined_folder, "enzymatic_reactions.csv"))
            michaelis_menten_constant = enzymatic_reactions_df.loc[
                (enzymatic_reactions_df["enzyme"] == "A"),
                "k_M"].item()
            concentration = calculate_average_concentration_of_X_in_inner_compartment_with_enzyme_A(combined_folder)
            data[index] = (internal_radius, flux, michaelis_menten_constant, concentration)
    return data

def plot_data(folder):
    fig, ax = plt.subplots(2, 2, figsize = (7,7))
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=['internal_radius', 'flux', 'michaelis_menten_constant', 'concentration'])
    K_m_values = np.sort(df["michaelis_menten_constant"].unique())
    min_internal_radius = np.min(df['internal_radius'].unique())
    for _, Km in enumerate(K_m_values):
        df_current = df[df["michaelis_menten_constant"]==Km]
        ax[0][0].scatter(df_current["internal_radius"], df_current["flux"], label=Km)
        ax[1][0].set_xlabel("relative radius of most \n inner membrane r*/R")
        ax[0][0].set_ylabel("flux")
        ax[1][0].scatter(df_current["internal_radius"], df_current["concentration"])
        ax[1][0].set_ylabel("average X concentration \n in volume with enzyme")
        flux_when_internal_radius_is_smallest = df_current.loc[df_current['internal_radius'] == min_internal_radius, 'flux'].item()
        ax[0][1].scatter(df_current["internal_radius"], df_current["flux"]/flux_when_internal_radius_is_smallest)
        ax[0][1].set_xlabel("relative radius of most \n inner membrane r*/R")
        ax[0][1].set_ylabel("flux / flux found with smallest radius")
        ax[1][1].axis("off")

    ax[0][0].legend(title="Km")
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "result.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)
    # python data/02d_enzymaticXtoY_2InnerBoundaries_modifyingPositionOfEnzymeRegion_modifyingKm/analysis.py data/02d_enzymaticXtoY_2InnerBoundaries_modifyingPositionOfEnzymeRegion_modifyingKm