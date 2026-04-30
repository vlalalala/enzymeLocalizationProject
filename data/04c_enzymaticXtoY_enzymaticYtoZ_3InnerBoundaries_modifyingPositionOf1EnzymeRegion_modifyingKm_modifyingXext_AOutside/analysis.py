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

def calculate_average_concentration_of_Y_in_inner_compartment_with_enzyme_B(combined_folder):
    system_geometry_file = os.path.join(combined_folder, "system_geometry_for_convergence.json")
    concentrations_file = os.path.join(combined_folder, ".species_steady_state_concentrations.json")
    if os.path.isfile(system_geometry_file) and os.path.isfile(concentrations_file):
        system_geometry = load_json(system_geometry_file)
        concentrations = load_json(concentrations_file)
    else:
        return None

    radii_in_region_2 = system_geometry["geometry_config"]["mesh_points_in_regions"][2]
    Y_concentrations_in_region_2 = []
    for radius_idx in range(len(radii_in_region_2)):
        Y_concentrations_in_region_2.append(concentrations[2][radius_idx]["Y"])
    
    concentration = 0
    for i in range(len(radii_in_region_2) - 1):
        r0 = radii_in_region_2[i]
        r1 = radii_in_region_2[i+1]
        c0 = Y_concentrations_in_region_2[i]
        c1 = Y_concentrations_in_region_2[i+1]
        
        # Volume-weighted average concentration in shell (assuming linear c(r))
        c_avg = (c0 * (3*r0**2 + 2*r0*r1 + r1**2) + c1 * (r0**2 + 2*r0*r1 + 3*r1**2)) \
                / (4 * (r0**2 + r0*r1 + r1**2))
        
        shell_volume = r1**3 - r0**3  # 4π/3 cancels with denominator
        concentration += c_avg * shell_volume

    concentration /= (radii_in_region_2[-1]**3 - radii_in_region_2[0]**3)
    return concentration

def get_data(folder):
    data = {}
    for folder_name in os.listdir(folder):
        match = re.match(r'^combined_(\d{6})$', folder_name)
        combined_folder = os.path.join(folder, folder_name)
        if match:
            index = match.group(1)  # keeps it as a string, preserving leading zeros
            geometry = read_yaml_file(os.path.join(combined_folder, "parameters_geometry.yaml"))
            pos_membrane = geometry["geometry_config"]["internal_membrane_relative_radii"][1]
            species_df = pd.read_csv(os.path.join(combined_folder, "species.csv"))
            X_external_concentration = species_df.loc[
                (species_df["name"] == "X"),
                "external_concentration"].item()
            enzymatic_reactions_df = pd.read_csv(os.path.join(combined_folder, "enzymatic_reactions.csv"))
            michaelis_menten_constant = enzymatic_reactions_df.loc[
                (enzymatic_reactions_df["enzyme"] == "B"),
                "k_M"].item()

            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            if os.path.isfile(fluxes_file):
                flux = load_json(fluxes_file)["Z"]
            else:
                flux = None

            average_Y_concentration = calculate_average_concentration_of_Y_in_inner_compartment_with_enzyme_B(combined_folder)
            data[index] = (pos_membrane, flux, X_external_concentration, average_Y_concentration, michaelis_menten_constant)
    return data


def plot_data(folder):
    cmap = plt.cm.viridis  # or any colormap you like
    
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=['pos_membrane', 'flux', 'X_external_concentration', 'average_Y_concentration', 'michaelis_menten_constant'])
    
    concentrations = np.sort(df['X_external_concentration'].unique())
    unique_k = sorted(df['michaelis_menten_constant'].unique())
    fig, ax = plt.subplots(len(concentrations), 2, figsize = (8,3*len(concentrations)))

    for row, X_external_concentration in enumerate(concentrations):
        current_df = df[df["X_external_concentration"]==X_external_concentration]
        for i, (k_val, group) in enumerate(current_df.groupby('michaelis_menten_constant')):
            color = cmap(i / len(unique_k))  # pick consistent color
            ax[row][0].scatter(
                group["pos_membrane"],
                group["flux"],
                label=str(k_val),
                norm=LogNorm(),
                color=color
            )
            ax[row][1].scatter(
                group["average_Y_concentration"],
                group["flux"],
                label=str(k_val),
                norm=LogNorm(),
                color=color
            )
            ax[row][1].axvline(x = k_val, color=color)
        ax[row][0].set_ylabel("flux")
        ax[row][0].set_title(f"X_ext = {X_external_concentration}")

    ax[-1][0].legend(title="K_m")
    ax[-1][0].set_xlabel("relative radius of most inner membrane r*/R")
    ax[-1][1].set_xlabel("average Y concentration where B enzyme located")
    #ax[1].axvline(x = 0.025)
    #ax[0].set_yscale("log")


    fig.tight_layout()
    fig.savefig(os.path.join(folder, "result.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)
    # python data/04d_enzymaticXtoY_enzymaticYtoZ_3InnerBoundaries_modifyingPositionOf1EnzymeRegion_modifyingKm_AOutside/analysis.py data/04d_enzymaticXtoY_enzymaticYtoZ_3InnerBoundaries_modifyingPositionOf1EnzymeRegion_modifyingKm_AOutside