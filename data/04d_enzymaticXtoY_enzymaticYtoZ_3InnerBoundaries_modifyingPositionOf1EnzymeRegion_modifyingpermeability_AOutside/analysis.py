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
from decimal import Decimal
def calculate_average_concentration_of_Y_in_inner_compartment_with_enzyme_B(combined_folder):
    system_geometry_file = os.path.join(combined_folder, "system_geometry_for_convergence.json")
    concentrations_file = os.path.join(combined_folder, ".species_steady_state_concentrations.json")
    if os.path.isfile(system_geometry_file) and os.path.isfile(concentrations_file):
        system_geometry = load_json(system_geometry_file)
        concentrations = load_json(concentrations_file)
    else:
        return None

    radii_in_region_1 = system_geometry["geometry_config"]["mesh_points_in_regions"][1]
    Y_concentrations_in_region_1 = []
    for radius_idx in range(len(radii_in_region_1)):
        Y_concentrations_in_region_1.append(concentrations[1][radius_idx]["Y"])
    
    concentration = 0
    for i in range(len(radii_in_region_1) - 1):
        r0 = radii_in_region_1[i]
        r1 = radii_in_region_1[i+1]
        c0 = Y_concentrations_in_region_1[i]
        c1 = Y_concentrations_in_region_1[i+1]
        
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
            pos_membrane = geometry["geometry_config"]["internal_membrane_relative_radii"][1]
            species_df = pd.read_csv(os.path.join(combined_folder, "species.csv"))
            X_permeability = species_df.loc[
                (species_df["name"] == "X"),
                "permeability_constant"].item()
            Y_permeability = species_df.loc[
                (species_df["name"] == "Y"),
                "permeability_constant"].item()
            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            if os.path.isfile(fluxes_file):
                flux = load_json(fluxes_file)["Z"]
            else:
                flux = None
            concentration = calculate_average_concentration_of_Y_in_inner_compartment_with_enzyme_B(combined_folder)
            data[index] = (pos_membrane, flux, X_permeability, Y_permeability, concentration, index)
    return data


def plot_data(folder):
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=['pos_membrane', 'flux', 'X_permeability', 'Y_permeability', 'concentration', 'index'])
    
    X_permeabilities = np.sort(df['X_permeability'].unique())
    Y_permeabilities = np.sort(df['Y_permeability'].unique())
    
    fig, ax = plt.subplots(len(X_permeabilities), len(Y_permeabilities), figsize = (4*len(X_permeabilities),3*len(Y_permeabilities)))

    for row, X_permeability in enumerate(X_permeabilities):
        for column, Y_permeability in enumerate(Y_permeabilities):
            current_df = df[(df["X_permeability"]==X_permeability) & (df["Y_permeability"]==Y_permeability)]
            ax[row][column].scatter(
                    current_df["pos_membrane"],
                    current_df["flux"],
                    norm=LogNorm(),
                    c = current_df["concentration"]
                )
            for _, point in current_df.iterrows():
                ax[row][column].annotate("{:.1e}".format(point['concentration']) + f"\n {point['index'].lstrip('0')}", (point['pos_membrane'], point['flux']))
            ax[row][column].set_ylabel("flux")
            ax[row][column].set_title(f"X_permeability= {X_permeability}, \n Y_permeability= {Y_permeability}")

        #ax[][0].set_xlabel("relative radius of most inner membrane r*/R")
        #ax[-1][1].set_xlabel("average Y concentration where B enzyme located")
        #ax[1].axvline(x = 0.025)
        #ax[0].set_yscale("log")


    fig.tight_layout()
    fig.savefig(os.path.join(folder, "result.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)
    # python data/04d_enzymaticXtoY_enzymaticYtoZ_3InnerBoundaries_modifyingPositionOf1EnzymeRegion_modifyingpermeability_AOutside/analysis.py data/04d_enzymaticXtoY_enzymaticYtoZ_3InnerBoundaries_modifyingPositionOf1EnzymeRegion_modifyingpermeability_AOutside