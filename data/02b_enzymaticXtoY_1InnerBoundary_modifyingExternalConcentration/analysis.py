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
from pathlib import Path
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})

def calculate_average_concentration_of_X_in_inner_compartment(combined_folder):
    system_geometry = load_json(os.path.join(combined_folder, "system_geometry_for_convergence.json"))
    concentrations = load_json(os.path.join(combined_folder, ".species_steady_state_concentrations.json"))
    
    radii_in_region_0 = system_geometry["geometry_config"]["mesh_points_in_regions"][0]
    X_concentrations_in_region_0 = []
    for radius_idx in range(len(radii_in_region_0)):
        X_concentrations_in_region_0.append(concentrations[0][radius_idx]["X"])
    
    concentration = 0
    for i in range(len(radii_in_region_0) - 1):
        r0 = radii_in_region_0[i]
        r1 = radii_in_region_0[i+1]
        c0 = X_concentrations_in_region_0[i]
        c1 = X_concentrations_in_region_0[i+1]
        
        # Volume-weighted average concentration in shell (assuming linear c(r))
        c_avg = (c0 * (3*r0**2 + 2*r0*r1 + r1**2) + c1 * (r0**2 + 2*r0*r1 + 3*r1**2)) \
                / (4 * (r0**2 + r0*r1 + r1**2))
        
        shell_volume = r1**3 - r0**3  # 4π/3 cancels with denominator
        concentration += c_avg * shell_volume

    concentration /= (radii_in_region_0[-1]**3 - radii_in_region_0[0]**3)
    return concentration

def get_data(folder):
    data = {}

    for folder_name in os.listdir(folder):
        match = re.match(r'^combined_(\d{6})$', folder_name)
        combined_folder = os.path.join(folder, folder_name)
        if match:
            index = match.group(1)  # keeps it as a string, preserving leading zeros
            species_df = pd.read_csv(os.path.join(combined_folder, "species.csv"))
            x_axis_value = species_df.loc[(species_df["name"] == "X"), "external_concentration"].item()
            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            if os.path.isfile(fluxes_file):
                flux = load_json(fluxes_file)["Y"]
                y_axis_value_1 = flux
                convergence_file = Path(combined_folder) / "concentration_convergence"
                if os.path.isfile(convergence_file):
                    y_axis_value_2 = calculate_average_concentration_of_X_in_inner_compartment(combined_folder) #[0][0]["X"]# this should be the average value but this gives a rough idea
                else:
                    y_axis_value_2 = None
                y_axis_value_3 = flux / x_axis_value
            else:
                y_axis_value_1 = None
                y_axis_value_2 = None
                y_axis_value_3 = None


            data[index] = (x_axis_value, y_axis_value_1, y_axis_value_2, y_axis_value_3)

    return data


def plot_data(folder):
    fig, ax = plt.subplots(3, 1, figsize = (3.5,8))
    data = get_data(folder)
    
    df = pd.DataFrame(data.values(), columns=['x', 'flux', 'X_conc_at_0', 'flux_per_X'])
    df = df.sort_values("x")
    ax[0].plot(df["x"], df["flux"])
    ax[1].plot(df["x"], df["X_conc_at_0"])
    ax[1].axvline(0.025, label = "MM constant", ls = ":")
    ax[2].plot(df["x"], df["flux_per_X"])
    ax[1].legend()

    # Set log scale on whichever axes need it
    ax[2].set_xlabel(r"external concentration of X / mol $\cdot$ m$^{-3}$")
    ax[0].set_ylabel(r"flux of Y / mol $\cdot$ s$^{-1}$")
    ax[1].set_ylabel(r'average concentration of X'+ "\n"+ r'in innermost region /'+ "\n"+ r"mol $\cdot$ m$^{-3}$")
    ax[2].set_ylabel(r'flux of Y divided by'+ "\n"+ r'external concentration of X')
    ax[0].set_box_aspect(1)
    ax[1].set_box_aspect(1)
    ax[2].set_box_aspect(1)

    ax[0].set_xscale('log')
    ax[1].set_xscale('log')
    ax[2].set_xscale('log')
    ax[1].set_yscale("log")

    fig.tight_layout()
    fig.savefig(os.path.join(folder, "fluxes.png"), dpi = 300)#, bbox_inches="tight")

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)
    # python data/02b_enzymaticXtoY_1InnerBoundary_modifyingExternalConcentration/analysis.py data/02b_enzymaticXtoY_1InnerBoundary_modifyingExternalConcentration