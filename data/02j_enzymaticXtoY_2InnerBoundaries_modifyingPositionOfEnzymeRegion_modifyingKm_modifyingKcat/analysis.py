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
            catalytic_rate = enzymatic_reactions_df.loc[
                (enzymatic_reactions_df["enzyme"] == "A"),
                "k_cat"].item()
            concentration = calculate_average_concentration_of_X_in_inner_compartment_with_enzyme_A(combined_folder)
            data[index] = (internal_radius, flux, michaelis_menten_constant, catalytic_rate, concentration)
    return data


def plot_data(folder):
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=['internal_radius', 'flux', 'michaelis_menten_constant', 'catalytic_rate', 'concentration'])
    
    internal_radii = np.sort(df['internal_radius'].unique())
    fig, ax = plt.subplots(len(internal_radii), 2, figsize = (4,3*len(internal_radii)), gridspec_kw={'width_ratios': [1, 0.1]})

    ############
    # Plot fluxes
    ############
    for internal_radius_idx, internal_radius in enumerate(internal_radii):
        current_df = df[df['internal_radius']==internal_radius]
        z_grid = current_df.pivot(index='catalytic_rate', columns='michaelis_menten_constant', values='flux')
        x_points = z_grid.columns.values
        y_points = z_grid.index.values
        z_points = z_grid.values  # 2D array of shape (len(y), len(x))
        
        mesh0 = ax[internal_radius_idx][0].pcolormesh(x_points, y_points, z_points, cmap='viridis', shading='auto',
                                #norm = Normalize(vmin=0, vmax=1)
                                #norm=LogNorm(vmin=np.nanmin(z_points), vmax=np.nanmax(z_points))
        )
        ax[internal_radius_idx][0].scatter(current_df['michaelis_menten_constant'], current_df['catalytic_rate'], color='red', s=5, zorder=5)
        #for _, point in df.iterrows():
        #    ax[0].annotate(f"{point['index'].lstrip('0')}", (point['x'], point['y']))

        fig.colorbar(mesh0, cax=ax[internal_radius_idx][1], label='flux')
        ax[internal_radius_idx][0].set_title(f"r* = {internal_radius}")
        ax[internal_radius_idx][0].set_xscale("log")
        ax[internal_radius_idx][0].set_yscale("log")
        ax[internal_radius_idx][1].set_box_aspect(10)
        ax[internal_radius_idx][0].set_box_aspect(1)

    # Set log scale on whichever axes need it
    ax[0][0].set_xlabel("K_m")
    ax[0][0].set_ylabel("k_cat")
    

    #mesh1 = ax[1][0].pcolormesh(x_points, y_points, z_points, cmap='viridis', shading='auto',
    #    norm = Normalize(vmin=3.5e-7, vmax=3.90e-7)
    #    #norm=LogNorm(vmin=np.nanmin(z_points), vmax=np.nanmax(z_points))
    #)
    #fig.colorbar(mesh1, cax=ax[1][1], label='flux')
    #ax[1][1].set_box_aspect(10)
    #ax[1][0].set_box_aspect(1)
    
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "fluxes.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)
    # python data/02j_enzymaticXtoY_2InnerBoundaries_modifyingPositionOfEnzymeRegion_modifyingKm_modifyingKcat/analysis.py data/02j_enzymaticXtoY_2InnerBoundaries_modifyingPositionOfEnzymeRegion_modifyingKm_modifyingKcat