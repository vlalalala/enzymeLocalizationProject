import sys
import matplotlib.pyplot as plt
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src'))
from auxiliary_functions_using_standard_library import load_json
import re
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
from auxiliary_functions import read_yaml_file
import ast
import matplotlib.image as mpimg
import numpy as np
from find_matching_parameter_value_combinations import filter_combined_folders
from pathlib import Path

def get_data(folder):
    data = {}
    for folder_name in os.listdir(folder):
        match = re.match(r'^combined_(\d{6})$', folder_name)
        combined_folder = os.path.join(folder, folder_name)
        if match:
            index = match.group(1)  # keeps it as a string, preserving leading zeros
            geometry = read_yaml_file(os.path.join(combined_folder, "parameters_geometry.yaml"))
            inner_radius = geometry["geometry_config"]["internal_membrane_relative_radii"][0]
            outer_inner_radius = geometry["geometry_config"]["internal_membrane_relative_radii"][1]

            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            try:
                if os.path.isfile(fluxes_file):
                    flux = load_json(fluxes_file)["Y"]
                else:
                    flux = None
            except:
                flux = None
            
            enzymatic_reactions_df = pd.read_csv(os.path.join(combined_folder, "enzymatic_reactions.csv"))
            michaelis_menten_constant = enzymatic_reactions_df.loc[
                (enzymatic_reactions_df["enzyme"] == "A"),
                "k_M"].item()
            catalytic_rate = enzymatic_reactions_df.loc[
                (enzymatic_reactions_df["enzyme"] == "A"),
                "k_cat"].item()
            if flux is None:
                print(index, inner_radius, outer_inner_radius, catalytic_rate, michaelis_menten_constant )
            data[index] = (#index,
                           inner_radius, outer_inner_radius, flux,
                           #k,
                           michaelis_menten_constant, catalytic_rate)
    return data

def find_optimal_radius(df):
    """Returns a dataframe with columns Km, Kcat, k, best_inner_radius, best_flux"""
    def get_best(g):
        idx = g['flux'].idxmax()
        return pd.Series({
            'best_inner_radius': g.loc[idx, 'inner_radius'],
            'best_flux': g.loc[idx, 'flux']
        })

    best_radius = (
        df.groupby(['Km', 'Kcat', 'k'])
        .apply(get_best)
        .reset_index()
    )
    return best_radius

def plot_data(folder):
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=[
        #'index',
        'inner_radius', "outer_inner_radius", 'flux',
        "Km", "Kcat"])
    #print(df)
    best_radius_df = find_optimal_radius(df)
    #print(best_radius_df)
    fig, ax = plt.subplots(1, 2, figsize = (4,3 * 1), gridspec_kw={'width_ratios': [1, 0.1]})

    pivot = best_radius_df.pivot(index='Km', columns='Kcat', values='best_inner_radius')
    Kcat_vals = pivot.columns.values
    Km_vals = pivot.index.values
    Z = pivot.values  # 2D array
    #print(Z.shape)

    mesh = ax[0].pcolormesh(
        Kcat_vals, Km_vals, Z,
        cmap='viridis', shading='auto'
    )
    ax[0].set_xlabel("k cat")
    ax[0].set_ylabel("Km")
    #mesh = ax[Km_idx][0].pcolormesh(group["Kcat"], group["k"], group["best_inner_radius"],
    #    cmap='viridis', shading='auto',
        #norm = Normalize(vmin=0, vmax=1)
        #norm=LogNorm(vmin=np.nanmin(z_points), vmax=np.nanmax(z_points)
    #)
    ax[0].scatter(best_radius_df["Kcat"], best_radius_df['Km'], color='red', s=5, zorder=5)
    #ax[0].set_title(f"Km = {Km_value}")
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    #for _, point in df.iterrows():
    #    ax[0].annotate(f"{point['index'].lstrip('0')}", (point['x'], point['y']))
    fig.colorbar(mesh, cax=ax[1], label="best relative radius \n of most inner membrane r*/R")
    ax[0].set_box_aspect(1)
    ax[1].set_box_aspect(10)
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "result.png"), dpi = 300)

def plot_steady_states(folder):
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=[
        'inner_radius', "outer_inner_radius", 'flux',
        "Km", "Kcat"])
    # one file per Km
    Kcat = np.sort(df["Kcat"].unique())
    inner_radii = np.sort(df["inner_radius"].unique())
    Km = np.sort(df["Km"].unique())
    for Km_value in Km:
        fig, ax = plt.subplots(len(Kcat), len(inner_radii), figsize = (4*len(inner_radii), 3*len(Kcat)))
        fig.suptitle(f"Km = {Km_value}")
        for Kcat_idx, Kcat_value in enumerate(Kcat):
            for inner_radius_idx, inner_radius in enumerate(inner_radii):
                current_df = df[(df["Kcat"]==Kcat_value)&(df["inner_radius"]==inner_radius)&(df['Km']==Km_value)]
                outer_inner_radius = list(current_df["outer_inner_radius"])[0]
                combinations = filter_combined_folders(
                    combined_root=folder,
                    criteria_yaml={
                        "options_parameters_geometry": {
                            "geometry_config": {
                                "internal_membrane_relative_radii": [inner_radius, outer_inner_radius]
                            }
                        }
                    },
                    criteria_csv=
                        {
                            "options_enzymatic_reactions": {"enzyme": {"A": {"k_cat": f"{Kcat_value}", "k_M": f"{Km_value}"}}}
                        }
                )
                if len(combinations)!=1:
                    raise ValueError(combinations)
                combination = combinations[0]#only one available
                file_to_plot = combination / "solver_iteration_data" / "interpolation_iteration_nr_0_final_concentrations.png"
                if os.path.isfile(file_to_plot):
                    img = mpimg.imread(str(file_to_plot))
                ax[Kcat_idx][inner_radius_idx].imshow(img)
                ax[Kcat_idx][inner_radius_idx].axis("off")
                ax[Kcat_idx][inner_radius_idx].set_title(f"Kcat = {Kcat_value}, \n inner radius = {inner_radius} \n {Path(combination).name}")
        fig.savefig(os.path.join(folder, f"complete_steadyStates_Km_{Km_value}.png"), dpi = 300)
            


if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)
    #plot_steady_states(FOLDER_TO_SOLVE)
    # python data/02l_enzymaticXtoY_2InnerBoundaries_modifyingPositionOfEnzymeRegion_modifyingKm_modifyingKcat/analysis.py data/02l_enzymaticXtoY_2InnerBoundaries_modifyingPositionOfEnzymeRegion_modifyingKm_modifyingKcat

# The smaller Km, the steeper the fall in the concentration of the reactant is within the shell where there is enzyme.