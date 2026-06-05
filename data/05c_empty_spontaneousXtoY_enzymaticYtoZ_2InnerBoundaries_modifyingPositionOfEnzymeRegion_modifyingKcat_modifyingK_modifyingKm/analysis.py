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

def get_data(folder):
    data = {}
    for folder_name in os.listdir(folder):
        match = re.match(r'^combined_(\d{6})$', folder_name)
        combined_folder = os.path.join(folder, folder_name)
        if match:
            index = match.group(1)  # keeps it as a string, preserving leading zeros
            geometry = read_yaml_file(os.path.join(combined_folder, "parameters_geometry.yaml"))
            inner_radius = geometry["geometry_config"]["internal_membrane_relative_radii"][0]
            spontaneous_reactions_df = pd.read_csv(os.path.join(combined_folder, "spontaneous_reactions.csv"))
            k = spontaneous_reactions_df.loc[
                (spontaneous_reactions_df["start_species"] == "X"),
                "k"].item()
            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            try:
                if os.path.isfile(fluxes_file):
                    flux = load_json(fluxes_file)["Z"]
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
            data[index] = (inner_radius, flux, k, michaelis_menten_constant, catalytic_rate)
    return data

def find_optimal_radius(df):
    """Returns a dataframe with columns Km, Kcat, k, best_inner_radius"""
    try:
        best_radius = (
            df.groupby(['Km', 'Kcat', 'k'])
            .apply(lambda g: g.loc[g['flux'].idxmax(), 'inner_radius'])
            .reset_index(name='best_inner_radius')
        )
    except:
        best_radius = pd.DataFrame(columns=['Km', 'Kcat', 'k', 'best_inner_radius'])
    return best_radius

def plot_data(folder):
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=['inner_radius', 'flux', 'k', "Km", "Kcat"])
    best_radius_df = find_optimal_radius(df)
    Km_num = len(df["Km"].unique())
    
    fig, ax = plt.subplots(Km_num, 2, figsize = (4,3 * Km_num), gridspec_kw={'width_ratios': [1, 0.1]})
    if Km_num == 1:
        ax = [ax] # so that indexing works

    for Km_idx, (Km_value, group) in enumerate(best_radius_df.groupby('Km')):
        pivot = group.pivot(index='k', columns='Kcat', values='best_inner_radius')
        Kcat_vals = pivot.columns.values
        k_vals = pivot.index.values
        Z = pivot.values  # 2D array

        mesh = ax[Km_idx][0].pcolormesh(
            Kcat_vals, k_vals, Z,
            cmap='viridis', shading='auto'
    )
        ax[Km_idx][0].set_xlabel("k cat")
        ax[Km_idx][0].set_ylabel("k")
        #mesh = ax[Km_idx][0].pcolormesh(group["Kcat"], group["k"], group["best_inner_radius"],
        #    cmap='viridis', shading='auto',
            #norm = Normalize(vmin=0, vmax=1)
            #norm=LogNorm(vmin=np.nanmin(z_points), vmax=np.nanmax(z_points)
        #)
        ax[Km_idx][0].scatter(group["Kcat"], group["k"], color='red', s=5, zorder=5)
        ax[Km_idx][0].set_title(f"Km = {Km_value}")
        ax[Km_idx][0].set_xscale("log")
        ax[Km_idx][0].set_yscale("log")
        #for _, point in df.iterrows():
        #    ax[0].annotate(f"{point['index'].lstrip('0')}", (point['x'], point['y']))
        fig.colorbar(mesh, cax=ax[Km_idx][1], label="best relative radius \n of most inner membrane r*/R")
        ax[Km_idx][0].set_box_aspect(1)
        ax[Km_idx][1].set_box_aspect(10)
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "result.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)
    # python data/05c_spontaneousXtoY_enzymaticYtoZ_2InnerBoundaries_modifyingPositionOfEnzymeRegion_modifyingKcat_modifyingK_modifyingKm/analysis.py data/05c_spontaneousXtoY_enzymaticYtoZ_2InnerBoundaries_modifyingPositionOfEnzymeRegion_modifyingKcat_modifyingK_modifyingKm