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
            if os.path.isfile(fluxes_file):
                flux = load_json(fluxes_file)["Z"]
            else:
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
    df = pd.DataFrame(data.values(), columns=['inner_radius', 'flux', 'k', "Km", "Kcat"])
    best_radius_df = find_optimal_radius(df)
    # only one value of Km here
    Km = list(best_radius_df["Km"].unique())[0]
    
    fig, ax = plt.subplots(2, 2, figsize = (5,3 * 2), gridspec_kw={'width_ratios': [1, 0.1]})

    pivot = best_radius_df.pivot(index='k', columns='Kcat', values='best_inner_radius')
    Kcat_vals = pivot.columns.values
    k_vals = pivot.index.values
    Z = pivot.values  # 2D array

    mesh = ax[0][0].pcolormesh(
        Kcat_vals, k_vals, Z,
        cmap='viridis', shading='auto'
    )

    ax[0][0].set_xlabel("k cat")
    ax[0][0].set_ylabel("k")
    ax[0][0].scatter(best_radius_df["Kcat"], best_radius_df["k"], color='red', s=5, zorder=5)
    ax[0][0].set_title(f"Km = {Km}")
    ax[0][0].set_xscale("log")
    ax[0][0].set_yscale("log")
    fig.colorbar(mesh, cax=ax[0][1], label="best relative radius \n of most inner membrane r*/R")
    ax[0][0].set_box_aspect(1)
    ax[0][1].set_box_aspect(10)

    pivot = best_radius_df.pivot(index='k', columns='Kcat', values='best_flux')
    Kcat_vals = pivot.columns.values
    k_vals = pivot.index.values
    Z = pivot.values  # 2D array

    mesh = ax[1][0].pcolormesh(
        Kcat_vals, k_vals, Z,
        cmap='viridis', shading='auto',
        norm=LogNorm(vmin=Z.min(), vmax=Z.max()),
    )

    ax[1][0].set_xlabel("k cat")
    ax[1][0].set_ylabel("k")
    ax[1][0].scatter(best_radius_df["Kcat"], best_radius_df["k"], color='red', s=5, zorder=5)
    ax[1][0].set_xscale("log")
    ax[1][0].set_yscale("log")

    #for _, point in df.iterrows():
    #    ax[0].annotate(f"{point['index'].lstrip('0')}", (point['x'], point['y']))
    fig.colorbar(mesh, cax=ax[1][1], label="best flux")
    ax[1][0].set_box_aspect(1)
    ax[1][1].set_box_aspect(10)

    fig.tight_layout()
    fig.savefig(os.path.join(folder, "result.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)
    # python data/05d_spontaneousXtoY_enzymaticYtoZ_2InnerBoundaries_modifyingPositionOfEnzymeRegion_modifyingKcat_lowPermeability/analysis.py data/05d_spontaneousXtoY_enzymaticYtoZ_2InnerBoundaries_modifyingPositionOfEnzymeRegion_modifyingKcat_lowPermeability