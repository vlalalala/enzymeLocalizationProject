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
            data[index] = (inner_radius, flux, k)
    return data


def plot_data(folder):
    fig, ax = plt.subplots(2, 1, figsize = (4,3))
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=['inner_radius', 'flux', 'k'])
    smallest_inner_radius = np.min(df["inner_radius"])
    for k_val, group in df.groupby('k'):
        group = group.sort_values("inner_radius")
        ax[0].plot(
            group["inner_radius"],
            group["flux"],
            label=str(k_val)
        )
        flux_for_innermost_inner_radius = group.loc[df["inner_radius"]==smallest_inner_radius, 'flux'].item()
        ax[1].scatter(group["inner_radius"], group["flux"]/flux_for_innermost_inner_radius, alpha = 0.5)
        ax[1].plot(group["inner_radius"], group["flux"]/flux_for_innermost_inner_radius)
    ax[0].legend(title="k")
    ax[0].set_xlabel("relative radius of most inner membrane r*/R")
    ax[0].set_ylabel("flux")
    ax[0].legend(title = "k")
    ax[1].axhline(1, c = "k", alpha = 0.2, ls = ":")
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "result.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)
    # python data/05_spontaneousXtoY_enzymaticYtoZ_2InnerBoundaries_modifyingPositionOfEnzymeRegion/analysis.py data/05_spontaneousXtoY_enzymaticYtoZ_2InnerBoundaries_modifyingPositionOfEnzymeRegion