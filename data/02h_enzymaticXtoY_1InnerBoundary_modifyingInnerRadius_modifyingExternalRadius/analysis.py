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


def get_data(folder):
    data = {}
    for folder_name in os.listdir(folder):
        match = re.match(r'^combined_(\d{6})$', folder_name)
        combined_folder = os.path.join(folder, folder_name)
        if match:
            index = match.group(1)  # keeps it as a string, preserving leading zeros
            geometry = read_yaml_file(os.path.join(combined_folder, "parameters_geometry.yaml"))
            inner_radius = geometry["geometry_config"]["internal_membrane_relative_radii"][0]
            external_radius = geometry["geometry_config"]["outer_membrane_radius"]
            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            if os.path.isfile(fluxes_file):
                flux = load_json(fluxes_file)["Y"]
            else:
                flux = None
            data[index] = (inner_radius, flux, external_radius, index)
    return data


def plot_data(folder):
    fig, ax = plt.subplots(1, 2, figsize = (8,3))
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=['inner_radius', 'flux', 'external_radius','index'])
    external_radii = np.sort(df['external_radius'].unique())
    smallest_inner_radius = np.min(df["inner_radius"])
    for radius in external_radii:
        current_df = df[df["external_radius"]==radius]
        current_df = current_df.sort_values("inner_radius")
        ax[0].scatter(current_df["inner_radius"], current_df["flux"], label=radius)
        
        flux_for_innermost_inner_radius = current_df.loc[df["inner_radius"]==smallest_inner_radius, 'flux'].item()
        ax[1].scatter(current_df["inner_radius"], current_df["flux"]/flux_for_innermost_inner_radius, label=radius, alpha = 0.5)
        ax[1].plot(current_df["inner_radius"], current_df["flux"]/flux_for_innermost_inner_radius, label=radius)
    ax[0].set_xlabel("relative membrane radius")
    ax[0].set_ylabel("flux")
    ax[1].set_xlabel("relative membrane radius")
    ax[1].set_ylabel("flux / flux at smallest inner radius")
    ax[0].legend(title = "R")
    ax[0].set_yscale("log")
    ax[1].axhline(1, c = "k", ls = ":")
    #ax[1].set_yscale("log")
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "fluxes.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)
    # python data/02h_enzymaticXtoY_1InnerBoundary_modifyingInnerRadius_modifyingExternalRadius/analysis.py data/02h_enzymaticXtoY_1InnerBoundary_modifyingInnerRadius_modifyingExternalRadius