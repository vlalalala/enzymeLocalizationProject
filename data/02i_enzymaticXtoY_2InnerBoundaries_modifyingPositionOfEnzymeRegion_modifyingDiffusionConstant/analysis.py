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
            species_df = pd.read_csv(os.path.join(combined_folder, "species.csv"))
            diffusion = species_df.loc[(species_df["name"] == "X"), "diffusion_constant"].item()
            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            if os.path.isfile(fluxes_file):
                flux = load_json(fluxes_file)["Y"]
            else:
                flux = None
            data[index] = (inner_radius, diffusion, flux)
    return data


def plot_data(folder):
    fig, ax = plt.subplots(1, 1, figsize = (4,3))
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=['inner_radius','diffusion', 'flux'])
    diffusion_constants = df["diffusion"].unique()
    for diffusion_constant in diffusion_constants:
        current_df = df[df["diffusion"]==diffusion_constant]
        current_df = current_df.sort_values("inner_radius")
        ax.plot(current_df["inner_radius"], current_df["flux"], label = diffusion_constant)
        ax.scatter(current_df["inner_radius"], current_df["flux"], alpha = 0.5)
    ax.set_xlabel("relative radius of most inner membrane r*/R")
    ax.set_ylabel("flux")
    ax.legend(title = "D")
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "result.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)
    # python data/02i_enzymaticXtoY_2InnerBoundaries_modifyingPositionOfEnzymeRegion_modifyingDiffusionConstant/analysis.py data/02i_enzymaticXtoY_2InnerBoundaries_modifyingPositionOfEnzymeRegion_modifyingDiffusionConstant