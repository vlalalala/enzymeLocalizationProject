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
import numpy as np
from find_matching_parameter_value_combinations import filter_combined_folders
import matplotlib.image as mpimg

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
            spontaneous_reactions_df = pd.read_csv(os.path.join(combined_folder, "spontaneous_reactions.csv"))
            k = spontaneous_reactions_df.loc[
                (spontaneous_reactions_df["start_species"] == "X"),
                "k"].item()
            if os.path.isfile(fluxes_file):
                flux = load_json(fluxes_file)["Y"]
            else:
                flux = None
                print(f"could not find flux file for {fluxes_file}")
            data[index] = (inner_radius, outer_inner_radius, k, flux)
    return data

def plot_data(folder):
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=['inner_radius', 'outer_inner_radius', 'k', 'flux'])
    # only one value of Km here
    fig, ax = plt.subplots(2, 1, figsize = (5,6))
    smallest_inner_radius = np.min(df["inner_radius"])
    for k in np.sort(df["k"].unique()):
        current_df = df[df["k"]==k]
        current_df = current_df.sort_values("inner_radius")
        ax[0].plot(current_df["inner_radius"], current_df["flux"], alpha = 1, label = k)
        flux_for_innermost_inner_radius = current_df.loc[df["inner_radius"]==smallest_inner_radius, 'flux'].item()
        ax[1].scatter(current_df["inner_radius"], current_df["flux"]/flux_for_innermost_inner_radius, alpha = 0.5)
        ax[1].plot(current_df["inner_radius"], current_df["flux"]/flux_for_innermost_inner_radius)
    ax[0].set_xlabel("inner radius")
    ax[1].set_xlabel("inner radius")
    ax[0].set_ylabel("flux")
    ax[1].set_ylabel("flux")
    ax[0].set_yscale("log")
    ax[0].legend(title="k")
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "result.png"), dpi = 300)

def plot_steady_states(folder):
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=['inner_radius','outer_inner_radius', 'k', 'flux'])
    # only one value of Km here
    inner_radii = np.sort(df["inner_radius"].unique())
    k_values =np.sort(df["k"].unique())
    fig, ax = plt.subplots(len(k_values), len(inner_radii), figsize = (4*len(inner_radii), 3*len(k_values)))
    
    for k_idx, k in enumerate(k_values):
        for inner_radius_idx, inner_radius in enumerate(inner_radii):
            current_df = df[(df["inner_radius"]==inner_radius)]
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
                criteria_csv={
                                "options_spontaneous_reactions": {"start_species": {"X": {"k": f"{k}"}}}

                    
                    }
            )
            if len(combinations)!=1:
                raise ValueError("more than one combination found")
            combination = combinations[0]#only one available
            file_to_plot = combination / "solver_iteration_data" / "interpolation_iteration_nr_0_final_concentrations.png"
            if os.path.isfile(file_to_plot):
                img = mpimg.imread(str(file_to_plot))
            ax[k_idx][inner_radius_idx].imshow(img)
            ax[k_idx][inner_radius_idx].axis("off")
            ax[k_idx][inner_radius_idx].set_title(f"inner radius = {inner_radius}, k = {k}")
        fig.savefig(os.path.join(folder, "complete_steadyStates.png"), dpi = 300)
            

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    #plot_steady_states(FOLDER_TO_SOLVE)
    plot_data(FOLDER_TO_SOLVE)
    # python data/01c_spontaneousXtoY_2InnerBoundaries_modifyingPositionOfRegion/analysis.py data/01c_spontaneousXtoY_2InnerBoundaries_modifyingPositionOfRegion