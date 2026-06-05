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
            if os.path.isfile(fluxes_file):
                flux = load_json(fluxes_file)["Z"]
            else:
                flux = None
                print(f"could not find flux file for {fluxes_file}")
            enzymatic_reactions_df = pd.read_csv(os.path.join(combined_folder, "enzymatic_reactions.csv"))
            catalytic_rate = enzymatic_reactions_df.loc[
                (enzymatic_reactions_df["enzyme"] == "A"),
                "k_cat"].item()
            data[index] = (inner_radius, outer_inner_radius, flux, catalytic_rate)
    return data

def plot_data(folder):
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=['inner_radius', 'outer_inner_radius', 'flux', "Kcat"])
    # only one value of Km here
    fig, ax = plt.subplots(1, 1, figsize = (4,3))
    Kcat = np.sort(df["Kcat"].unique())
    for Kcat_idx, Kcat_value in enumerate(Kcat):
        current_df = df[df["Kcat"]==Kcat_value]
        current_df = current_df.sort_values("inner_radius")
        ax.plot(current_df["inner_radius"], current_df["flux"], label = "{:.2E}".format(Kcat_value), ls = ":", alpha = 0.2)
    
    ax.set_xlabel("inner radius")
    ax.set_ylabel("flux")
    ax.legend(title="catalytic rate",loc='center left', bbox_to_anchor=(1, 0.5))
    #ax[0][0].set_yscale("log")

    fig.tight_layout()
    fig.savefig(os.path.join(folder, "result.png"), dpi = 300)

def plot_steady_states(folder):
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=['inner_radius','outer_inner_radius', 'flux', "Kcat"])
    # only one value of Km here
    Kcat = np.sort(df["Kcat"].unique())
    inner_radii = np.sort(df["inner_radius"].unique())
    fig, ax = plt.subplots(len(Kcat), len(inner_radii), figsize = (4*len(inner_radii), 3*len(Kcat)))
    for Kcat_idx, Kcat_value in enumerate(Kcat):
        for inner_radius_idx, inner_radius in enumerate(inner_radii):
            current_df = df[(df["Kcat"]==Kcat_value)&(df["inner_radius"]==inner_radius)]
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
                        "options_enzymatic_reactions": {"enzyme": {"A": {"k_cat": f"{Kcat_value}"}}}
                    }
            )
            combination = combinations[0]#only one available
            file_to_plot = combination / "solver_iteration_data" / "interpolation_iteration_nr_0_final_concentrations.png"
            if os.path.isfile(file_to_plot):
                img = mpimg.imread(str(file_to_plot))
            ax[Kcat_idx][inner_radius_idx].imshow(img)
            ax[Kcat_idx][inner_radius_idx].axis("off")
            ax[Kcat_idx][inner_radius_idx].set_title(f"Kcat = {Kcat_value}, \n inner radius = {inner_radius}")
    fig.savefig(os.path.join(folder, "complete_steadyStates.png"), dpi = 300)
            

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    #plot_steady_states(FOLDER_TO_SOLVE)
    plot_data(FOLDER_TO_SOLVE)
    # python data/05g_spontaneousXtoY_enzymaticYtoZ_2InnerBoundaries_modifyingPositionOfEnzymeRegion_modifyingKcat_lowPermeability/analysis.py data/05g_spontaneousXtoY_enzymaticYtoZ_2InnerBoundaries_modifyingPositionOfEnzymeRegion_modifyingKcat_lowPermeability