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
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})

def get_data(folder):
    data = {}

    for folder_name in os.listdir(folder):
        match = re.match(r'^combined_(\d{6})$', folder_name)
        combined_folder = os.path.join(folder, folder_name)
        if match:
            index = match.group(1)  # keeps it as a string, preserving leading zeros
            geometry = read_yaml_file(os.path.join(combined_folder, "parameters_geometry.yaml"))
            species_df = pd.read_csv(os.path.join(combined_folder, "species.csv"))
            X_permeability = species_df.loc[(species_df["name"] == "X"), "permeability_constant"].item()
            try:
                relative_membrane_value = geometry["geometry_config"]["internal_membrane_relative_radii"][0]
            except: # in case no inner boundary
                relative_membrane_value = None 
            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            if os.path.isfile(fluxes_file):
                y_flux = load_json(fluxes_file)["Y"]
            else:
                y_flux = None
            data[index] = (relative_membrane_value, X_permeability, y_flux)

    return data


def plot_data(folder):
    fig, ax = plt.subplots(2, 1, figsize = (4,6), sharex=False)
    data = get_data(folder)

    df = pd.DataFrame(data.values(), columns=['relative_membrane_value', 'X_permeability','flux'])
    for X_permeability, group in df.groupby('X_permeability'):
        # find flux when there is no inner membrane
        flux_when_no_inner_boundary = group.loc[df["relative_membrane_value"].isna(), 'flux'].item()
        # plot all others
        group = group[group["relative_membrane_value"]!=None]
        group = group.sort_values("relative_membrane_value")
        ax[0].plot(
            group["relative_membrane_value"],
            group["flux"],
            label='{:.0e}'.format(X_permeability),
        )
        ax[1].plot(
            group["relative_membrane_value"],
            group["flux"]/flux_when_no_inner_boundary,
            label='{:.0e}'.format(X_permeability),
        )
    ax[0].set_yscale("log")
    ax[1].set_xlabel(r"relative radius of innermost membrane $r_1/R$")
    ax[0].set_ylabel("species Y outward flux \n with 1 inner membrane \n"+ r"$j_X / \mathrm{mol} \cdot \mathrm{s}^{-1}$")
    ax[1].set_ylabel("species Y outward flux with 1 inner \n membrane divided by flux \n without membrane \n" + r"$j_X(r_1=r_1) / j_X(r_1=0)$")
    ax[0].legend(title=r"permeability $p_X$")
    
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "result.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)
    # python data/01d_spontaneousXtoY_1InnerBoundary_modifyingInnerRadius_modifyingXPermeability/analysis.py data/01d_spontaneousXtoY_1InnerBoundary_modifyingInnerRadius_modifyingXPermeability