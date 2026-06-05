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
            enzymes_df = pd.read_csv(os.path.join(combined_folder, "enzymes.csv"))
            allocationA_str = enzymes_df.loc[(enzymes_df["name"] == "A"), "allocation"].item()
            allocationA = ast.literal_eval(allocationA_str)
            allocationB_str = enzymes_df.loc[(enzymes_df["name"] == "B"), "allocation"].item()
            allocationB = ast.literal_eval(allocationB_str)
            allocationA_1 = allocationA[1]            
            allocationB_1 = allocationB[1]
            species_df = pd.read_csv(os.path.join(combined_folder, "species.csv"))
            X_external_concentration = species_df.loc[
                (species_df["name"] == "X"),
                "external_concentration"].item()
            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            if os.path.isfile(fluxes_file):
                Z_flux = load_json(fluxes_file)["Z"]
            else:
                Z_flux = None
            data[index] = (allocationA_1, allocationB_1, X_external_concentration, Z_flux)

    return data


def plot_data(folder):
    fig, ax = plt.subplots(1, 1, figsize = (5,3))
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=[
        'allocationA_1', 'allocationB_1', 'X_external_concentration', 'flux'])

    X_external_concentrations = sorted(df['X_external_concentration'].unique())
    for external_concentration in X_external_concentrations:
        current_df = df[df["X_external_concentration"]==external_concentration]
        current_df = current_df.sort_values("allocationA_1")
        ax.plot(current_df["allocationA_1"], current_df["flux"]/np.amax(current_df["flux"]),
                   label = external_concentration, alpha = 0.5)

    ax.set_xlabel("Proportion of enzyme A (X->Y) within penultimate section ")
    ax.set_ylabel("flux of Z")
    ax.legend(title="X ext")


    ############
    # Plot fluxes
    ############
    #z_grid = df.pivot(index='y', columns='x', values='flux')
    #x_points = z_grid.columns.values
    #y_points = z_grid.index.values
    #z_points = z_grid.values  # 2D array of shape (len(y), len(x))
    # 
    #mesh = ax[0].pcolormesh(x_points, y_points, z_points, cmap='viridis', shading='auto',
    #                        #norm = Normalize(vmin=0, vmax=1)
    #                        #norm=LogNorm(vmin=np.nanmin(z_points), vmax=np.nanmax(z_points)
    #)
    #ax[0].scatter(df['x'], df['y'], color='red', s=5, zorder=5)
    #fig.colorbar(mesh, cax=ax[1], label='flux')

    # Set log scale on whichever axes need it
    #ax[0].set_xlabel("proportion of \n enzyme A quantity \n in outermost region")
    #ax[0].set_ylabel("proportion of \n enzyme B quantity \n in outermost region")
    #ax[1].set_box_aspect(10)
    #ax[0].set_box_aspect(1)
    
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "result.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)

# python data/04k_enzymaticXtoY_enzymaticYtoZ_2InnerBoundaries_modifyingRatiosAandBwithinOutermostRegion/analysis.py data/04k_enzymaticXtoY_enzymaticYtoZ_2InnerBoundaries_modifyingRatiosAandBwithinOutermostRegion