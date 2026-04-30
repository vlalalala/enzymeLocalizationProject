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
            x_axis_value = allocationA[1]            
            y_axis_value = allocationB[1]
            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            if os.path.isfile(fluxes_file):
                z_axis_value = load_json(fluxes_file)["Z"]
            else:
                z_axis_value = None
            data[index] = (x_axis_value, y_axis_value, z_axis_value)

    return data


def plot_data(folder):
    fig, ax = plt.subplots(1, 2, figsize = (5,3), gridspec_kw={'width_ratios': [1, 0.1]})
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=['x', 'y', 'flux'])
    
    ############
    # Plot fluxes
    ############
    z_grid = df.pivot(index='y', columns='x', values='flux')
    x_points = z_grid.columns.values
    y_points = z_grid.index.values
    z_points = z_grid.values  # 2D array of shape (len(y), len(x))
    
    mesh = ax[0].pcolormesh(x_points, y_points, z_points, cmap='viridis', shading='auto',
                            #norm = Normalize(vmin=0, vmax=1)
                            #norm=LogNorm(vmin=np.nanmin(z_points), vmax=np.nanmax(z_points)
    )
    ax[0].scatter(df['x'], df['y'], color='red', s=5, zorder=5)
    fig.colorbar(mesh, cax=ax[1], label='flux')

    # Set log scale on whichever axes need it
    ax[0].set_xlabel("proportion of \n enzyme A quantity \n in outermost region")
    ax[0].set_ylabel("proportion of \n enzyme B quantity \n in outermost region")
    ax[1].set_box_aspect(10)
    ax[0].set_box_aspect(1)
    
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "fluxes.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)