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
            x_axis_value = geometry["geometry_config"]["internal_membrane_relative_radii"][0]
            enzymes_df = pd.read_csv(os.path.join(combined_folder, "enzymes.csv"))
            allocation_str = enzymes_df.loc[(enzymes_df["name"] == "A"), "allocation"].item()
            allocation = ast.literal_eval(allocation_str)
            y_axis_value = allocation[1]
            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            if os.path.isfile(fluxes_file):
                z_axis_value = load_json(fluxes_file)["Y"]
            else:
                z_axis_value = None
            data[index] = (x_axis_value, y_axis_value, z_axis_value, index)
    return data


def plot_data(folder):
    fig, ax = plt.subplots(2, 2, figsize = (5,6), gridspec_kw={'width_ratios': [1, 0.1]})
    data = get_data(folder)
    print(data)

    df = pd.DataFrame(data.values(), columns=['x', 'y', 'flux', 'index'])
    
    ############
    # Plot fluxes
    ############
    z_grid = df.pivot(index='y', columns='x', values='flux')
    x_points = z_grid.columns.values
    y_points = z_grid.index.values
    z_points = z_grid.values  # 2D array of shape (len(y), len(x))
    
    mesh0 = ax[0][0].pcolormesh(x_points, y_points, z_points, cmap='viridis', shading='auto',
                            #norm = Normalize(vmin=0, vmax=1)
                            #norm=LogNorm(vmin=np.nanmin(z_points), vmax=np.nanmax(z_points))
    )
    #ax[0][0].scatter(df['x'], df['y'], color='red', s=5, zorder=5)
    #for _, point in df.iterrows():
    #    ax[0].annotate(f"{point['index'].lstrip('0')}", (point['x'], point['y']))

    fig.colorbar(mesh0, cax=ax[0][1], label='flux')

    # Set log scale on whichever axes need it
    ax[0][0].set_xlabel("relative radius of inner membrane r*/R")
    ax[0][0].set_ylabel("proportion of enzyme \n in outermost region")
    ax[0][1].set_box_aspect(10)
    ax[0][0].set_box_aspect(1)

    mesh1 = ax[1][0].pcolormesh(x_points, y_points, z_points, cmap='viridis', shading='auto',
        norm = Normalize(vmin=3.5e-7, vmax=3.90e-7)
        #norm=LogNorm(vmin=np.nanmin(z_points), vmax=np.nanmax(z_points))
    )
    fig.colorbar(mesh1, cax=ax[1][1], label='flux')
    ax[1][1].set_box_aspect(10)
    ax[1][0].set_box_aspect(1)
    
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "fluxes.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)
    # python data/02A_enzymaticXtoY_1InnerBoundary_modifyingInnerRadius_modifyingEnzymeAllocation/analysis.py data/02A_enzymaticXtoY_1InnerBoundary_modifyingInnerRadius_modifyingEnzymeAllocation