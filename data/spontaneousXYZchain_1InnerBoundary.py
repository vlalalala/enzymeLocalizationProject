import sys
import matplotlib.pyplot as plt
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from auxiliary_functions_using_standard_library import load_json
import re
import pandas as pd
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize


def get_data(folder):
    """
    Load the best results from the combinations folders.
    """
    data = {}

    for folder_name in os.listdir(folder):
        match = re.match(r'^combined_(\d{6})$', folder_name)
        combined_folder = os.path.join(folder, folder_name)
        if match:
            index = match.group(1)  # keeps it as a string, preserving leading zeros
            # Get modified value
            spontaneous_reactions_df = pd.read_csv(
                os.path.join(combined_folder, "spontaneous_reactions.csv")
            )
            x_axis_value = spontaneous_reactions_df.loc[
                (spontaneous_reactions_df["start_species"] == "X")
                 & (spontaneous_reactions_df["end_species"] == "Y"),
                "k"].item()
            y_axis_value = spontaneous_reactions_df.loc[
                (spontaneous_reactions_df["start_species"] == "Y")
                 & (spontaneous_reactions_df["end_species"] == "Z"),
                "k"].item()
            json_path = os.path.join(combined_folder, 'best_result.json')
            if os.path.isfile(json_path):
                best_results = load_json(json_path)
                # Get optimized value
                z_axis_value1 = best_results["best_inner_membrane_radii"][0]
                z_axis_value2 = best_results["best_value"]
            else:
                z_axis_value1 = None
                z_axis_value2 = None
            data[index] = (x_axis_value, y_axis_value, z_axis_value1, z_axis_value2)

    return data


def plot_data(folder):
    fig, ax = plt.subplots(2, 2, figsize = (4,6), gridspec_kw={'width_ratios': [1, 0.1]})
    data = get_data(folder)
    
    df = pd.DataFrame(data.values(), columns=['x', 'y', 'membrane_pos', 'flux'])
    
    ############
    # Plot membrane positions
    ############
    z_grid_membrane_pos = df.pivot(index='y', columns='x', values='membrane_pos')
    x_points = z_grid_membrane_pos.columns.values
    y_points = z_grid_membrane_pos.index.values
    z_points_membrane_pos = z_grid_membrane_pos.values  # 2D array of shape (len(y), len(x))
    
    mesh = ax[0][0].pcolormesh(x_points, y_points, z_points_membrane_pos, cmap='viridis', shading='auto',
                            norm = Normalize(vmin=0, vmax=1)
                            #norm=LogNorm(vmin=np.nanmin(z_points), vmax=np.nanmax(z_points))
    )
    ax[0][0].scatter(df['x'], df['y'], color='red', s=5, zorder=5)
    fig.colorbar(mesh, cax=ax[0][1], label='optimal r*/R')

    # Set log scale on whichever axes need it
    ax[0][0].set_xscale('log')
    ax[0][0].set_yscale('log')
    ax[0][0].set_xlabel("reaction rate k for X->Y reaction")
    ax[0][0].set_ylabel("reaction rate k for Y->Z reaction")
    ax[0][1].set_box_aspect(10)
    ax[0][0].set_box_aspect(1)
    
    ############
    # Plot fluxes
    ############
    z_grid_flux = df.pivot(index='y', columns='x', values='flux')
    x_points = z_grid_flux.columns.values
    y_points = z_grid_flux.index.values
    z_points_flux = z_grid_flux.values  # 2D array of shape (len(y), len(x))
    
    mesh = ax[1][0].pcolormesh(x_points, y_points, z_points_flux, cmap='viridis', shading='auto',
                            #norm = Normalize(vmin=0, vmax=1)
                            norm=LogNorm(vmin=np.nanmin(z_points_flux), vmax=np.nanmax(z_points_flux))
    )
    ax[1][0].scatter(df['x'], df['y'], color='red', s=5, zorder=5)
    fig.colorbar(mesh, cax=ax[1][1], label='maximum flux')

    # Set log scale on whichever axes need it
    ax[1][0].set_xscale('log')
    ax[1][0].set_yscale('log')
    ax[1][0].set_xlabel("reaction rate k for X->Y reaction")
    ax[1][0].set_ylabel("reaction rate k for Y->Z reaction")
    ax[1][1].set_box_aspect(10)
    ax[1][0].set_box_aspect(1)


    fig.tight_layout()
    fig.savefig(os.path.join(folder, "optimization_result.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)