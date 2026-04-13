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
                z_axis_value = best_results["best_inner_membrane_radii"][0]
            else:
                z_axis_value = None
            data[index] = (x_axis_value, y_axis_value, z_axis_value)

    return data


def plot_data(folder):
    fig, ax = plt.subplots(1, 2, figsize = (4,3), gridspec_kw={'width_ratios': [1, 0.1]})
    data = get_data(folder)
    
    df = pd.DataFrame(data.values(), columns=['x', 'y', 'z'])
    # Pivot into a 2D grid
    z_grid = df.pivot(index='y', columns='x', values='z')
    x_points = z_grid.columns.values
    y_points = z_grid.index.values
    z_points = z_grid.values  # 2D array of shape (len(y), len(x))
    mesh = ax[0].pcolormesh(x_points, y_points, z_points, cmap='viridis', shading='auto',
                            norm = Normalize(vmin=0, vmax=1)
                            #norm=LogNorm(vmin=np.nanmin(z_points), vmax=np.nanmax(z_points))
    )
    ax[0].scatter(df['x'], df['y'], color='red', s=10, zorder=5)
    fig.colorbar(mesh, cax=ax[1], label='optimal r*/R')

    # Set log scale on whichever axes need it
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[0].set_xlabel("reaction rate k for X->Y reaction")
    ax[0].set_ylabel("reaction rate k for Y->Z reaction")
    ax[1].set_box_aspect(10)
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "optimization_result.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)