import sys
import matplotlib.pyplot as plt
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src'))
from auxiliary_functions_using_standard_library import load_json
import re
import pandas as pd
import numpy as np
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
            region_A = [region for region, val in allocationA.items() if val == 1][0]         
            region_B = [region for region, val in allocationB.items() if val == 1][0]  
            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            if os.path.isfile(fluxes_file):
                z_axis_value = load_json(fluxes_file)["Z"]
            else:
                z_axis_value = None
            species_df = pd.read_csv(os.path.join(combined_folder, "species.csv"))
            X_permeability = species_df.loc[
                (species_df["name"] == "X"),
                "permeability_constant"].item()
            Y_permeability = species_df.loc[
                (species_df["name"] == "Y"),
                "permeability_constant"].item()
            
            data[index] = (region_A, region_B, z_axis_value, index, X_permeability, Y_permeability)

    return data


def plot_data(folder):
    
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=['region_A', 'region_B', 'flux', 'index', 'X_permeability', 'Y_permeability'])
    X_permeabilities = np.sort(df['X_permeability'].unique())
    Y_permeabilities = np.sort(df['Y_permeability'].unique())

    fig, ax = plt.subplots(len(X_permeabilities),len(Y_permeabilities), figsize = (4*len(X_permeabilities),4*len(X_permeabilities)))
    for row, X_permeability in enumerate(X_permeabilities):
        for column, Y_permeability in enumerate(Y_permeabilities):
            current_df = df[(df["X_permeability"]==X_permeability) & (df["Y_permeability"]==Y_permeability)]
            
            z_grid = current_df.pivot(index='region_B', columns='region_A', values='flux')
            x_points = z_grid.columns.values
            y_points = z_grid.index.values
            z_points = z_grid.values  # 2D array of shape (len(y), len(x))
            mesh = ax[row][column].pcolormesh(x_points, y_points, z_points, cmap='viridis', shading='auto',
                                    #norm = Normalize(vmin=0, vmax=1)
                                    #norm=LogNorm(vmin=np.nanmin(z_points), vmax=np.nanmax(z_points)
            )
            ax[row][column].scatter(current_df['region_B'], current_df['region_A'], color='red', s=5, zorder=5)
            #fig.colorbar(mesh, cax=ax[row][1], label='flux')
            for _, point in current_df.iterrows():
                ax[row][column].annotate("{:.1e} \n".format(point['flux']) + str(point['index']).lstrip('0'), (point['region_B'], point['region_A']))

            # Set log scale on whichever axes need it
            ax[row][column].set_xlabel("region in which enzyme A is")
            ax[row][column].set_ylabel("region in which enzyme B is")
            #ax[row][1].set_box_aspect(10)
            ax[row][column].set_box_aspect(1)
            ax[row][column].set_title(f"p_X = {X_permeability}, p_Y = {Y_permeability}")
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "fluxes.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)
    # python data/04h_enzymaticXtoY_enzymaticYtoZ_4InnerBoundaries_modifyingEnzymeAllocations_modifyingPermeability/analysis.py data/04h_enzymaticXtoY_enzymaticYtoZ_4InnerBoundaries_modifyingEnzymeAllocations_modifyingPermeability