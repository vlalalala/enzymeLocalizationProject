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
            allocationC_str = enzymes_df.loc[(enzymes_df["name"] == "C"), "allocation"].item()
            allocationC = ast.literal_eval(allocationC_str)
            proportionAinner = allocationA[0]
            proportionBinner = allocationB[0]
            proportionCinner = allocationC[0]          
            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            if os.path.isfile(fluxes_file):
                flux = load_json(fluxes_file)["Z"]
            else:
                flux = None
            data[index] = (proportionAinner, proportionBinner, proportionCinner, flux)

    return data


def plot_data(folder):
    
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=['proportionAinner', 'proportionBinner', 'proportionCinner', 'flux'])

    fig, ax = plt.subplots(2, 2, figsize = (6, 6))
    for row in range(2):
        for column in range(2):
            current_df = df[(df["proportionAinner"]==row) & (df["proportionCinner"]==column)]
            ax[row][column].scatter(current_df["proportionBinner"], current_df["flux"])
            ax[row][column].set_xlabel("proportion of \n enzyme B quantity \n in innermost region")
            ax[row][column].set_ylabel("flux")
            ax[row][column].set_title(f"proportionAinner = {row},\n proportionCinner = {column}")
    
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "fluxes.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)
    # python data/06a_enzymaticWtoX_enzymaticXtoY_enzymaticYtoZ_1InnerBoundary_modifyingAllocation2ndEnzyme/analysis.py data/06a_enzymaticWtoX_enzymaticXtoY_enzymaticYtoZ_1InnerBoundary_modifyingAllocation2ndEnzyme