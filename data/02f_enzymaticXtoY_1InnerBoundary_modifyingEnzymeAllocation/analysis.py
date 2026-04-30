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
            allocation_str = enzymes_df.loc[(enzymes_df["name"] == "A"), "allocation"].item()
            allocation = ast.literal_eval(allocation_str)
            allocation_in_outermost_region = allocation[1]
            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            if os.path.isfile(fluxes_file):
                flux = load_json(fluxes_file)["Y"]
            else:
                flux = None
            data[index] = (allocation_in_outermost_region, flux, index)

    return data


def plot_data(folder):
    fig, ax = plt.subplots(1, 1, figsize = (4,3))
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=['allocation_in_outermost_region', 'flux', 'index'])
    
    ax.scatter(df["allocation_in_outermost_region"], df["flux"])
    ax.set_xlabel("proportion of enzyme in outermost region")
    ax.set_ylabel("flux")

    fig.tight_layout()
    fig.savefig(os.path.join(folder, "fluxes.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)
    # python data/02f_enzymaticXtoY_1InnerBoundary_modifyingEnzymeAllocation/analysis.py data/02f_enzymaticXtoY_1InnerBoundary_modifyingEnzymeAllocation