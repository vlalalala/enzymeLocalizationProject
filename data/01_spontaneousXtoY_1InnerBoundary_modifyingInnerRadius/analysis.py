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
            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            if os.path.isfile(fluxes_file):
                y_axis_value = load_json(fluxes_file)["Y"]
            else:
                y_axis_value = None
            data[index] = (x_axis_value, y_axis_value)

    return data


def plot_data(folder):
    fig, ax = plt.subplots(1, 1, figsize = (4,3))
    data = get_data(folder)

    df = pd.DataFrame(data.values(), columns=['x', 'flux'])

    ax.scatter(df["x"], df["flux"])

    ax.set_xlabel("relative radius of inner membrane r*/R")
    ax.set_ylabel("flux")

    
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "result.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)