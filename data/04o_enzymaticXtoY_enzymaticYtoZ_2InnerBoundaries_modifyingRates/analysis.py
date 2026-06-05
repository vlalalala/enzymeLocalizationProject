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
from decimal import Decimal

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
            enzymatic_reactions_df = pd.read_csv(os.path.join(combined_folder, "enzymatic_reactions.csv"))
            kcatA = enzymatic_reactions_df.loc[
                (enzymatic_reactions_df["enzyme"] == "A"),
                "k_cat"].item()
            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            if os.path.isfile(fluxes_file):
                Z_flux = load_json(fluxes_file)["Z"]
            else:
                Z_flux = None
            data[index] = (allocationA_1, allocationB_1, kcatA, X_external_concentration, Z_flux)

    return data


def plot_data(folder):
    fig, ax = plt.subplots(1, 1, figsize = (7,5))
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=[
        'allocationA_1', 'allocationB_1', 'kcatA', 'X_external_concentration', 'flux'])

    X_external_concentrations = sorted(df['X_external_concentration'].unique())
    kcatAs = sorted(df['kcatA'].unique())
    for kcatA in kcatAs:
        current_df = df[df["kcatA"]==kcatA]
        current_df = current_df[current_df["allocationA_1"]!=0]
        current_df = current_df.sort_values("allocationA_1")
        idx = current_df["flux"].idxmax()  # skipna=True by default
        x_of_max_flux = current_df.loc[idx, "allocationA_1"]
        #print(x_of_max_flux)
        #print(current_df[current_df["allocationA_1"]==x_of_max_flux]["flux"])
        ax.scatter([x_of_max_flux], [1])
        #print(x_of_max_flux)
        #ax.scatter(current_df["allocationA_1"], current_df["flux"]/np.amax(current_df["flux"]), alpha = 0.2)
        ax.plot(current_df["allocationA_1"], current_df["flux"]/np.amax(current_df["flux"]),
                   label = "{:.2E}".format(Decimal(kcatA)), alpha = 0.5)

    ax.set_xlabel("Proportion of enzyme A (X->Y) within penultimate section ")
    ax.set_ylabel("flux of Z / max flux of Z")
    ax.legend(title="catalytic rate of A",)#bbox_to_anchor=(1.5, 0.5))

    fig.tight_layout()
    fig.savefig(os.path.join(folder, "result.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)

# python data/04o_enzymaticXtoY_enzymaticYtoZ_2InnerBoundaries_modifyingRates/analysis.py data/04o_enzymaticXtoY_enzymaticYtoZ_2InnerBoundaries_modifyingRates