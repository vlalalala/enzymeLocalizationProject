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
            kMA = enzymatic_reactions_df.loc[
                (enzymatic_reactions_df["enzyme"] == "A"),
                "k_M"].item()
            kMB = enzymatic_reactions_df.loc[
                (enzymatic_reactions_df["enzyme"] == "B"),
                "k_M"].item()
            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            if os.path.isfile(fluxes_file):
                Z_flux = load_json(fluxes_file)["Z"]
            else:
                Z_flux = None
            data[index] = (allocationA_1, allocationB_1, kMA, kMB, X_external_concentration, Z_flux)

    return data

def plot_data(folder):
    from matplotlib.lines import Line2D
    fig, ax = plt.subplots(1, 1, figsize = (7,5))
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=[
        'allocationA_1', 'allocationB_1', 'kMA', 'kMB', 'X_external_concentration', 'flux'])
    X_external_concentrations = sorted(df['X_external_concentration'].unique())
    kMAs = sorted(df['kMA'].unique())
    kMBs = sorted(df['kMB'].unique())
    fig, axs = plt.subplots(
        1, len(kMBs),
        figsize=(5*len(kMBs), 4),
        sharex=True,
        sharey=True
    )

    colors = plt.cm.tab10(np.linspace(0, 1, len(kMAs)))
    color_map = dict(zip(kMAs, colors))

    for ax, kMB in zip(axs, kMBs):

        for kMA in kMAs:

            current_df = df[
                (df["kMA"] == kMA) &
                (df["kMB"] == kMB)
            ]

            current_df = current_df[current_df["allocationA_1"] != 0]
            current_df = current_df.sort_values("allocationA_1")

            y = current_df["flux"] / np.amax(current_df["flux"])

            idx = current_df["flux"].idxmax()
            x_of_max_flux = current_df.loc[idx, "allocationA_1"]

            ax.plot(
                current_df["allocationA_1"],
                y,
                color=color_map[kMA],
                label=f"{Decimal(kMA):.2E}"
            )

            ax.scatter(
                [x_of_max_flux],
                [1],
                color=color_map[kMA],
                zorder=5
            )

        ax.set_title(f"kM_ B = {Decimal(kMB):.2E}")
        ax.set_xlabel("allocationA_1")

    axs[0].set_ylabel("normalized flux")

    # one legend for all panels
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        title="kM A",
        loc="upper center",
        ncol=len(kMAs)
    )

    plt.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(os.path.join(folder, "result.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)

# python data/04q_enzymaticXtoY_enzymaticYtoZ_2InnerBoundaries_modifyingRates/analysis.py data/04q_enzymaticXtoY_enzymaticYtoZ_2InnerBoundaries_modifyingRates