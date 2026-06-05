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
            kcatB = enzymatic_reactions_df.loc[
                (enzymatic_reactions_df["enzyme"] == "B"),
                "k_cat"].item()
            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            if os.path.isfile(fluxes_file):
                Z_flux = load_json(fluxes_file)["Z"]
            else:
                Z_flux = None
            data[index] = (allocationA_1, allocationB_1, kcatA, kcatB, X_external_concentration, Z_flux)

    return data

def plot_data(folder):
    from matplotlib.lines import Line2D
    fig, ax = plt.subplots(1, 1, figsize = (7,5))
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=[
        'allocationA_1', 'allocationB_1', 'kcatA', 'kcatB', 'X_external_concentration', 'flux'])
    X_external_concentrations = sorted(df['X_external_concentration'].unique())
    kcatAs = sorted(df['kcatA'].unique())
    kcatBs = sorted(df['kcatB'].unique())
    fig, axs = plt.subplots(
        1, len(kcatBs),
        figsize=(5*len(kcatBs), 4),
        sharex=True,
        sharey=True
    )

    colors = plt.cm.tab10(np.linspace(0, 1, len(kcatAs)))
    color_map = dict(zip(kcatAs, colors))

    for ax, kcatB in zip(axs, kcatBs):

        for kcatA in kcatAs:

            current_df = df[
                (df["kcatA"] == kcatA) &
                (df["kcatB"] == kcatB)
            ]

            current_df = current_df[current_df["allocationA_1"] != 0]
            current_df = current_df.sort_values("allocationA_1")

            y = current_df["flux"] / np.amax(current_df["flux"])

            idx = current_df["flux"].idxmax()
            x_of_max_flux = current_df.loc[idx, "allocationA_1"]

            ax.plot(
                current_df["allocationA_1"],
                y,
                color=color_map[kcatA],
                label=f"{Decimal(kcatA):.2E}"
            )

            ax.scatter(
                [x_of_max_flux],
                [1],
                color=color_map[kcatA],
                zorder=5
            )

        ax.set_title(f"kcatB = {Decimal(kcatB):.2E}")
        ax.set_xlabel("allocationA_1")

    axs[0].set_ylabel("normalized flux")

    # one legend for all panels
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        title="kcatA",
        loc="upper center",
        ncol=len(kcatAs)
    )

    plt.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(os.path.join(folder, "result.png"), dpi = 300)

def plot_data_old(folder):
    from matplotlib.lines import Line2D
    fig, ax = plt.subplots(1, 1, figsize = (7,5))
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=[
        'allocationA_1', 'allocationB_1', 'kcatA', 'kcatB', 'X_external_concentration', 'flux'])

    X_external_concentrations = sorted(df['X_external_concentration'].unique())
    kcatAs = sorted(df['kcatA'].unique())
    kcatBs = sorted(df['kcatB'].unique())
    # Define styles
    colors = plt.cm.tab10(np.linspace(0, 1, len(kcatAs)))
    color_map = dict(zip(kcatAs, colors))

    linestyles = ["-", "--", ":"]
    ls_map = dict(zip(kcatBs, linestyles))
    markers = {
        kcatBs[0]: "o",
        kcatBs[1]: "s",
        kcatBs[2]: "^",
    }

    for kcatA in kcatAs:
        for kcatB in kcatBs:
            current_df = df[(df["kcatA"] == kcatA) &
                            (df["kcatB"] == kcatB)]

            current_df = current_df[current_df["allocationA_1"] != 0]
            current_df = current_df.sort_values("allocationA_1")

            idx = current_df["flux"].idxmax()
            x_of_max_flux = current_df.loc[idx, "allocationA_1"]

            ax.scatter([x_of_max_flux], [1],
                    color=color_map[kcatA], marker=markers[kcatB],)

            ax.plot(
                current_df["allocationA_1"],
                current_df["flux"] / np.amax(current_df["flux"]),
                color=color_map[kcatA],
                linestyle=ls_map[kcatB],
                alpha=0.5,
            )
    ax.set_xlabel("Proportion of enzyme A (X->Y) within penultimate section ")
    ax.set_ylabel("flux of Z / max flux of Z")
    ax.legend(title="catalytic rate of A",)#bbox_to_anchor=(1.5, 0.5))

    color_handles = [
    Line2D([0], [0],
            color=color_map[kcatA],
            lw=2,
            label=f"{Decimal(kcatA):.2E}")
        for kcatA in kcatAs
    ]

    leg1 = ax.legend(
        handles=color_handles,
        title="kcatA",
        loc="lower left"
    )

    # Legend for kcatB (linestyles)
    style_handles = [
        Line2D([0], [0],
            color="black",
            lw=2,
            linestyle=ls_map[kcatB],
            label=f"{Decimal(kcatB):.2E}")
        for kcatB in kcatBs
    ]

    leg2 = ax.legend(
        handles=style_handles,
        title="kcatB",
        loc="lower right"
    )

    ax.add_artist(leg1)


    fig.tight_layout()
    fig.savefig(os.path.join(folder, "result.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)

# python data/04p_enzymaticXtoY_enzymaticYtoZ_2InnerBoundaries_modifyingRates/analysis.py data/04p_enzymaticXtoY_enzymaticYtoZ_2InnerBoundaries_modifyingRates