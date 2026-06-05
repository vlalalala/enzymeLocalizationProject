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
            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            if os.path.isfile(fluxes_file):
                flux = load_json(fluxes_file)["Z"]
            else:
                flux = None
            enzymatic_reactions_df = pd.read_csv(os.path.join(combined_folder, "enzymatic_reactions.csv"))
            A_catalytic_rate = enzymatic_reactions_df.loc[
                (enzymatic_reactions_df["enzyme"] == "A"),
                "k_cat"].item()
            B_catalytic_rate = enzymatic_reactions_df.loc[
                (enzymatic_reactions_df["enzyme"] == "B"),
                "k_cat"].item()
            C_catalytic_rate = enzymatic_reactions_df.loc[
                (enzymatic_reactions_df["enzyme"] == "C"),
                "k_cat"].item()
            data[index] = (allocationA, allocationB, allocationC,
                           A_catalytic_rate, B_catalytic_rate, C_catalytic_rate,
                           flux
            )

    return data


def plot_data(folder):
    
    data = get_data(folder)
    df = pd.DataFrame(data.values(), columns=[
        'allocationA', 'allocationB', 'allocationC',
        'A_catalytic_rate', 'B_catalytic_rate', 'C_catalytic_rate',        
        'flux']
    )
    A_catalytic_rates = np.sort(df["A_catalytic_rate"].unique())
    B_catalytic_rates = np.sort(df["B_catalytic_rate"].unique())

    orders = {0:"ABC", 1:"BCA", 2:"CAB", 3:"CBA", 4:"BAC", 5:"ACB"}
    ls = {0:":", 1:"dashed", 2:"dashdot", 3:"dashdot", 4:"dashed", 5:"dotted"}
    fig, ax = plt.subplots(len(A_catalytic_rates), len(B_catalytic_rates),
                           figsize = (4*len(B_catalytic_rates), 3*len(A_catalytic_rates)))
    for row, A_catalytic_rate in enumerate(A_catalytic_rates):
        for column, B_catalytic_rate in enumerate(B_catalytic_rates):
            for order_idx, order in orders.items():
                A_loc = order.index('A')
                B_loc = order.index('B')
                C_loc = order.index('C')
                A_dict = {key: 0 for key in range(4)}
                B_dict = {key: 0 for key in range(4)}
                C_dict = {key: 0 for key in range(4)}
                A_dict[A_loc + 1] = 1
                B_dict[B_loc + 1] = 1
                C_dict[C_loc + 1] = 1
                current_df = df[(df["allocationA"]==A_dict)
                            & (df["allocationB"]==B_dict)
                            & (df["allocationC"]==C_dict)
                            & (df["A_catalytic_rate"]==A_catalytic_rate)
                            & (df["B_catalytic_rate"]==B_catalytic_rate)
                            ]
                current_df = current_df.sort_values("C_catalytic_rate")
                ax[row][column].plot(current_df["C_catalytic_rate"], current_df["flux"], label = order,
                                     ls = ls[order_idx], alpha = 0.6)
                ax[row][column].set_xlabel("C_catalytic_rate")
                ax[row][column].set_ylabel("flux")
                ax[row][column].set_title(f"A_catalytic_rate = {A_catalytic_rate},\n B_catalytic_rate = {B_catalytic_rate}")
                if row == 0 and column == 0:
                    ax[row][column].legend(title = "order")
                ax[row][column].set_xscale("log")
                ax[row][column].set_yscale("log")
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "fluxes.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    plot_data(FOLDER_TO_SOLVE)
    # python data/06b_enzymaticWtoX_enzymaticXtoY_enzymaticYtoZ_3InnerBoundaries_modifyingOrder/analysis.py data/06b_enzymaticWtoX_enzymaticXtoY_enzymaticYtoZ_3InnerBoundaries_modifyingOrder


# ABC always better than ACB
# (this makes sense, since a lot of Y will escape without being converted to Z, which is done further back in the center)

# ABC and BAC tend to be the best 2
# (since it's good that the last step (C) is most towards the outside)
# EXCEPT if the catalytic rate of A is very large (this effect is not seen if the catalytic rate of B is very large)

# ABC is better than BAC if the catalytic rate of B is much larger than the catalytic rate of A
# BAC is better than ABC if the catalytic rate of A is much larger than the catalytic rate of B

# CBA is always better than CAB