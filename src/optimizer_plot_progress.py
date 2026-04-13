import optuna
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from optimizer_suggest_trials import params_to_physical
from auxiliary_functions import read_yaml_file
import pandas as pd
from auxiliary_functions_using_standard_library import load_json
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np
from pathlib import Path
import re

markers = ['o', 's', '^', 'D', 'v', 'P', '*', 'X']  # circle, square, triangle, diamond, etc.

def plot_optimization_progress(
    folder_to_solve,
    data,
    enzymes_df,
    n_regions,
    n_rounds_to_plot,
    n_trials
):
    n_enzymes = len(enzymes_df)
    enzyme_names = enzymes_df["name"].tolist()
    num_rows = 1 + n_enzymes + 2
    fig, ax = plt.subplots(num_rows, 2, figsize = (4, 3 * num_rows), gridspec_kw={"width_ratios": [20, 1]})
    #print(data)
    all_fluxes = [
        data[round_idx][trial_idx]["flux_to_maximize"]
        for round_idx in range(n_rounds_to_plot)
        for trial_idx in range(n_trials)
    ] # all_fluxes is of length N_TRIALS * N_ROUNDS
    #print(all_fluxes)
    non_pruned_fluxes = [flux for flux in all_fluxes if flux is not None]
    if len(non_pruned_fluxes)>1:
        norm = mcolors.Normalize(vmin=min(non_pruned_fluxes), vmax=max(non_pruned_fluxes))
    cmap = cm.viridis

    for round_idx in range(n_rounds_to_plot):
        for trial_idx in range(n_trials):
            trial_data = data[round_idx][trial_idx]
            flux = trial_data["flux_to_maximize"]
            if flux is None:
                style = {'edgecolors': 'black', 'facecolors': 'none', 'linewidths': 1.5}
            else:
                color = cmap(norm(flux)) # if flux is not None, then the variable norm exists
                style = {'facecolors': color}
            ###################
            # Row 0 : plot membrane positions
            ###################
            radii = trial_data["inner_membrane_radii"]
            #print(round_idx, trial_idx, radii, trial_data["pruned"])
            for membrane_idx, radius in enumerate(radii):
                marker = markers[membrane_idx % len(markers)]
                ax[0][0].scatter(
                    round_idx, radius,
                    #color=color,
                    marker=marker,
                    s=50,
                    alpha=0.4,
                    **style
                )
                
            ###################
            # Row 1 : plot enzyme allocation
            ###################
            enzyme_allocations_list = trial_data["enzyme_allocations"]
            if enzyme_allocations_list is not None: # in case any enzyme has been allocated
                for enzyme_allocation_idx, enzyme_allocation in enumerate(enzyme_allocations_list):
                    marker = markers[enzyme_allocation_idx % len(markers)]
                    ax[1][0].scatter(
                        round_idx, enzyme_allocation,
                        color=color,
                        marker=marker,
                        s=50,
                        alpha=0.7
                    )

            #########################
            # Other rows: 
            ##########################
            for enzyme_idx in range(n_enzymes):
                enzyme_regional_alloc = trial_data["regional_alloc"][enzyme_idx]
                ax[2+enzyme_idx][0].set_title(f"enzyme {enzyme_names[enzyme_idx]}")
                for region, percentage in enzyme_regional_alloc.items():
                    marker = markers[region % len(markers)]
                    ax[2+enzyme_idx][0].scatter(
                        round_idx, percentage,
                        color=color,
                        marker=marker,
                        s=50,
                        alpha=0.7
                    )
    
    ####################################
    # Create the legends
    ####################################

    # Legend for membrane index -> shape 
    for membrane_idx in range(len(data[0][0]["inner_membrane_radii"])):
        ax[0][0].scatter([], [], marker=markers[membrane_idx % len(markers)],
                color="k", label=f"membrane {membrane_idx}")
    ax[0][0].legend(frameon=True)
    
    # Legend for enzyme -> shape 
    for enzyme_idx in range(n_enzymes):
        ax[1][0].scatter([], [], marker=markers[enzyme_idx % len(markers)],
                color="k", label=f"enzyme {enzyme_names[enzyme_idx]}")
    ax[1][0].legend(frameon=True)
    
    # Legend for enzyme -> shape 
    for region_idx in range(n_regions):
        ax[2][0].scatter([], [], marker=markers[region_idx % len(markers)],
                color="k", label=f"region {region_idx}")
    ax[2][0].legend(frameon=True)

    ####################################################
    # Add axes labels
    ####################################################
    for row in range(num_rows):
        ax[row][0].set_xlabel("round")
        ax[row][0].set_xticks(range(n_rounds_to_plot))
    ax[0][0].set_ylabel("normalized membrane radius r/R")
    ax[0][0].set_ylim(-0.05, 1.05)
    ax[0][0].set_yticks([0,0.25,0.5,0.75,1])
    ax[1][0].set_ylabel("allocation to enzyme")
    for enzyme_idx in range(n_enzymes):
        ax[2+enzyme_idx][0].set_ylabel("allocation of enzyme to region")
    
    ####################################################
    # Add colorbar
    ####################################################
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax = ax[0][1]
    plt.colorbar(sm, cax=cax, label="flux to maximize")
    for row in range(1, num_rows):
        ax[row][1].axis('off')

    fig.tight_layout()
    fig.savefig(os.path.join(folder_to_solve, "optimization_progress.png"), dpi=600)

def load_existing_data(folder_to_solve, study, n_rounds_to_load, n_trials, n_enzymes):
    # Find out total enzyme quantity and relative distance 
    conditions_info = read_yaml_file(os.path.join(folder_to_solve, "parameters_value_conditions.yaml"))
    total_enzyme_quantity = conditions_info["enzyme_total_fixed_quantity"]
    enzyme_maximum_concentration =  conditions_info["enzyme_maximum_concentration"]
    enzyme_maximum_concentration = 1 ####################################################################################  
    # Find out number of regions and minimum distance between membranes and the origin
    geometry_info = read_yaml_file(os.path.join(folder_to_solve, "parameters_geometry.yaml"))
    n_regions = len(geometry_info["geometry_config"]["internal_membrane_relative_radii"])+1
    external_radius = geometry_info["geometry_config"]["outer_membrane_radius"]
    data = {
        round_idx: {
            trial_idx: {
                "flux_to_maximize": None,
                "inner_membrane_radii": None,
                "regional_alloc": None,
                "enzyme_allocations": None
            }
            for trial_idx in range(n_trials)  # inner loop stays inside
        }
        for round_idx in range(n_rounds_to_load)      # outer loop for the outer dict
    }
    for round_idx in range(n_rounds_to_load):
        for trial_idx in range(n_trials):
            trial = next(
                (t for t in study.trials
                if t.user_attrs.get("round") == round_idx
                and t.user_attrs.get("trial_idx") == trial_idx),
                None
            )
            if trial is None:
                raise ValueError(f"There is something weird going on with round {round_idx}, trial {trial_idx}.")
            enzyme_allocations, regional_alloc, inner_membrane_radii, _ = params_to_physical(
                trial.params,
                n_enzymes=n_enzymes,
                n_regions=n_regions,
                total_enzyme_quantity=total_enzyme_quantity,
                enzyme_maximum_concentration=enzyme_maximum_concentration,
                external_radius=external_radius
            )  
            data[round_idx][trial_idx]["flux_to_maximize"] = trial.value
            data[round_idx][trial_idx]["inner_membrane_radii"] = inner_membrane_radii
            data[round_idx][trial_idx]["regional_alloc"] = regional_alloc
            data[round_idx][trial_idx]["enzyme_allocations"] = enzyme_allocations
    return data

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    largest_round_to_plot = int(sys.argv[2])
    # this is the number/name of the last one to plot!, not the index!

    storage_path = os.path.join(FOLDER_TO_SOLVE, "optuna_study.db")
    if not os.path.isfile(storage_path):
        raise FileNotFoundError(f"The file {storage_path} does not exist.")
    
    storage = f"sqlite:///{FOLDER_TO_SOLVE}/optuna_study.db"
    summaries = optuna.get_all_study_summaries(storage=storage)
    if len(summaries) == 0:
        raise ValueError(f"Within {storage} there are no studies...")

    study_name = "resource_allocation"
    study = optuna.load_study(
        study_name=study_name,
        storage=storage
    )

    optimization_config = load_json(os.path.join(Path(FOLDER_TO_SOLVE).parent, "optimization_config.lock.json"))
    n_trials = optimization_config["N_TRIALS"]
    n_rounds = optimization_config["N_ROUNDS"]
    # Find out number of enzymes
    enzymes_df = pd.read_csv(os.path.join(FOLDER_TO_SOLVE, "enzymes.csv"))
    n_enzymes = len(enzymes_df) # the first row is the header
    
    potential_optimization_convergence_file = os.path.join(FOLDER_TO_SOLVE, "optimization_convergence.txt")
    if os.path.isfile(potential_optimization_convergence_file):
        with open(potential_optimization_convergence_file) as f:
            info = f.read()
        final_round = int(re.search(r"round (\d+)", info).group(1))
        largest_round_to_plot -= 1 # largest round to plot is the number associated, so we have to decrease 1


    # if largest_round_to_plot "is the string 5", then n_rounds is actually 6 (since round 0 also counts)
    data = load_existing_data(
        FOLDER_TO_SOLVE, study,
        n_rounds_to_load = largest_round_to_plot+1,
        n_trials=n_trials,
        n_enzymes=n_enzymes
    )
    #print("data_loaded", data)
    geometry_info = read_yaml_file(os.path.join(FOLDER_TO_SOLVE, "parameters_geometry.yaml"))
    n_regions = len(geometry_info["geometry_config"]["internal_membrane_relative_radii"])+1
    
    

    plot_optimization_progress(
        FOLDER_TO_SOLVE,
        data,
        enzymes_df,
        n_regions,
        n_rounds_to_plot = largest_round_to_plot+1,
        n_trials=n_trials
    )