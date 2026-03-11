import optuna
import json
import numpy as np
import pandas as pd
import os
import csv
import shutil
import argparse
from pathlib import Path
from auxiliary_functions_using_standard_library import load_json
from auxiliary_functions import read_yaml_file, dump_in_yaml_file

def check_or_create_optimization_lock(df, n_trials, n_rounds):
    lock_path = os.path.join(df, "optimization_config.lock.json")
    os.makedirs(df, exist_ok=True)
    if os.path.exists(lock_path):
        with open(lock_path) as f:
            locked = json.load(f)
        if locked["N_TRIALS"] != n_trials or locked["N_ROUNDS"] != n_rounds:
            raise ValueError(
                f"N_TRIALS or N_ROUNDS changed since optimization started in {df}!\n"
                f"Locked:  N_TRIALS={locked['N_TRIALS']}, N_ROUNDS={locked['N_ROUNDS']}\n"
                f"Current: N_TRIALS={n_trials}, N_ROUNDS={n_rounds}\n"
                f"Delete {lock_path} only if you know what you are doing."
            )
    else:
        with open(lock_path, "w") as f:
            json.dump({"N_TRIALS": n_trials, "N_ROUNDS": n_rounds}, f, indent=2)
        print(f"Created optimization lock in {lock_path}")


def softmax(z):
    e = np.exp(z - np.max(z))
    return e / e.sum()

if __name__ == "__main__":
    # Parse arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--folder_to_solve", type=str, required=True)
    parser.add_argument("--round", type=int, required=True)
    parser.add_argument("--n_trials", type=int, required=True)
    parser.add_argument("--n_rounds", type=int, required=True)
    args = parser.parse_args()


    args = parser.parse_args()
    FOLDER_TO_SOLVE = args.folder_to_solve
    round_idx = args.round
    n_trials = args.n_trials
    n_rounds = args.n_rounds

    # Create lock of N_TRIALS and N_ROUNDS
    df = str(Path(args.folder_to_solve).parent)  # e.g. data_private/optuna_test
    check_or_create_optimization_lock(df, n_trials, n_rounds)

    # Find out number of enzymes
    enzymes_df = pd.read_csv(os.path.join(FOLDER_TO_SOLVE, "enzymes.csv"))
    n_enzymes = len(enzymes_df) # the first row is the header
    # Find out total enzyme quantity and relative distance 
    conditions_info = read_yaml_file(os.path.join(FOLDER_TO_SOLVE, "parameters_value_conditions.yaml"))
    total_enzyme_quantity = conditions_info["enzyme_total_fixed_quantity"]
    minimum_distance_between_membranes =  conditions_info["minimum_distance_between_membranes"]  
    # Find out number of regions and minimum distance between membranes and the origin
    geometry_info = read_yaml_file(os.path.join(FOLDER_TO_SOLVE, "parameters_geometry.yaml"))
    n_regions = len(geometry_info["geometry_config"]["internal_membrane_relative_radii"])+1
    external_radius = geometry_info["geometry_config"]["outer_membrane_radius"]
    n_inner_membranes = n_regions - 1
    minimum_normalized_distance_between_membranes = minimum_distance_between_membranes / external_radius
    
    
    

    storage = f"sqlite:///{FOLDER_TO_SOLVE}/optuna_study.db"
    # Create study on round 0, load on subsequent rounds
    study = optuna.create_study(
        study_name="resource_allocation",
        storage=storage,
        direction="maximize",
        load_if_exists=True
    )

    # Ask Optuna for a batch of N_TRIALS suggestions
    for trial_idx in range(n_trials):
        trial = study.ask()  # get a suggested trial without evaluating it
        trial.set_user_attr("round", round_idx)
        # --- Allocation to different enzymes: (n_enzymes - 1) free params ---
        #softmax([-5, 0, 0]) → [0.007, 0.497, 0.497]  # nearly equal split, one tiny
        #softmax([ 5, 0, 0]) → [0.987, 0.007, 0.007]  # almost all in first type
        #softmax([-5, 5, 0]) → [0.000, 0.993, 0.007]  # extreme concentration difference
        # ranges from -5 to 5 makes it possible to represent allocations
        # from nearly 0% to ~99% for any enzyme
        if n_enzymes > 0:
            z_enzymes = [trial.suggest_float(f"z_type_{i}", -5, 5)
                    for i in range(n_enzymes - 1)]
            z_enzymes.append(0.0) # last enzyme
            enzyme_allocations = (softmax(np.array(z_enzymes)) * total_enzyme_quantity).tolist()

        # --- Allocation of enzymes to different regions: n_enzymes x (n_regions - 1) free params ---
        regional_alloc = []
        for t in range(n_enzymes):
            z_regions = [trial.suggest_float(f"z_region_{t}_{r}", -5, 5)
                        for r in range(n_regions - 1)]
            z_regions.append(0.0)
            regional_alloc_list = softmax(np.array(z_regions)).tolist()
            regional_alloc.append({
                region: allocation
                for region, allocation in enumerate(regional_alloc_list)
            })

        # --- Allocation of membrane positions ---
        # Total "slack" available after reserving minimum gaps
        # Each region must be >= minimum_normalized_distance_between_membranes
        slack = 1.0 - n_regions * minimum_normalized_distance_between_membranes
        assert slack > 0, "distance between membranes too large for the number of thresholds"
        # Suggest unconstrained positive values (regions), normalize to slack
        raw_region_sizes = [trial.suggest_float(f"region_size_{i}", 0, 1)
            for i in range(n_regions)]
        # Normalize so region sizes sum to slack
        total = sum(raw_region_sizes)
        region_sizes = [g / total * slack for g in raw_region_sizes]
        inner_membrane_radii = []
        current_membrane_radius = 0.0
        for i in range(n_inner_membranes):
            current_membrane_radius += minimum_normalized_distance_between_membranes + region_sizes[i]
            inner_membrane_radii.append(current_membrane_radius)

        # Need constraint: distance between adjacent membranes, maximum concentration

        # Copy species.csv, spontaneous_reactions.csv, enzymatic_reactions.csv,
        # parameters_discretization.yaml, parameters_solver_input.yaml,
        # parameters_solver_output.yaml, parameters_value_conditions.yaml
        # onto each trial
        for file in [
            "species.csv", "spontaneous_reactions.csv", "enzymatic_reactions.csv",
            "parameters_discretization.yaml", "parameters_solver_input.yaml",
            "parameters_solver_output.yaml", "parameters_value_conditions.yaml"
        ]:
            shutil.copy(
                os.path.join(FOLDER_TO_SOLVE, f"{file}"),
                os.path.join(FOLDER_TO_SOLVE, f"optimization_round_{round_idx}/trial_{trial_idx}/{file}"),
            )
        
        # Create a modified parameters_geometry.yaml with the correct membrane radii
        geometry_info["geometry_config"]["internal_membrane_relative_radii"] = inner_membrane_radii
        dump_in_yaml_file(os.path.join(
                FOLDER_TO_SOLVE,
                f"optimization_round_{round_idx}/trial_{trial_idx}/parameters_geometry.yaml"),
                geometry_info
        )
        if n_enzymes > 0:
            # Create a modified enzymes.csv with the correct enzyme allocation
            enzymes_df["quantity"] = enzyme_allocations
            enzymes_df["allocation"] = regional_alloc
            enzymes_df.to_csv(
                os.path.join(
                    FOLDER_TO_SOLVE,
                    f"optimization_round_{round_idx}/trial_{trial_idx}/enzymes.csv"),
                index=False)
        else:
            shutil.copy(
                os.path.join(FOLDER_TO_SOLVE, "enzymes.csv"),
                os.path.join(FOLDER_TO_SOLVE, f"optimization_round_{round_idx}/trial_{trial_idx}/enzymes.csv"),
            )