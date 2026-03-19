import optuna
import os
import argparse
from auxiliary_functions import dump_json
from optimizer_suggest_trials import params_to_physical
import pandas as pd
from auxiliary_functions import read_yaml_file
from optuna.trial import TrialState
from pathlib import Path
from auxiliary_functions_using_standard_library import load_json


if __name__ == "__main__":
    # Parse arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--folder_to_solve", type=str, required=True)
    parser.add_argument("--round", type=int, required=True)
    args = parser.parse_args()
    FOLDER_TO_SOLVE = args.folder_to_solve
    round_idx =  args.round

    optimization_config = load_json(os.path.join(Path(FOLDER_TO_SOLVE).parent, "optimization_config.lock.json"))
    n_trials = optimization_config["N_TRIALS"]
    n_rounds = optimization_config["N_ROUNDS"]
    
    storage = f"sqlite:///{FOLDER_TO_SOLVE}/optuna_study.db"
    study = optuna.load_study(
        study_name="resource_allocation",
        storage=storage
    )

    enzymes_df = pd.read_csv(os.path.join(FOLDER_TO_SOLVE, "enzymes.csv"))
    n_enzymes = len(enzymes_df) # the first row is the header
    # Find out total enzyme quantity and relative distance 
    conditions_info = read_yaml_file(os.path.join(FOLDER_TO_SOLVE, "parameters_value_conditions.yaml"))
    total_enzyme_quantity = conditions_info["enzyme_total_fixed_quantity"]
    enzyme_maximum_concentration =  conditions_info["enzyme_maximum_concentration"]  
    # Find out number of regions and size
    geometry_info = read_yaml_file(os.path.join(FOLDER_TO_SOLVE, "parameters_geometry.yaml"))
    n_regions = len(geometry_info["geometry_config"]["internal_membrane_relative_radii"])+1
    external_radius = geometry_info["geometry_config"]["outer_membrane_radius"]

    eligible_trials = [
        t for t in study.trials
        if t.state == TrialState.COMPLETE
        and t.user_attrs.get("round", 0) <= round_idx
    ]

    if not eligible_trials:
        raise ValueError(f"No completed trials found up to round {round_idx}")
    
    print("Eligible trials: ", eligible_trials)

    best = max(eligible_trials, key=lambda t: t.value)
    best_enzyme_allocations, best_regional_alloc, best_inner_membrane_radii, _ = params_to_physical(
        best.params,
        n_enzymes=n_enzymes,
        n_regions=n_regions,
        total_enzyme_quantity=total_enzyme_quantity,
        enzyme_maximum_concentration=enzyme_maximum_concentration,
        external_radius=external_radius
    )

    result = {
        "best_value": best.value,
        "best_params": best.params,
        "best_enzyme_allocations": best_enzyme_allocations,
        "best_regional_alloc": best_regional_alloc,
        "best_inner_membrane_radii": best_inner_membrane_radii

    }
    dump_json(FOLDER_TO_SOLVE, f"optimization_round_{round_idx}_best", result)
    
    if round_idx == n_rounds-1:
        dump_json(FOLDER_TO_SOLVE, "best_result", result)