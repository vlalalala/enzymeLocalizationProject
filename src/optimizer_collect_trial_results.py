import optuna
import json
from pathlib import Path
import argparse
import os
from auxiliary_functions_using_standard_library import load_json
from optuna.trial import TrialState
import re
from auxiliary_functions import read_yaml_file

if __name__ == "__main__":
    # Parse arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--folder_to_solve", type=str, required=True)
    parser.add_argument("--round", type=int, required=True)
    
    args = parser.parse_args()
    FOLDER_TO_SOLVE = args.folder_to_solve
    round_idx = args.round

    optimization_params = read_yaml_file(os.path.join(FOLDER_TO_SOLVE, "parameters_optimization.yaml"))
    product_to_maximize = optimization_params["species_to_maximize"]

    round_folder = os.path.join(
        FOLDER_TO_SOLVE, f"optimization_round_{round_idx}"
    )

    storage = f"sqlite:///{FOLDER_TO_SOLVE}/optuna_study.db"

    study = optuna.load_study(
        study_name="resource_allocation",
        storage=storage
    )
    # Find number of trials used
    n_trials = load_json(Path(FOLDER_TO_SOLVE).parent / "optimization_config.lock.json")["N_TRIALS"]
    trial_dirs = [os.path.join(round_folder, f"trial_{trial_idx}") for trial_idx in range(n_trials)]
          
    # Feed results back to Optuna
    for trial_dir in trial_dirs:
        # Extract trial_idx from name
        trial_idx = int(re.search(r"trial_(\d+)", os.path.basename(trial_dir)).group(1))
        matching_trials = [
            t for t in study.trials
            if t.user_attrs.get("round") == round_idx
            and t.user_attrs.get("trial_idx") == trial_idx
        ]
        if len(matching_trials)!=1:
            raise ValueError("The number of matching trials to round and trial_idx is not equal to 1.")
        trial = matching_trials[-1]
        # skip pruned trials
        if trial.state == TrialState.PRUNED:
            print(f"Trial {trial_idx} was pruned, skipping.")
            continue
        
        # Non-pruned trials:
        fluxes = load_json(os.path.join(trial_dir, "fluxes.json"))
        try:
            scalar = fluxes[product_to_maximize]
        except:
            raise ValueError(f"Could not find species {product_to_maximize}.")
        
        if trial.state == TrialState.WAITING or trial.state == TrialState.RUNNING:
            # TrialState WAITING if the trial was created with ask()
            # but tell() was never called
            # RUNNING. the trial is still marked as running (this can
            # happen if the previous script run crashed before calling tell()
            # ) ... I do not really understand the difference
            study.tell(trial.number, scalar)
        elif trial.state == TrialState.COMPLETE:
            # tell() was already called for this trial
            print(f"Trial {trial_idx} already complete (value={trial.value}), skipping tell()")
        else:
            raise ValueError(f"The trialState is neither waiting, running, complete nor pruned: {trial.state}")
