import optuna
import json
from pathlib import Path
import argparse
import os
from auxiliary_functions_using_standard_library import load_json
from optuna.trial import TrialState

if __name__ == "__main__":
    # Parse arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--folder_to_solve", type=str, required=True)
    parser.add_argument("--round", type=int, required=True)
    parser.add_argument("--product_to_maximize", type=str, required=True)
    
    args = parser.parse_args()
    FOLDER_TO_SOLVE = args.folder_to_solve
    round_idx = args.round
    product_to_maximize = args.product_to_maximize

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
    for trial_idx, trial_dir in enumerate(trial_dirs):
        fluxes = load_json(os.path.join(trial_dir, "fluxes.json"))
        try:
            scalar = fluxes[product_to_maximize]
        except:
            raise ValueError(f"Could not find species {product_to_maximize}.")
        # Informs Optuna about the result of the trial
        trial = study.trials[trial_idx]
        if trial.state == TrialState.WAITING:
            study.tell(trial_idx, scalar)
        elif trial.state == TrialState.COMPLETE:
            print(f"Trial {trial_idx} already complete (value={trial.value}), skipping tell()")
        else:
            # RUNNING or PRUNED - safe to tell
            study.tell(trial_idx, scalar)
