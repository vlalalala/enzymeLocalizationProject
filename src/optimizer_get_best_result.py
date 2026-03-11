import optuna
import json
import argparse
from auxiliary_functions import dump_json

if __name__ == "__main__":
    # Parse arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--folder_to_solve", type=str, required=True)
    args = parser.parse_args()
    FOLDER_TO_SOLVE = args.folder_to_solve
    
    storage = f"sqlite:///{FOLDER_TO_SOLVE}/optuna_study.db"
    study = optuna.load_study(
        study_name="resource_allocation",
        storage=storage
    )

    best = study.best_trial
    result = {
        "best_value": best.value,
        "best_params": best.params
    }
    dump_json(FOLDER_TO_SOLVE, "best_results", result)

    print(f"Best value: {best.value}")
    print(f"Best params: {best.params}")