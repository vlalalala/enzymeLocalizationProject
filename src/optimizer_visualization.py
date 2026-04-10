
import os
import argparse
import optuna
import re
from optuna.trial import TrialState


if __name__ == "__main__":
    # Parse arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--folder_to_solve", type=str, required=True)
    args = parser.parse_args()
    FOLDER_TO_SOLVE = args.folder_to_solve
    
    # First, check whether convergence has already taken place.
    # If so, need to remove further trials in order to not affect the data
    potential_optimization_convergence_file = os.path.join(
        FOLDER_TO_SOLVE, "optimization_convergence.txt")
    if os.path.isfile(potential_optimization_convergence_file):
        with open(potential_optimization_convergence_file) as f:
            info = f.read()
        final_round = int(re.search(r"round (\d+)", info).group(1))
    else:
        final_round = None

    storage = f"sqlite:///{FOLDER_TO_SOLVE}/optuna_study.db"
    study = optuna.load_study(
        study_name="resource_allocation",
        storage=storage
    )

    if final_round is not None:
        trials = [
            t for t in study.trials
            if t.user_attrs.get("round", float("inf")) <= final_round
        ]
    else:
        trials = study.trials

    # Create an alias of the study but only with the filtered trials
    shadow_study = optuna.create_study(direction=study.direction)
    shadow_study.add_trials(trials)
    
    fig1 = optuna.visualization.plot_intermediate_values(shadow_study)
    #fig1.write_html(os.path.join(FOLDER_TO_SOLVE, "intermediate_values.html"))
    fig1.show()
    fig2 = optuna.visualization.plot_slice(shadow_study)
    fig2.show()
    fig3 = optuna.visualization.plot_edf(shadow_study)
    fig3.show()

    # Most visualization functions accept a trials argument
    fig4 = optuna.visualization.plot_optimization_history(shadow_study)
    fig4.show()
    
    try:
        importances = optuna.importance.get_param_importances(shadow_study)
        for param, importance in importances.items():
            print(f"{param}: {importance:.4f}")
        fig5 = optuna.visualization.plot_param_importances(shadow_study)
        fig5.show()
    except RuntimeError as e:
        print(f"Could not compute importances: {e}")
        print(f"Number of completed trials: {len([t for t in shadow_study.trials if t.state == TrialState.COMPLETE])}")
    
    fig6 = optuna.visualization.plot_parallel_coordinate(shadow_study)
    fig6.show()
    fig7 = optuna.visualization.plot_contour(shadow_study)
    fig7.show()
    

    #importances = optuna.importance.get_param_importances(shadow_study)
    #print(importances)