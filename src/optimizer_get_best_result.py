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
from optimizer_plot_progress import load_existing_data, plot_optimization_progress
from auxiliary_functions_framework_organization_using_standard_library import DelayedSignals
import numpy as np

def get_best_result_up_until_round_specified(
        folder_to_solve, round_idx, study
):
    """
    Creates a file named optimization_round_{round_idx}_best which tracks
    the best value up until that round (inclusive).
    
    :param folder_to_solve: Description
    """
    enzymes_df = pd.read_csv(os.path.join(folder_to_solve, "enzymes.csv"))
    n_enzymes = len(enzymes_df) # the first row is the header
    # Find out total enzyme quantity and relative distance 
    conditions_info = read_yaml_file(os.path.join(folder_to_solve, "parameters_value_conditions.yaml"))
    total_enzyme_quantity = conditions_info["enzyme_total_fixed_quantity"]
    enzyme_maximum_concentration =  conditions_info["enzyme_maximum_concentration"]  
    # Find out number of regions and size
    geometry_info = read_yaml_file(os.path.join(folder_to_solve, "parameters_geometry.yaml"))
    n_regions = len(geometry_info["geometry_config"]["internal_membrane_relative_radii"])+1
    external_radius = geometry_info["geometry_config"]["outer_membrane_radius"]

    eligible_trials = [
        t for t in study.trials
        if t.state == TrialState.COMPLETE
        and t.user_attrs.get("round", 0) <= round_idx
    ]

    if not eligible_trials:
        raise ValueError(f"No completed trials found up to round {round_idx}")
    
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
        "trial_round_best": best.user_attrs.get("round"),
        "trial_idx_best": best.user_attrs.get("trial_idx"),
        "best_params": best.params,
        "best_enzyme_allocations": best_enzyme_allocations,
        "best_regional_alloc": best_regional_alloc,
        "best_inner_membrane_radii": best_inner_membrane_radii,
        "non-pruned trials run": len(eligible_trials)
    }
    dump_json(folder_to_solve, f"optimization_round_{round_idx}_best", result)

    return result

#def track_relative_changes(dictionary):
#    """
#    Gives the relative difference of the round relative to the previous one.
    #"""#
#    relative_changes = {}
#    sorted_rounds = sorted(dictionary.keys())
#    for r in sorted_rounds:
#        if r == 0:
#            continue
#        current_value = dictionary[r]
#        previous_value = dictionary[r-1]
#        if previous_value == 0 or current_value is None:
#            relative_change = None
#        else:
#            relative_change = (current_value - previous_value) / abs(previous_value)#

#       relative_changes[r] = {
#            "relative_improvement": relative_change,
#        }
#    return relative_changes


def get_convergence(
        folder_to_solve,
        round_idx,
        number_of_trials_to_run_before_stopping,
        value_negligible_relative_change,
        n_startup_trials
    ):
    """
    Only completed trials count toward n_startup_trials — pruned and failed trials do not count.
    """
    # First, load all of the data from the best trials up until (and including) round_idx
    best_data = {}
    for round in range(round_idx+1):
        best_round_data = load_json(os.path.join(folder_to_solve, f"optimization_round_{round}_best.json"))
        best_data.update({round: best_round_data})
    
    rounds_cumulative_eligible_trials = {round: round_data["non-pruned trials run"] for round, round_data in best_data.items()}
    #print("rounds_cumulative_eligible_trials", rounds_cumulative_eligible_trials)
    rounds_number_of_trials_missing_until_round_idx = {round: rounds_cumulative_eligible_trials[round_idx]-rounds_cumulative_eligible_trials[round]
                                                       for round in rounds_cumulative_eligible_trials.keys()}
    #print("rounds_number_of_trials_missing_until_round_idx", rounds_number_of_trials_missing_until_round_idx)
    last_round_before_at_least_number_of_trials_to_run_before_stopping = max(
        (k for k, v in rounds_number_of_trials_missing_until_round_idx.items() if v >= number_of_trials_to_run_before_stopping),
        default=None)
    
    # Calculate how much difference/spread in the flux is expected from random trials
    eligible_trials = [
        t for t in study.trials
        if t.state == TrialState.COMPLETE
        and t.user_attrs.get("round", 0) <= round_idx
    ]
    n_startup_trials_values = [trial.value for trial in eligible_trials[:n_startup_trials]]
    flux_standard_devation_among_startup_trials = np.std(n_startup_trials_values)
    
    #print(f"The round we are using for comparing the change is round number {last_round_before_at_least_number_of_trials_to_run_before_stopping}. ")
    
    convergence = True
    # The point is to check that within the (at least) last number_of_trials_to_run_before_stopping number of trials,
    # the result hasn't changed a lot and the values of the parameters haven't changed a lot

    # For convergence, the relative flux increase within the last number_of_trials_to_run_before_stopping number of trials
    #  must be smaller than value_negligible_relative_change
    # For convergence, the position of the inner membranes must be the same (since these are anyways discretized, this is the best we can do)
    #
    ##### IMPORTANT!: once enzymes are added to the mix, the position of the inner membranes must be similar, but not necessarily equal
    #
    rounds_best_values = {round: round_data["best_value"] for round, round_data in best_data.items()}
    rounds_best_inner_membrane_positions = {round: round_data["best_inner_membrane_radii"] for round, round_data in best_data.items()}
    rounds_best_enzyme_allocations = {round: round_data["best_enzyme_allocations"] for round, round_data in best_data.items()}
    rounds_best_regional_alloc = {round: round_data["best_regional_alloc"] for round, round_data in best_data.items()}

    if last_round_before_at_least_number_of_trials_to_run_before_stopping is not None:
        # Values to compare
        previous_best_data_value = rounds_best_values[last_round_before_at_least_number_of_trials_to_run_before_stopping]
        current_best_data_value = rounds_best_values[round_idx]
        # Locations of membranes to compare
        previous_best_inner_membrane_positions = sorted(rounds_best_inner_membrane_positions[last_round_before_at_least_number_of_trials_to_run_before_stopping])
        current_best_inner_membrane_positions = sorted(rounds_best_inner_membrane_positions[round_idx])
        # Allocations of total enzyme quantity to specific enzymes to compare
        previous_best_enzyme_allocations = rounds_best_enzyme_allocations[last_round_before_at_least_number_of_trials_to_run_before_stopping]
        current_best_enzyme_allocations = rounds_best_enzyme_allocations[round_idx]
        # Allocations of enzyme quantities to specific regions to compare
        previous_best_regional_alloc = rounds_best_regional_alloc[last_round_before_at_least_number_of_trials_to_run_before_stopping]
        current_best_regional_alloc = rounds_best_regional_alloc[round_idx]


        # since we are maximizing the flux, current_best_data_value will be equal or larger than previous_best_data_value
        #if (current_best_data_value - previous_best_data_value) / previous_best_data_value > value_negligible_relative_change:
        if current_best_data_value - previous_best_data_value > value_negligible_relative_change * flux_standard_devation_among_startup_trials:
            print("The value has changed too much. Not converged yet.")
            convergence = False
        if previous_best_inner_membrane_positions != current_best_inner_membrane_positions:
            print("The radii are different.")
            convergence = False
        # Within enzyme allocations, the elements in the vectors must be similar element-wise
        # (each element corresponds to the allocation for one enzyme)
        if current_best_enzyme_allocations is not None:
            for enzyme_idx, allocation in enumerate(current_best_enzyme_allocations):
                if abs(allocation - previous_best_enzyme_allocations[enzyme_idx])/allocation > value_negligible_relative_change:
                    convergence = False
        # Within the allocation of each enzyme within each region, the elements of the vectors
        # must be similar element-wise
        for enzyme_idx, enzyme_info in current_best_regional_alloc:
            regions = enzyme_info.keys()
            for region in regions:
                if abs(enzyme_info[region] - previous_best_regional_alloc[enzyme_idx][region]) / enzyme_info[region] > value_negligible_relative_change:
                    convergence = False
        
    else:
        convergence = False
    
    return convergence, last_round_before_at_least_number_of_trials_to_run_before_stopping


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
    optimization_params = read_yaml_file(os.path.join(FOLDER_TO_SOLVE, "parameters_optimization.yaml"))
    if (optimization_params["sampler"]["type"] == "TPESampler"
        and optimization_params["sampler"]["sampler_parameters"]["TPESampler"]["n_startup_trials"] is not None):
        n_startup_trials = optimization_params["sampler"]["sampler_parameters"]["TPESampler"]["n_startup_trials"]
    elif (optimization_params["sampler"]["type"] == "CmaEsSampler"
        and optimization_params["sampler"]["sampler_parameters"]["CmaEsSampler"]["n_startup_trials"] is not None):
        n_startup_trials = optimization_params["sampler"]["sampler_parameters"]["CmaEsSampler"]["n_startup_trials"]
    else:
        n_startup_trials = 20
    
    storage = f"sqlite:///{FOLDER_TO_SOLVE}/optuna_study.db"
    study = optuna.load_study(
        study_name="resource_allocation",
        storage=storage
    )
    # since get_best_result_up_until_round_specified creates the .json file with the best
    # data until then, I want the convergence check to have to happen
    with DelayedSignals():
        result = get_best_result_up_until_round_specified(FOLDER_TO_SOLVE, round_idx, study)
            
        convergence, comparison_round = get_convergence(
            FOLDER_TO_SOLVE,
            round_idx,
            number_of_trials_to_run_before_stopping=optimization_params["convergence_params"]["number_of_trials_to_run_before_stopping"],
            value_negligible_relative_change=optimization_params["convergence_params"]["value_negligible_relative_change"],
            n_startup_trials=n_startup_trials
        )
        convergence_info_path = os.path.join(FOLDER_TO_SOLVE, "optimization_convergence.txt")
        # only create the analysis and plot through here if they have not already been done
        if convergence and not os.path.isfile(convergence_info_path):
            with open(convergence_info_path, "w") as f:
                f.write(f"Optimization converged! Comparing round {round_idx} to round {comparison_round}.\n")
            round_of_best = result["trial_round_best"]
            trial_idx_of_best = result["trial_idx_best"]

            #plot progress up to optimization (so snakemake is happy)
            enzymes_df = pd.read_csv(os.path.join(FOLDER_TO_SOLVE, "enzymes.csv"))
            n_enzymes = len(enzymes_df) # the first row is the header

            data = load_existing_data(FOLDER_TO_SOLVE, study, round_idx, n_trials, n_enzymes)
            geometry_info = read_yaml_file(os.path.join(FOLDER_TO_SOLVE, "parameters_geometry.yaml"))
            n_regions = len(geometry_info["geometry_config"]["internal_membrane_relative_radii"])+1
            #print("loaded_data", data)
            plot_optimization_progress(
                FOLDER_TO_SOLVE,
                data,
                enzymes_df,
                n_regions,
                n_rounds_to_plot=round_idx,
                n_trials = n_trials
            )
            dump_json(FOLDER_TO_SOLVE, "best_result", result)

        if round_idx == n_rounds-1:
            dump_json(FOLDER_TO_SOLVE, "best_result", result)

    