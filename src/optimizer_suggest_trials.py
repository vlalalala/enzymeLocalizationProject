import optuna
import json
import numpy as np
import pandas as pd
import os
import shutil
import argparse
from pathlib import Path
from auxiliary_functions import read_yaml_file, dump_in_yaml_file, dump_json
from optuna.trial import TrialState
from auxiliary_functions_framework_organization_using_standard_library import DelayedKeyboardInterrupt

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

def params_to_physical(
    params, n_enzymes, n_regions,
    total_enzyme_quantity,
    enzyme_maximum_concentration,
    external_radius,
    relative_delta_r = None
):
    n_inner_membranes = n_regions - 1
    prune = {"prune": False}
    # --- Allocation to different enzymes: (n_enzymes - 1) free params ---
    #softmax([-5, 0, 0]) → [0.007, 0.497, 0.497]  # nearly equal split, one tiny
    #softmax([ 5, 0, 0]) → [0.987, 0.007, 0.007]  # almost all in first type
    #softmax([-5, 5, 0]) → [0.000, 0.993, 0.007]  # extreme concentration difference
    # ranges from -5 to 5 makes it possible to represent allocations
    # from nearly 0% to ~99% for any enzyme
    if n_enzymes > 0:
        z_enzymes = [params[f"z_enzyme_types_{i}"] for i in range(n_enzymes - 1)]
        z_enzymes.append(0.0) # last enzyme, removes redundancy when searching space
        enzyme_allocations = (softmax(np.array(z_enzymes)) * total_enzyme_quantity).tolist()
    else:
        enzyme_allocations = None

    # --- Allocation of enzymes to different regions: n_enzymes x (n_regions - 1) free params ---
    regional_alloc = []
    for t in range(n_enzymes):
        z_regions = [params[f"z_region_{t}_{r}"] for r in range(n_regions - 1)]
        z_regions.append(0.0)
        regional_alloc_list = softmax(np.array(z_regions)).tolist()
        regional_alloc.append({
            region: allocation
            for region, allocation in enumerate(regional_alloc_list)
        })

    """
    First thing to check: if the total enzyme amount (sum of all enzyme quantities)
    exceeds the maximum concentration
    times the total volume of the vesicle -> raise an error
    """
    actual_volume = 4/3 * np.pi * external_radius**3
    if total_enzyme_quantity is not None and enzyme_maximum_concentration is not None:
        minimum_required_volume = total_enzyme_quantity / enzyme_maximum_concentration 
        if minimum_required_volume > actual_volume:
            raise ValueError(f"The maximum enzyme concentration is too low. The minimum required volume is {minimum_required_volume} and the actual volume is {actual_volume}")
    else:
        minimum_required_volume = 0
    """
    Calculate the volume "slack" to be distributed between the different regions
    """
    volume_slack = actual_volume - minimum_required_volume
    #print("volumes", volume_slack, actual_volume, minimum_required_volume)
    
    """
    Distribute this volume slack among the different regions
    """
    fraction_extra_volume = [params[f"fraction_extra_volume_{i}"] for i in range(n_regions-1)]
    last_frac = 1.0 - sum(fraction_extra_volume)
    if last_frac <= 0:
        prune["prune"] = True
        prune.update({"reason": "The volumes to allocate between the regions 1 to N-1 already account to more than 100%"})
    # Even if it is pruned, continue
    volume_slack_distribution = np.array(fraction_extra_volume + [last_frac])
    #print("volume_slack_distribution", volume_slack_distribution)
    extra_volume_per_region = (volume_slack_distribution * volume_slack).tolist()
    #print("extra_volume_per_region", extra_volume_per_region)
    
    """
    Find out the minimum volume for each region
    """
    # First calculate the total enzyme amount per region
    total_enzyme_in_region = [0.0] * n_regions
    for enzyme_index in range(n_enzymes):
        for region in range(n_regions):
            total_enzyme_in_region[region] += enzyme_allocations[enzyme_index] * regional_alloc[enzyme_index][region]

    # Then compute the minimum volume per region
    if enzyme_maximum_concentration is not None:
        minimum_volume_per_region = [
            total_enzyme_in_region[region] / enzyme_maximum_concentration
            for region in range(n_regions)
        ]
    else:
        minimum_volume_per_region = [0 for region in range(n_regions)]
    #print("minimum_volume_per_region", minimum_volume_per_region)
    """
    Compute inner membrane radii
    """
    inner_membrane_radii = []
    current_allocated_volume = 0
    for i in range(n_inner_membranes):
        current_allocated_volume += minimum_volume_per_region[i] + extra_volume_per_region[i]
        # (current_allocated_volume*3/(4*np.pi)) ** (1/3) means inner_membrane_radii is in units of meter
        # dividing by external_radius brings inner_membrane_radii into the range between 0 and 1
        inner_membrane_radii.append(
            (current_allocated_volume*3/(4*np.pi)) ** (1/3) / external_radius
        )
    #print("inner_membrane_radii", inner_membrane_radii)
    #print("relative_delta_r", relative_delta_r)
    # Check (only necessary when suggesting the trials. In other cases, do not pass relative_delta_r)
    
    if len(inner_membrane_radii) != 0 and relative_delta_r is not None:
        if inner_membrane_radii[0] < 1.5 * relative_delta_r or inner_membrane_radii[-1] > 1 - relative_delta_r * 1.5:
            prune["prune"] = True
            prune.update({"reason": "The distance between the first or last membrane radii is too close to the limits of 0,1"})
        if not all(b - a >= relative_delta_r * 2 for a, b in zip(inner_membrane_radii[:-1], inner_membrane_radii[1:])):
            prune["prune"] = True
            prune.update({"reason": "A distance between inner membranes is too small"})
        #if not all(0<inner_membrane_radius<1 for inner_membrane_radius in inner_membrane_radii):
        #    raise ValueError("An inner_membrane_radius is not between 0 and 1.")
    return enzyme_allocations, regional_alloc, inner_membrane_radii, prune

def create_files(
        folder_to_solve,
        round_idx,
        trial_idx,
        enzyme_allocations,
        regional_alloc,
        inner_membrane_radii,
        geometry_info,
    ):
    # Copy species.csv, spontaneous_reactions.csv, enzymatic_reactions.csv,
    # parameters_discretization.yaml, parameters_solver_input.yaml,
    # parameters_solver_output.yaml, parameters_value_conditions.yaml
    # onto each trial
    for file in [
        "species.csv", "spontaneous_reactions.csv", "enzymatic_reactions.csv",
        "parameters_discretization.yaml", "parameters_solver_input.yaml",
        "parameters_solver_output.yaml", "parameters_value_conditions.yaml"
    ]:  
        src = os.path.join(folder_to_solve, f"{file}")
        dst = os.path.join(folder_to_solve, f"optimization_round_{round_idx}/trial_{trial_idx}/{file}")
        if not os.path.isfile(src):
            raise FileNotFoundError(f"Source file not found: {src}")
        shutil.copy(src, dst)
        #print(f"Created {dst}")
        if not os.path.isfile(src):
            raise FileNotFoundError(f"Destination file not found: {dst}")
        
    
    # Create a modified parameters_geometry.yaml with the correct membrane radii
    geometry_info["geometry_config"]["internal_membrane_relative_radii"] = inner_membrane_radii
    dump_in_yaml_file(os.path.join(
            folder_to_solve,
            f"optimization_round_{round_idx}/trial_{trial_idx}/parameters_geometry.yaml"),
            geometry_info
    )
    #print(f"Created {os.path.join(folder_to_solve, f"optimization_round_{round_idx}/trial_{trial_idx}/parameters_geometry.yaml")}")
    if enzyme_allocations is not None:
        # Create a modified enzymes.csv with the correct enzyme allocation
        enzymes_df["quantity"] = enzyme_allocations
        enzymes_df["allocation"] = regional_alloc
        enzymes_df.to_csv(
            os.path.join(
                folder_to_solve,
                f"optimization_round_{round_idx}/trial_{trial_idx}/enzymes.csv"),
            index=False)
    else:
        shutil.copy(
            os.path.join(folder_to_solve, "enzymes.csv"),
            os.path.join(folder_to_solve, f"optimization_round_{round_idx}/trial_{trial_idx}/enzymes.csv"),
        )


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
    enzyme_maximum_concentration =  conditions_info["enzyme_maximum_concentration"]  
    # Find out number of regions and minimum distance between membranes and the origin
    geometry_info = read_yaml_file(os.path.join(FOLDER_TO_SOLVE, "parameters_geometry.yaml"))
    n_regions = len(geometry_info["geometry_config"]["internal_membrane_relative_radii"])+1
    external_radius = geometry_info["geometry_config"]["outer_membrane_radius"]
    n_inner_membranes = n_regions - 1
    discretization_info = read_yaml_file(os.path.join(FOLDER_TO_SOLVE, "parameters_discretization.yaml"))
    relative_delta_r = 1/(discretization_info["discretization_parameters"]["min_num_mesh_points"]-1)
    # relative_delta_r is in units of R
    
    optimization_params = read_yaml_file(os.path.join(FOLDER_TO_SOLVE, "parameters_optimization.yaml"))

    storage = f"sqlite:///{FOLDER_TO_SOLVE}/optuna_study.db"
    # Create study on round 0, load on subsequent rounds
    # Sampler in order to have reproducibility on suggestions
    if optimization_params["sampler"]["type"] == "TPESampler":
        sampler = optuna.samplers.TPESampler(
            seed=optimization_params["sampler"]["seed"]+round_idx,
            n_startup_trials=optimization_params["sampler"]["sampler_parameters"]["TPESampler"]["n_startup_trials"],
            gamma=lambda x: int(optimization_params["sampler"]["sampler_parameters"]["TPESampler"]["gamma"] * x)
        )# gamma is the exploration vs exploitation parameter -> the lower gamma, the more exploitation,
        # aka sampling close to the current optimal value. Default is 0.25
    elif optimization_params["sampler"]["type"] == "CmaEsSampler":
        sampler = optuna.samplers.CmaEsSampler(
            seed=optimization_params["sampler"]["seed"]+round_idx,
            sigma0 = optimization_params["sampler"]["sampler_parameters"]["CmaEsSampler"]["sigma0"],
            restart_strategy = optimization_params["sampler"]["sampler_parameters"]["CmaEsSampler"]["restart_strategy"],
            n_startup_trials = optimization_params["sampler"]["sampler_parameters"]["CmaEsSampler"]["n_startup_trials"],
        )
    # Claude says this is better if only 1 parameter to tune. Need to test...
    else:
        raise ValueError("optimization sampler not specified")
    # If I only do seed=42, all of the trials are exactly the same... :(
    study = optuna.create_study(
        study_name="resource_allocation",
        storage=storage,
        direction="maximize",
        load_if_exists=True,
        sampler=sampler
    )
    
    #existing = [
    #        t for t in study.trials
    #        if t.user_attrs.get("round") == 0
    #        and t.user_attrs.get("trial_idx") == 2
    #    ]
    #print("previous params: ", existing[0].params)
    # Ask Optuna for a batch of N_TRIALS suggestions
    for trial_idx in range(n_trials):
        existing = [
            t for t in study.trials
            if t.user_attrs.get("round") == round_idx
            and t.user_attrs.get("trial_idx") == trial_idx
        ]
        if len(existing) > 1:
            raise ValueError(f"Something went wrong. More than one trials have round {round_idx} and trial {trial_idx}.")
        if existing:
            trial = existing[-1]
        else:
            ########## IMPORTANT: if some keyboard interrupt occurs between
            # study.ask and the last .suggest_float, downstream stuff won't work...
            # Therefore: DelayedKeyboardInterrupt wrap
            with DelayedKeyboardInterrupt():
                trial = study.ask()  # get a suggested trial without evaluating it
                # access through study.trials[trial_idx]
                trial.set_user_attr("round", round_idx)
                trial.set_user_attr("trial_idx", trial_idx)
                
                # --- Allocation of total enzyme quantity to different enzymes: (n_enzymes - 1) free params ---
                if n_enzymes > 0:
                    z_enzymes = [trial.suggest_float(f"z_enzyme_types_{i}", -5, 5)
                            for i in range(n_enzymes - 1)]

                # --- Allocation of enzymes to different regions: n_enzymes x (n_regions - 1) free params ---
                for t in range(n_enzymes):
                    z_regions = [trial.suggest_float(f"z_region_{t}_{r}", -5, 5)
                                for r in range(n_regions - 1)]
                
                # --- Allocation of volume that is left from most packed distribution of enzymes     
                # z_i taken from Uniform(-5,5). After applying the fixed logit 0 and the softmax,
                # e^-5 = 0.0067 and e^5 = 148
                # p_i = e^{z_i} / sum_j e^{z_j}
                # for volume, best to pick fractions of volume freely (without drawing from the distribution
                # and then applying the softmax), which would require extreme logit values.
                # The number of free parameters is n_regions-1, since the extra volume should add up to 1
                fraction_extra_volume = [trial.suggest_float(f"fraction_extra_volume_{i}", 0.0, 1.0)
                    for i in range(n_regions-1)]
                #print("z_extra_volume:", z_extra_volume)
        #print(trial.params, round_idx, trial_idx)
        enzyme_allocations, regional_alloc, inner_membrane_radii, prune = params_to_physical(
            trial.params,
            n_enzymes=n_enzymes,
            n_regions=n_regions,
            total_enzyme_quantity=total_enzyme_quantity,
            enzyme_maximum_concentration=enzyme_maximum_concentration,
            external_radius=external_radius,
            relative_delta_r = relative_delta_r
        )
        if os.path.isfile(os.path.join(FOLDER_TO_SOLVE, "optimization_convergence.txt")):
            prune["prune"] = True
            prune.update({"reason": "The optimization procedure has already converged."})

        result = {
            "enzyme_allocations": enzyme_allocations,
            "regional_alloc": regional_alloc,
            "inner_membrane_radii": inner_membrane_radii,
            "prune": prune
        }
    
        if prune["prune"]:
            trial_dir = os.path.join(FOLDER_TO_SOLVE, f"optimization_round_{round_idx}/trial_{trial_idx}")
            #print(f"Infeasable {trial_dir} or system already optimized:  {prune["reason"]}, pruning.")
            current_trial = study.trials[trial.number]
            if current_trial.state == TrialState.RUNNING:
                study.tell(trial.number, state=TrialState.PRUNED)
            # already create output files of rules
            dump_json(trial_dir, "pruned", prune)
            # Define result placeholder for snakemake

        # Create files regardless of pruning so Snakemake is satisfied
        dump_json(
            os.path.join(FOLDER_TO_SOLVE, f"optimization_round_{round_idx}/trial_{trial_idx}"), "trial_parameters",
            result
        )
        create_files(
            folder_to_solve=FOLDER_TO_SOLVE,
            round_idx=round_idx,
            trial_idx=trial_idx,
            enzyme_allocations=enzyme_allocations,
            regional_alloc=regional_alloc,
            inner_membrane_radii=inner_membrane_radii,
            geometry_info=geometry_info
        )
    
    print(f"Created trial files for {os.path.join(FOLDER_TO_SOLVE, f"optimization_round_{round_idx}")}")







