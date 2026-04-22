import argparse
import os
from auxiliary_functions_using_standard_library import load_json
from auxiliary_functions import dump_json
from pathlib import Path
import shutil
import bisect
from auxiliary_functions import read_yaml_file, dump_in_yaml_file
import copy
import pandas as pd
import ast

def get_user_input_file_names(folder):
    file_names = [file
        for file in os.listdir(folder)
        if (os.path.isfile(os.path.join(folder, file))
            and (file.endswith(".csv") or file.endswith(".yaml"))
            and "progress" not in file)
    ]
    return file_names

def copy_file_names(file_names, source, destination):
    for file in file_names:
        shutil.copy(os.path.join(source, file),
                    os.path.join(destination, file))
        if not os.path.isfile(os.path.join(destination, file)):
            raise FileNotFoundError(f"Shutil failed silently.")

def copy_best_result_onto_check_folder(folder_to_solve, optimization_check_folder):
    os.makedirs(optimization_check_folder, exist_ok=True)
    # Find the folder that produced the best results
    best_results = load_json(os.path.join(folder_to_solve, "best_result.json"))
    path_to_best_result = Path(folder_to_solve) / f"optimization_round_{best_results["trial_round_best"]}" / f"trial_{best_results["trial_idx_best"]}"
    # copy all the .csv files and .yaml
    # (except the .csv files tracking progress) onto the optimization check folder,
    # as well as .system_geometry.json
    file_names_to_copy = get_user_input_file_names(path_to_best_result)
    file_names_to_copy += ["system_geometry_for_convergence.json"]
    copy_file_names(
        file_names = file_names_to_copy,
        source=path_to_best_result,
        destination=optimization_check_folder
    )

def neighbors(ordered_list, central_value):
    """Gets the immediately smaller and immediately larger values
    than a central value within an ordered list.
    """
    i = bisect.bisect_left(ordered_list, central_value)
    lower = ordered_list[i-1] if i > 0 else None
    upper = ordered_list[i+1] if i < len(ordered_list)-1 and ordered_list[i] == central_value else (
        ordered_list[i] if i < len(ordered_list) else None
    )
    return lower, upper

def create_optimization_combinations_differing_membrane_location(
    optimization_check_folder,
    optimization_combinations_up_until_now
    ):
    # Calculate new radii for each radius
    optimal_parameters_geometry = read_yaml_file(os.path.join(optimization_check_folder, "parameters_geometry.yaml"))
    system_geometry = load_json(os.path.join(optimization_check_folder, "system_geometry_for_convergence.json"))
    optimal_external_radius = system_geometry["geometry_config"]["outer_membrane_radius"]
    optimal_internal_membrane_radii = system_geometry["geometry_config"]["membrane_radii"][:-1]
    baseline_mesh_points = system_geometry["geometry_config"]["baseline_mesh_points"]
    for optimal_internal_membrane_radius_idx, optimal_internal_membrane_radius in enumerate(
            optimal_internal_membrane_radii):
        left_radius, right_radius = neighbors(baseline_mesh_points, optimal_internal_membrane_radius)
        for new_radius in [left_radius, right_radius]:
            new_relative_radius = new_radius/optimal_external_radius
            parameters_geometry = copy.deepcopy(optimal_parameters_geometry)
            parameters_geometry["geometry_config"]["internal_membrane_relative_radii"][optimal_internal_membrane_radius_idx] = new_relative_radius
            # Once the new parameters_geometry is done, copy all of the .csv and .yaml files (except the geometry one)
            file_names_to_copy = get_user_input_file_names(optimization_check_folder)
            file_names_to_copy.remove("parameters_geometry.yaml")
            new_modification_folder = os.path.join(optimization_check_folder, f"modification_{optimization_combinations_up_until_now}")
            os.makedirs(new_modification_folder, exist_ok=True)
            copy_file_names(
                file_names = file_names_to_copy,
                source=optimization_check_folder,
                destination=new_modification_folder
            )
            # Dump modified geometry file
            dump_in_yaml_file(os.path.join(new_modification_folder, "parameters_geometry.yaml"), parameters_geometry)
            if new_radius == left_radius:
                change = "decreased"
            else:
                change = "increased"
            with open(os.path.join(new_modification_folder, "info_on_modification.txt"), "w") as f:
                f.write(f"The radius of the inner membrane with index {optimal_internal_membrane_radius_idx} (starting from index 0, from left to right) was {change}.\n")
            # In order to keep folder names clear, track how many modification folders have been created
            optimization_combinations_up_until_now += 1
    return optimization_combinations_up_until_now

def create_optimization_combinations_differing_allocation_of_total_enzyme_quantity_to_different_enzymes(
    optimization_check_folder,
    optimization_combinations_up_until_now
    ):
    """
    Through this function we create modifications of the quantity of each enzyme, maintaining
    the total quantity
    """
    enzymes_df = pd.read_csv(os.path.join(optimization_check_folder, "enzymes.csv"))
    total_enzyme_quantity = sum(enzymes_df["quantity"])
    for enzyme_to_modify in enzymes_df["name"].unique():
        enzyme_quantity_among_the_other_enzymes = (
            total_enzyme_quantity
            - enzymes_df.loc[enzymes_df["name"] == enzyme_to_modify, "quantity"])
        for modification in ["increase", "decrease"]:
            enzymes_df_to_modify = copy.deepcopy(enzymes_df)
            if modification == "increase":
                factor = 1.01
                change = "increased"
            else:
                factor = 0.99
                change = "decreased"
            # Modify the value for that enzyme
            enzymes_df_to_modify.loc[enzymes_df_to_modify["name"] == enzyme_to_modify, "quantity"] *= factor
            # Now modify the value for the other enzymes so that the total quantity is maintained
            new_enzyme_quantity_among_the_other_enzymes = (
                total_enzyme_quantity
                - enzymes_df_to_modify.loc[enzymes_df_to_modify["name"] == enzyme_to_modify, "quantity"])
            pruned = False
            pruning_reason = ""
            #print("new", new_enzyme_quantity_among_the_other_enzymes.item()) # without the .item() it is a Series
            if (new_enzyme_quantity_among_the_other_enzymes.item() < 0 
                # aka the modified quantity for the enzyme is larger than the total amount
                or enzymes_df_to_modify.loc[enzymes_df_to_modify["name"] == enzyme_to_modify, "quantity"].item() < 0
                # aka the modified quantity of the enzyme is smaller than 0
                or new_enzyme_quantity_among_the_other_enzymes.item() > 0 and len(enzymes_df) == 1
                # aka the amount of enzyme was decreased and there are no other enzymes to "pick up" the missing quantity
                ):
                pruned = True
                pruning_reason = f"Enzyme {enzyme_to_modify} has a quantity below 0 or larger than the total enzyme quantity or has a smaller amount than the total enzyme quantity and there is only one enzyme."

            # The ratio among the other enzymes has to be maintained. Therefore, we change
            # the quantity of those enzymes by the same ratio, so as to "fill in" the new quantity
            # that is left for these other enzymes
            ratio = new_enzyme_quantity_among_the_other_enzymes / enzyme_quantity_among_the_other_enzymes
            for other_enzyme in enzymes_df_to_modify["name"].unique():
                if other_enzyme == enzyme_to_modify:
                    continue
                enzymes_df_to_modify.loc[enzymes_df_to_modify["name"] == other_enzyme, "quantity"] *= ratio
            
            # Once the new enzymes_df is done, copy all of the .csv and .yaml files (except the enzymes one)
            file_names_to_copy = get_user_input_file_names(optimization_check_folder)
            file_names_to_copy.remove("enzymes.csv")
            new_modification_folder = os.path.join(optimization_check_folder, f"modification_{optimization_combinations_up_until_now}")
            os.makedirs(new_modification_folder, exist_ok=True)
            copy_file_names(
                file_names = file_names_to_copy,
                source=optimization_check_folder,
                destination=new_modification_folder
            )
            # Dump modified enzymes file
            enzymes_df_to_modify.to_csv(os.path.join(new_modification_folder, "enzymes.csv"), encoding='utf-8-sig', index=False)
            if pruned:
                dump_json(new_modification_folder, "pruned", {"pruned": True, "reason": pruning_reason})
            with open(os.path.join(new_modification_folder, "info_on_modification.txt"), "w") as f:
                f.write(f"The quantity of enzyme {enzyme_to_modify} was {change}.\n")
            # In order to keep folder names clear, track how many modification folders have been created
            optimization_combinations_up_until_now += 1
    return optimization_combinations_up_until_now

################
# Probably need to add code to make it possible for pruning of case due to allocation not between 0 and 1
# without breaking the code
#################

def create_optimization_combinations_differing_allocation_of_enzyme_quantities_to_different_regions(
    optimization_check_folder,
    optimization_combinations_up_until_now
    ):
    """
    Through this function we create modifications of the quantity of each enzyme, maintaining
    the total quantity
    """
    enzymes_df = pd.read_csv(os.path.join(optimization_check_folder, "enzymes.csv"))
    for enzyme_to_modify in enzymes_df["name"].unique():
        allocation_dict = ast.literal_eval(enzymes_df.loc[enzymes_df["name"] == enzyme_to_modify, "allocation"].values[0])
        for region_to_modify in allocation_dict.keys():
            for modification in ["increase", "decrease"]:
                allocation_dict_to_modify = copy.deepcopy(allocation_dict)
                enzymes_df_to_modify = copy.deepcopy(enzymes_df)
                if modification == "increase":
                    factor = 1.01
                    change = "increased"
                else:
                    factor = 0.99
                    change = "decreased"

                # Modify the target region
                allocation_dict_to_modify[region_to_modify] *= factor

                # Compute ratio to renormalize the other regions
                previous_allocation_among_other_regions = 1 - allocation_dict[region_to_modify]
                target_allocation_among_other_regions = 1 - allocation_dict_to_modify[region_to_modify]
                ratio = target_allocation_among_other_regions / previous_allocation_among_other_regions

                # Renormalize the other regions
                for other_region in allocation_dict_to_modify.keys():
                    if other_region == region_to_modify:
                        continue
                    allocation_dict_to_modify[other_region] *= ratio

                if sum(allocation_dict_to_modify.values()) != 1.0:
                    raise ValueError(f"The sum is unequal to 1: {sum(allocation_dict_to_modify.values())}")                # Very important: if any of the values is not between 0 and 1 (0 and 1 inclusive), then state pruned
                pruned = False
                pruning_reason = ""
                for region, value in allocation_dict_to_modify.items():
                    if value < 0 or value > 1:
                        pruned = True
                        pruning_reason = f"Region {region} has a non-valid quantity (i.e. not within [0,1])."
                # Insert changed dictionary into dataframe
                enzymes_df_to_modify.loc[enzymes_df_to_modify["name"] == enzyme_to_modify, "allocation"] = str(allocation_dict_to_modify)

                # Once the new allocations dictionary is done, copy all of the .csv and .yaml files (except the enzymes one)
                file_names_to_copy = get_user_input_file_names(optimization_check_folder)
                file_names_to_copy.remove("enzymes.csv")
                new_modification_folder = os.path.join(optimization_check_folder, f"modification_{optimization_combinations_up_until_now}")
                os.makedirs(new_modification_folder, exist_ok=True)
                copy_file_names(
                    file_names = file_names_to_copy,
                    source=optimization_check_folder,
                    destination=new_modification_folder
                )
                # Dump modified enzymes file
                enzymes_df_to_modify.to_csv(os.path.join(new_modification_folder, "enzymes.csv"), encoding='utf-8-sig', index=False)
                if pruned:
                    dump_json(new_modification_folder, "pruned", {"pruned": True, "reason": pruning_reason})
                with open(os.path.join(new_modification_folder, "info_on_modification.txt"), "w") as f:
                    f.write(f"The quantity of enzyme {enzyme_to_modify} in region {region_to_modify} was {change}. {pruning_reason}\n")
                # In order to keep folder names clear, track how many modification folders have been created
                optimization_combinations_up_until_now += 1
    return optimization_combinations_up_until_now










if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--folder", type=str, required=True)
    parser.add_argument("--expected_number_modifications", type=int, required=True)
    args = parser.parse_args()
    FOLDER_TO_SOLVE = args.folder
    OPTIMIZATION_CHECK_FOLDER = os.path.join(FOLDER_TO_SOLVE, "optimization_check")
    # here, in case this code is run by snakemake, the optimization check folder should be 
    # created from scratch. Delete if it already exists
    if os.path.isdir(OPTIMIZATION_CHECK_FOLDER):
        shutil.rmtree(OPTIMIZATION_CHECK_FOLDER)
    
    EXPECTED_NUMBER_MODIFICATIONS = args.expected_number_modifications
    
    copy_best_result_onto_check_folder(FOLDER_TO_SOLVE, OPTIMIZATION_CHECK_FOLDER)

    optimization_combinations_up_until_now = 0

    optimization_combinations_up_until_now = create_optimization_combinations_differing_membrane_location(
        OPTIMIZATION_CHECK_FOLDER,
        optimization_combinations_up_until_now
    )

    optimization_combinations_up_until_now = create_optimization_combinations_differing_allocation_of_total_enzyme_quantity_to_different_enzymes(
        OPTIMIZATION_CHECK_FOLDER,
        optimization_combinations_up_until_now
    )

    optimization_combinations_up_until_now = create_optimization_combinations_differing_allocation_of_enzyme_quantities_to_different_regions(
        OPTIMIZATION_CHECK_FOLDER,
        optimization_combinations_up_until_now
    )

    if EXPECTED_NUMBER_MODIFICATIONS != optimization_combinations_up_until_now:
        shutil.rmtree(OPTIMIZATION_CHECK_FOLDER)
        raise ValueError(f"The expected number of modifications {EXPECTED_NUMBER_MODIFICATIONS} does not match the number of created folders {optimization_combinations_up_until_now}. Stopping snakemake from continuing.")
    
    


