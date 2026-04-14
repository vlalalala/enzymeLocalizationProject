import argparse
import os
from auxiliary_functions_using_standard_library import load_json, dump_json_base
from pathlib import Path
from auxiliary_functions import read_yaml_file

def get_fluxes(optimization_check_folder):
    fluxes = {}
    for folder in sorted(Path(optimization_check_folder).glob("modification_*")):
        if folder.is_dir():
            index = int(folder.name.split("_")[1])
            fluxes[index] = load_json(folder / "fluxes.json")
    return fluxes

def categorize_flux_changes(data: dict) -> dict:
    decrease = []
    increase = []
    same = []
    null = []

    for key, value in data.items():
        flux = value[0]
        if flux is None:
            null.append(key)
        elif flux < 0:
            decrease.append(key)
        elif flux > 0:
            increase.append(key)
        else:
            same.append(key)

    return {
        "flux decrease under modification": [len(decrease), decrease],
        "flux increase under modification": [len(increase), increase],
        "flux stays the same": [len(same), same],
        "null flux": [len(null), null]
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--folder", type=str, required=True)
    args = parser.parse_args()
    FOLDER_TO_SOLVE = args.folder
    OPTIMIZATION_CHECK_FOLDER = os.path.join(FOLDER_TO_SOLVE, "optimization_check")

    parameters_optimization = read_yaml_file(os.path.join(FOLDER_TO_SOLVE, "parameters_optimization.yaml"))
    species_to_maximize = parameters_optimization["species_to_maximize"]
    
    modified_fluxes = get_fluxes(OPTIMIZATION_CHECK_FOLDER)
    best_result_data = load_json(os.path.join(FOLDER_TO_SOLVE, "best_result.json"))
    original_flux = best_result_data["best_value"]

    ####################
    # Show comparison
    ####################
    comparison = {}
    for modification_idx, flux_data in modified_fluxes.items():
        if "pruned" in flux_data.keys():
            relative_change_in_flux = None
        else:
            modified_flux = flux_data[species_to_maximize]
            relative_change_in_flux = (modified_flux - original_flux) / original_flux
        with open(Path(OPTIMIZATION_CHECK_FOLDER) / f"modification_{modification_idx}" / "info_on_modification.txt") as f:
            info = f.read()
        comparison.update({modification_idx: (relative_change_in_flux, info)})

    dump_json_base(FOLDER_TO_SOLVE, "optimization_modifications_analysis", comparison)

    ##################
    # Show summary of comparison
    ##################
    categorization = categorize_flux_changes(comparison)
    dump_json_base(FOLDER_TO_SOLVE, "optimization_modifications_summary", categorization)

    




    
    



