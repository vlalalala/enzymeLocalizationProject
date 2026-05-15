import yaml
from pathlib import Path
import pandas as pd

def filter_combined_folders_yaml_criteria(combined_root: str | Path, criteria_per_option_yaml: dict):
    """
    combined_root: path containing options_* and combined_* folders
    criteria: dict like {"options_parameters_geometry": {"geometry_config": {"outer_membrane_radius": 0.5}},
                          "options_parameters_solver_input": {"some_param": 10}}
    Returns: list of Path objects to combined folders satisfying all criteria (combined_*)

    Example:
    filter_combined_folders(
        "../data_private/slurm_test2",
        {
            "options_parameters_geometry": {
                "geometry_config": {
                    "outer_membrane_radius": 1.0e-5
                }
            },
            "options_parameters_solver_input": {
                "geometry_parameters" : {
                    "num_mesh_points":400
                }
            }
        }
    )
    """
    combined_root = Path(combined_root)
    result = []

    for folder in sorted(combined_root.glob("combined_*")):
        match = True
        for options_folder_name, param_dict in criteria_per_option_yaml.items():

            # Extract the '*' from 'options_*'
            if not options_folder_name.startswith("options_"):
                raise ValueError(f"{options_folder_name} must start with 'options_'")

            basename = options_folder_name[len("options_"):]
            opt_file = folder / f"{basename}.yaml"
            
            if not opt_file.exists():
                match = False
                break
            with open(opt_file) as f:
                data = yaml.safe_load(f)
            # recursive check
            def check_params(d, p):
                for k, v in p.items():
                    if isinstance(v, dict):
                        if not check_params(d.get(k, {}), v):
                            return False
                    else:
                        if d.get(k) != v:
                            return False
                return True
            if not check_params(data, param_dict):
                match = False
                break
        if match:
            result.append(folder)
    return result

def csv_file_matches_criteria(csv_path: str | Path, conditions: dict) -> bool:
    """
    Example:
    csv_file_matches_criteria(
        "../data_private/slurm_test2/options_species/combo_000001.csv",
        {"name": {"Z": {"external_concentration": 0}}}
    )

    Returns True if the file satisfies all conditional constraints.
    (if there are no conditional constraints, it returns True)
    """
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"The file {csv_path} does not exist.")
    df = pd.read_csv(csv_path, dtype=str)  # keep everything as string
    for trigger_column, trigger_map in conditions.items():
        if trigger_column not in df.columns:
            return False
        for trigger_value, required_dict in trigger_map.items():
            # Select rows where trigger_column == trigger_value
            mask = df[trigger_column] == str(trigger_value)
            if not mask.any():
                # No such rows → nothing to check
                continue
            for required_column, required_value in required_dict.items():
                if required_column not in df.columns:
                    return False
                # Check that ALL selected rows satisfy the requirement
                N = 5  # number of leading characters to compare
                if not (
                    df.loc[mask, required_column]
                    .astype(str)
                    .str[:N]
                    .eq(str(required_value)[:N])
                    .all()
                ):
                    return False
    return True

def filter_combined_folders_csv_criteria(
    combined_root: str | Path,
    criteria_per_option_csv: dict,
):
    """
    Example:
    filter_combined_folders_csv_criteria(
        "../data_private/slurm_test2",
        {
            "options_enzymes": {"name": {"A": {"regions": "[0]"}}}
        }
    )
    """
    combined_root = Path(combined_root)
    result = []
    for folder in sorted(combined_root.glob("combined_*")):
        match = True
        for options_folder_name, conditions in criteria_per_option_csv.items():
            if not options_folder_name.startswith("options_"):
                raise ValueError(
                    f"{options_folder_name} must start with 'options_'"
                )
            # Extract the * part
            basename = options_folder_name[len("options_"):]
            csv_file = folder / f"{basename}.csv"
            if not csv_file.exists():
                match = False
                break
            if not csv_file_matches_criteria(csv_file, conditions):
                match = False
                break
        if match:
            result.append(folder)
    return result

######################################################################################
# Main function to use
def filter_combined_folders(
    combined_root: str | Path,
    criteria_yaml,
    criteria_csv,
):
    csv_matches = filter_combined_folders_csv_criteria(
        combined_root, criteria_csv
    )
    yaml_matches = filter_combined_folders_yaml_criteria(
        combined_root, criteria_yaml
    )
    return list(set(csv_matches) & set(yaml_matches))
#########################################################################################
# %%
####################### obsolete ####################################
def get_options_dict(folder: str | Path, file_ext=".yaml"):
    #index_options_folder
    """
    Returns a (nested) dictionary for the given folder
    where the first key is the name of each file with the given extension within that folder
    e.g.
    geometry_dict = get_options_dict("options_parameters_geometry")
    geometry_dict["combo_000002.yaml"]["geometry_config"]["outer_membrane_radius"] → 0.5
    """
    folder = Path(folder)
    if not folder.exists():
        raise FileNotFoundError(f"The folder {folder} does not exist.")
    options_info = {}
    for f in folder.glob(f"*{file_ext}"):
        with open(f) as fd:
            data = yaml.safe_load(fd)
        options_info[f.name] = data
    return options_info

def filter_files_by_param(options_dict, param_path, target_value):
    """
    Returns the file names (without the path) of the files where the target value is found
    within a specific file within a specific path of nested dictionaries.

    options_dict is the output of the function get_options_dict
    param_path: list of nested keys, e.g., ["geometry_config", "outer_membrane_radius"]
    target_value: the value we want to filter by
    
    e.g. filter_files_by_param(geometry_dict, ["geometry_config", "outer_membrane_radius"], 0.5)
    can return e.g. ["combo_000002"]
    """
    result = []
    for fname, data in options_dict.items():
        value = data
        try:
            for key in param_path:
                value = value[key]
        except KeyError:
            continue
        if value == target_value:
            result.append(fname)
    return result


def filter_option_csv_folder(options_folder: str, conditions: dict, combo_prefix="combo_"):
    """ Returns a list of all csv files within the options folder
    that conform to the conditions.
    Example:
    filter_option_csv_folder(
        "../data_private/slurm_test2/options_enzymes",
        {"name": {"A": {"regions": "[0]"}}}
    )
    """
    folder = Path(options_folder)
    matches = []
    for f in folder.glob(f"{combo_prefix}*.csv"):
        if csv_file_matches_criteria(f, conditions):
            matches.append(f.name)
    return matches
