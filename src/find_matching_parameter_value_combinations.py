import yaml
from pathlib import Path

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


def filter_combined_folders(combined_root: str | Path, criteria: dict):
    """
    combined_root: path containing options_* and combined_* folders
    criteria: dict like {"options_parameters_geometry": {"geometry_config": {"outer_membrane_radius": 0.5}},
                          "options_parameters_solver_input": {"some_param": 10}}
    Returns: list of Path objects to combined folders satisfying all criteria
    """
    combined_root = Path(combined_root)
    result = []

    for folder in sorted(combined_root.glob("combined_*")):
        match = True
        for opt_basename, param_dict in criteria.items():
            opt_file = folder / f"{opt_basename}.yaml"
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



def csv_file_matches_criteria(csv_path: str, conditions: dict) -> bool:
    """
    csv_path: path to CSV file
    conditions: dict mapping column to dict of condition value → required value
        Example: {
            "name": {"A": "R_alpha", "B": "R_beta"},
        }
        Meaning: if row[name] == "A" → row["regions"] == "R_alpha"
                 if row[name] == "B" → row["regions"] == "R_beta"
    Example:
    conditions = {
        "name": {
            "A": "R_alpha",
            "B": "R_beta",
        }
    }

    matches = csv_file_matches_criteria("options_foo/combo_00001.csv", conditions)
    print(matches)  # True/False
    """
    df = pd.read_csv(csv_path, dtype=str)  # all as strings
    for col, mapping in conditions.items():
        for val, required_region in mapping.items():
            mask = df[col] == val
            if not df.loc[mask, "regions"].eq(required_region).all():
                return False
    return True

def filter_option_csv_folder(folder: str, conditions: dict, combo_prefix="combo_"):
    folder = Path(folder)
    matches = []
    for f in folder.glob(f"{combo_prefix}*.csv"):
        if csv_file_matches_criteria(f, conditions):
            matches.append(f.name)
    return matches



def filter_combined_csv_folders(combined_root: str, criteria_per_option_csv: dict):
    """
    criteria_per_option_csv: dict mapping option folder basename → filter dict
    Example:
        {
            "foo": {"name": {"A": "R_alpha", "B": "R_beta"}},
            "bar": {"some_col": {"X": "Y"}}
        }

        criteria_per_option_csv = {
        "foo": {"name": {"A": "R_alpha", "B": "R_beta"}},
        "bar": {"id": {"X": "1"}},  # could mix CSV and YAML similarly
    }

    matching_folders = filter_combined_csv_folders("all_combos", criteria_per_option_csv)
    print(matching_folders)
    """
    combined_root = Path(combined_root)
    result = []

    for folder in sorted(combined_root.glob("combined_*")):
        match = True
        for option_basename, conditions in criteria_per_option_csv.items():
            csv_file = folder / f"{option_basename}.csv"
            if not csv_file.exists() or not csv_file_matches_criteria(csv_file, conditions):
                match = False
                break
        if match:
            result.append(folder)
    return result

"""
criteria = {
    "alpha": {"geometry_config": {"outer_membrane_radius": 0.5}},
    "beta": {"some_param": 10},
}

matching_folders = filter_combined_folders("all_combos", criteria)
print("Folders that match criteria:", matching_folders)
"""