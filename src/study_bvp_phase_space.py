#%%
from find_matching_parameter_value_combinations import filter_combined_folders
from auxiliary_functions_using_standard_library import load_json
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import numpy as np

columns_to_specify = {
    "options_species": ["name"],
    "options_spontaneous_reactions": ["start_species", "end_species"],
    "options_enzymes": ["name"],
    "options_enzymatic_reactions": ["start_species", "end_species", "enzyme"]
}

def get_2D_phase_space(combined_root, yaml_criteria, csv_criteria, x_axis, y_axis, file_name_to_plot):
    """
    combined_root, yaml_criteria, csv_criteria specify the (joint) criteria of all
    points that go into the diagram.
    x_axis and y_axis take the form
        ["options_species", {"name": "X"}, "k"]
    file_name_to_plot is the name of the .json file within each combined 
    """
    # Filter folders
    folders_to_plot = filter_combined_folders(
        combined_root, yaml_criteria, csv_criteria)
    # Save results to plot as a list of dictionaries of (folder, x_value, y_value, dict_within_file_name_to_plot)
    results_to_plot = []
    # Validate axis structure (existence of column to plot not checked yet)
    for axis in (x_axis, y_axis):
        option_folder, selector_dict, value_column = axis
        if option_folder not in columns_to_specify:
            raise ValueError(f"{option_folder} not in columns_to_specify")
        required_keys = columns_to_specify[option_folder]
        if set(selector_dict.keys()) != set(required_keys):
            raise ValueError(
                f"For {option_folder}, selector_dict must contain exactly "
                f"{required_keys}, but got {list(selector_dict.keys())}"
            )
    # Iterate over folders
    for folder in folders_to_plot:
        folder = Path(folder)
        
        def extract_axis_value(axis):
            option_folder, selector_dict, value_column = axis
            basename = option_folder[len("options_"):]
            csv_path = folder / f"{basename}.csv"
            if not csv_path.exists():
                raise FileNotFoundError(f"{csv_path} not found")
            df = pd.read_csv(csv_path, dtype=str)
            # build mask from selector_dict
            mask = pd.Series(True, index=df.index)
            for col, val in selector_dict.items():
                mask &= df[col] == str(val)
            matching_rows = df.loc[mask]
            if len(matching_rows) != 1:
                raise ValueError(
                    f"Expected exactly 1 matching row in {csv_path}, "
                    f"but found {len(matching_rows)}"
                )
            return matching_rows.iloc[0][value_column]
        
        x_value = extract_axis_value(x_axis)
        y_value = extract_axis_value(y_axis)

        # Open json file
        json_path = folder / file_name_to_plot
        try:
            result = load_json(json_path)
            results_to_plot.append({
                    "folder": folder,
                    "x": x_value,
                    "y": y_value,
                    "json": result,
                }
            )
        except:
            print(f"Could not get {json_path} info.")

    return results_to_plot    

def build_phase_space_matrix(
    phase_points,
    json_key,
    duplicate_mode="error",
):
    """
    phase_points: output of get_2D_phase_space()
    json_key: key inside each JSON dict to extract
    duplicate_mode:
        "error"  -> raise error if duplicate (x,y)
        "first"  -> keep first occurrence
        "mean"   -> average duplicates
    """
    # Extract (x, y, z)
    raw_points = []
    for entry in phase_points:
        x = float(entry["x"])
        y = float(entry["y"])
        json_data = entry["json"]
        if json_key not in json_data:
            raise KeyError(
                f"{json_key} not found in JSON for folder {entry['folder']}"
            )
        z = float(json_data[json_key])
        raw_points.append((x, y, z))

    # Handle duplicates
    point_dict = {}
    for x, y, z in raw_points:
        key = (x, y)
        point_dict.setdefault(key, []).append(z)
    processed_points = []
    for (x, y), z_list in point_dict.items():
        if len(z_list) > 1:
            if duplicate_mode == "error":
                raise ValueError(f"Duplicate point at {(x, y)}")
            elif duplicate_mode == "first":
                z_val = z_list[0]
            elif duplicate_mode == "mean":
                z_val = np.mean(z_list)
            else:
                raise ValueError("Invalid duplicate_mode")
        else:
            z_val = z_list[0]
        processed_points.append((x, y, z_val))

    # Sort points
    processed_points.sort(key=lambda t: (t[0], t[1]))

    # Build structured grid
    x_values = sorted(set(p[0] for p in processed_points))
    y_values = sorted(set(p[1] for p in processed_points))
    X_unique = np.array(x_values)
    Y_unique = np.array(y_values)
    Z_grid = np.full((len(Y_unique), len(X_unique)), np.nan)

    x_index = {v: i for i, v in enumerate(X_unique)}
    y_index = {v: i for i, v in enumerate(Y_unique)}

    for x, y, z in processed_points:
        xi = x_index[x]
        yi = y_index[y]
        Z_grid[yi, xi] = z

    return {
        "points": processed_points,   # sorted list (x,y,z)
        "X_unique": X_unique,
        "Y_unique": Y_unique,
        "Z_grid": Z_grid,
    }

def plot_2D_phase_space(
    phase_matrix,
    path_to_save,
    x_label=None, y_label=None,
    log_x=False, log_y=False,
    cmap = "viridis"
    ):
    fig, ax = plt.subplots(1,1, figsize = (3,3))
    
    X = phase_matrix["X_unique"]
    Y = phase_matrix["Y_unique"]
    Z = phase_matrix["Z_grid"]
    X_mesh, Y_mesh = np.meshgrid(X, Y)

    c = ax.pcolormesh(
        X_mesh,
        Y_mesh,
        Z,
        shading="auto",
        cmap=cmap,
    )
    # Scatterplot for actual values
    x_points = [float(p[0]) for p in phase_matrix["points"]]
    y_points = [float(p[1]) for p in phase_matrix["points"]]
    ax.scatter(x_points, y_points, color="red", s=10, zorder=3)

    fig.colorbar(c, ax=ax)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")

    fig.tight_layout()
    fig.savefig(path_to_save, dpi = 300, bbox_inches = "tight")


# %%
phase_space = get_2D_phase_space(
    combined_root = "../data_private/slurm_test2",
    yaml_criteria = {},
    csv_criteria = {
        "options_enzymes": {"name": {
            "A": {"regions": "[0]"},
            "B": {"regions": "[0]"},
            "C": {"regions": "[0]"},
        }},
        "options_spontaneous_reactions": {"start_species": {
            "Z": {"k": 0.1}
        }}
    },
    x_axis = [
        "options_spontaneous_reactions",
        {"start_species": "X", "end_species": "R1"},
        "k"],
    y_axis = ["options_spontaneous_reactions",
        {"start_species": "Y", "end_species": "R1"},
        "k"],
    file_name_to_plot = "fluxes.json"
)
matrix = build_phase_space_matrix(
    phase_space, "R"
)
plot_2D_phase_space(matrix, "../data_private/slurm_test2/phase_space1.png",
    "k1", "k2", log_x=True, log_y=True)

