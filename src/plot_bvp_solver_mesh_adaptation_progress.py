from pathlib import Path
import re
import os
import matplotlib.pyplot as plt
import sys
from auxiliary_functions_using_standard_library import pickle_load_binary, find_max_in_nested_dict
import pandas as pd
import numpy as np
from auxiliary_functions_using_standard_library import (pickle_load_binary,
    load_json)
from tqdm import tqdm
from run_bvp_solver import plot_steady_state_concentrations
import imageio
from auxiliary_functions_framework_organization import get_dict_with_correct_key_types_from_json_file
import matplotlib.pyplot as plt

def plot_convergence_progress(
    folder_to_solve,
    reaction_network
):  
    files = list(Path(folder_to_solve).glob(".progress_log_interpolating_*_times.csv"))
    files.sort()
    fig, ax = plt.subplots(4, 1, figsize = (3,9))
    colors = {species.name: color for species, color in zip(reaction_network.species, plt.cm.tab10.colors)}
    
    iterations_offset = 0
    for file_idx, file in enumerate(files):
        df = pd.read_csv(file)
        x = df["iteration"] + iterations_offset
        ax[0].plot(x, df["F_vector_norm"])
        ax[1].plot(x, df["tau"])
        for species in reaction_network.species:
            if file_idx == 0:
                label = species.name
            else:
                label = None
            ax[2].plot(x, df[f"{species.name}_absolute"], label=label, color=colors[species.name])
            ax[3].plot(x, df[f"{species.name}_relative"], label=label, color=colors[species.name])

        iterations_offset = x.iloc[-1] # next file starts after last iteration
        for i in range(4):
            ax[i].axvline(iterations_offset, ls = ":", color = "k", alpha = 0.6)

    # Find y value at x=10 across all data in ax[2]
    for line in ax[2].get_lines():
        x_data = line.get_xdata()
        y_data = line.get_ydata()
        idx = np.searchsorted(x_data, 100)
        if idx < len(y_data):
            y_at_10 = y_data[idx]
            break

    #ax[2].set_ylim(top=y_at_10)

    ax[0].set_ylabel("residual vector norm")
    ax[0].set_xlabel("iteration")
    ax[0].set_yscale('log')

    ax[1].set_ylabel("step size")
    ax[1].set_xlabel("iteration")
    ax[1].set_yscale('log')

    
    ax[2].set_ylabel("absolute difference between \n reaction and boundary flux")
    ax[2].set_xlabel("iteration")
    ax[2].legend()
    ax[2].set_yscale('log')

    ax[3].set_ylabel("relative difference between \n reaction and boundary flux")
    ax[3].set_xlabel("iteration")
    ax[3].legend()
    ax[3].set_yscale('log')
    
    fig.tight_layout()
    fig.savefig(os.path.join(folder_to_solve, "convergence.png"), bbox_inches = "tight", dpi=300)
    plt.close()

def make_newton_iterations_gif(
        folder_to_plot,
        reaction_network,
        species_lookup_dict,
    ):
    """
    Important: delete any .json files from previous simulations that may not be overwritten.
    (Should already automatically have been done by snakemake through cleanup_old_iterations rule)
    """
    file_to_create = os.path.join(folder_to_plot, "newton_iterations.gif")
    iteration_data_folder = os.path.join(folder_to_plot, "solver_iteration_data")
    # Get concentration files from which to make the gif
    progress_logs = sorted(Path(folder_to_plot).glob(".progress_log_interpolating_*_times.csv"))
    
    # First step: figure out the maximum y-value for all interpolation iterations
    max_y = 0
    for log_file in progress_logs:
        # Extract the interpolation number from the log filename
        match = re.search(r"\.progress_log_interpolating_(\d+)_times\.csv", log_file.name)
        if match:
            interp_nr = match.group(1)
            # Find all matching json files for this interpolation number
            json_files = list(Path(iteration_data_folder).glob(
                f"interpolation_iteration_nr_{interp_nr}_Newton_iteration_nr_*_concentrations.json"
            ))
            json_files.sort()
            for file in json_files:
                concentration_dict = load_json(file)
                max_value = find_max_in_nested_dict(concentration_dict)
                max_y = max(max_y, max_value)
            for species in reaction_network.species:
                max_y = max(max_y, species.external_concentration)

    max_y = max_y * 1.1 # make space for some vertical padding
    
    # Create concentration files
    png_files_created = []
    for log_file in progress_logs:
        # Extract the interpolation number from the log filename
        match = re.search(r"\.progress_log_interpolating_(\d+)_times\.csv", log_file.name)
        if match:
            interp_nr = match.group(1)
            # Find all matching json files for this interpolation number
            json_files = list(Path(iteration_data_folder).glob(
                f"interpolation_iteration_nr_{interp_nr}_Newton_iteration_nr_*_concentrations.json"
            ))
            json_files.sort()
            # add last one
            json_files.append(Path(os.path.join(iteration_data_folder, 
                f"interpolation_iteration_nr_{interp_nr}_final_concentrations.json")))
            system_geometry = load_json(
                os.path.join(
                    iteration_data_folder, f".system_geometry_interpolating_{interp_nr}_times.json"
                )
            )
            system_mesh = load_json(
                os.path.join(
                    iteration_data_folder, f".expanded_system_mesh_interpolating_{interp_nr}_times.json"
                )
            )

            for file in tqdm(json_files, file=sys.stderr):
                png_file = os.path.splitext(file)[0] + ".png" # remove .json and add .png
                if not os.path.isfile(png_file):
                    # create it
                    basename = os.path.basename(png_file)
                    matches = re.findall(r"\d+", basename)
                    interp_nr = int(matches[0])
                    if "final" in basename:
                        newton_nr = "final"
                    elif len(matches) > 1:
                        newton_nr = int(matches[1]) + 1
                    else:
                        newton_nr = "limit"
                    print(file, f"interpolation round #{interp_nr} iteration #{newton_nr}")
                    species_concentrations_to_plot_dict_with_strings = load_json(file)
                    species_concentrations_to_plot_dict = get_dict_with_correct_key_types_from_json_file(
                        species_concentrations_to_plot_dict_with_strings, species_lookup_dict)
                    fig, _ = plot_steady_state_concentrations(
                        reaction_network=reaction_network,
                        num_regions=system_geometry["geometry_config"]["num_regions"],
                        num_mesh_points_in_regions=system_geometry["geometry_config"]["num_mesh_points_in_regions"],
                        radii=system_mesh["radii"],
                        membrane_radii=system_geometry["geometry_config"]["membrane_radii"],
                        output_file_name=None,
                        species_concentrations_to_plot=species_concentrations_to_plot_dict,
                        system_geometry_dict=system_geometry,
                        title=f"interpolation round #{interp_nr} iteration #{newton_nr}",
                        ymax=max_y)

                    fig.savefig(png_file, dpi=300)
                    plt.close(fig)

                png_files_created.append(png_file)
    # Put all the pngs together
    print("creating gif", file_to_create)
    with imageio.get_writer(file_to_create, mode="I", loop=0, duration=0.1) as writer:
        for filename in tqdm(png_files_created, file=sys.stderr):
            writer.append_data(imageio.imread(filename))

    print("gif created!")
    # Create file with maximum concentration (to know whether previous png files
    #can be reused)
    with open(os.path.join(folder_to_plot, "max_y"), "w") as f:
        f.write(str(max_y))

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    REACTION_NETWORK = pickle_load_binary(os.path.join(FOLDER_TO_SOLVE, ".pickled_reaction_network"))
    plot_convergence_progress(FOLDER_TO_SOLVE, REACTION_NETWORK)
    make_newton_iterations_gif(
        FOLDER_TO_SOLVE, REACTION_NETWORK, {sp.name: sp for sp in REACTION_NETWORK.species}
    )
