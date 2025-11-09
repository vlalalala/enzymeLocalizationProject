import sys
import os
import re
from contextlib import redirect_stdout
import argparse
from tqdm import tqdm
import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_matrix
import matplotlib.pyplot as plt
from auxiliary_functions_using_standard_library import (rename_iteration_files,
    nested_max, all_non_negative, format_sci, pickle_load_binary,
    dump_json, find_sorted_unique_files_with_max_digits_and_max_value, load_json,
    find_max_in_nested_dict)
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme
from auxiliary_functions import save_matrix_as_sparse_txt
from auxiliary_functions_framework_organization import get_species_concentrations_from_json_file

import imageio.v3 as iio


def plot_steady_state_concentrations(output_file_name, species_concentrations_to_plot, title = None, ymax = None):
    x_values = []
    y_values = {}
    for species_idx, species in enumerate(REACTION_NETWORK.species):
        species_y_values = []
        for region in range(NUM_REGIONS):
            for n in range(NUM_MESH_POINTS_IN_REGIONS[region]):
                if species_idx == 0:
                    x_values.append(RADII[region][n])
                species_y_values.append(species_concentrations_to_plot[region][n][species])
        y_values[species] = species_y_values

    fig, ax = plt.subplots(1,1, figsize = (5,3))
    for x_value in x_values:
        ax.axvline(x_value/max(x_values), ymin = 0.95, ymax = 1, color="k")
    for species in REACTION_NETWORK.species:
        curve, = ax.plot(x_values/max(x_values), y_values[species], label=species.name)
        color = curve.get_color()
        ax.hlines(species.external_concentration, xmin=1, xmax = 1.1, color = color)
    ax.set_ylabel("concentration / M")
    ax.set_xlabel("relative distance to origin / r/R")
    ax.legend(
        loc='upper center',      # anchor point of legend
        bbox_to_anchor=(0.5, -0.25),  # (x, y) position in figure coordinates
        ncol=3,                  # number of columns
        frameon=False
    )
    for membrane_radius in MEMBRANE_RADII:
        ax.axvline(membrane_radius/max(MEMBRANE_RADII), linestyle = "--", alpha = 0.5, c = "k")

    max_value = max(max(y_values[species]) for species in REACTION_NETWORK.species)
    if ymax == None:
        ymax = max_value * 1.05
    if title != None:
        ax.set_title(title, loc="left")
    ax.set_ylim(ymin=0, ymax = ymax)
    ax.set_xlim(xmin=0, xmax = 1.1)
    fig.savefig(output_file_name, dpi = 300, bbox_inches='tight')
    plt.close(fig)

def make_newton_iterations_gif(iteration_data_folder, gif_output_folder):
    """
    Important: delete any .json files from previous simulations that may not be overwritten.
    """
    file_to_create = os.path.join(gif_output_folder, "newton_iterations.gif")
    sorted_files, max_digits = find_sorted_unique_files_with_max_digits_and_max_value(iteration_data_folder, ".iteration_nr_*_concentrations.json", max_iteration_value = MAX_NUM_NEWTON_ITERATIONS)
    max_concentration_value = 0
    for file in sorted_files:
        concentration_dict = load_json(file)
        max_value = find_max_in_nested_dict(concentration_dict)
        if max_value > max_concentration_value:
            max_concentration_value = max_value
    for species in REACTION_NETWORK.species:
        if species.external_concentration > max_concentration_value:
            max_concentration_value = species.external_concentration
    max_concentration_value *= 1.1
    png_files_created = []
    print("creating files for gif", file_to_create)
    for file in tqdm(sorted_files, file=sys.stderr):
        png_file = os.path.splitext(file)[0] + ".png" # remove .json
        number = int(re.findall(r"\d+", os.path.basename(png_file))[0])+1
        number_with_max_digits = f"{number:0{max_digits}d}"
        title_lines = [f"iteration #{number_with_max_digits}"]
        species_concentrations_to_plot_dict_with_strings = load_json(file)
        if SOLVER_PARAMS["VARIABLES_TO_SAVE"]["save_F_vector_norm"]:
            residual_norm_file= file.replace("_concentrations.json", "_F_vector_norm.json")
            residual_norm = load_json(residual_norm_file)
            title_lines.append(f"residual norm: {format_sci(residual_norm)}")
        species_concentrations_to_plot_dict = get_species_concentrations_from_json_file(
            species_concentrations_to_plot_dict_with_strings)
        plot_steady_state_concentrations(
            png_file, species_concentrations_to_plot_dict,
            title = "\n".join(title_lines),
            ymax = max_concentration_value)
        png_files_created.append(png_file)
    
    print("creating gif", file_to_create)
    with iio.imopen(file_to_create, "w") as writer:
        for filename in tqdm(png_files_created, file=sys.stderr):
            image = iio.imread(filename)
            writer.write(image)

if __name__ == "__main__":
    # Parse arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("folder_to_solve", type=str, help="Path to folder with system info")
    # use int(float()) to be able to pass scientific notation
    parser.add_argument("--max-iterations", type=lambda x: int(float(x)), help="Maximum Newton iterations") 
    parser.add_argument("--previous-solution", type=str, default=None,
                        help="Optional path to previous iteration solution file.")
    args = parser.parse_args()

    # Load all the passed information
    FOLDER_TO_SOLVE = args.folder_to_solve
    MAX_NUM_NEWTON_ITERATIONS = args.max_iterations
    PREVIOUS_SOLUTION = args.previous_solution

    # Load inputs and define global parameters
    REACTION_NETWORK = pickle_load_binary(os.path.join(FOLDER_TO_SOLVE, ".REACTION_NETWORK_pickle"))
    # Lookup to be able to match the species from the species name saved in .json files  
    SPECIES_LOOKUP = {sp.name: sp for sp in REACTION_NETWORK.species}

    # Step 1: Define all geometry variables
    R = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["outer_membrane_radius"]
    MESH_POINTS_IN_REGIONS = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["MESH_POINTS_IN_REGIONS"]
    NUM_MESH_POINTS_IN_REGIONS = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["NUM_MESH_POINTS_IN_REGIONS"]
    NUM_REGIONS = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["NUM_REGIONS"]
    MEMBRANE_RADII = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["MEMBRANE_RADII"]

    # Step 0: Get all solver parameters
    SOLVER_PARAMS = pickle_load_binary(os.path.join(FOLDER_TO_SOLVE, ".solver_info_pickle"))
    if SOLVER_PARAMS["OUTPUT_OPTIONS"]["create_gif_with_saved_data"] is True and SOLVER_PARAMS["VARIABLES_TO_SAVE"]["save_concentrations"] is False:
        raise ValueError("Cannot make the gif if the concentrations are not saved.")