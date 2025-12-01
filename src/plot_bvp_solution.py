import sys
import os
import re
import glob
import argparse
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from auxiliary_functions_using_standard_library import (
    format_sci, pickle_load_binary,
    load_json, find_max_in_nested_dict)
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme
from auxiliary_functions_framework_organization import get_species_concentrations_from_json_file
import imageio.v3 as iio


def plot_steady_state_concentrations(
        reaction_network,
        num_regions,
        num_mesh_points_in_regions,
        radii,
        membrane_radii,
        output_file_name, species_concentrations_to_plot, title = None, ymax = None):
    """ Plots the concentrations at a specific time point.
    output_file_name: should have .png or any other extension
    species_concentrations_to_plot: dictionary with keys region, n, species,
    title: string with title to put. Default to no title written.
    ymax: maximum value on y axis. If None (default), automatically computed from data.
    """
    x_values = []
    y_values = {}
    for species_idx, species in enumerate(reaction_network.species):
        species_y_values = []
        for region in range(num_regions):
            for n in range(num_mesh_points_in_regions[region]):
                if species_idx == 0:
                    x_values.append(radii[region][n])
                species_y_values.append(species_concentrations_to_plot[region][n][species])
        y_values[species] = species_y_values

    fig, ax = plt.subplots(1,1, figsize = (5,3))
    for x_value in x_values:
        ax.axvline(x_value/max(x_values), ymin = 0.95, ymax = 1, color="k")
    for species in reaction_network.species:
        curve, = ax.plot(np.array(x_values)/max(x_values), y_values[species], label=species.name)
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
    for membrane_radius in membrane_radii:
        ax.axvline(membrane_radius/max(membrane_radii), linestyle = "--", alpha = 0.5, c = "k")

    max_value = max(max(y_values[species]) for species in reaction_network.species)
    if ymax == None:
        ymax = max_value * 1.05
    if title != None:
        ax.set_title(title, loc="left")
    ax.set_ylim(ymin=0, ymax = ymax)
    ax.set_xlim(xmin=0, xmax = 1.1)
    fig.savefig(output_file_name, dpi = 300, bbox_inches='tight')
    plt.close(fig)

def make_newton_iterations_gif(
        reaction_network,
        num_regions,
        num_mesh_points_in_regions,
        radii,
        membrane_radii,
        iteration_data_folder, gif_output_folder, species_lookup_dict):
    """
    Important: delete any .json files from previous simulations that may not be overwritten.
    (Should already automatically have been done by snakemake through cleanup_old_iterations rule)
    """
    file_to_create = os.path.join(gif_output_folder, "newton_iterations.gif")
    # Get concentration files from which to make the gif
    files = glob.glob(os.path.join(iteration_data_folder, ".iteration_nr_*_concentrations.json"))
    files.sort()
    if not files:
        print(f"No iterations in {iteration_data_folder} from which to create gif")
        return   
    # Find out number of leading zeros within folder
    match = re.search(r'\.iteration_nr_(\d+)_concentrations\.json', os.path.basename(files[0]))
    if match:
        digit_count = len(match.group(1))
        print(f"The filenames have {digit_count} digits.")
    else:
        raise ValueError(f"Could not find a number in the filename {os.path.basename(files[0])}.")

    # Get the maximum concentration within the files (internal or external concentration)
    max_concentration_value = 0
    for file in files:
        concentration_dict = load_json(file)
        max_value = find_max_in_nested_dict(concentration_dict)
        max_concentration_value = max(max_concentration_value, max_value)
    for species in reaction_network.species:
        max_concentration_value = max(max_concentration_value, species.external_concentration)
    max_y = max_concentration_value * 1.1 # make space for some vertical padding
    
    # Read maximum y value previously used to create a gif
    # If the new highest concentration is higher, previous png files have to be
    # removed to then be rewritten with the higher y max
    max_y_filename = os.path.join(gif_output_folder, ".gif_ymax")
    if os.path.isfile(max_y_filename):
        with open(max_y_filename, "r") as f:
            previous_max_y = float(f.read().strip())
        if previous_max_y > max_y:
            max_y = previous_max_y
        else:
            for file in files:
                png_file = png_file = os.path.splitext(file)[0] + ".png"
                os.remove(png_file)

    # Create concentration files
    png_files_created = []
    print("creating files for gif", file_to_create)
    for file in tqdm(files, file=sys.stderr):
        png_file = os.path.splitext(file)[0] + ".png" # remove .json
        if os.path.isfile(png_file):
            continue
        number = int(re.findall(r"\d+", os.path.basename(png_file))[0])+1
        number_with_max_digits = f"{number:0{digit_count}d}"
        title_lines = [f"iteration #{number_with_max_digits}"]
        species_concentrations_to_plot_dict_with_strings = load_json(file)
        residual_norm_file= file.replace("_concentrations.json", "_F_vector_norm.json")
        if os.path.isfile(residual_norm_file):
            residual_norm = load_json(residual_norm_file)
            title_lines.append(f"residual norm: {format_sci(residual_norm)}")
        species_concentrations_to_plot_dict = get_species_concentrations_from_json_file(
            species_concentrations_to_plot_dict_with_strings, species_lookup_dict)
        plot_steady_state_concentrations(
            reaction_network=reaction_network,
            num_regions=num_regions,
            num_mesh_points_in_regions=num_mesh_points_in_regions,
            radii=radii,
            membrane_radii=membrane_radii,
            output_file_name=png_file,
            species_concentrations_to_plot=species_concentrations_to_plot_dict,
            title = "\n".join(title_lines),
            ymax = max_y)
        png_files_created.append(png_file)
    # Put all the pngs together
    print("creating gif", file_to_create)
    with iio.imopen(file_to_create, "w") as writer:
        for filename in tqdm(png_files_created, file=sys.stderr):
            image = iio.imread(filename)
            writer.write(image)
    print("gif created!")
    # Create file with maximum concentration (to know whether previous png files
    #can be reused)
    with open(max_y_filename, "w") as f:
        f.write(str(max_y))

if __name__ == "__main__":
    # Parse arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("folder_to_solve", type=str, help="Path to folder with system info")
    parser.add_argument("--plot_iteration", type=str, default=None,
                        help="Optional. Iteration number of which to plot the concentration.")
    args = parser.parse_args()

    # Load all the passed information
    FOLDER_TO_SOLVE = args.folder_to_solve
    PLOT_ITERATION = args.plot_iteration

    # Load inputs and define global parameters
    REACTION_NETWORK = pickle_load_binary(os.path.join(FOLDER_TO_SOLVE, ".pickled_reaction_network"))
    SYSTEM_GEOMETRY_DICT = load_json(os.path.join(FOLDER_TO_SOLVE, ".expanded_system_geometry"))
    SYSTEM_MESH_DICT= load_json(os.path.join(FOLDER_TO_SOLVE, ".expanded_system_mesh"))

    SPECIES_LOOKUP = {sp.name: sp for sp in REACTION_NETWORK.species}

    R = SYSTEM_GEOMETRY_DICT["geometry_config"]["outer_membrane_radius"]
    NUM_MESH_POINTS_IN_REGIONS = SYSTEM_GEOMETRY_DICT["geometry_config"]["num_mesh_points_in_regions"]
    NUM_REGIONS = SYSTEM_GEOMETRY_DICT["geometry_config"]["num_regions"]
    MEMBRANE_RADII = SYSTEM_GEOMETRY_DICT["geometry_config"]["membrane_radii"]

    RADII = SYSTEM_MESH_DICT["radii"]

    SOLVER_INPUT = load_json(os.path.join(FOLDER_TO_SOLVE, "solver_input.json"))
    if SOLVER_INPUT["output_options"]["create_gif_with_saved_data"] is True and SOLVER_INPUT["variables_to_save"]["save_concentrations"] is False:
        raise ValueError("Cannot make the gif if the concentrations are not saved.")
    
    ITERATIONS_FOLDER = os.path.join(FOLDER_TO_SOLVE, "solver_iteration_data")
    if PLOT_ITERATION is None:
        make_newton_iterations_gif(
            reaction_network=REACTION_NETWORK,
            num_regions=NUM_REGIONS,
            num_mesh_points_in_regions=NUM_MESH_POINTS_IN_REGIONS,
            radii=RADII,
            membrane_radii=MEMBRANE_RADII,
            iteration_data_folder=ITERATIONS_FOLDER,
            gif_output_folder=FOLDER_TO_SOLVE,
            species_lookup_dict=SPECIES_LOOKUP)
    else:
        # Find out concentrations file from the iteration number
        pattern = os.path.join(ITERATIONS_FOLDER, ".iteration_nr_*_concentrations.json")
        files = glob.glob(pattern)
        regex = re.compile(rf"\.iteration_nr_*{PLOT_ITERATION}_concentrations\.json$")
        matches = [f for f in files if regex.search(os.path.basename(f))]
        match_file= matches[0]
        species_concentrations_to_plot_with_strings = load_json(match_file)
        species_concentrations_to_plot_dict = get_species_concentrations_from_json_file(
            species_concentrations_to_plot_with_strings, SPECIES_LOOKUP)
        file_to_create = os.path.splitext(match_file)[0] + "_individual.png" # remove .json
        plot_steady_state_concentrations(
            reaction_network=REACTION_NETWORK,
            num_regions=NUM_REGIONS,
            num_mesh_points_in_regions=NUM_MESH_POINTS_IN_REGIONS,
            radii=RADII,
            membrane_radii=MEMBRANE_RADII,
            output_file_name=file_to_create,
            species_concentrations_to_plot=species_concentrations_to_plot_dict)