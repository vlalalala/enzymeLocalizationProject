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
from create_system_mesh import build_system_mesh

def make_newton_iterations_gif(
        folder_to_solve,
        creeping_reaction_simulation_folder,
    ):
    """
    Important: delete any .json files from previous simulations that may not be overwritten.
    (Should already automatically have been done by snakemake through cleanup_old_iterations rule)
    """
    file_to_create = os.path.join(creeping_reaction_simulation_folder, "newton_iterations.gif")
    iteration_data_folder = os.path.join(creeping_reaction_simulation_folder, "solver_iteration_data")
    # Get concentration files from which to make the gif
    progress_log = os.path.join(creeping_reaction_simulation_folder, ".progress_log.csv")
    reaction_network = pickle_load_binary(os.path.join(creeping_reaction_simulation_folder, ".pickled_reaction_network"))
    species_lookup = {sp.name: sp for sp in reaction_network.species}
    # First step: figure out the maximum y-value for all interpolation iterations
    max_y = 0
    # Extract the interpolation number from the log filename
    # Find all matching json files for this interpolation number
    json_files = list(Path(iteration_data_folder).glob(
        f"interpolation_iteration_nr_0_Newton_iteration_nr_*_concentrations.json"
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
    json_files = list(Path(iteration_data_folder).glob(
        f"interpolation_iteration_nr_0_Newton_iteration_nr_*_concentrations.json"
    ))
    json_files.sort()

    system_geometry = load_json(
        os.path.join(
            folder_to_solve, f".system_geometry.json"
        )
    )
    expanded_system_geometry_dict, expanded_system_mesh_dict = build_system_mesh(
                system_geometry, 
                reaction_network,
                0
            )
    png_files_created = []
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
                species_concentrations_to_plot_dict_with_strings, species_lookup)
            fig, _ = plot_steady_state_concentrations(
                reaction_network=reaction_network,
                num_regions=system_geometry["geometry_config"]["num_regions"],
                num_mesh_points_in_regions=system_geometry["geometry_config"]["num_mesh_points_in_regions"],
                radii=expanded_system_mesh_dict["radii"],
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
    try:
        with imageio.get_writer(file_to_create, mode="I", loop=0, duration=0.1) as writer:
            for filename in tqdm(png_files_created, file=sys.stderr):
                writer.append_data(imageio.imread(filename))
    except:
        print(f"Could not create {file_to_create}")

    print(f"gif created! {file_to_create}")
    # Create file with maximum concentration (to know whether previous png files
    #can be reused)
    with open(os.path.join(creeping_reaction_simulation_folder, "max_y"), "w") as f:
        f.write(str(max_y))