import sys
import os
import re
import glob
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from auxiliary_functions import read_yaml_file
from auxiliary_functions_using_standard_library import (
    format_sci, pickle_load_binary,
    load_json, find_max_in_nested_dict)
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme
from auxiliary_functions_framework_organization import get_species_concentrations_from_json_file
import imageio.v3 as iio

#plt.rcParams['text.usetex'] = True

def plot_steady_state_concentrations(
        reaction_network,
        num_regions,
        num_mesh_points_in_regions,
        radii,
        membrane_radii,
        output_file_name,
        species_concentrations_to_plot,
        system_geometry_dict,
        title = None,
        ymax = None
    ):
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

    fig, ax = add_theory_curve_to_ax(
                fig,
                ax,
                reaction_network,
                num_regions,
                membrane_radii,
                system_geometry_dict
            )
    
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
    fig.tight_layout()
    if output_file_name != None:
        fig.savefig(output_file_name, dpi = 300, bbox_inches='tight')
    return fig, ax

def make_newton_iterations_gif(
        reaction_network,
        num_regions,
        num_mesh_points_in_regions,
        radii,
        membrane_radii,
        iteration_data_folder,
        gif_output_folder,
        species_lookup_dict,
        system_geometry_dict
    ):
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
        species_concentrations_to_plot_dict_with_strings = load_json(file)
        species_concentrations_to_plot_dict = get_species_concentrations_from_json_file(
            species_concentrations_to_plot_dict_with_strings, species_lookup_dict)
        fig, ax = plot_steady_state_concentrations(
            reaction_network=reaction_network,
            num_regions=num_regions,
            num_mesh_points_in_regions=num_mesh_points_in_regions,
            radii=radii,
            membrane_radii=membrane_radii,
            output_file_name=None,
            species_concentrations_to_plot=species_concentrations_to_plot_dict,
            system_geometry_dict=system_geometry_dict,
            title = f"iteration #{number_with_max_digits}",
            ymax = max_y)
        fig.savefig(png_file, dpi = 300)
        plt.close(fig)
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

def plot_convergence_progress(
        folder_to_solve,
        reaction_network
):
    dataframe = pd.read_csv(os.path.join(folder_to_solve, ".convergence_logger.csv"))
    fig, ax = plt.subplots(2,2, figsize=(6,6))
    ax = np.ndarray.flatten(ax)
    for species in reaction_network.species:
        ax[0].plot(dataframe["iteration"], dataframe[species.name], label=species.name)
    ax[1].plot(dataframe["iteration"], dataframe["max_relative_change"])
    ax[2].plot(dataframe["iteration"], dataframe["max_Delta_u"])
    ax[3].plot(dataframe["iteration"], dataframe["F_vector_norm"])
    ax[2].set_xlabel("iteration")
    ax[3].set_xlabel("iteration")
    ax[0].set_ylabel("relative difference between \n reaction and boundary flux")
    #ax[0].set_ylabel(
    #    r"Relative difference $\frac{\lvert \Phi_{\text{react}} - \Phi_{\text{bound}} \rvert}"
    #    r"{\max\!\left(\lvert \Phi_{\text{react}} \rvert, \lvert \Phi_{\text{bound}} \rvert\right)}$"
    #)
    ax[1].set_ylabel("max_relative_change")
    ax[2].set_ylabel("max_Delta_u")
    ax[3].set_ylabel("F_vector_norm")
    ax[2].set_yscale('log')
    ax[3].set_yscale('log')
    ax[0].legend()
    fig.tight_layout()
    fig.savefig(os.path.join(folder_to_solve, "convergence.png"), bbox_inches = "tight", dpi=300)
    plt.close()

def add_theory_curve_to_ax(
        fig,
        ax,
        reaction_network,
        num_regions,
        membrane_radii,
        system_geometry_dict
    ):
    # Case of one spontaneous reaction with no inner boundaries
    if (len(reaction_network.spontaneous_reactions) == 1
        and len(reaction_network.enzymatic_reactions) == 0
        and num_regions == 1
    ):  
        s_reaction = reaction_network.spontaneous_reactions[0]
        s = s_reaction.start_species
        s_lambda = np.sqrt(s_reaction.k / s.diffusion_constant)
        external_radius = membrane_radii[-1]
        A = - ((s.permeability_constant / s.diffusion_constant * s.external_concentration * external_radius**2)
            / (
                np.exp(s_lambda * external_radius)*(s_lambda * external_radius - 1 + (s.permeability_constant * external_radius)/s.diffusion_constant)
               + np.exp(-s_lambda * external_radius) * (s_lambda * external_radius + 1 - (s.permeability_constant * external_radius)/s.diffusion_constant) 
              )
        )
        c = lambda r: 1/r * A *(np.exp(-s_lambda*r)-np.exp(s_lambda*r))
        r_to_plot = np.linspace(external_radius*0.01, external_radius, num = 100)
        ax.plot([r/external_radius for r in r_to_plot], [c(r) for r in r_to_plot], linestyle= "--",
                label = f"theory for {s.name}", zorder = -1, 
                linewidth = 1,
                color = "k",
                alpha = 0.5
        )
        ax.legend()
    
    # Case of one spontaneous reaction with one inner boundary
    elif (len(reaction_network.spontaneous_reactions) == 1
        and len(reaction_network.enzymatic_reactions) == 0
        and num_regions == 2
    ):   
        reaction = reaction_network.spontaneous_reactions[0]
        X = reaction.start_species
        X_lambda = np.sqrt(reaction.k / X.diffusion_constant)
        r_inner = membrane_radii[0]
        external_radius = membrane_radii[-1]
        s = np.sinh(X_lambda * r_inner)
        c = np.cosh(X_lambda * r_inner)
        beta = X.permeability_constant / X.diffusion_constant
        alpha = X.permeability_constant * external_radius / X.diffusion_constant
        
        rho = np.exp(-2 * X_lambda * r_inner) * (
            X.diffusion_constant * (r_inner**2 * X_lambda**2 * c + r_inner * X_lambda * (c - s) - s)
            + X.permeability_constant * r_inner**2 * X_lambda * (s + c)
        ) / (
            X.diffusion_constant * (r_inner**2 * X_lambda**2 * c - r_inner * X_lambda * (s + c) + s)
            + X.permeability_constant * r_inner**2 * X_lambda * (s - c)
        )

        A = (beta * X.external_concentration * external_radius**2)/(
            np.exp(-X_lambda * external_radius)*(alpha - X_lambda * external_radius - 1) + rho * np.exp(X_lambda * external_radius)*(alpha + X_lambda*external_radius -1)
        )
        
        B = rho * A

        S = (
            beta * r_inner * (A * np.exp(-X_lambda * r_inner) + B * np.exp(X_lambda * r_inner))
        ) / (
            X_lambda * r_inner * np.cosh(X_lambda * r_inner)
            - np.sinh(X_lambda * r_inner)
            + beta * r_inner * np.sinh(X_lambda * r_inner)
        )

        mesh_points_in_regions = system_geometry_dict["mesh_points_in_regions"]
        c_1 = lambda r : S * np.sinh(X_lambda * r) / r
        c_2 = lambda r : (A * np.exp(-X_lambda * r) + B * np.exp(X_lambda * r))/r
        c = [c_1, c_2]
        for region_idx in [0,1]:
            region_radii = mesh_points_in_regions[region_idx]
            # region_radii[0] skipped to avoid division by 0
            r_to_plot = np.linspace(region_radii[0]+external_radius*0.01, region_radii[-1], num = 100)
            if region_idx == 0:
                label = f"theory for {X.name}"
            else:
                label = None
            ax.plot([r/external_radius for r in r_to_plot],
                [c[region_idx](r) for r in r_to_plot], linestyle= "--",
                label = label, zorder = -1, 
                linewidth = 1,
                color = "k",
                alpha = 0.5
            )
        ax.legend()
    
    # Case of one spontaneous reaction with two inner boundaries
    elif (len(reaction_network.spontaneous_reactions) == 1
        and len(reaction_network.enzymatic_reactions) == 0
        and num_regions == 3
    ):  
        """This below was written by ChatGPT."""
        reaction = reaction_network.spontaneous_reactions[0]
        X = reaction.start_species

        lam = np.sqrt(reaction.k / X.diffusion_constant)  # lambda
        D = X.diffusion_constant
        p = X.permeability_constant
        beta = p / D

        # Two inner membranes
        R1 = membrane_radii[0]
        R2 = membrane_radii[1]
        R  = membrane_radii[2]

        c_ext = X.external_concentration

        # --- Helpers for evaluating c and c' for the basis solutions ---
        # For region 1: c1(r) = S1*sinh(lam r)/r
        def c1_val(r, S1):
            return S1 * np.sinh(lam * r) / r

        def c1_der(r, S1):
            # d/dr [ S1*sinh(lam r)/r ] = S1*(lam*cosh(lam r)/r - sinh(lam r)/r^2)
            return S1 * (lam * np.cosh(lam * r) / r - np.sinh(lam * r) / (r**2))

        # For region j>=2: cj(r) = (A*e^{-lam r} + B*e^{lam r})/r
        def cAB_val(r, A, B):
            return (A * np.exp(-lam * r) + B * np.exp(lam * r)) / r

        def cAB_der(r, A, B):
            # derivative of (f(r)/r) with f=A e^{-lam r}+B e^{lam r}
            f  = A * np.exp(-lam * r) + B * np.exp(lam * r)
            fp = -lam * A * np.exp(-lam * r) + lam * B * np.exp(lam * r)
            return fp / r - f / (r**2)

        # --- Build the 5x5 linear system M x = b for x=[S1, A2, B2, A3, B3] ---
        # Conditions:
        # (1) c1'(R1) = c2'(R1)
        # (2) c1'(R1) = beta*(c2(R1) - c1(R1))
        # (3) c2'(R2) = c3'(R2)
        # (4) c2'(R2) = beta*(c3(R2) - c2(R2))
        # (5) c3'(R)  = beta*(c_ext - c3(R))

        M = np.zeros((5, 5), dtype=float)
        b = np.zeros(5, dtype=float)

        # Row 0: c1'(R1) - c2'(R1) = 0
        # coefficients for S1, A2, B2, A3, B3
        M[0, 0] = c1_der(R1, 1.0)
        # -c2'(R1) contributions:
        # c2'(R1) is linear in A2,B2 so put negatives
        # We'll get coeffs by evaluating derivative with A=1,B=0 and A=0,B=1
        M[0, 1] = -cAB_der(R1, 1.0, 0.0)
        M[0, 2] = -cAB_der(R1, 0.0, 1.0)
        # A3,B3 not in this equation
        b[0] = 0.0

        # Row 1: c1'(R1) - beta*(c2(R1) - c1(R1)) = 0
        # => c1'(R1) + beta*c1(R1) - beta*c2(R1) = 0
        M[1, 0] = c1_der(R1, 1.0) + beta * c1_val(R1, 1.0)
        M[1, 1] = -beta * cAB_val(R1, 1.0, 0.0)
        M[1, 2] = -beta * cAB_val(R1, 0.0, 1.0)
        b[1] = 0.0

        # Row 2: c2'(R2) - c3'(R2) = 0
        M[2, 1] = cAB_der(R2, 1.0, 0.0)
        M[2, 2] = cAB_der(R2, 0.0, 1.0)
        M[2, 3] = -cAB_der(R2, 1.0, 0.0)
        M[2, 4] = -cAB_der(R2, 0.0, 1.0)
        b[2] = 0.0

        # Row 3: c2'(R2) - beta*(c3(R2) - c2(R2)) = 0
        # => c2'(R2) + beta*c2(R2) - beta*c3(R2) = 0
        M[3, 1] = cAB_der(R2, 1.0, 0.0) + beta * cAB_val(R2, 1.0, 0.0)
        M[3, 2] = cAB_der(R2, 0.0, 1.0) + beta * cAB_val(R2, 0.0, 1.0)
        M[3, 3] = -beta * cAB_val(R2, 1.0, 0.0)
        M[3, 4] = -beta * cAB_val(R2, 0.0, 1.0)
        b[3] = 0.0

        # Row 4: c3'(R) = beta*(c_ext - c3(R))
        # => c3'(R) + beta*c3(R) = beta*c_ext
        M[4, 3] = cAB_der(R, 1.0, 0.0) + beta * cAB_val(R, 1.0, 0.0)
        M[4, 4] = cAB_der(R, 0.0, 1.0) + beta * cAB_val(R, 0.0, 1.0)
        b[4] = beta * c_ext

        # Solve
        S1, A2, B2, A3, B3 = np.linalg.solve(M, b)

        # --- Define concentration functions for each region (compatible with your plotting pattern) ---
        def c_1(r):
            r = np.asarray(r, dtype=float)
            out = np.empty_like(r)
            # safe near 0: limit S1*sinh(lam r)/r -> S1*lam
            small = np.isclose(r, 0.0)
            out[small] = S1 * lam
            out[~small] = S1 * np.sinh(lam * r[~small]) / r[~small]
            return out

        def c_2(r):
            r = np.asarray(r, dtype=float)
            return (A2 * np.exp(-lam * r) + B2 * np.exp(lam * r)) / r

        def c_3(r):
            r = np.asarray(r, dtype=float)
            return (A3 * np.exp(-lam * r) + B3 * np.exp(lam * r)) / r

        c = [c_1, c_2, c_3]

        # --- Plot using your mesh_points_in_regions structure (now 3 regions) ---
        mesh_points_in_regions = system_geometry_dict["mesh_points_in_regions"]

        for region_idx in [0, 1, 2]:
            region_radii = mesh_points_in_regions[region_idx]
            # skip first point in region 0 to avoid r=0 if present
            start_idx = 1 if (region_idx == 0 and np.isclose(region_radii[0], 0.0)) else 0

            r_to_plot = np.linspace(region_radii[start_idx], region_radii[-1], num=100)

            label = f"theory for {X.name}" if region_idx == 0 else None
            ax.plot([rr / R for rr in r_to_plot],
                    [float(c[region_idx](rr)) for rr in r_to_plot],
                    linestyle="--",
                    label=label,
                    zorder=-1,
                    linewidth=1,
                    color = "k",
                    alpha = 0.5
            )
    
    
    
    return fig, ax

       













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
    SYSTEM_GEOMETRY_DICT = load_json(os.path.join(FOLDER_TO_SOLVE, ".expanded_system_geometry.json"))
    SYSTEM_MESH_DICT= load_json(os.path.join(FOLDER_TO_SOLVE, ".expanded_system_mesh.json"))

    SPECIES_LOOKUP = {sp.name: sp for sp in REACTION_NETWORK.species}

    R = SYSTEM_GEOMETRY_DICT["geometry_config"]["outer_membrane_radius"]
    NUM_MESH_POINTS_IN_REGIONS = SYSTEM_GEOMETRY_DICT["geometry_config"]["num_mesh_points_in_regions"]
    NUM_REGIONS = SYSTEM_GEOMETRY_DICT["geometry_config"]["num_regions"]
    MEMBRANE_RADII = SYSTEM_GEOMETRY_DICT["geometry_config"]["membrane_radii"]

    RADII = SYSTEM_MESH_DICT["radii"]

    SOLVER_INPUT = read_yaml_file(os.path.join(FOLDER_TO_SOLVE, "parameters_solver_input.yaml"))
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
            species_lookup_dict=SPECIES_LOOKUP,
            )
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
            species_concentrations_to_plot=species_concentrations_to_plot_dict,
            system_geometry_dict=SYSTEM_GEOMETRY_DICT["geometry_config"]
        )
    
    if SOLVER_INPUT["output_options"]["log_convergence_progress"]:
        plot_convergence_progress(FOLDER_TO_SOLVE, REACTION_NETWORK)
