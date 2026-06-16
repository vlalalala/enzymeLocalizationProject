import sys
import os
import re
import glob
import argparse
import pathlib
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src'))
from auxiliary_functions_using_standard_library import load_json
from solve_linear_reactions_analytically import SystemParams, solve, evaluate_solution
#from wrong_solve_linear_reactions_analytically import SystemParams, solve, evaluate_solution
#from reaction_diffusion_solver import SystemParams, solve, evaluate_solution
from plot_bvp_solution import plot_steady_state_concentrations
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from auxiliary_functions import read_yaml_file
from auxiliary_functions_using_standard_library import (
    format_sci, pickle_load_binary,
    load_json, find_max_in_nested_dict)
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme
from auxiliary_functions_framework_organization import get_dict_with_correct_key_types_from_json_file
from create_system_mesh import build_system_mesh
from pathlib import Path

def construct_analytical_solution(folder):
    return



def plot_numerics_and_analytical_solution(folder):
    """
    For now, the numerical iteration plotted is the 0th
    """
    file_to_create = os.path.join(folder, "plot.png")
    # load numerical solution
    reaction_network = pickle_load_binary(os.path.join(folder, ".pickled_reaction_network"))
    species_lookup = {sp.name: sp for sp in reaction_network.species}
    system_geometry_dict, expanded_system_mesh_dict = build_system_mesh(
                system_geometry_dict = load_json(os.path.join(folder, ".system_geometry.json")), 
                reaction_network = reaction_network,
                mesh_points_duplication_times = 0 # iteration 0, not final one
            )
    species_concentrations_to_plot_dict_with_strings = load_json(Path(folder) / "solver_iteration_data" / "interpolation_iteration_nr_0_final_concentrations.json")
    species_concentrations_to_plot_dict = get_dict_with_correct_key_types_from_json_file(
        species_concentrations_to_plot_dict_with_strings, species_lookup)

    fig, ax = plot_steady_state_concentrations(
        reaction_network = reaction_network,
        num_regions = system_geometry_dict["geometry_config"]["num_regions"],
        num_mesh_points_in_regions = system_geometry_dict["geometry_config"]["num_mesh_points_in_regions"],
        radii = expanded_system_mesh_dict["radii"],
        membrane_radii = system_geometry_dict["geometry_config"]["membrane_radii"],
        output_file_name= file_to_create ,
        species_concentrations_to_plot= species_concentrations_to_plot_dict,
        system_geometry_dict = system_geometry_dict["geometry_config"],
        title = None,
        ymax = None#1e-1,
    )

    # construct analytical solution
    geometry = read_yaml_file(os.path.join(folder, "parameters_geometry.yaml"))
    external_radius = geometry["geometry_config"]["outer_membrane_radius"]

    species_df = pd.read_csv(os.path.join(folder, "species.csv"))
    species_names = list(species_df["name"].unique())
    species_idx = {
        s: i for i, s in enumerate(species_names)
    }
    diffusions = [species_df.loc[(species_df["name"] == species_name),"diffusion_constant"].item()
                  for species_name in species_names]
    external_concentrations = [species_df.loc[(species_df["name"] == species_name),"external_concentration"].item()
                  for species_name in species_names]
    permeabilities = [species_df.loc[(species_df["name"] == species_name),"permeability_constant"].item()
                  for species_name in species_names]
    
    
    ##### Create file with enzyme concentrations
    enzymatic_reactions_df = pd.read_csv(os.path.join(folder, "enzymatic_reactions.csv"))
    enzyme_concentrations = load_json(os.path.join(folder, "enzyme_concentrations.json"))
    spontaneous_reactions_df = pd.read_csv(os.path.join(folder, "spontaneous_reactions.csv"))

    ########################################################################################
    # VERY IMPORTANT: THE PROGRAM ASSUMES THAT THE ORDER IN WHICH THE SPECIES ARE 
    # IN SPECIES.CSV IS THE ORDER OF THE CHEMICAL SPECIES WITHIN THE REACTION CHAIN;
    # THIS WILL NOT WORK IF SPECIES.CSV IS OUTSIDE OF ITS ORDER!
    # the number of reactions is assumed to be 1 less than the number of different species
    ########################################################################################
    k = np.zeros((system_geometry_dict["geometry_config"]["num_regions"], len(species_names), len(species_names)))
    
    for _, row in spontaneous_reactions_df.iterrows():
        i = species_idx[row["start_species"]]
        j = species_idx[row["end_species"]]
        rate = float(row["k"])
        k[:, i, i] -= rate
        k[:, j, i] += rate
    
    for region_idx in range(len(k)):
        for _, row in enzymatic_reactions_df.iterrows():
            i = species_idx[row["start_species"]]
            j = species_idx[row["end_species"]]
            kcat = float(row["k_cat"])
            KM   = float(row["k_M"])
            rate = kcat / KM * enzyme_concentrations[row["enzyme"]][region_idx]
            k[region_idx, i, i] -= rate
            k[region_idx, j, i] += rate
    print(k)
    ##### Compute analytical solution
    params = SystemParams(
        radii = np.array(geometry["geometry_config"]["internal_membrane_relative_radii"] + [1])*external_radius,   # R_1, R_2, R_3
        D     = np.array(diffusions),          # D_1, D_2
        K     = k,
        P     = np.tile(permeabilities, (len(geometry["geometry_config"]["internal_membrane_relative_radii"]), 1)), #np.array(permeabilities),          # internal permeabilities
        P_out = np.array(permeabilities),          # outer boundary permeabilities
        q_inf = np.array(external_concentrations),          # X_1=1, X_2=0 outside
    )

    sol = solve(params)
    # evaluate on a fine grid
    num_grid_points = 100
    r = np.linspace(params.radii[-1]*0.0001, params.radii[-1], num_grid_points)
    X = evaluate_solution(sol, r)

    for species_idx, species in enumerate(species_names):
        ax.plot(r/params.radii[-1], X[species_idx], label=str(species), ls = ":")
    ax.legend()

    fig.savefig(file_to_create, dpi = 300)





if __name__ == "__main__":
    # get all combined folders
    subfolders = [ f.path for f in os.scandir(pathlib.Path(__file__).parent.resolve()) if f.is_dir() ]
    combined_folders = [folder for folder in subfolders if "combined" in folder]
    for folder in combined_folders:
        #if "000002" not in folder:
        #   continue
        print(folder)
        #try:
        plot_numerics_and_analytical_solution(folder)
        #except:
        #    print(Path(folder).name, "failed")

    # load information
    # plot numerical solution as normal
    # construct input needed by analytical "solver"
    # add solution of analytical solver 
    
    # python data/00b_comparison_numerical_solution_to_analytical_linear_solutions/comparisons.py data/00b_comparison_numerical_solution_to_analytical_linear_solutions



"""
    k = np.zeros((system_geometry_dict["geometry_config"]["num_regions"], len(species_names), len(species_names)))
    for region_idx in range(len(k)):
        for reaction_idx in range(len(species_names)-1):
            start_species = species_names[reaction_idx]

            if start_species in spontaneous_reactions_df["start_species"].unique():
                k[region_idx][reaction_idx] = spontaneous_reactions_df.loc[(spontaneous_reactions_df["start_species"] == start_species), "k"].item()
            elif start_species in enzymatic_reactions_df["start_species"].unique():
                k_cat = enzymatic_reactions_df.loc[(enzymatic_reactions_df["start_species"] == start_species), "k_cat"].item()
                k_M = enzymatic_reactions_df.loc[(enzymatic_reactions_df["start_species"] == start_species), "k_M"].item()
                enzyme = enzymatic_reactions_df.loc[(enzymatic_reactions_df["start_species"] == start_species), "enzyme"].item()
                k[region_idx][reaction_idx] = k_cat / k_M * enzyme_concentrations[enzyme][region_idx]
            else:
                raise ValueError("start_species", start_species, "not found")


"""