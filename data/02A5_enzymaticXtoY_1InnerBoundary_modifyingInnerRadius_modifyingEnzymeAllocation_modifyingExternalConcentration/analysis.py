import sys
import matplotlib.pyplot as plt
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src'))
from auxiliary_functions_using_standard_library import load_json
import re
import pandas as pd
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
from auxiliary_functions import read_yaml_file
import ast
from find_matching_parameter_value_combinations import filter_combined_folders
import matplotlib.image as mpimg
from pathlib import Path


def calculate_average_concentration_of_X_weighted_by_enzyme_allocation(combined_folder):
    system_geometry_file = os.path.join(combined_folder, "system_geometry_for_convergence.json")
    if not os.path.isfile(system_geometry_file):
        return None
    try:
        system_geometry = load_json(system_geometry_file)
        concentrations = load_json(os.path.join(combined_folder, ".species_steady_state_concentrations.json"))
        enzymes_df = pd.read_csv(os.path.join(combined_folder, "enzymes.csv"))
    except:
        return None
    enzyme_allocation = ast.literal_eval(enzymes_df.loc[enzymes_df['name'] == "A", 'allocation'].iloc[0])
    average_concentration = 0
    for region_idx in range(2):
        radii_in_region = system_geometry["geometry_config"]["mesh_points_in_regions"][region_idx]
        X_concentrations_in_region = []
        for radius_idx in range(len(radii_in_region)):
            X_concentrations_in_region.append(concentrations[region_idx][radius_idx]["X"])
        concentration = 0
        for i in range(len(radii_in_region) - 1):
            r0 = radii_in_region[i]
            r1 = radii_in_region[i+1]
            c0 = X_concentrations_in_region[i]
            c1 = X_concentrations_in_region[i+1]
            
            # Volume-weighted average concentration in shell (assuming linear c(r))
            c_avg = (c0 * (3*r0**2 + 2*r0*r1 + r1**2) + c1 * (r0**2 + 2*r0*r1 + 3*r1**2)) \
                    / (4 * (r0**2 + r0*r1 + r1**2))
            
            shell_volume = r1**3 - r0**3  # 4π/3 cancels with denominator
            concentration += c_avg * shell_volume

        concentration /= (radii_in_region[-1]**3 - radii_in_region[0]**3)
        average_concentration += concentration * enzyme_allocation[region_idx]

    return average_concentration

def get_data(folder):
    data = {}
    for folder_name in os.listdir(folder):
        match = re.match(r'^combined_(\d{6})$', folder_name)
        combined_folder = os.path.join(folder, folder_name)
        if match:
            index = match.group(1)  # keeps it as a string, preserving leading zeros
            geometry = read_yaml_file(os.path.join(combined_folder, "parameters_geometry.yaml"))
            internal_relative_radius = geometry["geometry_config"]["internal_membrane_relative_radii"][0]
            species_df = pd.read_csv(os.path.join(combined_folder, "species.csv"))
            X_external_concentration = species_df.loc[
                (species_df["name"] == "X"),
                "external_concentration"].item()
            enzymes_df = pd.read_csv(os.path.join(combined_folder, "enzymes.csv"))
            allocation_str = enzymes_df.loc[(enzymes_df["name"] == "A"), "allocation"].item()
            allocation = ast.literal_eval(allocation_str)
            allocation_in_external = allocation[1]
            fluxes_file = os.path.join(combined_folder, "fluxes.json")
            if os.path.isfile(fluxes_file):
                Y_flux = load_json(fluxes_file)["Y"]
            else:
                Y_flux = None
            average_X_concentration_weighted = calculate_average_concentration_of_X_weighted_by_enzyme_allocation(combined_folder)
            data[index] = (internal_relative_radius, allocation_in_external, Y_flux, X_external_concentration,index, average_X_concentration_weighted)
    return data


def plot_steady_states(folder):
    data = get_data(folder)
    df = pd.DataFrame(data.values(),
                      columns=['internal_relative_radius', 'allocation_in_external', 'flux', 'X_external_concentration', 'index'])
    # one file per X_external_concentration
    inner_radii = np.sort(df["internal_relative_radius"].unique())
    allocation_in_external = np.sort(df["allocation_in_external"].unique())
    X_external_concentrations = np.sort(df["X_external_concentration"].unique())
    for X_external_concentration in X_external_concentrations:
        fig, ax = plt.subplots(len(allocation_in_external), len(inner_radii),
                               figsize = (4*len(inner_radii), 3*len(allocation_in_external)))
        fig.suptitle(f"X_ext = {X_external_concentration}")
        for allocation_idx, allocation in enumerate(allocation_in_external):
            for inner_radius_idx, inner_radius in enumerate(inner_radii):
                #current_df = df[(df["internal_relative_radius"]==inner_radius)&(df["allocation_in_external"]==allocation)&(df['X_external_concentration']==X_external_concentration)]
                combinations = filter_combined_folders(
                    combined_root=folder,
                    criteria_yaml={
                        "options_parameters_geometry": {
                            "geometry_config": {
                                "internal_membrane_relative_radii": [inner_radius]
                            }
                        }
                    },
                    criteria_csv=
                        {
                            "options_species": {"name": {"X": {"external_concentration": f"{X_external_concentration}"}}},
                            #"options_enzymes": {"name": {"A": {"allocation": "{0:"+ str(1-allocation) +", 1:"+str(allocation)+"}"}}},
                        }
                )
                # now select from combinations the one with the correct allocation
                enzyme_allocation_strings = {}
                for path_folder in combinations:
                    enzyme_df = pd.read_csv(path_folder / "enzymes.csv")
                    allocation_str = enzyme_df.loc[enzyme_df["name"] == "A", 'allocation'].iloc[0]
                    allocation_dict = ast.literal_eval(allocation_str)
                    allocation_outermost = allocation_dict[1]
                    enzyme_allocation_strings.update({path_folder: allocation_outermost})
                # find the path for which the allocation is closest
                combination, _ = min(enzyme_allocation_strings.items(), key=lambda x: abs(allocation - x[1]))
                
                #if len(combinations)!=1:
                #    raise ValueError(combinations)
                #combination = combinations[0]#only one available
                file_to_plot = combination / "solver_iteration_data" / "interpolation_iteration_nr_0_final_concentrations.png"
                if os.path.isfile(file_to_plot):
                    img = mpimg.imread(str(file_to_plot))
                ax[allocation_idx][inner_radius_idx].imshow(img)
                ax[allocation_idx][inner_radius_idx].axis("off")
                ax[allocation_idx][inner_radius_idx].set_title(f"Xext = {X_external_concentration}, \n inner radius = {inner_radius} \n {Path(combination).name}")
        fig.savefig(os.path.join(folder, f"complete_steadyStates_Xext_{X_external_concentration}.png"), dpi = 300)


def plot_data(folder):
    data = get_data(folder)
    print(data)

    df = pd.DataFrame(data.values(), columns=['x', 'y', 'flux', 'X_external_concentration', 'index'])
    external_concentrations = sorted(df["X_external_concentration"].unique())
    fig, ax = plt.subplots(1, 2, figsize = (5,6*len(external_concentrations)), gridspec_kw={'width_ratios': [1, 0.1]})
    
    ############
    # Plot fluxes
    ############
    for external_concentration_idx, external_concentration in enumerate(external_concentrations):
        current_df = df[df["X_external_concentration"]==external_concentration]
        z_grid = current_df.pivot(index='y', columns='x', values='flux')
        x_points = z_grid.columns.values
        y_points = z_grid.index.values
        z_points = z_grid.values  # 2D array of shape (len(y), len(x))
    
        mesh0 = ax[0].pcolormesh(x_points, y_points, z_points, cmap='viridis', shading='auto',
                                #norm = Normalize(vmin=0, vmax=1)
                                #norm=LogNorm(vmin=np.nanmin(z_points), vmax=np.nanmax(z_points))
        )
        #ax[external_concentration_idx][0].scatter(current_df['x'], current_df['y'], color='red', s=5, zorder=5)
        #for _, point in df.iterrows():
        #    ax[0].annotate(f"{point['index'].lstrip('0')}", (point['x'], point['y']))

        fig.colorbar(mesh0, cax=ax[1], label='flux')

        # Set log scale on whichever axes need it
        ax[0].set_xlabel("relative radius of inner membrane r*/R")
        ax[0].set_ylabel("proportion of enzyme \n in outermost region")
        ax[1].set_box_aspect(10)
        ax[0].set_box_aspect(1)
        ax[0].set_title(f"X_ext = {external_concentration} mol")
    
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "fluxes.png"), dpi = 300)

def plot_data_with_concentration_weighted(folder):
    data = get_data(folder)
    df = pd.DataFrame(data.values(),
                      columns=['x', 'y', 'flux',
                               'X_external_concentration', 'index',
                               'average_X_concentration_weighted']
    )
    external_concentrations = sorted(df["X_external_concentration"].unique())
    fig, ax = plt.subplots(len(external_concentrations), 4, figsize = (10,3*len(external_concentrations)),
                           gridspec_kw={'width_ratios': [1, 0.1]*2})
    
    ############
    # Plot fluxes
    ############
    for external_concentration_idx, external_concentration in enumerate(external_concentrations):
        current_df = df[df["X_external_concentration"]==external_concentration]
        z_grid = current_df.pivot(index='y', columns='x', values='flux')
        x_points = z_grid.columns.values
        y_points = z_grid.index.values
        z_points = z_grid.values  # 2D array of shape (len(y), len(x))
    
        mesh0 = ax[external_concentration_idx][0].pcolormesh(x_points, y_points, z_points, cmap='viridis', shading='auto',
                                #norm = Normalize(vmin=0, vmax=1)
                                #norm=LogNorm(vmin=np.nanmin(z_points), vmax=np.nanmax(z_points))
        )
        #ax[external_concentration_idx][0].scatter(current_df['x'], current_df['y'], color='red', s=5, zorder=5)
        #for _, point in df.iterrows():
        #    ax[0].annotate(f"{point['index'].lstrip('0')}", (point['x'], point['y']))

        fig.colorbar(mesh0, cax=ax[external_concentration_idx][1], label='flux')

        # Set log scale on whichever axes need it
        ax[external_concentration_idx][0].set_xlabel("relative radius of inner membrane r*/R")
        ax[external_concentration_idx][0].set_ylabel("proportion of enzyme \n in outermost region")
        ax[external_concentration_idx][1].set_box_aspect(10)
        ax[external_concentration_idx][0].set_box_aspect(1)
        ax[external_concentration_idx][0].set_title("X_ext = {:.1e}".format(external_concentration) +r" mol $\cdot$ m$^{-3}$")

        ##################
        # concentrations
        ##################
        z_grid = current_df.pivot(index='y', columns='x', values='average_X_concentration_weighted')
        z_points = z_grid.values  # 2D array of shape (len(y), len(x))
    
        mesh1 = ax[external_concentration_idx][2].pcolormesh(x_points, y_points, z_points, cmap='viridis', shading='auto',
                                #norm = Normalize(vmin=0, vmax=1)
                                #norm=LogNorm(vmin=np.nanmin(z_points), vmax=np.nanmax(z_points))
        )
        #ax[external_concentration_idx][0].scatter(current_df['x'], current_df['y'], color='red', s=5, zorder=5)
        #for _, point in df.iterrows():
        #    ax[0].annotate(f"{point['index'].lstrip('0')}", (point['x'], point['y']))

        fig.colorbar(mesh0, cax=ax[external_concentration_idx][3], label='average (weighted by \n volume and enzyme allocation) \n substrate concentration')

        # Set log scale on whichever axes need it
        ax[external_concentration_idx][2].set_xlabel("relative radius of inner membrane r*/R")
        ax[external_concentration_idx][2].set_ylabel("proportion of enzyme \n in outermost region")
        ax[external_concentration_idx][3].set_box_aspect(10)
        ax[external_concentration_idx][2].set_box_aspect(1)
        #ax[external_concentration_idx][0].set_title(f"X_ext = {external_concentration} mol")        
    
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "fluxes2.png"), dpi = 300)


if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    #plot_steady_states(FOLDER_TO_SOLVE)
    plot_data_with_concentration_weighted(FOLDER_TO_SOLVE)
    # python data/02A5_enzymaticXtoY_1InnerBoundary_modifyingInnerRadius_modifyingEnzymeAllocation_modifyingExternalConcentration/analysis.py data/02A5_enzymaticXtoY_1InnerBoundary_modifyingInnerRadius_modifyingEnzymeAllocation_modifyingExternalConcentration
    #print(filter_combined_folders(
    #    combined_root=FOLDER_TO_SOLVE,
    #    criteria_yaml={},
    #    criteria_csv=
    #        {
    #            "options_species": {"name": {"X": {"external_concentration": "10.0"}}},
    #            #"options_enzymes": {"name": {"A": {"allocation": "{0:"+ str(1-allocation) +", 1:"+str(allocation)+"}"}}},
    #        }
    #))