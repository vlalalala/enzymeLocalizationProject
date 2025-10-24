
#%%
import sys
import os
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy import integrate
from auxiliary_functions_using_standard_library import pickle_load_binary
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme

def create_bvp_result_plots(bvp_result, bvp_variable_indices, reaction_network):
    
    x_values = np.linspace(0, max(bvp_result.x), num = 100)
    integrals = {}
    num_enzymes = len(reaction_network.enzymes)

    enzyme_max_concentration = max(enzyme.concentration for enzyme in reaction_network.enzymes)
    norm = mcolors.Normalize(vmin=0, vmax=enzyme_max_concentration)
    cmap = matplotlib.colormaps['Oranges']
    scalar_map = cm.ScalarMappable(norm=norm, cmap=cmap)

    # MAIN PLOT
    fig, ax = plt.subplots(num_enzymes + 1,1, figsize = (5,3), sharex = True,
                        gridspec_kw={'height_ratios': [1]*num_enzymes + [num_enzymes*2]})
    fig.subplots_adjust(hspace=0)
    ax[0].set_title(f"steady state distribution, computed with {len(bvp_result.x)} nodes")
    # Main plot enzyme concentration
    for enzyme_idx, enzyme in enumerate(reaction_network.enzymes):
        ax[enzyme_idx].set_ylabel(enzyme.name, rotation=0, labelpad=20)
        # Get the label Text object and adjust its position
        label = ax[enzyme_idx].yaxis.get_label()
        # Move it vertically to center (around 0.5 in axes coords) and keep horizontal offset from labelpad
        label.set_verticalalignment('center')  # ensure vertical alignment
        # Manually set position: x controls horizontal offset, y controls vertical position (0 bottom, 0.5 center, 1 top)
        label.set_position((label.get_position()[0], 0.5))
        # Colorbar color
        color = scalar_map.to_rgba(enzyme.concentration)
        # Fill x-ranges where enzymes are at
        for localizationTuple in enzyme.localization:
            min_range = localizationTuple.minMaxLoc[0]
            max_range = localizationTuple.minMaxLoc[1]
            ax[enzyme_idx].fill_between(
                x_values/max(x_values), 0, 1, where=(
                    x_values/max(x_values) >= min_range)
                    & (x_values/max(x_values) <= max_range),
                color=color)
        ax[enzyme_idx].set_yticks([])
    
    # Main plot species concentration
    for species in reaction_network.species:
        # bvp_result.sol is the found solution for y as scipy.interpolate.PPoly instance, a C1 continuous
        # cubic spline.
        integrand = lambda r: bvp_result.sol(r)[bvp_variable_indices[species]["prim"]]* 10**3 * 4 * np.pi * r**2
        integral, _ = integrate.quad(integrand, 0, max(bvp_result.x))
        integrals[species.name] = integral
    
    # Find smallest exponent (in base 10) across all values
    if integrals:
        min_exp = min(int(np.floor(np.log10(abs(val)))) if val != 0 else 0 for val in integrals.values())
    else:
        min_exp = 0  # fallback

    # Format legend labels
    labels_with_integrals = {}
    for species in reaction_network.species:
        val = integrals[species.name]
        scaled_val = val / (10 ** min_exp)  # Adjust to common exponent
        label = f"{species.name}: {scaled_val:.0f}e{min_exp} mol"
        labels_with_integrals[species.name] = label
    
    # Plot and assign labels
    for species in reaction_network.species:    
        solution_y_values_to_plot = bvp_result.sol(x_values)[bvp_variable_indices[species]["prim"]]
        ax[-1].plot(x_values/max(x_values), solution_y_values_to_plot, label=labels_with_integrals[species.name])
    
    # Main plot mesh points visualization
    for node_x in bvp_result.x:
        ax[-1].axvline(node_x/max(x_values), ymin = 0.95, ymax = 1, c = "k", linewidth = 1)
    # Main plot labels
    ax[-1].set_xlabel("r/R")
    ax[-1].set_xlim([0, 1])
    ax[-1].set_ylabel("concentration / M")

    # ENZYME CONCENTRATION COLORBAR
    # Assume scalar_map is your ScalarMappable (from cmap and norm)
    fig_colorbar = plt.figure(figsize=(2, 4))  # tall figure for vertical bar
    # Add axes for colorbar (full figure area)
    cbar_ax = fig_colorbar.add_axes([0.2, 0.05, 0.3, 0.9])  # [left, bottom, width, height]
    # Create the colorbar
    cbar = fig_colorbar.colorbar(scalar_map, cax=cbar_ax, orientation='vertical')
    cbar.set_label("enzyme concentration")
    # SPECIES CONCENTRATION LEGEND
    fig_legend = plt.figure(figsize=(3, 2))
    ax_legend = fig_legend.add_subplot(111)
    ax_legend.axis('off')
    handles, labels = ax[-1].get_legend_handles_labels()
    ax_legend.legend(handles, labels, loc='center')

    return fig, fig_colorbar, fig_legend

#%%
if __name__ == "__main__":
    # Load all the information
    folder_to_plot = sys.argv[1]
    solved_reaction_network = pickle_load_binary(os.path.join(folder_to_plot, ".REACTION_NETWORK_pickle"))
    boundary_value_problem_result = pickle_load_binary(os.path.join(folder_to_plot, ".BOUNDARY_VALUE_PROBLEM_RESULT_pickle"))
    boundary_value_problem_variable_indices = pickle_load_binary(os.path.join(folder_to_plot, ".BOUNDARY_VALUE_PROBLEM_VARIABLE_INDICES_pickle"))


    # Run
    success = boundary_value_problem_result.success
    #print("success", success)
    #if not success:
    #    raise ImportError(f"The imported result for the bvp in folder {folder_to_plot} was not successful.")
    
    main_fig, enzyme_concentration_colorbar_fig, species_legend_fig =  create_bvp_result_plots(
        boundary_value_problem_result, boundary_value_problem_variable_indices,
        solved_reaction_network
    )

    # Save plots
    main_fig.savefig(os.path.join(folder_to_plot, "bvp_result.png"), dpi = 300,
        bbox_inches = "tight")
    enzyme_concentration_colorbar_fig.savefig(
        os.path.join(folder_to_plot, "bvp_enzyme_colorbar.png"), bbox_inches='tight', dpi= 300)
    species_legend_fig.savefig(
        os.path.join(folder_to_plot, "bvp_species_legend.png"), bbox_inches='tight',
        dpi= 300
    )


#%%
#species_integral_sum = np.sum([
#    integrals[species.name]
#    for species in reaction_network.species])
#annotation_list = [
#    f"{species.name}: {int(np.round(integrals[species.name]/species_integral_sum * 100))}%"
#    for species in reaction_network.species]
#annotation = ", ".join(annotation_list)

#cbar = fig.colorbar(scalar_map, ax=ax[:num_enzymes], orientation='vertical', fraction=0.02, pad=0.04)
#cbar.set_label("concentration")
#ax[-1].annotate(annotation, (0.05, 0.9), xycoords = "axes fraction", va = "top", ha = "left", fontsize = 9)

