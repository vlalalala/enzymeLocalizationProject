
#%%
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from scipy import integrate

final_mesh = res_a.x
solution_values_at_mesh = res_a.y
solution_derivatives_at_mesh = res_a.yp

success = res_a.success

integrals = {species.name: 0 for species in reaction_network.species}
x_values = np.linspace(0, radius, num = 100)
num_enzymes = len(reaction_network.enzymes)

enzyme_max_concentration = max(enzyme.concentration for enzyme in reaction_network.enzymes)

norm = mcolors.Normalize(vmin=0, vmax=enzyme_max_concentration)
cmap = cm.get_cmap('Oranges')
scalar_map = cm.ScalarMappable(norm=norm, cmap=cmap)

fig, ax = plt.subplots(num_enzymes + 1,1, figsize = (5,3), sharex = True,
                       gridspec_kw={'height_ratios': [1]*num_enzymes + [num_enzymes*2]})
fig.subplots_adjust(hspace=0)
ax[0].set_title(f"steady state distribution, computed with {len(res_a.x)} nodes")

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
    for localizationTuple in enzyme.localization:
        min_range = localizationTuple.minMaxLoc[0]
        max_range = localizationTuple.minMaxLoc[1]
        ax[enzyme_idx].fill_between(
            x_values, 0, 1, where=(x_values >= min_range) & (x_values <= max_range),
            color=color)
    ax[enzyme_idx].set_yticks([])

for species in reaction_network.species:
    sol = res_a.sol(x_values)[species_variable_indices[species]["prim"]]
    # Found solution for y as scipy.interpolate.PPoly instance, a C1 continuous
    # cubic spline.
    integral = integrate.quad(lambda r: res_a.sol(r)[species_variable_indices[species]["prim"]] * 4 * np.pi * r**2, 0, radius)
    integrals[species.name] = integral[0] # [0] is value [1] is error
    ax[-1].plot(x_values, sol, label=f"{species.name}: {int(np.round(integral[0]))} nM")

for node_x in res_a.x:
    ax[-1].axvline(node_x, ymin = 0.95, ymax = 1, c = "k", linewidth = 1)

ax[-1].set_xlabel("radius")
ax[-1].set_xlabel("distance to origin")
ax[-1].set_xlim([0, radius])
ax[-1].set_ylabel("concentration")

fig.savefig(os.path.join(violacein_folder, "test.png"), dpi = 300,
    bbox_inches = "tight")


# Assume scalar_map is your ScalarMappable (from cmap and norm)
fig_colorbar = plt.figure(figsize=(2, 4))  # tall figure for vertical bar
# Add axes for colorbar (full figure area)
cbar_ax = fig_colorbar.add_axes([0.2, 0.05, 0.3, 0.9])  # [left, bottom, width, height]
# Create the colorbar
cbar = fig_colorbar.colorbar(scalar_map, cax=cbar_ax, orientation='vertical')
cbar.set_label("enzyme concentration")
# Save the colorbar figure
fig_colorbar.savefig(os.path.join(violacein_folder, "test_colorbar.png"), bbox_inches='tight', dpi= 300)


fig_legend = plt.figure(figsize=(3, 2))
ax_legend = fig_legend.add_subplot(111)
ax_legend.axis('off')

handles, labels = ax[-1].get_legend_handles_labels()
legend = ax_legend.legend(handles, labels, loc='center')

fig_legend.savefig(os.path.join(violacein_folder, "test_colorbar.png"), bbox_inches='tight', dpi= 300)



#%%


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

