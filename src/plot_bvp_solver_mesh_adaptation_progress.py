from pathlib import Path
import os
import matplotlib.pyplot as plt
import sys
from auxiliary_functions_using_standard_library import pickle_load_binary
import pandas as pd
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

def plot_convergence_progress(
    folder_to_solve,
    reaction_network
):  
    files = list(Path(folder_to_solve).glob(".progress_log_interpolating_*_times.csv"))
    fig, ax = plt.subplots(4, 1, figsize = (3,9))
    files.sort()
    colors = {species.name: color for species, color in zip(reaction_network.species, plt.cm.tab10.colors)}
    
    iterations_offset = 0
    for file_idx, file in enumerate(files):
        print(file_idx, iterations_offset)
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

    

    

    #ax[0].set_ylabel(
    #    r"Relative difference $\frac{\lvert \Phi_{\text{react}} - \Phi_{\text{bound}} \rvert}"
    #    r"{\max\!\left(\lvert \Phi_{\text{react}} \rvert, \lvert \Phi_{\text{bound}} \rvert\right)}$"
    #)
    
    fig.tight_layout()
    fig.savefig(os.path.join(folder_to_solve, "convergence.png"), bbox_inches = "tight", dpi=300)
    plt.close()

if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    REACTION_NETWORK = pickle_load_binary(os.path.join(FOLDER_TO_SOLVE, ".pickled_reaction_network"))
    plot_convergence_progress(FOLDER_TO_SOLVE, REACTION_NETWORK)
