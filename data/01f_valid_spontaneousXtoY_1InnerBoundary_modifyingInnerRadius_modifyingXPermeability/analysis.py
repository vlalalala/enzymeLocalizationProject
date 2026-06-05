import sys
import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src'))
from decimal import Decimal
from plot_bvp_solution import (
    calculate_analytical_solution_1_spontaneous_reaction_1_region,
    calculate_analytical_solution_1_spontaneous_reaction_2_regions
)
import ast
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})

def plot_curve(folder):
    R = 1e-5
    k = 1
    D = 6.6e-11
    ext = 90e-8
    p = 25e-6
    r_I = 0.5
    c_with_inner = calculate_analytical_solution_1_spontaneous_reaction_2_regions(
                k = k, D = D, R = R, p = p, ext = ext, r_inner=r_I
            )
    fig, ax = plt.subplots(1,1, figsize = (4,3))
    radii = np.linspace(R*0.01, R, num = 100)
    ax.plot(radii, [c_with_inner(r) for r in radii])
    ax.set_ylim([-ext*0.1, ext*1.1])
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "test.png"), dpi = 300)

def test():
    R = 1e-5
    k = 1
    D = 6.6e-11
    ext = 90e-8
    p = 25e-6
    for r_inner in [1e-9, 1e-6, 1e-3]:
        c1 = calculate_analytical_solution_1_spontaneous_reaction_1_region(
            k,D,R,p,ext
        )
        c2 = calculate_analytical_solution_1_spontaneous_reaction_2_regions(
            k,D,R,p,ext,r_inner
        )

        print(
            r_inner,
            c1(R),
            c2(R),
            c2(R)/c1(R)
        )

def plot_data(folder):
    R = 1e-5
    k = 1
    D = 6.6e-11
    ext = 90e-8
    #p = 25e-6
    R_to_test = [1e-6, 1e-5, 1e-4]
    k_to_test = [0.1, 1, 10, 100]
    fig, ax = plt.subplots(1, 2, figsize = (4,3), gridspec_kw={'width_ratios': [1, 0.1]})
    tol = 1e-5 * ext
    P = np.zeros((100, 100))
    inner_radii = np.zeros((100, 100))
    Z = np.zeros((100, 100))
    for i, p in enumerate(np.logspace(-10, -1, num = 100)):
        c_without_inner = calculate_analytical_solution_1_spontaneous_reaction_1_region(
            k = k, D = D, R = R, p = p, ext = ext
        )
        flux_without = ext - c_without_inner(R)
        #print("without", flux_without)
        for j, r_I in enumerate(np.linspace(0.001*R, 0.999*R, num = 100)):
            try:
                c_with_inner = calculate_analytical_solution_1_spontaneous_reaction_2_regions(
                    k = k, D = D, R = R, p = p, ext = ext, r_inner=r_I
                )
                flux_with = ext - c_with_inner(R) # positive quantity
                #print("with", flux_with)
                P[i,j] = p
                inner_radii[i,j] = r_I / R
                #if abs(flux_without) < tol:
                #    Z[i,j] = np.nan
                if flux_with < 0:
                    raise ValueError
                else:
                    relative_flux = flux_with / flux_without
                    #if relative_flux < 0:
                    #Z[i,j] = np.nan
                    #else:
                    Z[i,j] = relative_flux
            except:
                print("didnt work")

    print("min:", np.nanmin(Z))
    print("max:", np.nanmax(Z))
    print("any <= 0:", np.any(Z <= 0))
    print("any NaN:", np.any(np.isnan(Z)))
   
    mesh0 = ax[0].pcolormesh(P, inner_radii, Z, cmap='viridis', shading='auto',
                             #norm=LogNorm(vmin=np.nanmin(Z), vmax=np.nanmax(Z))
                             #vmin = 0.0, vmax = 1
    )
    fig.colorbar(mesh0, cax=ax[1],
        label='flux with one inner membrane \n divided by flux without an inner membrane')

    # Set log scale on whichever axes need it
    ax[0].set_xlabel(r"permeability of species X $p_X$ / m$\cdot$s$^{-1}$")
    ax[0].set_ylabel(r"relative inner radius $r_I/R$")
    ax[1].set_box_aspect(10)
    ax[0].set_box_aspect(1)
    ax[0].set_xscale("log")
    
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "fluxes.png"), dpi = 300)

# since this is a first order reaction, it's the same to change the diffusion constant than
# to change the rate

def plot_data_mesh(folder):
    D = 6.6e-11
    ext = 90e-8
    R_to_test = [1e-6, 1e-5, 1e-4]
    k_to_test = [0.1, 1, 10, 100]
    fig, ax = plt.subplots(len(k_to_test), len(R_to_test) + 1,
        figsize = (4*len(R_to_test),3*len(k_to_test)), gridspec_kw={'width_ratios': [1]*len(R_to_test)+[0.1]})
    for R_idx, R in enumerate(R_to_test):
        for k_idx, k in enumerate(k_to_test):
            P = np.zeros((100, 100))
            inner_radii = np.zeros((100, 100))
            Z = np.zeros((100, 100))
            for i, p in enumerate(np.logspace(-10, -1, num = 100)):
                c_without_inner = calculate_analytical_solution_1_spontaneous_reaction_1_region(
                    k = k, D = D, R = R, p = p, ext = ext
                )
                flux_without = ext - c_without_inner(R)
                for j, r_I in enumerate(np.linspace(0.001*R, 0.999*R, num = 100)):
                    c_with_inner = calculate_analytical_solution_1_spontaneous_reaction_2_regions(
                        k = k, D = D, R = R, p = p, ext = ext, r_inner=r_I
                    )
                    flux_with = ext - c_with_inner(R) # positive quantity
                    #print("with", flux_with)
                    P[i,j] = p
                    inner_radii[i,j] = r_I / R
                    if flux_with < 0 or flux_without <0:
                        raise ValueError
                    else:
                        relative_flux = flux_with / flux_without
                        Z[i,j] = relative_flux

        
            mesh0 = ax[k_idx][R_idx].pcolormesh(P, inner_radii, Z, cmap='viridis', shading='auto',
                                    #norm=LogNorm(vmin=np.nanmin(Z), vmax=np.nanmax(Z))
                                    vmin = 0.0, vmax = 1
            )
            fig.colorbar(mesh0, cax=ax[k_idx][-1], label='flux with one inner membrane \n divided by flux without an inner membrane')

            if k_idx == len(k_to_test)-1:
            # Set log scale on whichever axes need it
                ax[k_idx][R_idx].set_xlabel(r"permeability of species X $p_X$ / m$\cdot$s$^{-1}$")
            if R_idx == 0:
                ax[k_idx][R_idx].set_ylabel(r"relative inner radius $r_I/R$")
            #if k_idx == 0:
            ax[k_idx][R_idx].set_title(r"external radius $R=$"+" {:.1e}".format(R)+"\n"+r"reaction rate $k=$"+" {:.1e}".format(k))
            ax[k_idx][-1].set_box_aspect(10)
            ax[k_idx][R_idx].set_box_aspect(1)
            ax[k_idx][R_idx].set_xscale("log")
    
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "fluxes_mesh.png"), dpi = 300)


if __name__ == "__main__":
    # Load all the passed information
    FOLDER_TO_SOLVE = sys.argv[1]
    test()
    plot_data(FOLDER_TO_SOLVE)
    plot_data_mesh(FOLDER_TO_SOLVE)
    #plot_curve(FOLDER_TO_SOLVE)
    # python data/01f_valid_spontaneousXtoY_1InnerBoundary_modifyingInnerRadius_modifyingXPermeability/analysis.py data/01f_valid_spontaneousXtoY_1InnerBoundary_modifyingInnerRadius_modifyingXPermeability


