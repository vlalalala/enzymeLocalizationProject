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
from pathlib import Path

def get_analytical_solution_function_linear(c_ext, R, p, D, k):
    l = np.sqrt(D/k)
    return lambda x: (p/D) * c_ext * np.cosh(x/l) / (1/l * np.sinh(R/l) + p/D * np.cosh(R/l))

def get_analytical_solution_function_spherical(c_ext, R, p, D, k):
    Lambda = np.sqrt(k/D)
    A =-p/D * c_ext * R**2 / (np.exp(Lambda * R) * (Lambda*R-1 + p*R/D) + np.exp(-Lambda*R) * (Lambda * R + 1 - p*R/D))
    return lambda r: (A * np.exp(-Lambda * r ) - A*np.exp(Lambda * r)) / r
    #return lambda r: (p/D) * c_ext * R**2 * np.sinh(Lambda * r) / r * 1 / ((Lambda * R - 1 + p * R / D) * np.sinh(Lambda * R) + p * R / D * np.cosh(Lambda * R) - np.sinh(Lambda * R)) 

def plot_data(folder):
    R = 1e-5
    c_ext = 25e-3
    p_list = list(np.logspace(-6, -2, num = 5))# + [1e20]#25e-6
    D = 6.6e-11
    k = 1e1
    fig, ax = plt.subplots(len(p_list), 1, figsize = (4,3*len(p_list)))
    for i, p in enumerate(p_list):
        linear_fct = get_analytical_solution_function_linear(c_ext=c_ext, R=R, p=p, D=D, k=k)
        spherical_fct = get_analytical_solution_function_spherical(c_ext=c_ext, R=R, p=p, D=D, k=k)

        x_values = np.linspace(0, R, num = 1000)
        ax[i].plot(x_values/R, [linear_fct(x_value) for x_value in x_values], label = "linear case")
        ax[i].plot(x_values/R, [spherical_fct(x_value) for x_value in x_values], label = "spherical case")
        ax[i].plot([1, 1.1], [c_ext, c_ext])
        ax[i].axvline(1, ls = ":", c = "k")
        ax[i].set_ylabel("concentration")
        ax[i].set_title(f"p = {p}")
    ax[-1].set_xlabel("location x/R")
    ax[0].legend()
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "result.png"), dpi = 300)

if __name__ == "__main__":
    # Load all the passed information
    folder = Path(__file__).resolve().parent
    plot_data(folder)
    # python data/00_comparison_diffusion_linearVsSpherical/analysis.py
