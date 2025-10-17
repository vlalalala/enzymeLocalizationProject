

def write_wl_script(params, output_wl_path, output_csv_path):
    wl_template = f"""

ClearAll["Global`*"];

(* --- Parameters --- *)
diff = {params['diff']};
p = {params['p']};
uExt = {params['uExt']};
vExt = {params['vExt']};
k = {params['k']};
c = {params['c']};

(* --- Reaction functions --- *)
f[u_, v_] := k*u;
g[u_, v_] := -k*u;

(* --- Left domain solver: x in [0,c] *)
leftSolver = ParametricNDSolveValue[
  {
   -diff*u''[x] + f[u[x], v[x]] == 0,
   -diff*v''[x] + g[u[x], v[x]] == 0,
   u[0] == u0, v[0] == v0,
   u'[0] == 0, v'[0] == 0
  },
  {u, v},
  {x, 0, c},
  {u0, v0}
];

(* --- Right domain solver: x in [c,1] *)
rightSolver = ParametricNDSolveValue[
  {
   -diff*u''[x] + f[u[x], v[x]] == 0,
   -diff*v''[x] + g[u[x], v[x]] == 0,
   u[c] == uC, v[c] == vC,
   u'[1] == (p/diff)*(uExt - u[1]),
   v'[1] == (p/diff)*(vExt - v[1])
  },
  {u, v},
  {x, c, 1},
  {uC, vC}
];

(* --- Matching function for FindRoot --- *)
matchingResiduals[{u0_?NumericQ, v0_?NumericQ, alpha_?NumericQ, beta_?NumericQ}] := Module[
  {left, right, uLc, vLc, duLc, dvLc, uRC, vRC, duRc, dvRc, uR1, vR1, duR1, dvR1},
  
  (* Solve left *)
  left = leftSolver[u0, v0];
  uLc = left[[1]][c] // N;
  vLc = left[[2]][c] // N;
  duLc = left[[1]]'[c] // N;
  dvLc = left[[2]]'[c] // N;
  
  (* Right initial values using alpha/beta *)
  uRC = alpha*uLc;
  vRC = beta*vLc;
  
  (* Solve right *)
  right = rightSolver[uRC, vRC];
  duRc = right[[1]]'[c] // N;
  dvRc = right[[2]]'[c] // N;
  uR1 = right[[1]][1] // N;
  vR1 = right[[2]][1] // N;
  duR1 = right[[1]]'[1] // N;
  dvR1 = right[[2]]'[1] // N;
  
  (* Residuals: flux continuity + Robin at x=1 *)
  {
    diff*duLc - diff*duRc,          (* u flux continuity *)
    diff*dvLc - diff*dvRc,          (* v flux continuity *)
    duR1 - (p/diff)*(uExt - uR1),  (* u Robin condition *)
    dvR1 - (p/diff)*(vExt - vR1)   (* v Robin condition *)
  }
];

(* --- Initial guesses for u0, v0, alpha, beta --- *)
initGuess = {1.0, 0.5, 1.0, 1.0};

(* --- Solve unknowns --- *)
solMatch = FindRoot[
  matchingResiduals[{u0, v0, alpha, beta}] == {0, 0, 0, 0},
  {{u0, initGuess[[1]]}, {v0, initGuess[[2]]}, {alpha, initGuess[[3]]}, {beta, initGuess[[4]]}},
  MaxIterations -> 200,
  WorkingPrecision -> 15
];

{u0Sol, v0Sol, alphaSol, betaSol} = {u0, v0, alpha, beta} /. solMatch;

Print["Solved parameters:"];
Print["u0 = ", u0Sol, ", v0 = ", v0Sol, ", alpha = ", alphaSol, ", beta = ", betaSol];

(* --- Reconstruct left/right solutions *)
leftFinal = leftSolver[u0Sol, v0Sol];
uLeftFunc = leftFinal[[1]];
vLeftFunc = leftFinal[[2]];

uRC = alphaSol*uLeftFunc[c];
vRC = betaSol*vLeftFunc[c];
rightFinal = rightSolver[uRC, vRC];
uRightFunc = rightFinal[[1]];
vRightFunc = rightFinal[[2]];

(* --- Piecewise full-domain solutions *)
uSol[x_] := Piecewise[{{uLeftFunc[x], x <= c}, {uRightFunc[x], x > c}}];
vSol[x_] := Piecewise[{{vLeftFunc[x], x <= c}, {vRightFunc[x], x > c}}];

(* --- Plot solutions --- *)
Plot[
  Evaluate[{uSol[x], vSol[x]}],
  {x, 0, 1},
  PlotLegends -> {"u(x)", "v(x)"},
  AxesLabel -> {"x", "Concentration"},
  PlotTheme -> "Detailed",
  PlotRange -> All
]

(* --- Export data to CSV --- *)
data = Table[{{x, N[uSol[x]], N[vSol[x]]}}, {{x, 0, 1, 0.01}}];
Export["{output_csv_path}", Prepend[data, {{"x", "u", "v"}}]];
Print["CSV written to {output_csv_path}"];
"""
    with open(output_wl_path, "w") as f:
        f.write(wl_template)


if __name__ == "__main__":
    folder_to_check_validity = sys.argv[1]
    write_wl_script(params, output_wl_path, output_csv_path)

"""
Need information about:
1
""" 
#%%
import os
import pickle
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from create_reaction_network import *
from auxiliary_functions_using_standard_library import pickle_load_binary

# Get the absolute path to this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Go up to the project root (one or more levels as needed)
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))

# Access other folders relative to the root
DATA_PATH = os.path.join(PROJECT_ROOT, "data")
DATA_PATH

violacein_example_folder = os.path.join(DATA_PATH, "violacein_2")

reaction_network = pickle_load_binary(os.path.join(violacein_example_folder, ".REACTION_NETWORK_pickle"))
reaction_network
#%%
reaction_network.species["Trp"].as_reactant_in[0]

# %%
reaction_term_arguments_list = []
for species in reaction_network.species:
    reaction_term_arguments_list.append(f"{species.name}_")
reaction_term_arguments_list.append("x_")
reaction_term_arguments = ", ".join(reaction_term_arguments_list)
reaction_term_arguments

#%%

for species in reaction_network.species:
    reaction_terms_string = ""
    for reaction in species.as_reactant_in + species.as_product_in:
        if reaction in species.as_reactant_in: # if acts as reactant, diminishes
            prefactor = "-"
        else:
            prefactor = "+"
        if isinstance(reaction, SpontaneousReaction):
            term = f"{prefactor}{reaction.k} * {reaction.start_species.name}"
        elif isinstance(reaction, EnzymaticReaction):
            term = f"{prefactor}{reaction.k_cat} * {reaction.enzyme} * {reaction.start_species.name} / ({reaction.k_M} + {reaction.start_species.name})"
        reaction_terms_string += term
    
    complete_reaction_string = f"reaction{species.name}[{}]"


"""
k[x_] := 1 + Sin[Pi x]  (* example *)

f[u_, v_, x_] := k[x] * u
g[u_, v_, x_] := -k[x] * v

 -diff*u''[x] + f[u[x], v[x], x] == 0,


k[x_] := Piecewise[{{1, x < 0.5}, {2, x >= 0.5}}]

Piecewise[{{val1, cond1}, {val2, cond2}, …}]


Piecewise[{{1, 0 < x < 0.5}, {2, x >= 0.5}}, 0]
the last value is the default value
"""