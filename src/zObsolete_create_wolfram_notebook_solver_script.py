#%%
import os
import pickle
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from create_reaction_network import Participant, Enzyme, Species, Reaction, EnzymaticReaction, SpontaneousReaction, Collection, System
from auxiliary_functions_using_standard_library import pickle_load_binary

# Get the absolute path to this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Go up to the project root (one or more levels as needed)
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))

# Access other folders relative to the root
DATA_PATH = os.path.join(PROJECT_ROOT, "data")

violacein_example_folder = os.path.join(DATA_PATH, "violacein_2")

reaction_network = pickle_load_binary(os.path.join(violacein_example_folder, ".REACTION_NETWORK_pickle"))

#%%
def write_wl_script(reaction_network_object, geometry_params, output_wl_path, output_csv_path):
    # ClearAll at the beginning
    clearAll_template = f"""
(* --- Clear all variables --- *)    
ClearAll["Global`*"];
"""
    # define parameters
    parameters_template = f"""
(* --- Parameters --- *)
R = {geometry_params["R"]};
r = {geometry_params["r"]};
"""
    # define species string
    species_string = ", ".join([species.name for species in reaction_network_object.species])
    unknown_concentration_left_species_string = ", ".join([species.name+"0" for species in reaction_network_object.species])
    unknown_concentration_right_species_string = ", ".join([species.name+"R" for species in reaction_network_object.species])
    # define arguments for the definitions of the reactions
    reaction_term_arguments_list = []
    for species in reaction_network_object.species:
        reaction_term_arguments_list.append(f"{species.name}_")
    reaction_term_arguments_list.append("x_")
    reaction_term_arguments = ", ".join(reaction_term_arguments_list)
    # define arguments for accessing the reactions within the differential equation
    reaction_term_arguments_for_diff_eq_list = []
    for species in reaction_network_object.species:
        reaction_term_arguments_for_diff_eq_list.append(f"{species.name}[x]")
    reaction_term_arguments_for_diff_eq_list.append("x")
    reaction_term_arguments_for_diff_eq = ", ".join(reaction_term_arguments_for_diff_eq_list)
    # define reactions
    for species in reaction_network_object.species:
        reaction_terms_string = ""
        for reaction in species.as_reactant_in + species.as_product_in:
            if reaction in species.as_reactant_in: # if acts as reactant, diminishes
                prefactor = "-"
            else:
                prefactor = "+"
            if isinstance(reaction, SpontaneousReaction):
                term = f"{prefactor}{reaction.k} * {reaction.start_species.name}"
            elif isinstance(reaction, EnzymaticReaction):
                term = f"{prefactor}{reaction.k_cat} * {reaction.enzyme}[x] * {reaction.start_species.name} / ({reaction.k_M} + {reaction.start_species.name})"
            reaction_terms_string += term
        species.reaction_string = f"reaction{species.name}[{reaction_term_arguments}] := {reaction_terms_string};"
    # write reactions in a multiline string
    reaction_multiline_string_list = []
    for species in reaction_network_object.species:
        reaction_multiline_string_list.append(species.reaction_string)
    reaction_multiline_string_template = "(* --- Reaction functions --- *)\n" + "\n".join(reaction_multiline_string_list)
    # define location of enzymes
    for enzyme in reaction_network_object.enzymes:
        individual_localization_info = []
        for localization in enzyme.localization:
            tuple_info = f"{{ {enzyme.concentration}, {localization.minMaxLoc[0]} < x < {localization.minMaxLoc[1]} }}"
            individual_localization_info.append(tuple_info)
        enzyme.concentration_info = f"{enzyme.name}[x_] := Piecewise[{{{", ".join(individual_localization_info)}}}]"
    # write location of enzymes in a multiline string
    enzyme_location_multiline_string_list = []
    for enzyme in reaction_network_object.enzymes:
        enzyme_location_multiline_string_list.append(enzyme.concentration_info)
    enzyme_location_multiline_string_template = "(* --- Enzyme concentration and location --- *)\n" + "\n".join(enzyme_location_multiline_string_list)
    # define information about system
    for species in reaction_network_object.species:
        # define differential equations
        species.steady_state_equation = f"-{species.diffusion_constant}*{species.name}''[x] + reaction{species.name}[{reaction_term_arguments_for_diff_eq}] == 0,"
        # define unknowns of the concentrations on the left
        species.unknown_concentration_left = f"{species.name}[0] == {species.name}0,"
        # define Neumann boundary condition on the left
        species.Neumann_condition_left = f"{species.name}'[0] == 0,"
        # define unknowns of the concentration on the right
        species.unknown_concentration_right = f"{species.name}[R] == {species.name}R,"
        # define flux boundary condition on the right
        species.flux_condition_right = f"{species.name}'[R] == {species.permeability_constant}/{species.diffusion_constant} * ({species.external_concentration} - {species.name}[R]),"
    # make template for left domain solver: write equations
    left_domain_info_list = []
    for species in reaction_network_object.species:
        left_domain_info_list.append(species.steady_state_equation)
        left_domain_info_list.append(species.unknown_concentration_left)
        left_domain_info_list.append(species.Neumann_condition_left)
    left_domain_info_multiline_string_template = "\n".join(left_domain_info_list)
    # make template for left domain solver: put everything together
    left_domain_template = f"""
(* --- Left domain solver: x in [0,r] *)
leftSolver = ParametricNDSolveValue[
{{
{left_domain_info_multiline_string_template}
}},
{{ {species_string} }},
{{ x, 0, r }},
{{ {unknown_concentration_left_species_string}}}
],
"""    
    # make template for right domain solver: write equations
    right_domain_info_list = []
    for species in reaction_network_object.species:
        right_domain_info_list.append(species.steady_state_equation)
        right_domain_info_list.append(species.unknown_concentration_right)
        right_domain_info_list.append(species.flux_condition_right)
    right_domain_info_multiline_string_template = "\n".join(right_domain_info_list)
    # make template for left domain solver: put everything together
    right_domain_template = f"""
(* --- Right domain solver: x in [r,R] *)
rightSolver = ParametricNDSolveValue[
{{
{right_domain_info_multiline_string_template}
}},
{{ {species_string} }},
{{ x, r, R }},
{{ {unknown_concentration_right_species_string}}}
],
"""



    wl_content = "\n".join([
    clearAll_template,
    parameters_template,
    reaction_multiline_string_template,
    enzyme_location_multiline_string_template,
    left_domain_template,
    right_domain_template
    ])
    with open(output_wl_path, "w") as f:
        f.write(wl_content)

#%%
write_wl_script(reaction_network, {"r": 0.5, "R": 1}, "test.wl", 
                "")




#%%
def write_wl_script(params, output_wl_path, output_csv_path):
    clearAll_template = f"ClearAll["Global`*"];"


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
