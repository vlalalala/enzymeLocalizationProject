import sys
import os
import numpy as np
from scipy.integrate import solve_bvp
from auxiliary_functions_using_standard_library import pickle_load_binary, pickle_dump_binary
from auxiliary_functions import LocationTuple
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme

def michaelis_menten_term(k_cat, k_M, c_enzyme, c_substrate, hill=1):
    return k_cat * c_enzyme * c_substrate**hill / (k_M**hill + c_substrate**hill)

def spontaneous_reaction_term(k, c_reactant):
    return k*c_reactant


def reaction_diffusion_system(r_mesh, variable_values_2d_array):
    """ Returns right-hand side.
    Return array with shape (n,m), n = number of variables, m = number of nodes
    r_mesh has shape (m,)
    value has shape (n, m) # for each substance n at each node m a specific value (float)

    THE VALUES IN R_MESH CAN CHANGE, so the enzymes concentration per mesh
    point have to be reevaluated each time according to the r_mesh being used...
    """
    # Get the information from the 2d array for each species
    species_conc = { # concentration
        species.name : variable_values_2d_array[species_variable_indices[species]["prim"]]
        for species in reaction_network.species
    }
    species_conc_der = { # 1st derivative of concentration
        species.name : variable_values_2d_array[species_variable_indices[species]["1stDer"]]
        for species in reaction_network.species
    }

    # Get the concentration of enzymes at the r_mesh nodes
    enzyme_conc = {enzyme: np.zeros(len(r_mesh))
        for enzyme in reaction_network.enzymes
    }
    for enzyme in reaction_network.enzymes:
        for locTuple in enzyme.localization:
            r1 = locTuple.minMaxLoc[0] * max(r_mesh)
            r2 = locTuple.minMaxLoc[1] * max(r_mesh)
            for index, radius in enumerate(r_mesh):
                if r1 <= radius <= r2:
                    enzyme_conc[enzyme][index] = enzyme.concentration

    system_right_hand_side = { # includes no reaction terms here
        species: {"prim": species_conc_der[species.name],
                  "1stDer": -2/r_mesh * species_conc_der[species.name]}
        for species in reaction_network.species
    }

    # Calculate the reaction terms
    for species in reaction_network.species:
        reaction_term_sum = np.zeros(len(r_mesh))
        for reaction in species.as_reactant_in + species.as_product_in:
            if isinstance(reaction, SpontaneousReaction):
                term = reaction.k * species_conc[reaction.start_species.name] #################################
            elif isinstance(reaction, EnzymaticReaction):
                term = reaction.k_cat * enzyme_conc[reaction.enzyme] * species_conc[reaction.start_species.name] / (reaction.k_M + species_conc[reaction.start_species.name])
            if reaction in species.as_reactant_in: # if acts as reactant, diminishes
                term *= -1
            reaction_term_sum += term
        # Add reaction term to diffusion term
        system_right_hand_side[species]["1stDer"] += -1/species.diffusion_constant * reaction_term_sum

    # Create the right-hand side 2d array
    values_2d_array = np.zeros((2*len(reaction_network.species), len(r_mesh)))
    for species in reaction_network.species:
        for variable_type in species_variable_indices[species].keys():
            # writes the species_variable_indices[species][variable_type]th row of the 2d array
            values_2d_array[species_variable_indices[species][variable_type]] = system_right_hand_side[species][variable_type]

    return values_2d_array

def boundary_conditions(y_a, y_b):
    """ At origin a we have reflexion, so the derivative variable is 0.
    At r = R a flux is given.
    y_a and y_b are the vectors with shape (n, ) (n=2*num_species)
    def bc(ya, yb):
        return np.array([
            ya[0] - 1,  # y0(a) = 1
            yb[1] - 0   # y1(b) = 0
        ])
    """
    # origin
    origin_conditions = [
        y_a[species_variable_indices[species]["1stDer"]] - 0
        for species in reaction_network.species
    ]
    # edge
    border_conditions = [(
        y_b[species_variable_indices[species]["1stDer"]]
        - species.permeability_constant / species.diffusion_constant
        * (species.external_concentration - y_b[species_variable_indices[species]["prim"]]))
        for species in reaction_network.species
    ]
    conditions = np.array(origin_conditions + border_conditions)
    return conditions

if __name__ == "__main__":
    # Load all the information
    folder_to_solve = sys.argv[1]
    reaction_network = pickle_load_binary(os.path.join(folder_to_solve, ".REACTION_NETWORK_pickle"))
    system_geometry_dict = pickle_load_binary(os.path.join(folder_to_solve, ".SYSTEM_GEOMETRY_pickle"))

    # Step 1: Define all geometry variables
    radius = system_geometry_dict["GEOMETRY_CONFIG"]["radius"]
    num_mesh_points = system_geometry_dict["GEOMETRY_CONFIG"]["num_mesh_points"]
    initial_r_mesh = np.linspace(radius*1e-4, radius, num = num_mesh_points) #nm # doesn't work within singularity

    # Step 2: Do preliminary work to make reaction diffusion system solvable
    """We need to define 2*n variables, where n is the number of species, since we
    have to convert the system with second derivatives to one with first derivatives
    with double the variables
    """
    species_variable_indices = { # each variable is mapped to an integer
        species: {"prim": species_idx*2,
                "1stDer": species_idx*2 + 1}
        for species_idx, species in enumerate(reaction_network.species)
    }
    # Set some initial values within the nodes for each variable 
    species_variables_inital_values_dict = {
        species: {key: np.zeros(num_mesh_points) for key in inner_dict}
        for species, inner_dict in species_variable_indices.items()
    }
    # Initial guesses can be modified like this, if needed
    species_variables_inital_values_dict[reaction_network.species["Trp"]]["prim"] = np.zeros(num_mesh_points)
    # Create the 2d array map
    initial_values_2d_array = np.zeros((2*len(reaction_network.species), num_mesh_points))
    for species in reaction_network.species:
        for variable_type in species_variable_indices[species].keys():
            # writes the species_variable_indices[species][variable_type]th row of the 2d array
            initial_values_2d_array[species_variable_indices[species][variable_type]] = species_variables_inital_values_dict[species][variable_type]

    # Step 3: run solver
    print(f"Solving bvp in {folder_to_solve}")
    bvp_result = solve_bvp(reaction_diffusion_system, boundary_conditions,
                    initial_r_mesh, initial_values_2d_array,
                    max_nodes=1000)
    print(f"Solved bvp in {folder_to_solve}")
    
    # Step 4: pickle
    pickle_dump_binary(
        os.path.join(folder_to_solve, ".BOUNDARY_VALUE_PROBLEM_RESULT_pickle"),
        bvp_result
    )
    pickle_dump_binary(
        os.path.join(folder_to_solve, ".BOUNDARY_VALUE_PROBLEM_VARIABLE_INDICES_pickle"),
        species_variable_indices
    )

    