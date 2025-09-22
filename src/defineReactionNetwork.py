#%% Imports
import os
import numpy as np
import numpy.typing as npt
from typing import Tuple
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from auxFcts import Ratio, define_ratio_from_string
import pickle

#%%

class PairTuple:
    def __init__(self, element1: float, element2: float ):
        self.pair = (min(element1, element2), max(element1, element2))

#%%
class Participant:
    def __init__(self, name: str) -> None:
        self.name = name

    def __str__(self): # short, e.g. for plotting with networkx
        return f"{self.name}"
    
    def __repr__(self):
        return f"{self.name}"

class Enzyme(Participant):
    def __init__(self, name: str, quantity: float, localization:npt.NDArray[PairTuple]) -> None:
        super().__init__(name)
        self.quantity = quantity
        self.localization = localization

    def __str__(self):
        return super().__str__()
    
    def __repr__(self):
        return f"enzyme {super().__repr__()}"

class Species(Participant):
    def __init__(self, name: str,
                 diffusion_constant: float = 1000e-9 * 60,
                 permeability_constant: float = 1e-9 * 60) -> None:
        super().__init__(name)
        self.diffusion_constant = diffusion_constant
        self.permeability_constant = permeability_constant
        self.as_reactant_in = []
        self.as_product_in = []
        self.first_time_derivative_terms = []
        self.second_time_derivative_terms = []

    def __str__(self):
        return super().__str__()
    
    def __repr__(self):
        return f"species {super().__repr__()}"
    
#%%
class Reaction:
    """ Base class that defines reactions. """
    def __init__(self, start_species: Species, end_species: Species, ratio_endtostart_species: Ratio) -> None:
        self.start_species = start_species
        self.end_species = end_species
        self.ratio_endtostart_species = ratio_endtostart_species
    def __repr__(self):
        return f"reaction from {self.start_species.name} to {self.start_species.name}"

class EnzymaticReaction(Reaction):
    """ Derived class from Reaction that defines enzymatic reactions. """
    def __init__(self, start_species: Species, end_species: Species, ratio_endtostart_species: Ratio, enzyme: Enzyme,
                 k_cat: float, k_M: float, hill: float) -> None:
        super().__init__(start_species, end_species, ratio_endtostart_species)
        self.enzyme = enzyme
        self.k_cat = k_cat
        self.k_M = k_M
        self.hill = hill
        self.k = k_cat/k_M

    def __repr__(self):
        return f"enzymatic {super().__repr__()} catalyzed by {self.enzyme}, with k_cat={self.k_cat}, k_M={self.k_M}, and k = {self.k}"

class SpontaneousReaction(Reaction):
    """ Derived class from Reaction that defines spontaneous reactions. """
    def __init__(self, start_species: Species, end_species: Species, ratio_endtostart_species: Ratio,
                 k: float):
        super().__init__(start_species, end_species, ratio_endtostart_species)
        self.k = k
        

    def __repr__(self):
        return f"spontaneous {super().__repr__()}, with k_forth={self.k_forth} and k_back={self.k_back}"

#%%
class System:
    def __init__(self, radius: float, mesh_n: int) -> None:
        self.radius = radius
        self.mesh_n = mesh_n
        self.network = nx.DiGraph() # initialize the reaction network as a directed graph
        self.reactions = []
        self.species = []
        self.enzymes = []
        

    def add_reaction(self, reaction: Reaction) -> None:
        """Adds one edge that represents the reaction. The attributes of the
        reaction are saved in a dictionary.
        Access through:
            for u, v, edge_data in self.network.edges(data=True):
                print(f"Edge from {u} to {v}; Edge attributes dict: {edge_data}")
        """
        self.network.add_edge(
            reaction.start_species, reaction.end_species,
            **vars(reaction)
        )
        # TODO: when I add a reaction, it should modify the species' attribute first derivative terms
        # Calls add_species
    
    def add_species(self, *species: Species) -> None:
        """Add species to the system. Put any number of them as parameters.
        As attributes of the network node we put the attributes of the object as a dict.
        """
        self.species += species
        self.network.add_nodes_from((s, vars(s)) for s in species)

    def draw_network(self) -> None:
        nx.draw_networkx(self.network)

#%% 
def generateReactionNetwork(
        enzymatic_reactions_csv_file, spontaneous_reactions_csv_file, enzymes_csv_file, species_csv_file):
    # pickle System
    return


"""
1. Extract data from .csv files
2. Define all reaction and participant objects.
3. Add the different reactions to the system object.
    With each reaction, add to the list of the attribute the species object and the enzyme object (in case not already there).
    Update the attribute of the species as_reactant_in, as_product_in
4. Check that when I modify an object defined outside of system, the object in the attribute of system also changes.
5. Construct the differential equations

The question is whether it is passing the object itself or a copy
"""
