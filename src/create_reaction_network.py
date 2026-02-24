#%% Imports
import os
import pickle
import numpy as np
import numpy.typing as npt
from typing import Tuple
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from auxiliary_functions_using_standard_library import load_json, pickle_load_binary, pickle_dump_binary
from auxiliary_functions import Ratio, define_ratio_from_string, print_network_info

import sys
# path of enzymeLocalizationProject
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(project_root) # to be able to access definitions of classes in auxFcts.py

#%%
class Participant:
    def __init__(self, name: str) -> None:
        self.name = name

    def __str__(self): # short, e.g. for plotting with networkx
        return f"{self.name}"
    
    def __repr__(self):
        return f"{self.name}"

class Enzyme(Participant):
    def __init__(self, name: str, quantity: float, regions: list) -> None:
        super().__init__(name)
        self.quantity = quantity
        self.regions = regions

    def __str__(self):
        return super().__str__()
    
    def __repr__(self):
        return f"enzyme {super().__repr__()}"

class Species(Participant):
    def __init__(self, name: str,
                 diffusion_constant: float,
                 external_concentration:float,
                 **kwargs) -> None:
        super().__init__(name)
        self.diffusion_constant = diffusion_constant
        self.external_concentration = external_concentration
        for key, value in kwargs.items():
            setattr(self, key, value)
        self.as_reactant_in = []
        self.as_product_in = []

    def __str__(self):
        return super().__str__()
    
    def __repr__(self):
        return f"species {super().__repr__()}"
    
    def __eq__(self, other): # to be able to use object as key in dictionary, when the dictionary has been pickled and unpicklend
        return isinstance(other, Species) and self.name == other.name

    def __hash__(self):
        return hash(self.name)

class Reaction:
    """ Base class that defines reactions. """
    def __init__(self, start_species: Species, end_species: Species, ratio_endtostart_species: Ratio) -> None:
        self.start_species = start_species
        self.end_species = end_species
        self.ratio_endtostart_species = ratio_endtostart_species
        self.name = ""

    def __repr__(self):
        return f"reaction from {self.start_species.name} to {self.end_species.name}"

class EnzymaticReaction(Reaction):
    """ Derived class from Reaction that defines enzymatic reactions. """
    def __init__(self, start_species, end_species, ratio_endtostart_species: Ratio, enzyme,
                 k_cat: float, k_M: float, hill: float) -> None:
        super().__init__(start_species, end_species, ratio_endtostart_species)
        self.enzyme = enzyme
        self.k_cat = k_cat
        self.k_M = k_M
        self.hill = hill
        self.k = k_cat/k_M

    def __str__(self):
        return f"enzymatic {super().__repr__()} catalyzed by {self.enzyme}"

    def __repr__(self):
        return f"{self.__str__()}, with k_cat={self.k_cat}, k_M={self.k_M}, and k = {self.k}"

class SpontaneousReaction(Reaction):
    """ Derived class from Reaction that defines spontaneous reactions. """
    def __init__(self, start_species: Species, end_species: Species, ratio_endtostart_species: Ratio,
                 k: float):
        super().__init__(start_species, end_species, ratio_endtostart_species)
        self.k = k
    
    def __str__(self):
        return f"spontaneous {super().__repr__()}"

    def __repr__(self):
        return f"{self.__str__()}, with k={self.k}"

#%%
class Collection:
    """Defined in order to be able to get object from list by name and not (solely) by index
    (overfill [] operator)
    """
    def __init__(self, network_object_list):
        self._network_object_list = network_object_list  # list of objects
        self._by_name = {network_object.name: network_object for network_object in network_object_list}  # dict for fast lookup

    def __getitem__(self, key):
        if isinstance(key, int):
            return self._network_object_list[key]
        elif isinstance(key, str):
            return self._by_name[key]
        else:
            raise TypeError("Key must be int (index) or str (network object name)")

    def __iter__(self):
        return iter(self._network_object_list)

    def __len__(self):
        return len(self._network_object_list)

    def append(self, network_object):
        self._network_object_list.append(network_object)
        self._by_name[network_object.name] = network_object
    
    def __add__(self, other):
        if not isinstance(other, Collection):
            return NotImplemented
        combined_list = self._network_object_list + other._network_object_list
        return Collection(combined_list)

    def __repr__(self):
        return repr(self._network_object_list)

#%%
class System:
    def __init__(self) -> None:
        self.enzymatic_reactions = []
        self.spontaneous_reactions = []
        self.species = []
        self.enzymes = []
        self.network = nx.DiGraph() # initialize the reaction network as a directed graph

    def inform_species_about_reactions(self):
        for reaction in self.enzymatic_reactions + self.spontaneous_reactions:
            reaction.start_species.as_reactant_in.append(reaction)
            reaction.end_species.as_product_in.append(reaction)

    def fill_network(self):
        """Adds all necessary information inside the networkx graph."""
        self.network.add_nodes_from((s, vars(s)) for s in self.species)
        for reaction in self.enzymatic_reactions + self.spontaneous_reactions:
            self.network.add_edge(reaction.start_species, reaction.end_species, **vars(reaction))

    def create_simplified_network(self):
        """Only uses/has the most basic information that must be plotted."""
        simple_G = nx.DiGraph()
        # Add nodes with simple names
        for node in self.network.nodes():
            simple_G.add_node(node.name)
        # Add edges using names
        for start_species, end_species, reaction_data in self.network.edges(data=True):
            start_species_name = start_species.name
            end_species_name = end_species.name
            # Extract relevant info from edge attributes
            enzyme = reaction_data.get("enzyme", None)
            weight = None
            label = None
            if enzyme:
                weight = reaction_data.get("k_cat") # raises an error if k_cat is not found
                label = getattr(enzyme, "name") # raises an error if name does not exist
            else:
                weight = reaction_data.get("k")

            simple_G.add_edge(start_species_name, end_species_name, weight = weight)
            if label:
                simple_G[start_species_name][end_species_name]["label"] = label
        
        return simple_G

    def draw_network(self, case_folder) -> None:
        # Uses simplified graph to plot
        plotting_network = self.create_simplified_network()
        # Define plot size
        fig, ax = plt.subplots(figsize=(5, 5))  # Width x Height in inches
        pos = nx.nx_agraph.graphviz_layout(plotting_network, prog="dot")
        plotting_network.graph["graph"] = {"rankdir": "LR"}
        #if nx.check_planarity(plotting_network)[0]:
        #    pos = nx.planar_layout(plotting_network)
        #else:
        #pos = nx.spring_layout(plotting_network, k=2.0, iterations=200, weight=None)
        #pos = nx.kamada_kawai_layout(plotting_network, weight=None, scale = 100.0) # other layout options: shell, circular, kamada_kawai...
        #pos = nx.spring_layout(plotting_network, k=2.0, iterations=100, weight=None, scale=6.0)
        # Draw nodes
        nx.draw_networkx_nodes(
            plotting_network,
            pos,
            ax = ax,
            node_size=1000,        # Size of the nodes. area of node marker in points squared
            node_color="skyblue",
            edgecolors="black"    # Border color around nodes
        )
        # Draw node labels
        nx.draw_networkx_labels(
            plotting_network,
            pos,
            ax = ax,
            font_size=12,
            font_weight="bold"
        )
        # Draw edges
        nx.draw_networkx_edges(
            plotting_network, pos,
            ax=ax,
            arrows=True,
            arrowstyle='-|>',
            arrowsize=20,
            connectionstyle='arc3,rad=0.0',
            min_source_margin=15,  # shifts arrow start away from node center
            min_target_margin=15   # shifts arrow end before hitting node border
        )
        # Draw edge labels (like reaction names or weights)
        edge_labels = nx.get_edge_attributes(plotting_network, 'label')
        nx.draw_networkx_edge_labels(
            plotting_network,
            pos,
            ax = ax,
            edge_labels=edge_labels,
            font_size=10,
            label_pos=0.5  # Position along the edge: 0=start, 1=end, 0.5=middle
        )
        # Turn off axes
        ax.set_axis_off()
        # Show the plot
        fig.tight_layout()
        fig.savefig(os.path.join(case_folder,"reaction_network_graph.png"), dpi = 300, bbox_inches='tight')
        plt.close(fig)

def create_reaction_network(case_folder, csv_file_names):
    """
    """
    # Step 0: initialize object where all the information will be saved
    system = System()

    # Step 1: Import dataframes
    dataframes = {}
    for csv_file_name in csv_file_names:
        pickle_path = os.path.join(case_folder, f".pickled_dataframe_{csv_file_name}")
        dataframes[csv_file_name] = pickle_load_binary(pickle_path)
    
    # Step 2: read rows for each dataframe and create an object for each; save each object in System
    enzymes = [Enzyme(**row) for _, row in dataframes["enzymes"].iterrows()]
    system.enzymes = Collection(enzymes)
    species = [Species(**row) for _, row in dataframes["species"].iterrows()]
    system.species = Collection(species)
    # For the reactions, creates the .name attribute
    enzymatic_reactions = []
    for _, row in dataframes["enzymatic_reactions"].iterrows():
        enzymatic_reaction = EnzymaticReaction(**row)
        enzymatic_reaction.name = f"{enzymatic_reaction.start_species}->{enzymatic_reaction.end_species} {enzymatic_reaction.enzyme}"
        enzymatic_reactions.append(enzymatic_reaction)
    system.enzymatic_reactions = Collection(enzymatic_reactions)
    
    spontaneous_reactions = []
    for _, row in dataframes["spontaneous_reactions"].iterrows():
        spontaneous_reaction = SpontaneousReaction(**row)
        spontaneous_reaction.name = f"{spontaneous_reaction.start_species}->{spontaneous_reaction.end_species}"
        spontaneous_reactions.append(spontaneous_reaction)
    system.spontaneous_reactions = Collection(spontaneous_reactions)

    # Step 3: convert each enzymatic_reaction.enzyme, reaction.start_species, reaction.end_species
    # from a string to the object itself
    for enzymatic_reaction in system.enzymatic_reactions:
        # find the enzyme object that has a name that matches the enzyme (name) in the reaction
        enzyme_object = next((enzyme for enzyme in system.enzymes
                              if enzyme.name == enzymatic_reaction.enzyme), None)
        enzymatic_reaction.enzyme = enzyme_object

    for reaction in system.enzymatic_reactions + system.spontaneous_reactions:
        start_species_object = next((species for species in system.species
                                     if species.name == reaction.start_species), None)
        end_species_object = next((species for species in system.species
                                     if species.name == reaction.end_species), None)
        reaction.start_species = start_species_object
        reaction.end_species = end_species_object  

    # Step 3: fill as_reactant_in and as_product_in for all species objects
    system.inform_species_about_reactions()

    # Step 4: fill the networkx graph
    system.fill_network()

    # Step 5: draw the network graph and save a png file of it
    #try:
    #    system.draw_network(case_folder)
    #except:
    #    print("Unable to draw network.")
    pickle_dump_binary(os.path.join(case_folder, ".pickled_reaction_network"), system)


if __name__ == "__main__":
    folder_to_check_validity = sys.argv[1]
    reaction_network_info_file_names = load_json("src/_template_reaction_network.json").keys()
    create_reaction_network(folder_to_check_validity, reaction_network_info_file_names)
