import sys
import os
from itertools import count
from auxiliary_functions_using_standard_library import pickle_load_binary, load_json
from auxiliary_functions import dump_json
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme

def build_point_ids_dict(reaction_network: System, num_mesh_points_in_regions: dict) -> dict:
    """Builds a nested dict mapping (region, mesh_point, species) to unique IDs.
    Access values of unique IDs through point_ids_dict[region_idx][mesh_point_idx][species]
    (remider that the key is the species object itself, not its .name attribute)
    reaction_network is of type System (see file src/create_reaction_network.py)
    num_mesh_points_in_regions is a dictionary where the key is an integer (which region)
    and the value is the number of mesh points that are within that region.
    """
    counter = count()  # local counter — resets each time you call the function
    point_ids_dict = {
        region_idx: {
            mesh_point_idx: {
                species: next(counter)
                for species in reaction_network.species
            }
            for mesh_point_idx in range(num_mesh_points_in_regions[region_idx])
        }
        for region_idx in range(len(num_mesh_points_in_regions))
    }
    return point_ids_dict

def build_reverse_point_ids_dict(point_ids_dict) -> dict:
    """ Takes the return dict from build_point_ids_dict and constructs a
    inverse dictionary.
    The key is the index of the node and the value is (region, n, species).
    """
    reverse_point_ids_dict = {
        value: (region_idx, mesh_point_idx, species)
        for region_idx, mesh_points in point_ids_dict.items()
        for mesh_point_idx, species_map in mesh_points.items()
        for species, value in species_map.items()
    }
    return reverse_point_ids_dict

def build_radii_dict(mesh_points_in_regions) -> dict:
    """ Returns a dictionary dict_name[region][n] : radius_to_origin
    """
    radii_dict = {
        region_idx : {
            mesh_point_idx : mesh_points_in_regions[region_idx][mesh_point_idx]
            for mesh_point_idx in range(len(mesh_points_in_regions[region_idx]))
        }
        for region_idx in range(len(mesh_points_in_regions))
    }
    return radii_dict

def build_point_infos_dict(num_mesh_points_in_regions) -> dict:
    """ Gives information on whether the node is within the bulk of the region,
    the left-most node within the region or the right-most node within the region.
    dict_name[region][n] : "i" or "l" or "r", respectively
    """
    point_infos_dict = {
        region_idx : {
            mesh_point_idx : "l" if mesh_point_idx==0 else ("r" if mesh_point_idx==num_mesh_points_in_regions[region_idx]-1 else "i")
            for mesh_point_idx in range(num_mesh_points_in_regions[region_idx])
        }
        for region_idx in range(len(num_mesh_points_in_regions))
    }
    return point_infos_dict

def build_point_neighbor_dict(num_mesh_points_in_regions) -> dict:
    """ Gives a list of (region, n) tuples for each [region][n] which specifies
    the node information of spatial neighbors (and itself).
    dict_name[region][n] : [(region, n-1), (region, n), (region, n+1)] e.g. for
    "i" nodes. for "l" and "r" nodes it depends on whether the node is at an inner
    boundary or at r=0 or r=R. Always goes from left to right.
    """
    neighbors_dict = {}
    point_infos = build_point_infos_dict(num_mesh_points_in_regions)
    for region_idx, mesh_points in point_infos.items():
        for n, kind in mesh_points.items():
            if kind == "i":
                # Internal: previous, self, next
                neighbors_dict[(region_idx, n)] = [
                    (region_idx, n - 1),
                    (region_idx, n),
                    (region_idx, n + 1),
                ]
            elif kind == "l":
                # Left boundary: connect to previous region's rightmost point (if it exists)
                if region_idx > 0:
                    prev_region_last = num_mesh_points_in_regions[region_idx - 1] - 1
                    neighbors_dict[(region_idx, n)] = [
                        (region_idx - 1, prev_region_last),
                    ]
                else:
                    neighbors_dict[(region_idx, n)] = []
                neighbors_dict[(region_idx, n)].append( (region_idx, 0) )
                neighbors_dict[(region_idx, n)].append( (region_idx, 1) )
            elif kind == "r":
                # Right boundary: self’s region last, then 0 and 1 in same region
                last_in_region = num_mesh_points_in_regions[region_idx] - 1
                neighbors_dict[(region_idx, n)] = [
                    (region_idx, last_in_region-1),
                    (region_idx, last_in_region),
                ]
                if region_idx < len(num_mesh_points_in_regions)-1:
                    neighbors_dict[(region_idx, n)].append( (region_idx+1, 0) )
    return neighbors_dict

if __name__ == "__main__":
    FOLDER_TO_SOLVE = sys.argv[1]    
    # Step 0: Load inputs and define global parameters
    REACTION_NETWORK = pickle_load_binary(os.path.join(FOLDER_TO_SOLVE, ".pickled_REACTION_NETWORK"))
    SYSTEM_GEOMETRY_DICT = load_json(os.path.join(FOLDER_TO_SOLVE, ".expanded_system_geometry.json"))
    SOLVER_INPUT = load_json(os.path.join(FOLDER_TO_SOLVE, "solver_input.json"))

    # Step 1: Define all geometry variables
    MESH_POINTS_IN_REGIONS = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["MESH_POINTS_IN_REGIONS"]
    NUM_MESH_POINTS_IN_REGIONS = SYSTEM_GEOMETRY_DICT["GEOMETRY_CONFIG"]["NUM_MESH_POINTS_IN_REGIONS"]
    
    # Step 2: Define structures to access geometry information
    POINT_IDS = build_point_ids_dict(REACTION_NETWORK, NUM_MESH_POINTS_IN_REGIONS)
    REVERSE_POINT_IDS = build_reverse_point_ids_dict(POINT_IDS)
    RADII = build_radii_dict(MESH_POINTS_IN_REGIONS)
    DELTA_R = RADII[0][1]-RADII[0][0] # the different points within a region are equally spaced
    NUM_POINTS = len(REVERSE_POINT_IDS) # each point saves the concentration for one species at one node
    POINT_INFOS = build_point_infos_dict(NUM_MESH_POINTS_IN_REGIONS)
    NEIGHBORS = build_point_neighbor_dict(NUM_MESH_POINTS_IN_REGIONS)
    
    # Check that each region has at least 3 points
    for region, radii in RADII.items():
        if len(radii)<3:
            raise ValueError(f"Region {region} has less than 3 points, such that the diffusion term does not work.")

    # Save dictionary in .json file for readability
    dict_to_dump = {
        "POINT_IDS": POINT_IDS,
        "REVERSE_POINT_IDS": REVERSE_POINT_IDS,
        "RADII": RADII,
        "DELTA_R": DELTA_R,
        "NUM_POINTS": NUM_POINTS,
        "POINT_INFOS": POINT_INFOS,
        "NEIGHBORS": NEIGHBORS,

    }
    dump_json(FOLDER_TO_SOLVE, ".expanded_system_mesh", dict_to_dump)


    