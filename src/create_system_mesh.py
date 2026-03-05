import sys
import os
from itertools import count
from auxiliary_functions_using_standard_library import pickle_load_binary, load_json
from auxiliary_functions import dump_json, read_yaml_file
from create_reaction_network import System, Collection, EnzymaticReaction, Species, SpontaneousReaction, Enzyme

def interpolate_midpoints(values, interpolation_times):
    """
    Inserts midpoints between neighboring values,
    doing this a number interpolation_times of times
    and returns the new list.
    values : list of numbers (ordered)
    """
    if interpolation_times < 0:
        raise ValueError("times must be non-negative")
    result = list(values)
    for _ in range(interpolation_times):
        new_list = []
        for i in range(len(result) - 1):
            new_list.append(result[i])
            midpoint = (result[i] + result[i + 1]) / 2
            new_list.append(midpoint)
        new_list.append(result[-1])
        result = new_list
    return result

def define_mesh_points(
    system_geometry_dict,
    mesh_points_duplication_times
):
    """
    system_geometry_dict has the baseline geometry
    """
    baseline_mesh_points = system_geometry_dict["geometry_config"]["baseline_mesh_points"]
    modified_mesh_points = interpolate_midpoints(baseline_mesh_points, mesh_points_duplication_times)
    system_geometry_dict["geometry_config"]["modified_mesh_points"] = modified_mesh_points
    boundary_radii = system_geometry_dict["geometry_config"]["boundary_radii"]
    
    mesh_points_in_regions = {
        region_idx : [mesh_point for mesh_point in modified_mesh_points if boundary_radii[region_idx]<=mesh_point<=boundary_radii[region_idx+1]]
        for region_idx in range(system_geometry_dict["geometry_config"]["num_regions"])
    }
    system_geometry_dict["geometry_config"]["mesh_points_in_regions"]  = mesh_points_in_regions
    num_mesh_points_in_regions = {
        region_idx: len(mesh_points_in_region)
        for region_idx, mesh_points_in_region in mesh_points_in_regions.items()
    }
    system_geometry_dict["geometry_config"]["num_mesh_points_in_regions"]  = num_mesh_points_in_regions

    return system_geometry_dict

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

def build_system_mesh(system_geometry_dict, reaction_network, mesh_points_duplication_times):
    # Define mesh points for given number of interpolation iterations
    system_geometry_nested_dict = define_mesh_points(system_geometry_dict, mesh_points_duplication_times)
    
    # Define aliases for easier access
    num_mesh_points_in_regions = system_geometry_nested_dict["geometry_config"]["num_mesh_points_in_regions"]
    mesh_points_in_regions = system_geometry_nested_dict["geometry_config"]["mesh_points_in_regions"]
    
    # Define structures to access geometry information
    point_ids = build_point_ids_dict(reaction_network, num_mesh_points_in_regions)
    reverse_point_ids = build_reverse_point_ids_dict(point_ids)
    radii = build_radii_dict(mesh_points_in_regions)
    delta_r = radii[0][1]-radii[0][0] # the different points within a region are equally spaced
    num_points = len(reverse_point_ids) # each point saves the concentration for one species at one node
    point_infos = build_point_infos_dict(num_mesh_points_in_regions)
    neighbors = build_point_neighbor_dict(num_mesh_points_in_regions)
    
    # Check that each region has at least 3 points
    for region, region_radii in radii.items():
        if len(region_radii)<3:
            raise ValueError(f"Region {region} has less than 3 points, such that the diffusion term does not work.")

    # Save dictionary in .json file for readability
    expanded_system_mesh = {
        "point_ids": point_ids,
        "reverse_point_ids": reverse_point_ids,
        "radii": radii,
        "delta_r": delta_r,
        "num_points": num_points,
        "point_infos": point_infos,
        "neighbors": neighbors
    }
    return system_geometry_nested_dict, expanded_system_mesh


    