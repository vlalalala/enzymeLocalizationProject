import os
import numpy as np
import ast
from auxiliary_functions import save_matrix_as_sparse_txt, dump_json

def get_species_concentrations_from_json_file(
        imported_concentrations_dict_from_json, species_lookup_dict) -> dict:
    """Make keys that got converted to strings instead of integers into integers again
    """
    dict_in_correct_concentrations_format = {
            int(region_idx): {
                int(mesh_point_idx): {
                    species_lookup_dict[species_name]: data
                    for species_name, data in mesh_point_info.items()}
                for mesh_point_idx, mesh_point_info in region_info.items()}
            for region_idx, region_info in imported_concentrations_dict_from_json.items()
        }
    return dict_in_correct_concentrations_format

def get_correct_point_ids_dict(imported_point_ids_dict, species_lookup_dict
    ):
    """Works exactly the same as get_species_concentrations_from_json_file, just that every
    [region][n][species] maps to the point_id and not to the concentration
    """
    dict_in_correct_concentrations_format = get_species_concentrations_from_json_file(imported_point_ids_dict, species_lookup_dict)
    return dict_in_correct_concentrations_format

def get_correct_reverse_point_ids_dict(imported_reverse_point_ids_dict,
        species_lookup_dict
    ):
    dict_inc_correct_format = {
        int(k): [*v[:-1], species_lookup_dict[v[-1]]]
        for k, v in imported_reverse_point_ids_dict.items()
    }
    return dict_inc_correct_format

def get_correct_neighbors_dict(imported_neighbors_dict):
    """From the json file, the neighbors tuple gets saved as a list. Convert back.
    """
    new_dict = {}
    for k, v in imported_neighbors_dict.items():
        # Parse the string "[0, 0]" into a real list [0, 0]
        key_as_list = ast.literal_eval(k)
        # Convert it to a tuple
        key_as_tuple = tuple(key_as_list)
        new_dict[key_as_tuple] = v
    return new_dict

def save_newton_iteration_data(
    folder_to_save_in, iter_string, J_to_save, F_to_save,
    species_concentrations_to_save, du_to_save, variables_to_save_dictionary):
    """Function to save data as needed, based on the dictionary variables_to_save_dictionary.
    """
    if variables_to_save_dictionary["save_F_vector"]:
        np.savetxt(os.path.join(folder_to_save_in, f".iteration_nr_{iter_string}_F.txt"), F_to_save, fmt="%.15e", delimiter="\n")
    if variables_to_save_dictionary["save_F_vector_norm"]:
        dump_json(folder_to_save_in, f".iteration_nr_{iter_string}_F_vector_norm",
            np.linalg.norm(F_to_save))
    if variables_to_save_dictionary["save_J_matrix"]:
        save_matrix_as_sparse_txt(J_to_save, os.path.join(folder_to_save_in, f".iteration_nr_{iter_string}_J_matrix"))
    if variables_to_save_dictionary["save_du_vector"]:
        np.savetxt(os.path.join(folder_to_save_in, f".iteration_nr_{iter_string}_du_vector.txt"), du_to_save, fmt="%.15e", delimiter="\n")
    if variables_to_save_dictionary["save_du_vector_max"]:
        dump_json(folder_to_save_in, f".iteration_nr_{iter_string}_du_vector_max",
            max(du_to_save))
    if variables_to_save_dictionary["save_concentrations"]:    
        dump_json(folder_to_save_in, f".iteration_nr_{iter_string}_concentrations",
            species_concentrations_to_save)