import os
import numpy as np
from auxiliary_functions_using_standard_library import dump_json
from auxiliary_functions import save_matrix_as_sparse_txt

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