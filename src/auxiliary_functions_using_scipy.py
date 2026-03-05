from scipy.sparse import csr_matrix, coo_matrix
import numpy as np
import os
from auxiliary_functions import dump_json


def save_matrix_as_sparse_txt(matrix: np.ndarray, filepath: str):
    """
    Save a 2D NumPy array as a sparse matrix in .txt format (row, col, value).
    Only nonzero entries are stored.
    """
    # Convert to CSR sparse matrix
    sparse_mat = csr_matrix(matrix)
    coo = sparse_mat.tocoo()
    # Stack row, col, data
    data = np.column_stack((coo.row, coo.col, coo.data))
    # Save to .txt
    np.savetxt(filepath+".txt", data, fmt=["%d", "%d", "%.15e"],
               header="row\tcol\tvalue", delimiter="\t", comments='')

def save_newton_iteration_data(
    folder_to_save_in, filename, J_to_save, F_to_save,
    species_concentrations_to_save, du_to_save, variables_to_save_dictionary):
    """Function to save data as needed, based on the dictionary variables_to_save_dictionary.
    """
    if variables_to_save_dictionary["save_F_vector"]:
        np.savetxt(os.path.join(folder_to_save_in, f"{filename}_F.txt"), F_to_save, fmt="%.15e", delimiter="\n")
    if variables_to_save_dictionary["save_F_vector_norm"]:
        dump_json(folder_to_save_in, f"{filename}_F_vector_norm",
            np.linalg.norm(F_to_save))
    if variables_to_save_dictionary["save_J_matrix"]:
        save_matrix_as_sparse_txt(J_to_save, os.path.join(folder_to_save_in, f"{filename}_J_matrix"))
    if variables_to_save_dictionary["save_du_vector"]:
        np.savetxt(os.path.join(folder_to_save_in, f"{filename}_du_vector.txt"), du_to_save, fmt="%.15e", delimiter="\n")
    if variables_to_save_dictionary["save_du_vector_max"]:
        dump_json(folder_to_save_in, f"{filename}_du_vector_max",
            max(du_to_save))
    if variables_to_save_dictionary["save_concentrations"]:    
        dump_json(folder_to_save_in, f"{filename}_concentrations",
            species_concentrations_to_save)