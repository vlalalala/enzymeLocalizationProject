import os
import json
import numpy as np
from math import gcd
import pandas as pd
import ast
from scipy.sparse import csr_matrix, coo_matrix
import csv

def dump_json(dump_directory: str, file_basename: str, dict_to_dump: dict):
    """
    Dump a dictionary (possibly nested) as JSON, converting Species objects to their .name.
    Handles Species in keys, values, lists, sets, tuples, etc.
    Set file_basename without the .json ending.
    """
    
    def convert_to_saveable_type(obj):
        # Convert Species objects (keys or values)
        if hasattr(obj, "__class__") and obj.__class__.__name__ == "Species":
            return obj.name

        # Handle NumPy arrays and scalars
        elif isinstance(obj, np.ndarray):
            return obj.tolist()  # convert ndarray → list
        elif isinstance(obj, (np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, (np.int32, np.int64)):
            return int(obj)
        
        # Handle dictionaries (convert both keys and values)
        elif isinstance(obj, dict):
            return {str(convert_to_saveable_type(k)): convert_to_saveable_type(v) for k, v in obj.items()}

        # Handle lists, tuples, sets
        elif isinstance(obj, (list, tuple, set)):
            return [convert_to_saveable_type(i) for i in obj]

        # Base case: leave unchanged
        else:
            return obj

    converted = convert_to_saveable_type(dict_to_dump)
    os.makedirs(dump_directory, exist_ok=True)
    path = os.path.join(dump_directory, f"{file_basename}.json")
    with open(path, "w") as f:
        json.dump(converted, f, indent=4)

class Ratio:
    """Auxiliary class. Defines a ratio of a:b """
    def __init__(self, a: int, b: int) -> None:
        if b == 0:
            raise ValueError("Denominator of a ratio cannot be zero")
        g = gcd(a, b) # get greatest common denominator to have simplified ratios
        self.numerator = a // g
        self.denominator = b // g

    def __repr__(self) -> str:
        return f"{self.numerator}:{self.denominator}"

    def to_float(self) -> float:
        return self.numerator / self.denominator

    def apply(self, value: float) -> float:
        """How much of a do we get if a:b and we have an amount "value" of b?"""
        return value * self.numerator / self.denominator

    def inverse(self):
        return Ratio(self.denominator, self.numerator)

def define_ratio_from_string(ratio_string: str) -> Ratio:
    """Checks that the string ratio is correct (raises error if not)
    and returns the corresponding Ratio object."""
    try:
        left, right = ratio_string.split(':')
        if not left.isdigit() or not right.isdigit():
            raise ValueError("The numbers left and right of the colon are not integers.")
        return Ratio(int(left), int(right))
    except ValueError:
        raise ValueError("Input must be in the format 'int:int', e.g. '1:5'.")

class LocationTuple:
    def __init__(self, element1: float, element2: float):
        if (not isinstance(element1, (int, float))
            or not isinstance(element2, (int, float))
            or not 0 <= element1 <= 1
            or not 0 <= element2 <= 1
            or not element1 < element2
        ):
            raise TypeError(
                f"The location information from {element1} and {element2} are not valid.")
        self.minMaxLoc = (element1, element2)
    
    def return_within_tuple(self, number):
        if self.minMaxLoc[0] <= number <= self.minMaxLoc[1]:
            return True
        else:
            return False

def define_region_list(region_string:str, num_regions: int):
    list_of_regions = ast.literal_eval(region_string)
    if all(0<=region<num_regions for region in list_of_regions):
        return list_of_regions
    raise ValueError("The region string", region_string, "does not have a valid region.")

def define_tuple_of_locationPairTuples_from_string(location_string: str):
    list_of_location_tuples = ast.literal_eval(location_string)
    # Converts list of tuples to list of objects of class LocationTuple
    list_of_LocationTuples = []
    for tuple_idx, location_tuple in enumerate(list_of_location_tuples):
        list_of_LocationTuples.append(LocationTuple(*location_tuple))
    # Checks that the second element in each tuple object is smaller than the first
    # element in the next tuple (no overlap between tuples,
    # written from closest to 0 to closest to 1)
    for tuple_idx in range(len(list_of_LocationTuples[:-1])):
        if list_of_LocationTuples[tuple_idx].minMaxLoc[1] >= list_of_LocationTuples[tuple_idx+1].minMaxLoc[0]:
            raise ValueError(f"There are issues in the location_string {location_string}")
    
    return list_of_LocationTuples

def no_empty_cells(dataframe: pd.DataFrame):
    """Returns True if all the values in the cells exist. Else False. """
    return not dataframe.isna().any().any()

def no_repeated_rows_in_csv_file(csv_file: str):
    """
    Returns True if each row of the csv file is unique.
    Raises ValueError with duplicated rows otherwise.
    It's much easier to first check the csv files than to check the modified
    dataframe, pandas has issues if the types of the elements are not string.
    """
    df = pd.read_csv(csv_file)
    # Check for duplicated rows
    duplicates = df[df.duplicated()]

    if duplicates.empty:
        return True
    else:
        raise ValueError(f"Found {len(duplicates)} duplicated row(s):", duplicates)

def check_correct_type(dataframe: pd.DataFrame, columns_of_interest: list, expected_type) -> bool:
    """ Returns True or raises Error (good for asserting) if all the values within
    the columns listed in columns_of_interest have the type expected_type.
    """
    for col in columns_of_interest:
        if not dataframe[col].map(lambda x: isinstance(x, expected_type)).all():
            raise TypeError(f"Column '{col}' contains values not of type {expected_type}")
    return True

def checks_lack_of_repetitions(array)-> bool:
    """Returns True only if there are no repetitions in the array. 
    Raises an error if these are repeated (good for asserting).
    """
    unique_values, values_counts = np.unique(array, return_counts=True)
    if any(count for count in values_counts if count!=1):
        repetitions_dict = {value: count for value, count in zip(unique_values, values_counts)
                            if count!=1}
        raise ValueError("There are repeated values:", repetitions_dict)
    return True

def print_network_info(graph):
    """Prints all the graph info."""
    print("NODES:")
    for node in graph.nodes(data=True):
        print(f"  - {node}")
        print(f"    ↪ type: {type(node[0])}, attrs: {node[1]}")
    print("\nEDGES:")
    for u, v, attrs in graph.edges(data=True):
        print(f"  - from {u} to {v}")
        print(f"    ↪ types: {type(u)}, {type(v)}")
        print(f"    ↪ attrs: {attrs}")

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

def convert_sparse_txt_to_csv(filepath_txt: str, filepath_csv: str):
    """
    Read a sparse .txt file saved in (row, col, value) format
    and create a full CSV of the dense matrix.
    """
    # Load the sparse txt file, skipping header
    coo_data = np.loadtxt(filepath_txt, skiprows=1)
    if coo_data.ndim == 1:  # handle single nonzero element
        coo_data = coo_data.reshape(1,3)
    
    rows, cols, vals = coo_data.T
    rows = rows.astype(int)
    cols = cols.astype(int)
    
    # Determine matrix size
    nrows = rows.max() + 1
    ncols = cols.max() + 1
    
    # Build dense matrix
    dense_mat = np.zeros((nrows, ncols))
    dense_mat[rows, cols] = vals
    
    # Save as CSV
    np.savetxt(filepath_csv+".csv", dense_mat, delimiter=",", fmt="%.6e")
