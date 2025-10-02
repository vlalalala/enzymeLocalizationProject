import numpy as np
from math import gcd
import pandas as pd
import pickle
import yaml
from pathlib import Path
import ast

def load_yaml_as_dict(yaml_file_path: str):
    return yaml.safe_load(Path(yaml_file_path).read_text())

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

def define_list_of_locationPairTuples_from_string(location_string: str):
    list_of_location_tuples = ast.literal_eval(location_string)
    # Converts list of tuples to list of objects of class LocationTuple
    for tuple_idx, location_tuple in enumerate(list_of_location_tuples):
        list_of_location_tuples[tuple_idx] = LocationTuple(*location_tuple)
    # Checks that the second element in each tuple object is smaller than the first
    # element in the next tuple (no overlap between tuples,
    # written from closest to 0 to closest to 1)
    for tuple_idx in range(list_of_location_tuples[:-1]):
        if list_of_location_tuples[idx].minMaxLoc[1] >= list_of_location_tuples[idx+1].minMaxLoc[0]:
            raise ValueError(f"There are issues in the location_string {location_string}")
    return list_of_location_tuples

def no_empty_cells(dataframe: pd.DataFrame):
    """Returns True if all the values in the cells exist. Else False. """
    return not dataframe.isna().any().any()

def no_repeated_rows(dataframe: pd.DataFrame):
    """Returns True if each row of the dataframe is only once. Else False."""
    duplicate_rows = dataframe.duplicated().to_dict()
    duplicate_rows = {row: is_dup for row, is_dup in duplicate_rows.items() if is_dup}
    if len(duplicate_rows) == 0:
        return True
    else:
        raise ValueError("The rows", duplicate_rows.keys(), "are found multiple times.")

def check_correct_type(dataframe: pd.DataFrame, columns_of_interest: list, expected_type) -> bool:
    for col in columns_of_interest:
        if not dataframe[col].map(lambda x: isinstance(x, expected_type)).all():
            raise TypeError(f"Column '{col}' contains values not of type {expected_type}")
    return True

def checks_lack_of_repetitions(array)-> bool:
    """Returns True only if there are no repetitions in the array.
    Raises an error if these are repeated"""
    unique_values, values_counts = np.unique(array, return_counts=True)
    if any(count for count in values_counts if count!=1):
        repetitions_dict = {value: count for value, count in zip(unique_values, values_counts)
                            if count!=1}
        raise ValueError("There are repeated values:", repetitions_dict)
    else:
        return True

def pickle_dump_binary(path, variable):
    with open(path, 'wb') as f:
        pickle.dump(variable, f)

def pickle_load_binary(path):
    with open(path, 'rb') as f:
        loaded_variable = pickle.load(f)
    return loaded_variable

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

