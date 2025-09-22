import numpy as np
from math import gcd
import pandas as pd

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

def no_empty_cells(dataframe: pd.DataFrame):
    """Returns True if all the values in the cells exist. Else False. """
    return not dataframe.isna().any().any()

def check_correct_type(dataframe, columns_of_interest, expected_type):
    for col in columns_of_interest:
        print(col, expected_type)
        if not dataframe[col].map(lambda x: isinstance(x, expected_type)).all():
            raise TypeError(f"Column '{col}' contains values not of type {expected_type}")
    return True

#%%
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

