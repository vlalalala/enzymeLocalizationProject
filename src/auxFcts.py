from math import gcd

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
    try:
        left, right = ratio_string.split(':')
        if not left.isdigit() or not right.isdigit():
            raise ValueError("The numbers left and right of the colon are not integers.")
        return Ratio(int(left), int(right))
    except ValueError:
        raise ValueError("Input must be in the format 'int:int', e.g. '1:5'.")