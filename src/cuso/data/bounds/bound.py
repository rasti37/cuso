"""Bound class"""

import math
from typing import Optional

from cuso.data.types import Variable


class Bound:
    """Implementation of a bound.

    It may be bounded above or bounded below or both.
    """

    def __init__(self, x: Variable, lower: Optional[int], upper: Optional[int]):
        """Initialize the Bound

        Args:
            x (Variable): Variable being bounded.
            lower (Optional[int]): lower bound, if one exists, else None
            upper (Optional[int]): upper bound, if one exists, else None
        """
        self.x = x
        self.lower = lower
        self.upper = upper

    def __contains__(self, value):
        if self.lower is not None and value <= self.lower:
            return False
        if self.upper is not None and value >= self.upper:
            return False
        return True

    def __iter__(self):
        return iter((self.lower, self.upper))

    def _bound_as_str(self, boundval) -> str:
        if abs(boundval) < 1000:
            return str(boundval)
        sgn = boundval // abs(boundval)
        lg = math.log2(abs(boundval))
        return f"{sgn * 2}^{lg:.1f}"

    def __repr__(self):
        s = ""
        if self.lower is not None:
            s += self._bound_as_str(self.lower) + " < "
        s += str(self.x)
        if self.upper is not None:
            s += " < " + self._bound_as_str(self.upper)
        return s
