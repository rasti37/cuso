"""This module implements shift polynomial selection strategies.

Given a multivariate Coppersmith problem, generate different
shift polynomial sets to try to solve with Coppersmith's method.
"""

from .shift_poly_selection import ShiftPolyStrategy

from .graph import GraphShiftPolys
from .optimal import OptimalShiftPolys
