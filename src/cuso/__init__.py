"""
This module provides an automated method for finding bounded roots of systems of polynomials.
"""

from . import exceptions, utils
from .data import (
    bounds,
    problem,
    relations,
    solutions,
    Relation,
    RelationSet,
    Solution,
    SolutionSet,
    BoundSet,
)
from .wrapper import find_small_roots
