"""This module implements Solutions to multivarite Coppersmith problems,
a SolutionSet, partial variants, and converters between different solutions
(for example, in the case where there is a change of variables).
"""

from .converter import SolutionConverter, RenameSolutionConverter
from .solutions import Solution, SolutionSet, PartialSolution, PartialSolutionSet
