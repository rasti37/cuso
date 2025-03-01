"""Base class for a PartialSolver"""

from abc import abstractmethod

from cuso.data import PartialSolutionSet
from ..multivariate_solver import MultivariateSolver


class PartialSolver(MultivariateSolver):
    """Same as Solver, but may only return a partial solution.

    This means that in a multivariate system (x, y), we may
    have found the value of x in a small root, but not the value of y.
    """

    @abstractmethod
    def solve(self) -> PartialSolutionSet:
        """Find and return a set of partial solutions which solve the Coppersmith problem.

        Returns:
            PartialSolutionSet: The set of recovered solutions.
        """
