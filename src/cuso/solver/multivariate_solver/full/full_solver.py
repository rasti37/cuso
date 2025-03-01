"""Base class for a FullSolver"""

from abc import abstractmethod

from cuso.data import SolutionSet
from ..multivariate_solver import MultivariateSolver


class FullSolver(MultivariateSolver):
    """Same as Solver, but may only return a full solution.

    This means that in a multivariate system (x, y), we must
    find values of both x and y.
    """

    @abstractmethod
    def solve(self) -> SolutionSet:
        """Find and return a set of full solutions which solve the Coppersmith problem.

        Returns:
            SolutionSet: The set of recovered solutions.
        """
