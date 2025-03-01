"""Base class for a solver of a multivariate Coppersmith problem"""

from abc import abstractmethod
from typing import Union

from cuso.solver.solver import Solver
from cuso.data.problem import MultivariateCoppersmithProblem
from cuso.data.solutions import PartialSolutionSet, SolutionSet


class MultivariateSolver(Solver):
    """Base class for a solver of a multivariate Coppersmith problem"""

    def __init__(self, problem: MultivariateCoppersmithProblem):
        """Construct a solver for a multivariate Coppersmith problem

        Args:
            problem (MultivariateCoppersmithProblem): the problem to be solved
        """
        super().__init__()
        self.problem = problem

    def check(self, solution: SolutionSet) -> bool:
        """Check whether a SolutionSet satisfies the Coppersmith problem

        Args:
            solution (SolutionSet): Set of solutions

        Returns:
            bool: True if the solutions satisfy the problem.
        """
        return self.problem.check(solution)

    @abstractmethod
    def solve(self) -> Union[PartialSolutionSet, SolutionSet]:
        """Find and return a set of solutions which solve the Coppersmith problem.

        Returns:
            SolutionSet: The set of recovered solutions.
        """
