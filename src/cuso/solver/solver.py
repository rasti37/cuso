"""Abstract solver class."""

from abc import abstractmethod
import logging
from typing import Any

class Solver:
    """Base class to represent a Solver.

    A Solver represents an instance of a problem with an easily checked
    solution. For debugging purposes, it's possible to set the expected
    solution to check intermediate results.
    """
    def __init__(self):
        self.logger = logging.getLogger("cuso." + self.__class__.__name__)

        self.has_solution = False
        self.expected = None

    @abstractmethod
    def check(self, solution: Any) -> bool:
        """Check if a solution satisfies the problem.

        Args:
            solution (Any): The solution being checked.

        Returns:
            bool: True if the solution is valid.
        """

    def set_expected(self, solution: Any):
        """Set the expected solution of the problem.

        Args:
            solution (Solution): Expected solution.
        """
        if not self.check(solution):
            self.logger.warning("Expected solution does not satisfy problem.")
        else:
            self.logger.debug("Expected solution satisfies problem.")
            self.has_solution = True
            self.expected = solution

    @abstractmethod
    def solve(self) -> Any:
        """Generate a solution."""