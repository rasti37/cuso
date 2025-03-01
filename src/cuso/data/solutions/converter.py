"""This file contains classes which convert Solution and SolutionSet types.

It is used, for example, when we have a change of varibles and need to convert
the solution which uses the new variables to one which uses the old variables.
"""

from abc import abstractmethod
from typing import Any, List

from .solutions import SolutionSet, Solution


class SolutionConverter:
    """Generic class for converting Solution and SolutionSet types."""

    @abstractmethod
    def convert_solution_to_new(self, solution: Solution) -> Solution:
        """Convert a Solution of the old type to the new type.

        Args:
            solution (Solution): Solution of the old problem.

        Returns:
            Solution: Solution of the new problem.
        """

    def convert_to_new(self, solutions: SolutionSet) -> SolutionSet:
        """Convert a SolutionSet of the old type to the new type.

        Args:
            solutions (SolutionSet): SolutionSet of the old problem.

        Returns:
            SolutionSet: SolutionSet of the new problem.
        """
        return SolutionSet([self.convert_solution_to_new(soln) for soln in solutions])

    @abstractmethod
    def convert_solution_to_old(self, solution: Solution) -> Solution:
        """Convert a Solution of the new type to the old type.

        Args:
            solution (Solution): Solution of the new problem.

        Returns:
            Solution: Solution of the old problem.
        """

    def convert_to_old(self, solutions: SolutionSet) -> SolutionSet:
        """Convert a SolutionSet of the new type to the old type.

        Args:
            solutions (SolutionSet): SolutionSet of the new problem.

        Returns:
            SolutionSet: SolutionSet of the old problem.
        """
        return SolutionSet([self.convert_solution_to_old(soln) for soln in solutions])


class RenameSolutionConverter(SolutionConverter):
    """Implementation that just renames variables in the solution."""

    def __init__(self, orig: List[Any], new: List[Any]):
        """Construct this type.

        Args:
            orig (List[Any]): old variables
            new (List[Any]): new variables

        Raises:
            ValueError: invalid parameters
        """
        if len(orig) != len(new):
            raise ValueError("Old and new must be same length")
        self.orig = orig
        self.new = new

    def convert_solution_to_new(self, solution: Solution) -> Solution:
        """Convert a Solution of the old type to the new type.

        Args:
            solution (Solution): Solution of the old problem.

        Returns:
            Solution: Solution of the new problem.
        """
        to_subs = {o: n for o, n in zip(self.orig, self.new)}
        new_soln = {to_subs.get(k, k): v for k, v in solution.items()}
        return Solution(new_soln)

    def convert_solution_to_old(self, solution: Solution) -> Solution:
        """Convert a Solution of the new type to the old type.

        Args:
            solution (Solution): Solution of the new problem.

        Returns:
            Solution: Solution of the old problem.
        """
        to_subs = {n: o for o, n in zip(self.orig, self.new)}
        old_soln = {to_subs.get(k, k): v for k, v in solution.items()}
        return Solution(old_soln)
