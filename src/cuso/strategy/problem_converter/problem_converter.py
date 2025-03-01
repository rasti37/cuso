"""Convert one multivariate Coppersmith problem to another
"""

from abc import abstractmethod
from typing import Tuple

from cuso.data import (
    MultivariateCoppersmithProblem,
    SolutionConverter,
    RelationConverter,
)
from cuso.strategy.strategy import Strategy


class MultivariateProblemConverter(Strategy):
    """Abstract class that describes a shift polynomial selection strategy."""

    @abstractmethod
    def run(
        self, problem: MultivariateCoppersmithProblem
    ) -> Tuple[MultivariateCoppersmithProblem, RelationConverter, SolutionConverter]:
        """Convert a multivariate Coppersmith problem to another.

        Args:
            problem (MultivariateCoppersmithProblem): Problem to convert.

        Returns:
            Tuple[MultivariateCoppersmithProblem, RelationConverter, SolutionConverter]:
                Tuple of new problem, converter that maps relations from the original
                problem to the new problem, and converter that map solutions
                from the original problem to the new problem.
        """
