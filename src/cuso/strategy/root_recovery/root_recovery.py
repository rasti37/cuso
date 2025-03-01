"""Define a Root Recovery Strategy
"""

from abc import abstractmethod
from typing import Dict, Optional

from cuso.data import RelationSet, BoundSet, Lattice, PartialSolutionSet, SolutionSet
from cuso.strategy.strategy import Strategy


class RootRecovery(Strategy):
    """Abstract class that describes a way to recover roots from a reduced lattice."""

    @abstractmethod
    def run(
        self,
        reduced_lattice: Lattice,
        input_relations: RelationSet,
        bounds: BoundSet,
        solver_kwargs: Dict = None,
        expected: Optional[SolutionSet] = None,
    ) -> PartialSolutionSet:
        """Find a set of partial solutions based on a reduced lattice basis.

        Args:
            reduced_lattice (Lattice): Reduced lattice.
            input_relations (RelationSet): original relations with shared root.
            bounds (BoundSet): original bounds
            solver_kwargs (Dict, optional): Arguments used to initialize the
                parent multivariate Coppersmith solver.
            expected (Optional[SolutionSet], optional): Intended solution.

        Returns:
            PartialSolutionSet: Set of partial solutions
        """
