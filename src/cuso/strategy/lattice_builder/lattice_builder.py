"""Strategy to build Coppersmith lattices."""

from abc import abstractmethod

from cuso.data import RelationSet, BoundSet, Lattice
from ..strategy import Strategy


class LatticeBuilder(Strategy):
    """Lattice builder strategy."""

    @abstractmethod
    def run(self, relations: RelationSet, bounds: BoundSet) -> Lattice:
        """Return the Coppersmith lattice for a set of relations and bounds.

        Args:
            relations (RelationSet): shift polynomial relations
            bounds (BoundSet): bounds on the desired root

        Returns:
            Lattice: Coppersmith lattice
        """
