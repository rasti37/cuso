"""Define a Shift Polynomial Strategy
"""

from abc import abstractmethod
from typing import Iterator

from cuso.data import RelationSet, BoundSet
from cuso.strategy.strategy import Strategy


class ShiftPolyStrategy(Strategy):
    """Abstract class that describes a shift polynomial selection strategy."""

    @abstractmethod
    def run(
        self, input_relations: RelationSet, bounds: BoundSet
    ) -> Iterator[RelationSet]:
        """Generate shift RelationSets from input RelationSets

        Args:
            input_relations (RelationSet): input relations
            bounds (BoundSet): root bounds

        Yields:
            Iterator[RelationSet]: Yields sets of relations that share the
                same bounded roots as the input relations.
        """
