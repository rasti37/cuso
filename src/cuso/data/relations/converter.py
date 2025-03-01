"""This file contains classes which convert Relation and RelationSet types.

It is used, for example, when we have a change of varibles and need to convert
the relation which uses the old variables to one which uses the new variables.
"""

from abc import abstractmethod
from typing import List, Any

from cuso.data.types import Polynomial
from .relation import Relation
from .relation_set import RelationSet


class RelationConverter:
    """Generic class for converting Relation and RelationSet types."""

    @abstractmethod
    def convert_polynomial_to_new(self, poly: Polynomial) -> Polynomial:
        """Convert a Polynomial of the old type to the new type.

        Args:
            poly (Polynomial): Polynomial in the old problem.

        Returns:
            Polynomial: Polynomial in the new problem.
        """

    def convert_relation_to_new(self, relation: Relation) -> Relation:
        """Convert a Relation of the old type to the new type.

        Args:
            relation (Relation): Relation in the old problem.

        Returns:
            Relation: Relation in the new problem.
        """
        return Relation(
            self.convert_polynomial_to_new(relation.polynomial), relation.modulus
        )

    def convert_to_new(self, relations: RelationSet) -> RelationSet:
        """Convert a RelationSet of the old type to the new type.

        Args:
            relations (RelationSet): Relations in the old problem.

        Returns:
            RelationSet: Relations in the new problem.
        """
        return RelationSet([self.convert_relation_to_new(rel) for rel in relations])

    @abstractmethod
    def convert_polynomial_to_old(self, poly: Polynomial) -> Polynomial:
        """Convert a Polynomial of the new type to the old type.

        Args:
            poly (Polynomial): Polynomial in the new problem.

        Returns:
            Polynomial: Polynomial in the old problem.
        """

    def convert_relation_to_old(self, relation: Relation) -> Relation:
        """Convert a Relation of the new type to the old type.

        Args:
            relation (Relation): Relation in the new problem.

        Returns:
            Relation: Relation in the old problem.
        """
        return Relation(
            self.convert_polynomial_to_old(relation.polynomial), relation.modulus
        )

    def convert_to_old(self, relations: RelationSet) -> RelationSet:
        """Convert a RelationSet of the new type to the old type.

        Args:
            relations (RelationSet): Relations in the new problem.

        Returns:
            RelationSet: Relations in the old problem.
        """
        return RelationSet([self.convert_relation_to_old(rel) for rel in relations])


class RenameRelationConverter(RelationConverter):
    """Implementation that just renames variables in the relations."""

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
        if len(orig) == 0:
            raise ValueError("Must have at least one named element")
        self.orig = orig
        self.new = new

        self.orig_ring = orig[0].parent()
        self.new_ring = new[0].parent()

    def convert_polynomial_to_new(self, poly: Polynomial) -> Polynomial:
        """Convert a Polynomial of the old type to the new type.

        Args:
            poly (Polynomial): Polynomial in the old problem.

        Returns:
            Polynomial: Polynomial in the new problem.
        """
        to_subs = {o: n for o, n in zip(self.orig, self.new)}
        new_poly = self.new_ring(poly.subs(to_subs))
        return new_poly

    def convert_polynomial_to_old(self, poly: Polynomial) -> Polynomial:
        """Convert a Polynomial of the new type to the old type.

        Args:
            poly (Polynomial): Polynomial in the new problem.

        Returns:
            Polynomial: Polynomial in the old problem.
        """
        to_subs = {n: o for o, n in zip(self.orig, self.new)}
        old_poly = self.orig_ring(poly.subs(to_subs))
        return old_poly
