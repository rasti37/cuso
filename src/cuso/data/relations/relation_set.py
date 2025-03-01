"""This module implements the RelationSet, which describes a set of relations."""

from typing import List, Union

from sage.all import PolynomialRing

from cuso.data.solutions import Solution
from cuso.data.types import Variable
from .relation import Relation


class RelationSet:
    """Represent a collection of Relations over a common set of variables."""

    _relations: List[Relation] = None
    _ring: PolynomialRing = None
    _unknown_moduli: List[Variable] = None

    def __init__(self, relations: Union[List[Relation], "RelationSet"]):
        if isinstance(relations, RelationSet):
            self._relations = relations._relations[:]
            self._ring = relations._ring
            self._unknown_moduli = relations._unknown_moduli
            return

        # Check types
        if not isinstance(relations, list):
            raise TypeError("RelationSet input must be list of Relations")
        if not all(isinstance(rel, Relation) for rel in relations):
            raise TypeError("RelationSet input must be list of Relations")

        # All relations must share the same ring
        if len(relations) > 0:
            ring = relations[0].ring()
            if not all(rel.ring() == ring for rel in relations):
                raise ValueError(
                    "All Relations in RelationSet must share a polynomial ring"
                )
        else:
            raise ValueError("List of relations cannot be empty.")

        unknown_moduli = []
        for rel in relations:
            unknown_moduli += rel.unknown_moduli()
        unknown_moduli = list(set(unknown_moduli))

        self._relations = relations[:]
        self._ring = ring
        self._unknown_moduli = unknown_moduli

    def ring(self) -> PolynomialRing:
        """Returns the polynomial ring.

        Returns:
            PolynomialRing: Common polynomial ring of relations
        """
        return self._ring

    def variables(self) -> List:
        """Return the list of variables in the polynomial ring.

        Returns:
            List[]: Generators of polynomial ring.
        """
        return self.ring().gens()

    def unknown_moduli(self):
        """Return the list of symbolic moduli.

        Returns:
            List[]: List of symbols that appear in moduli.
        """
        return self._unknown_moduli

    def check(self, solution: Solution) -> bool:
        """Check if the solution satisfies all relations.

        Args:
            solution (Solution): Dictionary of variable values.

        Returns:
            bool: The solution satisfies all relations.
        """
        for relation in self:
            if not relation.check(solution):
                return False
        return True

    def __getitem__(self, key):
        """Allow iteration over and accessing Relations in a RelationSet"""
        if isinstance(key, slice):
            return RelationSet(self._relations[key])
        return self._relations[key]

    def __len__(self) -> int:
        """Number of relations.

        Returns:
            int: Number of relations.
        """
        return len(self._relations)

    def __add__(self, other: "RelationSet") -> "RelationSet":
        """Allow adding two RelationSet objects together

        Args:
            other (RelationSet): RelationSet to add

        Returns:
            RelationSet: RelationSet containing all relations
        """
        other = RelationSet(other)
        return RelationSet(self._relations + other._relations)

    def __repr__(self):
        s = "RelationSet with "
        if len(self) == 1:
            s += f"1 relation in {self._ring.gens()[0]}"
        else:
            s += f"{len(self)} relations in {self._ring.gens()}"
        return s
