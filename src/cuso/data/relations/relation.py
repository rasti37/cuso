"""Represent a constraint, or an integer polynomial with either a modular or integer constraint."""

from typing import Optional, Union, List

from sage.all import Polynomial, Expression, ZZ, Integer, PolynomialRing
from sage.rings.polynomial.multi_polynomial import MPolynomial

from cuso.data.solutions import Solution
from cuso.data.types import Variable, Modulus


class Relation:
    """Coppersmith relation of a polynomial and a constraint."""

    def __init__(self, polynomial: Polynomial, modulus: Optional[Modulus] = None):
        """Initializes the Relation.

        Args:
            polynomial (Polynomial): Polynomial defined in ZZ[x]
            modulus (Optional[Modulus], optional): Modulus the Relation satisfies.
                If None, the relation is 0 at the small root. Defaults to None.
        """
        # Check type
        if not isinstance(polynomial, (Polynomial, MPolynomial)):
            raise TypeError("Relation polynomial must be SageMath Polynomial.")
        if not polynomial.base_ring() == ZZ:
            raise TypeError("Relation polynomial must be defined in ZZ[]")
        if modulus is not None and not isinstance(modulus, (int, Integer, Expression)):
            raise TypeError("Modulus must be integer, expression, or None")

        if isinstance(modulus, Integer):
            modulus = int(modulus)

        self.polynomial = polynomial
        self.modulus = modulus

    def ring(self) -> PolynomialRing:
        """Get the Polynomial Ring for this relation.

        Returns:
            PolynomialRing: The ring
        """
        return self.polynomial.parent()

    def variables(self) -> List[Variable]:
        """Get the variables in the polynomial ring

        Returns:
            List[Variable]: List of polynomial ring generators
        """
        return self.ring().gens()

    def unknown_moduli(self) -> List[Variable]:
        """Get the list of unknown moduli"""
        if self.modulus is None or isinstance(self.modulus, int):
            return []
        return self.modulus.variables()

    def check(self, solution: Solution) -> bool:
        """Check if the solution satisfies the constrained relation.

        Args:
            solution (Solution): Values of unknowns.

        Returns:
            bool: True if the values satisfy the relation.
        """
        root = tuple(solution[xi] for xi in self.variables())
        value = self.polynomial(*root)
        if self.modulus is None:
            return value == 0
        if isinstance(self.modulus, int):
            mod_value = self.modulus
        else:
            mod_value = self.modulus.subs(dict(solution))
            mod_value = int(mod_value)
        return value % mod_value == 0

    def __mul__(self, other: Union["Relation", Polynomial]) -> "Relation":
        """Multiply two Relations together.

        The modulus of the product is the product of the individual moduli.

        Args:
            other (Union[Relation;, Polynomial]): The polynomial to multiply this relation by.

        Raises:
            ValueError: The polynomial rings of each term must match

        Returns:
            Relation: the product of both relations.
        """
        if isinstance(other, Relation):
            if self.ring() != other.ring():
                raise ValueError("Rings do not match.")
            new_polynomial = self.polynomial * other.polynomial
            if self.modulus is None:
                new_modulus = other.modulus
            elif other.modulus is None:
                new_modulus = self.modulus
            else:
                new_modulus = self.modulus * other.modulus
        else:
            new_polynomial = self.polynomial * self.ring()(other)
            new_modulus = self.modulus
        return Relation(new_polynomial, new_modulus)

    def __repr__(self):
        s = "Relation("
        s += str(self.polynomial)
        if self.modulus is None:
            s += " == 0"
        else:
            s += " == 0 modulo " + str(self.modulus)
        s += ")"
        return s
