"""Implementation of ideals of relations"""

from typing import List

from sage.all import gcd, Integer, PolynomialRing, Ideal

from cuso.data.types import Polynomial, Modulus


class RelationIdeal:
    """Implement ideals of relations.

    A RelationIdeal is a polynomial ideal in ZZ[x1,...] where all polynomials
    in the ideal share a constraint. That is, at a bounded solution $r$ for
    a Coppersmith problem, $f(r) == 0 (mod p)$ for all $f$ in the relation ideal.
    (The constraint could also be integer, so $f(r) = 0$). Relation ideals can
    be added and multiplied together to form new relation ideals.
    """

    def __init__(
        self, polynomials: List[Polynomial], ring: PolynomialRing, modulus: Modulus
    ):
        """Construct a RelationIdeal

        Args:
            polynomials (List[Polynomial]): Polynomial generators of the ideal
            ring (PolynomialRing): parent ring of the polynomials
            modulus (Modulus): Shared modulus constraint of all generators
        """
        self.ideal: Ideal = ring.ideal(polynomials)
        self._ring: PolynomialRing = ring
        self.modulus: Modulus = modulus

    def ring(self) -> PolynomialRing:
        """Get the parent ring

        Returns:
            PolynomialRing: ring containing this ideal
        """
        return self._ring

    def groebner_basis(self) -> List[Polynomial]:
        """Return the Groebner basis which generates this ideal.

        Returns:
            List[Polynomial]: Groebner basis
        """
        # Check if univariate
        if len(self._ring.gens()) == 1:
            # Add an unused variable so that the ring is multivariate
            new_ring = self._ring.extend_variables("tmp_unused")
            new_ideal = new_ring.ideal(self.ideal)
            new_groebner = new_ideal.groebner_basis()
            groebner = [self._ring(f) for f in new_groebner]
            if groebner == [0]:
                return []
            return groebner
        gb = self.ideal.groebner_basis()
        if gb == [0]:
            return []
        return gb

    def __add__(self, other: "RelationIdeal") -> "RelationIdeal":
        """Add two relation ideals together.

        If all polynomials in J1 have modulus p1, and all
        polynomials in J2 have modulus p2, then all polynomials
        in J1 + J2 have modulus lcm(p1, p2).

        Args:
            other (RelationIdeal): The other relation ideal to add

        Raises:
            TypeError: Incompatible relation ideals
            ValueError: Incompatible relation ideals

        Returns:
            RelationIdeal: The sum of ideals
        """
        if not isinstance(other, RelationIdeal):
            raise TypeError("Can only add relation ideal to relation ideal")
        if self.ring() != other.ring():
            raise ValueError("Relation ideals must have the same ring")
        new_ideal = self.ideal + other.ideal
        new_ring = self.ring()
        if self.modulus is None:
            new_modulus = other.modulus
        elif other.modulus is None:
            new_modulus = self.modulus
        else:
            new_modulus = gcd(self.modulus, other.modulus)
            if isinstance(new_modulus, Integer):
                new_modulus = int(new_modulus)
        return RelationIdeal(new_ideal.gens(), new_ring, new_modulus)

    def __mul__(self, other: "RelationIdeal") -> "RelationIdeal":
        """Multiply two relation ideals together.

        If all polynomials in J1 have modulus p1, and all
        polynomials in J2 have modulus p2, then all polynomials
        in J1 * J2 have modulus p1 * p2.

        Args:
            other (RelationIdeal): The other relation ideal to multiply

        Raises:
            TypeError: Incompatible relation ideals
            ValueError: Incompatible relation ideals

        Returns:
            RelationIdeal: The product of ideals
        """
        if not isinstance(other, RelationIdeal):
            raise TypeError("Can only multiply relation ideal by relation ideal")
        if self.ring() != other.ring():
            raise ValueError("Relation ideals must have the same ring")
        new_ideal = self.ideal * other.ideal
        new_ring = self.ring()
        if self.modulus is None:
            new_modulus = other.modulus
        elif other.modulus is None:
            new_modulus = self.modulus
        else:
            new_modulus = self.modulus * other.modulus
        return RelationIdeal(new_ideal.gens(), new_ring, new_modulus)
