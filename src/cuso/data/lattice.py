"""Implement a generic Lattice as is used in Coppersmith's method.

The lattice may use the primal or the dual Coppersmith construction.
"""

from typing import List, Optional

from sage.all import vector, Infinity, Integer, Matrix

from .types import Modulus, Monomial
from .relations import Relation


class Lattice:
    """Class representing a Coppersmith lattice."""

    def __init__(
        self,
        basis: Matrix,
        monomials: List[Monomial],
        is_primal: bool = None,
        modulus: Optional[Modulus] = None,
        scale_factors: Optional[List[Integer]] = None,
    ):
        """Coppersmith lattice.

        Args:
            basis (Matrix): Unscaled basis matrix.
            monomials (List[Monomial]): Ordered list of monomials.
            is_primal (bool, optional): True if the lattice vectors represent solutions,
                false if the lattice vectors represent polynomials.
            modulus (Modulus, optional): Shared modulus of all relations, if there is one.
            scale_factors (List[Integer], optional): Amount to scale each column by.
                An entry in this list can be Infinity, so lattice reduction infinitely
                penalizes any nonzero value in the corresponding column.
        """
        self.basis: Matrix = basis
        self.monomials: List[Monomial] = monomials
        self.modulus: Modulus = modulus
        self.scale_factors: List[Integer] = scale_factors
        self.is_primal: bool = is_primal

    def get_scaled_vector(self, index: int) -> vector:
        """Return the scaled basis vector in row INDEX.

        Args:
            index (int): The row to return.

        Raises:
            ValueError: The scaled vector is infinitely large.

        Returns:
            vector: scaled basis vector
        """
        vec = self.get_vector(index)
        if self.scale_factors is None:
            scaled_vec = vec
        else:
            scaled_vec = []
            for s_i, v_i in zip(self.scale_factors, vec):
                if s_i == Infinity:
                    if v_i == 0:
                        scaled_vec += [0]
                    else:
                        raise ValueError("Vector is infinitely large")
                else:
                    scaled_vec += [s_i * v_i]
            scaled_vec = vector(scaled_vec)
        return scaled_vec

    def get_vector(self, index: int) -> vector:
        """Return the unscaled basis vector in row INDEX

        Args:
            index (int): The row to return

        Raises:
            TypeError: Index must be integer
            IndexError: Index out of range

        Returns:
            vector: unscaled basis vector
        """
        if not isinstance(index, int):
            raise TypeError("Index must be integer")
        if index < 0 or index >= self.basis.nrows():
            raise IndexError

        return self.basis[index]

    def get_relation(self, index: int) -> Relation:
        """Return the relation corresponding to basis vector INDEX.

        In a dual lattice, each vector corresponds to a relation.
        
        Args:
            index (int): Index of the relation to return

        Raises:
            TypeError: Lattice is not dual.

        Returns:
            Relation: basis relation at the specified index.
        """
        if self.is_primal is not False:
            raise TypeError("Can only get relation from dual lattice.")
        vec = self.get_vector(index)
        poly = 0
        for v_i, m_i in zip(vec, self.monomials):
            poly += v_i * m_i
        return Relation(poly, self.modulus)

    def rank(self) -> int:
        """Get the rank of the lattice.

        Returns:
            int: rank
        """
        return self.basis.nrows()

    def dimension(self) -> int:
        """Get the dimension of the lattice.

        Returns:
            int: dimension
        """
        return self.basis.ncols()
