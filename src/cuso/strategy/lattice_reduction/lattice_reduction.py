"""Define a Lattice Reduction Strategy
"""

from abc import abstractmethod
from typing import List

from sage.all import copy, lcm, QQ, Infinity, Integer, Matrix

from cuso.data import Lattice
from cuso.strategy.strategy import Strategy


class LatticeReduction(Strategy):
    """Abstract class that describes a shift polynomial selection strategy.

    This class allows us to work with "infinite" scale factors, so if a
    column is scaled by +Infinity, then the shortest vectors in the reduced
    lattice basis will have 0 in the corresponding column. This is achieved
    by replacing the infinite scale factor with increasingly large scale
    factors.
    """

    @abstractmethod
    def reduce_integer_basis(self, basis: Matrix) -> Matrix:
        """Perform lattice basis reduction on an integer matrix.

        Args:
            basis (Matrix): input basis

        Returns:
            Matrix: reduced basis
        """
        raise NotImplementedError

    def scale_basis(self, basis: Matrix, scale_factors: List[Integer]) -> Matrix:
        """Take an unscaled basis and scale it by finite scale factors

        Args:
            basis (Matrix): unscaled basis
            scale_factors (List[Integer]): scale factors

        Returns:
            Matrix: scaled basis
        """
        basis = copy(basis)
        for i, s_i in enumerate(scale_factors):
            basis[:, i] *= s_i
        return basis

    def unscale_basis(self, basis: Matrix, scale_factors: List[Integer]) -> Matrix:
        """Take a scaled basis and unscale it by finite scale factors

        Args:
            basis (Matrix): scaled basis
            scale_factors (List[Integer]): scale factors

        Returns:
            Matrix: unscaled basis
        """
        basis = copy(basis)
        for i, s_i in enumerate(scale_factors):
            basis[:, i] /= s_i
        return basis

    def run(self, lattice: Lattice) -> Lattice:
        """Perform lattice reduction on a Coppersmith lattice.

        Args:
            lattice (Lattice): input lattice

        Returns:
            Lattice: reduced lattice
        """
        basis = lattice.basis
        scale_factors = lattice.scale_factors
        rank = lattice.rank()

        self.logger.debug("Beginning lattice reduction.")

        # Get common denominator
        denom = lcm([QQ(x).denominator() for x in scale_factors if x < Infinity])
        inf_mask = [x == Infinity for x in scale_factors]
        inf_factor = 1 << 1000  # Large value to scale by

        # Get integer scale factors
        int_scale_factors = [
            int(denom * s_i) if s_i < Infinity else inf_factor for s_i in scale_factors
        ]

        scaled_basis = self.scale_basis(basis, int_scale_factors)
        while True:
            red_scaled_basis = self.reduce_integer_basis(scaled_basis)

            # We're done if the inf-scaled entries are all 0 for the short vectors
            expected_zero_rows = rank - sum(inf_mask)
            insufficient_reduction = False
            for j, mask_j in enumerate(inf_mask):
                if mask_j:
                    if red_scaled_basis[:expected_zero_rows, j] != 0:
                        insufficient_reduction = True
                        break
            if insufficient_reduction:
                # Need to scale up the inf-scaled columns to penalize nonzeros even more
                inf_factor *= inf_factor
                scaled_basis = red_scaled_basis
                for j, mask_j in enumerate(inf_mask):
                    if mask_j:
                        int_scale_factors[j] *= inf_factor
                        scaled_basis[:, j] *= inf_factor
            else:
                break
        red_basis = self.unscale_basis(red_scaled_basis, int_scale_factors)

        self.logger.debug("Done.")
        red_lattice = Lattice(
            red_basis,
            lattice.monomials,
            modulus=lattice.modulus,
            scale_factors=scale_factors,
            is_primal=lattice.is_primal,
        )
        return red_lattice
