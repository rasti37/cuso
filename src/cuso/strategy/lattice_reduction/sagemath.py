"""Wrapper around Sage's lattice reduction methods"""

from sage.all import Matrix

from .lattice_reduction import LatticeReduction


class SageLatticeReduction(LatticeReduction):
    """Strategy that uses Sage's default lattice reduction method."""

    def reduce_integer_basis(self, basis: Matrix) -> Matrix:
        """Perform lattice basis reduction on an integer matrix.

        Args:
            basis (Matrix): input basis

        Returns:
            Matrix: reduced basis
        """
        self.logger.debug("Beginning lattice reduction.")
        red_basis = basis.LLL()
        self.logger.debug("Done.")
        return red_basis
