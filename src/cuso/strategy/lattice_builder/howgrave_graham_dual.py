"""This file contains the strategy to build a Howgrave-Graham dual lattice."""

from sage.all import Matrix, ZZ

from cuso.data import RelationSet, BoundSet, Lattice
from .lattice_builder import LatticeBuilder


class DualLatticeBuilder(LatticeBuilder):
    """Strategy to build a Howgrave-Graham dual lattice."""

    def run(self, relations: RelationSet, bounds: BoundSet) -> Lattice:
        if not isinstance(relations, RelationSet):
            raise TypeError("LatticeBuilder requires RelationSet as input")
        if not isinstance(bounds, BoundSet):
            raise TypeError("LatticeBuilder requires BoundSet as input")

        # Ensure all relations share a modulus
        if len(relations) == 0:
            raise ValueError("Must specify at least one relation")
        modulus = relations[0].modulus
        if not all(rel.modulus == modulus for rel in relations):
            raise ValueError("All shift relations must share the same modulus")

        monomials = []
        for rel in relations:
            monomials += rel.polynomial.monomials()
        monomials = sorted(list(set(monomials)))

        rank = len(relations)
        dimension = len(monomials)
        self.logger.debug(
            "Building lattice of rank %d and dimension %d", rank, dimension
        )

        M = Matrix(ZZ, rank, dimension)
        for i, rel in enumerate(relations):
            f = rel.polynomial
            for m in f.monomials():
                j = monomials.index(m)
                cij = f.monomial_coefficient(m)
                M[i, j] = cij

        # scale
        scale_factors = []
        for m in monomials:
            Xj = bounds.get_abs_bound(m)
            scale_factors += [Xj]

        L = Lattice(
            M, monomials, modulus=modulus, scale_factors=scale_factors, is_primal=False
        )
        self.logger.info(
            "Built a dual lattice basis of rank %u and dimension %u",
            L.rank(),
            L.dimension(),
        )
        return L
