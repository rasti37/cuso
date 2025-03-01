"""This file contains the strategy to build a Coppersmith primal lattice."""

from sage.all import Matrix, ZZ, QQ, Infinity

from cuso.data import RelationSet, BoundSet, Lattice
from .lattice_builder import LatticeBuilder


class PrimalLatticeBuilder(LatticeBuilder):
    """Strategy to build a Coppersmith primal lattice."""

    def run(self, relations: RelationSet, bounds: BoundSet) -> Lattice:
        if not isinstance(relations, RelationSet):
            raise TypeError("LatticeBuilder requires RelationSet as input")
        if not isinstance(bounds, BoundSet):
            raise TypeError("LatticeBuilder requires BoundSet as input")

        monomials = []
        for rel in relations:
            monomials += rel.polynomial.monomials()
        monomials = sorted(list(set(monomials)))

        # Get number of relations with a modular constraint
        num_mod_rels = 0
        for rel in relations:
            if rel.modulus is None:
                continue
            if isinstance(rel.modulus, int):
                num_mod_rels += 1
            else:
                raise TypeError("Cannot build Primal lattice from symbolic modulus")

        num_monoms = len(monomials)
        rank = num_monoms + num_mod_rels
        dimension = num_monoms + len(relations)
        self.logger.debug(
            "Building lattice of rank %d and dimension %d", rank, dimension
        )

        M = Matrix(ZZ, rank, dimension)

        for i in range(rank):
            M[i, i] = 1

        mod_rel_ind = 0
        for j, rel in enumerate(relations):
            # Column num_monoms + j encodes this relation
            f = rel.polynomial
            for m in f.monomials():
                i = monomials.index(m)
                cij = f.monomial_coefficient(m)
                M[i, num_monoms + j] = cij
            if rel.modulus is not None:
                # Add the modulus to the lower right corner
                M[num_monoms + mod_rel_ind, num_monoms + j] = rel.modulus
                mod_rel_ind += 1

        # scale
        scale_factors = []
        for m in monomials:
            Xj = bounds.get_abs_bound(m)
            scale_factors += [1 / QQ(Xj)]
        scale_factors += [Infinity] * len(relations)

        L = Lattice(
            M,
            monomials,
            scale_factors=scale_factors,
            is_primal=True,
        )
        self.logger.info(
            "Built a primal lattice basis of rank %u and dimension %u",
            L.rank(),
            L.dimension(),
        )
        return L
