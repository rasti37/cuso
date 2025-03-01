"""This file implements the root recover strategy for primal lattices"""

from typing import Dict, List, Optional

from sage.all import Infinity, vector

from cuso.data import (
    Lattice,
    RelationSet,
    BoundSet,
    PartialSolutionSet,
    PartialSolution,
    SolutionSet,
)
from cuso.exceptions import SolveFailureError
from cuso.data.types import Monomial, Variable
from .root_recovery import RootRecovery


class PrimalRecovery(RootRecovery):
    """Root recovery for primal lattices"""

    def recover_short_vectors(self, lattice: Lattice) -> List[int]:
        """Return a list of indices where the lattice vector may correspond
        to a bounded root.

        Args:
            lattice (Lattice): reduced primal Coppersmith lattice

        Returns:
            List[int]: List of indices where the basis vector is sufficiently short.
        """
        # Look for vectors with infinity norm <= 1
        short_inds = []
        for i in range(lattice.rank()):
            try:
                vec = lattice.get_scaled_vector(i)
            except ValueError:
                # Could happen if the lattice is scaled by Infinity
                continue
            if vec.norm(Infinity) <= 1:
                short_inds += [i]
        return short_inds

    def get_solution_from_single_vector(
        self, vec: vector, monomials: List[Monomial], variables: List[Variable]
    ) -> PartialSolutionSet:
        """Given an unscaled vector and list of monomials, get the values
        of the variables.

        Args:
            vec (vector): unscaled vector where entry i is the assignment of
                the monomial monomials[i]
            monomials (List[Monomial]): list of monomials
            variables (List[Variable]): list of variables to recover

        Returns:
            PartialSolutionSet: Set of partial solutions.
        """
        if 1 not in monomials:
            raise SolveFailureError("Unable to confirm sign of monomials")

        one_ind = monomials.index(1)
        if vec[one_ind] < 0:
            # Flip sign
            vec *= -1

        root = {}
        for x_i in variables:
            if x_i in monomials:
                x_ind = monomials.index(x_i)
                root[x_i] = int(vec[x_ind])

        if len(root) == 0:
            return PartialSolutionSet()

        # Check for consistency
        for v_i, m_i in zip(vec, monomials):
            v_rec = m_i.subs(root)
            if v_i != v_rec:
                return PartialSolutionSet([])

        return PartialSolutionSet([PartialSolution(root)])

    def run(
        self,
        reduced_lattice: Lattice,
        input_relations: RelationSet,
        bounds: BoundSet,
        solver_kwargs: Dict = None,
        expected: Optional[SolutionSet] = None,
    ) -> PartialSolutionSet:
        # Get new integer relations
        short_vec_inds = self.recover_short_vectors(reduced_lattice)

        self.logger.debug(
            "Found %d short vectors in the primal lattice", len(short_vec_inds)
        )
        if len(short_vec_inds) == 0:
            raise SolveFailureError("Lattice reduction found no short vectors.")

        monoms = reduced_lattice.monomials
        if len(short_vec_inds) == 1 and 1 in monoms:
            xs = input_relations.variables()
            if set(xs).issubset(set(monoms)):
                short_vec = reduced_lattice.get_vector(short_vec_inds[0])
                return self.get_solution_from_single_vector(short_vec, monoms, xs)

        raise SolveFailureError("Unable to extract solution from primal lattice")
