"""Coppersmith solver for when input relations are linear."""

from typing import Tuple

from sage.all import Infinity, lcm, sqrt, ceil
from fpylll import IntegerMatrix
from fpylll.fplll.gso import MatGSO
from fpylll.fplll.enumeration import Enumeration, EnumerationError

from cuso.exceptions import SolveFailureError
from cuso.strategy.lattice_builder import PrimalLatticeBuilder
from cuso.strategy.lattice_reduction import Flatter
from cuso.strategy.problem_converter import RecenterConverter
from cuso.data import (
    BoundSet,
    Lattice,
    MultivariateCoppersmithProblem,
    PartialSolution,
    PartialSolutionSet,
    Relation,
    RelationSet,
    SolutionConverter,
)
from .partial_solver import PartialSolver


class LinearSolver(PartialSolver):
    """Solver for the subset of Coppersmith problem where all polynomials are linear.

    We can just use normal lattice reduction (and lattice enumeration of small vectors)
    to recover the set of small solutions.
    """

    def _get_linear_relations(self, rels: RelationSet):
        """Get the subset of input relations that are linear"""

        linear_rels = []
        for rel in rels:
            if rel.modulus is not None and not isinstance(rel.modulus, int):
                # No symbolic moduli here
                continue
            if len(rel.polynomial.parent().gens()) == 1:
                tot_degree = rel.polynomial.degree()
            else:
                tot_degree = rel.polynomial.degree(std_grading=True)
            if tot_degree <= 1:
                linear_rels += [rel]
        if len(linear_rels) == 0:
            raise SolveFailureError("No linear relations")
        rels = []
        for rel in linear_rels:
            if rel.modulus is None:
                rels += [rel]
            else:
                rels += [Relation(rel.polynomial % rel.modulus, rel.modulus)]
        return RelationSet(rels)

    def _build_lattice(self, rels: RelationSet, bounds: BoundSet) -> Lattice:
        builder = PrimalLatticeBuilder()
        lattice = builder.run(rels, bounds)
        return lattice

    def _reduce_lattice(self, lattice: Lattice) -> Lattice:
        reducer = Flatter()
        red_lat = reducer.run(lattice)
        return red_lat

    def _recover_short_vectors(
        self, lattice: Lattice, bounds: BoundSet
    ) -> PartialSolutionSet:
        """Look for all basis vectors with infinity norm <= 1"""
        short_inds = []
        denom = 1
        for i in range(lattice.rank()):
            try:
                vec = lattice.get_scaled_vector(i)
            except ValueError:
                # Could happen if the lattice is scaled by Infinity
                continue
            if vec.norm(Infinity) <= 1:
                short_inds += [i]
                denom = lcm(denom, vec.denominator())
        if len(short_inds) == 0:
            raise SolveFailureError("No sufficiently short vectors")
        if len(short_inds) == 1:
            vec = lattice.get_vector(short_inds[0])
            partial_soln = {}
            if 1 in lattice.monomials:
                one_ind = lattice.monomials.index(1)
                if vec[one_ind] == -1:
                    vec = tuple(-v_i for v_i in vec)
                elif vec[one_ind] != 1:
                    return PartialSolutionSet()
            for i, m_i in enumerate(lattice.monomials):
                if m_i != 1:
                    partial_soln[m_i] = vec[i]
            if not bounds.check(partial_soln):
                return PartialSolutionSet()
            return PartialSolutionSet([PartialSolution(partial_soln)])
        if len(short_inds) >= 2:
            # Use length of shortest vector to see if running the lattice enumerator
            # would return too many values to be useful
            shvec = lattice.get_scaled_vector(short_inds[0])
            shvec_len = shvec.norm(Infinity)
            if shvec_len < 1 / 100:
                raise SolveFailureError("Too many candidate solutions")

        A = IntegerMatrix(len(short_inds), len(lattice.monomials))
        A_unscaled = IntegerMatrix(len(short_inds), len(lattice.monomials))
        for i, ii in enumerate(short_inds):
            vec = lattice.get_scaled_vector(ii)
            vec_unscaled = lattice.get_vector(ii)
            for j in range(len(lattice.monomials)):
                A[i, j] = int(vec[j] * denom)
                A_unscaled[i, j] = int(vec_unscaled[j])
        dim = len(lattice.monomials)
        bound_len = int(ceil(denom * sqrt(dim)))
        M = MatGSO(A)
        M_mpfr = MatGSO(A, float_type="mpfr")
        _ = M.update_gso()
        M_mpfr.update_gso()

        partial_solns = []

        def callbackf(v):
            nonlocal partial_solns
            v = tuple(round(v_i) for v_i in v)
            vec = A_unscaled.multiply_left(v)

            partial_soln = {}
            if 1 in lattice.monomials:
                one_ind = lattice.monomials.index(1)
                if vec[one_ind] == -1:
                    vec = tuple(-v_i for v_i in vec)
                elif vec[one_ind] != 1:
                    return False
            for i, m_i in enumerate(lattice.monomials):
                if m_i != 1:
                    partial_soln[m_i] = vec[i]
            if not bounds.check(partial_soln):
                return False

            partial_solns += [partial_soln]
            if len(partial_solns) > 100:
                # Too many possible solutions!
                # Return True here to halt enumeration
                return True
            return False

        enum = Enumeration(M, callbackf=callbackf)
        target = [0] * len(lattice.monomials)
        if 1 in lattice.monomials:
            target[lattice.monomials.index(1)] = int(denom)
        target = M_mpfr.from_canonical(target)
        exp = bound_len.bit_length()
        bound_len /= 1 << exp
        self.logger.debug("Trying lattice enumeration of linear relations.")
        try:
            enum.enumerate(0, A.nrows, bound_len**2, 2 * exp, target=target)
            # If we got here, then that means over 100 solutions were found.
            raise SolveFailureError("Too many candidates")
        except EnumerationError:
            # Expected.
            pass
        return PartialSolutionSet([PartialSolution(soln) for soln in partial_solns])

    def _get_centered_problem(self) -> Tuple[MultivariateCoppersmithProblem, SolutionConverter]:
        old_prob = self.problem
        new_prob, _, soln_conv = RecenterConverter().run(old_prob)
        return new_prob, soln_conv

    def solve(self) -> PartialSolutionSet:
        prob, soln_conv = self._get_centered_problem()
        rels = self._get_linear_relations(prob.relations)
        bounds = prob.bounds
        lattice = self._build_lattice(rels, bounds)
        red_lat = self._reduce_lattice(lattice)
        roots = self._recover_short_vectors(red_lat, bounds)
        self.logger.info("Found %u possible (partial) root(s)", len(roots))
        new_roots = PartialSolutionSet(
            [PartialSolution(root) for root in roots if bounds.check(root)]
        )
        return soln_conv.convert_to_old(new_roots)
