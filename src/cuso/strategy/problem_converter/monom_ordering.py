"""Implement problem converter to set weight term order."""

import math
from typing import Tuple

from sage.all import PolynomialRing, ZZ, TermOrder

from cuso.data import (
    BoundSet,
    MultivariateCoppersmithProblem,
    SolutionConverter,
    RelationConverter,
)

from cuso.data.relations.converter import RenameRelationConverter
from cuso.data.solutions.converter import RenameSolutionConverter
from .problem_converter import MultivariateProblemConverter


class BoundedMonomialOrderConverter(MultivariateProblemConverter):
    """MultivariateProblemConverter that changes the TermOrder to a weight order
    based on the bounds.

    The weight of term x_i ^ e_i is e_i * log X_i.
    """

    def _get_new_ring(
        self, orig_ring: PolynomialRing, orig_bounds: BoundSet
    ) -> PolynomialRing:
        """Get the polynomial ring with the new term order.

        Args:
            orig_ring (PolynomialRing): original polynomial ring
            orig_bounds (BoundSet): original bounds

        Returns:
            PolynomialRing: Weight order polynomial ring.
        """
        # Need to set the term order
        xs = orig_ring.gens()
        if len(xs) == 1:
            # Univariate
            return orig_ring
        bound_vals = [orig_bounds.get_abs_bound(x) for x in xs]
        lg_bounds = [math.log2(b) if b > 1 else 1 for b in bound_vals]
        order = TermOrder("wdeglex", tuple(lg_bounds))
        ring = PolynomialRing(ZZ, orig_ring.variable_names(), order=order)
        return ring

    def _convert_bounds(
        self,
        bounds: BoundSet,
        orig_ring: PolynomialRing,
        new_ring: PolynomialRing,
    ) -> BoundSet:
        """Get new bounds.

        These are the same as the old bounds, but the variables are in a
        different polynomial ring.

        Args:
            bounds (BoundSet): original bounds
            orig_ring (PolynomialRing): original ring
            new_ring (PolynomialRing): new ring

        Returns:
            BoundSet: new bounds
        """
        new_bounds = {}
        for xi in bounds:
            if xi in orig_ring.gens():
                ind = orig_ring.gens().index(xi)
                new_xi = new_ring.gens()[ind]
                new_bounds[new_xi] = bounds[xi]
            else:
                new_bounds[xi] = bounds[xi]
        new_bounds = BoundSet(new_bounds)
        return new_bounds

    def run(
        self, problem: MultivariateCoppersmithProblem
    ) -> Tuple[MultivariateCoppersmithProblem, RelationConverter, SolutionConverter]:
        orig_ring = problem.relations.ring()
        new_ring = self._get_new_ring(orig_ring, problem.bounds)
        rel_converter = RenameRelationConverter(orig_ring.gens(), new_ring.gens())

        new_rels = rel_converter.convert_to_new(problem.relations)
        new_bounds = self._convert_bounds(problem.bounds, orig_ring, new_ring)

        new_problem = MultivariateCoppersmithProblem(new_rels, new_bounds)
        soln_converter = RenameSolutionConverter(orig_ring.gens(), new_ring.gens())
        return new_problem, rel_converter, soln_converter
