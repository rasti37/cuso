"""This file contains tools for centering the unknowns around zero."""

from typing import Tuple, List

from sage.all import PolynomialRing, ZZ

from cuso.data import (
    BoundSet,
    MultivariateCoppersmithProblem,
    SolutionConverter,
    Solution,
    RelationConverter,
)

from .problem_converter import MultivariateProblemConverter


class RecenterSolution(SolutionConverter):
    """Implementation of SolutionConverter which transforms variables
    from x to x_c so that the bounds of x_c are [-B, B].
    """

    def __init__(
        self, old_ring: PolynomialRing, new_ring: PolynomialRing, centers: List[int]
    ):
        """Initialize RecenterSolution converter

        Args:
            old_ring (PolynomialRing): original ring
            new_ring (PolynomialRing): new ring with centered variables
            centers (List[int]): how much to shift each variable by
        """
        self.old_ring: PolynomialRing = old_ring
        self.new_ring: PolynomialRing = new_ring
        self.centers: List[int] = centers

    def convert_solution_to_old(self, solution: Solution) -> Solution:
        old_soln = {}
        for xi in solution:
            if xi in self.new_ring.gens():
                ind = self.new_ring.gens().index(xi)
                old_xi = self.old_ring.gens()[ind]
                old_soln[old_xi] = solution[xi] + self.centers[ind]
            else:
                old_soln[xi] = solution[xi]
        return Solution(old_soln)

    def convert_solution_to_new(self, solution: Solution) -> Solution:
        new_soln = {}
        for xi in solution:
            if xi in self.old_ring.gens():
                ind = self.old_ring.gens().index(xi)
                new_xi = self.new_ring.gens()[ind]
                new_soln[new_xi] = solution[xi] - self.centers[ind]
            else:
                new_soln[xi] = solution[xi]
        return Solution(new_soln)


class RecenterRelation(RelationConverter):
    """Implementation of RelationConverter which transforms variables
    from x to x_c so that the bounds of x_c are [-B, B]"""

    def __init__(
        self, orig_ring: PolynomialRing, new_ring: PolynomialRing, centers: List[int]
    ):
        """Initialize RecenterRelation converter

        Args:
            old_ring (PolynomialRing): original ring
            new_ring (PolynomialRing): new ring with centered variables
            centers (List[int]): how much to shift each variable by
        """
        self.orig_ring = orig_ring
        self.new_ring = new_ring
        self.centers = centers

    def convert_polynomial_to_new(self, poly):
        orig_ring = self.orig_ring
        new_ring = self.new_ring
        centers = self.centers

        varsub = {
            old_xi: new_xi + center_i
            for old_xi, new_xi, center_i in zip(
                orig_ring.gens(), new_ring.gens(), centers
            )
        }
        new_poly = new_ring(poly.subs(varsub))
        return new_poly

    def convert_polynomial_to_old(self, poly):
        orig_ring = self.orig_ring
        new_ring = self.new_ring
        centers = self.centers

        varsub = {
            new_xi: old_xi - center_i
            for old_xi, new_xi, center_i in zip(
                orig_ring.gens(), new_ring.gens(), centers
            )
        }
        orig_poly = orig_ring(poly.subs(varsub))
        return orig_poly


class RecenterConverter(MultivariateProblemConverter):
    """Implementation of MultivariateProblem Converter which transforms variables
    from x to x_c so that the bounds of x_c are [-B, B]"""

    def _convert_bounds(
        self,
        bounds: BoundSet,
        centers: List[int],
        orig_ring: PolynomialRing,
        new_ring: PolynomialRing,
    ) -> BoundSet:
        new_bounds = {}
        for xi in bounds:
            if xi in orig_ring.gens():
                lbound, ubound = bounds[xi]
                ind = orig_ring.gens().index(xi)
                new_xi = new_ring.gens()[ind]
                new_lbound = lbound - centers[ind]
                new_ubound = ubound - centers[ind]
                new_bounds[new_xi] = (new_lbound, new_ubound)
            else:
                new_bounds[xi] = bounds[xi]
        new_bounds = BoundSet(new_bounds)
        return new_bounds

    def run(
        self, problem: MultivariateCoppersmithProblem
    ) -> Tuple[MultivariateCoppersmithProblem, RelationConverter, SolutionConverter]:
        ring = problem.relations.ring()
        bounds = problem.bounds

        centers = []
        new_varnames = []
        for xi in ring.gens():
            # Get old bounds
            lbound = bounds.get_lower_bound(xi)
            ubound = bounds.get_upper_bound(xi)
            center = (lbound + ubound) // 2
            centers += [center]
            if center == 0:
                new_varnames += [str(xi)]
            else:
                new_name = f"{str(xi)}_c"
                new_varnames += [new_name]

        # Compute new problem
        orig_ring = problem.relations.ring()
        new_ring = PolynomialRing(ZZ, ",".join(new_varnames))
        rel_converter = RecenterRelation(orig_ring, new_ring, centers)
        new_rels = rel_converter.convert_to_new(problem.relations)
        new_bounds = self._convert_bounds(problem.bounds, centers, orig_ring, new_ring)

        for old_xi, new_xi, center in zip(orig_ring.gens(), new_ring.gens(), centers):
            if center == 0:
                continue
            self.logger.info(
                "Recentering variable %s from %s to %s",
                str(old_xi),
                bounds[old_xi],
                new_bounds[new_xi],
            )
            self.logger.debug("%s = %s - %d", str(new_xi), str(old_xi), center)

        new_problem = MultivariateCoppersmithProblem(new_rels, new_bounds)
        soln_converter = RecenterSolution(orig_ring, new_ring, centers)
        return new_problem, rel_converter, soln_converter
