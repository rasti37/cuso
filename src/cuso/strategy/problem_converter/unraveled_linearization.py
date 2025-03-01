"""Implement unraveled linearization via augmenting multivariate Coppersmith problems."""

import math
from operator import mul
from typing import Tuple, List

from sage.all import PolynomialRing, ZZ, TermOrder

from cuso.data import (
    BoundSet,
    MultivariateCoppersmithProblem,
    SolutionConverter,
    Solution,
    RelationConverter,
    Relation,
    RelationSet,
)
from cuso.data.types import Polynomial, Variable
from cuso.data.relations.converter import RenameRelationConverter

from .problem_converter import MultivariateProblemConverter


class UnravelSolution(SolutionConverter):
    """Convert a solution based on the given unraveled linearization terms."""

    def __init__(self, old_ring, new_ring, ul_terms):
        self.old_ring = old_ring
        self.new_ring = new_ring
        self.ul_terms = ul_terms

    def convert_solution_to_old(self, solution: Solution) -> Solution:
        old_xs = self.old_ring.gens()
        new_xs = self.new_ring.gens()
        old_soln = {}

        for new_xi in solution:
            if new_xi not in new_xs:
                # Modulus
                old_soln[new_xi] = solution[new_xi]
                continue

            ind = new_xs.index(new_xi)
            if ind < len(old_xs):
                old_soln[old_xs[ind]] = solution[new_xi]
            else:
                # This is u_i. Discard.
                pass
        return Solution(old_soln)

    def convert_solution_to_new(self, solution: Solution) -> Solution:
        new_soln = {}
        nvars_orig = len(self.old_ring.gens())
        nvars_new = len(self.new_ring.gens())
        old_xs = self.old_ring.gens()
        new_xs = self.new_ring.gens()

        # Add the x_s from the solution
        old_x_vals = []
        for i, old_xi in enumerate(self.old_ring.gens()):
            xi = new_xs[i]
            if old_xi in solution:
                old_x_vals += [solution[old_xi]]
            else:
                old_x_vals += [None]

        for xi in solution:
            if xi in old_xs:
                ind = old_xs.index(xi)
                new_soln[new_xs[ind]] = solution[xi]
            else:
                new_soln[xi] = solution[xi]

        # Add the u_s from the solution
        for i in range(nvars_new - nvars_orig):
            ul_term = self.ul_terms[i]
            new_val = int(ul_term(*old_x_vals))
            xi = new_xs[nvars_orig + i]
            new_soln[xi] = new_val

        return Solution(new_soln)


class UnraveledLinearization(MultivariateProblemConverter):
    """Implement the problem converter for unraveled linearization.

    This works by introducing new variables. If we have
        f(x, y) = x*y + Ax + 1 (mod N)
    and want to use unraveled linearization for x*y + 1, then this
    class introduces a new variable u, a new relation
        f_ul(x, y, u) = x*y - u + 1
    so that the ideal <f> + <f_ul> contains the unraveled
        g(x, y, u) = u + A*x
    """

    def __init__(self, unrav_lin_terms: List[Polynomial]):
        super().__init__()
        self.unrav_lin_terms = unrav_lin_terms

    def _get_new_ring(self, orig_problem):
        orig_relations = orig_problem.relations
        orig_bounds = orig_problem.bounds

        # Need to set the term order
        orig_ring = orig_relations.ring()
        xs = orig_ring.gens()

        # get new variable names
        ul_names = [f"u_{i+1}" for i in range(len(self.unrav_lin_terms))]
        for ul_name in ul_names:
            assert ul_name not in orig_ring.variable_names()
        ul_names = list(orig_ring.variable_names()) + ul_names

        orig_order = orig_ring.term_order()
        if orig_order.is_weighted_degree_order():
            weights = orig_order.weights()
        else:
            bound_vals = [orig_bounds.get_abs_bound(x) for x in xs]
            weights = tuple(math.log2(b) for b in bound_vals)

        ul_weights = list(weights)
        for ul_rel in self.unrav_lin_terms:
            degs = [monom.degrees() for monom in ul_rel.monomials()]
            weight = max(sum(map(mul, deg, weights)) for deg in degs)
            ul_weights += [weight]

        order = TermOrder("wdeglex", tuple(ul_weights))
        ring = PolynomialRing(ZZ, ul_names, order=order)
        return ring

    def _convert_bounds(self, bounds, orig_ring, new_ring):
        new_bounds = dict()
        for k, v in bounds.items():
            if k in orig_ring.gens():
                new_bounds[new_ring(k)] = v
            else:
                new_bounds[k] = v

        # Add UL bounds
        u_s = new_ring.gens()[-len(self.unrav_lin_terms) :]
        for ui, ul_term in zip(u_s, self.unrav_lin_terms):
            lbound = bounds.get_lower_bound(ul_term)
            ubound = bounds.get_upper_bound(ul_term)
            new_bounds[ui] = (lbound, ubound)

        new_bounds = BoundSet(new_bounds)
        return new_bounds

    def _convert_relations(
        self,
        old_relations: RelationSet,
        rel_converter: RelationConverter,
        u_s: List[Variable],
    ) -> RelationSet:
        # Convert old relations first
        rel_list = []
        for old_rel in old_relations:
            rel_list += [rel_converter.convert_relation_to_new(old_rel)]

        # Add new UL relations with integer constraints
        for u_i, ul_term in zip(u_s, self.unrav_lin_terms):
            new_poly = u_i - rel_converter.convert_polynomial_to_new(ul_term)
            rel_list += [Relation(new_poly)]

        return RelationSet(rel_list)

    def run(
        self, problem: MultivariateCoppersmithProblem
    ) -> Tuple[MultivariateCoppersmithProblem, RelationConverter, SolutionConverter]:
        orig_ring = problem.relations.ring()
        new_ring = self._get_new_ring(problem)
        new_vars = new_ring.gens()
        rel_converter = RenameRelationConverter(
            orig_ring.gens(), new_vars[: len(orig_ring.gens())]
        )
        u_s = new_vars[len(orig_ring.gens()) :]

        new_rels = self._convert_relations(problem.relations, rel_converter, u_s)

        new_bounds = self._convert_bounds(problem.bounds, orig_ring, new_ring)

        new_problem = MultivariateCoppersmithProblem(new_rels, new_bounds)
        soln_converter = UnravelSolution(orig_ring, new_ring, self.unrav_lin_terms)
        return new_problem, rel_converter, soln_converter
