"""This module implements finding partial roots using Groebner bases."""

import warnings

from sage.all import QQ, Ideal

from cuso.data import RelationSet, PartialSolution, PartialSolutionSet
from cuso.exceptions import SolveFailureError
from .partial_solver import PartialSolver


class GroebnerSolver(PartialSolver):
    """Recover partial roots using Groebner bases.

    If there are relations with integer constraints in the input, then
    they share small roots over the rationals. We can use Groebner bases
    to compute the variety corresponding to the ideal, recovering the
    root as rationals. The integer roots are a subset, which we return.
    """

    def _get_variety(self, ideal: Ideal):
        try:
            return ideal.variety(algorithm="msolve", proof=False)
        except (ValueError, TypeError) as e:
            if "positive-dimensional ideal" in str(e):
                # Nothing we can do here
                raise
            # Maybe this is an old version of Sage. This might be recoverable.
            warnings.warn(
                "Variety calculation with msolve failed, it's possible that you need to upgrade "
                "SageMath to a newer version (9.6+) that allows calculating varieties with "
                "other algorithms, or you may need to install msolve (https://msolve.lip6.fr/). "
                "Attempting to use the old, slower implementation.",
                UserWarning
            )
        return ideal.variety()

    def _get_partial_roots(self, relations: RelationSet):
        int_rels = [rel for rel in relations if rel.modulus is None]
        if len(int_rels) == 0:
            raise SolveFailureError("No integer relations for Groebner solver.")
        self.logger.debug("Using Groebner to try to solve integer relations")

        polys = [rel.polynomial for rel in int_rels]
        ring = polys[0].parent()
        if len(ring.gens()) == 1:
            # Univariate case
            # Get roots of first polynomial
            roots = [x_val for x_val, _ in polys[0].roots()]

            # Check to see which roots satisfy all polynomials
            roots = [root for root in roots if all(f(root) == 0 for f in polys)]
            x = ring.gens()[0]
            solns = [{x: root} for root in roots]
            return solns
        else:
            # See if there are any variables that don't appear in
            # the polynomials. If there are, add a dummy equation
            # to the ideal so we get a 0-dimensional variety.
            unused_gens = []
            for xi in ring.gens():
                deg = max(f.degree(xi) for f in polys)
                if deg == 0:
                    unused_gens += [xi]
            for xi in unused_gens:
                polys += [xi]

            ring_Q = ring.change_ring(QQ, order="degrevlex")
            ideal = ring_Q.ideal(polys)
            try:
                variety = self._get_variety(ideal)
            except ValueError as exc:
                # raise SolveFailureError("Positive dimension ideal.") from exc
                # Positive dimension ideal
                # Can we get any partial solutions from the GB?
                gb = ideal.groebner_basis()
                soln = {}
                improved = True
                while improved:
                    improved = False
                    for g in gb:
                        g = g.subs(soln)
                        if g == 0:
                            continue
                        if sum(g.degrees()) == 1:
                            # Then this equation is a * x + b == 0, so we can solve for x
                            (x,) = g.variables()
                            if x in unused_gens:
                                continue
                            a = g.monomial_coefficient(x)
                            b = g.constant_coefficient()
                            val_x = int(-b // a)
                            soln[x] = val_x
                            improved = True
                if len(soln) == 0:
                    # Unable to find any solutions
                    raise SolveFailureError("Positive dimension ideal.") from exc
                soln = {ring(x_i): v_i for x_i, v_i in soln.items()}
                solns = [soln]
                return solns
            if len(variety) == 0:
                raise SolveFailureError("Variety is empty.")
            solns = []
            for point in variety:
                soln = {}
                for xi in point.keys():
                    if ring(xi) not in unused_gens:
                        soln[ring(xi)] = int(point[xi])
                solns += [soln]
            return solns

    def solve(self) -> PartialSolutionSet:
        prob = self.problem
        rels = prob.relations
        roots = self._get_partial_roots(rels)
        self.logger.info("Found %u possible (partial) root(s)", len(roots))
        return PartialSolutionSet(
            [PartialSolution(root) for root in roots if prob.bounds.check(root)]
        )
