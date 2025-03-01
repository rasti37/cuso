"""This module implements the AutomatedSolver"""

from functools import reduce

from sage.all import Expression, gcd, PolynomialRing, ZZ

from cuso.data import (
    MultivariateCoppersmithProblem,
    Relation,
    RelationSet,
    BoundSet,
    Solution,
    SolutionSet,
    PartialSolution,
    PartialSolutionSet,
)
from ..partial import AutomatedPartialSolver
from .full_solver import FullSolver


class AutomatedSolver(FullSolver):
    """Implementation of the Automated Coppersmith solver.

    At a high level, use the AutomatedPartial solver to find partial solutions
    to the Coppersmith problem. If some variable values were not recovered,
    use the recovered values to set up a new Coppersmith problem to recover
    the remaining values.
    """

    def __init__(self, problem: MultivariateCoppersmithProblem, *args, **kwargs):
        super().__init__(problem)
        self.my_args = args
        self.my_kwargs = kwargs

    def _get_unknown_moduli(
        self, relations: RelationSet, soln: PartialSolution
    ) -> PartialSolution:
        soln = soln.copy()

        evaluations = []
        moduli = []
        root = {x_i: v for x_i, v in soln.items() if x_i in relations.ring().gens()}

        for rel in relations:
            evaluation = rel.polynomial.subs(root)
            try:
                evaluation = int(evaluation)
                evaluations += [evaluation]
                moduli += [rel.modulus]
            except TypeError:
                # Could fail if we only have solutions for some of the variables
                continue

        for p_i in relations.unknown_moduli():
            # If p_i is in the modulus, compute GCD
            eval_p = []
            mod_p = []
            for eval_i, mod_i in zip(evaluations, moduli):
                if isinstance(mod_i, Expression) and mod_i.degree(p_i) > 0:
                    eval_p += [eval_i]
                    mod_p += [mod_i]
            eval_p = reduce(gcd, eval_p, 0)
            mod_p = reduce(gcd, mod_p, 0)
            if mod_p != p_i:
                raise NotImplementedError

            # Could still have problems if working mod p * q and
            # finding multiples of p and of q
            soln[p_i] = int(eval_p)

        return soln

    def _get_remaining_unknowns(
        self, relations: RelationSet, solution: PartialSolution, bounds: BoundSet
    ) -> SolutionSet:
        old_ring = relations.ring()
        old_vars = [xi for xi in old_ring.gens() if xi not in solution]
        if len(old_vars) == 0:
            solns = SolutionSet([Solution(solution)])
            if self.problem.check(solns):
                return solns
            else:
                return SolutionSet()

        new_varnames = ",".join([str(xi) for xi in old_vars])
        new_ring = PolynomialRing(ZZ, new_varnames)
        varsubs = {old_xi: new_xi for old_xi, new_xi in zip(old_vars, new_ring.gens())}
        varsubs.update(solution)

        new_rels = []
        for rel in relations:
            new_poly = rel.polynomial.subs(varsubs)
            new_poly = new_ring(str(new_poly))
            if new_poly == 0:
                continue
            if isinstance(rel.modulus, int):
                new_mod = rel.modulus
            else:
                if rel.modulus is not None:
                    new_mod = rel.modulus.subs(solution)
                else:
                    new_mod = None
            new_rels += [Relation(new_poly, new_mod)]
        new_relations = RelationSet(new_rels)

        new_bounds = {}
        for k, v in bounds.items():
            if k in solution:
                continue
            if k in varsubs:
                new_bounds[varsubs[k]] = v
            else:
                new_bounds[k] = v
        new_bounds = BoundSet(new_bounds)

        prob = MultivariateCoppersmithProblem(new_relations, new_bounds)
        solver = AutomatedSolver(prob, *self.my_args, **self.my_kwargs)
        if self.has_solution:
            new_expected = []
            for soln in self.expected:
                # Convert from old ring to new ring
                new_soln = {}
                for xi in soln:
                    if xi in old_ring.gens():
                        if xi not in solution:
                            new_soln[varsubs[xi]] = soln[xi]
                    else:
                        new_soln[xi] = soln[xi]
                new_expected += [new_soln]
            new_expected = SolutionSet(new_expected)
            solver.set_expected(new_expected)

        solns = []
        new_solns = solver.solve()
        for new_soln in new_solns:
            # Convert from new ring to old ring
            orig_soln = {}
            for xi in new_soln:
                if xi in new_ring.gens():
                    orig_soln[old_ring(xi)] = new_soln[xi]
                else:
                    orig_soln[xi] = new_soln[xi]
            orig_soln.update(solution)
            solns += [orig_soln]
        return [Solution(soln) for soln in solns]

    def _expand_solution(
        self,
        partial_solns: PartialSolutionSet,
        relations: RelationSet,
        bounds: BoundSet,
    ) -> SolutionSet:
        solns = [
            self._get_unknown_moduli(relations, partial_soln)
            for partial_soln in partial_solns
        ]
        ext_solns = []
        for soln in solns:
            self.logger.info("Found partial solution")
            for xi, vi in soln.items():
                self.logger.info("\t%s = %s", str(xi), str(vi))
            ext_solns += self._get_remaining_unknowns(relations, soln, bounds)
        return SolutionSet(ext_solns)

    def solve(self) -> SolutionSet:
        prob = self.problem

        rels = prob.relations

        auto_partial_solver = AutomatedPartialSolver(
            self.problem, *self.my_args, **self.my_kwargs
        )
        if self.has_solution:
            auto_partial_solver.set_expected(self.expected)
        partial_solns = auto_partial_solver.solve()

        solns = self._expand_solution(partial_solns, rels, prob.bounds)
        return solns
