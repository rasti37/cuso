"""This module implements automatically finding partial roots that solve Coppersmith problems."""

from cuso.exceptions import SolveFailureError
from cuso.data.problem import MultivariateCoppersmithProblem
from cuso.data.solutions import PartialSolutionSet
from .partial_solver import PartialSolver
from .coppersmith import CoppersmithSolver
from .groebner import GroebnerSolver
from .linear import LinearSolver


class AutomatedPartialSolver(PartialSolver):
    """Automatic recovery of partial roots in multivariate Coppersmith problems

    We try multiple different solvers, and return the partial solution if any
    of them succeed. First, we try the linear solver, which uses lattice reduction
    when there are enough linear relations to find a small number of solutions.
    Second, we try the Groebner solver, which uses Groebner bases to find the
    variety corresponding to relations with integer constraints. Third, we try
    the Coppersmith solver, which generates shift ideals and Coppersmith lattices.
    """

    linear_solver = None
    newton_solver = None
    groebner_solver = None
    coppersmith_solver = None

    def __init__(
        self,
        problem: MultivariateCoppersmithProblem,
        unraveled_linearization_relations=None,
        use_intermediate_sizes=None,
        use_graph_optimization=None,
        use_primal_strategy=None,
    ):
        super().__init__(problem)

        self.linear_solver = LinearSolver(
            problem,
        )
        self.groebner_solver = GroebnerSolver(
            problem,
        )
        self.coppersmith_solver = CoppersmithSolver(
            problem,
            unraveled_linearization_relations=unraveled_linearization_relations,
            use_intermediate_sizes=use_intermediate_sizes,
            use_graph_optimization=use_graph_optimization,
            use_primal_strategy=use_primal_strategy,
        )

    def solve(self) -> PartialSolutionSet:
        solvers = [self.linear_solver, self.groebner_solver, self.coppersmith_solver]
        # Try to solve using different methods
        for solver in solvers:
            try:
                if self.has_solution:
                    solver.set_expected(self.expected)
                return solver.solve()
            except SolveFailureError:
                continue
        raise SolveFailureError("All multivariate solvers failed.")
