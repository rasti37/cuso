"""This module implements finding partial roots using shift polynomials."""

from typing import Dict

from cuso.data.problem import MultivariateCoppersmithProblem
from cuso.exceptions import SolveFailureError
from cuso.strategy.problem_converter import (
    MultivariateProblemConverter,
    ChainConverter,
)
from cuso.strategy.shift_polynomial_selection import (
    ShiftPolyStrategy,
    GraphShiftPolys,
    OptimalShiftPolys,
)
from cuso.strategy.lattice_reduction import LatticeReduction, Flatter
from cuso.strategy.lattice_builder import (
    LatticeBuilder,
    DualLatticeBuilder,
    PrimalLatticeBuilder,
)
from cuso.strategy.root_recovery import (
    RootRecovery,
    HastadHowgraveGraham,
    PrimalRecovery,
)

from .partial_solver import PartialSolver


class CoppersmithSolver(PartialSolver):
    """Recover partial roots using shift polynomials.

    We use the input relations to generate shift relations, then
    for each set of shift relations we build and reduce the
    Coppersmith lattice. We try to recover roots. If this succeeds,
    return the partial roots. If it fails, try the next set of
    shift relations.
    """

    problem_converter: MultivariateProblemConverter = None
    shift_poly_strategy: ShiftPolyStrategy = None
    lattice_building_strategy: LatticeBuilder = None
    lattice_reduction_strategy: LatticeReduction = None
    root_recovery_strategy: RootRecovery = None
    use_primal_strategy: bool = None
    use_graph_optimization: bool = None

    def __init__(
        self,
        problem: MultivariateCoppersmithProblem,
        unraveled_linearization_relations=None,
        use_intermediate_sizes=None,
        use_graph_optimization=None,
        use_primal_strategy=None,
    ):
        super().__init__(problem=problem)
        self._ul = unraveled_linearization_relations
        self._solver_kwargs = {
            "use_intermediate_sizes": use_intermediate_sizes,
            "use_graph_optimization": use_graph_optimization,
            "use_primal_strategy": use_primal_strategy,
        }

        mod_mul_known, has_symbolic_mod = self._get_configuration()

        if use_primal_strategy is None:
            # Determine whether or not to use the primal strategy
            if mod_mul_known or has_symbolic_mod:
                use_primal_strategy = False
            else:
                use_primal_strategy = True

        if use_graph_optimization is None:
            use_graph_optimization = self._get_use_graph_optimization()

        if use_intermediate_sizes is None:
            use_intermediate_sizes = True

        # Check to see if the configuration is valid
        self.problem_converter = ChainConverter(
            unrav_lin_relations=unraveled_linearization_relations
        )

        shift_poly_strategy = OptimalShiftPolys(
            use_intermediate_sizes=use_intermediate_sizes
        )
        if use_graph_optimization:
            if not mod_mul_known:
                raise ValueError(
                    "Graph optimization requires a known multiple of the modulus"
                )
            shift_poly_strategy = GraphShiftPolys(shift_poly_strategy)
        self.shift_poly_strategy = shift_poly_strategy

        self.lattice_reduction_strategy = Flatter()

        if use_primal_strategy:
            self.lattice_building_strategy = PrimalLatticeBuilder()
            self.root_recovery_strategy = PrimalRecovery()
        else:
            self.lattice_building_strategy = DualLatticeBuilder()
            self.root_recovery_strategy = HastadHowgraveGraham(
                use_l1_of_vectors=True,
            )

    def _get_use_graph_optimization(self) -> bool:
        mod_mul_known, _ = self._get_configuration()
        if not mod_mul_known:
            # Graph optimization would fail here.
            return False

        # If we're using unraveled linearization, that's probably
        # a hint to use graph optimization.
        if self._ul is not None and len(self._ul) > 0:
            return True

        # Look for small coefficients.
        for rel in self.problem.relations:
            coefs = rel.polynomial.coefficients()
            num_small_coefs = len([c for c in coefs if abs(c) < 100])
            if num_small_coefs > 1:
                self.logger.info(
                    "Small coefficients detected, disabling graph optimization."
                )
                self.logger.info(
                    "Graph optimization can be enabled by passing use_graph_optimization=True"
                    " to find_small_roots()"
                )
                return False

        # It's probably OK to use the graph optimization
        return True

    def _get_configuration(self):
        # Check if a multiple of the modulus is known
        modulus_multiple_known = False
        for rel in self.problem.relations:
            if rel.modulus is not None:
                if rel.polynomial.is_constant():
                    # Constant polynomial with an unknown modulus
                    modulus_multiple_known = True
                    break
                if isinstance(rel.modulus, int):
                    # Known modulus
                    modulus_multiple_known = True
                    break

        # Check if any of the relations have a symbolic modulus
        has_symbolic_moduli = False
        for rel in self.problem.relations:
            if rel.modulus is not None and not isinstance(rel.modulus, int):
                has_symbolic_moduli = True

        return modulus_multiple_known, has_symbolic_moduli

    def solve(self) -> Dict:
        shift_poly_strat = self.shift_poly_strategy
        lattice_building_strategy = self.lattice_building_strategy
        latred_strategy = self.lattice_reduction_strategy
        root_recovery = self.root_recovery_strategy

        relations = self.problem.relations
        num_modular = len([r for r in relations if r.modulus is not None])
        num_integer = len(relations) - num_modular
        num_variables = len(relations.ring().gens())
        self.logger.info(
            "Attempting the automated multivariate Coppersmith method with %u variable(s)",
            num_variables,
        )
        self.logger.info(
            "Multivariate Coppersmith problem has %d modular relation(s) "\
            "and %d integer relation(s)",
            num_modular,
            num_integer,
        )
        new_prob, _, soln_converter = self.problem_converter.run(self.problem)

        if self.has_solution:
            _expected = soln_converter.convert_to_new(self.expected)
        else:
            _expected = None

        input_rels = new_prob.relations
        input_bounds = new_prob.bounds

        for shift_rels in shift_poly_strat.run(input_rels, input_bounds):
            self.logger.debug("Generated %d shift relations", len(shift_rels))
            if self.has_solution:
                for soln in _expected:
                    if shift_rels.check(soln):
                        self.logger.debug("Expected value satisfies shift relations.")
                    else:
                        self.logger.warning(
                            "Expected value does not satisfy shift relations."
                        )

            try:
                lattice = lattice_building_strategy.run(shift_rels, input_bounds)
                reduced_lattice = latred_strategy.run(lattice)
                solutions = root_recovery.run(
                    reduced_lattice,
                    input_rels,
                    input_bounds,
                    solver_kwargs=self._solver_kwargs,
                    expected=_expected,
                )
                solutions = soln_converter.convert_to_old(solutions)
                return solutions
            except SolveFailureError:
                self.logger.debug("Coppersmith method failed with given shift polys")
