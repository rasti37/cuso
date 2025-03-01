"""This file implements the root recover strategy for dual lattices"""

from typing import Dict, Optional

from sage.all import sqrt

from cuso.data import (
    Lattice,
    Relation,
    RelationSet,
    BoundSet,
    PartialSolutionSet,
    SolutionSet,
    MultivariateCoppersmithProblem,
)
from cuso.solver import multivariate_solver
from cuso.exceptions import SolveFailureError
from .root_recovery import RootRecovery


class HastadHowgraveGraham(RootRecovery):
    """Root recovery for dual lattices"""

    def __init__(self, use_l1_of_vectors: Optional[bool] = None):
        """Construct the strategy.

        Args:
            use_l1_of_vectors (bool, optional): Use the L1 norm to
                determine if vectors are short enough rather than
                the L2 norm.
        """
        super().__init__()
        if use_l1_of_vectors is None:
            use_l1_of_vectors = True
        self.use_l1_of_vectors: bool = use_l1_of_vectors

    def recover_int_relations(self, lattice: Lattice, bounds: BoundSet) -> RelationSet:
        """Given a Coppersmith dual lattice, recover the basis relations such that
        the evaluation at the bounded root is exactly 0.

        Args:
            lattice (Lattice): Reduced dual lattice
            bounds (BoundSet): Bounds of Coppersmith problem

        Raises:
            ValueError: Invalid parameters
            SolveFailureError: No short vectors in lattice.

        Returns:
            RelationSet: Integer relations which share a bounded root.
        """
        modulus = lattice.modulus
        if modulus is None:
            raise ValueError("HG root recovery requires modular relation ideal.")

        # Get Howgrave-Graham bound
        max_l1_norm = bounds.get_lower_bound(modulus)

        short_vector_inds = []
        if self.use_l1_of_vectors:
            # Get the vectors which satisfy the L1 norm bound
            for i in range(lattice.rank()):
                vec = lattice.get_scaled_vector(i)
                if vec.norm(1) < max_l1_norm:
                    short_vector_inds += [i]
        else:
            # Use the L2 norm to bound the L1 norm
            dim = lattice.dimension()
            for i in range(lattice.rank()):
                vec = lattice.get_scaled_vector(i)
                if vec.norm(2) / sqrt(dim) < max_l1_norm:
                    short_vector_inds += [i]

        self.logger.info(
            "Found %d integer relations in the lattice", len(short_vector_inds)
        )
        if len(short_vector_inds) == 0:
            raise SolveFailureError("Lattice reduction found no short vectors.")

        # Convert these vector to integer relations
        self.logger.debug("Found new integer relations:")
        rels = []
        for ind in short_vector_inds:
            rel = lattice.get_relation(ind)
            # Set modulus to None
            rel = Relation(rel.polynomial, None)
            rels += [rel]
            self.logger.debug("%s", str(rel.polynomial))

        return RelationSet(rels)

    def run(
        self,
        reduced_lattice: Lattice,
        input_relations: RelationSet,
        bounds: BoundSet,
        solver_kwargs: Dict = None,
        expected: Optional[SolutionSet] = None,
    ) -> PartialSolutionSet:
        # Get new integer relations
        int_rels = self.recover_int_relations(reduced_lattice, bounds)
        new_relations = RelationSet(input_relations + int_rels)

        new_prob = MultivariateCoppersmithProblem(new_relations, bounds)

        solver_kwargs = solver_kwargs if solver_kwargs is not None else {}
        solver = multivariate_solver.AutomatedPartialSolver(new_prob, **solver_kwargs)
        if expected is not None:
            solver.set_expected(expected)
        return solver.solve()
