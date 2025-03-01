"""MultivariateCoppersmithProblem class"""

from typing import List
import logging

from sage.all import Integer

from .relations import RelationSet
from .bounds import BoundSet
from .solutions import SolutionSet
from .types import Variable

__all__ = ["MultivariateCoppersmithProblem"]

logger = logging.getLogger("cuso.MultivariateCoppersmithProblem")


class MultivariateCoppersmithProblem:
    """Representation of a multivariate Coppersmith problem"""

    def __init__(self, relations: RelationSet, bounds: BoundSet):
        """Initialize.

        Args:
            relations (RelationSet): Input relations.
            bounds (BoundSet): Bounds on the desired solutions.
        """
        self.relations: RelationSet = relations
        self.bounds: BoundSet = bounds

    def variables(self) -> List[Variable]:
        """Get variables appearing in the relations.

        Returns:
            List[Variables]: List of variables
        """
        return self.relations.variables()

    def unknown_moduli(self) -> List[Variable]:
        """Get variables appearing in the moduli.

        Returns:
            List[Variable]: List of symbols in the moduli.
        """
        return self.relations.unknown_moduli()

    def check(self, solutions: SolutionSet) -> bool:
        """Check whether the solutions satisfy the relations.

        Args:
            solutions (SolutionSet): List of solutions.

        Raises:
            TypeError: Solution is in invalid format.

        Returns:
            bool: True if the solutions satisfy the problem.
        """
        if not isinstance(solutions, SolutionSet):
            raise TypeError("Expected solutions to be SolutionSet")

        # Check that all unknown appear in the solution
        variables = self.variables()
        moduli_variables = self.unknown_moduli()
        for soln in solutions:
            for xi in variables:
                if xi not in soln:
                    logger.warning("Value of %s missing from solution %s", xi, soln)
                    return False
                if not isinstance(soln[xi], (int, Integer)):
                    raise TypeError(f"Solution value {soln[xi]} must be integer.")
            for pi in moduli_variables:
                if pi not in soln:
                    logger.warning("Value of modulus %s missing from %s", pi, soln)
                    return False
                if not isinstance(soln[pi], (int, Integer)):
                    raise TypeError(f"Solution value {soln[pi]} must be integer.")

        for soln in solutions:
            # Check that the solution satisfies the bounds
            if not self.bounds.check(soln):
                return False
            if not self.relations.check(soln):
                return False

        return True

    def __repr__(self):
        variables = self.variables()
        nvars = len(variables)
        nrels = len(self.relations)
        s = "Multivariate Coppersmith Problem in "
        if nvars == 1:
            s += f"{nvars} variable {variables[0]}"
        else:
            s += f"{nvars} variables {variables}"
        if nrels == 1:
            s += " with 1 relation"
        else:
            s += f" with {len(self.relations)} relations"
        return s
