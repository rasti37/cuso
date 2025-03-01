"""BoundSet class"""

from collections import UserDict
from typing import Dict

from sage.all import Expression, Integer
from sage.all import Polynomial as SagePolynomial
from sage.rings.polynomial.multi_polynomial import MPolynomial

from cuso.data.solutions import Solution
from cuso.data.types import Polynomial
from .bound import Bound


class BoundSet(UserDict):
    """Represent a collection of upper and lower bounds on variables."""

    def get_lower_bound(self, expr: Expression) -> int:
        """Return the lower bound for the given expression.

        Args:
            expr (Expression): Expression to evaluate

        Raises:
            ValueError: The lower bound is not defined for one of the symbols

        Returns:
            int: Integer lower bound for the expression
        """
        if isinstance(expr, int):
            return expr
        if expr in self:
            lbound = self[expr].lower
            if lbound is None:
                raise ValueError(f"Lower bound not set for {expr}")
            return lbound
        lower_bounds = {x: bound.lower for x, bound in self.items()}
        try:
            if isinstance(expr, Expression):
                return int(expr.subs(lower_bounds))
            return -self.get_poly_max_bound(-expr)
        except TypeError as exc:
            raise ValueError(
                f"Could not find lower bound for {expr}. Are all lower bounds specified?"
            ) from exc

    def get_upper_bound(self, expr: Expression) -> int:
        """Return the upper bound for the given expression.

        Args:
            expr (Expression): Expression to evaluate

        Raises:
            ValueError: The upper bound is not defined for one of the symbols

        Returns:
            int: Integer upper bound for the expression
        """
        if isinstance(expr, int):
            return expr
        if expr in self:
            ubound = self[expr].upper
            if ubound is None:
                raise ValueError(f"Upper bound not set for {expr}")
            return ubound
        upper_bounds = {x: bound.upper for x, bound in self.items()}
        try:
            if isinstance(expr, Expression):
                return int(expr.subs(upper_bounds))
            return self.get_poly_max_bound(expr)
        except TypeError as exc:
            err_str = f"Could not find upper bound for {expr}. Are all upper bounds specified?"
            raise ValueError(err_str) from exc

    def get_abs_bound(self, expr: Expression) -> int:
        """Return the upper bound for the absolute value of the given expression.

        Args:
            expr (Expression): Expression to evaluate

        Raises:
            ValueError: The absolute value bound is not defined for one of the symbols

        Returns:
            int: Integer upper bound for the absolute value of the expression
        """
        upper_bound = self.get_upper_bound(expr)
        lower_bound = self.get_lower_bound(expr)
        return max(abs(upper_bound), abs(lower_bound))

    def get_poly_max_bound(self, poly: Polynomial) -> int:
        """Return a bound on the absolute value of the evaluation of a polynomial.

        Args:
            poly (Polynomial): Polynomial being evaluated

        Returns:
            int: maximum absolute value of the polynomial within the bounds
        """
        maxval = 0
        ring = poly.parent()
        for ci, mi in zip(poly.coefficients(), poly.monomials()):
            maxabs = [max(map(abs, self[xi])) for xi in ring.gens()]
            maxterm = abs(ci) * int(mi(*maxabs))
            maxval += maxterm
        return int(maxval)

    def check(self, solution: Solution) -> bool:
        """Check whether the solution satisfies the bounds.

        Args:
            solution (Dict): Dictionary of values

        Returns:
            bool: All supplied values match bounds.
        """
        for value in solution.values():
            if not isinstance(value, (int, Integer)):
                raise TypeError(f"Solution value {value} is not integer.")

        for xi, bound in self.items():
            if xi not in solution:
                continue
            val = solution[xi]
            if val not in bound:
                return False
        return True

    def __setitem__(self, key, value):
        # Check keys
        type_error_s = "Keys must be symbols or generators of a polynomial ring"
        if isinstance(key, (MPolynomial, SagePolynomial)):
            if isinstance(key, MPolynomial) and not key.is_generator():
                raise TypeError(type_error_s)
            if isinstance(key, SagePolynomial) and not key.is_gen():
                raise TypeError(type_error_s)
        elif isinstance(key, Expression):
            if not key.is_symbol():
                raise TypeError(type_error_s)
        else:
            raise TypeError(type_error_s)
        if not isinstance(value, Bound):
            value = Bound(key, *value)
        super().__setitem__(key, value)

    def __repr__(self):
        s = "BoundSet for "
        xs = tuple(self.keys())
        if len(xs) == 1:
            s += str(xs[0])
        else:
            s += str(xs)

        return s

    def __str__(self):
        s = "Multivariate Coppersmith Bounds(\n"
        for bnd in self.values():
            s += "\t" + str(bnd) + ",\n"
        s += ")"
        return s
