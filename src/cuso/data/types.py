"""Type hints for cuso objects"""

from typing import Union

from sage.all import Expression
from sage.all import Polynomial as SagePolynomial
from sage.rings.polynomial.multi_polynomial import MPolynomial

__all__ = ["Modulus, Monomial, Polynomial, Variable"]

Modulus = Union[int, Expression]
Polynomial = Union[SagePolynomial, MPolynomial]
Monomial = Polynomial
Variable = Union[str, Expression]
