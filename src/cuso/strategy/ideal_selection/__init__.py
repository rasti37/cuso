"""This module implements ideal selection strategies.

Given a multivariate Coppersmith problem, generate different ideals
such that shift polynomials in ideal J_i have modulus N_i. Multiple
ideals are returned with increasing modulus, so if a problem is not
solved with low multiplicity, we can try it with high multiplicity.
"""

from .ideal_selection import RelationIdealGenerator
