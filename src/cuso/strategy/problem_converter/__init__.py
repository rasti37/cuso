"""This module implements problem converter strategies.

Given a multivariate Coppersmith problem, return a new
multivariate Coppersmith problem such that solving the
new problem allows us to solve the original problem.
This concept applies to variable recentering (so solutions
are bounded by (-B, B) instead of (0, 2B), for example)
and unraveled linearization.
"""

from .problem_converter import MultivariateProblemConverter

from .chain import ChainConverter
from .monom_ordering import BoundedMonomialOrderConverter
from .recenter import RecenterConverter
from .unraveled_linearization import UnraveledLinearization
