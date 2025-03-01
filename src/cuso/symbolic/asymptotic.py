import logging

from sage.all import lcm

from .ideal import get_ideal
from .limit import get_delta_star
from .optimal import get_optimal_shift_polys
from .polytope import get_monomials_from_vertices
from .polynomial import get_polynomials
from .problem import SymbolicCoppersmithProblem
from .shift_set_properties import ShiftProperties
from .utils import lc_no_N

logger = logging.getLogger("cuso.symbolic.Asymptotic")


def get_asymptotic_bounds(
    coppersmith_problem: SymbolicCoppersmithProblem,
    bounds_guess,
    monomial_vertices=None,
    precompute_multiplicity=1,
    tshifts_to_include=None,
    skip_reduce=False,
):
    logger.info("Getting asymptotic bounds for multivariate Coppersmith problem")
    for line in str(coppersmith_problem).split("\n"):
        logger.debug(line)

    # Get the ideal we are searching for symbolic polynomials in
    logger.info("Computing J = J^k_pre + J_inf")
    J, J_inf = get_ideal(coppersmith_problem, bounds_guess, precompute_multiplicity)
    G_inf = J_inf.groebner_basis()
    # Check if G_inf is empty
    if G_inf == [0]:
        G_inf = []

    logger.info("Computing M_1 from M_vert")
    M_1 = get_monomials_from_vertices(J.ring(), monomial_vertices)

    # Get optimal shift polynomials for these monomials
    logger.info("Getting optimal shift polynomials for M_1")
    # Need to skip reduction because otherwise sage exceeds recursion depth in template evaluation
    S_1 = get_optimal_shift_polys(M_1, J, skip_reduce=skip_reduce)

    S_props = ShiftProperties(S_1, G_inf)

    logger.info("Recovering polynomials.")
    nvars = len(J.ring().gens()[1:])
    if tshifts_to_include is None:
        tshifts_to_include = [True] * nvars
    elif tshifts_to_include in [True, False]:
        tshifts_to_include = [tshifts_to_include] * nvars
    s_dim, s_xs, s_C = get_polynomials(S_props, tshifts_to_include)

    delta_star, tau = get_delta_star(
        coppersmith_problem,
        precompute_multiplicity,
        s_dim,
        s_xs,
        s_C,
        tshifts_to_include,
    )

    # The following must be coprime with N for the precomputation substitutions to be valid:
    g = 1
    coprime_condition = []
    for f in S_1:
        coef_denoms = [c.denominator() for c in f.coefficients()]
        for c in coef_denoms:
            if c != 1 and c not in coprime_condition:
                coprime_condition += [c]
            g = lcm(g, c)

    return delta_star, tau, coprime_condition
