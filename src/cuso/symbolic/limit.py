import logging

from sage.all import QQ, RR, PolynomialRing, minimize_constrained

from .problem import SymbolicCoppersmithProblem

logger = logging.getLogger("cuso.symbolic.Limit")


def limit_univariate(to_optimize):
    num = to_optimize.numerator()
    den = to_optimize.denominator()

    # Limit is the same if we discard low order terms
    lim_num = num.lt()
    lim_den = den.lt()

    if lim_num.degree() != lim_den.degree():
        logger.warning("Limit is either zero or infinite")
        raise ValueError()
    
    return QQ(lim_num / lim_den)


def limit_multivariate(to_optimize):
    """We are given a multivariate expression in terms of k, t_1, ...

    Given a choice of tau, let t_i = tau_i.
    We can compute the limit for that choice of tau as k -> infinity.
    Derive an expression for tau |-> lim_{k->infty} to_optimize(k, k*tau)
    Then optimize tau
    """
    # Replace t_i with tau_i * k
    ring = to_optimize.parent()
    k, *ts = ring.gens()
    tau_names = [str(t_i).replace("t_", "tau_") for t_i in ts]
    Qt = PolynomialRing(QQ, ",".join(tau_names)).fraction_field()
    taus = Qt.gens()
    new_ring = PolynomialRing(Qt, "k").fraction_field()

    (k_new,) = new_ring.gens()
    subs = {t_i: tau_i * k_new for t_i, tau_i in zip(ts, taus)}
    subs[k] = k_new

    to_optimize = new_ring(to_optimize.subs(subs))

    num = to_optimize.numerator()
    den = to_optimize.denominator()

    # Limit is the same if we discard low order terms
    lim_num = num.lt()
    lim_den = den.lt()

    if lim_num.degree() != lim_den.degree():
        logger.warning("Limit is either zero or infinite")
        raise ValueError()

    lim = Qt(lim_num / lim_den)
    Rt = PolynomialRing(RR, ",".join(tau_names)).fraction_field()
    lim = Rt(lim)

    # maximize lim
    # Maximizing upper bound by minimizing f
    f = lambda p: -float(lim(*p))
    # All values of tau are >= 0
    constraints = [
        (0, 1000) for _ in range(len(taus))
    ]
    logger.debug("Trying to maximize %s", lim)
    maximizer = minimize_constrained(f, constraints, [0]*len(taus))
    logger.debug("Bound maximized by (%s) = (%s)", ",".join(map(str, taus)), ",".join(map(str, maximizer)))
    delta_star = -f(maximizer)
    return delta_star, maximizer


def get_delta_star(
    problem: SymbolicCoppersmithProblem,
    precompute_multiplicity,
    s_dim,
    s_xs,
    s_C,
    tshifts_to_include,
):
    """Compute expression 3 of Lemma 7 of eprint 2024/1577"""

    # We have to be careful here, since S_1 has shared root modulo p^k_pre
    # so the actual bound on the X_i values should be divided by k_pre
    k_pre = QQ(precompute_multiplicity)
    ab_s = [
        (bound.constant_term / k_pre, bound.linear_term / k_pre)
        for bound in problem.bounds
    ]
    # Similarly, C_j will be N, so s_C reports the number of copies of N.
    # thus log_{p^k_pre} N = (log N) / (log p) / k_pre
    log_p_C = (1 / QQ(problem.divisor_bound)) / k_pre

    ring = s_dim.parent()
    k = ring.gens()[0]
    numerator = k * s_dim
    for (a_i, _), s_xi in zip(ab_s, s_xs):
        numerator -= a_i * s_xi
    numerator -= log_p_C * s_C

    denominator = sum(b_i * s_xi for (_, b_i), s_xi in zip(ab_s, s_xs))

    to_optimize = numerator / denominator

    if sum(tshifts_to_include) == 0:
        delta_star = limit_univariate(to_optimize)
        taus = (0,) * len(tshifts_to_include)
        return delta_star, taus

    delta_star, tau_subset = limit_multivariate(to_optimize)
    taus = [0] * len(tshifts_to_include)
    i = 0
    for j, tshift_included in enumerate(tshifts_to_include):
        if tshift_included:
            taus[j] = tau_subset[i]
            i += 1

    return delta_star, tuple(taus)
