import itertools
import logging
from typing import List
import warnings

from sage.all import PolynomialRing, QQ

from .shift_set_properties import ShiftProperties

logger = logging.getLogger("cuso.symbolic.Polynomial")

def get_eval_points(dim, tshifts_to_include):
    multiplicities = list(range(1, dim + 4))

    max_ts = tuple(2 if t_i else 0 for t_i in tshifts_to_include)

    eval_pts = []
    for k in multiplicities:
        for ts in itertools.product(*[range(0, max_t+1) for max_t in max_ts]):
            eval_pts += [(k, ts)]
    return eval_pts

def get_polynomial(xs, ys, tshifts_to_include):
    nvars = 1 + sum(tshifts_to_include)
    if nvars == 1:
        R = PolynomialRing(QQ, "k")
        xs = [x[0] for x in xs]
        points = list(zip(xs, ys))
        p = R.lagrange_polynomial(points)
        return p
    
    # Multivariate
    t_names = [f"t_{i}" for i in range(len(tshifts_to_include)) if tshifts_to_include[i]]
    varnames = ("k", *t_names)
    R = PolynomialRing(QQ, varnames)
    try:
        p = R.interpolation(nvars + 2, xs, ys)
    except AttributeError as exc:
        warnings.warn(
            "Multivariate polynomial interpolation failed. It's likely that you need to "
            "upgrade SageMath to a newer version (9.8+) that has this feature."
            UserWarning
        )
        raise RuntimeError("SageMath version is too old.")

    return p

def get_polynomials(S_prop: ShiftProperties, tshifts_to_include: List[bool]):
    dim = S_prop.dim()
    points = get_eval_points(dim, tshifts_to_include)

    evaluations = []
    xs = []
    logger.debug("Querying %d points", len(points))
    for k, ts in points:
        logger.debug("Determine lattice construction for k = %d, t = %s", k, ts)
        evaluation = S_prop.s_ktul(k, ts)
        evaluations += [evaluation]

    xs += [(k, *ts) for k, ts in points]
    y_dim = [y[0] for y in evaluations]
    y_Xs  = list(zip(*[y[1] for y in evaluations]))
    y_C  = [y[2] for y in evaluations]

    logger.debug("Interpolating polynomials")
    s_dim = get_polynomial(xs, y_dim, tshifts_to_include)
    s_xs = [
        get_polynomial(xs, y_Xi, tshifts_to_include) for y_Xi in y_Xs
    ]
    s_C = get_polynomial(xs, y_C, tshifts_to_include)

    if any(tshifts_to_include):
        if s_dim.total_degree() > dim:
            s_dim = None
        s_xs = [
            None if s_xi.total_degree() > dim + 1 else s_xi for s_xi in s_xs
        ]
        if s_C.total_degree() > dim + 1:
            s_C = None
    else:
        if s_dim.degree() > dim:
            s_dim = None
        s_xs = [
            None if s_xi.degree() > dim + 1 else s_xi for s_xi in s_xs
        ]
        if s_C.degree() > dim + 1:
            s_C = None

    if s_dim is None or any([s_xi is None for s_xi in s_xs]) or s_C is None:
        raise ValueError("Polynomials do not have expected total degree")

    return s_dim, s_xs, s_C