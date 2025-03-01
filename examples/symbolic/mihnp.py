from functools import reduce
import itertools
import operator

from sage.all import PolynomialRing, QQ

from common import print_results
from cuso.symbolic import SymbolicCoppersmithProblem, SymbolicBounds, get_asymptotic_bounds

def mihnp(num_samples=3, simplified=True):
    # unknowns alpha, eps_i
    # known p, x_i, b_i
    if simplified:
        known_list = ",".join(
            [f"c_{i}" for i in range(num_samples)] + [f"b_{i}" for i in range(num_samples)] + [f"d_{i}" for i in range(num_samples)]
        )
    else:
        known_list = ",".join(
            [f"c_{i}" for i in range(num_samples)] + [f"b_{i}" for i in range(num_samples)]
        )
    CoeffField = PolynomialRing(QQ, known_list).fraction_field()
    knowns = CoeffField.gens()
    c_s = knowns[:num_samples]
    b_s = knowns[num_samples:2*num_samples]
    if simplified:
        d_s = knowns[2*num_samples:]

    unknown_list = ",".join(
        ["alpha"] + [f"eps_{i}" for i in range(num_samples)]
    )
    R = PolynomialRing(CoeffField, unknown_list)
    alpha, *eps = R.gens()


    # List of relations that hold modulo p
    mod_rels = []
    for i in range(num_samples):
        if simplified:
            mod_rel = alpha * eps[i] + c_s[i] * eps[i] + b_s[i] * alpha + d_s[i]
        else:
            mod_rel = (alpha + c_s[i]) * (eps[i] + b_s[i]) - 1
        mod_rels += [
            mod_rel
        ]

    bounds = [SymbolicBounds(1, 0)] + [SymbolicBounds(0, 1)] * num_samples
    prob = SymbolicCoppersmithProblem(mod_rels, bounds, modulus_name="p")
    bounds_guess = [1] + [0.1] * num_samples
    precompute_multiplicity = 2
    if num_samples in [4, 5]:
        power = 1
    else:
        power = 2
    cube = itertools.product(*[[1,xi**power] for xi in eps])
    monomial_vertices = [reduce(operator.mul, vert, 1) for vert in cube]
    tshifts_to_include = [False] * (num_samples + 1)

    delta, tau, coprime_cond = get_asymptotic_bounds(
        prob,
        bounds_guess,
        monomial_vertices,
        precompute_multiplicity=precompute_multiplicity,
        tshifts_to_include=tshifts_to_include,
    )
    print_results(
        prob,
        bounds_guess,
        monomial_vertices,
        precompute_multiplicity,
        tshifts_to_include,
        delta,
        tau,
        coprime_cond,
    )

def run_symbolic():
    mihnp(num_samples=3, simplified=False)
    mihnp(num_samples=4, simplified=False)
    mihnp(num_samples=5, simplified=False)


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)
    run_symbolic()
