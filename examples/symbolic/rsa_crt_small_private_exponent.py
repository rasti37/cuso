from sage.all import PolynomialRing, QQ

from common import print_results
from cuso.symbolic import (
    SymbolicCoppersmithProblem,
    SymbolicBounds,
    convert_with_unraveled_linearization,
    get_asymptotic_bounds,
)

def run_symbolic():
    CoeffField = PolynomialRing(QQ, "N").fraction_field()
    R = PolynomialRing(CoeffField, "p, q, k, l")

    N, = CoeffField.gens()
    p,q,k,l = R.gens()

    # List of relations that hold modulo p
    mod_rels = [
        -1 - k * (p - 1),
        -1 - l * (q - 1),
    ]
    int_rels = [
        p * q - N,
    ]
    ul = [
        p + q,
        k * p + l * q - k - l + 2,
        k + l - 1,
        k * l,
        k * l * p + k * l * q - k * p - l * q + k + l - 1,
    ]

    bounds = [
        SymbolicBounds(QQ(1)/2, 0),
        SymbolicBounds(QQ(1)/2, 0),
        SymbolicBounds(QQ(1)/2, 1),
        SymbolicBounds(QQ(1)/2, 1),
    ]
    bounds_guess = [0.5, 0.5, 0.65, 0.65]

    prob = SymbolicCoppersmithProblem(mod_rels, bounds, int_rels=int_rels, modulus_name="e")
    new_vars, prob, bounds_guess = convert_with_unraveled_linearization(prob, ul, bounds_guess)

    p, q, k, l, u, v, w, x, y = new_vars
    monomial_vertices = [1, u, v**2, w**2, x**2, y**2, u*w**2, u*y**2, v*w, v*x, v*y]
    precompute_multiplicity = 4
    tshifts_to_include = [False] * len(new_vars)

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


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)
    run_symbolic()
