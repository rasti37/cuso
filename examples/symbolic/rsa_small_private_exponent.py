from sage.all import PolynomialRing, QQ

from common import print_results
from cuso.symbolic import (
    SymbolicCoppersmithProblem,
    SymbolicBounds,
    convert_with_unraveled_linearization,
    get_asymptotic_bounds,
)


def run_symbolic():
    CoeffField = PolynomialRing(QQ, "a").fraction_field()
    R = PolynomialRing(CoeffField, "x,y")

    (a,) = CoeffField.gens()
    x, y = R.gens()

    # List of relations that hold modulo p
    mod_rels = [
        x * y + a * x - 1,
    ]
    ul = [
        x * y - 1,
    ]

    bounds = [
        SymbolicBounds(0, 1),  # delta
        SymbolicBounds(QQ(1) / 2, 0),  # 1/2
    ]
    bounds_guess = [0.1, 0.5]

    prob = SymbolicCoppersmithProblem(mod_rels, bounds)
    (x, y, u), prob, bounds_guess = convert_with_unraveled_linearization(prob, ul, bounds_guess)

    monomial_vertices = [
        1,
        x**3,
        u**3,
        u**3 * y,
    ]
    precompute_multiplicity = 3
    tshifts_to_include = [True,True,True]

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
