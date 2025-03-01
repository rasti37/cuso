from sage.all import PolynomialRing, QQ

from common import print_results
from cuso.symbolic import SymbolicCoppersmithProblem, get_asymptotic_bounds


def run_symbolic():
    CoeffField = PolynomialRing(QQ, "a, c").fraction_field()
    R = PolynomialRing(CoeffField, "x")

    a, c = CoeffField.gens()
    (x,) = R.gens()

    # List of relations that hold modulo p
    mod_rels = [
        (a + x) ** 3 - c,
    ]
    bounds = [1]
    prob = SymbolicCoppersmithProblem(mod_rels, bounds)
    bounds_guess = [0.1]
    monomial_vertices = [1, x**3]
    precompute_multiplicity = 1
    tshifts_to_include = [False]

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
