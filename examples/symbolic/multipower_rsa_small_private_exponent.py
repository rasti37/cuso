from sage.all import PolynomialRing, QQ, TermOrder

from common import print_results
from cuso.symbolic import (
    SymbolicCoppersmithProblem,
    SymbolicBounds,
    get_asymptotic_bounds,
)


def run_symbolic():
    CoeffField = PolynomialRing(QQ, "e").fraction_field()

    r = 3
    p_len = QQ(1)/(r + 1)
    bounds_guess = [0.1]

    R = PolynomialRing(CoeffField, "x")
    e, = CoeffField.gens()
    x, = R.gens()
    
    mod_rels = [
        e*x - 1,
    ]
    rel_mod_deg = [
        2 # ex - 1 = 0 mod p^2
    ]

    bounds = [1]
    prob = SymbolicCoppersmithProblem(
        mod_rels,
        bounds,
        relation_modulus_degrees=rel_mod_deg,
        modulus_name="p",
        modulus_multiple_name="N",
        mod_mult_modulus_degree=3, # N is 0 mod p^3
        divisor_bound=p_len,
    )
    monomial_vertices = [1, x**4]
    precompute_multiplicity = 6
    tshifts_to_include = [True]

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
