from sage.all import PolynomialRing, QQ, TermOrder

from common import print_results
from cuso.symbolic import (
    SymbolicCoppersmithProblem,
    SymbolicBounds,
    get_asymptotic_bounds,
)


def run_symbolic():
    CoeffField = PolynomialRing(QQ, "E, N").fraction_field()

    gamma = QQ(45) / 100
    R = PolynomialRing(CoeffField, "x,y")
    E, N = CoeffField.gens()
    x, y = R.gens()

    mod_rels = [
        E - x,
        N - y,
    ]
    rel_mod_deg = [
        1, # E - x is 0 mod g
        2, # N - y is 0 mod g^2
    ]
    bounds = [SymbolicBounds(0, 1), SymbolicBounds(1 / QQ(2) / gamma, 0)]
    # N minus 1 (Nm1) is multiple of modulus squared
    prob = SymbolicCoppersmithProblem(
        mod_rels,
        bounds,
        relation_modulus_degrees=rel_mod_deg,
        modulus_name="g",
        modulus_multiple_name="Nm1",
        mod_mult_modulus_degree=1, # N - 1 is 0 mod g
        divisor_bound=gamma,
    )

    bounds_guess = [0.1, 0.5]
    monomial_vertices = [1, x**7, y**5]
    precompute_multiplicity = 6
    tshifts_to_include = [True, True]

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
