from sage.all import PolynomialRing, QQ

from common import print_results
from cuso.symbolic import SymbolicCoppersmithProblem, get_asymptotic_bounds


def run_symbolic():
    CoeffField = PolynomialRing(QQ, "A_msb, B_msb, C_msb").fraction_field()
    R = PolynomialRing(CoeffField, "x,y,z")

    A_msb, B_msb, C_msb = CoeffField.gens()
    x,y,z = R.gens()

    # List of relations that hold modulo p
    mod_rels = [
        (A_msb+x)**2 + 12*(A_msb+x) - 4*(A_msb+x)*(B_msb+y)**2 - 8*(B_msb+y)**2 + 36,
        (C_msb+z)**2 + 12*(C_msb+z) - 4*(C_msb+z)*(A_msb+x)**2 - 8*(A_msb+x)**2 + 36,
    ]
    bounds = [1, 1, 1]
    bounds_guess = [0.1, 0.1, 0.1]
    prob = SymbolicCoppersmithProblem(mod_rels, bounds, modulus_name="p")

    monomial_vertices = [
        1, x**4, x**4*z**2, z**4,
        y**4, x**2*y**4, x**2*y**4*z**2, y**4*z**2,
    ]
    tshifts_to_include = [False, False, False]
    precompute_multiplicity = 2

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
