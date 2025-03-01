from sage.all import PolynomialRing, QQ

from common import print_results
from cuso.symbolic import SymbolicCoppersmithProblem, get_asymptotic_bounds


def pacd(num_samples=1):
    known_list = ",".join([f"c_{i}" for i in range(num_samples)])
    CoeffField = PolynomialRing(QQ, known_list).fraction_field()
    knowns = CoeffField.gens()
    c_s = knowns

    unknown_list = ",".join([f"x_{i}" for i in range(num_samples)])
    R = PolynomialRing(CoeffField, unknown_list)
    xs = R.gens()

    # List of relations that hold modulo p
    mod_rels = []
    for i in range(num_samples):
        mod_rel = -xs[i] + c_s[i]
        mod_rels += [mod_rel]

    bounds = [1] * num_samples
    prob = SymbolicCoppersmithProblem(
        mod_rels,
        bounds,
        modulus_name="p",
        modulus_multiple_name="N",
        divisor_bound=QQ(4) / 10,
    )
    bounds_guess = [0.1] * num_samples
    if num_samples == 1:
        k, deg = 1, 1
    if num_samples == 2:
        k, deg = 3, 4
    if num_samples == 3:
        k, deg = 3, 4
    monomial_vertices = [1] + [xi**deg for xi in xs]
    precompute_multiplicity = k
    tshifts_to_include = [True] * num_samples

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
    pacd(1)
    pacd(2)
    pacd(3)


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)
    run_symbolic()
