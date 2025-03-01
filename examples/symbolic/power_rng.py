from sage.all import PolynomialRing, QQ

from common import print_results
from cuso.symbolic import (
    SymbolicCoppersmithProblem,
    SymbolicBounds,
    convert_with_unraveled_linearization,
    get_asymptotic_bounds,
)

def power_rng(num_samples=2):
    nrels = num_samples - 1
    a_names = [f"a_{i}" for i in range(nrels)]
    b_names = [f"b_{i}" for i in range(nrels)]
    all_names = ",".join(a_names + b_names)
    CoeffField = PolynomialRing(QQ, all_names).fraction_field()
    known_vars = CoeffField.gens()
    a_s = known_vars[:nrels]
    b_s = known_vars[nrels:]
    
    unk_names = ",".join([f"x_{i}" for i in range(num_samples)])
    R = PolynomialRing(CoeffField, unk_names)

    x_s = R.gens()

    # List of relations that hold modulo p
    mod_rels = []
    ul = []
    for i in range(nrels):
        mod_rel = x_s[i]**2 + a_s[i]*x_s[i] + b_s[i] - x_s[i+1]
        mod_rels += [mod_rel]
        ul += [x_s[i]**2 - x_s[i+1]]

    bounds = [1] * num_samples
    bounds_guess = [0.1] * num_samples
    prob = SymbolicCoppersmithProblem(mod_rels, bounds)
    new_vars, prob, bounds_guess = convert_with_unraveled_linearization(prob, ul, bounds_guess)

    monom_vertices = None
    new_xs = new_vars[:num_samples]
    new_us = new_vars[num_samples:]
    monom_vertices = [1]
    for i in range(num_samples):
        monom_vertices += [new_xs[i]**(2**(num_samples - i - 1))]
    for i in range(nrels):
        monom_vertices += [new_us[i]**(2**(nrels - i))]

    if num_samples == 2:
        k = 2
    elif num_samples == 3:
        k = 4
    else:
        raise

    monomial_vertices = monom_vertices
    precompute_multiplicity = k
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

def run_symbolic():
    power_rng(2)
    power_rng(3)


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)
    run_symbolic()
