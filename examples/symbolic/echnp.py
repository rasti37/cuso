from functools import reduce
import operator

from sage.all import PolynomialRing, QQ

from common import print_results
from cuso.symbolic import SymbolicCoppersmithProblem, get_asymptotic_bounds

def echnp(num_samples=3):
    # unknowns alpha, eps_i
    # known p, x_i, b_i
    known_list = ",".join(
        ["h0", "a", "b"] + [f"h_{i},xQ_{i}" for i in range(num_samples)]
    )
    CoeffField = PolynomialRing(QQ, known_list).fraction_field()
    knowns = CoeffField.gens()
    h0, a, b, *_ = knowns
    h_s = [knowns[i] for i in range(3, len(knowns), 2)]
    xQ_s = [knowns[i] for i in range(4, len(knowns), 2)]

    unknown_list = ",".join(
        ["x_0"] + [f"y_{i}" for i in range(num_samples)]
    )
    R = PolynomialRing(CoeffField, unknown_list)
    x_0, *ys = R.gens()

    # List of relations that hold modulo p
    # From XSWH22 section 2.4
    mod_rels = []
    for i in range(num_samples):
        y_i = ys[i]
        hi = h_s[i]
        xQi = xQ_s[i]

        A_i = hi*(h0 - xQi)**2 - 2*h0**2*xQi - 2*(a + xQi**2)*h0 - 2*a*xQi - 4*b
        B_i = 2*(hi*(h0 - xQi) - 2*h0*xQi - a - xQi**2)
        C_i = hi - 2*xQi
        D_i = (h0 - xQi)**2
        E_i = 2*(h0 - xQi)
        mod_rel = A_i + B_i*x_0 + C_i*x_0**2 + D_i*y_i + E_i*x_0*y_i + x_0**2*y_i
        mod_rels += [
            mod_rel
        ]

    bounds = [1] * (num_samples + 1)
    prob = SymbolicCoppersmithProblem(mod_rels, bounds, modulus_name="p")
    bounds_guess = [0.1] * (num_samples + 1)

    precompute_multiplicity = 1
    monomial_vertices = [1, x_0**2] + [y_i**2 for y_i in ys] + [x_0**2*y_i for y_i in ys]
    monomial_vertices += [reduce(operator.mul, ys, 1)]
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
    echnp(num_samples=3)


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)
    run_symbolic()
