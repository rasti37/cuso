from sage.all import PolynomialRing

from .problem import SymbolicCoppersmithProblem, SymbolicBounds


def convert_with_unraveled_linearization(
    prob: SymbolicCoppersmithProblem, unrav_lins, bounds_guess
):
    R = prob.ring
    if len(unrav_lins) == 0:
        return R.gens(), prob, bounds_guess
    assert max(prob.rel_mod_deg) == 1

    u_vars = [f"u_{i}" for i in range(1, len(unrav_lins) + 1)]
    new_vars = (*R.variable_names(), *u_vars)
    new_R = PolynomialRing(R.base_ring(), new_vars)

    mod_rels = [new_R(f) for f in prob.mod_rels]
    int_rels = [new_R(f) for f in prob.int_rels]
    u_s = new_R.gens()[-len(u_vars) :]
    for i in range(len(u_s)):
        int_rels += [u_s[i] - unrav_lins[i]]

    new_bounds = prob.bounds[:]
    new_bounds_guess = bounds_guess[:]
    for i, ul_i in enumerate(unrav_lins):
        v, f = 0, 0
        b = 0
        max_term = max(
            ul_i.monomials(),
            key=lambda x: sum(bgi * ei for bgi, ei in zip(bounds_guess, x.degrees())),
        )
        degs = max_term.degrees()
        for j, sj in enumerate(degs):
            v += sj * prob.bounds[j].constant_term
            f += sj * prob.bounds[j].linear_term
            b += sj * bounds_guess[j]
        new_bounds += [SymbolicBounds(v, f)]
        new_bounds_guess += [b]
    new_prob = SymbolicCoppersmithProblem(
        mod_rels=mod_rels,
        bounds=new_bounds,
        int_rels=int_rels,
        modulus_multiple_name=prob.modulus_multiple_name,
        modulus_name=prob.modulus_name,
        divisor_bound=prob.divisor_bound,
    )
    return new_R.gens(), new_prob, new_bounds_guess
