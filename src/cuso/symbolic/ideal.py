import functools
import itertools
import logging

from sage.all import TermOrder, PolynomialRing

from .problem import SymbolicCoppersmithProblem

logger = logging.getLogger("cuso.symbolic.Ideal")


"""
Compute the asymptotic Coppersmith bound symbolically
"""
class IdealFinder:
    def __init__(self, Js, J_inf):
        self.Js = Js
        self.J_inf = J_inf

    @functools.lru_cache()
    def get_ideal(self, multiplicity):
        if min(multiplicity) < 0:
            # Invalid ideal
            return None
            
        if sum(multiplicity) == 0:
            R = self.J_inf.ring()
            return R.ideal([1])

        J = self.J_inf

        # For each of the ideals J1 in mod_rels, see if there's a previous
        # ideal J2 such that the modulus of J1 * J2 is a multiple of the
        # modulus of J and J2 is "smaller" than this multiplicity
        mod_rels = self.Js
        for mults, J1 in mod_rels.items():
            diffgen = itertools.product(*[range(d+1) for d in mults])
            for diff in diffgen:
                if diff == mults:
                    continue

                J2_mults = tuple(ai - bi + di for ai, bi, di in zip(multiplicity, mults, diff))
                J2 = self.get_ideal(J2_mults)
                if J2 is None:
                    continue
                
                J_prod = J1 * J2
                J = J + J_prod

        #lg_mod = self.get_modulus_lbound_size(multiplicity)
        logger.debug("Generated ideal for multiplicity %d", multiplicity)
        return J

def get_base_ideal(coppersmith_problem: SymbolicCoppersmithProblem, bounds_guess):
    # Given a multivariate Coppersmith problem and a guess of the
    # optimal bounds (so we are working with a useful monomial order),
    # return the ideals J_p and J_inf that represent the input polynomials
    # Also return p_weight, the weight of p w.r.t. the term order.
    # We would prefer to normalize so p_weight = 1., but TermOrder does
    # not allow noninteger weights, so we scale and round.

    base_R = coppersmith_problem.mod_rels[0].parent()
    coeffRing = base_R.base_ring()
    xs = base_R.variable_names()

    # Introduce the multiple of the modulus as a new variable
    new_xs = (coppersmith_problem.modulus_multiple_name, *xs)

    # Get the monomial order
    weights = (1 / coppersmith_problem.divisor_bound, *bounds_guess)
    p_weight = 1000
    weights_as_ints = tuple(round(w * p_weight) for w in weights)
    order = TermOrder("wdeglex", weights_as_ints)

    # Define the new ring
    new_R = PolynomialRing(coeffRing, new_xs, order=order)

    # Get all polynomials that are 0 mod p in this new ring
    N = new_R.gens()[0]
    new_mod_rels = [new_R(f) for f in coppersmith_problem.mod_rels] + [N]
    new_rel_mod_degs = coppersmith_problem.rel_mod_deg + [coppersmith_problem.mod_mult_modulus_degree]
    new_int_rels = [new_R(f) for f in coppersmith_problem.int_rels]

    Js = {(e,): new_R.ideal([f for f, p_deg in zip(new_mod_rels, new_rel_mod_degs) if p_deg == e]) for e in set(new_rel_mod_degs)}
    J_inf = new_R.ideal(new_int_rels)
    return Js, J_inf, p_weight


def get_ideal(
    coppersmith_problem: SymbolicCoppersmithProblem,
    bounds_guess,
    precompute_multiplicity,
):
    logger.debug("Computing ideals J_p and J_inf")
    Js, J_inf, p_weight = get_base_ideal(coppersmith_problem, bounds_guess)

    # Get the ideal we are searching for symbolic polynomials in
    logger.debug("Precomputing to multiplicity k_pre = %d", precompute_multiplicity)
    logger.debug("Computing ideal J = J_p**k_pre + J_inf")
    ideal_finder = IdealFinder(Js, J_inf)
    J = ideal_finder.get_ideal((precompute_multiplicity,))
    return J, J_inf
