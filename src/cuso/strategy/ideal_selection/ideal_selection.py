"""Implement constructing ideals from relations"""

from typing import Iterator
import functools
import itertools
from math import log2

from sage.all import Integer, gcd, Expression

from cuso.strategy.strategy import Strategy
from cuso.data import RelationSet, BoundSet, RelationIdeal, Relation
from cuso.utils import weighted_combinations


class RelationIdealGenerator(Strategy):
    """Construct relation ideals for a given set of relations and bounds.

    As long as one of the relations has a modular constraint, we have at least
    one modular relation ideal. By taking products and sums of relation
    ideals, we get a new modular relation ideal. We yield relation ideals
    in increasing modulus size, corresponding to increasing multiplicity.
    """

    def __init__(self):
        super().__init__()
        self._factors = None
        self._J_ps = None
        self._J_inf = None
        self._ring = None

    def _get_list_of_factors(self, moduli):
        # Return list of pairwise coprime factors that divide all moduli
        # It would be more efficient to use a product tree, but this is easier.
        factors = []
        if 0 in moduli:
            factors += [0]

        # Get list of integers
        modset = []
        for mod in moduli:
            if isinstance(mod, (int, Integer)):
                modset += [mod]
            else:
                # Get leading coefficient and ensure the modulus
                # is a single term
                for p_i in mod.variables():
                    if p_i not in factors:
                        factors += [p_i]
                    assert len(mod.coefficients(p_i)) == 1
                int_part = mod.coefficients()[0][0]
                modset += [int_part]

        modset = list(set([mod for mod in modset if mod not in [0, 1]]))
        while True:
            # Check to see if all pairs have 0 gcd
            inds_shared_factor = None
            for i, mod_i in enumerate(modset):
                for j, mod_j in enumerate(modset):
                    if j <= i:
                        continue
                    if gcd(mod_i, mod_j) != 1:
                        inds_shared_factor = i, j
                        break
            if inds_shared_factor is None:
                break

            i, j = inds_shared_factor
            g = int(gcd(modset[i], modset[j]))
            new_modset = [g, modset[i] // g, modset[j] // g]
            new_modset += [
                m_i for i, m_i in enumerate(modset) if i not in inds_shared_factor
            ]
            new_modset = self._get_list_of_factors(new_modset)
            modset = new_modset

        return modset + factors

    def _get_base_ideals(self, relations, ring):
        # Cancel out any shared factors in the moduli
        new_rels = []
        for rel in relations:
            if rel.modulus is None:
                new_rels += [rel]
            else:
                g = functools.reduce(
                    gcd, rel.polynomial.coefficients() + [rel.modulus], 0
                )
                if g == 1:
                    new_rels += [rel]
                else:
                    ring = rel.ring()
                    new_rels += [Relation(rel.polynomial // g, rel.modulus // g)]
        relations = RelationSet(new_rels)

        # Get set of moduli
        moduli = [rel.modulus for rel in relations]
        # If there is no modulus, ignore
        moduli = [mod for mod in moduli if mod]
        moduli = list(set(moduli))
        factors = self._get_list_of_factors(moduli)

        # Get exponents
        exponents = []
        for mod in moduli:
            exponent = [0] * len(factors)
            for j, q_j in enumerate(factors):
                if isinstance(q_j, int):
                    while mod % q_j == 0:
                        # The modulus is divisible by this factor
                        mod //= q_j
                        exponent[j] += 1
                else:
                    if isinstance(mod, Expression):
                        exponent[j] = mod.degree(q_j)
                        mod /= q_j ** exponent[j]
            assert mod == 1
            exponents += [tuple(exponent)]

        # Split relations by their modulus
        rel_inf = []
        rel_mod = {exponent: [] for exponent in exponents}

        for rel in relations:
            mod = rel.modulus
            poly = rel.polynomial
            if mod is None:
                rel_inf += [poly]
            else:
                ind = moduli.index(mod)
                exp = exponents[ind]
                rel_mod[exp] += [poly]

        # Ensure the factors get their own ideal
        for i, q_i in enumerate(factors):
            if isinstance(q_i, int):
                exp = [0] * len(factors)
                exp[i] = 1
                exp = tuple(exp)
                if exp not in exponents:
                    exponents += [exp]
                    moduli += [q_i]
                    rel_mod[exp] = []
                rel_mod[exp] += [ring(q_i)]

        J_inf = RelationIdeal(rel_inf, ring, None)
        J_ps = {
            exp: RelationIdeal(rel_mod[exp], ring, mod)
            for exp, mod in zip(exponents, moduli)
        }
        return factors, J_ps, J_inf

    @functools.lru_cache()
    def _get_ideal(self, multiplicity):
        if min(multiplicity) < 0:
            # Invalid ideal
            return None

        if sum(multiplicity) == 0:
            return RelationIdeal([1], self._ring, 1)

        J = self._J_inf

        # For each of the ideals J1 in mod_rels, see if there's a previous
        # ideal J2 such that the modulus of J1 * J2 is a multiple of the
        # modulus of J and J2 is "smaller" than this multiplicity
        mod_ideals = self._J_ps
        for exp, ideal1 in mod_ideals.items():
            diffgen = itertools.product(*[range(d + 1) for d in exp])
            applied_diffs = []
            for diff in diffgen:
                if diff == exp:
                    continue

                is_dominated = False
                for applied_diff in applied_diffs:
                    # Check to see if the applied_diff dominates this diff.
                    # If so, there's no use using this diff.
                    if all(ai <= bi for ai, bi in zip(applied_diff, diff)):
                        is_dominated = True
                if is_dominated:
                    break

                exp_smaller = tuple(
                    ai - bi + di for ai, bi, di in zip(multiplicity, exp, diff)
                )
                ideal2 = self._get_ideal(exp_smaller)
                if ideal2 is None:
                    continue

                applied_diffs += [diff]
                J = J + ideal1 * ideal2

        self.logger.debug("Generated ideal for multiplicity %s", multiplicity)
        return J

    def run(self, relations: RelationSet, bounds: BoundSet) -> Iterator[RelationIdeal]:
        """Generate relation ideals for input relations.

        Args:
            relations (RelationSet): input relations
            bounds (BoundSet): bounds on the variables

        Yields:
            Iterator[RelationIdeal]: RelationIdeals to find shift relations in
        """
        ring = relations.ring()
        factors, J_ps, J_inf = self._get_base_ideals(relations, ring)

        self.logger.info("Recovered %u modular ideal(s)", len(J_ps))

        self._factors = factors
        self._J_ps = J_ps
        self._J_inf = J_inf
        self._ring = J_inf.ring()

        if len(J_ps) == 0:
            # Only have integer relations.
            yield J_inf, J_inf
            return

        weights = [log2(bounds.get_lower_bound(q)) for q in self._factors]
        for exp, _ in weighted_combinations(weights):
            if sum(exp) == 0:
                continue
            J = self._get_ideal(exp)
            mod_bound = bounds.get_lower_bound(J.modulus)
            s_exp = str(exp[0]) if len(exp) == 1 else str(exp)
            self.logger.info(
                "Generated ideal of multiplicity %s with %f-bit modulus",
                s_exp,
                log2(mod_bound),
            )

            yield J, J_inf
