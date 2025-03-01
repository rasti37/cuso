"""This file implements the optimal shift polynomial selection strategy."""

from typing import Iterator
import math

from sage.all import is_power_of_two

from cuso.data import Relation, RelationSet, BoundSet
from cuso.strategy.ideal_selection import RelationIdealGenerator
from cuso.utils import weighted_combinations

from .shift_poly_selection import ShiftPolyStrategy


class OptimalShiftPolys(ShiftPolyStrategy):
    """Optimal shift polynomial selection strategy with provable guarantees."""

    def __init__(self, *args, use_intermediate_sizes=True, **kwargs):
        super().__init__(*args, **kwargs)
        if not isinstance(use_intermediate_sizes, bool):
            raise TypeError("use_intermediate_sizes should be True or False")
        self.use_intermediate_sizes = use_intermediate_sizes

    def _get_monomial_set(self, modulus, ring, ideal_inf, bounds: BoundSet):
        bound_vals = [bounds.get_abs_bound(x) for x in ring.gens()]
        lg_bounds = list(map(math.log2, bound_vals))
        if modulus is None:
            # Working without a modulus
            inf_lms = []
            max_weight = float("inf")
        else:
            inf_lms = [g.lm() for g in ideal_inf.groebner_basis()]
            max_weight = math.log2(bounds.get_upper_bound(modulus))

        for exps, weight in weighted_combinations(lg_bounds):
            if weight >= max_weight:
                break

            monom = 1
            for xi, ei in zip(ring.gens(), exps):
                monom *= xi**ei

            for g_lm in inf_lms:
                if monom % g_lm == 0:
                    break
            else:
                yield monom

    def _get_shift_poly_with_lm(self, monom, groebner_basis, ideal):
        T = [g for g in groebner_basis if monom % g.lm() == 0]
        if len(T) == 0:
            return None

        g = min(T, key=lambda g: abs(g.lc()))
        h = g * monom // g.lm()
        h_prime = h.lt() + ideal.ideal.reduce(h - h.lt())
        return Relation(h_prime, ideal.modulus)

    def _is_intermediate_output_size(self, sz):
        # Return true if sz >= 10 and sz has binary form
        # 0b1000... or
        # 0b1100...
        if sz < 2:
            return False
        if is_power_of_two(sz):
            return True
        if sz % 3 == 0 and is_power_of_two(sz // 3):
            return True
        return False

    def _get_shift_polys_for_ideal(
        self, ideal, ideal_inf, bounds
    ) -> Iterator[RelationSet]:
        self.logger.debug("Calculating monomials within modulus bound.")
        monomials = self._get_monomial_set(
            ideal.modulus, ideal.ring(), ideal_inf, bounds
        )
        self.logger.debug("Computing Groebner basis.")
        groebner_basis = ideal.groebner_basis()
        shift_polys = []
        for monom in monomials:
            shift_poly = self._get_shift_poly_with_lm(monom, groebner_basis, ideal)
            if shift_poly:
                shift_polys += [shift_poly]
            if self.use_intermediate_sizes and self._is_intermediate_output_size(
                len(shift_polys)
            ):
                yield RelationSet(shift_polys)
        if not self.use_intermediate_sizes or not self._is_intermediate_output_size(
            len(shift_polys)
        ):
            if len(shift_polys) != 0:
                # We exited the loop and did not yield on the last loop iteration
                yield RelationSet(shift_polys)
        self.logger.info(
            "Computed %u shift polynomials with monomials smaller than the modulus",
            len(shift_polys),
        )

    def run(
        self, input_relations: RelationSet, bounds: BoundSet
    ) -> Iterator[RelationSet]:
        ideal_generator = RelationIdealGenerator()
        for ideal, ideal_inf in ideal_generator.run(input_relations, bounds):
            self.logger.debug("Generating optimal shift polynomials for all monomials")
            for shift_polys in self._get_shift_polys_for_ideal(
                ideal, ideal_inf, bounds
            ):
                yield shift_polys
