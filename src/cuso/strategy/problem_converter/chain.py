"""This file contains tools for chaining together multiple problem converters."""

from typing import Tuple, List

from cuso.data import (
    MultivariateCoppersmithProblem,
    SolutionConverter,
    RelationConverter,
)
from cuso.data.types import Polynomial

from .problem_converter import MultivariateProblemConverter
from .monom_ordering import BoundedMonomialOrderConverter
from .recenter import RecenterConverter
from .unraveled_linearization import UnraveledLinearization


class ChainSolnConverter(SolutionConverter):
    """Chain together multiple SolutionConverters

    Args:
        SolutionConverter (List[SolutionConverter]): List of
            solution converters to chain together.
    """

    def __init__(self, soln_converters: List[SolutionConverter]):
        self.soln_converters: List[SolutionConverter] = soln_converters

    def convert_solution_to_old(self, solution):
        for soln_converter in self.soln_converters[::-1]:
            solution = soln_converter.convert_solution_to_old(solution)
        return solution

    def convert_solution_to_new(self, solution):
        for soln_converter in self.soln_converters:
            solution = soln_converter.convert_solution_to_new(solution)
        return solution


class ChainRelConverter(RelationConverter):
    """Chain together multiple RelationConverters

    Args:
        RelationConverter (List[RelationConverter]): List of
            relation converters to chain together.
    """

    def __init__(self, rel_converters: List[RelationConverter]):
        self.rel_converters: List[RelationConverter] = rel_converters

    def convert_polynomial_to_old(self, poly):
        for rel_converter in self.rel_converters[::-1]:
            poly = rel_converter.convert_polynomial_to_old(poly)
        return poly

    def convert_polynomial_to_new(self, poly):
        for rel_converter in self.rel_converters:
            poly = rel_converter.convert_polynomial_to_new(poly)
        return poly


class ChainConverter(MultivariateProblemConverter):
    """Chain together multiple MultivariateProblemConverters."""

    def __init__(
        self, do_recentering: bool = True, unrav_lin_relations: List[Polynomial] = None
    ):
        """Initialize the ChainConverter.

        Args:
            do_recentering (bool, optional): If True, recenter variables so the center of the
                bounded range is at 0.
            unrav_lin_relations (List[Polynomial], optional): If nonempty, also perform unraveled
                linearization using the specified polynomials
        """
        super().__init__()

        self.unrav_lin_relations = unrav_lin_relations
        self.do_recentering = do_recentering

    def run(
        self, problem: MultivariateCoppersmithProblem
    ) -> Tuple[MultivariateCoppersmithProblem, RelationConverter, SolutionConverter]:

        # Converters to use
        converters = []
        if self.do_recentering:
            converters += [RecenterConverter()]
        converters += [BoundedMonomialOrderConverter()]

        soln_converters = []
        rel_converters = []
        cur_prob = problem
        for conv in converters:
            new_prob, new_rel_converter, new_soln_converter = conv.run(cur_prob)
            cur_prob = new_prob
            soln_converters += [new_soln_converter]
            rel_converters += [new_rel_converter]

        rel_converter = ChainRelConverter(rel_converters)

        if self.unrav_lin_relations and len(self.unrav_lin_relations) != 0:
            new_ul_rels = [
                rel_converter.convert_polynomial_to_new(f_ul)
                for f_ul in self.unrav_lin_relations
            ]
            ul_conv = UnraveledLinearization(new_ul_rels)
            converters += [ul_conv]

            new_prob, new_rel_converter, new_soln_converter = ul_conv.run(cur_prob)
            cur_prob = new_prob
            soln_converters += [new_soln_converter]
            rel_converters += [new_rel_converter]
            rel_converter = ChainRelConverter(rel_converters)

        soln_converter = ChainSolnConverter(soln_converters)

        return cur_prob, rel_converter, soln_converter
