""" This module provides a helpful wrapper function for using cuso.
"""

import logging
import operator
from collections import namedtuple
from typing import Dict, List, Optional, Tuple, Union

from sage.all import Polynomial, Expression, Integer, ZZ, var, PolynomialRing
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.abc import IntegerModRing

from cuso.data import (
    BoundSet,
    Relation,
    RelationSet,
    Solution,
    SolutionSet,
    PartialSolution,
    PartialSolutionSet,
    MultivariateCoppersmithProblem,
)
from cuso.solver import AutomatedSolver, AutomatedPartialSolver

logger = logging.getLogger("cuso.Wrapper")

# Define types that can be converted to cuso's internal types
BoundsLike = Union[BoundSet, Dict, List, int, Integer]
ModulusLike = Union[int, Integer, str, Expression]
ModuliLike = Union[ModulusLike, List[ModulusLike]]
PolynomialLike = Union[Polynomial, MPolynomial, Expression]
RelationLike = Union[Relation, PolynomialLike]
RelationSetLike = Union[RelationSet, List[RelationLike], RelationLike]
SolutionSetLike = Union[SolutionSet, List[Dict], Dict]

ModulusInformation = namedtuple(
    "ModulusInformation",
    ["is_symbolic", "value", "arg_value", "multiple", "lower_bound", "upper_bound"],
)


def do_varsub(item, varsub: dict):
    """Substitute variables.

    Args:
        item (_type_): object with old variables
        varsub (dict): mapping of old variables to new variables

    Raises:
        TypeError: item not of recognized type.

    Returns:
        _type_: item with old variables replaced
    """
    if item is None:
        return None
    if isinstance(item, dict):
        new_dict = {}
        for xi in item:
            if xi in varsub:
                new_xi = varsub[xi]
            else:
                new_xi = xi
            new_dict[new_xi] = item[xi]
        return new_dict
    if isinstance(item, list):
        return [do_varsub(item_i, varsub) for item_i in item]
    if isinstance(item, Solution):
        return Solution(do_varsub(dict(item), varsub))
    if isinstance(item, SolutionSet):
        return SolutionSet(do_varsub(list(item), varsub))
    if isinstance(item, PartialSolution):
        return PartialSolution(do_varsub(dict(item), varsub))
    if isinstance(item, PartialSolutionSet):
        return PartialSolutionSet(do_varsub(list(item), varsub))
    if isinstance(item, (Polynomial, MPolynomial, Expression)):
        return item.subs(varsub)
    raise TypeError


def parse_modulus(
    modulus: ModulusLike = None,
    modulus_multiple: Optional[int] = None,
    modulus_lower_bound: Optional[int] = None,
    modulus_upper_bound: Optional[int] = None,
) -> Optional[ModulusInformation]:
    """Return a consistent description of the modulus, if there is one.

    Args:
        modulus (int, str, Expression, optional): symbolic or integer modulus value. Defaults to None.
        modulus_multiple (int, optional): multiple of modulus. Defaults to None.
        modulus_lower_bound (int, optional): lower bound of modulus. Defaults to None.
        modulus_upper_bound (int, optional): upper bound of modulus. Defaults to None.

    Raises:
        ValueError: input is not consistent
        TypeError: input is invalid

    Returns:
        ModulusInformation: tuple of (is_symbolic, value, multiple, lbound, ubound)
    """
    if not (modulus or modulus_multiple or modulus_lower_bound or modulus_upper_bound):
        return None

    # Return (is_symbolic, value, multiple, lbound, ubound)
    is_symbolic = None
    if isinstance(modulus, (int, Integer)):
        is_symbolic = False
        value = int(modulus)
        if modulus_multiple is not None:
            raise ValueError("Do not set both integer modulus and modulus multiple")
        if modulus_lower_bound is not None:
            raise ValueError("Do not set both integer modulus and modulus lower bound")
        if modulus_upper_bound is not None:
            raise ValueError("Do not set both integer modulus and modulus upper bound")
        return ModulusInformation(is_symbolic, value, value, value, value, value)
    if modulus is None or isinstance(modulus, (str, Expression)):
        is_symbolic = True
        arg_value = modulus
        if modulus is None:
            value = None
        elif isinstance(modulus, str):
            value = var(modulus)
        else:
            value = modulus
        if modulus_lower_bound is None:
            raise ValueError("Must specify modulus_lower_bound for symbolic modulus")
        lbound = int(modulus_lower_bound)
        multiple = int(modulus_multiple) if modulus_multiple else None
        if modulus_upper_bound:
            ubound = int(modulus_upper_bound)
        else:
            if modulus_multiple:
                ubound = multiple
            else:
                ubound = None
        return ModulusInformation(
            is_symbolic, value, arg_value, multiple, lbound, ubound
        )
    raise TypeError("Unrecognized modulus configuration")


def parse_moduli(
    modulus: Optional[ModuliLike] = None,
    modulus_multiple: Optional[Union[int, List[int]]] = None,
    modulus_lower_bound: Optional[Union[int, List[int]]] = None,
    modulus_upper_bound: Optional[Union[int, List[int]]] = None,
) -> Optional[Union[ModulusInformation, List[ModulusInformation]]]:
    """The moduli may be specified as either a single, shared modulus, or as a list.

    Detect which is the case, and return either the shared modulus or a list of the
    individual moduli.

    Args:
        modulus (Optional[ModuliLike], optional): [List of] modulus values.
        modulus_multiple (Optional[Union[int, List[int]]], optional): [List of] modulus multiples.
        modulus_lower_bound (Optional[Union[int, List[int]]], optional): [List of] modulus lower bounds.
        modulus_upper_bound (Optional[Union[int, List[int]]], optional): [List of] modulus upper bounds.

    Returns:
        Optional[Union[ModulusInformation, List[ModulusInformation]]]: [List of] individual moduli
    """
    args = [modulus, modulus_multiple, modulus_lower_bound, modulus_upper_bound]
    if any(isinstance(arg, list) for arg in args):
        num_moduli = None
        for arg in args:
            if isinstance(arg, list):
                num_moduli = len(arg)
            elif arg is not None:
                err = (
                    "If one of modulus, modulus_multiple, modulus_lower_bound, modulus_upper_bound"
                    " is a list, all must be either a list or None."
                )
                raise ValueError(err)
        # Check that all have the same length
        for i, arg_i in enumerate(args):
            if arg_i is None:
                args[i] = [None] * num_moduli
            else:
                if len(arg_i) != num_moduli:
                    err = "If using lists, all modulus.* arguments must have the same length."
                    raise ValueError(err)
        moduli, moduli_multiple, moduli_lower_bound, moduli_upper_bound = args
        return [
            parse_modulus(mod, mod_mul, mod_lb, mod_ub)
            for mod, mod_mul, mod_lb, mod_ub in zip(
                moduli, moduli_multiple, moduli_lower_bound, moduli_upper_bound
            )
        ]
    return parse_modulus(
        modulus, modulus_multiple, modulus_lower_bound, modulus_upper_bound
    )


def parse_relations(
    relations: RelationSetLike,
    modulus: ModuliLike = None,
    modulus_multiple: Optional[int] = None,
    modulus_lower_bound: Optional[int] = None,
    modulus_upper_bound: Optional[int] = None,
) -> Tuple[RelationSet, List[ModulusInformation], PolynomialRing, Dict]:
    """Take user-provided relations as input and convert to a consistent type.

    Args:
        relations (RelationSetLike): Input relation(s). Can be a RelationSet, a single Relation,
                                     a Polynomial, or an Expression
        modulus (ModuliLike, optional): Moduli for the relations, if there is one.
        modulus_multiple (int, optional): Integer multiple of the moduli, if it is known.
        modulus_lower_bound (int, optional): Lower bound of the moduli, if it is known.
        modulus_upper_bound (int, optional): Upper bound of the moduli, if it is known.

    Returns:
        Tuple[RelationSet, List[ModulusInformation], PolynomialRing, Dict]:
            the cuso RelationSet, as well as other useful information.
    """
    inferred_moduli = parse_moduli(
        modulus=modulus,
        modulus_multiple=modulus_multiple,
        modulus_lower_bound=modulus_lower_bound,
        modulus_upper_bound=modulus_upper_bound,
    )

    expect_modulus = False
    as_polys = False
    varsub = {}
    if isinstance(relations, RelationSet):
        expect_modulus = False
    else:
        if not isinstance(relations, list):
            if isinstance(inferred_moduli, list) and len(inferred_moduli) != 1:
                raise ValueError(
                    "Cannot specify moduli as list if there is only one relation."
                )
            relations = [relations]
        if len(relations) == 0:
            raise ValueError("Must specify at least one relation.")
        if isinstance(inferred_moduli, list):
            if len(inferred_moduli) != len(relations):
                raise ValueError("Must have as many moduli as relations.")
        else:
            inferred_moduli = [inferred_moduli for _ in range(len(relations))]
        rel_types = set(type(rel) for rel in relations)
        if len(rel_types) != 1:
            raise TypeError("All relations must have the same type!")
        rel_type = list(rel_types)[0]
        if issubclass(rel_type, Relation):
            expect_modulus = False
            relations = RelationSet(relations)
        elif issubclass(rel_type, (Polynomial, MPolynomial)):
            # Make sure all polynomials are defined over ZZ[x].
            # If any are defined over ZZ_N[x], infer the modulus is N
            polys = []
            moduli = []
            parent_ring = None
            for poly in relations:
                poly_ring = poly.base_ring()
                if poly_ring == ZZ:
                    polys += [poly]
                    moduli += [None]
                elif isinstance(poly_ring, IntegerModRing):
                    polys += [poly.change_ring(ZZ)]
                    moduli += [int(poly_ring.characteristic())]
                else:
                    raise TypeError(
                        "Input relation is a polynomial not in ZZ or Integers(N)"
                    )
                if parent_ring:
                    if polys[-1].parent() != parent_ring:
                        raise TypeError(
                            "Input polynomials not defined over same polynomial ring."
                        )
                else:
                    parent_ring = polys[-1].parent()
            if any(moduli):
                # Some of the polys are defined over Z_N
                expect_modulus = False
            else:
                # All polys are defined over ZZ
                expect_modulus = True
            as_polys = True
        elif issubclass(rel_type, Expression):
            # Convert expressions to elements of polynomial ring ZZ[x]
            unk_vars = []
            for rel in relations:
                for xi in rel.variables():
                    if xi not in unk_vars:
                        unk_vars += [xi]
            unk_var_names = [str(xi) for xi in unk_vars]
            parent_ring = PolynomialRing(ZZ, unk_var_names)
            for orig_xi, new_xi in zip(unk_vars, parent_ring.gens()):
                varsub[orig_xi] = new_xi
            polys = []
            moduli = []

            for rel in relations:
                if rel.is_relational():
                    if rel.operator() != operator.eq:
                        raise ValueError(
                            "Only relational expressions of the form A == B are supported."
                        )
                    rel = rel.left_hand_side() - rel.right_hand_side()
                poly = parent_ring(rel)
                polys += [poly]
                moduli += [None]

            expect_modulus = True
            as_polys = True
        else:
            raise TypeError(f"Unrecognized relation type {rel_type}")

    # if expect_modulus and inferred_modulus is None:
    #    raise ValueError("Modulus information is expected.")
    no_inferred_moduli = (inferred_moduli is None) or (
        isinstance(inferred_moduli, list)
        and all(inferred_modulus is None for inferred_modulus in inferred_moduli)
    )
    if not expect_modulus and not no_inferred_moduli:
        raise ValueError(
            "Modulus should not be specified if it can be inferred from the relations"
        )

    if expect_modulus:
        # polys is a list of polynomials
        assert as_polys
        assert isinstance(inferred_moduli, list)
        new_polys = []
        moduli = []
        for i, (poly, inferred_modulus) in enumerate(zip(polys, inferred_moduli)):
            if inferred_modulus is not None:
                is_symbolic, value, arg_value, multiple, lbound, ubound = (
                    inferred_modulus
                )
                if not is_symbolic:
                    modulus = value
                else:
                    if value is None:
                        # Need to pick modulus name
                        ring_vars = parent_ring.variable_names()
                        if "p" not in ring_vars and len(inferred_moduli) == 1:
                            modname = "p"
                        elif "q" not in ring_vars and len(inferred_moduli) == 1:
                            modname = "q"
                        else:
                            i = 0
                            while True:
                                modname = f"p_{i}"
                                if modname not in ring_vars:
                                    break
                        new_modulus = var(modname)
                        modulus = new_modulus
                    else:
                        modulus = value
                    if arg_value is not None:
                        varsub[arg_value] = modulus
                    inferred_modulus = ModulusInformation(
                        is_symbolic, modulus, arg_value, multiple, lbound, ubound
                    )
                    inferred_moduli[i] = inferred_modulus
            else:
                is_symbolic = False
                modulus = None
            new_polys += [poly]
            moduli += [modulus]
            if is_symbolic and multiple:
                # Add polynomial for the multiple of p
                new_polys += [parent_ring(multiple)]
                moduli += [modulus]
        polys = new_polys
    if as_polys:
        relations = []
        for poly, mod in zip(polys, moduli):
            relations += [Relation(poly, mod)]
        relations = RelationSet(relations)
    parent_ring = relations[0].polynomial.parent()

    return relations, inferred_moduli, parent_ring, varsub


def parse_bounds(
    bounds: BoundsLike,
    ring: PolynomialRing,
    inferred_moduli: List[ModulusInformation],
    varsub: Dict,
) -> BoundSet:
    """Return a consistent representation of the Coppersmith bounds.

    Args:
        bounds (BoundsLike): Description of the bounds.
        ring (PolynomialRing): polynomial ring.
        inferred_modulus (ModulusInformation): common modulus information.
        varsub (Dict): Variable substitution

    Returns:
        BoundSet: the cuso BoundSet for this problem
    """
    # Start by turning bounds into a dict
    do_replace = False
    if isinstance(bounds, BoundSet):
        do_replace = True
    elif isinstance(bounds, dict):
        do_replace = True
    elif isinstance(bounds, list):
        unknowns = ring.gens()
        if len(bounds) != len(unknowns):
            raise ValueError(
                f"If specifying a list of bounds, it must correspond to unknowns {ring.gens()}"
            )
        if any(int(bound_i) <= 0 for bound_i in bounds):
            raise ValueError("Bounds must be positive")
        bounds = {xi: (-bound_i, bound_i) for xi, bound_i in zip(ring.gens(), bounds)}
        do_replace = False
    elif isinstance(bounds, (int, Integer)):
        bounds = int(bounds)
        if bounds <= 0:
            raise ValueError("Bounds must be positive")
        bounds = {xi: (-bounds, bounds) for xi in ring.gens()}
        do_replace = False
    else:
        raise TypeError("Unexpected bounds type")

    if do_replace:
        bounds = do_varsub(bounds, varsub)

    if inferred_moduli is not None:
        for inferred_modulus in inferred_moduli:
            if inferred_modulus is None:
                continue
            is_symbolic, value, _, __, lbound, ubound = inferred_modulus
            if is_symbolic:
                bounds[value] = (lbound, ubound)

    # bounds is dict. If any values are integer, convert to tuple
    bounds = {
        k: (-v, v) if isinstance(v, (int, Integer)) else v for k, v in bounds.items()
    }
    bounds = BoundSet(bounds)
    return bounds


def parse_solution(solutions: SolutionSetLike, varsub: Dict) -> SolutionSet:
    """Return a consistent representation of the expected solution.

    Args:
        solutions (SolutionSetLike): Expected solution.
        varsub (Dict): Variable substitution information.

    Raises:
        TypeError: solutions are the wrong type.

    Returns:
        SolutionSet: The cuso SolutionSet for this problem
    """
    if solutions is None:
        return None
    if isinstance(solutions, dict):
        solutions = [solutions]
    if not isinstance(solutions, list):
        raise TypeError("Expected solution should be dict or list of dicts.")
    # If variables in solution are represented by strings, convert to variables
    str_lookup = {str(x_i): x_i for x_i in varsub.keys()}
    solutions = [
        {str_lookup.get(x_i, x_i): v_i for x_i, v_i in solution.items()}
        for solution in solutions
    ]
    return SolutionSet(do_varsub(solutions, varsub))


def find_small_roots(
    relations: RelationSetLike,
    bounds: BoundsLike,
    modulus: Optional[ModuliLike] = None,
    modulus_multiple: Optional[Union[int, List[int]]] = None,
    modulus_lower_bound: Optional[Union[int, List[int]]] = None,
    modulus_upper_bound: Optional[Union[int, List[int]]] = None,
    unraveled_linearization_relations: Optional[List[PolynomialLike]] = None,
    use_graph_optimization: Optional[bool] = None,
    use_intermediate_sizes: bool = True,
    allow_partial_solutions: bool = False,
    expected_solution: Optional[SolutionSetLike] = None,
) -> List[Dict]:
    """Find bounded roots of a system of polynomial equations.

    This method is intended to be permissive in how the problem is specified. It will
    attempt to convert the problem to cuso's internal types, inferring the RelationSet
    and BoundSet intended by the user.

    Args:
        relations (RelationSetLike): The input relations.
            This can be a cuso RelationSet, a single relation, or a list of relations.
            Each relation is either a cuso Relation, a Sage Polynomial, a Sage
            MPolynomial, or a Sage Expression.
        bounds (BoundsLike): The variable bounds.
            This can be a cuso BoundSet, a dictionary of bounds, a list of bounds on
            the absolute value, or a bound on the absolute value. If the bounds are
            given as a dictionary, the key is a variable and the value is either a
            tuple of (lower bound, upper bound) or a bound on the absolute value.
        modulus (ModulusLike, optional): The modulus for the input relations.
            If there is just one value here, then the modulus is shared by all
            relations. Otherwise it is a list of moduli, one for each relation.
            Each modulus can be an integer if the modulus is known, or a Sage variable
            or string if the modulus is unknown.
        modulus_multiple (Union[int, List[int]], optional): Integer multiple of the modulus.
        modulus_lower_bound (Union[int, List[int]], optional): Integer lower bound of the
            modulus, or moduli if given as a list.
        modulus_upper_bound (Union[int, List[int]], optional): Integer upper bound of the
            modulus, or moduli if given as a list.
        unraveled_linearization_relations (List[PolynomialLike], optional):
            List of unraveled linearization polynomials.
        use_graph_optimization (bool): If True, use the graph-based shift polynomial
            strategy. In some cases, such as when the input polynomials have small
            coefficients, the graph-based method may fail. If this is the case, it
            is sometimes worth disabling this strategy and falling back to the provable
            strategy, even though it is slower. Enabled by default, unless the input
            polynomials have small coefficients.
        use_intermediate_sizes (bool): If True, attempt to use fewer shift polynomials.
            Rather than build the entire set of optimal shift polynomials all at once,
            try the first 2, then the first 3, then 4, 6, 8, 12, 16, 24, 32, and so on.
            This can decrease the cost of calculating all shift polynomials if only a
            few are required. Defaults to True.
        allow_partial_solutions (bool): If True, return as soon as some variable value
            is known, even if some values are unknown. Defaults to False.
        expected_solution (SolutionSetLike, optional): The expected solution. When
            testing and debugging Coppersmith's method, it is sometimes useful to
            provide the intended root to ensure that intermediate results are computed
            correctly. This value can be a dictionary mapping variables to their value,
            a list of dictionaries, or a cuso SolutionSet.

    Returns:
        List[Dict]: A list of bounded solutions of the input polynomials, represented
            as dictionaries.
    """
    relations, inferred_moduli, ring, varsub = parse_relations(
        relations,
        modulus=modulus,
        modulus_multiple=modulus_multiple,
        modulus_lower_bound=modulus_lower_bound,
        modulus_upper_bound=modulus_upper_bound,
    )
    bounds = parse_bounds(bounds, ring, inferred_moduli, varsub)
    expected = parse_solution(expected_solution, varsub)

    ul_rels = do_varsub(unraveled_linearization_relations, varsub)

    problem = MultivariateCoppersmithProblem(relations, bounds)
    if allow_partial_solutions:
        solver_type = AutomatedPartialSolver
    else:
        solver_type = AutomatedSolver
    solver = solver_type(
        problem,
        unraveled_linearization_relations=ul_rels,
        use_intermediate_sizes=use_intermediate_sizes,
        use_graph_optimization=use_graph_optimization,
    )
    if expected:
        solver.set_expected(expected)
    solns = solver.solve()

    # Convert solns back using varsub
    varsub_inv = {v: k for k, v in varsub.items()}
    new_solns = do_varsub(solns, varsub_inv)

    # Return as list of dictionaries
    return [dict(soln) for soln in new_solns]


__all__ = ["find_small_roots"]
