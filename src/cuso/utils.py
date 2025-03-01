"""This file contains common mathematical utilities."""

from typing import Iterator, Tuple, List
from heapq import heappush, heappop

from cuso.data.types import Polynomial


def weighted_combinations(weights: Tuple[float]) -> Iterator[Tuple[List[int], float]]:
    """Return integer combinations of increasing weight.

    Given a weight vector (w_1, ..., w_k), this method returns
    a tuple of nonnegative integers (e_1, ..., e_k) and a weight
    W = e_1*w_1 + ... + e_k*w_k such W increases monotonically.

    Args:
        weights (Tuple[float]): vector of weights

    Yields:
        Tuple[List[int], float]: tuple of integers and combined weight
    """
    n = len(weights)

    heap = []
    heappush(heap, (0, (0,) * n))

    while True:
        score, exps = heappop(heap)
        for i in range(n):
            new_exps = list(exps)
            new_exps[i] += 1
            new_exps = tuple(new_exps)
            new_score = sum(ai * bi for ai, bi in zip(weights, new_exps))
            if (new_score, new_exps) not in heap:
                heappush(heap, (new_score, new_exps))
        yield exps, score


def is_suitable(polys: List[Polynomial]) -> bool:
    """Checks whether a set of polynomials is suitable.

    Given S, first we compute M = support(S), then we check if
    |M| = |S| and M = LM(S).

    Args:
        polys (List[Polynomial]): collection of polynomials to check

    Returns:
        bool: True if S is (M, <)-suitable
    """
    all_monoms = []
    for poly in polys:
        all_monoms += poly.monomials()
    all_monoms = set(all_monoms)

    lms = [f.lm() for f in polys]
    return len(lms) == len(all_monoms) and set(lms) == all_monoms
