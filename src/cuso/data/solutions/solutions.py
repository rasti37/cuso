"""Define the Solution and SolutionSet types.

These are just wrappers around dictionaries and lists.
"""

__all__ = [
    "PartialSolution",
    "PartialSolutionSet",
    "Solution",
    "SolutionSet",
]


class PartialSolution(dict):
    """Representation of a partial solution to a system of relations"""


class PartialSolutionSet(list):
    """Represent a collection of partial solutions"""


class Solution(PartialSolution):
    """Representation of a solution to a system of relations"""


class SolutionSet(PartialSolutionSet):
    """Represent a collection of Solutions"""
