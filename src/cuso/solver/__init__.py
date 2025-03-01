"""This module implements cuse "Solvers"

A Solver is used to solve problems whent there is a concise, easily checked solution.
For example, recovering the set of roots for a multivariate polynomial system.
"""

from .multivariate_solver import AutomatedSolver, AutomatedPartialSolver
