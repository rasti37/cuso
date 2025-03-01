"""This module implements lattice building strategies.

Given a set of shift polynomials, build a lattice such that
reducing the lattice helps recover the shared root of the
shift polynomials. Usually we use the Howgrave-Graham dual
construction, but in some cases we also use Coppersmith's
primal construction.
"""

from .lattice_builder import LatticeBuilder

from .howgrave_graham_dual import DualLatticeBuilder
from .primal import PrimalLatticeBuilder
