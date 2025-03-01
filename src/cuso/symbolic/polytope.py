import functools
import logging
import operator

from sage.all import Polyhedron, PolynomialRing

logger = logging.getLogger("cuso.symbolic.Polytope")

def get_monomials_from_vertices(R: PolynomialRing, m_vertices):
    """Return the convex hull of a set of monomials.

    Args:
        R (PolynomialRing): Parent ring of the monomials
        m_vertices (List[Polynomial]): The set M_vert

    Returns:
        List[Polynomial]: List of monomials in the convex hull
    """
    logger.debug("M_vert has %d monomials", len(m_vertices))
    if isinstance(m_vertices[0], tuple):
        vertices = m_vertices
    else:
        vertices = [
            R(v).degrees() for v in m_vertices
        ]
    P = Polyhedron(
        vertices=vertices
    )
    xs = R.gens()
    M = []
    for pt in P.integral_points():
        m = functools.reduce(operator.mul, [xi**ei for (xi, ei) in zip(xs, pt)], 1)
        M += [m]
    logger.debug("M_1 has %d monomials", len(M))
    return M
