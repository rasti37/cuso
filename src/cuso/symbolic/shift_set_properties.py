from collections import namedtuple
import itertools
from functools import cache

from sage.all import Polyhedron

from .utils import lm_no_N

LatticeProperty = namedtuple("LatticeProperty", "dimension, X_s, LC")

class AugPolyhedron:
    def __init__(self, points, pointvals):
        self.points = points
        self.pointvals = pointvals

    def ambient_dim(self):
        return len(self.points[0])
    
    def dim(self):
        return Polyhedron(vertices=self.points).dim()

    @staticmethod
    def from_polytope(polytope, lc_vals=None):
        # Fill in interior points
        points = polytope.integral_points()

        origin = (0,) * len(points[0])
        points_with_origin = Polyhedron(vertices = [origin] + list(points)).integral_points()
        if lc_vals is None:
            lc_vals = {(0,)*len(points[0]): 0}
        aug_points = {}
        for point in points_with_origin:
            minval = None
            for k, v in lc_vals.items():
                # Is this monomial less than ours?
                if all(ai <= bi for ai, bi in zip(k, point)):
                    if minval is None:
                        minval = v
                    minval = min(v, minval)
            aug_points[tuple(point)] = minval
        points = [tuple(point) for point in points]
        return AugPolyhedron(points, aug_points)
    
    def value(self, point):
        if point in self.pointvals:
            return self.pointvals[point]
        
        v = None
        for x in self.pointvals.keys():
            if all(ai <= bi for ai, bi in zip(x, point)):
                if v is None:
                    v = self.pointvals[x]
                else:
                    v = min(v, self.pointvals[x])
        self.pointvals[point] = v
        return v

    def __add__(self, other):
        assert isinstance(other, AugPolyhedron)

        P1 = self.pointvals
        P2 = other.pointvals

        new_pointvals = {}
        new_pointvals_arr = []
        for p1 in P1:
            for p2 in P2:
                v = P1[p1] + P2[p2]
                new_p = tuple(ai + bi for ai, bi in zip(p1, p2))
                new_pointvals_arr += [(new_p, v)]
        new_pointvals_arr = sorted(new_pointvals_arr)
        prev_k = None
        for k, v in new_pointvals_arr:
            if k != prev_k:
                new_pointvals[k] = v
                prev_k = k


        new_points = []
        for p1 in self.points:
            for p2 in other.points:
                new_p = tuple(ai + bi for ai, bi in zip(p1, p2))
                new_points += [new_p]
        new_points = list(set(new_points))
        return AugPolyhedron(new_points, new_pointvals)
    
    def t_shift(self, max_displacement):
        # Return an augmented polynomial that is this one, shifted by all
        # amounts up to displacement. For example, if displacement is (2,),
        # then return this polytope union (polytope shifted by 1) union
        # (polytope shifted by 2)
        disp_iter = itertools.product(*[range(1+max_t) for max_t in max_displacement])
        new_pts = []
        new_vals = self.pointvals.copy()
        for disp in disp_iter:
            for p_orig in self.pointvals:
                p_shift = tuple(ai + bi for ai, bi in zip(p_orig, disp))
                if p_shift not in new_vals:
                    new_vals[p_shift] = self.pointvals[p_orig]
                else:
                    new_vals[p_shift] = min(new_vals[p_shift], self.pointvals[p_orig])
                if p_orig in self.points and p_shift not in new_pts:
                    new_pts += [p_shift]
        return AugPolyhedron(new_pts, new_vals)
    
    def unravel(self, ul):
        if len(ul) == 0:
            return AugPolyhedron(self.points, self.pointvals)
        
        points = self.points
        while True:
            did_unravel = False
            new_points = []
            for point in points:
                for f in ul:
                    if all(ai >= bi for ai, bi in zip(point, f)):
                        # To be unraveled
                        break
                else:
                    # Was not unravelled
                    new_points += [point]

            points = new_points
            if not did_unravel:
                break
        return AugPolyhedron(points, self.pointvals)
    
    def __str__(self):
        return str(self.points)

class ShiftProperties:
    def __init__(self, S_1, G_inf):
        self.S_1 = S_1
        self.G_inf = G_inf

        self.G_inf_LM = [lm_no_N(g) for g in G_inf]
        self.ul = [g.degrees()[1:] for g in self.G_inf_LM]

    def dim(self):
        S_1 = self._S_k(1)
        return S_1.dim()

    @cache
    def _S_k(self, k):
        if k > 1:
            S_prev = self._S_k(k - 1)
            S_1 = self._S_k(1)
            S_k = S_prev + S_1
            return S_k

        vertices = []
        for f in self.S_1:
            for m in f.monomials():
                # This monomial is in the polytope,
                # ignore the degree of N.
                vertices += [m.degrees()[1:]]
        P = Polyhedron(vertices=vertices)

        # Augmented polytope which includes the
        # contribution of N
        lc_vals = {}
        for f in self.S_1:
            lt = f.lm()
            lgN = lt.degrees()[0]
            monom = lt.degrees()[1:]
            lc_vals[monom] = lgN

        S_1 = AugPolyhedron.from_polytope(P, lc_vals=lc_vals)
        return S_1
    
    @cache
    def _S_kt(self, k, ts):
        if sum(ts) == 0:
            S_kt = self._S_k(k)
            return S_kt
        
        first_nz = 0
        while ts[first_nz] == 0:
            first_nz += 1

        new_ts = list(ts)[:]
        t_disp = [0] * len(ts)
        new_ts[first_nz] -= 1
        t_disp[first_nz] = 1
        S_kt_prev = self._S_kt(k, tuple(new_ts))
        S_kt = S_kt_prev.t_shift(t_disp)
        return S_kt

    @cache
    def _S_ktul(self, k, ts):
        S_kt = self._S_kt(k, ts)
        S_ktul = S_kt.unravel(self.ul)
        return S_ktul

    def _aug_poly_to_property(self, S) -> LatticeProperty:
        polytope_points = {
            pt: S.value(pt) for pt in S.points
        }

        lattice_dim = len(polytope_points)

        per_coordinate_sum = [0] * S.ambient_dim()
        for point in polytope_points:
            for i, x_i in enumerate(point):
                per_coordinate_sum[i] += x_i

        lc_term = 0
        for pt in polytope_points:
            lc_term += polytope_points[pt]
        
        return LatticeProperty(lattice_dim, per_coordinate_sum, lc_term)

    def s_ktul(self, k, ts) -> LatticeProperty:
        S_ktul = self._S_ktul(k, tuple(ts))
        return self._aug_poly_to_property(S_ktul)

    def s_kt(self, k, ts) -> LatticeProperty:
        S_kt = self._S_kt(k, tuple(ts))
        return self._aug_poly_to_property(S_kt)

    def s_k(self, k) -> LatticeProperty:
        S_k = self._S_k(k)
        return self._aug_poly_to_property(S_k)

