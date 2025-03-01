"""This file implements the graph-based shift polynomial selection strategy."""

from typing import Iterator, Optional
import math

import igraph as ig

from cuso.data import BoundSet, RelationSet, Relation
from cuso.utils import is_suitable

from .shift_poly_selection import ShiftPolyStrategy


class GraphShiftPolys(ShiftPolyStrategy):
    """Heuristic shift polynomial selection strategy based on graph optimization."""

    def __init__(self, sub_shift_polys: ShiftPolyStrategy):
        super().__init__()
        self.sub_shift_polys: ShiftPolyStrategy = sub_shift_polys

    def _approx_shvec_bound(self, shift_polys, bounds: BoundSet):
        logdet = sum(math.log2(bounds.get_abs_bound(f.lt())) for f in shift_polys)
        return logdet / len(shift_polys)

    def _maximum_closure(self, G):
        # Use Picard's algorithm to find maximum closure
        G_Picard = G.copy()
        weights = G.vs["weight"]

        # Add source and sink
        G_Picard.add_vertices(2)
        n = len(weights)
        source = n
        sink = n + 1

        capacities = []
        for _ in G_Picard.get_edgelist():
            # Existing edges get infinite capacity
            capacities += [float("inf")]
        for i, w_i in enumerate(weights):
            if w_i > 0:
                G_Picard.add_edges([(source, i)])
                capacities += [w_i]
            else:
                G_Picard.add_edges([(i, sink)])
                capacities += [-w_i]

        G_Picard.es["capacity"] = capacities

        cut = G_Picard.mincut(source, sink, capacity=G_Picard.es["capacity"])
        part = cut.partition
        # partition that includes source
        part = [p for p in part if source in p][0]
        return [node for node in part if node != source]

    def _refine_once(self, polys, bounds: BoundSet) -> Optional[RelationSet]:
        monoms = sorted([f.lm() for f in polys])

        # Ensure polys are in the same order as monomials
        polys = sorted(polys)

        # Build weighted graph
        edges = []
        for f in polys:
            f_monoms = f.monomials()
            f_lm = max(f_monoms)
            for m in f_monoms:
                if m == f_lm:
                    continue
                i = monoms.index(f_lm)
                j = monoms.index(m)
                edges += [(i, j)]
        G = ig.Graph(len(monoms), edges, directed=True)

        # Set weight of vertex. Weight is based on leading term of f
        unnormalized_weights = [math.log2(bounds.get_abs_bound(f.lt())) for f in polys]
        avg = sum(unnormalized_weights) / len(unnormalized_weights)
        weights = [-w + avg for w in unnormalized_weights]
        G.vs["weight"] = weights

        subset = self._maximum_closure(G)
        closure_total_weight = sum([weights[i] for i in subset])
        if abs(closure_total_weight) < 1e-6:
            return None
        else:
            return [polys[i] for i in subset]

    def _refine_shift_polys(self, shift_polys, bounds) -> Optional[RelationSet]:
        self.logger.debug("Running graph optimizer on %d polynomials.", len(shift_polys))

        modulus = shift_polys[0].modulus
        polys = []
        for rel in shift_polys:
            if modulus != rel.modulus:
                raise ValueError("All shift polynomials must share the same modulus")
            polys += [rel.polynomial]
        if not is_suitable(polys):
            raise ValueError("Shift polynomials must be (M,<)-suitable")

        num_iters = 0
        while True:
            num_iters += 1
            result = self._refine_once(polys, bounds)
            if result is None:
                break
            else:
                polys = result
        self.logger.debug("Converged after %d iterations.", num_iters)
        self.logger.debug("Found subset of %d polynomials.", len(polys))

        # Look at expected length of shortest vector
        log_shvec_len = self._approx_shvec_bound(polys, bounds)
        log_mod = math.log2(bounds.get_lower_bound(modulus))
        if log_shvec_len < log_mod:
            # If determinant bound is beneficial, return these relations
            relations = [Relation(f, modulus) for f in polys]
            relations = RelationSet(relations)
            self.logger.info(
                "Graph optimization found a subset of %u shift relations (out of %u)",
                len(relations),
                len(shift_polys),
            )
            return relations
        else:
            self.logger.debug(
                "Expected shortest vector is not smaller than modulus. Skipping."
            )
            return None

    def run(
        self, input_relations: RelationSet, bounds: BoundSet
    ) -> Iterator[RelationSet]:
        for shift_polys in self.sub_shift_polys.run(input_relations, bounds):
            refined = self._refine_shift_polys(shift_polys, bounds)
            if refined is not None:
                yield refined
