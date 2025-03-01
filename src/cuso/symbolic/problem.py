class SymbolicBounds:
    def __init__(self, a, b=None):
        if b is None:
            if isinstance(a, SymbolicBounds):
                a, b = a.ab
            elif isinstance(a, tuple) and len(a) == 2:
                a, b = a
            elif isinstance(a, list) and len(a) == 2:
                a, b = a
            else:
                a, b = 0, a
        # Bound is (a + b*delta) log_p
        self.ab = a, b
        self.linear_term = b
        self.constant_term = a

    def __str__(self):
        a, b = self.ab
        if b == 0:
            bd = ""
        elif b == 1:
            bd = "delta"
        else:
            bd = f"{b} * delta"
        if a == 0:
            return bd
        else:
            if b == 0:
                return f"{a}"
            return f"{a} + " + bd

class SymbolicCoppersmithProblem:
    def __init__(
        self,
        mod_rels,
        bounds,
        int_rels=None,
        modulus_multiple_name="N",
        modulus_name="p",
        divisor_bound=1,
        relation_modulus_degrees=None,
        mod_mult_modulus_degree=None,
    ):

        self.mod_rels = mod_rels
        if relation_modulus_degrees is None:
            self.rel_mod_deg = [1] * len(self.mod_rels)
        elif len(relation_modulus_degrees) != len(self.mod_rels):
            raise ValueError()
        else:
            self.rel_mod_deg = relation_modulus_degrees
        if not isinstance(bounds, list):
            bounds = [bounds]
        self.bounds = [SymbolicBounds(bound) for bound in bounds]
        self.int_rels = int_rels if int_rels is not None else []

        self.modulus_name = modulus_name
        if divisor_bound == 1:
            self.modulus_multiple_name = modulus_name
        else:
            self.modulus_multiple_name = modulus_multiple_name
        if mod_mult_modulus_degree is None:
            mod_mult_modulus_degree = 1
        self.mod_mult_modulus_degree = mod_mult_modulus_degree
        self.divisor_bound = divisor_bound
        self.ring = self.mod_rels[0].parent()

    def __str__(self):
        ring = self.mod_rels[0].parent()
        s = f"Multivariate Coppersmith problem in variables {ring.gens()}"
        s += " subject to constrained polynomials"
        for f, p_deg in zip(self.mod_rels, self.rel_mod_deg):
            if p_deg == 1:
                s += f"\n\t{f} == 0 (mod {self.modulus_name})"
            else:
                s += f"\n\t{f} == 0 (mod {self.modulus_name}^{p_deg})"
        for f in self.int_rels:
            s += f"\n\t{f} == 0"
        if self.divisor_bound == 1:
            s += f"\nwhere {self.modulus_name} is known,"
        else:
            if self.mod_mult_modulus_degree == 1:
                mods = f"{self.modulus_name}"
            else:
                mods = f"{self.modulus_name}^{self.mod_mult_modulus_degree}"
            s += f"\nwhere {self.modulus_multiple_name} is a known multiple of {mods}"
            s += f" and {self.modulus_name} > {self.modulus_multiple_name}^{self.divisor_bound},"
        s += " and the small roots satisfy"
        for xi, bi in zip(ring.gens(), self.bounds):
            s += f"\n\t|{xi}| < {self.modulus_name}^" + "{" + str(bi) + "}"
        return s
