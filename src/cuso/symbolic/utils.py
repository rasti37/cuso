def monom_no_N(m):
    N = m.parent().gens()[0]
    return m.subs({N: 1})

def lm_no_N(f):
    # f is in the ring K[N, x_1, x_2, ...]
    # return the leading monomial in terms of (x_1, x_2, ...)
    return lt_no_N(f)[1]

def lc_no_N(f):
    return lt_no_N(f)[0]

def lt_no_N(f):
    ring = f.parent()
    max_monom = 0
    max_coef = 0
    for coef, monom in zip(f.coefficients(), f.monomials()):
        m_no_N = monom_no_N(monom)
        coef *= monom / m_no_N

        if m_no_N > max_monom:
            max_monom = m_no_N
            max_coef = coef
    return (max_coef, ring(max_monom))

def coeff_no_N(f, m):
    # Get the coefficient of m in f, ignoring N
    c = 0
    for coef, monom in zip(f.coefficients(), f.monomials()):
        if monom.degrees()[1:] == m.degrees()[1:]:
            c += (monom / m) * coef
    return c