import logging
import math

from sage.all import randrange, var, inverse_mod

import common
import cuso


def generate_challenge(modulus_len, private_exponent_fraction, seed=None):
    theory_bound = 1 - math.sqrt(2) / 2
    if private_exponent_fraction > theory_bound:
        logging.warning(
            "Setting private_exponent_fraction to %f exceeds theoretical bound %f",
            private_exponent_fraction,
            theory_bound,
        )
    priv_exp_len = math.floor(modulus_len * private_exponent_fraction)
    priv, pub = common.generate_rsa_key(modulus_len, d_len=priv_exp_len, seed=seed)

    challenge = (pub, priv_exp_len)
    solution = (priv,)

    return solution, challenge


def solve_challenge_bd99(challenge, solution=None):
    (n, e), priv_exp_len = challenge

    # We represent unknowns with variables.
    # We could equivalently set k and s to the generators of
    # PolynomialRing(ZZ, 'k, s')
    k, s = var("k, s")

    # According to the RSA equations,
    # e * d == 1 + k * phi(n)
    # e * d == 1 + k * (n - (p + q) + 1)
    # 0 == 1 + k * (n - s + 1) (mod e)
    relation = 1 + k * (n - s + 1)

    # Specify bounds as a dictionary of (lower_bound, upper_bound) pairs
    # s = p + q
    prime_lbound = 1 << ((n.bit_length() // 2) - 1)
    prime_ubound = 1 << (n.bit_length() // 2)
    s_bound = (2 * prime_lbound, 2 * prime_ubound)
    # k ~= e * d / n ~= d
    bounds = {
        k: (1, 1 << priv_exp_len),
        s: s_bound,
    }

    # When testing, it is often helpful to provide the expected solution
    # as a list of dictionaries
    if solution is not None:
        ((_p, _q, _d),) = solution
        _s = _p + _q
        _k = (e * _d - 1) // (n - _s + 1)
        expected_solution = [
            {
                k: _k,
                s: _s,
            }
        ]
    else:
        expected_solution = None

    roots = cuso.find_small_roots(
        relations=[relation],
        bounds=bounds,
        modulus=e,
        expected_solution=expected_solution,
    )

    assert len(roots) > 0
    k_value = int(roots[0][k])
    s_value = int(roots[0][s])

    # e * d == 1 + k * (n - (p + q) + 1)
    d = (1 + k_value * (n - s_value + 1)) // e

    # Check that this is a valid decryption exponent
    r = randrange(1, n)
    assert pow(pow(r, e, n), d, n) == r
    print("Successfully recovered small RSA private exponent", d)
    return d


def solve_challenge_hm10(challenge, solution=None):
    (n, e), priv_exp_len = challenge

    # We represent unknowns with variables.
    # We could equivalently set k and s to the generators of
    # PolynomialRing(ZZ, 'k, s')
    k, s = var("k, s")

    # According to the RSA equations,
    # e * d == 1 + k * phi(n)
    # e * d == 1 + k * (n - (p + q) + 1)
    # 0 == 1 + k * (n - s + 1) (mod e)
    relation = 1 + k * (n - s + 1)

    # Herrmann and May's approach uses unraveled linearization
    ul = [1 - k * s]

    # Specify bounds as a dictionary of (lower_bound, upper_bound) pairs
    # s = p + q
    prime_lbound = 1 << ((n.bit_length() // 2) - 1)
    prime_ubound = 1 << (n.bit_length() // 2)
    s_bound = (2 * prime_lbound, 2 * prime_ubound)
    # k ~= e * d / n ~= d
    bounds = {
        k: (1, 1 << priv_exp_len),
        s: s_bound,
    }

    # When testing, it is often helpful to provide the expected solution
    # as a list of dictionaries
    if solution is not None:
        ((_p, _q, _d),) = solution
        _s = _p + _q
        _k = (e * _d - 1) // (n - _s + 1)
        expected_solution = [
            {
                k: _k,
                s: _s,
            }
        ]
    else:
        expected_solution = None

    roots = cuso.find_small_roots(
        relations=[relation],
        bounds=bounds,
        modulus=e,
        unraveled_linearization_relations=ul,
        expected_solution=expected_solution,
    )

    assert len(roots) > 0
    k_value = int(roots[0][k])
    s_value = int(roots[0][s])

    # e * d == 1 + k * (n - (p + q) + 1)
    d = (1 + k_value * (n - s_value + 1)) // e

    # Check that this is a valid decryption exponent
    r = randrange(1, n)
    assert pow(pow(r, e, n), d, n) == r
    print("Successfully recovered small RSA private exponent", d)
    return d


def solve_challenge_cuso(challenge, solution=None):
    (n, e), priv_exp_len = challenge

    # We represent unknowns with variables.
    k, p, q = var("k, p, q")

    # According to the RSA equations,
    # e * d == 1 + k * phi(n)
    # e * d == 1 + k * (p - 1) * (q - 1)
    # 0 == 1 + k * (p - 1) * (q - 1) (mod e)
    phi = (p - 1) * (q - 1)
    relation_exp = 1 + k * phi
    # We also have n == p * q
    relation_n = n - p * q
    relations = [
        relation_exp,
        relation_n,
    ]
    # The first relation is modulo e, the second doesn't have a modulus.
    modulus = [
        e,
        None,
    ]

    # Specify bounds as a dictionary of (lower_bound, upper_bound) pairs
    # s = p + q
    prime_lbound = 1 << ((n.bit_length() // 2) - 1)
    prime_ubound = 1 << (n.bit_length() // 2)
    # k ~= e * d / n ~= d
    bounds = {
        k: (1, 1 << priv_exp_len),
        p: (prime_lbound, prime_ubound),
        q: (prime_lbound, prime_ubound),
    }

    # When testing, it is often helpful to provide the expected solution
    # as a list of dictionaries
    if solution is not None:
        ((_p, _q, _d),) = solution
        _k = (e * _d - 1) // ((_p - 1) * (_q - 1))
        expected_solution = [
            {
                k: _k,
                p: _p,
                q: _q,
            }
        ]
    else:
        expected_solution = None

    roots = cuso.find_small_roots(
        relations=relations,
        bounds=bounds,
        modulus=modulus,
        expected_solution=expected_solution,
    )

    assert len(roots) > 0
    phi = phi.subs(roots[0])
    d = int(inverse_mod(e, phi))

    # Check that this is a valid decryption exponent
    r = randrange(1, n)
    assert pow(pow(r, e, n), d, n) == r
    print("Successfully recovered small RSA private exponent", d)
    return d


def main():
    logging.basicConfig(level=logging.INFO)
    soln, chal = generate_challenge(2048, 0.26)
    # We implement three approaches to solve this problem.
    # The first is Boneh and Durfee's approach from 1999.
    # The second is Herrmann and May's approach from 2010, which uses unraveled linearization.
    # The third is a more natural way of expressing the problem using cuso.

    # solve_challenge_bd99(challenge=chal, solution=soln)
    # solve_challenge_hm10(challenge=chal, solution=soln)
    solve_challenge_cuso(challenge=chal, solution=soln)


if __name__ == "__main__":
    main()
