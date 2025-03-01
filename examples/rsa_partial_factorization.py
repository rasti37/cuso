import logging
import math

from sage.all import var

import common
import cuso


def generate_challenge(modulus_len, unknown_fraction, seed=None):
    priv, pub = common.generate_rsa_key(modulus_len, seed=seed)

    n, e = pub
    p, q, *_ = priv

    if unknown_fraction > 0.25:
        logging.warning(
            "Setting unknown_fraction to %f exceeds theoretical bound 1 / 4",
            unknown_fraction,
        )

    # We leak some number of least significant bits of p.
    # This is to show that the library works even if the input polynomial is not monic.
    num_unknown_msbs = math.floor(modulus_len * unknown_fraction)
    num_known_lsbs = (modulus_len // 2) - num_unknown_msbs

    p_lsb_value = p % (2**num_known_lsbs)

    challenge = (pub, p_lsb_value, num_known_lsbs)
    solution = (priv,)

    return solution, challenge


def solve_challenge(challenge, solution=None):
    (n, e), p_lsb_value, num_known_lsbs = challenge

    # We represent unknowns with variables.
    # We could equivalently set x to the generator of
    # PolynomialRing(ZZ, 'x')
    x = var("x")
    p_expr = x * 2**num_known_lsbs + p_lsb_value
    p_len = n.bit_length() // 2

    # The library allows flexibility in specifying input relations.
    # Here, we could also set the relations to
    #     p = var('p')
    #     [
    #       cuso.relation.Relation(p_expr, p),
    #       cuso.relation.Relation(N, p)
    #     ]
    # which means there are two relations with symbolic modulus given by variable p
    relations = [p_expr]

    # Specify bounds on x as a dictionary of (lower_bound, upper_bound) pairs
    num_unk_msbs = p_len - num_known_lsbs
    bounds = {
        x: (0, 2**num_unk_msbs),
    }

    # When testing, it is often helpful to provide the expected solution
    # as a list of dictionaries
    if solution is not None:
        ((_p, *_),) = solution
        _msb_value = (_p - p_lsb_value) >> num_known_lsbs
        expected_solution = [{x: _msb_value, "p": _p}]
    else:
        expected_solution = None

    roots = cuso.find_small_roots(
        relations=relations,
        bounds=bounds,
        modulus="p",
        modulus_multiple=n,
        modulus_lower_bound=1 << (p_len - 1),
        modulus_upper_bound=1 << p_len,
        expected_solution=expected_solution,
    )

    assert len(roots) > 0
    p = int(roots[0]["p"])

    # Check that we got the correct plaintext
    assert 1 < p < n and n % p == 0
    print("Successfully recovered RSA factorization.")
    q = n // p
    return p, q


def main():
    logging.basicConfig(level=logging.INFO)
    soln, chal = generate_challenge(2048, 0.24)
    solve_challenge(challenge=chal, solution=soln)


if __name__ == "__main__":
    main()
