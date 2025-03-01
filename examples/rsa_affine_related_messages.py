import logging
import math

from sage.all import randrange, var

import common
import cuso


def generate_challenge(modulus_len, related_bits_fraction, seed=None):
    priv, pub = common.generate_rsa_key(modulus_len, seed=seed)

    n, e = pub

    if related_bits_fraction > 1 / e**2:
        logging.warning(
            "Setting related_bits_fraction to %f exceeds theoretical bound 1 / %d",
            related_bits_fraction,
            e**2,
        )

    pt_1 = int(randrange(n))
    r_bitlen = math.floor(modulus_len * related_bits_fraction)
    related_value = randrange(1 << r_bitlen)

    pt_2 = (pt_1 + related_value) % n

    ct_1 = pow(pt_1, e, n)
    ct_2 = pow(pt_2, e, n)

    challenge = (pub, ct_1, ct_2, r_bitlen)
    solution = (priv, pt_1, pt_2)

    return solution, challenge


def solve_challenge(challenge, solution=None):
    (n, e), ct_1, ct_2, r_bitlen = challenge

    # We represent unknowns with variables.
    # We could equivalently set m and r to the generators of
    # PolynomialRing(ZZ, 'm, r')
    m, r = var("m, r")
    pt_1 = m
    pt_2 = m + r

    # Set relations using RSA equation
    relations = [
        ct_1 == pt_1**e,
        ct_2 == pt_2**e,
    ]

    # Specify bounds as a dictionary of (lower_bound, upper_bound) pairs
    bounds = {
        m: (0, n),
        r: (0, 1 << r_bitlen),
    }

    # When testing, it is often helpful to provide the expected solution
    # as a list of dictionaries
    if solution is not None:
        _priv, _pt_1, _pt_2 = solution
        _m = _pt_1
        _r = (_pt_2 - _pt_1) % n
        expected_solution = [
            {
                m: _m,
                r: _r,
            }
        ]
    else:
        expected_solution = None

    roots = cuso.find_small_roots(
        relations=relations,
        bounds=bounds,
        modulus=n,
        expected_solution=expected_solution,
        allow_partial_solutions=False,
    )

    assert len(roots) > 0
    pt_1 = int(pt_1.subs(roots[0]))
    pt_2 = int(pt_2.subs(roots[0]))

    # Check that these are the correct plaintexts
    assert pow(pt_1, e, n) == ct_1
    assert pow(pt_2, e, n) == ct_2
    print("Successfully recovered affine related plaintexts.")
    return pt_1, pt_2


def main():
    logging.basicConfig(level=logging.INFO)
    soln, chal = generate_challenge(2048, 0.084)
    solve_challenge(challenge=chal, solution=soln)


if __name__ == "__main__":
    main()
