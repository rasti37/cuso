import logging
import math

from sage.all import randrange, var

import common
import cuso


def generate_challenge(modulus_len: int, unknown_fraction: float, seed=None) -> tuple:
    _, pub = common.generate_rsa_key(modulus_len, seed=seed)

    n, e = pub

    if unknown_fraction > 1 / e:
        logging.warning(
            "Setting unknown_fraction to %f exceeds theoretical bound 1 / %d",
            unknown_fraction,
            e,
        )

    plaintext_value = int(randrange(n))
    ciphertext_value = pow(plaintext_value, e, n)
    # We pick a random plaintext and leak some number of least significant bits.
    # This is to show that the library works even if the input polynomial is not monic.
    num_unknown_msbs = math.floor(modulus_len * unknown_fraction)
    num_known_lsbs = modulus_len - num_unknown_msbs

    lsb_value = plaintext_value % (2**num_known_lsbs)

    challenge = (pub, ciphertext_value, lsb_value, num_known_lsbs)
    solution = (plaintext_value,)

    return solution, challenge


def solve_challenge(challenge, solution=None):
    (n, e), ct, lsb, lsb_len = challenge

    # We represent unknowns with variables.
    # We could equivalently set x to the generator of
    # PolynomialRing(ZZ, 'x')
    x = var("x")
    pt = x * 2**lsb_len + lsb

    # The library allows flexibility in specifying input relations.
    # Here, all the following are equivalent
    #     pt**e == ct
    #     pt**e - ct
    #     [pt**e - ct]
    #     cuso.relation.Relation(pt**e - ct, n)
    #     cuso.relation_set.RelationSet([cuso.relation.Relation(pt**e - ct, n)])
    relations = pt**e == ct

    # Specify bounds as a dictionary of (lower_bound, upper_bound) pairs
    num_unk_msbs = n.bit_length() - lsb_len
    bounds = {
        x: (0, 2**num_unk_msbs),
    }

    # When testing, it is often helpful to provide the expected solution
    # as a list of dictionaries
    if solution is not None:
        (_pt,) = solution
        _msb_value = (_pt - lsb) >> lsb_len
        expected_solution = [{x: _msb_value}]
    else:
        expected_solution = None

    roots = cuso.find_small_roots(
        relations=relations,
        bounds=bounds,
        modulus=n,
        expected_solution=expected_solution,
    )

    assert len(roots) > 0
    # Substitute the recovered root back into the plaintext expression
    pt_value = int(pt.subs(roots[0]))

    # Check that we got the correct plaintext
    assert pow(pt_value, e, n) == ct
    print("Successfully recovered stereotyped RSA plaintext.")
    return pt_value


def main():
    logging.basicConfig(level=logging.INFO)
    soln, chal = generate_challenge(2048, 0.30)
    solve_challenge(challenge=chal, solution=soln)


if __name__ == "__main__":
    main()
