import logging
import math

from sage.all import var, randrange, set_random_seed, random_prime

import cuso


def generate_challenge(
    sample_len, num_samples, divisor_fraction, err_fraction, seed=None
):
    set_random_seed(seed)

    divisor_len = math.ceil(sample_len * divisor_fraction)
    err_len = math.floor(divisor_len * err_fraction)

    p = random_prime(1 << divisor_len, proof=False, lbound=1 << (divisor_len - 1))

    samples = []
    root = {"p": p}
    for i in range(num_samples):
        q_i = randrange(1 << (sample_len - divisor_len))

        r_i = randrange(-(1 << err_len), 1 << err_len)
        a_i = p * q_i + r_i

        root[f"r_{i}"] = r_i

        samples += [a_i]

    challenge = (samples, divisor_len, err_len)
    solution = (root,)

    return solution, challenge


def solve_challenge(challenge, solution=None):
    samples, divisor_len, err_len = challenge

    # We represent unknowns with variables.
    relations = []
    bounds = {}

    for i, a_i in enumerate(samples):
        r_i = var(f"r_{i}")
        rel = a_i - r_i

        relations += [rel]
        bounds[r_i] = (-(2**err_len), 2**err_len)

    # When testing, it is often helpful to provide the expected solution
    if solution is not None:
        (_root,) = solution
        expected_solution = [
            _root,
        ]
    else:
        expected_solution = None

    roots = cuso.find_small_roots(
        relations=relations,
        bounds=bounds,
        modulus="p",
        modulus_lower_bound=1 << (divisor_len - 1),
        modulus_upper_bound=1 << divisor_len,
        expected_solution=expected_solution,
    )

    assert len(roots) > 0
    p = roots[0]["p"]

    print("Successfully recovered approximate divisor.")
    return p


def main():
    logging.basicConfig(level=logging.INFO)
    soln, chal = generate_challenge(1024, 4, 0.3, 0.35)
    solve_challenge(challenge=chal, solution=soln)


if __name__ == "__main__":
    main()
