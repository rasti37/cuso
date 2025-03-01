import logging
import math

from sage.all import var, inverse_mod, randrange, set_random_seed, random_prime

import cuso


def generate_challenge(modulus_len, num_samples, known_fraction, seed=None):
    set_random_seed(seed)

    p = random_prime(1 << modulus_len, proof=False, lbound=1 << (modulus_len - 1))

    alpha = randrange(p)

    num_known_msbs = math.ceil(known_fraction * modulus_len)
    num_unknown_lsbs = modulus_len - num_known_msbs

    samples = []
    root = {"alpha": alpha}
    for i in range(num_samples):
        x_i = randrange(p)

        t = int(inverse_mod(x_i + alpha, p))
        lsb_i = t % (2**num_unknown_lsbs)
        b_i = (t >> num_unknown_lsbs) << num_unknown_lsbs

        root[f"eps_{i}"] = lsb_i

        samples += [(x_i, b_i, num_unknown_lsbs)]

    challenge = (p, samples)
    solution = (root,)

    return solution, challenge


def solve_challenge(challenge, solution=None):
    p, samples = challenge

    # We represent unknowns with variables.
    alpha = var("alpha")  # Represents hidden number
    relations = []
    bounds = {alpha: (0, p)}

    for i, (x_i, b_i, num_unknown_lsbs) in enumerate(samples):
        eps_i = var(f"eps_{i}")
        # ECDSA equation b_i = MSB(1 / (alpha + x_i))
        rel = (b_i + eps_i) * (alpha + x_i) == 1

        relations += [rel]
        bounds[eps_i] = (0, 2**num_unknown_lsbs)

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
        modulus=p,
        expected_solution=expected_solution,
    )

    assert len(roots) > 0
    alpha = roots[0][alpha]
    print("Successfully recovered hidden number.")
    return alpha


def main():
    logging.basicConfig(level=logging.INFO)
    soln, chal = generate_challenge(256, 5, 0.53)
    solve_challenge(challenge=chal, solution=soln)


if __name__ == "__main__":
    main()
