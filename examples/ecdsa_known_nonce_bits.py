from hashlib import sha256
import logging
import math
import os

from sage.all import var, inverse_mod
import ecdsa
from ecdsa.util import PRNG, sigdecode_string
from ecdsa.keys import _truncate_and_convert_digest

import cuso


def get_rsh(vk, msg: bytes, sig: bytes):
    """Recover the values of r, and s, and h from a message and signature pair.

    Args:
        vk: Verification key
        msg_i (bytes): message bytes
        sig_i (bytes): signature bytes
    """
    # Extract out r, s, and h
    r, s = sigdecode_string(sig, vk.curve.order)
    h = _truncate_and_convert_digest(vk.default_hashfunc(msg).digest(), vk.curve, True)
    return r, s, h


def generate_challenge(modulus_len, num_samples, known_fraction, seed=None):
    if num_samples * known_fraction < 1:
        logging.warning(
            "Leaking %f of the nonce bits for %u samples does not meet theoretical requirement",
            known_fraction,
            num_samples,
        )

    rng = PRNG(seed) if seed is not None else os.urandom
    curve = {
        112: ecdsa.SECP112r1,
        128: ecdsa.SECP128r1,
        160: ecdsa.SECP160r1,
        256: ecdsa.SECP256k1,
    }[modulus_len]
    hashfunc = sha256

    n = curve.order
    sk = ecdsa.SigningKey.generate(curve=curve, hashfunc=hashfunc, entropy=rng)
    vk = sk.get_verifying_key()

    x = sk.privkey.secret_multiplier

    num_known_lsbs = math.ceil(known_fraction * modulus_len)

    samples = []
    root = {"x": x}
    for i in range(num_samples):
        msg_i = rng(32)
        sig_i = sk.sign(msg_i, entropy=rng)

        r_i, s_i, h_i = get_rsh(vk, msg_i, sig_i)
        # Compute k_i using
        # s == k^-1 (h + rx)
        k_i = int(inverse_mod(s_i, n) * (h_i + r_i * x) % n)

        k_i_lsb = k_i % (2**num_known_lsbs)
        k_i_msb = (k_i - k_i_lsb) >> num_known_lsbs

        root[f"k_{i}_msb"] = k_i_msb

        samples += [(msg_i, sig_i, k_i_lsb, num_known_lsbs)]

    challenge = (vk, samples)
    solution = (sk, root)

    return solution, challenge


def solve_challenge(challenge, solution=None):
    vk, samples = challenge
    n = int(vk.curve.order)

    # We represent unknowns with variables.
    x = var("x")  # Represents secret key
    relations = []
    bounds = {x: (0, n)}

    for i, (msg_i, sig_i, k_i_lsb, num_known_lsbs) in enumerate(samples):
        assert vk.verify(sig_i, msg_i)

        r_i, s_i, h_i = get_rsh(vk, msg_i, sig_i)
        k_i_msb = var(f"k_{i}_msb")

        k_i = k_i_msb * 2**num_known_lsbs + k_i_lsb

        # ECDSA equation s == k^-1 (h + rx)
        rel = s_i * k_i == h_i + r_i * x

        relations += [rel]
        num_unknown_msbs = n.bit_length() - num_known_lsbs
        bounds[k_i_msb] = (0, 2**num_unknown_msbs)

    # When testing, it is often helpful to provide the expected solution
    if solution is not None:
        _sk, _root = solution
        expected_solution = [
            _root,
        ]
    else:
        expected_solution = None

    roots = cuso.find_small_roots(
        relations=relations,
        bounds=bounds,
        modulus=n,
        expected_solution=expected_solution,
    )

    assert len(roots) > 0
    # Check that we got the correct private key
    x = roots[0][x]
    sk = ecdsa.SigningKey.from_secret_exponent(
        x, vk.curve, hashfunc=vk.default_hashfunc
    )
    assert sk.verifying_key == vk
    print("Successfully recovered ECDSA signing key.")
    return sk


def main():
    logging.basicConfig(level=logging.INFO)
    soln, chal = generate_challenge(256, 25, 0.05)
    solve_challenge(challenge=chal, solution=soln)


if __name__ == "__main__":
    main()
