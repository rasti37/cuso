""" Common utilities for cuso example scripts
"""

from sage.all import random_prime, set_random_seed, randrange, gcd, inverse_mod


def generate_rsa_key(bit_length, e_value=None, e_len=None, d_len=None, seed=None):
    """Generate an RSA key

    Args:
        bit_length (int): bit length of the RSA modulus
        e_value (int, optional): public exponent. Defaults to None.
        e_len (int, optional): length of public exponent. Defaults to None.
        d_len (int, optional): length of private exponent. Defaults to None.
        seed (int, optional): RNG seed. Defaults to None.

    Returns:
        _type_: tuple of private key, public key
    """
    p_len = bit_length // 2
    q_len = bit_length - p_len

    if seed is not None:
        set_random_seed(seed)

    p_lbound = 3 << (p_len - 2)
    p_ubound = 1 << p_len
    q_lbound = 3 << (q_len - 2)
    q_ubound = 1 << q_len

    # Check exponent specification
    if e_value is not None:
        assert e_len is None and d_len is None
    elif e_len is not None:
        assert e_value is None and d_len is None
    elif d_len is not None:
        assert e_len is None and e_value is None

    while True:
        if e_value is not None:
            e = e_value
        elif e_len is not None:
            e = randrange(1 << (e_len - 1), 1 << e_len)
        elif d_len is not None:
            # We swap d and e later
            e = randrange(1 << (d_len - 1), 1 << d_len)
        else:
            e = 3

        p = int(random_prime(p_ubound, proof=False, lbound=p_lbound))
        q = int(random_prime(q_ubound, proof=False, lbound=q_lbound))

        n = p * q

        if (
            p.bit_length() == p_len
            and q.bit_length() == q_len
            and n.bit_length() == bit_length
            and gcd(e, p - 1) == 1
            and gcd(e, q - 1) == 1
        ):
            break

    phi = (p - 1) * (q - 1)
    d = int(inverse_mod(e, phi))
    if d_len is not None:
        # Swap d and e
        d, e = e, d

    pub = (n, e)
    priv = (p, q, d)
    return priv, pub
