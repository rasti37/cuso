# cuso #

cuso is a (Copper)smith (So)lver which automatically finds small roots of multivariate polynomials. It is intended for use by CTF players and researchers, and it abstracts away as many of the complex details about lattices and Coppersmith's method as possible. cuso is the corresponding implementation for the paper "[Solving Multivariate Coppersmith Problems with Known Moduli](https://eprint.iacr.org/2024/1577)."

## Multivariate Coppersmith Problems
Many cryptanalytic problems involve finding a bounded solution to a system of modular or integer polynomials. This includes factoring RSA moduli when 1/4 of the bits are known, recovering RSA plaintexts with fixed padding, and more. Formally, a multivariate Coppersmith problem is defined by a system of polynomials with integer coefficients. Each polynomial $f_i$ has the form

$$ 
\begin{align*}
f_i(x_1, \dots, x_\ell) &\equiv 0 \pmod{N_i} & \text{ where } N_i\text{ is known, or} \\
f_i(x_1, \dots, x_\ell) &\equiv 0 \pmod{p_i} & \text{ where } p_i\text{ is unknown, or} \\
f_i(x_1, \dots, x_\ell) &= 0 & \text {with no modulus}.
\end{align*}
$$

In some cases when $p_i$ is unknown, a multiple $N_i$ of $p_i$ is known. A Coppersmith problem also includes a collection of integer bounds $(X_1, \dots, X_\ell)$. The goal is to find all solutions $(r_1, \dots, r_\ell)$ that satisfy all of the polynomials and the bounds:

$$ |r_i| < X_i \quad 1 \le i \le \ell. $$

## Usage
Example scripts are found in the [examples](examples) directory. Let's say we are trying to recover the factorization of RSA modulus $N$ given the least significant bits of the factor $p$.
~~~python
import cuso
from sage.all import var

N = 7803960993055714764790127684076211863657988334421459147942569078964485679825779792612991909294331632362413082874625938556130231837335088680459104018893959
p_lsb = 779432292600303509528505569280402295375268661469
num_lsbs = 160

# x represents the unknown lsbs
x = var('x')
# f(x) == 0 (mod p) when x is the 96 most significant bits of p
f = 2**160 * x + p_lsb

# Only one polynomial in our system
relations = [f]
# We give lower and upper bounds of x
bounds = {x: (0, 2**96)}
roots = cuso.find_small_roots(
    relations,
    bounds,
    modulus = "p",
    modulus_multiple = N,
    modulus_lower_bound = 2**255,
)
assert len(roots) > 0
p = roots[0]["p"]
assert N % p == 0
print(f"Recovered factor p = {p}")
~~~

We can also handle the more complicated, multivariate problem of recovering a small RSA secret exponent.
~~~python
import cuso
from sage.all import var, inverse_mod

N = 7803960993055714764790127684076211863657988334421459147942569078964485679825779792612991909294331632362413082874625938556130231837335088680459104018893959
e = 7338455574804691890149651420174638680950093405994546367015621732049630822037090156191940077914112346956150999747345628657616947167917355899104335150736479
priv_exp_len = 136

p, q, k = var('p, q, k')
# By the RSA equation,
#     e * d = 1 + k*(p - 1)*(q - 1)
#         N = p * q
# where d, k, p, and q are unknown.
# The first relation is modulo e;
# the second doesn't have a modulus.
phi = (p - 1) * (q - 1)
relations = [
    0 == 1 + k * phi,
    N == p * q,
]
moduli = [
    e,
    None,
]

# We give lower and upper bounds of p, q, and k.
bounds = {
    p: (2**255, 2**256),
    q: (2**255, 2**256),
    k: (0, 2**priv_exp_len),
}
roots = cuso.find_small_roots(
    relations,
    bounds,
    modulus=moduli,
)
assert len(roots) > 0
phi = phi.subs(roots[0])
d = int(inverse_mod(e, phi))
print(f"Recovered small exponent d = {d}")
~~~

There are also [examples](examples/symbolic/) for computing asymptotic bounds automatically.

## Installation
Install dependencies:
* [SageMath](https://www.sagemath.org/) (version 9.8 or higher is recommended)
* [flatter](https://github.com/keeganryan/flatter.git)
* [msolve](https://github.com/algebraic-solving/msolve)

Install cuso:
```
git clone https://github.com/keeganryan/cuso.git
cd cuso/
pip install .
```

## Details
Coppersmith's [original 1996 work](https://doi.org/10.1007/3-540-68339-9_14) shows how to solve Coppersmith problems involving a single univariate polynomial $f$ modulo $N$, but the technique is notoriously hard to generalize to systems of multivariate polynomials. The main challenge is choosing good _shift polynomials_, or polynomial combinations of the input polynomials which involve overlapping sets of monomials. Previously, achieving the best results for many Coppersmith problems required a deep understanding of Coppersmith's method and handcrafted shift polynomial selection strategies.

cuso is based on the paper "[Solving Multivariate Coppersmith Problems with Known Moduli](https://eprint.iacr.org/2024/1577)," which presents three _automatic_ shift polynomial strategies that meet or exceed the performance of the existing manual strategies. There is a provably optimal strategy, a heuristic graph-based strategy, and a precomptuation strategy. The first two are implemented in cuso.

The provable strategy is based on the theory of Groebner bases over Euclidean domains. Using [Howgrave-Graham's dual construction](https://doi.org/10.1007/bfb0024458), Coppersmith's method heuristically succeeds when it finds enough integer linear combinations of shift polynomials that have small coefficients. The provable strategy produces shift polynomials with the guarantee that their integer linear span contains all of these small Howgrave-Graham polynomials. However, the downside to this is that the total number of shift polynomials may be large, and it might not be possible to efficiently find the small linear combinations using lattice reduction.

The graph-based strategy is an optimization that finds small subsets of the optimal shift polynomials, so lattice reduction becomes less expensive. This optimization is not proven to work in all cases, but it is often extremely useful in practice.

cuso applies these strategies automatically, eliminating the need to manually specify shift polynomial strategies or understand how to parameterize the complicated inner behavior of Coppersmith's method. cuso also implements the symbolic strategies described in the paper, which can be used to heuristically compute the asymptotic bounds for a given problem.

## License
cuso is licensed under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.en.html).
