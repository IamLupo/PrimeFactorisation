#!/usr/bin/env python3

import math
import random
import time

from sympy import Matrix


# ========================================================================
# START EXPERIMENT 128
# LLL / Coppersmith-style recovery from partial factor residue
# ========================================================================

print("=" * 72)
print("START EXPERIMENT 128")
print("LLL / Coppersmith-style recovery from partial factor residue")
print("=" * 72)


# ------------------------------------------------------------------------
# RADICES
#
# We deliberately use several small coprime radices.
#
# Their product M is the amount of known information about p.
#
# The central question:
#
#       p = A + M*x
#
# with x relatively small.
#
# Can LLL recover x without enumerating x?
# ------------------------------------------------------------------------

R1 = [17, 43, 47]

M = math.prod(R1)

print()
print("R1 =", R1)
print("M  =", M)


# ------------------------------------------------------------------------
# MILLER-RABIN
# ------------------------------------------------------------------------

def is_probable_prime(n):

    if n < 2:
        return False

    small = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    )

    for p in small:

        if n == p:
            return True

        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        d //= 2
        s += 1

    bases = (
        2,
        325,
        9375,
        28178,
        450775,
        9780504,
        1795265022,
    )

    for a in bases:

        if a % n == 0:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):

            x = (x * x) % n

            if x == n - 1:
                break

        else:
            return False

    return True


# ------------------------------------------------------------------------
# PRIME GENERATOR
# ------------------------------------------------------------------------

def random_prime(bits):

    while True:

        x = random.getrandbits(bits)

        x |= 1
        x |= 1 << (bits - 1)

        if is_probable_prime(x):
            return x


# ------------------------------------------------------------------------
# SEMIPRIME
# ------------------------------------------------------------------------

def random_semiprime(bits):

    pb = bits // 2
    qb = bits - pb

    while True:

        p = random_prime(pb)
        q = random_prime(qb)

        if p == q:
            continue

        p, q = sorted((p, q))

        n = p * q

        if n.bit_length() == bits:
            return n, p, q


# ------------------------------------------------------------------------
# CRT
# ------------------------------------------------------------------------

def crt_pairwise(residues, moduli):

    x = 0
    modulus = 1

    for a, m in zip(residues, moduli):

        inv = pow(modulus, -1, m)

        t = ((a - x) * inv) % m

        x += modulus * t
        modulus *= m
        x %= modulus

    return x, modulus


# ------------------------------------------------------------------------
# TRUE CRT RESIDUE OF p
# ------------------------------------------------------------------------

def factor_residue(p):

    residues = [
        p % r
        for r in R1
    ]

    A, modulus = crt_pairwise(
        residues,
        R1
    )

    return residues, A, modulus


# ------------------------------------------------------------------------
# BRUTE-FORCE x CONTROL
#
# p = A + M*x
#
# Search x directly.
# ------------------------------------------------------------------------

def brute_x_search(n, A, M):

    limit = math.isqrt(n)

    if A == 0:
        x0 = 1
    else:
        x0 = 0

    # First x giving p >= 2.
    if A + M * x0 < 2:

        x0 = (
            (2 - A + M - 1)
            // M
        )

    if A + M * x0 > limit:

        return None, 0

    x = x0
    tested = 0

    while True:

        p = A + M * x

        if p > limit:
            break

        tested += 1

        if n % p == 0:

            return p, tested

        x += 1

    return None, tested


# ------------------------------------------------------------------------
# LLL POLYNOMIAL BASIS
#
# We want a small root x of
#
#     f(x) = A + M*x
#
# modulo an unknown factor p of n.
#
# We construct the simplest degree-1 Coppersmith-style lattice from
#
#     f(x)
#     x*f(x)
#     n
#
# with coefficient scaling X.
#
# This is intentionally an exploratory construction rather than a claim
# that this is an optimal Coppersmith implementation.
#
# The resulting short vectors are interpreted as integer polynomials.
# ------------------------------------------------------------------------

def polynomial_coefficients_from_vector(
    vector,
    n,
    A,
    M,
):
    """
    Basis convention:

        row 0 = n
        row 1 = X*f(x)
        row 2 = X*x*f(x)

    where

        f(x) = A + M*x.

    Coefficients are returned in powers of x.
    """

    c0 = int(vector[0])
    c1 = int(vector[1])
    c2 = int(vector[2])

    # ------------------------------------------------------------
    # n contributes only to x^0.
    #
    # X*f(x) contributes:
    #
    #   X*A + X*M*x
    #
    # X*x*f(x) contributes:
    #
    #   X*A*x + X*M*x^2
    # ------------------------------------------------------------

    return [
        c0 * n + c1 * 0 + c2 * 0,
        c1 * 0,
        c2 * 0,
    ]


# ------------------------------------------------------------------------
# Build lattice.
#
# Rather than trying to interpret raw LLL vectors directly, we use a
# scaled coefficient representation explicitly.
# ------------------------------------------------------------------------

def build_lattice(
    n,
    A,
    M,
    X,
):

    # Polynomial rows:
    #
    # g0 = n
    # g1 = X*f(x)
    # g2 = X*x*f(x)
    #
    # coefficient columns are:
    #
    # x^0, x^1, x^2
    #

    B = Matrix([
        [n,           0,           0],
        [X * A,       X * M,       0],
        [0,           X * A,       X * M],
    ])

    return B


# ------------------------------------------------------------------------
# Evaluate a polynomial at x.
# ------------------------------------------------------------------------

def eval_poly(coeffs, x):

    result = 0

    for c in reversed(coeffs):

        result = result * x + c

    return result


# ------------------------------------------------------------------------
# LLL SEARCH
#
# We don't assume that the first vector is necessarily the solution.
# We inspect every reduced basis vector and small integer combinations
# of pairs of vectors.
# ------------------------------------------------------------------------

def lll_recover(
    n,
    A,
    M,
    X,
):

    B = build_lattice(
        n,
        A,
        M,
        X,
    )

    R = B.lll()

    vectors = []

    for i in range(R.rows):

        vectors.append(
            [
                int(R[i, j])
                for j in range(R.cols)
            ]
        )

    # ------------------------------------------------------------
    # Convert lattice vectors into polynomials.
    #
    # Our basis encoding:
    #
    #   v0*n
    #   v1*X*f
    #   v2*X*x*f
    #
    # therefore:
    #
    #   c0 = v0*n + v1*X*A
    #   c1 = v1*X*M + v2*X*A
    #   c2 = v2*X*M
    # ------------------------------------------------------------

    polys = []

    for v in vectors:

        v0, v1, v2 = v

        coeffs = [
            v0 * n + v1 * X * A,
            v1 * X * M + v2 * X * A,
            v2 * X * M,
        ]

        polys.append(coeffs)

    # ------------------------------------------------------------
    # Also try small combinations of reduced basis vectors.
    # ------------------------------------------------------------

    extra = []

    for i in range(len(vectors)):

        for j in range(i + 1, len(vectors)):

            for alpha in (-2, -1, 1, 2):

                for beta in (-2, -1, 1, 2):

                    v = [
                        alpha * vectors[i][k]
                        + beta * vectors[j][k]
                        for k in range(3)
                    ]

                    v0, v1, v2 = v

                    coeffs = [
                        v0 * n + v1 * X * A,
                        v1 * X * M + v2 * X * A,
                        v2 * X * M,
                    ]

                    extra.append(coeffs)

    polys.extend(extra)

    # ------------------------------------------------------------
    # Root search.
    #
    # This is intentionally restricted to integer roots in the
    # expected interval.
    #
    # For degree 1 we can solve exactly.
    # For degree 2 we use discriminant.
    # ------------------------------------------------------------

    candidates = []

    for coeffs in polys:

        c0, c1, c2 = coeffs

        # --------------------------------------------------------
        # Constant polynomial.
        # --------------------------------------------------------

        if c1 == 0 and c2 == 0:
            continue

        # --------------------------------------------------------
        # Linear polynomial.
        # --------------------------------------------------------

        if c2 == 0:

            if c1 == 0:
                continue

            if (-c0) % c1 == 0:

                x = (-c0) // c1

                candidates.append(x)

            continue

        # --------------------------------------------------------
        # Quadratic.
        # --------------------------------------------------------

        D = c1 * c1 - 4 * c2 * c0

        if D < 0:
            continue

        sD = math.isqrt(D)

        if sD * sD != D:
            continue

        for sign in (-1, 1):

            numerator = -c1 + sign * sD
            denominator = 2 * c2

            if denominator != 0:

                if numerator % denominator == 0:

                    x = numerator // denominator

                    candidates.append(x)

    # ------------------------------------------------------------
    # Deduplicate.
    # ------------------------------------------------------------

    candidates = sorted(set(candidates))

    return candidates, vectors, polys


# ------------------------------------------------------------------------
# VERIFY CANDIDATE
# ------------------------------------------------------------------------

def verify_x(
    n,
    A,
    M,
    x,
):

    if x < 0:
        return None

    p = A + M * x

    if p <= 1:
        return None

    if n % p != 0:
        return None

    q = n // p

    return p, q


# ------------------------------------------------------------------------
# OPTIONAL FULL CARRY CHECK
# ------------------------------------------------------------------------

def carry_state(
    p,
    q,
    r,
    s,
):

    k, a = divmod(p, r)
    ell, b = divmod(q, s)

    c1 = (k * b) // s
    c2 = (ell * a) // r

    beta = (k * b) % s
    alpha = (ell * a) % r

    c3 = (
        r * beta
        + s * alpha
        + a * b
    ) // (r * s)

    E = c1 + c2 + c3
    Q = (p * q) // (r * s)

    return {
        "k": k,
        "ell": ell,
        "a": a,
        "b": b,
        "E": E,
        "Q": Q,
        "K": Q - E,
    }


def verify_carry_grid(
    p,
    q,
):

    for r in R1:

        C = carry_state(
            p,
            q,
            r,
            19,
        )

        if C["K"] != C["k"] * C["ell"]:
            return False

    return True


# ------------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------------

random.seed(128)


# We want cases where M is meaningful but does NOT itself uniquely
# determine p below sqrt(n).
#
# With 48-bit n and M=34357, there can still be hundreds of x values.
TEST_BITS = [
    30,
    36,
    42,
    48,
]


for bits in TEST_BITS:

    print()
    print("=" * 72)
    print("GENERATING", bits, "-BIT SEMIPRIME")
    print("=" * 72)

    n, true_p, true_q = random_semiprime(bits)

    residues, A, modulus = factor_residue(true_p)

    true_x = (
        true_p - A
    ) // M

    x_bound = (
        math.isqrt(n) - A
    ) // M

    print()
    print("-" * 72)

    print("n bits =", n.bit_length())
    print("n      =", n)
    print("true p =", true_p)
    print("true q =", true_q)

    print()
    print("PARTIAL FACTOR INFORMATION")

    print("residues =", residues)
    print("A mod M  =", A)
    print("M        =", modulus)

    print()
    print("TRUE x")
    print("x =", true_x)
    print("x bound =", x_bound)

    # ------------------------------------------------------------
    # Direct x scan.
    # ------------------------------------------------------------

    print()
    print("BRUTE-FORCE x SCAN")

    t0 = time.perf_counter()

    p_bf, tested_bf = brute_x_search(
        n,
        A,
        M,
    )

    t1 = time.perf_counter()

    brute_time = t1 - t0

    print(
        "x tested =",
        tested_bf
    )

    print(
        "runtime =",
        f"{brute_time:.6f}",
        "s"
    )

    print(
        "found =",
        p_bf == true_p
    )

    # ------------------------------------------------------------
    # Try several lattice scales.
    #
    # X is the expected magnitude of x.
    # ------------------------------------------------------------

    print()
    print("LLL RECOVERY")

    scales = [
        max(1, x_bound),
        max(1, x_bound // 2),
        max(1, 2 * x_bound),
        max(1, 4 * x_bound),
    ]

    successful_scales = []

    total_lll_time = 0.0

    for X in scales:

        t0 = time.perf_counter()

        candidates, vectors, polys = lll_recover(
            n,
            A,
            M,
            X,
        )

        t1 = time.perf_counter()

        elapsed = t1 - t0
        total_lll_time += elapsed

        found = False
        recovered = None

        for x in candidates:

            result = verify_x(
                n,
                A,
                M,
                x,
            )

            if result is None:
                continue

            p, q = result

            if {
                p,
                q
            } == {
                true_p,
                true_q
            }:

                found = True
                recovered = (
                    x,
                    p,
                    q,
                )

                break

        print()
        print("X =", X)
        print("reduced basis vectors =", len(vectors))
        print("candidate roots =", len(candidates))
        print("runtime =", f"{elapsed:.6f}", "s")
        print("TRUE FACTOR RECOVERED =", found)

        if found:

            successful_scales.append(X)

            x, p, q = recovered

            print("recovered x =", x)
            print("recovered p =", p)
            print("recovered q =", q)

    # ------------------------------------------------------------
    # Summary.
    # ------------------------------------------------------------

    print()
    print("SUMMARY")

    print(
        "brute-force x candidates =",
        tested_bf
    )

    print(
        "LLL total runtime =",
        f"{total_lll_time:.6f}",
        "s"
    )

    print(
        "successful LLL scales =",
        successful_scales
    )

    # ------------------------------------------------------------
    # Carry verification when successful.
    # ------------------------------------------------------------

    if successful_scales:

        print()
        print("CARRY VERIFICATION")

        print(
            "carry grid verified =",
            verify_carry_grid(
                true_p,
                true_q,
            )
        )


# ========================================================================
# FINISHED EXPERIMENT 128
# ========================================================================

print()
print("=" * 72)
print("FINISHED EXPERIMENT 128")
print("=" * 72)
