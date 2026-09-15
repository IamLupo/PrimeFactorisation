#!/usr/bin/env python3

import math
import random
import time


# ========================================================================
# START EXPERIMENT 132
# 2x2 rank-1 determinant collapse
# ========================================================================

print("=" * 72)
print("START EXPERIMENT 132")
print("2x2 rank-1 determinant collapse")
print("=" * 72)


# ------------------------------------------------------------------------
# GRID
# ------------------------------------------------------------------------

R1 = [17, 43]
R2 = [19, 37]

print()
print("R1 =", R1)
print("R2 =", R2)

for r in R1:
    for s in R2:
        print(
            f"R[{r},{s}] = {r*s}"
        )


# ------------------------------------------------------------------------
# MILLER-RABIN
# ------------------------------------------------------------------------

def is_probable_prime(n):

    if n < 2:
        return False

    small_primes = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    )

    for p in small_primes:

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
# RANDOM PRIME
# ------------------------------------------------------------------------

def random_prime(bits):

    while True:

        x = random.getrandbits(bits)

        x |= 1
        x |= 1 << (bits - 1)

        if is_probable_prime(x):
            return x


# ------------------------------------------------------------------------
# RANDOM SEMIPRIME
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
# EXACT CELL
# ------------------------------------------------------------------------

def cell(n, p, q, r, s):

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

    Q = n // (r * s)

    K = Q - E

    # Direct algebraic identity.
    assert K == k * ell

    return {
        "Q": Q,
        "E": E,
        "K": K,
        "k": k,
        "ell": ell,
        "a": a,
        "b": b,
        "c1": c1,
        "c2": c2,
        "c3": c3,
        "alpha": alpha,
        "beta": beta,
    }


# ------------------------------------------------------------------------
# 2x2 GRID
# ------------------------------------------------------------------------

def true_grid(n, p, q):

    g = {}

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            g[(i, j)] = cell(
                n,
                p,
                q,
                r,
                s,
            )

    return g


# ------------------------------------------------------------------------
# DETERMINANT
# ------------------------------------------------------------------------

def determinant_values(g):

    Q00 = g[(0, 0)]["Q"]
    Q01 = g[(0, 1)]["Q"]
    Q10 = g[(1, 0)]["Q"]
    Q11 = g[(1, 1)]["Q"]

    E00 = g[(0, 0)]["E"]
    E01 = g[(0, 1)]["E"]
    E10 = g[(1, 0)]["E"]
    E11 = g[(1, 1)]["E"]

    K00 = g[(0, 0)]["K"]
    K01 = g[(0, 1)]["K"]
    K10 = g[(1, 0)]["K"]
    K11 = g[(1, 1)]["K"]

    DQ = (
        Q00 * Q11
        - Q01 * Q10
    )

    DK = (
        K00 * K11
        - K01 * K10
    )

    # ------------------------------------------------------------
    # Expanded carry-side determinant correction.
    #
    # DQ =
    #
    #   Q00*E11
    # + Q11*E00
    # - Q01*E10
    # - Q10*E01
    # + E00*E11
    # - E01*E10
    #
    # ------------------------------------------------------------

    correction = (
        Q00 * E11
        + Q11 * E00
        - Q01 * E10
        - Q10 * E01
        + E00 * E11
        - E01 * E10
    )

    return {
        "DQ": DQ,
        "DK": DK,
        "correction": correction,
        "expanded_identity": (
            DQ == correction
        ),
    }


# ------------------------------------------------------------------------
# NORMALIZED CARRY DETERMINANT
#
# The raw determinant can be enormous.
#
# We also separate:
#
#     L = linear E contribution
#
#     C = quadratic E contribution
#
# so that
#
#     DQ = L + C.
# ------------------------------------------------------------------------

def determinant_components(g):

    Q00 = g[(0, 0)]["Q"]
    Q01 = g[(0, 1)]["Q"]
    Q10 = g[(1, 0)]["Q"]
    Q11 = g[(1, 1)]["Q"]

    E00 = g[(0, 0)]["E"]
    E01 = g[(0, 1)]["E"]
    E10 = g[(1, 0)]["E"]
    E11 = g[(1, 1)]["E"]

    L = (
        Q00 * E11
        + Q11 * E00
        - Q01 * E10
        - Q10 * E01
    )

    C = (
        E00 * E11
        - E01 * E10
    )

    DQ = (
        Q00 * Q11
        - Q01 * Q10
    )

    return DQ, L, C


# ------------------------------------------------------------------------
# FIRST QUESTION:
#
# How large is the "collapsed" determinant compared with Q^2?
# ------------------------------------------------------------------------

def scaling_report(n, g):

    vals = determinant_values(g)

    Qs = [
        g[k]["Q"]
        for k in g
    ]

    max_Q = max(Qs)

    DQ = vals["DQ"]

    return {
        "DQ": DQ,
        "abs_DQ": abs(DQ),
        "max_Q": max_Q,
        "Q2": max_Q * max_Q,
        "ratio": (
            abs(DQ)
            / (max_Q * max_Q)
        ),
    }


# ------------------------------------------------------------------------
# SECOND QUESTION:
#
# Do the carry values have a low-rank pattern of their own?
#
# Calculate:
#
#     det(E) = E00 E11 - E01 E10
#
# and compare it against the linear contribution.
# ------------------------------------------------------------------------

def carry_determinant(g):

    E00 = g[(0, 0)]["E"]
    E01 = g[(0, 1)]["E"]
    E10 = g[(1, 0)]["E"]
    E11 = g[(1, 1)]["E"]

    return (
        E00 * E11
        - E01 * E10
    )


# ------------------------------------------------------------------------
# THIRD QUESTION:
#
# Do neighboring radix cells produce predictable carry differences?
#
# ------------------------------------------------------------------------

def carry_differences(g):

    return {
        "row0": (
            g[(0, 1)]["E"]
            - g[(0, 0)]["E"]
        ),

        "row1": (
            g[(1, 1)]["E"]
            - g[(1, 0)]["E"]
        ),

        "col0": (
            g[(1, 0)]["E"]
            - g[(0, 0)]["E"]
        ),

        "col1": (
            g[(1, 1)]["E"]
            - g[(0, 1)]["E"]
        ),
    }


# ------------------------------------------------------------------------
# FOURTH QUESTION:
#
# Compare the determinant identity with a "random carry" control.
#
# We keep the actual Q values but replace the four E values by random
# values having approximately the same magnitude.
#
# This is NOT intended as a probability proof. It simply tells us whether
# the observed determinant cancellation is much more constrained than
# arbitrary numbers of the same scale.
# ------------------------------------------------------------------------

def random_carry_control(g, samples, rng):

    Es = [
        g[(0, 0)]["E"],
        g[(0, 1)]["E"],
        g[(1, 0)]["E"],
        g[(1, 1)]["E"],
    ]

    EMAX = max(Es)

    Q00 = g[(0, 0)]["Q"]
    Q01 = g[(0, 1)]["Q"]
    Q10 = g[(1, 0)]["Q"]
    Q11 = g[(1, 1)]["Q"]

    DQ = (
        Q00 * Q11
        - Q01 * Q10
    )

    hits = 0

    closest = None

    for _ in range(samples):

        e00 = rng.randint(0, EMAX)
        e01 = rng.randint(0, EMAX)
        e10 = rng.randint(0, EMAX)
        e11 = rng.randint(0, EMAX)

        correction = (
            Q00 * e11
            + Q11 * e00
            - Q01 * e10
            - Q10 * e01
            + e00 * e11
            - e01 * e10
        )

        error = abs(
            correction - DQ
        )

        if closest is None or error < closest:
            closest = error

        if error == 0:
            hits += 1

    return {
        "samples": samples,
        "exact_hits": hits,
        "closest_error": closest,
    }


# ------------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------------

random.seed(132)
rng = random.Random(132)

TEST_BITS = [
    30,
    36,
    42,
    48,
    54,
]


for bits in TEST_BITS:

    print()
    print("=" * 72)
    print("GENERATING", bits, "-BIT SEMIPRIME")
    print("=" * 72)

    n, p, q = random_semiprime(bits)

    print()
    print("-" * 72)

    print("n bits =", n.bit_length())
    print("n      =", n)
    print("true p =", p)
    print("true q =", q)

    # ------------------------------------------------------------
    # Build true grid.
    # ------------------------------------------------------------

    g = true_grid(
        n,
        p,
        q,
    )

    # ------------------------------------------------------------
    # Print Q and E.
    # ------------------------------------------------------------

    print()
    print("Q MATRIX")

    print(
        "[",
        g[(0, 0)]["Q"],
        g[(0, 1)]["Q"],
        "]"
    )

    print(
        "[",
        g[(1, 0)]["Q"],
        g[(1, 1)]["Q"],
        "]"
    )

    print()
    print("E MATRIX")

    print(
        "[",
        g[(0, 0)]["E"],
        g[(0, 1)]["E"],
        "]"
    )

    print(
        "[",
        g[(1, 0)]["E"],
        g[(1, 1)]["E"],
        "]"
    )

    print()
    print("K MATRIX")

    print(
        "[",
        g[(0, 0)]["K"],
        g[(0, 1)]["K"],
        "]"
    )

    print(
        "[",
        g[(1, 0)]["K"],
        g[(1, 1)]["K"],
        "]"
    )

    # ------------------------------------------------------------
    # Determinant identity.
    # ------------------------------------------------------------

    print()
    print("RANK-1 CHECK")

    d = determinant_values(g)

    print(
        "det(Q) =",
        d["DQ"]
    )

    print(
        "det(K) =",
        d["DK"]
    )

    print(
        "carry correction =",
        d["correction"]
    )

    print(
        "det(Q) == carry correction =",
        d["expanded_identity"]
    )

    # ------------------------------------------------------------
    # Components.
    # ------------------------------------------------------------

    DQ, L, C = determinant_components(g)

    print()
    print("DETERMINANT COMPONENTS")

    print(
        "D_Q =",
        DQ
    )

    print(
        "linear E contribution L =",
        L
    )

    print(
        "quadratic E contribution C =",
        C
    )

    print(
        "L + C == D_Q =",
        L + C == DQ
    )

    # ------------------------------------------------------------
    # Scale.
    # ------------------------------------------------------------

    scale = scaling_report(
        n,
        g
    )

    print()
    print("SCALING")

    print(
        "|D_Q| =",
        scale["abs_DQ"]
    )

    print(
        "max(Q)^2 =",
        scale["Q2"]
    )

    print(
        "|D_Q| / max(Q)^2 =",
        f"{scale['ratio']:.12e}"
    )

    # ------------------------------------------------------------
    # Carry determinant.
    # ------------------------------------------------------------

    print()
    print("CARRY DETERMINANT")

    E_det = carry_determinant(g)

    print(
        "det(E) =",
        E_det
    )

    print(
        "|det(E)| =",
        abs(E_det)
    )

    # ------------------------------------------------------------
    # Carry differences.
    # ------------------------------------------------------------

    print()
    print("CARRY DIFFERENCES")

    diffs = carry_differences(g)

    for name, value in diffs.items():

        print(
            name,
            "=",
            value
        )

    # ------------------------------------------------------------
    # Direct true-state substitution.
    #
    # This checks that the four E values really satisfy the only
    # determinant equation available at this level.
    # ------------------------------------------------------------

    print()
    print("TRUE STATE SUBSTITUTION")

    E00 = g[(0, 0)]["E"]
    E01 = g[(0, 1)]["E"]
    E10 = g[(1, 0)]["E"]
    E11 = g[(1, 1)]["E"]

    print(
        "E00 =",
        E00
    )

    print(
        "E01 =",
        E01
    )

    print(
        "E10 =",
        E10
    )

    print(
        "E11 =",
        E11
    )

    print(
        "det(K) = 0 =",
        (
            (g[(0, 0)]["Q"] - E00)
            * (g[(1, 1)]["Q"] - E11)
            -
            (g[(0, 1)]["Q"] - E01)
            * (g[(1, 0)]["Q"] - E10)
        ) == 0
    )

    # ------------------------------------------------------------
    # Random carry control.
    #
    # Keep sample count modest so this remains quick.
    # ------------------------------------------------------------

    print()
    print("RANDOM CARRY CONTROL")

    control = random_carry_control(
        g,
        samples=20000,
        rng=rng,
    )

    print(
        "random samples =",
        control["samples"]
    )

    print(
        "exact determinant matches =",
        control["exact_hits"]
    )

    print(
        "closest random correction error =",
        control["closest_error"]
    )

    # ------------------------------------------------------------
    # Structural conclusion for this instance.
    # ------------------------------------------------------------

    print()
    print("INSTANCE SUMMARY")

    print(
        "carry determinant magnitude =",
        abs(E_det)
    )

    print(
        "D_Q magnitude =",
        abs(DQ)
    )

    print(
        "linear contribution magnitude =",
        abs(L)
    )

    print(
        "quadratic contribution magnitude =",
        abs(C)
    )


# ========================================================================
# FINISHED EXPERIMENT 132
# ========================================================================

print()
print("=" * 72)
print("FINISHED EXPERIMENT 132")
print("=" * 72)
