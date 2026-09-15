#!/usr/bin/env python3

import math
import random


# ========================================================================
# START EXPERIMENT 133
# Corrected determinant identity and carry-difference factorization
# ========================================================================

print("=" * 72)
print("START EXPERIMENT 133")
print("Corrected determinant identity and carry-difference factorization")
print("=" * 72)


# ------------------------------------------------------------------------
# RADICES
# ------------------------------------------------------------------------

R1 = [17, 43]
R2 = [19, 37]

print()
print("R1 =", R1)
print("R2 =", R2)


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
# PRIME
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
# CELL
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
# GRID
# ------------------------------------------------------------------------

def grid(n, p, q):

    out = {}

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            out[(i, j)] = cell(
                n,
                p,
                q,
                r,
                s,
            )

    return out


# ------------------------------------------------------------------------
# CORRECT DETERMINANT IDENTITY
# ------------------------------------------------------------------------

def determinant_identity(g):

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

    detQ = (
        Q00 * Q11
        - Q01 * Q10
    )

    detE = (
        E00 * E11
        - E01 * E10
    )

    L = (
        Q00 * E11
        + Q11 * E00
        - Q01 * E10
        - Q10 * E01
    )

    detK = (
        K00 * K11
        - K01 * K10
    )

    reconstructed = L - detE

    return {
        "detQ": detQ,
        "detE": detE,
        "L": L,
        "detK": detK,
        "reconstructed": reconstructed,
        "identity_ok": (
            detQ == reconstructed
            and detK == 0
        ),
    }


# ------------------------------------------------------------------------
# CARRY DIFFERENCE VARIABLES
# ------------------------------------------------------------------------

def difference_variables(g):

    E00 = g[(0, 0)]["E"]
    E01 = g[(0, 1)]["E"]
    E10 = g[(1, 0)]["E"]
    E11 = g[(1, 1)]["E"]

    Q00 = g[(0, 0)]["Q"]
    Q01 = g[(0, 1)]["Q"]
    Q10 = g[(1, 0)]["Q"]
    Q11 = g[(1, 1)]["Q"]

    # ------------------------------------------------------------
    # Horizontal carry differences.
    # ------------------------------------------------------------

    d0 = E01 - E00
    d1 = E11 - E10

    # ------------------------------------------------------------
    # Vertical carry differences.
    # ------------------------------------------------------------

    v0 = E10 - E00
    v1 = E11 - E01

    # ------------------------------------------------------------
    # Mixed difference.
    # ------------------------------------------------------------

    mixed = (
        E11
        - E10
        - E01
        + E00
    )

    # ------------------------------------------------------------
    # Same quantities on Q.
    # ------------------------------------------------------------

    qh0 = Q01 - Q00
    qh1 = Q11 - Q10

    qv0 = Q10 - Q00
    qv1 = Q11 - Q01

    qmixed = (
        Q11
        - Q10
        - Q01
        + Q00
    )

    return {
        "d0": d0,
        "d1": d1,
        "v0": v0,
        "v1": v1,
        "mixed": mixed,
        "qh0": qh0,
        "qh1": qh1,
        "qv0": qv0,
        "qv1": qv1,
        "qmixed": qmixed,
    }


# ------------------------------------------------------------------------
# EXPRESS DET(E) USING DIFFERENCES
#
# Put
#
#     E01 = E00 + d0
#     E10 = E00 + v0
#     E11 = E00 + d0 + v1
#
# Then:
#
#     det(E) = E00*E11 - E01*E10
#
# We can investigate whether the common E00 dependence cancels in
# useful combinations.
# ------------------------------------------------------------------------

def difference_det_expansion(g):

    E00 = g[(0, 0)]["E"]

    diffs = difference_variables(g)

    d0 = diffs["d0"]
    v0 = diffs["v0"]
    v1 = diffs["v1"]

    E01 = E00 + d0
    E10 = E00 + v0
    E11 = E00 + d0 + v1

    direct = (
        E00 * E11
        - E01 * E10
    )

    expanded = (
        E00 * (d0 + v1 - v0)
        - d0 * v0
    )

    return {
        "direct": direct,
        "expanded": expanded,
        "equal": direct == expanded,
    }


# ------------------------------------------------------------------------
# TEST WHETHER DETERMINANT DEPENDS PRIMARILY ON DIFFERENCES
# ------------------------------------------------------------------------

def normalized_measures(g):

    Q00 = g[(0, 0)]["Q"]
    Q01 = g[(0, 1)]["Q"]
    Q10 = g[(1, 0)]["Q"]
    Q11 = g[(1, 1)]["Q"]

    E00 = g[(0, 0)]["E"]
    E01 = g[(0, 1)]["E"]
    E10 = g[(1, 0)]["E"]
    E11 = g[(1, 1)]["E"]

    detQ = (
        Q00 * Q11
        - Q01 * Q10
    )

    maxQ = max(
        abs(Q00),
        abs(Q01),
        abs(Q10),
        abs(Q11),
    )

    maxE = max(
        abs(E00),
        abs(E01),
        abs(E10),
        abs(E11),
    )

    diffs = difference_variables(g)

    max_diff = max(
        abs(diffs["d0"]),
        abs(diffs["d1"]),
        abs(diffs["v0"]),
        abs(diffs["v1"]),
        abs(diffs["mixed"]),
    )

    return {
        "detQ": detQ,
        "Q2": maxQ * maxQ,
        "Emax": maxE,
        "diffmax": max_diff,
        "detQ_over_Q2": (
            abs(detQ) / (maxQ * maxQ)
        ),
        "diff_over_E": (
            max_diff / maxE
            if maxE
            else 0.0
        ),
    }


# ------------------------------------------------------------------------
# RANDOM CARRY MATRIX CONTROL
#
# Preserve the scale of E but randomize it.
#
# This tests whether the exact determinant identity is special compared
# with arbitrary carry-sized integers.
# ------------------------------------------------------------------------

def random_control(g, rng, samples=10000):

    Q00 = g[(0, 0)]["Q"]
    Q01 = g[(0, 1)]["Q"]
    Q10 = g[(1, 0)]["Q"]
    Q11 = g[(1, 1)]["Q"]

    Emax = max(
        g[(0, 0)]["E"],
        g[(0, 1)]["E"],
        g[(1, 0)]["E"],
        g[(1, 1)]["E"],
    )

    detQ = (
        Q00 * Q11
        - Q01 * Q10
    )

    exact = 0
    closest = None

    for _ in range(samples):

        e00 = rng.randint(0, Emax)
        e01 = rng.randint(0, Emax)
        e10 = rng.randint(0, Emax)
        e11 = rng.randint(0, Emax)

        L = (
            Q00 * e11
            + Q11 * e00
            - Q01 * e10
            - Q10 * e01
        )

        detE = (
            e00 * e11
            - e01 * e10
        )

        error = abs(
            (L - detE)
            - detQ
        )

        if error == 0:
            exact += 1

        if closest is None or error < closest:
            closest = error

    return exact, closest


# ------------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------------

random.seed(133)
rng = random.Random(133)

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

    g = grid(n, p, q)

    # ------------------------------------------------------------
    # Print matrices.
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

    # ------------------------------------------------------------
    # Correct determinant identity.
    # ------------------------------------------------------------

    d = determinant_identity(g)

    print()
    print("CORRECTED DETERMINANT IDENTITY")

    print(
        "det(Q) =",
        d["detQ"]
    )

    print(
        "L =",
        d["L"]
    )

    print(
        "det(E) =",
        d["detE"]
    )

    print(
        "L - det(E) =",
        d["reconstructed"]
    )

    print(
        "det(K) =",
        d["detK"]
    )

    print(
        "IDENTITY VERIFIED =",
        d["identity_ok"]
    )

    # ------------------------------------------------------------
    # Difference variables.
    # ------------------------------------------------------------

    diffs = difference_variables(g)

    print()
    print("CARRY DIFFERENCE VARIABLES")

    for name in (
        "d0",
        "d1",
        "v0",
        "v1",
        "mixed",
    ):

        print(
            name,
            "=",
            diffs[name]
        )

    print()
    print("Q DIFFERENCE VARIABLES")

    for name in (
        "qh0",
        "qh1",
        "qv0",
        "qv1",
        "qmixed",
    ):

        print(
            name,
            "=",
            diffs[name]
        )

    # ------------------------------------------------------------
    # Difference determinant expansion.
    # ------------------------------------------------------------

    de = difference_det_expansion(g)

    print()
    print("det(E) DIFFERENCE EXPANSION")

    print(
        "direct =",
        de["direct"]
    )

    print(
        "expanded =",
        de["expanded"]
    )

    print(
        "equal =",
        de["equal"]
    )

    # ------------------------------------------------------------
    # Normalized measurements.
    # ------------------------------------------------------------

    m = normalized_measures(g)

    print()
    print("NORMALIZED STRUCTURE")

    print(
        "|det(Q)| =",
        abs(m["detQ"])
    )

    print(
        "max(Q)^2 =",
        m["Q2"]
    )

    print(
        "|det(Q)| / max(Q)^2 =",
        f"{m['detQ_over_Q2']:.12e}"
    )

    print(
        "max(E) =",
        m["Emax"]
    )

    print(
        "max carry difference =",
        m["diffmax"]
    )

    print(
        "max difference / max(E) =",
        f"{m['diff_over_E']:.12e}"
    )

    # ------------------------------------------------------------
    # Random control.
    # ------------------------------------------------------------

    print()
    print("RANDOM CARRY CONTROL")

    exact, closest = random_control(
        g,
        rng,
        samples=10000,
    )

    print(
        "samples =",
        10000
    )

    print(
        "exact determinant matches =",
        exact
    )

    print(
        "closest error =",
        closest
    )

    print()
    print("=" * 72)


# ========================================================================
# FINISHED EXPERIMENT 133
# ========================================================================

print()
print("=" * 72)
print("FINISHED EXPERIMENT 133")
print("=" * 72)
