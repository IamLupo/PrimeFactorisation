#!/usr/bin/env python3

import math
import random


# ========================================================================
# START EXPERIMENT 138
# Exact cross-ratio cancellation of the rank-1 quotient matrix
# ========================================================================

print("=" * 72)
print("START EXPERIMENT 138")
print("Exact cross-ratio cancellation of the rank-1 quotient matrix")
print("=" * 72)


# ------------------------------------------------------------------------
# RADICES
# ------------------------------------------------------------------------

R1 = [17, 43, 59]
R2 = [19, 37, 61]

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
        "K": K,
        "E": E,
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

def build_grid(n, p, q):

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
# CROSS-RATIO EXPRESSION
#
# D =
#
#   Q10*Q01 - Q11*Q00
#
# Since:
#
#   K10*K01 = K11*K00
#
# the K*K terms disappear.
#
# ------------------------------------------------------------------------

def cross_expression(g):

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

    # Observable expression.
    DQ = (
        Q10 * Q01
        - Q11 * Q00
    )

    # The rank-1 quotient product cancels.
    DK = (
        K10 * K01
        - K11 * K00
    )

    # Expand the remaining terms explicitly.
    remainder = (
        K10 * E01
        + K01 * E10
        + E10 * E01
        - K11 * E00
        - K00 * E11
        - E11 * E00
    )

    return {
        "DQ": DQ,
        "DK": DK,
        "remainder": remainder,
        "identity_ok": (
            DK == 0
            and DQ == remainder
        ),
    }


# ------------------------------------------------------------------------
# FACTOR THE LINEAR K PART
#
# Linear K contribution:
#
#   K10 E01 + K01 E10 - K11 E00 - K00 E11
#
# Substitute:
#
#   K10 = k1*l0
#   K01 = k0*l1
#   K11 = k1*l1
#   K00 = k0*l0
#
# giving:
#
#   k1*l0*E01
# + k0*l1*E10
# - k1*l1*E00
# - k0*l0*E11.
#
# Group by k0/k1 or l0/l1.
# ------------------------------------------------------------------------

def factor_linear_K(g):

    k0 = g[(0, 0)]["k"]
    k1 = g[(1, 0)]["k"]

    l0 = g[(0, 0)]["ell"]
    l1 = g[(0, 1)]["ell"]

    E00 = g[(0, 0)]["E"]
    E01 = g[(0, 1)]["E"]
    E10 = g[(1, 0)]["E"]
    E11 = g[(1, 1)]["E"]

    direct = (
        k1 * l0 * E01
        + k0 * l1 * E10
        - k1 * l1 * E00
        - k0 * l0 * E11
    )

    # Group by k1 and k0:
    #
    # k1*(l0*E01-l1*E00)
    # +k0*(l1*E10-l0*E11)

    grouped_k = (
        k1 * (
            l0 * E01
            - l1 * E00
        )
        +
        k0 * (
            l1 * E10
            - l0 * E11
        )
    )

    # Group by l1 and l0:
    #
    # l1*(k0*E10-k1*E00)
    # +l0*(k1*E01-k0*E11)

    grouped_l = (
        l1 * (
            k0 * E10
            - k1 * E00
        )
        +
        l0 * (
            k1 * E01
            - k0 * E11
        )
    )

    return {
        "direct": direct,
        "grouped_k": grouped_k,
        "grouped_l": grouped_l,
        "k_equal": (
            direct == grouped_k
        ),
        "l_equal": (
            direct == grouped_l
        ),
    }


# ------------------------------------------------------------------------
# CARRY-ONLY QUADRATIC PART
# ------------------------------------------------------------------------

def carry_quadratic(g):

    E00 = g[(0, 0)]["E"]
    E01 = g[(0, 1)]["E"]
    E10 = g[(1, 0)]["E"]
    E11 = g[(1, 1)]["E"]

    return (
        E10 * E01
        - E11 * E00
    )


# ------------------------------------------------------------------------
# OBSERVABLE SCALE
# ------------------------------------------------------------------------

def scale_report(g):

    Q_values = [
        g[key]["Q"]
        for key in g
    ]

    E_values = [
        g[key]["E"]
        for key in g
    ]

    DQ = cross_expression(g)["DQ"]

    maxQ = max(
        abs(x)
        for x in Q_values
    )

    maxE = max(
        abs(x)
        for x in E_values
    )

    return {
        "DQ": DQ,
        "maxQ": maxQ,
        "maxE": maxE,
        "DQ_over_Q2": (
            abs(DQ)
            / (maxQ * maxQ)
        ),
    }


# ------------------------------------------------------------------------
# COMPARE WITH ORDINARY DETERMINANT
# ------------------------------------------------------------------------

def ordinary_determinant(g):

    Q00 = g[(0, 0)]["Q"]
    Q01 = g[(0, 1)]["Q"]
    Q10 = g[(1, 0)]["Q"]
    Q11 = g[(1, 1)]["Q"]

    return (
        Q00 * Q11
        - Q01 * Q10
    )


# ------------------------------------------------------------------------
# DIFFERENCE FORM
#
# Let:
#
#   dk = k1-k0
#   dl = l1-l0
#
# Because:
#
#   k0*r0 + a0 = k1*r1 + a1
#
# the differences are strongly connected to the residues.
# ------------------------------------------------------------------------

def quotient_differences(g):

    k0 = g[(0, 0)]["k"]
    k1 = g[(1, 0)]["k"]

    l0 = g[(0, 0)]["ell"]
    l1 = g[(0, 1)]["ell"]

    a0 = g[(0, 0)]["a"]
    a1 = g[(1, 0)]["a"]

    b0 = g[(0, 0)]["b"]
    b1 = g[(0, 1)]["b"]

    return {
        "dk": k1 - k0,
        "dl": l1 - l0,
        "da": a1 - a0,
        "db": b1 - b0,
        "k_relation": (
            R1[0] * k0
            + a0
            ==
            R1[1] * k1
            + a1
        ),
        "l_relation": (
            R2[0] * l0
            + b0
            ==
            R2[1] * l1
            + b1
        ),
    }


# ------------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------------

random.seed(138)

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

    g = build_grid(
        n,
        p,
        q,
    )

    # ------------------------------------------------------------
    # Display K.
    # ------------------------------------------------------------

    print()
    print("K MATRIX")

    for i in range(3):

        print(
            [
                g[(i, j)]["K"]
                for j in range(3)
            ]
        )

    print()
    print("E MATRIX")

    for i in range(3):

        print(
            [
                g[(i, j)]["E"]
                for j in range(3)
            ]
        )

    # ------------------------------------------------------------
    # Cross-expression.
    # ------------------------------------------------------------

    x = cross_expression(g)

    print()
    print("CROSS-RATIO CANCELLATION")

    print(
        "DQ =",
        x["DQ"]
    )

    print(
        "DK =",
        x["DK"]
    )

    print(
        "remainder =",
        x["remainder"]
    )

    print(
        "IDENTITY VERIFIED =",
        x["identity_ok"]
    )

    # ------------------------------------------------------------
    # Factor linear K contribution.
    # ------------------------------------------------------------

    f = factor_linear_K(g)

    print()
    print("LINEAR K FACTORIZATION")

    print(
        "direct =",
        f["direct"]
    )

    print(
        "grouped by k =",
        f["grouped_k"]
    )

    print(
        "grouped by l =",
        f["grouped_l"]
    )

    print(
        "k grouping exact =",
        f["k_equal"]
    )

    print(
        "l grouping exact =",
        f["l_equal"]
    )

    # ------------------------------------------------------------
    # Carry-only contribution.
    # ------------------------------------------------------------

    C = carry_quadratic(g)

    print()
    print(
        "CARRY-ONLY CROSS TERM =",
        C
    )

    # ------------------------------------------------------------
    # Quotient differences.
    # ------------------------------------------------------------

    d = quotient_differences(g)

    print()
    print("QUOTIENT / RESIDUE DIFFERENCES")

    print(
        "dk =",
        d["dk"]
    )

    print(
        "dl =",
        d["dl"]
    )

    print(
        "da =",
        d["da"]
    )

    print(
        "db =",
        d["db"]
    )

    print(
        "k radix relation =",
        d["k_relation"]
    )

    print(
        "ell radix relation =",
        d["l_relation"]
    )

    # ------------------------------------------------------------
    # Scale.
    # ------------------------------------------------------------

    scale = scale_report(g)

    print()
    print("SCALE")

    print(
        "|DQ| =",
        abs(scale["DQ"])
    )

    print(
        "max(Q) =",
        scale["maxQ"]
    )

    print(
        "max(E) =",
        scale["maxE"]
    )

    print(
        "|DQ| / max(Q)^2 =",
        f"{scale['DQ_over_Q2']:.12e}"
    )

    # ------------------------------------------------------------
    # Compare cross expression to ordinary determinant.
    # ------------------------------------------------------------

    OD = ordinary_determinant(g)

    print()
    print("COMPARISON")

    print(
        "ordinary det(Q) =",
        OD
    )

    print(
        "cross expression =",
        x["DQ"]
    )

    print(
        "ratio |cross| / |ordinary| =",
        (
            abs(x["DQ"]) / abs(OD)
            if OD != 0
            else float("inf")
        )
    )

    print()
    print("=" * 72)


# ========================================================================
# FINISHED EXPERIMENT 138
# ========================================================================

print()
print("=" * 72)
print("FINISHED EXPERIMENT 138")
print("=" * 72)
