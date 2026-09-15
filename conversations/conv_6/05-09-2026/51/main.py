#!/usr/bin/env python3

import math
import random


# ========================================================================
# START EXPERIMENT 134
# Eliminate the absolute carry level from the 2x2 determinant
# ========================================================================

print("=" * 72)
print("START EXPERIMENT 134")
print("Eliminate the absolute carry level from the 2x2 determinant")
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
# DIFFERENCES
# ------------------------------------------------------------------------

def carry_differences(g):

    E00 = g[(0, 0)]["E"]
    E01 = g[(0, 1)]["E"]
    E10 = g[(1, 0)]["E"]
    E11 = g[(1, 1)]["E"]

    return {
        "E00": E00,

        "d0": E01 - E00,

        "v0": E10 - E00,

        "v1": E11 - E01,

        "mixed": (
            E11
            - E10
            - E01
            + E00
        ),
    }


# ------------------------------------------------------------------------
# DIRECT DETERMINANT
# ------------------------------------------------------------------------

def detQ(g):

    Q00 = g[(0, 0)]["Q"]
    Q01 = g[(0, 1)]["Q"]
    Q10 = g[(1, 0)]["Q"]
    Q11 = g[(1, 1)]["Q"]

    return (
        Q00 * Q11
        - Q01 * Q10
    )


# ------------------------------------------------------------------------
# ELIMINATED DETERMINANT EQUATION
#
# Let
#
#   E01 = X + d0
#   E10 = X + v0
#   E11 = X + d0 + v1
#
# where X = E00.
#
# Then det(Q-E)=0 can be written:
#
#   DQ = X * A + B
#
# for coefficients A,B depending only on Q and the differences.
#
# We derive these explicitly and solve for X.
# ------------------------------------------------------------------------

def eliminated_equation(g):

    Q00 = g[(0, 0)]["Q"]
    Q01 = g[(0, 1)]["Q"]
    Q10 = g[(1, 0)]["Q"]
    Q11 = g[(1, 1)]["Q"]

    diffs = carry_differences(g)

    X = diffs["E00"]
    d0 = diffs["d0"]
    v0 = diffs["v0"]
    v1 = diffs["v1"]

    # ------------------------------------------------------------
    # E values represented using X.
    # ------------------------------------------------------------

    E00 = X
    E01 = X + d0
    E10 = X + v0
    E11 = X + d0 + v1

    # ------------------------------------------------------------
    # Compute linear coefficient symbolically:
    #
    # L =
    #
    # Q00*E11
    # + Q11*E00
    # - Q01*E10
    # - Q10*E01
    #
    # ------------------------------------------------------------

    A = (
        Q00
        + Q11
        - Q01
        - Q10
    )

    B = (
        Q00 * (d0 + v1)
        - Q01 * v0
        - Q10 * d0
    )

    # ------------------------------------------------------------
    # Correct carry determinant:
    #
    # det(E)
    # =
    # X*(v1-v0) - d0*v0
    #
    # Therefore:
    #
    # DQ = L - det(E)
    #
    #     = X*A + B
    #       - X*(v1-v0)
    #       + d0*v0
    #
    #     = X*COEFF + CONSTANT
    # ------------------------------------------------------------

    coeff = (
        A
        - v1
        + v0
    )

    constant = (
        B
        + d0 * v0
    )

    DQ = detQ(g)

    # Equation:
    #
    #     DQ = X*coeff + constant
    #
    return {
        "DQ": DQ,
        "A": A,
        "B": B,
        "coeff": coeff,
        "constant": constant,
        "X_true": X,
        "reconstructed": (
            X * coeff
            + constant
        ),
    }


# ------------------------------------------------------------------------
# SOLVE FOR E00
# ------------------------------------------------------------------------

def solve_E00(g):

    e = eliminated_equation(g)

    coeff = e["coeff"]
    rhs = e["DQ"] - e["constant"]

    if coeff == 0:

        return {
            **e,
            "solvable": False,
            "E00_recovered": None,
        }

    if rhs % coeff != 0:

        return {
            **e,
            "solvable": False,
            "E00_recovered": None,
        }

    X = rhs // coeff

    return {
        **e,
        "solvable": True,
        "E00_recovered": X,
    }


# ------------------------------------------------------------------------
# RECONSTRUCT COMPLETE E MATRIX
# ------------------------------------------------------------------------

def reconstruct_E(g, E00):

    diffs = carry_differences(g)

    d0 = diffs["d0"]
    v0 = diffs["v0"]
    v1 = diffs["v1"]

    return {
        "E00": E00,
        "E01": E00 + d0,
        "E10": E00 + v0,
        "E11": E00 + d0 + v1,
    }


# ------------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------------

random.seed(134)

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

    g = grid(
        n,
        p,
        q,
    )

    # ------------------------------------------------------------
    # Print Q.
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

    # ------------------------------------------------------------
    # Carry differences.
    # ------------------------------------------------------------

    d = carry_differences(g)

    print()
    print("CARRY DIFFERENCE VARIABLES")

    for name in (
        "E00",
        "d0",
        "v0",
        "v1",
        "mixed",
    ):

        print(
            name,
            "=",
            d[name]
        )

    # ------------------------------------------------------------
    # Direct determinant.
    # ------------------------------------------------------------

    DQ = detQ(g)

    print()
    print("DIRECT determinant(Q) =", DQ)

    # ------------------------------------------------------------
    # Eliminate E00.
    # ------------------------------------------------------------

    result = solve_E00(g)

    print()
    print("ELIMINATED EQUATION")

    print(
        "coefficient =",
        result["coeff"]
    )

    print(
        "constant =",
        result["constant"]
    )

    print(
        "true E00 =",
        result["X_true"]
    )

    print(
        "reconstructed DQ =",
        result["reconstructed"]
    )

    print(
        "equation verified =",
        result["reconstructed"] == DQ
    )

    print()
    print("SOLVING FOR E00")

    print(
        "solvable =",
        result["solvable"]
    )

    if result["solvable"]:

        E00_recovered = result[
            "E00_recovered"
        ]

        print(
            "recovered E00 =",
            E00_recovered
        )

        print(
            "true E00 =",
            result["X_true"]
        )

        print(
            "E00 EXACT =",
            E00_recovered
            == result["X_true"]
        )

        # --------------------------------------------------------
        # Reconstruct all E values.
        # --------------------------------------------------------

        recovered_E = reconstruct_E(
            g,
            E00_recovered,
        )

        true_E = {
            "E00": g[(0, 0)]["E"],
            "E01": g[(0, 1)]["E"],
            "E10": g[(1, 0)]["E"],
            "E11": g[(1, 1)]["E"],
        }

        print()
        print("RECONSTRUCTED E MATRIX")

        print(
            "[",
            recovered_E["E00"],
            recovered_E["E01"],
            "]"
        )

        print(
            "[",
            recovered_E["E10"],
            recovered_E["E11"],
            "]"
        )

        print()
        print(
            "FULL E MATRIX EXACT =",
            recovered_E == true_E
        )

    else:

        print(
            "E00 could not be recovered exactly."
        )

    # ------------------------------------------------------------
    # Compare the size of the unknowns.
    # ------------------------------------------------------------

    print()
    print("STATE SIZE")

    print(
        "|E00| =",
        abs(d["E00"])
    )

    print(
        "|d0| =",
        abs(d["d0"])
    )

    print(
        "|v0| =",
        abs(d["v0"])
    )

    print(
        "|v1| =",
        abs(d["v1"])
    )

    print(
        "|mixed| =",
        abs(d["mixed"])
    )

    print()
    print("=" * 72)


# ========================================================================
# FINISHED EXPERIMENT 134
# ========================================================================

print()
print("=" * 72)
print("FINISHED EXPERIMENT 134")
print("=" * 72)
