#!/usr/bin/env python3

import math
import random


# ========================================================================
# START EXPERIMENT 135
# 3x3 rank-1 grid: independent minor constraints
# ========================================================================

print("=" * 72)
print("START EXPERIMENT 135")
print("3x3 rank-1 grid: independent minor constraints")
print("=" * 72)


# ------------------------------------------------------------------------
# 3x3 RADIX GRID
#
# Keep the radices pairwise coprime within each direction.
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
# MATRIX HELPERS
# ------------------------------------------------------------------------

def matrix(g, key):

    return [
        [
            g[(i, j)][key]
            for j in range(3)
        ]
        for i in range(3)
    ]


def print_matrix(name, M):

    print()
    print(name)

    for row in M:
        print(row)


# ------------------------------------------------------------------------
# 2x2 MINOR
# ------------------------------------------------------------------------

def minor2(M, i, j, m, n):

    return (
        M[i][j] * M[m][n]
        - M[i][n] * M[m][j]
    )


# ------------------------------------------------------------------------
# ALL 2x2 MINORS
#
# There are 9 possible row-pair / column-pair minors in a 3x3 matrix.
# For rank 1, every one must vanish.
# ------------------------------------------------------------------------

def all_minors(M):

    out = []

    for i in range(3):

        for m in range(i + 1, 3):

            for j in range(3):

                for n in range(j + 1, 3):

                    value = (
                        M[i][j] * M[m][n]
                        - M[i][n] * M[m][j]
                    )

                    out.append(
                        {
                            "rows": (i, m),
                            "cols": (j, n),
                            "value": value,
                        }
                    )

    return out


# ------------------------------------------------------------------------
# Q MINOR CORRECTION
#
# Since:
#
#       K = Q-E
#
# and every 2x2 K-minor is zero,
#
#       det(Q-E) = 0
#
# for each 2x2 submatrix.
#
# We calculate:
#
#       DQ
#
# and the exact carry correction separately.
# ------------------------------------------------------------------------

def minor_correction(Q, E, i, m, j, n):

    q00 = Q[i][j]
    q01 = Q[i][n]
    q10 = Q[m][j]
    q11 = Q[m][n]

    e00 = E[i][j]
    e01 = E[i][n]
    e10 = E[m][j]
    e11 = E[m][n]

    DQ = (
        q00 * q11
        - q01 * q10
    )

    L = (
        q00 * e11
        + q11 * e00
        - q01 * e10
        - q10 * e01
    )

    DE = (
        e00 * e11
        - e01 * e10
    )

    # Correct identity:
    #
    # det(Q) = L - det(E)
    #
    reconstructed = L - DE

    return {
        "DQ": DQ,
        "L": L,
        "DE": DE,
        "reconstructed": reconstructed,
        "ok": (
            reconstructed == DQ
        ),
    }


# ------------------------------------------------------------------------
# NORMALIZED MINOR SIZE
# ------------------------------------------------------------------------

def minor_scale_report(Q, E):

    Q_minors = all_minors(Q)
    E_minors = all_minors(E)

    max_q_minor = max(
        abs(x["value"])
        for x in Q_minors
    )

    max_e_minor = max(
        abs(x["value"])
        for x in E_minors
    )

    return {
        "max_Q_minor": max_q_minor,
        "max_E_minor": max_e_minor,
    }


# ------------------------------------------------------------------------
# INDEPENDENT MINOR COUNT
#
# The nine minors are algebraically dependent.
#
# For a rank-1 3x3 matrix, the constraints are equivalent to the matrix
# having rank <= 1.
#
# We inspect how many minors are numerically zero for K and how large
# the corresponding observable Q minors are.
# ------------------------------------------------------------------------

def rank_of_integer_matrix(M):

    """
    Tiny exact rank routine for 3x3 integer matrices using Gaussian
    elimination over fractions represented as rational arithmetic.
    """

    from fractions import Fraction

    A = [
        [
            Fraction(x)
            for x in row
        ]
        for row in M
    ]

    rows = 3
    cols = 3

    rank = 0

    for col in range(cols):

        pivot = None

        for row in range(rank, rows):

            if A[row][col] != 0:

                pivot = row
                break

        if pivot is None:
            continue

        A[rank], A[pivot] = (
            A[pivot],
            A[rank],
        )

        pivot_value = A[rank][col]

        for j in range(col, cols):

            A[rank][j] /= pivot_value

        for row in range(rows):

            if row == rank:
                continue

            factor = A[row][col]

            if factor == 0:
                continue

            for j in range(col, cols):

                A[row][j] -= (
                    factor
                    * A[rank][j]
                )

        rank += 1

    return rank


# ------------------------------------------------------------------------
# 3x3 MINOR CORRECTION SUMMARY
# ------------------------------------------------------------------------

def correction_summary(Q, E):

    summary = []

    for i in range(3):

        for m in range(i + 1, 3):

            for j in range(3):

                for n in range(j + 1, 3):

                    info = minor_correction(
                        Q,
                        E,
                        i,
                        m,
                        j,
                        n,
                    )

                    summary.append(
                        {
                            "rows": (i, m),
                            "cols": (j, n),
                            **info,
                        }
                    )

    return summary


# ------------------------------------------------------------------------
# CARRY DIFFERENCES
#
# Select the (0,0) value as a reference and express every E entry
# as:
#
#       Eij = E00 + Dij
#
# This gives 8 difference variables.
#
# The experiment measures whether the 3x3 system gives more algebraic
# constraints than the 2x2 case.
# ------------------------------------------------------------------------

def carry_difference_matrix(E):

    E00 = E[0][0]

    D = [
        [
            E[i][j] - E00
            for j in range(3)
        ]
        for i in range(3)
    ]

    return D


# ------------------------------------------------------------------------
# MIXED DIFFERENCES
# ------------------------------------------------------------------------

def mixed_difference_matrix(E):

    out = []

    for i in range(2):

        row = []

        for j in range(2):

            value = (
                E[i + 1][j + 1]
                - E[i + 1][j]
                - E[i][j + 1]
                + E[i][j]
            )

            row.append(value)

        out.append(row)

    return out


# ------------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------------

random.seed(135)

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

    Q = matrix(g, "Q")
    E = matrix(g, "E")
    K = matrix(g, "K")

    # ------------------------------------------------------------
    # Print matrices.
    # ------------------------------------------------------------

    print_matrix(
        "Q MATRIX",
        Q,
    )

    print_matrix(
        "E MATRIX",
        E,
    )

    print_matrix(
        "K MATRIX",
        K,
    )

    # ------------------------------------------------------------
    # Rank.
    # ------------------------------------------------------------

    print()
    print("EXACT MATRIX RANK")

    print(
        "rank(Q) =",
        rank_of_integer_matrix(Q)
    )

    print(
        "rank(E) =",
        rank_of_integer_matrix(E)
    )

    print(
        "rank(K) =",
        rank_of_integer_matrix(K)
    )

    # ------------------------------------------------------------
    # K minors.
    # ------------------------------------------------------------

    K_minors = all_minors(K)

    zero_K_minors = sum(
        1
        for x in K_minors
        if x["value"] == 0
    )

    print()
    print("K MINORS")

    print(
        "total 2x2 minors =",
        len(K_minors)
    )

    print(
        "zero K minors =",
        zero_K_minors
    )

    # ------------------------------------------------------------
    # Q minors.
    # ------------------------------------------------------------

    Q_minors = all_minors(Q)

    print()
    print("Q MINORS")

    for item in Q_minors:

        print(
            "rows =",
            item["rows"],
            "cols =",
            item["cols"],
            "value =",
            item["value"],
        )

    # ------------------------------------------------------------
    # Carry correction for every minor.
    # ------------------------------------------------------------

    summary = correction_summary(
        Q,
        E,
    )

    print()
    print("MINOR CORRECTION IDENTITIES")

    all_ok = True

    for item in summary:

        if not item["ok"]:
            all_ok = False

        print(
            "rows =",
            item["rows"],
            "cols =",
            item["cols"],
            "DQ =",
            item["DQ"],
            "L =",
            item["L"],
            "DE =",
            item["DE"],
            "L-DE =",
            item["reconstructed"],
            "OK =",
            item["ok"],
        )

    print()
    print(
        "ALL 3x3 MINOR IDENTITIES VERIFIED =",
        all_ok
    )

    # ------------------------------------------------------------
    # Carry difference matrix.
    # ------------------------------------------------------------

    D = carry_difference_matrix(E)

    print_matrix(
        "E DIFFERENCE MATRIX (relative to E00)",
        D,
    )

    # ------------------------------------------------------------
    # Mixed differences.
    # ------------------------------------------------------------

    MD = mixed_difference_matrix(E)

    print_matrix(
        "E MIXED DIFFERENCE MATRIX",
        MD,
    )

    # ------------------------------------------------------------
    # Scale.
    # ------------------------------------------------------------

    scale = minor_scale_report(
        Q,
        E,
    )

    print()
    print("MINOR SCALE")

    print(
        "max |2x2 minor(Q)| =",
        scale["max_Q_minor"]
    )

    print(
        "max |2x2 minor(E)| =",
        scale["max_E_minor"]
    )

    # ------------------------------------------------------------
    # Compare E differences to E itself.
    # ------------------------------------------------------------

    E00 = E[0][0]

    max_difference = max(
        abs(E[i][j] - E00)
        for i in range(3)
        for j in range(3)
    )

    max_E = max(
        abs(E[i][j])
        for i in range(3)
        for j in range(3)
    )

    max_mixed = max(
        abs(x)
        for row in MD
        for x in row
    )

    print()
    print("DIFFERENCE SCALING")

    print(
        "max |Eij-E00| =",
        max_difference
    )

    print(
        "max |Eij| =",
        max_E
    )

    print(
        "max |mixed E difference| =",
        max_mixed
    )

    if max_E:

        print(
            "difference / E scale =",
            f"{max_difference / max_E:.12e}"
        )

    # ------------------------------------------------------------
    # Number of independent-looking minor magnitudes.
    # ------------------------------------------------------------

    distinct_abs_Q_minors = sorted(
        set(
            abs(x["value"])
            for x in Q_minors
        )
    )

    distinct_abs_E_minors = sorted(
        set(
            abs(x["value"])
            for x in all_minors(E)
        )
    )

    print()
    print("MINOR DIVERSITY")

    print(
        "distinct |Q minors| =",
        len(distinct_abs_Q_minors)
    )

    print(
        "distinct |E minors| =",
        len(distinct_abs_E_minors)
    )

    # ------------------------------------------------------------
    # Show whether the 3x3 grid gives additional exact rank
    # constraints beyond one 2x2 determinant.
    # ------------------------------------------------------------

    print()
    print("STRUCTURAL RESULT")

    if zero_K_minors == 9:

        print(
            "ALL 2x2 MINORS OF K VANISH"
        )

    if rank_of_integer_matrix(K) == 1:

        print(
            "K IS EXACTLY RANK 1"
        )

    print(
        "NUMBER OF CARRY VARIABLES =",
        9
    )

    print(
        "NUMBER OF E-DIFFERENCES RELATIVE TO E00 =",
        8
    )

    print(
        "NUMBER OF 2x2 MINOR EQUATIONS =",
        9
    )

    print()
    print("=" * 72)


# ========================================================================
# FINISHED EXPERIMENT 135
# ========================================================================

print()
print("=" * 72)
print("FINISHED EXPERIMENT 135")
print("=" * 72)
