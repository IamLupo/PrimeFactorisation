#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 143 — EXACT p-ADIC LOW-RANK LIFT / MINOR-VALUATION AUDIT
==============================================================================

Experiment 142 found striking finite-field rank collapse for the normalized
cofactor matrix R:

    rank(R mod 2) = 1
    rank(R mod 3) = 1
    rank(R mod 7) = 2

This experiment does NOT attempt a Smith reduction.

Instead it asks a more targeted question:

    Does the low-rank collapse survive to higher p-adic precision?

For a matrix A over Z, a determinantal rank bound

    rank <= r  modulo p^e

is tested by checking whether every (r+1)x(r+1) minor is divisible by p^e.

Therefore:

    rank <= 1 mod p^e
        iff every 2x2 minor is divisible by p^e;

    rank <= 2 mod p^e
        iff every 3x3 minor is divisible by p^e.

We compute the exact minimum p-adic valuation among those minors.

This distinguishes:

    accidental rank collapse modulo p

from

    genuine p-adic low-rank structure.

The experiment also records the first valuation at which the low-rank
condition fails, and the locations of the least-divisible minors.

No floating point.
No SymPy.
No extrapolation.
No claim about the original (p,q)-kernel.
"""


from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import gcd
import sys


# ============================================================================
# EXACT DATA
# ============================================================================

Q = {
    1: [
        -126258696,
        -11600759760,
        2668721436,
        1764373740,
        -1338089411,
        495451247,
    ],
    3: [
        9955176,
        -1263551016,
        -152369292,
        -128667196,
        421514439,
    ],
    5: [
        -62398,
        4771718,
        16027881,
    ],
    7: [
        1,
    ],
}

D = {
    1: 5,
    3: 4,
    5: 2,
    7: 0,
}

TRANSITIONS = [
    (1, 3),
    (3, 5),
    (5, 7),
]

DELTA_11 = 9512681472


# ============================================================================
# BASIC HELPERS
# ============================================================================

def digits(n: int) -> int:
    return len(str(abs(n))) if n else 1


def valuation(n: int, p: int) -> int:
    n = abs(int(n))

    if n == 0:
        return 10**9

    v = 0

    while n % p == 0:
        n //= p
        v += 1

    return v


def gcd_many(values) -> int:
    g = 0

    for x in values:
        g = gcd(g, abs(int(x)))

    return abs(g)


def q_value(p: int, r: int) -> int:
    if r < 0:
        return 0

    if r >= len(Q[p]):
        return 0

    return Q[p][r]


def monomial_basis(p: int, d: int):
    return [
        1,
        d,
        d * d,
        p,
        p * d,
        p * p,
    ]


# ============================================================================
# EXACT BAREISS DETERMINANT
# ============================================================================

def bareiss_det(A):
    if not A:
        return 1

    M = [
        [int(x) for x in row]
        for row in A
    ]

    n = len(M)

    if n == 1:
        return M[0][0]

    previous = 1
    sign = 1

    for k in range(n - 1):

        pivot_row = None

        for r in range(k, n):
            if M[r][k] != 0:
                pivot_row = r
                break

        if pivot_row is None:
            return 0

        if pivot_row != k:
            M[k], M[pivot_row] = M[pivot_row], M[k]
            sign *= -1

        pivot = M[k][k]

        for i in range(k + 1, n):

            for j in range(k + 1, n):

                value = (
                    M[i][j] * pivot
                    - M[i][k] * M[k][j]
                )

                if k > 0:

                    if value % previous != 0:
                        raise ArithmeticError(
                            "Bareiss division failed."
                        )

                    value //= previous

                M[i][j] = value

        for i in range(k + 1, n):
            M[i][k] = 0

        previous = pivot

    return sign * M[-1][-1]


# ============================================================================
# BUILD FULL 14x14 SYSTEM
# ============================================================================

def build_full_system():

    M = []

    for p, p_next in TRANSITIONS:

        for r in range(D[p] + 1):

            d = D[p] - r

            q0 = q_value(p, r)
            q1 = q_value(p, r + 1)

            basis = monomial_basis(p, d)

            row = [0] * 14

            for j in range(6):
                row[j] = basis[j] * q0

            for j in range(6):
                row[6 + j] = basis[j] * q1

            if d == 0:
                row[12] = 1
                row[13] = p

            M.append(row)

    return M


# ============================================================================
# BUILD 12x12 F-CORE
# ============================================================================

def build_core(M):

    boundary0 = M[5]
    boundary1 = M[10]

    H = [
        [boundary0[12], boundary0[13]],
        [boundary1[12], boundary1[13]],
    ]

    H_det = (
        H[0][0] * H[1][1]
        - H[0][1] * H[1][0]
    )

    if H_det == 0:
        raise ArithmeticError(
            "Boundary block singular."
        )

    H_inv = [
        [
            Fraction(H[1][1], H_det),
            Fraction(-H[0][1], H_det),
        ],
        [
            Fraction(-H[1][0], H_det),
            Fraction(H[0][0], H_det),
        ],
    ]

    rhs0 = q_value(3, 5)
    rhs1 = q_value(5, 4)

    h_const = [
        H_inv[i][0] * rhs0
        + H_inv[i][1] * rhs1
        for i in range(2)
    ]

    h_coeff = [
        [Fraction(0) for _ in range(12)]
        for _ in range(2)
    ]

    for i in range(2):

        for j in range(12):

            h_coeff[i][j] = -(
                H_inv[i][0] * boundary0[j]
                + H_inv[i][1] * boundary1[j]
            )

    core = []

    for idx in range(14):

        if idx in (5, 10):
            continue

        row = M[idx]

        reduced = []

        for j in range(12):

            value = (
                Fraction(row[j])
                + Fraction(row[12]) * h_coeff[0][j]
                + Fraction(row[13]) * h_coeff[1][j]
            )

            if value.denominator != 1:
                raise ArithmeticError(
                    "Non-integral core entry."
                )

            reduced.append(value.numerator)

        core.append(reduced)

    return core, H_det


# ============================================================================
# NORMALIZED COFACTOR MATRIX
# ============================================================================

def maximal_minor(A, omit_row, omit_col):

    rows = [
        i for i in range(len(A))
        if i != omit_row
    ]

    cols = [
        j for j in range(len(A[0]))
        if j != omit_col
    ]

    return [
        [A[i][j] for j in cols]
        for i in rows
    ]


def normalized_cofactors(A):

    n = len(A)

    R = []

    for i in range(n):

        row = []

        for j in range(n):

            det_minor = bareiss_det(
                maximal_minor(A, i, j)
            )

            if det_minor % DELTA_11 != 0:
                raise ArithmeticError(
                    f"Delta_11 does not divide minor ({i},{j})."
                )

            row.append(
                det_minor // DELTA_11
            )

        R.append(row)

    return R


# ============================================================================
# SMALL-MINOR DETERMINANTS
# ============================================================================

def det_2x2(a, b, c, d):
    return a * d - b * c


def det_3x3(M):

    return (
        M[0][0] * (
            M[1][1] * M[2][2]
            - M[1][2] * M[2][1]
        )
        - M[0][1] * (
            M[1][0] * M[2][2]
            - M[1][2] * M[2][0]
        )
        + M[0][2] * (
            M[1][0] * M[2][1]
            - M[1][1] * M[2][0]
        )
    )


# ============================================================================
# 2x2 MINOR VALUATION AUDIT
# ============================================================================

def audit_2x2(A, p):

    n = len(A)

    minimum = 10**9
    positions = []
    zero_count = 0
    total = 0

    for rows in combinations(range(n), 2):

        i, j = rows

        for cols in combinations(range(n), 2):

            a, b = cols

            value = det_2x2(
                A[i][a],
                A[i][b],
                A[j][a],
                A[j][b],
            )

            total += 1

            if value == 0:
                zero_count += 1
                continue

            v = valuation(value, p)

            if v < minimum:
                minimum = v
                positions = [(i, j, a, b)]

            elif v == minimum:
                positions.append(
                    (i, j, a, b)
                )

    return {
        "minimum": minimum,
        "positions": positions,
        "zero_count": zero_count,
        "total": total,
    }


# ============================================================================
# 3x3 MINOR VALUATION AUDIT
# ============================================================================

def audit_3x3(A, p):

    n = len(A)

    minimum = 10**9
    positions = []
    zero_count = 0
    total = 0

    for rows in combinations(range(n), 3):

        for cols in combinations(range(n), 3):

            i, j, k = rows
            a, b, c = cols

            value = det_3x3([
                [
                    A[i][a],
                    A[i][b],
                    A[i][c],
                ],
                [
                    A[j][a],
                    A[j][b],
                    A[j][c],
                ],
                [
                    A[k][a],
                    A[k][b],
                    A[k][c],
                ],
            ])

            total += 1

            if value == 0:
                zero_count += 1
                continue

            v = valuation(value, p)

            if v < minimum:
                minimum = v
                positions = [
                    (i, j, k, a, b, c)
                ]

            elif v == minimum:
                positions.append(
                    (i, j, k, a, b, c)
                )

    return {
        "minimum": minimum,
        "positions": positions,
        "zero_count": zero_count,
        "total": total,
    }


# ============================================================================
# MODULAR RANK
# ============================================================================

def rank_mod(A, p):

    M = [
        [x % p for x in row]
        for row in A
    ]

    rows = len(M)
    cols = len(M[0])
    rank = 0

    for col in range(cols):

        pivot = None

        for r in range(rank, rows):

            if M[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            continue

        M[rank], M[pivot] = M[pivot], M[rank]

        inv = pow(
            M[rank][col],
            -1,
            p,
        )

        for j in range(col, cols):
            M[rank][j] = (
                M[rank][j] * inv
            ) % p

        for r in range(rows):

            if r == rank:
                continue

            factor = M[r][col]

            if factor == 0:
                continue

            for j in range(col, cols):
                M[r][j] = (
                    M[r][j]
                    - factor * M[rank][j]
                ) % p

        rank += 1

        if rank == rows:
            break

    return rank


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 143 — EXACT p-ADIC LOW-RANK LIFT / "
        "MINOR-VALUATION AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. BUILD
    # ------------------------------------------------------------------

    M = build_full_system()
    core, H_det = build_core(M)
    R = normalized_cofactors(core)

    print()
    print("=" * 78)
    print("1. EXACT DATA OBJECT")
    print("=" * 78)

    print("  core_shape=(12,12)")
    print("  normalized_cofactor_shape=(12,12)")
    print(f"  Delta_11={DELTA_11}")
    print(f"  H_det={H_det}")

    # ------------------------------------------------------------------
    # 2. MOD p BASELINE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. FINITE-FIELD BASELINE")
    print("=" * 78)

    baseline = {}

    for p in (2, 3, 7):

        rank = rank_mod(R, p)

        baseline[p] = rank

        print(
            f"  prime={p}: rank_mod_p={rank}"
        )

    # ------------------------------------------------------------------
    # 3. RANK-1 p-ADIC LIFT TEST
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. RANK-1 p-ADIC LIFT TEST")
    print("=" * 78)

    rank1_ok = {}

    for p in (2, 3):

        audit = audit_2x2(R, p)

        min_v = audit["minimum"]

        rank1_ok[p] = min_v >= 1

        print()
        print(
            f"  prime={p}"
        )

        print(
            f"    minimum_v2x2={min_v}"
        )

        print(
            f"    rank<=1_mod_p="
            f"{min_v >= 1}"
        )

        print(
            f"    rank<=1_mod_p2="
            f"{min_v >= 2}"
        )

        print(
            f"    rank<=1_mod_p3="
            f"{min_v >= 3}"
        )

        print(
            f"    total_2x2_minors={audit['total']}"
        )

        print(
            f"    zero_2x2_minors={audit['zero_count']}"
        )

        print(
            f"    minimum_positions_sample="
            f"{audit['positions'][:10]}"
        )

    # ------------------------------------------------------------------
    # 4. RANK-2 p-ADIC LIFT TEST
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. RANK-2 p-ADIC LIFT TEST")
    print("=" * 78)

    rank2_ok = {}

    # For p=7 the observed rank is 2, so this is the natural test.
    for p in (7,):

        audit = audit_3x3(R, p)

        min_v = audit["minimum"]

        rank2_ok[p] = min_v >= 1

        print()
        print(
            f"  prime={p}"
        )

        print(
            f"    minimum_v3x3={min_v}"
        )

        print(
            f"    rank<=2_mod_p="
            f"{min_v >= 1}"
        )

        print(
            f"    rank<=2_mod_p2="
            f"{min_v >= 2}"
        )

        print(
            f"    rank<=2_mod_p3="
            f"{min_v >= 3}"
        )

        print(
            f"    total_3x3_minors={audit['total']}"
        )

        print(
            f"    zero_3x3_minors={audit['zero_count']}"
        )

        print(
            f"    minimum_positions_sample="
            f"{audit['positions'][:10]}"
        )

    # ------------------------------------------------------------------
    # 5. FIRST FAILURE PRECISION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. FIRST p-ADIC FAILURE PRECISION")
    print("=" * 78)

    for p in (2, 3):

        audit = audit_2x2(R, p)
        min_v = audit["minimum"]

        first_failure = min_v + 1

        print(
            f"  prime={p}: "
            f"rank-one condition holds through e={min_v}; "
            f"fails at e={first_failure}"
        )

    audit7 = audit_3x3(R, 7)

    print(
        f"  prime=7: "
        f"rank-two condition holds through e={audit7['minimum']}; "
        f"fails at e={audit7['minimum'] + 1}"
    )

    # ------------------------------------------------------------------
    # 6. LOW-RANK CHARACTER
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. EXACT LOW-RANK CHARACTER")
    print("=" * 78)

    print(
        """
The rank test is interpreted through minors, not through ordinary
field rank over a non-field ring.

For prime p:

    all 2x2 minors divisible by p
        <=> determinantal rank <= 1 modulo p;

    all 3x3 minors divisible by p
        <=> determinantal rank <= 2 modulo p.

The minimum p-adic valuation of those minors tells exactly how far
that divisibility persists into higher p-adic precision.

This gives a clean test of whether the rank-one/rank-two phenomena
seen in Experiment 142 are merely modulo-p accidents or survive
to higher p-adic order.
"""
    )

    # ------------------------------------------------------------------
    # 7. FINAL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. FINAL EXACTNESS")
    print("=" * 78)

    baseline_ok = (
        baseline[2] == 1
        and baseline[3] == 1
        and baseline[7] == 2
    )

    # The experiment is successful provided the audits completed and
    # the baseline rank facts are reproduced.
    final_ok = (
        baseline_ok
        and all(
            p in rank1_ok
            for p in (2, 3)
        )
        and 7 in rank2_ok
    )

    print(
        f"  baseline_rank_exact={baseline_ok}"
    )

    print(
        f"  rank1_padic_audit_completed="
        f"{set(rank1_ok) == {2,3}}"
    )

    print(
        f"  rank2_padic_audit_completed="
        f"{7 in rank2_ok}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 143 COMPLETE")


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:
        print("\nInterrupted.")
        sys.exit(130)

    except Exception as exc:

        print(
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )

        raise

