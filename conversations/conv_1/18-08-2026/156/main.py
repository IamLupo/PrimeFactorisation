#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 156 — EXACT MINIMAL-COFACTOR 2-ADIC PROFILE / LIFT AUDIT
==============================================================================

Experiment 155 established the exact identity of supports:

    minimal 11x11 cofactors at v2 = 24

are exactly

    (5,3), (5,5), (5,7), (5,9), (5,10), (5,11),

and these are exactly the nonzero positions of

    R = C / Delta_11

modulo 2.

Therefore the mod-2 rank-one layer is literally the parity residue
of the minimum-valuation cofactor layer.

Experiment 156 studies that six-entry layer directly.

For every cofactor in the distinguished row i=5 we compute:

    C[5,j],
    v2(C[5,j]),
    C[5,j] / 2^24,
    (C[5,j] / 2^24) mod 2^e

for several powers 2^e.

We then test:

    1. which columns attain the minimum valuation;
    2. whether the normalized odd values agree modulo 4, 8, 16, 32;
    3. gcd/content of the six normalized integers;
    4. pairwise ratios modulo 2^e where inverses exist;
    5. exact equality patterns among the six residues;
    6. whether higher-valuation cofactors in the same row are
       congruent to zero at the corresponding precision;
    7. the complete valuation profile of row 5;
    8. comparison with the corresponding normalized R row.

The point is to determine whether the leading cofactor strip is:

    merely six unrelated odd integers,

or

    a coherent 2-adic vector with additional structure.

No recurrence is fitted.
No SymPy.
No floating point.
No extrapolation.
No connection to the original (p,q)-kernel is asserted.
"""

from __future__ import annotations

from math import gcd
from functools import reduce
from fractions import Fraction
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
N = 12
MIN_V2 = 24
DISTINGUISHED_ROW = 5


# ============================================================================
# BASIC HELPERS
# ============================================================================

def q_value(p: int, r: int) -> int:
    if r < 0 or r >= len(Q[p]):
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


def v_p(n: int, p: int) -> int | float:
    n = abs(int(n))

    if n == 0:
        return float("inf")

    v = 0

    while n % p == 0:
        n //= p
        v += 1

    return v


def gcd_list(values):
    values = [
        abs(int(x))
        for x in values
        if int(x) != 0
    ]

    if not values:
        return 0

    return reduce(gcd, values)


# ============================================================================
# EXACT BAREISS
# ============================================================================

def bareiss_det(A):

    if not A:
        return 1

    M = [[int(x) for x in row] for row in A]
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
            M[k], M[pivot_row] = (
                M[pivot_row],
                M[k],
            )
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
                            "Bareiss exact division failed."
                        )

                    value //= previous

                M[i][j] = value

        for i in range(k + 1, n):
            M[i][k] = 0

        previous = pivot

    return sign * M[-1][-1]


# ============================================================================
# BUILD SYSTEM
# ============================================================================

def build_system():

    M = []

    for p, pnext in TRANSITIONS:

        for r in range(D[p] + 1):

            d = D[p] - r

            q0 = q_value(p, r)
            q1 = q_value(p, r + 1)

            b = monomial_basis(p, d)

            row = [0] * 14

            for j in range(6):
                row[j] = b[j] * q0

            for j in range(6):
                row[6 + j] = b[j] * q1

            if d == 0:
                row[12] = 1
                row[13] = p

            M.append(row)

    return M


# ============================================================================
# BUILD F-CORE
# ============================================================================

def build_core(M):

    b0 = M[5]
    b1 = M[10]

    H = [
        [b0[12], b0[13]],
        [b1[12], b1[13]],
    ]

    det_h = (
        H[0][0] * H[1][1]
        - H[0][1] * H[1][0]
    )

    if det_h == 0:
        raise ArithmeticError(
            "Boundary block singular."
        )

    H_inv = [
        [
            Fraction(H[1][1], det_h),
            Fraction(-H[0][1], det_h),
        ],
        [
            Fraction(-H[1][0], det_h),
            Fraction(H[0][0], det_h),
        ],
    ]

    h_coeff = [
        [Fraction(0) for _ in range(12)]
        for _ in range(2)
    ]

    for i in range(2):

        for j in range(12):

            h_coeff[i][j] = -(
                H_inv[i][0] * b0[j]
                + H_inv[i][1] * b1[j]
            )

    core = []

    for idx in range(14):

        if idx in (5, 10):
            continue

        row = M[idx]

        out = []

        for j in range(12):

            value = (
                Fraction(row[j])
                + Fraction(row[12]) * h_coeff[0][j]
                + Fraction(row[13]) * h_coeff[1][j]
            )

            if value.denominator != 1:
                raise ArithmeticError(
                    "Non-integral F-core entry."
                )

            out.append(value.numerator)

        core.append(out)

    return core


# ============================================================================
# COFACTOR MATRIX
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


def cofactor_matrix(A):

    n = len(A)
    C = []

    for i in range(n):

        row = []

        for j in range(n):

            minor = bareiss_det(
                maximal_minor(A, i, j)
            )

            if (i + j) & 1:
                minor = -minor

            row.append(minor)

        C.append(row)

    return C


# ============================================================================
# NORMALIZED MATRIX
# ============================================================================

def normalized_matrix(C):

    R = []

    for row in C:

        out = []

        for x in row:

            if x % DELTA_11 != 0:
                raise ArithmeticError(
                    "Delta_11 does not divide cofactor."
                )

            out.append(
                x // DELTA_11
            )

        R.append(out)

    return R


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 156 — EXACT MINIMAL-COFACTOR "
        "2-ADIC PROFILE / LIFT AUDIT"
    )
    print("=" * 78)

    M = build_system()
    A = build_core(M)
    C = cofactor_matrix(A)
    R = normalized_matrix(C)

    # ------------------------------------------------------------------
    # 1. BASIC VALIDATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT DATA VALIDATION")
    print("=" * 78)

    row = C[DISTINGUISHED_ROW]

    print(
        f"  core_shape={len(A)}x{len(A[0])}"
    )

    print(
        f"  cofactor_shape={len(C)}x{len(C[0])}"
    )

    print(
        f"  distinguished_row={DISTINGUISHED_ROW}"
    )

    print(
        f"  Delta_11={DELTA_11}"
    )

    print(
        f"  v2(Delta_11)={v_p(DELTA_11, 2)}"
    )

    # ------------------------------------------------------------------
    # 2. FULL ROW VALUATION PROFILE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. DISTINGUISHED COFACTOR ROW VALUATION PROFILE")
    print("=" * 78)

    row_v2 = [
        v_p(x, 2)
        for x in row
    ]

    for j, (x, v) in enumerate(
        zip(row, row_v2)
    ):

        print(
            f"  j={j}: "
            f"v2={v} "
            f"value={x}"
        )

    # ------------------------------------------------------------------
    # 3. MINIMAL LAYER
    # ------------------------------------------------------------------

    minimal_cols = [
        j
        for j, v in enumerate(row_v2)
        if v == MIN_V2
    ]

    print()
    print("=" * 78)
    print("3. MINIMAL v2=24 LAYER")
    print("=" * 78)

    print(
        f"  minimal_columns={minimal_cols}"
    )

    print(
        f"  count={len(minimal_cols)}"
    )

    normalized_leading = {}

    for j in minimal_cols:

        value = row[j]

        if value % (2 ** MIN_V2) != 0:
            raise ArithmeticError(
                "Minimal cofactor not divisible by 2^24."
            )

        u = value // (2 ** MIN_V2)

        normalized_leading[j] = u

        print(
            f"  j={j}: "
            f"C/2^24={u}"
        )

    # ------------------------------------------------------------------
    # 4. ODD CONTENT
    # ------------------------------------------------------------------

    leading_values = list(
        normalized_leading.values()
    )

    leading_gcd = gcd_list(
        leading_values
    )

    print()
    print("=" * 78)
    print("4. LEADING-LAYER INTEGER CONTENT")
    print("=" * 78)

    print(
        f"  gcd_of_six_leading_values={leading_gcd}"
    )

    print(
        f"  gcd_v2={v_p(leading_gcd, 2)}"
    )

    print(
        f"  gcd_v3={v_p(leading_gcd, 3)}"
    )

    print(
        f"  gcd_v5={v_p(leading_gcd, 5)}"
    )

    # ------------------------------------------------------------------
    # 5. RESIDUES MOD 2^e
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. EXACT 2-ADIC LIFT OF MINIMAL COFACTORS")
    print("=" * 78)

    for e in (1, 2, 3, 4, 5, 6):

        modulus = 2 ** e

        residues = [
            normalized_leading[j] % modulus
            for j in minimal_cols
        ]

        distinct = len(set(residues))

        print(
            f"  modulus=2^{e}={modulus}: "
            f"residues={residues} "
            f"distinct={distinct}"
        )

    # ------------------------------------------------------------------
    # 6. PAIRWISE RATIO TEST
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. PAIRWISE 2-ADIC RATIO AUDIT")
    print("=" * 78)

    base_col = minimal_cols[0]
    base_value = normalized_leading[base_col]

    print(
        f"  base_column={base_col}"
    )

    print(
        f"  base_value={base_value}"
    )

    for j in minimal_cols[1:]:

        value = normalized_leading[j]

        ratios = []

        for e in (2, 3, 4, 5, 6):

            modulus = 2 ** e

            try:
                inv = pow(
                    base_value % modulus,
                    -1,
                    modulus,
                )

                ratio = (
                    value * inv
                ) % modulus

            except ValueError:

                ratio = None

            ratios.append(
                (modulus, ratio)
            )

        print(
            f"  column={j}: ratios={ratios}"
        )

    # ------------------------------------------------------------------
    # 7. NORMALIZED R ROW
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. NORMALIZED R ROW")
    print("=" * 78)

    R_row = R[DISTINGUISHED_ROW]

    for j, x in enumerate(R_row):

        print(
            f"  j={j}: "
            f"R={x} "
            f"parity={x & 1}"
        )

    # ------------------------------------------------------------------
    # 8. COMPARE 2^24-STRIPPED COFACTORS WITH R
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. STRIPPED COFACTOR VS NORMALIZED R")
    print("=" * 78)

    print(
        "  Compare:"
    )

    print(
        "      C[5,j] / 2^24"
    )

    print(
        "  with"
    )

    print(
        "      R[5,j] = C[5,j] / Delta_11."
    )

    print(
        "  Delta_11 / 2^24 = "
        f"{DELTA_11 // (2 ** MIN_V2)}"
    )

    scale = DELTA_11 // (2 ** MIN_V2)

    stripped_exact = True

    for j in minimal_cols:

        lhs = (
            row[j]
            // (2 ** MIN_V2)
        )

        rhs = R_row[j] * scale

        equal = (
            lhs == rhs
        )

        stripped_exact = (
            stripped_exact
            and equal
        )

        print(
            f"  j={j}: "
            f"stripped={lhs} "
            f"R_scaled={rhs} "
            f"exact={equal}"
        )

    # ------------------------------------------------------------------
    # 9. HIGHER-VALUATION COLUMNS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. HIGHER-VALUATION ROW-5 COLUMNS")
    print("=" * 78)

    for j in range(N):

        if j in minimal_cols:
            continue

        x = row[j]
        v = row_v2[j]

        excess = (
            v - MIN_V2
            if v != float("inf")
            else None
        )

        print(
            f"  j={j}: "
            f"v2={v} "
            f"excess={excess} "
            f"R_parity={R_row[j] & 1}"
        )

    # ------------------------------------------------------------------
    # 10. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 155 established that the mod-2 support of R is exactly
the set of cofactors attaining the minimum valuation v2=24.

Experiment 156 now examines the actual six odd numbers underneath
that support.

There are several possible outcomes:

    * all six values agree modulo larger powers of 2;
    * they differ immediately after parity;
    * they have a small common odd factor;
    * their pairwise ratios stabilize 2-adically;
    * higher-valuation columns form a predictable next layer.

The exact identity

    C[5,j] = 2^24 * u_j

for the six minimal columns isolates the leading arithmetic layer
without any determinant normalization.

The comparison with R is purely algebraic:

    R[5,j]
      = C[5,j] / Delta_11
      = u_j / (Delta_11 / 2^24).

Thus no new fitting is involved.

The purpose is to see whether the six-entry rank-one leading layer
contains any additional 2-adic regularity beyond its support.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(A) == 12
        and len(C) == 12
        and len(R) == 12
        and len(minimal_cols) == 6
        and all(
            row_v2[j] == MIN_V2
            for j in minimal_cols
        )
        and stripped_exact
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  core_exact={len(A) == 12}"
    )

    print(
        f"  cofactor_exact={len(C) == 12}"
    )

    print(
        f"  six_minimal_columns_exact="
        f"{len(minimal_cols) == 6}"
    )

    print(
        f"  v2_minimal_layer_exact="
        f"{all(row_v2[j] == MIN_V2 for j in minimal_cols)}"
    )

    print(
        f"  stripped_cofactor_identity_exact="
        f"{stripped_exact}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 156 COMPLETE")


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

