#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 157 — EXACT MINIMAL-COFACTOR ODD-CONTENT / 3-7-ADIC AUDIT
==============================================================================

Experiment 156 established:

    C[5,j] = 2^24 * u_j

for

    j in {3,5,7,9,10,11},

with

    gcd(u_j) = 1701 = 3^5 * 7,

while

    Delta_11 / 2^24 = 567 = 3^4 * 7.

Hence the corresponding normalized cofactor entries have the exact
common factor

    gcd(R[5,j]) = 3.

Experiment 157 isolates this odd arithmetic content.

The experiment computes:

    1. the exact six-entry minimal cofactor vector u;
    2. gcd(u);
    3. gcd of the normalized R entries;
    4. quotient identities

           gcd(u) / (Delta_11 / 2^24),

    5. 3-adic valuations of u and R;
    6. 7-adic valuations of u and R;
    7. residues after dividing out the common factor 3;
    8. 3^e residue profiles;
    9. 7^e residue profiles;
   10. pairwise ratios in the 3-adic and 7-adic unit groups.

The goal is narrow:

    Is the common factor 3 merely an accidental gcd,
    or does it correspond to a coherent odd-prime layer left over
    after the 2^24 normalization and the Delta_11 normalization?

No polynomial fitting.
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
# DATA
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

ROW = 5
MIN_V2 = 24
MIN_COLS = [3, 5, 7, 9, 10, 11]


# ============================================================================
# HELPERS
# ============================================================================

def q_value(p, r):
    if r < 0 or r >= len(Q[p]):
        return 0
    return Q[p][r]


def basis(p, d):
    return [
        1,
        d,
        d * d,
        p,
        p * d,
        p * p,
    ]


def valuation(n, p):
    n = abs(int(n))

    if n == 0:
        return float("inf")

    out = 0

    while n % p == 0:
        n //= p
        out += 1

    return out


def gcd_list(values):
    values = [abs(int(x)) for x in values if x != 0]

    if not values:
        return 0

    return reduce(gcd, values)


# ============================================================================
# BAREISS
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
                            "Bareiss exact division failed."
                        )

                    value //= previous

                M[i][j] = value

        for i in range(k + 1, n):
            M[i][k] = 0

        previous = pivot

    return sign * M[-1][-1]


# ============================================================================
# SYSTEM
# ============================================================================

def build_system():

    M = []

    for p, pnext in TRANSITIONS:

        for r in range(D[p] + 1):

            d = D[p] - r

            q0 = q_value(p, r)
            q1 = q_value(p, r + 1)

            b = basis(p, d)

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
# F-CORE
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
        raise ArithmeticError("Boundary block singular.")

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
# COFACTORS
# ============================================================================

def maximal_minor(A, omit_row, omit_col):

    rows = [i for i in range(len(A)) if i != omit_row]
    cols = [j for j in range(len(A[0])) if j != omit_col]

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

            value = bareiss_det(
                maximal_minor(A, i, j)
            )

            if (i + j) & 1:
                value = -value

            row.append(value)

        C.append(row)

    return C


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 157 — EXACT MINIMAL-COFACTOR "
        "ODD-CONTENT / 3-7-ADIC AUDIT"
    )
    print("=" * 78)

    M = build_system()
    A = build_core(M)
    C = cofactor_matrix(A)

    row = C[ROW]

    u = {
        j: row[j] // (2 ** MIN_V2)
        for j in MIN_COLS
    }

    R_row = {
        j: row[j] // DELTA_11
        for j in MIN_COLS
    }

    # ------------------------------------------------------------------
    # 1. VALIDATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT VALIDATION")
    print("=" * 78)

    print(
        f"  distinguished_row={ROW}"
    )

    print(
        f"  minimal_columns={MIN_COLS}"
    )

    print(
        f"  v2_scale={MIN_V2}"
    )

    print(
        f"  Delta_11={DELTA_11}"
    )

    delta_odd = DELTA_11 // (2 ** MIN_V2)

    print(
        f"  Delta_11 / 2^24={delta_odd}"
    )

    # ------------------------------------------------------------------
    # 2. GCD CONTENT
    # ------------------------------------------------------------------

    gcd_u = gcd_list(u.values())
    gcd_R = gcd_list(R_row.values())

    print()
    print("=" * 78)
    print("2. ODD CONTENT")
    print("=" * 78)

    print(
        f"  gcd_u={gcd_u}"
    )

    print(
        f"  gcd_R={gcd_R}"
    )

    print(
        f"  gcd_u / Delta_odd="
        f"{Fraction(gcd_u, delta_odd)}"
    )

    print(
        f"  gcd_u_factorization="
        f"3^{valuation(gcd_u,3)} * "
        f"7^{valuation(gcd_u,7)}"
    )

    print(
        f"  gcd_R_v3={valuation(gcd_R,3)}"
    )

    print(
        f"  gcd_R_v7={valuation(gcd_R,7)}"
    )

    # ------------------------------------------------------------------
    # 3. ENTRYWISE 3/7 VALUATIONS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. ENTRYWISE 3-ADIC / 7-ADIC VALUATIONS")
    print("=" * 78)

    for j in MIN_COLS:

        print(
            f"  j={j}: "
            f"v3(u)={valuation(u[j],3)} "
            f"v7(u)={valuation(u[j],7)} "
            f"v3(R)={valuation(R_row[j],3)} "
            f"v7(R)={valuation(R_row[j],7)}"
        )

    # ------------------------------------------------------------------
    # 4. DIVIDE OUT COMMON FACTOR 3
    # ------------------------------------------------------------------

    if gcd_R != 0 and gcd_R % 3 != 0:
        reduced_factor = 1
    else:
        reduced_factor = 3

    reduced = {
        j: R_row[j] // reduced_factor
        for j in MIN_COLS
    }

    print()
    print("=" * 78)
    print("4. COMMON-FACTOR-REMOVED NORMALIZED LAYER")
    print("=" * 78)

    print(
        f"  factor_removed={reduced_factor}"
    )

    for j in MIN_COLS:

        print(
            f"  j={j}: "
            f"R/3={reduced[j]}"
        )

    print(
        f"  gcd_after_dividing_3="
        f"{gcd_list(reduced.values())}"
    )

    # ------------------------------------------------------------------
    # 5. 3-ADIC RESIDUES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. 3-ADIC RESIDUE PROFILE")
    print("=" * 78)

    for e in (1, 2, 3, 4, 5):

        modulus = 3 ** e

        residues = [
            reduced[j] % modulus
            for j in MIN_COLS
        ]

        print(
            f"  modulus=3^{e}={modulus}: "
            f"residues={residues} "
            f"distinct={len(set(residues))}"
        )

    # ------------------------------------------------------------------
    # 6. 7-ADIC RESIDUES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. 7-ADIC RESIDUE PROFILE")
    print("=" * 78)

    for e in (1, 2, 3):

        modulus = 7 ** e

        residues = [
            reduced[j] % modulus
            for j in MIN_COLS
        ]

        print(
            f"  modulus=7^{e}={modulus}: "
            f"residues={residues} "
            f"distinct={len(set(residues))}"
        )

    # ------------------------------------------------------------------
    # 7. PAIRWISE UNIT RATIOS AT 3
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. PAIRWISE 3-ADIC UNIT-RATIO AUDIT")
    print("=" * 78)

    base_col = MIN_COLS[0]
    base = reduced[base_col]

    print(
        f"  base_column={base_col}"
    )

    for j in MIN_COLS[1:]:

        ratios = []

        for e in (1, 2, 3, 4):

            modulus = 3 ** e

            try:
                inv = pow(
                    base % modulus,
                    -1,
                    modulus,
                )

                ratio = (
                    reduced[j]
                    * inv
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
    # 8. PAIRWISE UNIT RATIOS AT 7
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. PAIRWISE 7-ADIC UNIT-RATIO AUDIT")
    print("=" * 78)

    base = reduced[base_col]

    for j in MIN_COLS[1:]:

        ratios = []

        for e in (1, 2, 3):

            modulus = 7 ** e

            try:
                inv = pow(
                    base % modulus,
                    -1,
                    modulus,
                )

                ratio = (
                    reduced[j]
                    * inv
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
    # 9. STRUCTURAL RELATION TO DELTA_11
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. DELTA_11 ODD-CONTENT RELATION")
    print("=" * 78)

    print(
        f"  Delta_11=2^24 * {delta_odd}"
    )

    print(
        f"  gcd_u={gcd_u}"
    )

    print(
        f"  gcd_u = Delta_odd * {Fraction(gcd_u, delta_odd)}"
    )

    print(
        f"  gcd_R={gcd_R}"
    )

    print(
        f"  exact_identity="
        f"{gcd_u == delta_odd * gcd_R}"
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
Experiment 156 found

    gcd(u_3,u_5,u_7,u_9,u_10,u_11) = 1701,

while

    Delta_11 / 2^24 = 567.

Therefore the corresponding normalized cofactor entries have

    gcd(R_5,*) = 3.

Experiment 157 isolates this odd content.

The main exact question is whether the residual factor 3 is itself
part of a coherent odd-prime layer.

The experiment therefore separates:

    2-adic normalization:
        C = 2^24 * u;

    determinant normalization:
        R = u / 567;

    residual common content:
        gcd(R_minimal entries) = 3.

If division by 3 exposes a particularly simple 3-adic or 7-adic
pattern, that is evidence for an additional arithmetic layer.

If the residues immediately become irregular, then the factor 3
should be regarded simply as common content rather than as a new
structural coordinate.

No claim is made from finite residue coincidence alone.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(A) == 12
        and len(C) == 12
        and gcd_u == 1701
        and gcd_R == 3
        and delta_odd == 567
        and gcd_u == delta_odd * gcd_R
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  core_exact={len(A) == 12}"
    )

    print(
        f"  gcd_u_exact={gcd_u == 1701}"
    )

    print(
        f"  Delta_odd_exact={delta_odd == 567}"
    )

    print(
        f"  gcd_R_exact={gcd_R == 3}"
    )

    print(
        f"  content_identity_exact="
        f"{gcd_u == delta_odd * gcd_R}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 157 COMPLETE")


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

