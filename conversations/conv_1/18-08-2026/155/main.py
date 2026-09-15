#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 155 — EXACT MINIMAL-COFACTOR SUPPORT / MOD-2 DIRECTION AUDIT
==============================================================================

Experiment 154 established:

    min v2 of 1x1 minors  = 0
    min v2 of 2x2 minors  = 0
    min v2 of 3x3 minors  = 0

but

    min v2 of 11x11 minors = 24
    v2(det A)              = 33.

Thus the large 2-adic obstruction is concentrated at the cofactor
level rather than appearing gradually with minor size.

Experiment 155 studies the geometry of the 11x11 minors that attain
the minimum valuation 24.

For each cofactor position (i,j):

    C_ij = signed/minor cofactor,

we record

    v2(C_ij),
    C_ij / Delta_11 mod 2,

and whether that normalized cofactor is nonzero modulo 2.

The experiment then asks:

    1. Which exact cofactor positions attain v2 = 24?
    2. Do those positions coincide with the support of R mod 2?
    3. Is the support pattern determined by a small set of rows/columns?
    4. Do the six minimal cofactors generate exactly the rank-one
       direction of R mod 2?
    5. Are the minimal cofactors aligned with the first 2-adic defect
       E mod 2?

This directly connects the Smith/cofactor valuation layer to the
rank-one phenomenon already observed in R.

No recurrence is fitted.
No SymPy.
No floating point.
No extrapolation.
No connection to the original (p,q)-kernel is asserted.
"""

from __future__ import annotations

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


# ============================================================================
# HELPERS
# ============================================================================

def q_value(p, r):
    if r < 0 or r >= len(Q[p]):
        return 0
    return Q[p][r]


def monomial_basis(p, d):
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

    v = 0

    while n % p == 0:
        n //= p
        v += 1

    return v


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
# COFACTOR / NORMALIZED COFACTOR
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

            # Cofactor sign.
            if (i + j) & 1:
                minor = -minor

            row.append(minor)

        C.append(row)

    return C


# ============================================================================
# MOD-2 RANK
# ============================================================================

def rref_mod2(A):

    M = [
        [int(x) & 1 for x in row]
        for row in A
    ]

    rows = len(M)
    cols = len(M[0])

    pivot_row = 0
    pivots = []

    for col in range(cols):

        pivot = None

        for r in range(pivot_row, rows):

            if M[r][col]:
                pivot = r
                break

        if pivot is None:
            continue

        M[pivot_row], M[pivot] = (
            M[pivot],
            M[pivot_row],
        )

        for r in range(rows):

            if r == pivot_row:
                continue

            if M[r][col]:

                for j in range(col, cols):
                    M[r][j] ^= M[pivot_row][j]

        pivots.append(col)
        pivot_row += 1

        if pivot_row == rows:
            break

    return M, pivots


def rank_mod2(A):
    _, pivots = rref_mod2(A)
    return len(pivots)


# ============================================================================
# MOD-8 RANK-ONE / DEFECT
# ============================================================================

def build_R_and_defect(C):

    R = []

    for row in C:

        rrow = []

        for value in row:

            if value % DELTA_11 != 0:
                raise ArithmeticError(
                    "Delta_11 does not divide cofactor."
                )

            rrow.append(
                value // DELTA_11
            )

        R.append(rrow)

    # Exact mod-8 rank-one reconstruction.
    pivot = None
    inverse = None

    for i in range(N):

        for j in range(N):

            x = R[i][j] % 8

            if x == 0:
                continue

            try:
                inverse = pow(x, -1, 8)
            except ValueError:
                continue

            pivot = (i, j)
            break

        if pivot is not None:
            break

    if pivot is None:
        raise ArithmeticError(
            "No mod-8 unit pivot."
        )

    i0, j0 = pivot

    u = [
        R[i][j0] % 8
        for i in range(N)
    ]

    v = [
        R[i0][j] * inverse % 8
        for j in range(N)
    ]

    factor_exact = True

    for i in range(N):

        for j in range(N):

            if (
                u[i] * v[j] % 8
                != R[i][j] % 8
            ):

                factor_exact = False
                break

        if not factor_exact:
            break

    if not factor_exact:
        raise ArithmeticError(
            "Mod-8 rank-one reconstruction failed."
        )

    E = []

    for i in range(N):

        row = []

        for j in range(N):

            diff = (
                R[i][j]
                - u[i] * v[j]
            )

            if diff % 8 != 0:
                raise ArithmeticError(
                    "First defect is not integral."
                )

            row.append(
                (diff // 8) & 1
            )

        E.append(row)

    return R, E, pivot


# ============================================================================
# SUPPORT
# ============================================================================

def support_matrix(A):

    return [
        [
            (i, j)
            for j in range(len(A[0]))
            if A[i][j]
        ]
        for i in range(len(A))
    ]


def nonzero_positions(A):

    out = []

    for i in range(len(A)):

        for j in range(len(A[0])):

            if A[i][j]:
                out.append((i, j))

    return out


def row_support(A):

    return [
        i
        for i, row in enumerate(A)
        if any(row)
    ]


def col_support(A):

    return [
        j
        for j in range(len(A[0]))
        if any(A[i][j] for i in range(len(A)))
    ]


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 155 — EXACT MINIMAL-COFACTOR SUPPORT / "
        "MOD-2 DIRECTION AUDIT"
    )
    print("=" * 78)

    M = build_system()
    A = build_core(M)

    C = cofactor_matrix(A)

    # ------------------------------------------------------------------
    # 1. BASIC VALIDATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT CORE / COFACTOR VALIDATION")
    print("=" * 78)

    print(
        f"  core_shape={len(A)}x{len(A[0])}"
    )

    print(
        f"  cofactor_shape={len(C)}x{len(C[0])}"
    )

    det_A = bareiss_det(A)

    print(
        f"  det_v2={valuation(det_A, 2)}"
    )

    print(
        f"  Delta_11={DELTA_11}"
    )

    print(
        f"  Delta_11_v2={valuation(DELTA_11, 2)}"
    )

    # ------------------------------------------------------------------
    # 2. FULL COFACTOR VALUATION MATRIX
    # ------------------------------------------------------------------

    v2_matrix = [
        [
            valuation(C[i][j], 2)
            for j in range(N)
        ]
        for i in range(N)
    ]

    min_v2 = min(
        v2_matrix[i][j]
        for i in range(N)
        for j in range(N)
        if C[i][j] != 0
    )

    print()
    print("=" * 78)
    print("2. COFACTOR 2-ADIC VALUATION MATRIX")
    print("=" * 78)

    for row in v2_matrix:

        print(
            "  "
            + " ".join(
                str(x) if x != float("inf")
                else "inf"
                for x in row
            )
        )

    print(
        f"  minimum_v2={min_v2}"
    )

    # ------------------------------------------------------------------
    # 3. MINIMAL COFACTOR POSITIONS
    # ------------------------------------------------------------------

    minimal_positions = [
        (i, j)
        for i in range(N)
        for j in range(N)
        if C[i][j] != 0
        and valuation(C[i][j], 2) == min_v2
    ]

    print()
    print("=" * 78)
    print("3. MINIMAL 11x11 COFACTOR POSITIONS")
    print("=" * 78)

    print(
        f"  count={len(minimal_positions)}"
    )

    print(
        f"  positions={minimal_positions}"
    )

    # ------------------------------------------------------------------
    # 4. NORMALIZED COFACTOR SUPPORT MOD 2
    # ------------------------------------------------------------------

    R, E, defect_pivot = (
        build_R_and_defect(C)
    )

    R2 = [
        [
            R[i][j] & 1
            for j in range(N)
        ]
        for i in range(N)
    ]

    print()
    print("=" * 78)
    print("4. NORMALIZED COFACTOR MOD-2 SUPPORT")
    print("=" * 78)

    print(
        f"  rank_mod2={rank_mod2(R2)}"
    )

    print(
        f"  support_positions={nonzero_positions(R2)}"
    )

    print(
        f"  support_rows={row_support(R2)}"
    )

    print(
        f"  support_columns={col_support(R2)}"
    )

    # ------------------------------------------------------------------
    # 5. MINIMAL COFACTOR SUPPORT VS R MOD 2
    # ------------------------------------------------------------------

    minimal_mask = [
        [
            int((i, j) in set(minimal_positions))
            for j in range(N)
        ]
        for i in range(N)
    ]

    support_overlap = [
        (i, j)
        for i, j in minimal_positions
        if R2[i][j]
    ]

    print()
    print("=" * 78)
    print("5. MINIMAL-COFACTOR / R-MOD-2 OVERLAP")
    print("=" * 78)

    print(
        f"  minimal_positions={len(minimal_positions)}"
    )

    print(
        f"  positions_nonzero_after_normalization="
        f"{len(support_overlap)}"
    )

    print(
        f"  overlap={support_overlap}"
    )

    print(
        f"  all_minimal_cofactors_survive_mod2="
        f"{len(support_overlap) == len(minimal_positions)}"
    )

    # ------------------------------------------------------------------
    # 6. ROW/COLUMN PROFILE OF MINIMAL COFACTORS
    # ------------------------------------------------------------------

    minimal_rows = sorted(
        set(i for i, _ in minimal_positions)
    )

    minimal_cols = sorted(
        set(j for _, j in minimal_positions)
    )

    row_counts = {
        i: sum(
            1
            for r, _ in minimal_positions
            if r == i
        )
        for i in minimal_rows
    }

    col_counts = {
        j: sum(
            1
            for _, c in minimal_positions
            if c == j
        )
        for j in minimal_cols
    }

    print()
    print("=" * 78)
    print("6. MINIMAL-COFACTOR ROW/COLUMN PROFILE")
    print("=" * 78)

    print(
        f"  minimal_rows={minimal_rows}"
    )

    print(
        f"  minimal_columns={minimal_cols}"
    )

    print(
        f"  row_counts={row_counts}"
    )

    print(
        f"  column_counts={col_counts}"
    )

    # ------------------------------------------------------------------
    # 7. COMPARE WITH KNOWN R DIRECTION
    # ------------------------------------------------------------------

    rank1_positions = [
        (i, j)
        for i in range(N)
        for j in range(N)
        if R2[i][j]
    ]

    print()
    print("=" * 78)
    print("7. RANK-ONE DIRECTION VS MINIMAL COFACTORS")
    print("=" * 78)

    print(
        f"  R_mod2_support={rank1_positions}"
    )

    print(
        f"  minimal_cofactor_support="
        f"{minimal_positions}"
    )

    print(
        f"  support_difference_R_minus_minimal="
        f"{sorted(set(rank1_positions) - set(minimal_positions))}"
    )

    print(
        f"  support_difference_minimal_minus_R="
        f"{sorted(set(minimal_positions) - set(rank1_positions))}"
    )

    # ------------------------------------------------------------------
    # 8. FIRST DEFECT SUPPORT
    # ------------------------------------------------------------------

    E2 = [
        [
            E[i][j]
            for j in range(N)
        ]
        for i in range(N)
    ]

    print()
    print("=" * 78)
    print("8. FIRST 2-ADIC DEFECT SUPPORT")
    print("=" * 78)

    print(
        f"  defect_rank_mod2={rank_mod2(E2)}"
    )

    print(
        f"  defect_nonzero_positions="
        f"{nonzero_positions(E2)}"
    )

    print(
        f"  defect_rows={row_support(E2)}"
    )

    print(
        f"  defect_columns={col_support(E2)}"
    )

    # ------------------------------------------------------------------
    # 9. MINIMAL COFACTOR VS DEFECT PROJECTION
    # ------------------------------------------------------------------

    minimal_row_set = set(minimal_rows)
    minimal_col_set = set(minimal_cols)

    defect_on_minimal_rows = [
        (i, j)
        for i, j in nonzero_positions(E2)
        if i in minimal_row_set
    ]

    defect_on_minimal_cols = [
        (i, j)
        for i, j in nonzero_positions(E2)
        if j in minimal_col_set
    ]

    print()
    print("=" * 78)
    print("9. MINIMAL-COFACTOR / FIRST-DEFECT ALIGNMENT")
    print("=" * 78)

    print(
        f"  defect_positions_on_minimal_rows="
        f"{defect_on_minimal_rows}"
    )

    print(
        f"  defect_positions_on_minimal_columns="
        f"{defect_on_minimal_cols}"
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
Experiment 154 showed a sharp arithmetic transition:

    small minors:
        minimum v2 = 0;

    11x11 minors:
        minimum v2 = 24;

    determinant:
        v2 = 33.

Experiment 155 asks where the v2=24 layer actually lives.

The normalized matrix is

    R = C / Delta_11,

where C is the cofactor matrix.

If the v2=24 cofactors are exactly the entries that survive
modulo Delta_11 as odd integers, then the leading rank-one
matrix R mod 2 is literally the residue of the minimal cofactor
layer.

That would establish a precise arithmetic relationship between:

    minimal 11x11 minors
        and
    the rank-one normalized cofactor layer.

The first defect E is then the next 2-adic layer beyond that
minimal cofactor residue.

The experiment therefore separates:

    leading cofactor layer
        from
    first correction layer.

No claim is made that this is the source of the original kernel.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(A) == 12
        and len(C) == 12
        and min_v2 == 24
        and len(minimal_positions) == 6
        and rank_mod2(R2) == 1
        and rank_mod2(E2) == 3
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
        f"  minimum_cofactor_v2_match={min_v2 == 24}"
    )

    print(
        f"  six_minimizers_match="
        f"{len(minimal_positions) == 6}"
    )

    print(
        f"  R_rank_mod2_match={rank_mod2(R2) == 1}"
    )

    print(
        f"  defect_rank_mod2_match={rank_mod2(E2) == 3}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 155 COMPLETE")


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
