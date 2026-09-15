#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 161 — EXACT MOD-2 ADJUGATE / TRUE LEFT-RIGHT KERNEL AUDIT
==============================================================================

Experiment 160 corrected the cofactor/adjugate orientation and proved:

    adj(A) A = det(A) I
    A adj(A) = det(A) I

and after normalization:

    AdjR A = d12 I
    A AdjR = d12 I

with

    v2(d12) = 9.

However, the previous kernel-space test contained one remaining
orientation mistake:

    nullspace_mod2(A)

is ker(A), the RIGHT nullspace.

The LEFT nullspace is

    ker(A^T).

Experiment 161 corrects this final distinction.

It tests:

    1. rank(A mod 2);
    2. rank(AdjR mod 2);
    3. true left nullspace ker(A^T);
    4. true right nullspace ker(A);
    5. every nonzero row of AdjR mod 2 lies in ker(A^T);
    6. every nonzero column of AdjR mod 2 lies in ker(A);
    7. the adjugate row-space as a subspace of ker(A^T);
    8. the adjugate column-space as a subspace of ker(A);
    9. exact equality/dimension comparisons;
   10. support vectors of AdjR mod 2.

This is a pure linear-algebra correction.

No recurrence fitting.
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
# BUILD SYSTEM
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
# COFACTOR / ADJUGATE
# ============================================================================

def maximal_minor(A, omit_row, omit_col):

    rows = [i for i in range(len(A)) if i != omit_row]
    cols = [j for j in range(len(A[0])) if j != omit_col]

    return [
        [A[i][j] for j in cols]
        for i in rows
    ]


def cofactor_matrix(A):

    C = []

    for i in range(N):

        row = []

        for j in range(N):

            value = bareiss_det(
                maximal_minor(A, i, j)
            )

            if (i + j) & 1:
                value = -value

            row.append(value)

        C.append(row)

    return C


def normalize(C):

    R = []

    for row in C:

        out = []

        for x in row:

            if x % DELTA_11 != 0:
                raise ArithmeticError(
                    "Delta_11 divisibility failure."
                )

            out.append(
                x // DELTA_11
            )

        R.append(out)

    return R


# ============================================================================
# MOD-2 LINEAR ALGEBRA
# ============================================================================

def transpose(A):

    return [
        [
            A[i][j]
            for i in range(len(A))
        ]
        for j in range(len(A[0]))
    ]


def rref_mod2(A):

    M = [
        [int(x) & 1 for x in row]
        for row in A
    ]

    rows = len(M)
    cols = len(M[0])

    pivots = []
    pivot_row = 0

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


def row_space_mod2(A):

    R, _ = rref_mod2(A)

    return [
        row
        for row in R
        if any(row)
    ]


def nullspace_mod2(A):

    R, pivots = rref_mod2(A)

    n = len(A[0])

    pivot_set = set(pivots)

    free_cols = [
        j for j in range(n)
        if j not in pivot_set
    ]

    basis = []

    for f in free_cols:

        v = [0] * n
        v[f] = 1

        for i, pc in enumerate(pivots):

            if R[i][f]:
                v[pc] = 1

        basis.append(v)

    return basis


def in_span(v, basis):

    if not basis:
        return not any(v)

    return (
        rank_mod2(basis)
        == rank_mod2(basis + [v])
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 161 — EXACT MOD-2 ADJUGATE / "
        "TRUE LEFT-RIGHT KERNEL AUDIT"
    )
    print("=" * 78)

    M = build_system()
    A = build_core(M)
    C = cofactor_matrix(A)
    R = normalize(C)

    AdjR = transpose(R)

    A2 = [
        [x & 1 for x in row]
        for row in A
    ]

    AdjR2 = [
        [x & 1 for x in row]
        for row in AdjR
    ]

    # ------------------------------------------------------------------
    # 1. RANKS
    # ------------------------------------------------------------------

    rank_A = rank_mod2(A2)
    rank_AdjR = rank_mod2(AdjR2)

    print()
    print("=" * 78)
    print("1. EXACT MOD-2 RANKS")
    print("=" * 78)

    print(
        f"  rank_A_mod2={rank_A}"
    )

    print(
        f"  rank_AdjR_mod2={rank_AdjR}"
    )

    # ------------------------------------------------------------------
    # 2. TRUE LEFT / RIGHT NULLSPACES
    # ------------------------------------------------------------------

    right_kernel = nullspace_mod2(A2)
    left_kernel = nullspace_mod2(
        transpose(A2)
    )

    print()
    print("=" * 78)
    print("2. TRUE LEFT / RIGHT NULLSPACES")
    print("=" * 78)

    print(
        f"  right_nullity={len(right_kernel)}"
    )

    print(
        f"  left_nullity={len(left_kernel)}"
    )

    for i, v in enumerate(left_kernel):

        print(
            f"  left_basis_{i}={v}"
        )

    # ------------------------------------------------------------------
    # 3. ADJUGATE SUPPORT
    # ------------------------------------------------------------------

    active_rows = [
        i
        for i, row in enumerate(AdjR2)
        if any(row)
    ]

    active_cols = [
        j
        for j in range(N)
        if any(
            AdjR2[i][j]
            for i in range(N)
        )
    ]

    print()
    print("=" * 78)
    print("3. ADJUGATE SUPPORT MOD 2")
    print("=" * 78)

    print(
        f"  active_rows={active_rows}"
    )

    print(
        f"  active_columns={active_cols}"
    )

    for i in active_rows:

        print(
            f"  row_{i}={AdjR2[i]}"
        )

    # ------------------------------------------------------------------
    # 4. ROWS VS TRUE LEFT NULLSPACE
    # ------------------------------------------------------------------

    nonzero_rows = [
        row
        for row in AdjR2
        if any(row)
    ]

    row_inclusion = [
        in_span(
            row,
            left_kernel,
        )
        for row in nonzero_rows
    ]

    print()
    print("=" * 78)
    print("4. ADJUGATE ROWS VS TRUE LEFT NULLSPACE")
    print("=" * 78)

    print(
        f"  row_inclusion_results={row_inclusion}"
    )

    print(
        f"  all_rows_in_left_nullspace="
        f"{all(row_inclusion)}"
    )

    # ------------------------------------------------------------------
    # 5. COLUMNS VS TRUE RIGHT NULLSPACE
    # ------------------------------------------------------------------

    AdjRT = transpose(AdjR2)

    nonzero_cols = [
        [
            AdjR2[i][j]
            for i in range(N)
        ]
        for j in active_cols
    ]

    col_inclusion = [
        in_span(
            col,
            right_kernel,
        )
        for col in nonzero_cols
    ]

    print()
    print("=" * 78)
    print("5. ADJUGATE COLUMNS VS TRUE RIGHT NULLSPACE")
    print("=" * 78)

    print(
        f"  column_inclusion_results={col_inclusion}"
    )

    print(
        f"  all_columns_in_right_nullspace="
        f"{all(col_inclusion)}"
    )

    # ------------------------------------------------------------------
    # 6. SPACE DIMENSIONS
    # ------------------------------------------------------------------

    row_space = row_space_mod2(AdjR2)
    col_space = row_space_mod2(
        transpose(AdjR2)
    )

    print()
    print("=" * 78)
    print("6. SPACE DIMENSIONS")
    print("=" * 78)

    print(
        f"  adjugate_rowspace_rank={len(row_space)}"
    )

    print(
        f"  adjugate_columnspace_rank={len(col_space)}"
    )

    print(
        f"  true_left_nullity={len(left_kernel)}"
    )

    print(
        f"  true_right_nullity={len(right_kernel)}"
    )

    # ------------------------------------------------------------------
    # 7. DIRECT MATRIX TESTS
    # ------------------------------------------------------------------

    def matvec_mod2(A, v):

        out = []

        for row in A:

            s = 0

            for a, b in zip(row, v):

                s ^= (
                    (a & 1)
                    & (b & 1)
                )

            out.append(s)

        return out

    direct_left = [
        matvec_mod2(
            transpose(A2),
            row,
        )
        for row in nonzero_rows
    ]

    direct_right = [
        matvec_mod2(
            A2,
            col,
        )
        for col in nonzero_cols
    ]

    print()
    print("=" * 78)
    print("7. DIRECT KERNEL TEST")
    print("=" * 78)

    print(
        f"  Adj_row_annihilation_results={direct_left}"
    )

    print(
        f"  Adj_column_annihilation_results={direct_right}"
    )

    # ------------------------------------------------------------------
    # 8. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The exact adjugate identities imply:

    Adj(A) A = det(A) I,
    A Adj(A) = det(A) I.

Modulo 2, since 2 divides det(A), this gives:

    every row of Adj(A) mod 2
        lies in ker(A^T);

    every column of Adj(A) mod 2
        lies in ker(A).

The previous experiment used the wrong nullspace orientation when
testing the first statement.

Experiment 161 performs both tests against their actual nullspaces.

There are now three distinct objects:

    A mod 2:
        rank 3;

    Adj(A)/Delta_11 mod 2:
        rank 1;

    true left/right nullspaces of A:
        dimension 9 on each side.

The interesting structural question is therefore not whether the
adjugate is "large" or whether its coefficients are complicated.

It is whether the rank-one adjugate residue selects a distinguished
one-dimensional subspace of the nine-dimensional left/right kernels.

A positive answer is an exact linear-algebraic statement.
"""
    )

    # ------------------------------------------------------------------
    # 9. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        rank_A == 3
        and rank_AdjR == 1
        and len(left_kernel) == 9
        and len(right_kernel) == 9
        and all(row_inclusion)
        and all(col_inclusion)
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  rank_A_exact={rank_A == 3}"
    )

    print(
        f"  rank_adjugate_mod2_exact={rank_AdjR == 1}"
    )

    print(
        f"  true_left_nullity_exact={len(left_kernel) == 9}"
    )

    print(
        f"  true_right_nullity_exact={len(right_kernel) == 9}"
    )

    print(
        f"  adjugate_rows_in_true_left_kernel="
        f"{all(row_inclusion)}"
    )

    print(
        f"  adjugate_columns_in_true_right_kernel="
        f"{all(col_inclusion)}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 161 COMPLETE")


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

