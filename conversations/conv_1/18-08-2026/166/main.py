#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 166 — EXACT MOD-2 SCHUR QUOTIENT / ADJUGATE PAIRING AUDIT
==============================================================================

Established:

    rank_F2(F) = 4
    rank_F2(A) = 3
    rank_F2(K) = 1

with

    A = F + K.

Moreover,

    Row(F) = Row(A) ⊕ Row(K),
    Col(F) = Col(A) ⊕ Col(K)

modulo 2.

The normalized adjugate of A has rank one modulo 2, but its row and
column directions are not contained in the image spaces of F or A.

Experiment 166 asks whether the one-dimensional Schur quotient and the
one-dimensional adjugate kernel directions are canonically paired.

For a square matrix over a field, row-space annihilators and quotient
directions are naturally dual. Here we explicitly construct:

    Q_row  = Row(F) / Row(A),
    Q_col  = Col(F) / Col(A),

and compare them with:

    ker(A^T),
    ker(A).

The main test is whether the Schur correction direction is detected
nontrivially by the adjugate kernel direction through the standard
dot product over F_2.

This is an exact finite-dimensional pairing test.

No recurrence fitting.
No SymPy.
No floating point.
No extrapolation.
No claim about the original (p,q)-kernel.
"""

from __future__ import annotations

from fractions import Fraction
import sys


# ============================================================================
# DATA
# ============================================================================

QDATA = {
    1: [-126258696, -11600759760, 2668721436,
        1764373740, -1338089411, 495451247],
    3: [9955176, -1263551016, -152369292,
        -128667196, 421514439],
    5: [-62398, 4771718, 16027881],
    7: [1],
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

FULL_N = 14
CORE_N = 12

BOUNDARY_ROWS = (5, 10)
BOUNDARY_COLS = (12, 13)

DELTA_11 = 9512681472


# ============================================================================
# BUILDING
# ============================================================================

def q_value(p, r):
    if r < 0 or r >= len(QDATA[p]):
        return 0
    return QDATA[p][r]


def basis(p, d):
    return [1, d, d*d, p, p*d, p*p]


def bareiss_det(A):
    if not A:
        return 1

    M = [[int(x) for x in row] for row in A]
    n = len(M)

    if n == 1:
        return M[0][0]

    prev = 1
    sign = 1

    for k in range(n - 1):

        pivot_row = None

        for r in range(k, n):
            if M[r][k]:
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

                x = (
                    M[i][j] * pivot
                    - M[i][k] * M[k][j]
                )

                if k:
                    if x % prev:
                        raise ArithmeticError(
                            "Bareiss exact division failed."
                        )
                    x //= prev

                M[i][j] = x

        for i in range(k + 1, n):
            M[i][k] = 0

        prev = pivot

    return sign * M[-1][-1]


def build_full_system():
    M = []

    for p, _ in TRANSITIONS:

        for r in range(D[p] + 1):

            d = D[p] - r
            q0 = q_value(p, r)
            q1 = q_value(p, r + 1)

            b = basis(p, d)

            row = [0] * FULL_N

            for j in range(6):
                row[j] = b[j] * q0

            for j in range(6):
                row[6+j] = b[j] * q1

            if d == 0:
                row[12] = 1
                row[13] = p

            M.append(row)

    return M


def build_core(M):

    retained_rows = [
        i for i in range(FULL_N)
        if i not in BOUNDARY_ROWS
    ]

    boundary_rows = [
        M[i] for i in BOUNDARY_ROWS
    ]

    F = [
        [M[i][j] for j in range(CORE_N)]
        for i in retained_rows
    ]

    B = [
        [M[i][j] for j in BOUNDARY_COLS]
        for i in retained_rows
    ]

    H = [
        [row[j] for j in BOUNDARY_COLS]
        for row in boundary_rows
    ]

    C = [
        [row[j] for j in range(CORE_N)]
        for row in boundary_rows
    ]

    detH = H[0][0]*H[1][1] - H[0][1]*H[1][0]

    if detH == 0:
        raise ArithmeticError("det(H)=0")

    Hinv = [
        [Fraction(H[1][1], detH),
         Fraction(-H[0][1], detH)],
        [Fraction(-H[1][0], detH),
         Fraction(H[0][0], detH)],
    ]

    def mm(A, B):
        out = [
            [Fraction(0) for _ in range(len(B[0]))]
            for _ in range(len(A))
        ]
        for i in range(len(A)):
            for k in range(len(B)):
                if A[i][k] == 0:
                    continue
                for j in range(len(B[0])):
                    out[i][j] += (
                        Fraction(A[i][k]) * Fraction(B[k][j])
                    )
        return out

    corr = mm(B, mm(Hinv, C))

    K = [
        [-corr[i][j] for j in range(CORE_N)]
        for i in range(CORE_N)
    ]

    K_int = [
        [
            x.numerator
            if x.denominator == 1
            else (_raise_nonintegral())
            for x in row
        ]
        for row in K
    ]

    A = [
        [
            F[i][j] + K_int[i][j]
            for j in range(CORE_N)
        ]
        for i in range(CORE_N)
    ]

    return F, K_int, A, H


def _raise_nonintegral():
    raise ArithmeticError("Non-integral Schur correction.")


# ============================================================================
# F2 LINEAR ALGEBRA
# ============================================================================

def mod2(A):
    return [[int(x) & 1 for x in row] for row in A]


def transpose(A):
    return [list(x) for x in zip(*A)]


def rref2(A):
    M = mod2(A)

    if not M:
        return [], []

    rows = len(M)
    cols = len(M[0])

    pivots = []
    pr = 0

    for c in range(cols):

        p = None

        for r in range(pr, rows):
            if M[r][c]:
                p = r
                break

        if p is None:
            continue

        M[pr], M[p] = M[p], M[pr]

        for r in range(rows):
            if r != pr and M[r][c]:
                for j in range(cols):
                    M[r][j] ^= M[pr][j]

        pivots.append(c)
        pr += 1

        if pr == rows:
            break

    return M, pivots


def rank2(A):
    return len(rref2(A)[1])


def row_basis2(A):
    R, _ = rref2(A)
    return [r for r in R if any(r)]


def nullspace2(A):
    if not A:
        return []

    R, pivots = rref2(A)
    n = len(A[0])
    P = set(pivots)

    free = [j for j in range(n) if j not in P]

    out = []

    for f in free:

        v = [0] * n
        v[f] = 1

        for i, p in enumerate(pivots):
            if R[i][f]:
                v[p] = 1

        out.append(v)

    return out


def in_span2(v, B):
    if not B:
        return not any(v)

    if len(v) != len(B[0]):
        raise ValueError("dimension mismatch")

    return rank2(B) == rank2(B + [v])


def dot2(a, b):
    if len(a) != len(b):
        raise ValueError("dot-product dimension mismatch")
    s = 0
    for x, y in zip(a, b):
        s ^= (x & 1) & (y & 1)
    return s


def cofactor_matrix(A):
    n = len(A)
    C = []

    for i in range(n):

        row = []

        for j in range(n):

            rows = [
                r for r in range(n)
                if r != i
            ]

            cols = [
                c for c in range(n)
                if c != j
            ]

            minor = [
                [A[r][c] for c in cols]
                for r in rows
            ]

            x = bareiss_det(minor)

            if (i+j) & 1:
                x = -x

            row.append(x)

        C.append(row)

    return C


def normalize_cofactor(C):
    return [
        [
            x // DELTA_11
            if x % DELTA_11 == 0
            else (_raise_nondivisible())
            for x in row
        ]
        for row in C
    ]


def _raise_nondivisible():
    raise ArithmeticError(
        "Delta_11 does not divide cofactor."
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 166 — EXACT MOD-2 SCHUR QUOTIENT / "
        "ADJUGATE PAIRING AUDIT"
    )
    print("=" * 78)

    M = build_full_system()
    F, K, A, H = build_core(M)

    # ------------------------------------------------------------------
    # 1. BASE RANKS
    # ------------------------------------------------------------------

    M2 = mod2(M)
    F2 = mod2(F)
    K2 = mod2(K)
    A2 = mod2(A)

    rankM = rank2(M2)
    rankF = rank2(F2)
    rankK = rank2(K2)
    rankA = rank2(A2)

    print()
    print("=" * 78)
    print("1. EXACT MOD-2 RANKS")
    print("=" * 78)

    print(f"  rank(M)={rankM}")
    print(f"  rank(F_ret)={rankF}")
    print(f"  rank(K)={rankK}")
    print(f"  rank(A)={rankA}")
    print(f"  F_minus_A={rankF-rankA}")

    # ------------------------------------------------------------------
    # 2. ROW / COLUMN SPACES
    # ------------------------------------------------------------------

    F_rows = row_basis2(F2)
    A_rows = row_basis2(A2)
    K_rows = row_basis2(K2)

    F_cols = row_basis2(transpose(F2))
    A_cols = row_basis2(transpose(A2))
    K_cols = row_basis2(transpose(K2))

    print()
    print("=" * 78)
    print("2. SPACE DIMENSIONS")
    print("=" * 78)

    print(f"  Row(F)={len(F_rows)}")
    print(f"  Row(A)={len(A_rows)}")
    print(f"  Row(K)={len(K_rows)}")

    print(f"  Col(F)={len(F_cols)}")
    print(f"  Col(A)={len(A_cols)}")
    print(f"  Col(K)={len(K_cols)}")

    # ------------------------------------------------------------------
    # 3. QUOTIENT / INTERSECTION
    # ------------------------------------------------------------------

    row_FA_sum = rank2(F_rows + A_rows)
    row_FK_sum = rank2(F_rows + K_rows)
    row_AK_sum = rank2(A_rows + K_rows)

    col_FA_sum = rank2(F_cols + A_cols)
    col_FK_sum = rank2(F_cols + K_cols)
    col_AK_sum = rank2(A_cols + K_cols)

    print()
    print("=" * 78)
    print("3. EXACT QUOTIENT / INTERSECTION DATA")
    print("=" * 78)

    print(
        f"  dim(Row(F) cap Row(A))="
        f"{len(F_rows)+len(A_rows)-row_FA_sum}"
    )

    print(
        f"  dim(Row(F) cap Row(K))="
        f"{len(F_rows)+len(K_rows)-row_FK_sum}"
    )

    print(
        f"  dim(Row(A) cap Row(K))="
        f"{len(A_rows)+len(K_rows)-row_AK_sum}"
    )

    print(
        f"  dim(Col(F) cap Col(A))="
        f"{len(F_cols)+len(A_cols)-col_FA_sum}"
    )

    print(
        f"  dim(Col(F) cap Col(K))="
        f"{len(F_cols)+len(K_cols)-col_FK_sum}"
    )

    print(
        f"  dim(Col(A) cap Col(K))="
        f"{len(A_cols)+len(K_cols)-col_AK_sum}"
    )

    # ------------------------------------------------------------------
    # 4. KERNELS OF A
    # ------------------------------------------------------------------

    leftA = nullspace2(transpose(A2))
    rightA = nullspace2(A2)

    print()
    print("=" * 78)
    print("4. CORE KERNELS")
    print("=" * 78)

    print(
        f"  left_nullity(A)={len(leftA)}"
    )

    print(
        f"  right_nullity(A)={len(rightA)}"
    )

    # ------------------------------------------------------------------
    # 5. ADJUGATE DIRECTIONS
    # ------------------------------------------------------------------

    C = cofactor_matrix(A)
    R = normalize_cofactor(C)
    AdjR = transpose(R)
    AdjR2 = mod2(AdjR)

    adj_row_space = row_basis2(AdjR2)
    adj_col_space = row_basis2(
        transpose(AdjR2)
    )

    if len(adj_row_space) != 1:
        raise ArithmeticError(
            "Expected rank-one adjugate row space."
        )

    if len(adj_col_space) != 1:
        raise ArithmeticError(
            "Expected rank-one adjugate column space."
        )

    adj_row = adj_row_space[0]
    adj_col = adj_col_space[0]

    print()
    print("=" * 78)
    print("5. ADJUGATE DIRECTIONS")
    print("=" * 78)

    print(
        f"  adj_row={adj_row}"
    )

    print(
        f"  adj_col={adj_col}"
    )

    print(
        f"  adj_row_in_left_kernel="
        f"{in_span2(adj_row, leftA)}"
    )

    print(
        f"  adj_col_in_right_kernel="
        f"{in_span2(adj_col, rightA)}"
    )

    # ------------------------------------------------------------------
    # 6. QUOTIENT REPRESENTATIVES
    # ------------------------------------------------------------------

    row_quotient_dim = (
        len(F_rows) - len(A_rows)
    )

    col_quotient_dim = (
        len(F_cols) - len(A_cols)
    )

    row_F_not_A = [
        v for v in F_rows
        if not in_span2(v, A_rows)
    ]

    col_F_not_A = [
        v for v in F_cols
        if not in_span2(v, A_cols)
    ]

    print()
    print("=" * 78)
    print("6. EXPLICIT QUOTIENT REPRESENTATIVES")
    print("=" * 78)

    print(
        f"  Row quotient dimension="
        f"{row_quotient_dim}"
    )

    for i, v in enumerate(row_F_not_A):
        print(
            f"  row_quotient_candidate_{i}={v}"
        )

    print(
        f"  Col quotient dimension="
        f"{col_quotient_dim}"
    )

    for i, v in enumerate(col_F_not_A):
        print(
            f"  col_quotient_candidate_{i}={v}"
        )

    # ------------------------------------------------------------------
    # 7. PAIRING MATRIX
    # ------------------------------------------------------------------

    """
    The row-space quotient is a subspace of the dual ambient space
    relative to the right kernel.

    The column-space quotient is paired with the left kernel.

    We compute all pairings between candidate quotient vectors and
    bases of the complementary kernel spaces.

    The purpose is descriptive: determine whether the quotient direction
    is detected by the distinguished adjugate line.
    """

    row_pairings = []

    for q in row_F_not_A:

        row_pairings.append(
            [
                dot2(q, v)
                for v in leftA
            ]
        )

    col_pairings = []

    for q in col_F_not_A:

        col_pairings.append(
            [
                dot2(q, v)
                for v in rightA
            ]
        )

    print()
    print("=" * 78)
    print("7. QUOTIENT / KERNEL PAIRINGS")
    print("=" * 78)

    print(
        f"  row_quotient_pairings="
        f"{row_pairings}"
    )

    print(
        f"  col_quotient_pairings="
        f"{col_pairings}"
    )

    # Pair quotient representatives directly with adjugate directions.
    row_adj_pairing = [
        dot2(q, adj_col)
        for q in row_F_not_A
    ]

    col_adj_pairing = [
        dot2(q, adj_row)
        for q in col_F_not_A
    ]

    print(
        f"  row_quotient_x_adj_col="
        f"{row_adj_pairing}"
    )

    print(
        f"  col_quotient_x_adj_row="
        f"{col_adj_pairing}"
    )

    # ------------------------------------------------------------------
    # 8. SCHUR MATRIX CHECK
    # ------------------------------------------------------------------

    H_mod2 = mod2(H)

    print()
    print("=" * 78)
    print("8. SCHUR BLOCK")
    print("=" * 78)

    print(
        f"  H={H}"
    )

    print(
        f"  det(H)="
        f"{H[0][0]*H[1][1]-H[0][1]*H[1][0]}"
    )

    print(
        f"  H_mod2={H_mod2}"
    )

    # ------------------------------------------------------------------
    # 9. INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The previous experiment established the exact direct-sum geometry

    Row(F) = Row(A) ⊕ Row(K)

and

    Col(F) = Col(A) ⊕ Col(K)

modulo 2.

Experiment 166 asks what the missing one-dimensional direction means
from the dual/kernel point of view.

The core adjugate line is a canonical one-dimensional subspace of

    ker(A^T)

and

    ker(A).

The quotient directions

    Row(F)/Row(A),
    Col(F)/Col(A)

are tested against these kernel directions using the exact F_2 dot
product.

A nonzero pairing would identify the Schur-lost direction by a
canonical dual functional.

A zero pairing would show that the quotient direction is invisible
to the adjugate line and that the two structures are independent at
this level.

This is a finite exact linear-algebra question only.
"""
    )

    # ------------------------------------------------------------------
    # 10. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        rankM == 4
        and rankF == 4
        and rankK == 1
        and rankA == 3
        and row_quotient_dim == 1
        and col_quotient_dim == 1
        and len(adj_row_space) == 1
        and len(adj_col_space) == 1
        and in_span2(adj_row, leftA)
        and in_span2(adj_col, rightA)
    )

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  rank_M_4={rankM == 4}"
    )

    print(
        f"  rank_F_4={rankF == 4}"
    )

    print(
        f"  rank_K_1={rankK == 1}"
    )

    print(
        f"  rank_A_3={rankA == 3}"
    )

    print(
        f"  row_quotient_dimension_1="
        f"{row_quotient_dim == 1}"
    )

    print(
        f"  col_quotient_dimension_1="
        f"{col_quotient_dim == 1}"
    )

    print(
        f"  adjugate_rank_one="
        f"{len(adj_row_space) == 1 and len(adj_col_space) == 1}"
    )

    print(
        f"  adj_row_kernel_membership="
        f"{in_span2(adj_row, leftA)}"
    )

    print(
        f"  adj_col_kernel_membership="
        f"{in_span2(adj_col, rightA)}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 166 COMPLETE")


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

