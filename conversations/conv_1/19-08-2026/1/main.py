#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 170R — EXACT FIRST-LIFT FACTOR LOCALIZATION /
                SCHUR CORRECTION ROW-SELECTION AUDIT
==============================================================================

Goal:

    Experiment 169S showed

        L1 = (B adj(H))/2       rank 1 mod 2
        R1 = (adj(H) C)/2       rank 2 mod 2
        K  = -L1 C = -B R1      rank 1 mod 2

    with

        B : 12 x 2
        H : 2 x 2
        C : 2 x 12

This experiment determines whether the rank-one K mod 2 comes from
one selected row of C, and separately describes the rank-two R1 side.

No SymPy.
No floating point.
No extrapolation.
No recurrence fitting.
No connection to the original (p,q)-kernel.
"""

from __future__ import annotations

import sys
from fractions import Fraction


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

FULL_N = 14
CORE_N = 12

BOUNDARY_ROWS = (5, 10)
BOUNDARY_COLS = (12, 13)


# ============================================================================
# BASIC HELPERS
# ============================================================================

def q_value(p: int, r: int) -> int:
    if r < 0 or r >= len(Q[p]):
        return 0
    return Q[p][r]


def basis(p: int, d: int) -> list[int]:
    return [
        1,
        d,
        d * d,
        p,
        p * d,
        p * p,
    ]


def assert_rectangular(A, name="matrix"):
    if not A:
        return

    width = len(A[0])

    for i, row in enumerate(A):
        if len(row) != width:
            raise ValueError(
                f"{name}: row {i} width {len(row)}, expected {width}."
            )


def to_mod2(A):
    assert_rectangular(A)
    return [
        [int(x) & 1 for x in row]
        for row in A
    ]


# ============================================================================
# MATRIX OPERATIONS
# ============================================================================

def transpose(A):
    assert_rectangular(A, "transpose input")

    if not A:
        return []

    return [
        [A[i][j] for i in range(len(A))]
        for j in range(len(A[0]))
    ]


def matmul(A, B):
    assert_rectangular(A, "left matrix")
    assert_rectangular(B, "right matrix")

    if not A or not B:
        return []

    if len(A[0]) != len(B):
        raise ValueError(
            f"Matrix multiplication mismatch: "
            f"{len(A)}x{len(A[0])} times "
            f"{len(B)}x{len(B[0])}."
        )

    out = [
        [Fraction(0) for _ in range(len(B[0]))]
        for _ in range(len(A))
    ]

    for i in range(len(A)):
        for k in range(len(B)):

            a = Fraction(A[i][k])

            if a == 0:
                continue

            for j in range(len(B[0])):
                out[i][j] += a * Fraction(B[k][j])

    return out


def matrix_to_int(A, name):
    out = []

    for i, row in enumerate(A):

        new_row = []

        for j, x in enumerate(row):

            x = Fraction(x)

            if x.denominator != 1:
                raise ArithmeticError(
                    f"{name}[{i}][{j}] non-integral: {x}"
                )

            new_row.append(x.numerator)

        out.append(new_row)

    return out


def divide_exact(A, d, name):

    out = []

    for i, row in enumerate(A):

        new_row = []

        for j, x in enumerate(row):

            if x % d != 0:
                raise ArithmeticError(
                    f"{name}[{i}][{j}]={x} "
                    f"is not divisible by {d}."
                )

            new_row.append(x // d)

        out.append(new_row)

    return out


# ============================================================================
# FULL SYSTEM
# ============================================================================

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
                row[6 + j] = b[j] * q1

            if d == 0:
                row[12] = 1
                row[13] = p

            M.append(row)

    if len(M) != FULL_N:
        raise ArithmeticError(
            f"Expected {FULL_N} rows, got {len(M)}."
        )

    return M


def build_blocks(M):

    retained = [
        i
        for i in range(FULL_N)
        if i not in BOUNDARY_ROWS
    ]

    retained_rows = [
        M[i]
        for i in retained
    ]

    boundary_rows = [
        M[i]
        for i in BOUNDARY_ROWS
    ]

    # 12 x 12
    F = [
        [
            row[j]
            for j in range(CORE_N)
        ]
        for row in retained_rows
    ]

    # 12 x 2
    B = [
        [
            row[j]
            for j in BOUNDARY_COLS
        ]
        for row in retained_rows
    ]

    # 2 x 2
    H = [
        [
            row[j]
            for j in BOUNDARY_COLS
        ]
        for row in boundary_rows
    ]

    # 2 x 12
    C = [
        [
            row[j]
            for j in range(CORE_N)
        ]
        for row in boundary_rows
    ]

    det_h = (
        H[0][0] * H[1][1]
        - H[0][1] * H[1][0]
    )

    if det_h != 2:
        raise ArithmeticError(
            f"Expected det(H)=2, got {det_h}."
        )

    adj_h = [
        [
            H[1][1],
            -H[0][1],
        ],
        [
            -H[1][0],
            H[0][0],
        ],
    ]

    return F, B, H, adj_h, C, det_h


# ============================================================================
# MOD-2 LINEAR ALGEBRA
# ============================================================================

def rref_mod2(A):

    M = to_mod2(A)

    if not M:
        return [], []

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

                for j in range(cols):
                    M[r][j] ^= M[pivot_row][j]

        pivots.append(col)
        pivot_row += 1

        if pivot_row == rows:
            break

    return M, pivots


def rank_mod2(A):
    if not A:
        return 0

    return len(rref_mod2(A)[1])


def row_basis_mod2(A):

    R, _ = rref_mod2(A)

    return [
        row
        for row in R
        if any(row)
    ]


def matrix_support(A):

    return [
        (i, j)
        for i in range(len(A))
        for j in range(len(A[0]))
        if A[i][j]
    ]


def vector_support(v):

    return [
        i
        for i, x in enumerate(v)
        if x
    ]


# ============================================================================
# RANK-ONE FACTORIZATION FOR RECTANGULAR MATRICES
# ============================================================================

def rank_one_factorization_mod2(A):

    """
    For m x n matrix A over F_2, if rank(A)=1, find

        A = u v^T.

    Returns (u, v), where

        len(u) = m,
        len(v) = n.
    """

    A2 = to_mod2(A)

    if rank_mod2(A2) != 1:
        return None, None

    m = len(A2)
    n = len(A2[0])

    pivot_row = None
    pivot_col = None

    for i in range(m):
        for j in range(n):
            if A2[i][j]:
                pivot_row = i
                pivot_col = j
                break

        if pivot_row is not None:
            break

    if pivot_row is None:
        return None, None

    u = [
        A2[i][pivot_col]
        for i in range(m)
    ]

    v = [
        A2[pivot_row][j]
        for j in range(n)
    ]

    predicted = [
        [
            u[i] & v[j]
            for j in range(n)
        ]
        for i in range(m)
    ]

    if predicted != A2:
        raise ArithmeticError(
            "Rectangular rank-one factorization reconstruction failed."
        )

    return u, v


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 170R — EXACT FIRST-LIFT FACTOR LOCALIZATION / "
        "SCHUR CORRECTION ROW-SELECTION AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # BUILD
    # ------------------------------------------------------------------

    M = build_full_system()

    F, B, H, adj_h, C, det_h = (
        build_blocks(M)
    )

    B2 = to_mod2(B)
    C2 = to_mod2(C)
    H2 = to_mod2(H)
    adj_h2 = to_mod2(adj_h)

    # ------------------------------------------------------------------
    # 1. BASELINE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. BASELINE")
    print("=" * 78)

    print(
        f"  rank(B_mod2)={rank_mod2(B2)}"
    )

    print(
        f"  rank(C_mod2)={rank_mod2(C2)}"
    )

    print(
        f"  rank(H_mod2)={rank_mod2(H2)}"
    )

    print(
        f"  rank(adj_H_mod2)={rank_mod2(adj_h2)}"
    )

    print(
        f"  det_H={det_h}"
    )

    # ------------------------------------------------------------------
    # 2. FIRST LIFTS
    # ------------------------------------------------------------------

    L = matrix_to_int(
        matmul(B, adj_h),
        "L"
    )

    R = matrix_to_int(
        matmul(adj_h, C),
        "R"
    )

    L1 = divide_exact(
        L,
        2,
        "L"
    )

    R1 = divide_exact(
        R,
        2,
        "R"
    )

    L1_2 = to_mod2(L1)
    R1_2 = to_mod2(R1)

    rank_L1 = rank_mod2(L1_2)
    rank_R1 = rank_mod2(R1_2)

    print()
    print("=" * 78)
    print("2. FIRST 2-ADIC LIFTS")
    print("=" * 78)

    print(
        f"  L_shape=12x2"
    )

    print(
        f"  R_shape=2x12"
    )

    print(
        f"  L_even=True"
    )

    print(
        f"  R_even=True"
    )

    print(
        f"  rank_L1_mod2={rank_L1}"
    )

    print(
        f"  rank_R1_mod2={rank_R1}"
    )

    print("  L1_mod2:")

    for row in L1_2:
        print(
            "    "
            + " ".join(str(x) for x in row)
        )

    print("  R1_mod2:")

    for row in R1_2:
        print(
            "    "
            + " ".join(str(x) for x in row)
        )

    # ------------------------------------------------------------------
    # 3. K
    # ------------------------------------------------------------------

    K_left = matrix_to_int(
        matmul(L1, C),
        "L1*C"
    )

    K = [
        [
            -x
            for x in row
        ]
        for row in K_left
    ]

    K_right = matrix_to_int(
        matmul(B, R1),
        "B*R1"
    )

    K_right = [
        [
            -x
            for x in row
        ]
        for row in K_right
    ]

    K2 = to_mod2(K)

    left_right_exact = (
        K == K_right
    )

    print()
    print("=" * 78)
    print("3. EXACT K RECOVERY")
    print("=" * 78)

    print(
        f"  K_from_L1_exact=True"
    )

    print(
        f"  K_from_R1_exact=True"
    )

    print(
        f"  left_right_exact={left_right_exact}"
    )

    print(
        f"  rank_K_mod2={rank_mod2(K2)}"
    )

    print(
        f"  K_support={matrix_support(K2)}"
    )

    # ------------------------------------------------------------------
    # 4. L1 MATRIX-UNIT TEST
    # ------------------------------------------------------------------

    L1_support = matrix_support(
        L1_2
    )

    is_matrix_unit = (
        len(L1_support) == 1
    )

    selected_row = None
    selected_col = None

    if is_matrix_unit:

        selected_row, selected_col = (
            L1_support[0]
        )

    print()
    print("=" * 78)
    print("4. L1 MATRIX-UNIT TEST")
    print("=" * 78)

    print(
        f"  L1_support={L1_support}"
    )

    print(
        f"  L1_is_matrix_unit="
        f"{is_matrix_unit}"
    )

    if is_matrix_unit:

        print(
            f"  selected_row={selected_row}"
        )

        print(
            f"  selected_C_row={selected_col}"
        )

    # ------------------------------------------------------------------
    # 5. SELECTED C ROW
    # ------------------------------------------------------------------

    selected_C_row_matches = False
    selected_C_row = None

    if is_matrix_unit:

        selected_C_row = C2[selected_col]

        predicted = [
            [0 for _ in range(CORE_N)]
            for _ in range(CORE_N)
        ]

        for j in range(CORE_N):
            predicted[selected_row][j] = (
                selected_C_row[j]
            )

        selected_C_row_matches = (
            predicted == K2
        )

    print()
    print("=" * 78)
    print("5. SELECTED-C-ROW TEST")
    print("=" * 78)

    if selected_C_row is not None:

        print(
            f"  C_row_index={selected_col}"
        )

        print(
            f"  C_row_support="
            f"{vector_support(selected_C_row)}"
        )

    print(
        f"  selected_C_row_reconstructs_K="
        f"{selected_C_row_matches}"
    )

    # ------------------------------------------------------------------
    # 6. C STRUCTURE
    # ------------------------------------------------------------------

    C_rank = rank_mod2(C2)
    C_rows = row_basis_mod2(C2)
    C_cols = row_basis_mod2(
        transpose(C2)
    )

    K_rank = rank_mod2(K2)

    print()
    print("=" * 78)
    print("6. C / K STRUCTURE")
    print("=" * 78)

    print(
        f"  rank_C_mod2={C_rank}"
    )

    print(
        f"  C_row_basis={C_rows}"
    )

    print(
        f"  C_column_space_rank={len(C_cols)}"
    )

    print(
        f"  rank_K_mod2={K_rank}"
    )

    # ------------------------------------------------------------------
    # 7. RIGHT-LIFT ASYMMETRY
    # ------------------------------------------------------------------

    B_R1 = matrix_to_int(
        matmul(B, R1),
        "B*R1"
    )

    B_R1_2 = to_mod2(B_R1)

    print()
    print("=" * 78)
    print("7. RIGHT-LIFT ASYMMETRY")
    print("=" * 78)

    print(
        f"  rank_B_mod2={rank_mod2(B2)}"
    )

    print(
        f"  rank_R1_mod2={rank_R1}"
    )

    print(
        f"  rank_B_R1_mod2={rank_mod2(B_R1_2)}"
    )

    print(
        f"  B_R1_matches_K_mod2="
        f"{B_R1_2 == K2}"
    )

    print(
        "  B_R1_mod2:"
    )

    for row in B_R1_2:
        print(
            "    "
            + " ".join(str(x) for x in row)
        )

    # ------------------------------------------------------------------
    # 8. R1 ROW-SPACE
    # ------------------------------------------------------------------

    R1_rows = row_basis_mod2(R1_2)
    R1_cols = row_basis_mod2(
        transpose(R1_2)
    )

    B_rows = row_basis_mod2(B2)
    B_cols = row_basis_mod2(
        transpose(B2)
    )

    print()
    print("=" * 78)
    print("8. RIGHT-LIFT SPACE DATA")
    print("=" * 78)

    print(
        f"  R1_row_space_rank={len(R1_rows)}"
    )

    print(
        f"  R1_column_space_rank={len(R1_cols)}"
    )

    print(
        f"  B_row_space_rank={len(B_rows)}"
    )

    print(
        f"  B_column_space_rank={len(B_cols)}"
    )

    print(
        f"  R1_row_basis={R1_rows}"
    )

    print(
        f"  R1_column_basis={R1_cols}"
    )

    # ------------------------------------------------------------------
    # 9. CORE RANK
    # ------------------------------------------------------------------

    A = [
        [
            F[i][j] + K[i][j]
            for j in range(CORE_N)
        ]
        for i in range(CORE_N)
    ]

    rank_F = rank_mod2(F)
    rank_A = rank_mod2(A)

    print()
    print("=" * 78)
    print("9. CORE RANK TRANSITION")
    print("=" * 78)

    print(
        f"  rank_F_ret_mod2={rank_F}"
    )

    print(
        f"  rank_K_mod2={K_rank}"
    )

    print(
        f"  rank_A_mod2={rank_A}"
    )

    print(
        f"  rank_drop={rank_F-rank_A}"
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
Experiment 169S established:

    rank(L1 mod 2) = 1
    rank(R1 mod 2) = 2
    rank(K mod 2)  = 1.

The left factor has shape 12x2.

If its unique nonzero entry is at (i,j), then

    L1 mod 2 = e_i e_j^T

and therefore

    K mod 2 = e_i * C_j mod 2.

That would mean the three-entry K support is not an independent
phenomenon: it is exactly one selected row of C injected into one
selected row of the core.

The right side behaves differently. R1 has rank two, but

    B mod 2

has rank one, and the product

    B R1

still has rank one.

Thus the right lift contains an additional direction that disappears
when acted on by B.

This experiment separates:

    left-side rank-one localization

from

    right-side rank-two-then-collapse behavior.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    det_ok = (det_h == 2)
    L1_ok = (rank_L1 == 1)
    R1_ok = (rank_R1 == 2)
    B_ok = (rank_mod2(B2) == 1)
    C_ok = (C_rank == 1)
    K_ok = (K_rank == 1)
    F_ok = (rank_F == 4)
    A_ok = (rank_A == 3)

    final_ok = (
        det_ok
        and L1_ok
        and R1_ok
        and B_ok
        and C_ok
        and K_ok
        and F_ok
        and A_ok
        and left_right_exact
        and selected_C_row_matches
        and (B_R1_2 == K2)
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  det_H_equals_2={det_ok}"
    )

    print(
        f"  L1_rank_one={L1_ok}"
    )

    print(
        f"  R1_rank_two={R1_ok}"
    )

    print(
        f"  B_rank_one={B_ok}"
    )

    print(
        f"  C_rank_one={C_ok}"
    )

    print(
        f"  K_rank_one={K_ok}"
    )

    print(
        f"  F_rank_four={F_ok}"
    )

    print(
        f"  A_rank_three={A_ok}"
    )

    print(
        f"  left_right_K_recovery={left_right_exact}"
    )

    print(
        f"  selected_C_row_reconstructs_K="
        f"{selected_C_row_matches}"
    )

    print(
        f"  B_R1_reproduces_K="
        f"{B_R1_2 == K2}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 170R COMPLETE")


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