#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 167R — EXACT SCHUR-QUOTIENT / RANK-ONE-CORRECTION
CANONICALIZATION AUDIT
==============================================================================

Corrected replacement for Experiment 167.

The previous crash came from using Python list concatenation:

    row_q + K_row

instead of addition in F_2:

    row_q XOR K_row.

Everything below uses explicit mod-2 vector addition.

Established exact facts:

    rank_F2(F_ret) = 4
    rank_F2(A)     = 3
    rank_F2(K)     = 1

with

    A = F_ret + K.

The experiment asks whether the rank-one Schur correction K itself
represents the unique quotient directions

    Row(F_ret) / Row(A)

and

    Col(F_ret) / Col(A).

No mixed-dimensional comparisons.

No SymPy.
No floating point.
No extrapolation.
No recurrence fitting.
No connection to the original (p,q)-kernel.
"""

from __future__ import annotations

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

FULL_N = 14
CORE_N = 12

BOUNDARY_ROWS = (5, 10)
BOUNDARY_COLS = (12, 13)

DELTA_11 = 9512681472


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
                f"{name}: row {i} has width {len(row)}, "
                f"expected {width}."
            )


# ============================================================================
# F_2 VECTOR OPERATIONS
# ============================================================================

def add2(v, w):
    """Vector addition in F_2 = componentwise XOR."""

    if len(v) != len(w):
        raise ValueError(
            f"F2 vector dimension mismatch: "
            f"{len(v)} versus {len(w)}."
        )

    return [
        (a ^ b)
        for a, b in zip(v, w)
    ]


def scalar_mul2(c, v):
    c &= 1

    if c == 0:
        return [0] * len(v)

    return list(v)


# ============================================================================
# BAREISS
# ============================================================================

def bareiss_det(A):

    if not A:
        return 1

    assert_rectangular(A, "Bareiss input")

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
                            "Bareiss exact division failed."
                        )

                    value //= previous

                M[i][j] = value

        for i in range(k + 1, n):
            M[i][k] = 0

        previous = pivot

    return sign * M[-1][-1]


# ============================================================================
# MATRIX OPERATIONS
# ============================================================================

def transpose(A):

    assert_rectangular(A, "transpose input")

    if not A:
        return []

    return [
        [
            A[i][j]
            for i in range(len(A))
        ]
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
        [
            Fraction(0)
            for _ in range(len(B[0]))
        ]
        for _ in range(len(A))
    ]

    for i in range(len(A)):

        for k in range(len(B)):

            if A[i][k] == 0:
                continue

            a = Fraction(A[i][k])

            for j in range(len(B[0])):

                out[i][j] += (
                    a * Fraction(B[k][j])
                )

    return out


def matrix_to_integer(A, name):

    out = []

    for i, row in enumerate(A):

        new_row = []

        for j, x in enumerate(row):

            x = Fraction(x)

            if x.denominator != 1:
                raise ArithmeticError(
                    f"{name}[{i}][{j}] is non-integral: {x}"
                )

            new_row.append(x.numerator)

        out.append(new_row)

    return out


# ============================================================================
# BUILD FULL SYSTEM
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
            f"Expected {FULL_N} equations, got {len(M)}."
        )

    assert_rectangular(M, "full system")

    return M


# ============================================================================
# EXACT SCHUR COMPLEMENT
# ============================================================================

def build_schur_data(M):

    retained_indices = [
        i
        for i in range(FULL_N)
        if i not in BOUNDARY_ROWS
    ]

    retained_rows = [
        M[i]
        for i in retained_indices
    ]

    boundary_rows = [
        M[i]
        for i in BOUNDARY_ROWS
    ]

    F = [
        [
            row[j]
            for j in range(CORE_N)
        ]
        for row in retained_rows
    ]

    B = [
        [
            row[j]
            for j in BOUNDARY_COLS
        ]
        for row in retained_rows
    ]

    H = [
        [
            row[j]
            for j in BOUNDARY_COLS
        ]
        for row in boundary_rows
    ]

    C = [
        [
            row[j]
            for j in range(CORE_N)
        ]
        for row in boundary_rows
    ]

    det_H = (
        H[0][0] * H[1][1]
        - H[0][1] * H[1][0]
    )

    if det_H == 0:
        raise ArithmeticError(
            "Boundary block H is singular over Q."
        )

    H_inv = [
        [
            Fraction(H[1][1], det_H),
            Fraction(-H[0][1], det_H),
        ],
        [
            Fraction(-H[1][0], det_H),
            Fraction(H[0][0], det_H),
        ],
    ]

    B_Hinv = matmul(B, H_inv)
    B_Hinv_C = matmul(B_Hinv, C)

    K_fraction = [
        [
            -B_Hinv_C[i][j]
            for j in range(CORE_N)
        ]
        for i in range(CORE_N)
    ]

    K = matrix_to_integer(
        K_fraction,
        "Schur correction K",
    )

    A = [
        [
            F[i][j] + K[i][j]
            for j in range(CORE_N)
        ]
        for i in range(CORE_N)
    ]

    return (
        F,
        K,
        A,
        H,
        retained_indices,
    )


# ============================================================================
# MOD-2 LINEAR ALGEBRA
# ============================================================================

def to_mod2(A):

    assert_rectangular(A)

    return [
        [
            int(x) & 1
            for x in row
        ]
        for row in A
    ]


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

        for r in range(
            pivot_row,
            rows,
        ):

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

    if not A:
        return []

    R, _ = rref_mod2(A)

    return [
        row
        for row in R
        if any(row)
    ]


def nullspace_mod2(A):

    if not A:
        return []

    R, pivots = rref_mod2(A)

    n = len(A[0])

    pivot_set = set(pivots)

    free_cols = [
        j
        for j in range(n)
        if j not in pivot_set
    ]

    basis = []

    for free in free_cols:

        v = [0] * n
        v[free] = 1

        for i, pivot in enumerate(pivots):

            if R[i][free]:
                v[pivot] = 1

        basis.append(v)

    return basis


def in_span_mod2(v, basis):

    if not basis:
        return not any(v)

    width = len(basis[0])

    if len(v) != width:
        raise ValueError(
            "Dimension mismatch in span test."
        )

    if any(
        len(b) != width
        for b in basis
    ):
        raise ValueError(
            "Nonuniform basis dimension."
        )

    return (
        rank_mod2(basis)
        == rank_mod2(
            basis + [v]
        )
    )


def quotient_representative(
    super_basis,
    sub_basis,
):

    if not super_basis:
        return None

    width = len(super_basis[0])

    if any(
        len(v) != width
        for v in sub_basis
    ):
        raise ValueError(
            "Subspace dimension mismatch."
        )

    for v in super_basis:

        if not in_span_mod2(
            v,
            sub_basis,
        ):
            return v

    return None


# ============================================================================
# COFACTOR / ADJUGATE
# ============================================================================

def maximal_minor(
    A,
    omit_row,
    omit_col,
):

    rows = [
        i
        for i in range(len(A))
        if i != omit_row
    ]

    cols = [
        j
        for j in range(len(A[0]))
        if j != omit_col
    ]

    return [
        [
            A[i][j]
            for j in cols
        ]
        for i in rows
    ]


def cofactor_matrix(A):

    n = len(A)

    C = []

    for i in range(n):

        row = []

        for j in range(n):

            value = bareiss_det(
                maximal_minor(
                    A,
                    i,
                    j,
                )
            )

            if (i + j) & 1:
                value = -value

            row.append(value)

        C.append(row)

    return C


def normalize_cofactor(C):

    R = []

    for row in C:

        new_row = []

        for x in row:

            if x % DELTA_11 != 0:
                raise ArithmeticError(
                    "Delta_11 divisibility failure."
                )

            new_row.append(
                x // DELTA_11
            )

        R.append(new_row)

    return R


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 167R — EXACT SCHUR-QUOTIENT / "
        "RANK-ONE-CORRECTION CANONICALIZATION"
    )
    print("=" * 78)

    M = build_full_system()

    (
        F,
        K,
        A,
        H,
        retained_indices,
    ) = build_schur_data(M)

    F2 = to_mod2(F)
    K2 = to_mod2(K)
    A2 = to_mod2(A)
    M2 = to_mod2(M)

    # ------------------------------------------------------------------
    # 1. BASELINE
    # ------------------------------------------------------------------

    rank_M = rank_mod2(M2)
    rank_F = rank_mod2(F2)
    rank_K = rank_mod2(K2)
    rank_A = rank_mod2(A2)

    F_rows = row_basis_mod2(F2)
    A_rows = row_basis_mod2(A2)
    K_rows = row_basis_mod2(K2)

    F_cols = row_basis_mod2(
        transpose(F2)
    )

    A_cols = row_basis_mod2(
        transpose(A2)
    )

    K_cols = row_basis_mod2(
        transpose(K2)
    )

    print()
    print("=" * 78)
    print("1. EXACT BASELINE")
    print("=" * 78)

    print(
        f"  rank(M)={rank_M}"
    )

    print(
        f"  rank(F_ret)={rank_F}"
    )

    print(
        f"  rank(K)={rank_K}"
    )

    print(
        f"  rank(A)={rank_A}"
    )

    print(
        f"  rank_drop={rank_F-rank_A}"
    )

    # ------------------------------------------------------------------
    # 2. QUOTIENT REPRESENTATIVES
    # ------------------------------------------------------------------

    row_q = quotient_representative(
        F_rows,
        A_rows,
    )

    col_q = quotient_representative(
        F_cols,
        A_cols,
    )

    if row_q is None:
        raise ArithmeticError(
            "Could not find row quotient representative."
        )

    if col_q is None:
        raise ArithmeticError(
            "Could not find column quotient representative."
        )

    print()
    print("=" * 78)
    print("2. CANONICAL QUOTIENT REPRESENTATIVES")
    print("=" * 78)

    print(
        f"  dim(Row(F)/Row(A))="
        f"{len(F_rows)-len(A_rows)}"
    )

    print(
        f"  row_quotient_representative={row_q}"
    )

    print(
        f"  dim(Col(F)/Col(A))="
        f"{len(F_cols)-len(A_cols)}"
    )

    print(
        f"  col_quotient_representative={col_q}"
    )

    # ------------------------------------------------------------------
    # 3. K DIRECTIONS
    # ------------------------------------------------------------------

    if len(K_rows) != 1:
        raise ArithmeticError(
            "K does not have rank-one row space."
        )

    if len(K_cols) != 1:
        raise ArithmeticError(
            "K does not have rank-one column space."
        )

    K_row_dir = K_rows[0]
    K_col_dir = K_cols[0]

    print()
    print("=" * 78)
    print("3. RANK-ONE K DIRECTIONS")
    print("=" * 78)

    print(
        f"  K_row_direction={K_row_dir}"
    )

    print(
        f"  K_col_direction={K_col_dir}"
    )

    # ------------------------------------------------------------------
    # 4. K LIES IN F BUT NOT A?
    # ------------------------------------------------------------------

    K_row_in_F = in_span_mod2(
        K_row_dir,
        F_rows,
    )

    K_row_in_A = in_span_mod2(
        K_row_dir,
        A_rows,
    )

    K_col_in_F = in_span_mod2(
        K_col_dir,
        F_cols,
    )

    K_col_in_A = in_span_mod2(
        K_col_dir,
        A_cols,
    )

    print()
    print("=" * 78)
    print("4. K VS ROW/COLUMN SPACES")
    print("=" * 78)

    print(
        f"  K_row_in_Row(F)={K_row_in_F}"
    )

    print(
        f"  K_row_in_Row(A)={K_row_in_A}"
    )

    print(
        f"  K_col_in_Col(F)={K_col_in_F}"
    )

    print(
        f"  K_col_in_Col(A)={K_col_in_A}"
    )

    # ------------------------------------------------------------------
    # 5. SAME QUOTIENT CLASS
    # ------------------------------------------------------------------

    row_sum = add2(
        row_q,
        K_row_dir,
    )

    col_sum = add2(
        col_q,
        K_col_dir,
    )

    row_same_class = in_span_mod2(
        row_sum,
        A_rows,
    )

    col_same_class = in_span_mod2(
        col_sum,
        A_cols,
    )

    print()
    print("=" * 78)
    print("5. QUOTIENT-CLASS EQUIVALENCE")
    print("=" * 78)

    print(
        f"  row_q XOR K_row={row_sum}"
    )

    print(
        f"  same_row_quotient_class="
        f"{row_same_class}"
    )

    print(
        f"  col_q XOR K_col={col_sum}"
    )

    print(
        f"  same_col_quotient_class="
        f"{col_same_class}"
    )

    # ------------------------------------------------------------------
    # 6. ALL NONZERO K ROWS / COLUMNS
    # ------------------------------------------------------------------

    nonzero_rows = [
        row
        for row in K2
        if any(row)
    ]

    nonzero_cols = [
        col
        for col in transpose(K2)
        if any(col)
    ]

    all_rows_same_class = True

    for row in nonzero_rows:

        if not in_span_mod2(
            add2(row, row_q),
            A_rows,
        ):
            all_rows_same_class = False
            break

    all_cols_same_class = True

    for col in nonzero_cols:

        if not in_span_mod2(
            add2(col, col_q),
            A_cols,
        ):
            all_cols_same_class = False
            break

    print()
    print("=" * 78)
    print("6. ALL NONZERO K ROWS / COLUMNS")
    print("=" * 78)

    print(
        f"  nonzero_K_rows={len(nonzero_rows)}"
    )

    print(
        f"  nonzero_K_columns={len(nonzero_cols)}"
    )

    print(
        f"  all_K_rows_same_quotient_class="
        f"{all_rows_same_class}"
    )

    print(
        f"  all_K_columns_same_quotient_class="
        f"{all_cols_same_class}"
    )

    # ------------------------------------------------------------------
    # 7. ADJUGATE REFERENCE
    # ------------------------------------------------------------------

    cof = cofactor_matrix(A)
    R = normalize_cofactor(cof)
    AdjR = transpose(R)
    AdjR2 = to_mod2(AdjR)

    adj_rows = row_basis_mod2(AdjR2)
    adj_cols = row_basis_mod2(
        transpose(AdjR2)
    )

    if len(adj_rows) != 1:
        raise ArithmeticError(
            "Normalized adjugate row space is not rank one."
        )

    if len(adj_cols) != 1:
        raise ArithmeticError(
            "Normalized adjugate column space is not rank one."
        )

    adj_row = adj_rows[0]
    adj_col = adj_cols[0]

    print()
    print("=" * 78)
    print("7. ADJUGATE REFERENCE")
    print("=" * 78)

    print(
        f"  adj_row={adj_row}"
    )

    print(
        f"  adj_col={adj_col}"
    )

    print(
        f"  adj_row_in_Row(A)="
        f"{in_span_mod2(adj_row, A_rows)}"
    )

    print(
        f"  adj_col_in_Col(A)="
        f"{in_span_mod2(adj_col, A_cols)}"
    )

    # ------------------------------------------------------------------
    # 8. INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The rank-one Schur correction K is now compared with the unique
one-dimensional quotient classes.

The quotient dimensions are

    dim(Row(F)/Row(A)) = 1,
    dim(Col(F)/Col(A)) = 1.

Because K has rank one, it has unique nonzero row and column directions
over F_2.

The canonical question is:

    [K_row] = [q_row] in Row(F)/Row(A)?

and

    [K_col] = [q_col] in Col(F)/Col(A)?

Over F_2 this is tested exactly by

    K_row XOR q_row ∈ Row(A),

and similarly for columns.

The experiment also checks every nonzero row and column of K, although
rank one already implies that all nonzero rows lie in the same
one-dimensional row space, and likewise for columns.

A positive result would establish:

    K is not merely rank one;

    K is the canonical representative of the Schur-lost quotient
    direction.

A negative result would mean the rank-one correction is responsible
for the rank change but is not canonically aligned with the quotient
representative chosen from Row(F) or Col(F).
"""
    )

    # ------------------------------------------------------------------
    # 9. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        rank_M == 4
        and rank_F == 4
        and rank_K == 1
        and rank_A == 3
        and len(F_rows) - len(A_rows) == 1
        and len(F_cols) - len(A_cols) == 1
        and K_row_in_F
        and K_col_in_F
        and not K_row_in_A
        and not K_col_in_A
        and row_same_class
        and col_same_class
        and all_rows_same_class
        and all_cols_same_class
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  rank_M_4={rank_M == 4}"
    )

    print(
        f"  rank_F_4={rank_F == 4}"
    )

    print(
        f"  rank_K_1={rank_K == 1}"
    )

    print(
        f"  rank_A_3={rank_A == 3}"
    )

    print(
        f"  row_quotient_one_dimensional="
        f"{len(F_rows)-len(A_rows) == 1}"
    )

    print(
        f"  col_quotient_one_dimensional="
        f"{len(F_cols)-len(A_cols) == 1}"
    )

    print(
        f"  K_row_in_F={K_row_in_F}"
    )

    print(
        f"  K_col_in_F={K_col_in_F}"
    )

    print(
        f"  K_row_not_in_A={not K_row_in_A}"
    )

    print(
        f"  K_col_not_in_A={not K_col_in_A}"
    )

    print(
        f"  row_quotient_class_match="
        f"{row_same_class}"
    )

    print(
        f"  col_quotient_class_match="
        f"{col_same_class}"
    )

    print(
        f"  every_K_row_matches_quotient="
        f"{all_rows_same_class}"
    )

    print(
        f"  every_K_col_matches_quotient="
        f"{all_cols_same_class}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 167R COMPLETE")


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