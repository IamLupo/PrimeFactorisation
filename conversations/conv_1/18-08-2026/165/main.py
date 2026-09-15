#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 165R — EXACT SCHUR-CORRECTION / MOD-2 RANK-DROP AUDIT
==============================================================================

Corrected replacement for Experiment 165.

Established data:

    rank_F2(M)     = 4
    rank_F2(F_ret) = 4
    rank_F2(A)     = 3

where

    M     = full 14x14 system,
    F_ret = retained 12x12 operator block before Schur elimination,
    A     = exact 12x12 Schur complement.

Therefore the interesting rank change is

    F_ret  ->  A,

not

    M -> F_ret.

The boundary block is

    H = [[1,1],[1,3]]

with

    det(H) = 2.

Hence the Schur complement is invertible over Q but not modulo 2.

This experiment computes the exact correction

    K = A - F_ret = -B H^{-1} C

and studies its mod-2 row/column geometry.

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

FULL_SIZE = 14
CORE_SIZE = 12

BOUNDARY_ROWS = (5, 10)
BOUNDARY_COLS = (12, 13)

DELTA_11 = 9512681472


# ============================================================================
# HELPERS
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
# BAREISS DETERMINANT
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
# FULL 14x14 SYSTEM
# ============================================================================

def build_full_system():

    M = []
    metadata = []

    for p, p_next in TRANSITIONS:

        for r in range(D[p] + 1):

            d = D[p] - r

            q0 = q_value(p, r)
            q1 = q_value(p, r + 1)

            b = basis(p, d)

            row = [0] * FULL_SIZE

            for j in range(6):
                row[j] = b[j] * q0

            for j in range(6):
                row[6 + j] = b[j] * q1

            if d == 0:
                row[12] = 1
                row[13] = p

            M.append(row)

            metadata.append(
                {
                    "p": p,
                    "p_next": p_next,
                    "r": r,
                    "d": d,
                }
            )

    if len(M) != FULL_SIZE:
        raise ArithmeticError(
            f"Expected {FULL_SIZE} rows, got {len(M)}."
        )

    assert_rectangular(M, "full system")

    return M, metadata


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

    a_rows = len(A)
    a_cols = len(A[0])

    b_rows = len(B)
    b_cols = len(B[0])

    if a_cols != b_rows:
        raise ValueError(
            f"Matrix multiply mismatch: "
            f"{a_rows}x{a_cols} times {b_rows}x{b_cols}."
        )

    out = [
        [
            Fraction(0)
            for _ in range(b_cols)
        ]
        for _ in range(a_rows)
    ]

    for i in range(a_rows):

        for k in range(a_cols):

            aik = Fraction(A[i][k])

            if aik == 0:
                continue

            for j in range(b_cols):
                out[i][j] += (
                    aik
                    * Fraction(B[k][j])
                )

    return out


def matrix_add(A, B):

    assert_rectangular(A, "A")
    assert_rectangular(B, "B")

    if (
        len(A) != len(B)
        or len(A[0]) != len(B[0])
    ):
        raise ValueError(
            "Matrix addition dimension mismatch."
        )

    return [
        [
            Fraction(A[i][j])
            + Fraction(B[i][j])
            for j in range(len(A[0]))
        ]
        for i in range(len(A))
    ]


def matrix_sub(A, B):

    assert_rectangular(A, "A")
    assert_rectangular(B, "B")

    if (
        len(A) != len(B)
        or len(A[0]) != len(B[0])
    ):
        raise ValueError(
            "Matrix subtraction dimension mismatch."
        )

    return [
        [
            Fraction(A[i][j])
            - Fraction(B[i][j])
            for j in range(len(A[0]))
        ]
        for i in range(len(A))
    ]


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
# EXACT SCHUR DATA
# ============================================================================

def build_schur_objects(M):

    boundary_equation_rows = [
        M[idx]
        for idx in BOUNDARY_ROWS
    ]

    retained_equation_rows = [
        M[idx]
        for idx in range(FULL_SIZE)
        if idx not in BOUNDARY_ROWS
    ]

    # B: retained equations x boundary variables, 12x2
    B = [
        [
            row[col]
            for col in BOUNDARY_COLS
        ]
        for row in retained_equation_rows
    ]

    # H: boundary equations x boundary variables, 2x2
    H = [
        [
            row[col]
            for col in BOUNDARY_COLS
        ]
        for row in boundary_equation_rows
    ]

    # C: boundary equations x operator variables, 2x12
    C = [
        [
            row[col]
            for col in range(CORE_SIZE)
        ]
        for row in boundary_equation_rows
    ]

    det_H = (
        H[0][0] * H[1][1]
        - H[0][1] * H[1][0]
    )

    if det_H == 0:
        raise ArithmeticError(
            "H is singular over Q."
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

    # F_ret: 12x12 retained operator block
    F_ret = [
        [
            row[col]
            for col in range(CORE_SIZE)
        ]
        for row in retained_equation_rows
    ]

    # A = F_ret - B H^{-1} C
    correction = matmul(
        B,
        matmul(
            H_inv,
            C,
        )
    )

    K_fraction = [
        [
            -x
            for x in row
        ]
        for row in correction
    ]

    K = matrix_to_integer(
        K_fraction,
        "Schur correction K",
    )

    A_fraction = matrix_add(
        F_ret,
        K,
    )

    A = matrix_to_integer(
        A_fraction,
        "Schur core A",
    )

    return {
        "retained_rows": retained_equation_rows,
        "B": B,
        "H": H,
        "C": C,
        "H_inv": H_inv,
        "F_ret": F_ret,
        "K": K,
        "A": A,
        "det_H": det_H,
    }


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

    A = to_mod2(A)

    if not A:
        return [], []

    rows = len(A)
    cols = len(A[0])

    pivots = []
    pivot_row = 0

    for col in range(cols):

        pivot = None

        for r in range(
            pivot_row,
            rows,
        ):

            if A[r][col]:
                pivot = r
                break

        if pivot is None:
            continue

        A[pivot_row], A[pivot] = (
            A[pivot],
            A[pivot_row],
        )

        for r in range(rows):

            if r == pivot_row:
                continue

            if A[r][col]:

                for j in range(cols):
                    A[r][j] ^= A[pivot_row][j]

        pivots.append(col)
        pivot_row += 1

        if pivot_row == rows:
            break

    return A, pivots


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

    dimension = len(v)

    if any(
        len(b) != dimension
        for b in basis
    ):
        raise ValueError(
            "Dimension mismatch in span test."
        )

    return (
        rank_mod2(basis)
        == rank_mod2(
            basis + [v]
        )
    )


def intersection_dimension(
    basis_A,
    rank_A,
    basis_B,
    rank_B,
    combined_rank,
):

    return (
        rank_A
        + rank_B
        - combined_rank
    )


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

    for i, row in enumerate(C):

        new_row = []

        for j, x in enumerate(row):

            if x % DELTA_11 != 0:
                raise ArithmeticError(
                    f"Delta_11 does not divide C[{i}][{j}]."
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
        "EXPERIMENT 165R — EXACT SCHUR-CORRECTION / "
        "MOD-2 RANK-DROP AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # BUILD
    # ------------------------------------------------------------------

    M, metadata = build_full_system()

    S = build_schur_objects(M)

    F_ret = S["F_ret"]
    K = S["K"]
    A = S["A"]
    B = S["B"]
    H = S["H"]
    C = S["C"]
    det_H = S["det_H"]

    # Exact consistency.
    exact_A_equals_F_plus_K = (
        A == matrix_to_integer(
            matrix_add(F_ret, K),
            "F_ret + K",
        )
    )

    # ------------------------------------------------------------------
    # 1. BASIC RANKS
    # ------------------------------------------------------------------

    full_rank = rank_mod2(M)
    F_rank = rank_mod2(F_ret)
    A_rank = rank_mod2(A)
    K_rank = rank_mod2(K)

    print()
    print("=" * 78)
    print("1. EXACT RANK DATA")
    print("=" * 78)

    print(
        f"  M_shape={len(M)}x{len(M[0])}"
    )

    print(
        f"  F_ret_shape={len(F_ret)}x{len(F_ret[0])}"
    )

    print(
        f"  A_shape={len(A)}x{len(A[0])}"
    )

    print(
        f"  K_shape={len(K)}x{len(K[0])}"
    )

    print(
        f"  rank_M_mod2={full_rank}"
    )

    print(
        f"  rank_F_ret_mod2={F_rank}"
    )

    print(
        f"  rank_A_mod2={A_rank}"
    )

    print(
        f"  rank_K_mod2={K_rank}"
    )

    print(
        f"  exact_A_equals_F_plus_K="
        f"{exact_A_equals_F_plus_K}"
    )

    print(
        f"  rank_drop_F_to_A="
        f"{F_rank - A_rank}"
    )

    # ------------------------------------------------------------------
    # 2. SCHUR BLOCK
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. EXACT SCHUR BLOCK")
    print("=" * 78)

    print(
        f"  H={H}"
    )

    print(
        f"  det_H={det_H}"
    )

    print(
        "  H_mod2="
        f"{to_mod2(H)}"
    )

    print(
        f"  B_shape={len(B)}x{len(B[0])}"
    )

    print(
        f"  C_shape={len(C)}x{len(C[0])}"
    )

    # ------------------------------------------------------------------
    # 3. K MOD 2
    # ------------------------------------------------------------------

    K2 = to_mod2(K)

    print()
    print("=" * 78)
    print("3. SCHUR CORRECTION K MOD 2")
    print("=" * 78)

    for row in K2:
        print(
            "  "
            + " ".join(
                str(x)
                for x in row
            )
        )

    # ------------------------------------------------------------------
    # 4. ROW-SPACE GEOMETRY
    # ------------------------------------------------------------------

    F_row = row_basis_mod2(F_ret)
    A_row = row_basis_mod2(A)
    K_row = row_basis_mod2(K)

    F_row_rank = len(F_row)
    A_row_rank = len(A_row)
    K_row_rank = len(K_row)

    combined_FA = rank_mod2(
        F_row + A_row
    )

    combined_FK = rank_mod2(
        F_row + K_row
    )

    combined_AK = rank_mod2(
        A_row + K_row
    )

    int_FA = intersection_dimension(
        F_row,
        F_row_rank,
        A_row,
        A_row_rank,
        combined_FA,
    )

    int_FK = intersection_dimension(
        F_row,
        F_row_rank,
        K_row,
        K_row_rank,
        combined_FK,
    )

    int_AK = intersection_dimension(
        A_row,
        A_row_rank,
        K_row,
        K_row_rank,
        combined_AK,
    )

    print()
    print("=" * 78)
    print("4. ROW-SPACE GEOMETRY")
    print("=" * 78)

    print(
        f"  dim(Row(F))={F_row_rank}"
    )

    print(
        f"  dim(Row(A))={A_row_rank}"
    )

    print(
        f"  dim(Row(K))={K_row_rank}"
    )

    print(
        f"  dim(Row(F)+Row(A))="
        f"{combined_FA}"
    )

    print(
        f"  dim(Row(F)∩Row(A))="
        f"{int_FA}"
    )

    print(
        f"  dim(Row(F)+Row(K))="
        f"{combined_FK}"
    )

    print(
        f"  dim(Row(F)∩Row(K))="
        f"{int_FK}"
    )

    print(
        f"  dim(Row(A)+Row(K))="
        f"{combined_AK}"
    )

    print(
        f"  dim(Row(A)∩Row(K))="
        f"{int_AK}"
    )

    # ------------------------------------------------------------------
    # 5. COLUMN-SPACE GEOMETRY
    # ------------------------------------------------------------------

    F_col = row_basis_mod2(
        transpose(
            to_mod2(F_ret)
        )
    )

    A_col = row_basis_mod2(
        transpose(
            to_mod2(A)
        )
    )

    K_col = row_basis_mod2(
        transpose(K2)
    )

    F_col_rank = len(F_col)
    A_col_rank = len(A_col)
    K_col_rank = len(K_col)

    combined_FC = rank_mod2(
        F_col + A_col
    )

    combined_FK_col = rank_mod2(
        F_col + K_col
    )

    combined_AK_col = rank_mod2(
        A_col + K_col
    )

    int_FC = (
        F_col_rank
        + A_col_rank
        - combined_FC
    )

    int_FK_col = (
        F_col_rank
        + K_col_rank
        - combined_FK_col
    )

    int_AK_col = (
        A_col_rank
        + K_col_rank
        - combined_AK_col
    )

    print()
    print("=" * 78)
    print("5. COLUMN-SPACE GEOMETRY")
    print("=" * 78)

    print(
        f"  dim(Col(F))={F_col_rank}"
    )

    print(
        f"  dim(Col(A))={A_col_rank}"
    )

    print(
        f"  dim(Col(K))={K_col_rank}"
    )

    print(
        f"  dim(Col(F)+Col(A))="
        f"{combined_FC}"
    )

    print(
        f"  dim(Col(F)∩Col(A))="
        f"{int_FC}"
    )

    print(
        f"  dim(Col(F)+Col(K))="
        f"{combined_FK_col}"
    )

    print(
        f"  dim(Col(F)∩Col(K))="
        f"{int_FK_col}"
    )

    print(
        f"  dim(Col(A)+Col(K))="
        f"{combined_AK_col}"
    )

    print(
        f"  dim(Col(A)∩Col(K))="
        f"{int_AK_col}"
    )

    # ------------------------------------------------------------------
    # 6. ADJUGATE
    # ------------------------------------------------------------------

    cof = cofactor_matrix(A)

    R = normalize_cofactor(cof)

    AdjR = transpose(R)
    AdjR2 = to_mod2(AdjR)

    Adj_row = row_basis_mod2(AdjR2)
    Adj_col = row_basis_mod2(
        transpose(AdjR2)
    )

    if len(Adj_row) != 1:
        raise ArithmeticError(
            "Normalized adjugate row space is not rank one."
        )

    if len(Adj_col) != 1:
        raise ArithmeticError(
            "Normalized adjugate column space is not rank one."
        )

    adj_row_dir = Adj_row[0]
    adj_col_dir = Adj_col[0]

    adj_row_in_A = in_span_mod2(
        adj_row_dir,
        A_row,
    )

    adj_row_in_F = in_span_mod2(
        adj_row_dir,
        F_row,
    )

    adj_col_in_A = in_span_mod2(
        adj_col_dir,
        A_col,
    )

    adj_col_in_F = in_span_mod2(
        adj_col_dir,
        F_col,
    )

    print()
    print("=" * 78)
    print("6. ADJUGATE DIRECTION")
    print("=" * 78)

    print(
        f"  adj_row_direction={adj_row_dir}"
    )

    print(
        f"  adj_col_direction={adj_col_dir}"
    )

    print(
        f"  adj_row_in_Row(A)={adj_row_in_A}"
    )

    print(
        f"  adj_row_in_Row(F)={adj_row_in_F}"
    )

    print(
        f"  adj_col_in_Col(A)={adj_col_in_A}"
    )

    print(
        f"  adj_col_in_Col(F)={adj_col_in_F}"
    )

    # ------------------------------------------------------------------
    # 7. FULL SYSTEM REFERENCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. FULL SYSTEM REFERENCE")
    print("=" * 78)

    print(
        f"  rank_M_mod2={full_rank}"
    )

    print(
        f"  rank_F_ret_mod2={F_rank}"
    )

    print(
        f"  rank_A_mod2={A_rank}"
    )

    print(
        f"  M_minus_F={full_rank - F_rank}"
    )

    print(
        f"  F_minus_A={F_rank - A_rank}"
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
The current evidence isolates the rank change more sharply than the
earlier experiments.

We have:

    rank_F2(M)     = 4
    rank_F2(F_ret) = 4
    rank_F2(A)     = 3.

Therefore the full system and the retained operator block already
have the same mod-2 rank.

The one-dimensional rank loss occurs specifically in the rational
Schur correction:

    A = F_ret + K,

where

    K = -B H^{-1} C

and

    det(H)=2.

Experiment 165 measures the relative positions of Row(F_ret),
Row(A), and Row(K), together with the corresponding column spaces.

The strongest possible simple outcome would be:

    dim(Row(F_ret) ∩ Row(A)) = 3,

so that the Schur core is a codimension-one subspace of the retained
row space.

A similarly clean column statement would identify the same loss on
the column side.

The rank of K alone is not enough: a correction can have rank larger
than one while changing the final row space by only one dimension.

The normalized adjugate line is checked only against these legitimate
12-dimensional spaces.

No causal interpretation is imposed.
"""
    )

    # ------------------------------------------------------------------
    # 9. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        F_rank == 4
        and A_rank == 3
        and full_rank == 4
        and exact_A_equals_F_plus_K
        and F_rank - A_rank == 1
        and len(Adj_row) == 1
        and len(Adj_col) == 1
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  F_ret_exact={F_rank == 4}"
    )

    print(
        f"  A_exact={A_rank == 3}"
    )

    print(
        f"  M_exact={full_rank == 4}"
    )

    print(
        f"  A_equals_F_plus_K_exact="
        f"{exact_A_equals_F_plus_K}"
    )

    print(
        f"  Schur_rank_drop_one="
        f"{F_rank - A_rank == 1}"
    )

    print(
        f"  adjugate_rank_one="
        f"{len(Adj_row) == 1 and len(Adj_col) == 1}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 165R COMPLETE")


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