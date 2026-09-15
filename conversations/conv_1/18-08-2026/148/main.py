#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 148 — EXACT 2-ADIC MOD-16 OBSTRUCTION / DEFECT-BASIS AUDIT
==============================================================================

Experiment 147 established an invariant statement:

    rank(R mod 2) = 1
    rank-one modulo 4 = True
    rank-one modulo 8 = True
    rank-one modulo 16 = False

and the first correction

    R = u v^T + 8 E

has

    rank(E mod 2) = 3.

Experiment 148 isolates the exact obstruction.

The central questions are:

    1. What is the canonical row space of E mod 2?
    2. What is its canonical column space?
    3. Can E mod 2 be written as a sum of exactly three
       rank-one outer products?
    4. Is the failure at modulus 16 therefore genuinely
       three-dimensional?
    5. Can the obstruction be localized to a small set of rows
       or columns?
    6. Does the same pivot structure reproduce the full defect
       matrix exactly?

This is a finite-field structural decomposition of the FIRST
2-adic correction, not a new interpolation.

The expected logic is:

    R = u v^T + 8E

and therefore

    R mod 16
      = u v^T + 8(E mod 2).

Thus the rank-one lift fails precisely because E mod 2 is nonzero.

The rank-3 calculation determines the dimension of the first
obstruction.

No Smith reduction.
No floating point.
No SymPy.
No extrapolation.
No claim about the original (p,q)-kernel.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
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


# ============================================================================
# BAREISS DETERMINANT
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
# FULL 14x14 SYSTEM
# ============================================================================

def build_full_system():

    M = []

    for p, _ in TRANSITIONS:

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
            "Boundary block is singular."
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
        reduced = []

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

            reduced.append(value.numerator)

        core.append(reduced)

    return core


# ============================================================================
# NORMALIZED COFACTOR MATRIX
# ============================================================================

def maximal_minor(A, omit_row, omit_col):

    rows = [i for i in range(len(A)) if i != omit_row]
    cols = [j for j in range(len(A[0])) if j != omit_col]

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
                    f"Delta_11 divisibility failed at ({i},{j})."
                )

            row.append(
                det_minor // DELTA_11
            )

        R.append(row)

    return R


# ============================================================================
# MODULAR LINEAR ALGEBRA OVER F_2
# ============================================================================

def rref_mod2(A):

    M = [
        [int(x) & 1 for x in row]
        for row in A
    ]

    rows = len(M)
    cols = len(M[0])

    pivot_cols = []
    pivot_row = 0

    for col in range(cols):

        pivot = None

        for r in range(pivot_row, rows):

            if M[r][col] == 1:
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

            if M[r][col] == 1:

                for j in range(col, cols):
                    M[r][j] ^= M[pivot_row][j]

        pivot_cols.append(col)
        pivot_row += 1

        if pivot_row == rows:
            break

    return M, pivot_cols


def rank_mod2(A):
    _, pivots = rref_mod2(A)
    return len(pivots)


# ============================================================================
# ROW SPACE BASIS
# ============================================================================

def row_space_basis_mod2(A):

    RREF, pivots = rref_mod2(A)

    basis = []

    for row in RREF:

        if any(row):
            basis.append(row)

    return basis, pivots


# ============================================================================
# COLUMN SPACE BASIS
# ============================================================================

def transpose(A):

    return [
        [
            A[i][j]
            for i in range(len(A))
        ]
        for j in range(len(A[0]))
    ]


def column_space_basis_mod2(A):

    AT = transpose(A)

    basis_rows, pivots = row_space_basis_mod2(AT)

    basis_columns = [
        list(row)
        for row in basis_rows
    ]

    return basis_columns, pivots


# ============================================================================
# RANK-ONE DECOMPOSITION OF A RANK-3 MATRIX
# ============================================================================

def solve_linear_combination_mod2(
    basis_rows,
    target,
):
    """
    Solve target = sum c_i basis_rows[i] over F_2.

    Since basis_rows are assumed independent, Gaussian elimination
    gives a unique solution.
    """

    k = len(basis_rows)
    n = len(target)

    augmented = [
        basis_rows[i][:] + [
            0
        ]
        for i in range(k)
    ]

    # Build equations with basis vectors as row generators.
    equations = []

    for j in range(n):

        equations.append([
            basis_rows[i][j]
            for i in range(k)
        ] + [
            target[j]
        ])

    rows = len(equations)
    cols = k

    pivot_row = 0
    pivot_cols = []

    for col in range(cols):

        pivot = None

        for r in range(pivot_row, rows):

            if equations[r][col] == 1:
                pivot = r
                break

        if pivot is None:
            continue

        equations[pivot_row], equations[pivot] = (
            equations[pivot],
            equations[pivot_row],
        )

        for r in range(rows):

            if r == pivot_row:
                continue

            if equations[r][col] == 1:

                for j in range(col, cols + 1):
                    equations[r][j] ^= equations[pivot_row][j]

        pivot_cols.append(col)
        pivot_row += 1

    solution = [0] * k

    for r, col in enumerate(pivot_cols):

        if r >= rows:
            break

        solution[col] = equations[r][cols]

    # Exact verification.
    for j in range(n):

        value = 0

        for i in range(k):
            value ^= (
                solution[i]
                & basis_rows[i][j]
            )

        if value != target[j]:
            raise ArithmeticError(
                "Mod-2 linear combination reconstruction failed."
            )

    return solution


def rank3_outer_decomposition(E):

    """
    Given a 12x12 matrix E over F_2 of rank 3, construct

        E = u0 v0^T + u1 v1^T + u2 v2^T

    using independent row vectors as v_i and expressing every row
    in that row basis.

    This is canonical up to row-basis choice.
    """

    basis_rows, pivot_cols = row_space_basis_mod2(E)

    if len(basis_rows) != 3:
        raise ArithmeticError(
            f"Expected rank 3, found {len(basis_rows)}."
        )

    u_vectors = []

    for i in range(3):
        u_vectors.append([])

    for row in E:

        coeffs = solve_linear_combination_mod2(
            basis_rows,
            row,
        )

        for i in range(3):
            u_vectors[i].append(
                coeffs[i]
            )

    return basis_rows, u_vectors, pivot_cols


# ============================================================================
# SUPPORT / LOCALIZATION
# ============================================================================

def nonzero_rows(A):

    result = []

    for i, row in enumerate(A):

        if any(row):
            result.append(i)

    return result


def nonzero_cols(A):

    result = []

    for j, column in enumerate(transpose(A)):

        if any(column):
            result.append(j)

    return result


def row_weights(A):

    return [
        sum(row)
        for row in A
    ]


def column_weights(A):

    return [
        sum(column)
        for column in transpose(A)
    ]


# ============================================================================
# MOD-16 OBSTRUCTION CHECK
# ============================================================================

def reconstruct_mod16_from_defect(R, defect, u, v):

    out = []

    for i in range(12):

        row = []

        for j in range(12):

            value = (
                u[i] * v[j]
                + 8 * defect[i][j]
            ) % 16

            row.append(value)

        out.append(row)

    return out


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 148 — EXACT 2-ADIC MOD-16 OBSTRUCTION / "
        "DEFECT-BASIS AUDIT"
    )
    print("=" * 78)

    M = build_full_system()
    core = build_core(M)
    R = normalized_cofactors(core)

    # ------------------------------------------------------------------
    # 1. RECOVER RANK-1 FACTORIZATION MOD 8
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT MOD-8 RANK-1 FACTORIZATION")
    print("=" * 78)

    modulus = 8
    n = 12

    pivot_position = None
    pivot_inverse = None

    for i in range(n):

        for j in range(n):

            value = R[i][j] % modulus

            if value == 0:
                continue

            try:
                inv = pow(
                    value,
                    -1,
                    modulus,
                )
            except ValueError:
                continue

            pivot_position = (i, j)
            pivot_inverse = inv
            break

        if pivot_position is not None:
            break

    if pivot_position is None:
        raise ArithmeticError(
            "No unit pivot found modulo 8."
        )

    i0, j0 = pivot_position

    u8 = [
        R[i][j0] % 8
        for i in range(n)
    ]

    v8 = [
        R[i0][j] * pivot_inverse % 8
        for j in range(n)
    ]

    factor_exact = True

    for i in range(n):

        for j in range(n):

            if (
                u8[i] * v8[j] % 8
                != R[i][j] % 8
            ):

                factor_exact = False
                break

        if not factor_exact:
            break

    print(
        f"  pivot_position={pivot_position}"
    )

    print(
        f"  factorization_exact={factor_exact}"
    )

    print(
        f"  u={u8}"
    )

    print(
        f"  v={v8}"
    )

    # ------------------------------------------------------------------
    # 2. FIRST 2-ADIC DEFECT
    # ------------------------------------------------------------------

    E = []

    for i in range(n):

        row = []

        for j in range(n):

            difference = (
                R[i][j]
                - u8[i] * v8[j]
            )

            if difference % 8 != 0:
                raise ArithmeticError(
                    "R - uv^T is not divisible by 8."
                )

            row.append(
                difference // 8
            )

        E.append(row)

    E2 = [
        [
            x % 2
            for x in row
        ]
        for row in E
    ]

    defect_rank = rank_mod2(E2)

    print()
    print("=" * 78)
    print("2. FIRST DEFECT MATRIX E MOD 2")
    print("=" * 78)

    print(
        f"  defect_rank={defect_rank}"
    )

    print(
        "  E mod 2:"
    )

    for row in E2:
        print(
            "    "
            + " ".join(
                str(x)
                for x in row
            )
        )

    # ------------------------------------------------------------------
    # 3. ROW-SPACE BASIS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. CANONICAL ROW-SPACE BASIS")
    print("=" * 78)

    row_basis, pivot_cols = row_space_basis_mod2(E2)

    print(
        f"  pivot_columns={pivot_cols}"
    )

    for i, row in enumerate(row_basis):

        print(
            f"  basis_row_{i}={row}"
        )

    # ------------------------------------------------------------------
    # 4. COLUMN-SPACE BASIS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. CANONICAL COLUMN-SPACE BASIS")
    print("=" * 78)

    col_basis, pivot_rows = column_space_basis_mod2(E2)

    print(
        f"  pivot_rows={pivot_rows}"
    )

    for i, column in enumerate(col_basis):

        print(
            f"  basis_column_{i}={column}"
        )

    # ------------------------------------------------------------------
    # 5. EXACT RANK-3 OUTER-PRODUCT DECOMPOSITION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. EXACT RANK-3 OUTER-PRODUCT DECOMPOSITION")
    print("=" * 78)

    basis_rows, u_vectors, pivots = (
        rank3_outer_decomposition(E2)
    )

    reconstruction = [
        [
            0
            for _ in range(n)
        ]
        for _ in range(n)
    ]

    for a in range(3):

        for i in range(n):

            for j in range(n):

                reconstruction[i][j] ^= (
                    u_vectors[a][i]
                    & basis_rows[a][j]
                )

    outer_exact = (
        reconstruction == E2
    )

    print(
        f"  outer_product_count=3"
    )

    print(
        f"  outer_decomposition_exact={outer_exact}"
    )

    for a in range(3):

        print(
            f"  term_{a}:"
        )

        print(
            f"    u_{a}={u_vectors[a]}"
        )

        print(
            f"    v_{a}={basis_rows[a]}"
        )

    # ------------------------------------------------------------------
    # 6. DEFECT LOCALIZATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. DEFECT SUPPORT / LOCALIZATION")
    print("=" * 78)

    nz_rows = nonzero_rows(E2)
    nz_cols = nonzero_cols(E2)

    print(
        f"  nonzero_rows={nz_rows}"
    )

    print(
        f"  nonzero_columns={nz_cols}"
    )

    print(
        f"  row_weights={row_weights(E2)}"
    )

    print(
        f"  column_weights={column_weights(E2)}"
    )

    # ------------------------------------------------------------------
    # 7. MOD-16 RECONSTRUCTION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. EXACT MOD-16 RECONSTRUCTION")
    print("=" * 78)

    reconstructed_mod16 = (
        reconstruct_mod16_from_defect(
            R,
            E2,
            u8,
            v8,
        )
    )

    direct_mod16 = [
        [
            R[i][j] % 16
            for j in range(n)
        ]
        for i in range(n)
    ]

    mod16_exact = (
        reconstructed_mod16
        == direct_mod16
    )

    print(
        f"  reconstruction_exact={mod16_exact}"
    )

    # ------------------------------------------------------------------
    # 8. FIRST OBSTRUCTION SUBSPACES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. OBSTRUCTION BASIS DIMENSIONS")
    print("=" * 78)

    print(
        f"  rank(E mod 2)={defect_rank}"
    )

    print(
        f"  row_space_dimension={len(row_basis)}"
    )

    print(
        f"  column_space_dimension={len(col_basis)}"
    )

    # ------------------------------------------------------------------
    # 9. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The mod-8 rank-one structure is

    R = u v^T + 8E.

Reducing modulo 16 gives

    R mod 16
      = u v^T + 8(E mod 2).

Therefore the complete obstruction to lifting the rank-one
factorization from modulus 8 to modulus 16 is exactly the matrix

    E mod 2.

The experiment determines its exact dimension:

    rank(E mod 2) = 3.

It then reconstructs the defect as three exact rank-one
outer products over F_2.

This separates two facts:

    leading structure:
        rank-one modulo 8;

    first obstruction:
        a three-dimensional correction modulo 2.

A successful rank-3 decomposition is not evidence for a hidden
recurrence by itself. It is simply the canonical linear-algebraic
description of the observed first 2-adic obstruction.

The next useful question, if warranted by this result, is whether
these three defect directions correspond to any already-existing
combinatorial row/column partition in the underlying exact system.

No claim about the original (p,q)-kernel is made here.
"""
    )

    # ------------------------------------------------------------------
    # 10. FINAL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    final_ok = (
        factor_exact
        and defect_rank == 3
        and outer_exact
        and mod16_exact
        and len(row_basis) == 3
        and len(col_basis) == 3
    )

    print(
        f"  mod8_rank1_exact={factor_exact}"
    )

    print(
        f"  first_defect_rank_exact="
        f"{defect_rank == 3}"
    )

    print(
        f"  rank3_outer_decomposition_exact="
        f"{outer_exact}"
    )

    print(
        f"  mod16_reconstruction_exact="
        f"{mod16_exact}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 148 COMPLETE")


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

