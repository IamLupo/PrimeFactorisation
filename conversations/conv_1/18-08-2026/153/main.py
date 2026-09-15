#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 153R — EXACT 2-ADIC COFACTOR / F-CORE NULLSPACE AUDIT
==============================================================================

Corrected version of Experiment 153.

The previous script assumed that a particular pivot-based construction
would reproduce R mod 2. That assumption was too strong and caused an
unnecessary fatal error.

This version separates:

    rank(R mod 2)

from

    existence of the selected pivot reconstruction.

The exact audit proceeds safely:

    1. build the exact 12x12 F-core;
    2. build the normalized cofactor matrix R;
    3. compute rank(A mod 2);
    4. compute rank(R mod 2);
    5. recover row/column nullspaces of A mod 2;
    6. if R mod 2 has rank one, construct and independently verify
       an outer-product representation;
    7. otherwise report the actual rank without aborting;
    8. analyse the first 2-adic defect E.

No fitted recurrence.
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
# EXACT BAREISS DETERMINANT
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
# BUILD FULL SYSTEM
# ============================================================================

def build_system():

    M = []
    meta = []

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

            meta.append(
                {
                    "p": p,
                    "pnext": pnext,
                    "r": r,
                    "d": d,
                }
            )

    return M, meta


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

            mdet = bareiss_det(
                maximal_minor(A, i, j)
            )

            if mdet % DELTA_11 != 0:
                raise ArithmeticError(
                    "Delta_11 divisibility failure."
                )

            row.append(
                mdet // DELTA_11
            )

        R.append(row)

    return R


# ============================================================================
# MOD-2 LINEAR ALGEBRA
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

            if M[r][col]:

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


def transpose(A):

    return [
        [
            A[i][j]
            for i in range(len(A))
        ]
        for j in range(len(A[0]))
    ]


def row_space_basis_mod2(A):

    R, _ = rref_mod2(A)

    return [
        row
        for row in R
        if any(row)
    ]


def column_space_basis_mod2(A):

    return row_space_basis_mod2(
        transpose(A)
    )


# ============================================================================
# MOD-2 NULLSPACE
# ============================================================================

def nullspace_mod2(A):

    R, pivots = rref_mod2(A)

    ncols = len(A[0])
    pivot_set = set(pivots)

    free_cols = [
        j for j in range(ncols)
        if j not in pivot_set
    ]

    basis = []

    for free in free_cols:

        v = [0] * ncols
        v[free] = 1

        for row_index, pivot_col in enumerate(pivots):

            if R[row_index][free]:
                v[pivot_col] = 1

        basis.append(v)

    return basis, pivots


# ============================================================================
# VECTOR OPERATIONS MOD 2
# ============================================================================

def dot2(a, b):

    value = 0

    for x, y in zip(a, b):
        value ^= (
            (x & 1)
            & (y & 1)
        )

    return value


def mat_vec2(A, v):

    return [
        dot2(row, v)
        for row in A
    ]


def vector_in_span_mod2(v, basis):

    if not basis:
        return not any(v)

    before = rank_mod2(basis)
    after = rank_mod2(
        basis + [v]
    )

    return before == after


# ============================================================================
# SAFE RANK-ONE FACTORIZATION
# ============================================================================

def try_rank_one_factor_mod2(A):

    n = len(A)

    if rank_mod2(A) != 1:
        return {
            "rank": rank_mod2(A),
            "exists": False,
            "exact": False,
            "row": None,
            "column": None,
            "pivot": None,
        }

    pivot = None

    for i in range(n):

        for j in range(n):

            if A[i][j] & 1:
                pivot = (i, j)
                break

        if pivot is not None:
            break

    if pivot is None:
        return {
            "rank": 1,
            "exists": False,
            "exact": False,
            "row": None,
            "column": None,
            "pivot": None,
        }

    i0, j0 = pivot

    row = [
        A[i0][j] & 1
        for j in range(n)
    ]

    column = [
        A[i][j0] & 1
        for i in range(n)
    ]

    reconstruction = [
        [
            column[i] & row[j]
            for j in range(n)
        ]
        for i in range(n)
    ]

    exact = (
        reconstruction
        == [
            [
                A[i][j] & 1
                for j in range(n)
            ]
            for i in range(n)
        ]
    )

    return {
        "rank": 1,
        "exists": True,
        "exact": exact,
        "row": row,
        "column": column,
        "pivot": pivot,
    }


# ============================================================================
# FIRST 2-ADIC DEFECT
# ============================================================================

def make_defect(R):

    n = len(R)

    pivot = None
    inverse = None

    for i in range(n):

        for j in range(n):

            x = R[i][j] % 8

            if x == 0:
                continue

            try:
                inv = pow(x, -1, 8)
            except ValueError:
                continue

            pivot = (i, j)
            inverse = inv
            break

        if pivot is not None:
            break

    if pivot is None:
        raise ArithmeticError(
            "No unit pivot modulo 8."
        )

    i0, j0 = pivot

    u = [
        R[i][j0] % 8
        for i in range(n)
    ]

    v = [
        (
            R[i0][j] * inverse
        ) % 8
        for j in range(n)
    ]

    exact = True

    for i in range(n):

        for j in range(n):

            if (
                u[i] * v[j] % 8
                != R[i][j] % 8
            ):

                exact = False
                break

        if not exact:
            break

    if not exact:
        raise ArithmeticError(
            "Mod-8 factorization failed."
        )

    E = []

    for i in range(n):

        row = []

        for j in range(n):

            diff = (
                R[i][j]
                - u[i] * v[j]
            )

            if diff % 8 != 0:
                raise ArithmeticError(
                    "R - uv^T is not divisible by 8."
                )

            row.append(
                (diff // 8) & 1
            )

        E.append(row)

    return {
        "pivot": pivot,
        "u": u,
        "v": v,
        "E": E,
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 153R — EXACT 2-ADIC COFACTOR / "
        "F-CORE NULLSPACE AUDIT"
    )
    print("=" * 78)

    M, meta = build_system()
    A = build_core(M)
    R = normalized_cofactors(A)

    A2 = [
        [
            x & 1
            for x in row
        ]
        for row in A
    ]

    R2 = [
        [
            x & 1
            for x in row
        ]
        for row in R
    ]

    # ------------------------------------------------------------------
    # 1. BASIC VALIDATION
    # ------------------------------------------------------------------

    rank_A = rank_mod2(A2)
    rank_R = rank_mod2(R2)

    print()
    print("=" * 78)
    print("1. EXACT MOD-2 VALIDATION")
    print("=" * 78)

    print(
        f"  A_shape={len(A)}x{len(A[0])}"
    )

    print(
        f"  rank(A mod 2)={rank_A}"
    )

    print(
        f"  rank(R mod 2)={rank_R}"
    )

    # ------------------------------------------------------------------
    # 2. CORE NULLSPACES
    # ------------------------------------------------------------------

    right_null, right_pivots = nullspace_mod2(A2)
    left_null, left_pivots = nullspace_mod2(
        transpose(A2)
    )

    print()
    print("=" * 78)
    print("2. F-CORE NULLSPACES")
    print("=" * 78)

    print(
        f"  right_nullity={len(right_null)}"
    )

    print(
        f"  left_nullity={len(left_null)}"
    )

    print(
        f"  right_pivots={right_pivots}"
    )

    print(
        f"  left_pivots={left_pivots}"
    )

    for i, v in enumerate(right_null):
        print(
            f"  right_null_{i}={v}"
        )

    for i, v in enumerate(left_null):
        print(
            f"  left_null_{i}={v}"
        )

    # ------------------------------------------------------------------
    # 3. SAFE R-MOD-2 FACTORIZATION
    # ------------------------------------------------------------------

    factor = try_rank_one_factor_mod2(R2)

    print()
    print("=" * 78)
    print("3. NORMALIZED COFACTOR MOD-2 AUDIT")
    print("=" * 78)

    print(
        f"  rank_mod_2={factor['rank']}"
    )

    print(
        f"  rank_one_factorization_exists="
        f"{factor['exists']}"
    )

    print(
        f"  rank_one_factorization_exact="
        f"{factor['exact']}"
    )

    print(
        f"  pivot={factor['pivot']}"
    )

    if factor["exists"]:

        print(
            f"  row_direction={factor['row']}"
        )

        print(
            f"  column_direction={factor['column']}"
        )

    # ------------------------------------------------------------------
    # 4. NULLSPACE RELATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. NULLSPACE / COFACTOR-DIRECTION RELATION")
    print("=" * 78)

    if factor["exists"]:

        r_row = factor["row"]
        r_col = factor["column"]

        right_dots = [
            dot2(r_row, v)
            for v in right_null
        ]

        left_dots = [
            dot2(r_col, v)
            for v in left_null
        ]

        print(
            f"  row_direction_dot_right_null={right_dots}"
        )

        print(
            f"  column_direction_dot_left_null="
            f"{left_dots}"
        )

        print(
            "  row_direction_in_right_nullspan="
            f"{vector_in_span_mod2(r_row, right_null)}"
        )

        print(
            "  column_direction_in_left_nullspan="
            f"{vector_in_span_mod2(r_col, left_null)}"
        )

        core_rows = row_space_basis_mod2(A2)
        core_cols = column_space_basis_mod2(A2)

        print(
            "  row_direction_in_core_rowspace="
            f"{vector_in_span_mod2(r_row, core_rows)}"
        )

        print(
            "  column_direction_in_core_colspace="
            f"{vector_in_span_mod2(r_col, core_cols)}"
        )

    else:

        print(
            "  skipped_direction_tests="
            "R mod 2 is not rank one."
        )

    # ------------------------------------------------------------------
    # 5. FIRST 2-ADIC DEFECT
    # ------------------------------------------------------------------

    defect = make_defect(R)
    E = defect["E"]
    E2rank = rank_mod2(E)

    print()
    print("=" * 78)
    print("5. FIRST 2-ADIC DEFECT")
    print("=" * 78)

    print(
        f"  mod8_pivot={defect['pivot']}"
    )

    print(
        f"  mod8_factorization_exact=True"
    )

    print(
        f"  defect_rank_mod_2={E2rank}"
    )

    E_row_space = row_space_basis_mod2(E)
    E_col_space = column_space_basis_mod2(E)

    print(
        f"  defect_row_space_rank={len(E_row_space)}"
    )

    print(
        f"  defect_column_space_rank={len(E_col_space)}"
    )

    # ------------------------------------------------------------------
    # 6. CORE NULLSPACES VS DEFECT
    # ------------------------------------------------------------------

    right_images = [
        mat_vec2(E, v)
        for v in right_null
    ]

    left_images = [
        mat_vec2(transpose(E), v)
        for v in left_null
    ]

    right_defect_annihilates = all(
        not any(v)
        for v in right_images
    )

    left_defect_annihilates = all(
        not any(v)
        for v in left_images
    )

    print()
    print("=" * 78)
    print("6. CORE NULLSPACE / FIRST-DEFECT AUDIT")
    print("=" * 78)

    print(
        f"  all_right_core_nullvectors_killed_by_E="
        f"{right_defect_annihilates}"
    )

    print(
        f"  all_left_core_nullvectors_killed_by_E_transpose="
        f"{left_defect_annihilates}"
    )

    # ------------------------------------------------------------------
    # 7. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The corrected experiment separates two logically different tests.

First:

    what is the exact rank of R mod 2?

Second:

    if that rank is one, does its unique row/column direction
    have a recognizable relation to the nullspaces of the F-core A?

The first correction E is treated independently:

    R = u v^T + 8E.

Thus even if the normalized cofactor direction has no simple relation
to the mod-2 nullspaces of A, the first defect can still be tested
against those nullspaces.

A positive result would establish an exact linear-algebraic relation
between the F-core singular structure and the normalized cofactor
structure.

A negative result rules out the simplest nullspace explanation.

The experiment deliberately does not infer causation from rank alone.
"""
    )

    # ------------------------------------------------------------------
    # 8. FINAL
    # ------------------------------------------------------------------

    baseline_ok = (
        len(A) == 12
        and len(R) == 12
    )

    final_ok = (
        baseline_ok
        and rank_A == 3
        and E2rank == 3
    )

    print()
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  data_exact={baseline_ok}"
    )

    print(
        f"  core_rank_expected={rank_A == 3}"
    )

    print(
        f"  R_rank_computed={rank_R}"
    )

    print(
        f"  first_defect_rank_expected={E2rank == 3}"
    )

    print(
        f"  nullspace_audit_completed=True"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 153R COMPLETE")


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