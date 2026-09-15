#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 149 — EXACT 2-ADIC DEFECT BLOCK / SUPPORT / SUBSPACE AUDIT
==============================================================================

Experiment 148 established

    R = u v^T + 8E

with

    rank(E mod 2) = 3,

and an exact decomposition

    E = u_0 v_0^T + u_1 v_1^T + u_2 v_2^T

over F_2.

The 12 coordinates of the F-core are naturally ordered as

    columns 0..5   = F_0 coefficients,
    columns 6..11  = F_1 coefficients.

Experiment 149 asks whether the first 2-adic obstruction respects this
intrinsic block structure.

We test:

    1. E in the 4 block quadrants

           F0<-F0   F0<-F1
           F1<-F0   F1<-F1

    2. rank of each block over F_2;

    3. row/column support inside each block;

    4. whether the defect row-space is contained in, or decomposes
       relative to, the F0/F1 coordinate subspaces;

    5. whether the rank-3 defect splits naturally across the two
       coordinate blocks;

    6. whether deleting the distinguished boundary coordinates
       (the coordinates that appeared as sparse defect directions)
       lowers the defect rank.

No fitting is performed.

No new recurrence is searched for.

The goal is to determine whether the first 2-adic obstruction is aligned
with the already-existing algebraic decomposition of the operator.

No connection to the original (p,q)-kernel is asserted.
"""

from __future__ import annotations

from fractions import Fraction
import sys


# ============================================================================
# EXACT INPUT DATA
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
# HELPERS
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
# FULL SYSTEM
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
# NORMALIZED COFACTORS
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

            mdet = bareiss_det(
                maximal_minor(A, i, j)
            )

            if mdet % DELTA_11 != 0:
                raise ArithmeticError(
                    "Delta_11 divisibility failure."
                )

            row.append(mdet // DELTA_11)

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

    pivots = []
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

        pivots.append(col)
        pivot_row += 1

        if pivot_row == rows:
            break

    return M, pivots


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


# ============================================================================
# MOD-8 RANK-1 FACTORIZATION AND FIRST DEFECT
# ============================================================================

def rank1_factor_mod8(A):

    n = len(A)
    modulus = 8

    pivot_pos = None
    pivot_inverse = None

    for i in range(n):

        for j in range(n):

            value = A[i][j] % modulus

            if value == 0:
                continue

            try:
                inverse = pow(
                    value,
                    -1,
                    modulus,
                )
            except ValueError:
                continue

            pivot_pos = (i, j)
            pivot_inverse = inverse
            break

        if pivot_pos is not None:
            break

    if pivot_pos is None:
        return None

    i0, j0 = pivot_pos

    u = [
        A[i][j0] % modulus
        for i in range(n)
    ]

    v = [
        (
            A[i0][j] * pivot_inverse
        ) % modulus
        for j in range(n)
    ]

    exact = True

    for i in range(n):

        for j in range(n):

            if (
                u[i] * v[j] % modulus
                != A[i][j] % modulus
            ):

                exact = False
                break

        if not exact:
            break

    if not exact:
        return None

    E = []

    for i in range(n):

        row = []

        for j in range(n):

            value = (
                A[i][j]
                - u[i] * v[j]
            )

            if value % 8 != 0:
                raise ArithmeticError(
                    "Non-integral first defect."
                )

            row.append(
                (value // 8) & 1
            )

        E.append(row)

    return {
        "pivot": pivot_pos,
        "u": u,
        "v": v,
        "E": E,
    }


# ============================================================================
# BLOCK EXTRACTION
# ============================================================================

def submatrix(A, rows, cols):

    return [
        [
            A[i][j]
            for j in cols
        ]
        for i in rows
    ]


def matrix_support(A):

    rows = []
    cols = []

    for i, row in enumerate(A):
        if any(row):
            rows.append(i)

    for j in range(len(A[0])):
        if any(
            A[i][j]
            for i in range(len(A))
        ):
            cols.append(j)

    return rows, cols


def row_weights(A):

    return [
        sum(row)
        for row in A
    ]


def col_weights(A):

    return [
        sum(
            A[i][j]
            for i in range(len(A))
        )
        for j in range(len(A[0]))
    ]


# ============================================================================
# DELETION-RANK AUDIT
# ============================================================================

def deletion_rank(A, remove_rows=(), remove_cols=()):

    keep_rows = [
        i
        for i in range(len(A))
        if i not in set(remove_rows)
    ]

    keep_cols = [
        j
        for j in range(len(A[0]))
        if j not in set(remove_cols)
    ]

    B = submatrix(
        A,
        keep_rows,
        keep_cols,
    )

    if not B or not B[0]:
        return 0

    return rank_mod2(B)


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 149 — EXACT 2-ADIC DEFECT BLOCK / "
        "SUPPORT / SUBSPACE AUDIT"
    )
    print("=" * 78)

    M = build_full_system()
    core = build_core(M)
    R = normalized_cofactors(core)

    factor = rank1_factor_mod8(R)

    if factor is None:
        raise ArithmeticError(
            "Could not recover exact mod-8 rank-one factorization."
        )

    E = factor["E"]

    # ------------------------------------------------------------------
    # 1. DATA VALIDATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT DATA VALIDATION")
    print("=" * 78)

    print(
        "  R_shape=(12,12)"
    )

    print(
        f"  mod8_rank1_factorization=True"
    )

    print(
        f"  first_defect_rank={rank_mod2(E)}"
    )

    # ------------------------------------------------------------------
    # 2. NATURAL F0/F1 BLOCK DECOMPOSITION
    # ------------------------------------------------------------------

    F0 = list(range(0, 6))
    F1 = list(range(6, 12))

    blocks = {
        "F0<-F0": (F0, F0),
        "F0<-F1": (F0, F1),
        "F1<-F0": (F1, F0),
        "F1<-F1": (F1, F1),
    }

    print()
    print("=" * 78)
    print("2. F0/F1 BLOCK DECOMPOSITION OF FIRST DEFECT")
    print("=" * 78)

    block_ranks = {}

    for name, (rows, cols) in blocks.items():

        B = submatrix(
            E,
            rows,
            cols,
        )

        rank = rank_mod2(B)

        block_ranks[name] = rank

        print(
            f"  {name}: shape={len(rows)}x{len(cols)} "
            f"rank_mod_2={rank}"
        )

        print(
            f"    row_weights={row_weights(B)}"
        )

        print(
            f"    col_weights={col_weights(B)}"
        )

    # ------------------------------------------------------------------
    # 3. BLOCK SUPPORT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. BLOCK SUPPORT")
    print("=" * 78)

    for name, (rows, cols) in blocks.items():

        B = submatrix(
            E,
            rows,
            cols,
        )

        sr, sc = matrix_support(B)

        print(
            f"  {name}: nonzero_rows={sr}"
        )

        print(
            f"  {name}: nonzero_cols={sc}"
        )

    # ------------------------------------------------------------------
    # 4. GLOBAL DEFECT SUPPORT
    # ------------------------------------------------------------------

    support_rows, support_cols = matrix_support(E)

    print()
    print("=" * 78)
    print("4. GLOBAL DEFECT SUPPORT")
    print("=" * 78)

    print(
        f"  nonzero_rows={support_rows}"
    )

    print(
        f"  nonzero_columns={support_cols}"
    )

    print(
        f"  row_weights={row_weights(E)}"
    )

    print(
        f"  col_weights={col_weights(E)}"
    )

    # ------------------------------------------------------------------
    # 5. DISTINGUISHED-COORDINATE DELETION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. DISTINGUISHED-COORDINATE DELETION AUDIT")
    print("=" * 78)

    tests = {
        "remove_row_0": (
            (0,),
            (),
        ),
        "remove_row_5": (
            (5,),
            (),
        ),
        "remove_rows_0_5": (
            (0, 5),
            (),
        ),
        "remove_col_0": (
            (),
            (0,),
        ),
        "remove_col_3": (
            (),
            (3,),
        ),
        "remove_cols_0_3": (
            (),
            (0, 3),
        ),
        "remove_row5_col0": (
            (5,),
            (0,),
        ),
        "remove_rows0_5_cols0_3": (
            (0, 5),
            (0, 3),
        ),
    }

    deletion_results = {}

    for name, (rr, cc) in tests.items():

        rank = deletion_rank(
            E,
            rr,
            cc,
        )

        deletion_results[name] = rank

        print(
            f"  {name}: rank={rank}"
        )

    # ------------------------------------------------------------------
    # 6. BLOCK-COMPONENT RANKS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. F0/F1 COMPONENT RANK TEST")
    print("=" * 78)

    for component_name in (
        "F0_ROWS_ONLY",
        "F1_ROWS_ONLY",
        "F0_COLS_ONLY",
        "F1_COLS_ONLY",
    ):

        if component_name == "F0_ROWS_ONLY":
            rank = deletion_rank(
                E,
                F1,
                (),
            )

        elif component_name == "F1_ROWS_ONLY":
            rank = deletion_rank(
                E,
                F0,
                (),
            )

        elif component_name == "F0_COLS_ONLY":
            rank = deletion_rank(
                E,
                (),
                F1,
            )

        else:
            rank = deletion_rank(
                E,
                (),
                F0,
            )

        print(
            f"  {component_name}: rank={rank}"
        )

    # ------------------------------------------------------------------
    # 7. FIRST-DEFECT INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The first 2-adic defect is a 12x12 matrix attached to the F-core.

The natural coordinate split is

    F_0 : coordinates 0..5,
    F_1 : coordinates 6..11.

This experiment does not alter that coordinate system.

Instead it asks whether the exact rank-3 obstruction

    E mod 2

is concentrated in one block, split across blocks, or only appears
after mixing the two blocks.

The deletion audit is deliberately small:

    if removing one or two distinguished coordinates lowers the
    defect rank, then those coordinates are structurally involved;

    if no such deletion changes the rank, the obstruction is more
    genuinely distributed.

A low-rank individual block is not by itself evidence of a deeper
recurrence. The relevant observation is the relationship between
the four blocks and the global rank-3 space.

No source-level interpretation is assumed.
"""
    )

    # ------------------------------------------------------------------
    # 8. FINAL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    final_ok = (
        len(R) == 12
        and len(R[0]) == 12
        and rank_mod2(E) == 3
        and all(
            name in block_ranks
            for name in blocks
        )
        and len(deletion_results) == len(tests)
    )

    print(
        f"  data_exact={len(R) == 12 and len(R[0]) == 12}"
    )

    print(
        f"  first_defect_rank_exact={rank_mod2(E) == 3}"
    )

    print(
        f"  block_audit_completed={len(block_ranks) == 4}"
    )

    print(
        f"  deletion_audit_completed="
        f"{len(deletion_results) == len(tests)}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 149 COMPLETE")


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
