#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 164U — EXACT MOD-2 BOUNDARY-RANK / SCHUR-EXTENSION AUDIT
==============================================================================

HARDENED DIMENSION-CORRECT VERSION

The previous 164 variants repeatedly mixed different ambient vector spaces.

This version treats the spaces explicitly.

FULL SYSTEM
    M : 14 x 14

RETAINED EQUATIONS
    M_ret : 12 x 14

F-ONLY PROJECTION
    F_ret : 12 x 12
    F_all : 14 x 12

CORE
    A : 12 x 12

IMPORTANT:

    Row spaces of M_ret and M live in different ambient spaces if M_ret
    is kept as a 12-row matrix.

Therefore, before comparing their column spaces, M_ret is embedded into
the 14-row ambient space by inserting zero rows at the two removed
boundary-equation positions.

Likewise, all quotient comparisons explicitly validate dimensions.

The goal is to isolate the actual mod-2 rank extension:

    rank_F2(M) = 4
    rank_F2(A) = 3

and determine whether that +1 extension is attributable to the boundary
equations/columns in a dimensionally valid way.

No SymPy.
No floating point.
No recurrence fitting.
No extrapolation.
No connection to the original (p,q)-kernel.
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
                f"{name} is not rectangular at row {i}: "
                f"expected width {width}, got {len(row)}."
            )


def ambient_dimension(A):
    assert_rectangular(A)
    if not A:
        return (0, 0)
    return (len(A), len(A[0]))


# ============================================================================
# BAREISS
# ============================================================================

def bareiss_det(A):

    if not A:
        return 1

    assert_rectangular(A, "Bareiss input")

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
# BUILD FULL SYSTEM
# ============================================================================

def build_full_system():

    M = []
    metadata = []

    for p, pnext in TRANSITIONS:

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

            metadata.append(
                {
                    "p": p,
                    "pnext": pnext,
                    "r": r,
                    "d": d,
                }
            )

    if len(M) != FULL_N:
        raise ArithmeticError(
            f"Expected {FULL_N} full rows, got {len(M)}."
        )

    assert_rectangular(M, "M")

    return M, metadata


# ============================================================================
# SCHUR COMPLEMENT
# ============================================================================

def build_core(M):

    b0 = M[BOUNDARY_ROWS[0]]
    b1 = M[BOUNDARY_ROWS[1]]

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
            "Boundary block H is singular over Q."
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
        [Fraction(0) for _ in range(CORE_N)]
        for _ in range(2)
    ]

    for a in range(2):
        for j in range(CORE_N):

            h_coeff[a][j] = -(
                H_inv[a][0] * b0[j]
                + H_inv[a][1] * b1[j]
            )

    retained_rows = [
        i
        for i in range(FULL_N)
        if i not in BOUNDARY_ROWS
    ]

    core = []

    for i in retained_rows:

        row = M[i]
        out = []

        for j in range(CORE_N):

            value = (
                Fraction(row[j])
                + Fraction(row[12]) * h_coeff[0][j]
                + Fraction(row[13]) * h_coeff[1][j]
            )

            if value.denominator != 1:
                raise ArithmeticError(
                    "Schur complement entry is non-integral."
                )

            out.append(value.numerator)

        core.append(out)

    assert_rectangular(core, "core")

    return core, H, h_coeff, retained_rows


# ============================================================================
# MOD-2 LINEAR ALGEBRA
# ============================================================================

def to_mod2(A):
    assert_rectangular(A, "to_mod2 input")

    return [
        [int(x) & 1 for x in row]
        for row in A
    ]


def transpose(A):
    assert_rectangular(A, "transpose input")

    if not A:
        return []

    return [
        [A[i][j] for i in range(len(A))]
        for j in range(len(A[0]))
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


def row_space_basis_mod2(A):

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

    result = []

    for free in free_cols:

        v = [0] * n
        v[free] = 1

        for i, pivot_col in enumerate(pivots):

            if R[i][free]:
                v[pivot_col] = 1

        result.append(v)

    return result


def in_span_mod2(v, basis_vectors):

    if not basis_vectors:
        return not any(v)

    dim = len(v)

    for b in basis_vectors:
        if len(b) != dim:
            raise ValueError(
                "Dimension mismatch in span test: "
                f"vector dimension={dim}, "
                f"basis dimension={len(b)}."
            )

    before = rank_mod2(basis_vectors)
    after = rank_mod2(
        basis_vectors + [v]
    )

    return before == after


def quotient_witnesses(subspace_basis,
                       full_basis):

    if subspace_basis:

        dimension = len(subspace_basis[0])

        for v in subspace_basis:
            if len(v) != dimension:
                raise ValueError(
                    "Nonuniform subspace basis dimensions."
                )

        for v in full_basis:
            if len(v) != dimension:
                raise ValueError(
                    "Quotient comparison dimension mismatch."
                )

    elif full_basis:

        dimension = len(full_basis[0])

        for v in full_basis:
            if len(v) != dimension:
                raise ValueError(
                    "Nonuniform full basis dimensions."
                )

    current = list(subspace_basis)
    current_rank = rank_mod2(current)

    witnesses = []

    for v in full_basis:

        trial = current + [v]
        new_rank = rank_mod2(trial)

        if new_rank > current_rank:

            witnesses.append(v)
            current.append(v)
            current_rank = new_rank

    return witnesses


# ============================================================================
# COFACTOR / ADJUGATE
# ============================================================================

def maximal_minor(A, omit_row, omit_col):

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


def normalize(C):

    R = []

    for row in C:

        out = []

        for value in row:

            if value % DELTA_11 != 0:
                raise ArithmeticError(
                    "Delta_11 divisibility failure."
                )

            out.append(
                value // DELTA_11
            )

        R.append(out)

    return R


def get_adjugate_line(R):

    AdjR = transpose(R)
    AdjR2 = to_mod2(AdjR)

    row_basis = row_space_basis_mod2(AdjR2)
    col_basis = row_space_basis_mod2(
        transpose(AdjR2)
    )

    if len(row_basis) != 1:
        raise ArithmeticError(
            "Expected rank-one adjugate row space."
        )

    if len(col_basis) != 1:
        raise ArithmeticError(
            "Expected rank-one adjugate column space."
        )

    return (
        AdjR,
        AdjR2,
        row_basis[0],
        col_basis[0],
    )


# ============================================================================
# EMBEDDING HELPERS
# ============================================================================

def embed_retained_rows_zero(retained_matrix,
                             retained_indices,
                             total_rows):

    """
    Embed a 12 x n retained-equation matrix into the 14-row ambient space
    by placing zero rows at the removed boundary equations.

    Result: 14 x n.
    """

    assert_rectangular(
        retained_matrix,
        "retained_matrix"
    )

    if len(retained_matrix) != len(
        retained_indices
    ):
        raise ValueError(
            "Retained-row count mismatch."
        )

    if any(
        i < 0 or i >= total_rows
        for i in retained_indices
    ):
        raise ValueError(
            "Retained-row index out of range."
        )

    width = (
        len(retained_matrix[0])
        if retained_matrix
        else 0
    )

    result = [
        [0] * width
        for _ in range(total_rows)
    ]

    for local_i, global_i in enumerate(
        retained_indices
    ):

        result[global_i] = list(
            retained_matrix[local_i]
        )

    return result


def project_full_left_kernel(
    full_left_kernel,
    retained_indices,
):

    """
    Full left-kernel vectors have 14 coordinates, one per equation.

    Projection removes the two boundary-equation coordinates and therefore
    produces 12-dimensional vectors.
    """

    return [
        [
            v[i]
            for i in retained_indices
        ]
        for v in full_left_kernel
    ]


def project_full_right_kernel(
    full_right_kernel,
):

    """
    Full right-kernel vectors have 14 variable coordinates.

    Projection removes the two boundary-variable coordinates.
    """

    return [
        [
            v[j]
            for j in range(FULL_N)
            if j not in BOUNDARY_COLS
        ]
        for v in full_right_kernel
    ]


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 164U — EXACT MOD-2 BOUNDARY-RANK / "
        "SCHUR-EXTENSION AUDIT"
    )
    print("=" * 78)

    M, metadata = build_full_system()

    A, H, h_coeff, retained_indices = build_core(M)

    C = cofactor_matrix(A)
    R = normalize(C)

    AdjR, AdjR2, adj_row_dir, adj_col_dir = (
        get_adjugate_line(R)
    )

    M2 = to_mod2(M)
    A2 = to_mod2(A)

    # ------------------------------------------------------------------
    # 1. BASE SYSTEMS
    # ------------------------------------------------------------------

    retained_rows_full = [
        M[i]
        for i in retained_indices
    ]

    F_ret = [
        M[i][:CORE_N]
        for i in retained_indices
    ]

    F_all = [
        M[i][:CORE_N]
        for i in range(FULL_N)
    ]

    H_cols = [
        [
            M[i][j]
            for j in BOUNDARY_COLS
        ]
        for i in range(FULL_N)
    ]

    boundary_rows = [
        M[i]
        for i in BOUNDARY_ROWS
    ]

    full_rank = rank_mod2(M2)
    core_rank = rank_mod2(A2)
    f_ret_rank = rank_mod2(F_ret)
    f_all_rank = rank_mod2(F_all)
    h_col_rank = rank_mod2(H_cols)
    boundary_row_rank = rank_mod2(
        boundary_rows
    )

    print()
    print("=" * 78)
    print("1. EXACT MOD-2 RANKS")
    print("=" * 78)

    print(
        f"  M_shape=14x14"
    )

    print(
        f"  full_rank_mod2={full_rank}"
    )

    print(
        f"  A_shape=12x12"
    )

    print(
        f"  core_rank_mod2={core_rank}"
    )

    print(
        f"  F_ret_shape=12x12"
    )

    print(
        f"  F_ret_rank_mod2={f_ret_rank}"
    )

    print(
        f"  F_all_shape=14x12"
    )

    print(
        f"  F_all_rank_mod2={f_all_rank}"
    )

    print(
        f"  boundary_column_rank_mod2={h_col_rank}"
    )

    print(
        f"  boundary_row_rank_mod2={boundary_row_rank}"
    )

    # ------------------------------------------------------------------
    # 2. F-BLOCK QUOTIENT
    # ------------------------------------------------------------------

    f_ret_basis = row_space_basis_mod2(F_ret)
    f_all_basis = row_space_basis_mod2(F_all)

    f_quotient = quotient_witnesses(
        f_ret_basis,
        f_all_basis,
    )

    print()
    print("=" * 78)
    print("2. F-BLOCK QUOTIENT")
    print("=" * 78)

    print(
        f"  retained_F_rowspace_rank="
        f"{len(f_ret_basis)}"
    )

    print(
        f"  all_F_rowspace_rank="
        f"{len(f_all_basis)}"
    )

    print(
        f"  quotient_dimension="
        f"{len(f_quotient)}"
    )

    for i, w in enumerate(f_quotient):
        print(
            f"  quotient_witness_{i}={w}"
        )

    # ------------------------------------------------------------------
    # 3. FULL ROW-SPACE QUOTIENT
    # ------------------------------------------------------------------

    full_row_basis = row_space_basis_mod2(
        M2
    )

    retained_row_basis = row_space_basis_mod2(
        to_mod2(retained_rows_full)
    )

    full_row_quotient = quotient_witnesses(
        retained_row_basis,
        full_row_basis,
    )

    print()
    print("=" * 78)
    print("3. FULL ROW-SPACE QUOTIENT")
    print("=" * 78)

    print(
        f"  retained_rowspace_rank="
        f"{len(retained_row_basis)}"
    )

    print(
        f"  full_rowspace_rank="
        f"{len(full_row_basis)}"
    )

    print(
        f"  quotient_dimension="
        f"{len(full_row_quotient)}"
    )

    for i, w in enumerate(
        full_row_quotient
    ):
        print(
            f"  quotient_witness_{i}={w}"
        )

    # ------------------------------------------------------------------
    # 4. DIMENSION-CORRECT COLUMN-SPACE COMPARISON
    # ------------------------------------------------------------------

    """
    A 12x14 retained matrix has columns living in F_2^12.

    A 14x14 full matrix has columns living in F_2^14.

    We therefore embed the retained matrix into F_2^14 by adding zero
    rows in the two deleted boundary-equation positions.

    ONLY THEN do we compare column spaces.
    """

    retained_embedded = embed_retained_rows_zero(
        retained_rows_full,
        retained_indices,
        FULL_N,
    )

    retained_embedded_basis = row_space_basis_mod2(
        transpose(
            to_mod2(retained_embedded)
        )
    )

    full_column_basis = row_space_basis_mod2(
        transpose(M2)
    )

    full_column_quotient = quotient_witnesses(
        retained_embedded_basis,
        full_column_basis,
    )

    print()
    print("=" * 78)
    print("4. DIMENSION-CORRECT COLUMN-SPACE QUOTIENT")
    print("=" * 78)

    print(
        f"  retained_embedded_shape=14x14"
    )

    print(
        f"  retained_embedded_columnspace_rank="
        f"{len(retained_embedded_basis)}"
    )

    print(
        f"  full_columnspace_rank="
        f"{len(full_column_basis)}"
    )

    print(
        f"  quotient_dimension="
        f"{len(full_column_quotient)}"
    )

    for i, w in enumerate(
        full_column_quotient
    ):
        print(
            f"  column_quotient_witness_{i}={w}"
        )

    # ------------------------------------------------------------------
    # 5. KERNELS
    # ------------------------------------------------------------------

    full_right_kernel = nullspace_mod2(M2)
    full_left_kernel = nullspace_mod2(
        transpose(M2)
    )

    core_right_kernel = nullspace_mod2(A2)
    core_left_kernel = nullspace_mod2(
        transpose(A2)
    )

    projected_right = project_full_right_kernel(
        full_right_kernel
    )

    projected_left = project_full_left_kernel(
        full_left_kernel,
        retained_indices,
    )

    projected_right_rank = (
        rank_mod2(projected_right)
        if projected_right
        else 0
    )

    projected_left_rank = (
        rank_mod2(projected_left)
        if projected_left
        else 0
    )

    print()
    print("=" * 78)
    print("5. KERNEL PROJECTION")
    print("=" * 78)

    print(
        f"  full_right_nullity="
        f"{len(full_right_kernel)}"
    )

    print(
        f"  projected_right_rank="
        f"{projected_right_rank}"
    )

    print(
        f"  core_right_nullity="
        f"{len(core_right_kernel)}"
    )

    print(
        f"  full_left_nullity="
        f"{len(full_left_kernel)}"
    )

    print(
        f"  projected_left_rank="
        f"{projected_left_rank}"
    )

    print(
        f"  core_left_nullity="
        f"{len(core_left_kernel)}"
    )

    # ------------------------------------------------------------------
    # 6. ADJUGATE LINE — ONLY 12D SPACES
    # ------------------------------------------------------------------

    adj_row_in_core = in_span_mod2(
        adj_row_dir,
        core_left_kernel,
    )

    adj_col_in_core = in_span_mod2(
        adj_col_dir,
        core_right_kernel,
    )

    adj_row_in_projected_left = in_span_mod2(
        adj_row_dir,
        projected_left,
    )

    adj_col_in_projected_right = in_span_mod2(
        adj_col_dir,
        projected_right,
    )

    print()
    print("=" * 78)
    print("6. ADJUGATE-LINE MEMBERSHIP")
    print("=" * 78)

    print(
        f"  adj_row_direction={adj_row_dir}"
    )

    print(
        f"  adj_col_direction={adj_col_dir}"
    )

    print(
        f"  adj_row_in_core_left_kernel="
        f"{adj_row_in_core}"
    )

    print(
        f"  adj_col_in_core_right_kernel="
        f"{adj_col_in_core}"
    )

    print(
        f"  adj_row_in_projected_full_left="
        f"{adj_row_in_projected_left}"
    )

    print(
        f"  adj_col_in_projected_full_right="
        f"{adj_col_in_projected_right}"
    )

    # ------------------------------------------------------------------
    # 7. NULLITY GAPS
    # ------------------------------------------------------------------

    right_nullity_gap = (
        len(full_right_kernel)
        - len(core_right_kernel)
    )

    left_nullity_gap = (
        len(full_left_kernel)
        - len(core_left_kernel)
    )

    rank_gap = (
        full_rank - core_rank
    )

    print()
    print("=" * 78)
    print("7. RANK / NULLITY GAPS")
    print("=" * 78)

    print(
        f"  full_rank={full_rank}"
    )

    print(
        f"  core_rank={core_rank}"
    )

    print(
        f"  rank_gap={rank_gap}"
    )

    print(
        f"  right_nullity_gap="
        f"{right_nullity_gap}"
    )

    print(
        f"  left_nullity_gap="
        f"{left_nullity_gap}"
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
The previous 164 scripts failed because they compared vectors from
different ambient spaces.

This version fixes that at the construction level.

The correct spaces are:

    F_ret:
        12-dimensional row vectors;

    F_all:
        12-dimensional row vectors;

    retained_embedded column space:
        14-dimensional column vectors;

    full column space:
        14-dimensional column vectors;

    projected full kernels:
        12-dimensional vectors;

    core kernels:
        12-dimensional vectors.

The rational Schur block remains

    H = [[1,1],[1,3]],

with

    det(H)=2.

Thus H is invertible over Q but singular over F_2.

The experiment therefore does not treat the Schur complement as a
mod-2 unimodular transformation.

The main numerical relation being tested is

    rank_F2(M) - rank_F2(A) = 1,

together with the corresponding kernel projection and quotient
dimensions.

The distinguished adjugate line is tested only inside 12-dimensional
kernel spaces, so these membership statements are dimensionally exact.
"""
    )

    # ------------------------------------------------------------------
    # 9. FINAL EXACTNESS
    # ------------------------------------------------------------------

    rank_gap_exact = (
        rank_gap == 1
    )

    f_gap_exact = (
        f_all_rank - f_ret_rank == 1
    )

    f_quotient_exact = (
        len(f_quotient)
        == f_all_rank - f_ret_rank
    )

    full_row_quotient_exact = (
        len(full_row_quotient)
        == full_rank
        - len(retained_row_basis)
    )

    full_col_quotient_exact = (
        len(full_column_quotient)
        == full_rank
        - len(retained_embedded_basis)
    )

    final_ok = (
        full_rank == 4
        and core_rank == 3
        and f_ret_rank == 3
        and f_all_rank == 4
        and rank_gap_exact
        and f_gap_exact
        and f_quotient_exact
        and full_row_quotient_exact
        and full_col_quotient_exact
        and right_nullity_gap == 1
        and left_nullity_gap == 1
        and projected_right_rank == CORE_N - core_rank
        and projected_left_rank == CORE_N - core_rank
        and adj_row_in_core
        and adj_col_in_core
        and adj_row_in_projected_left
        and adj_col_in_projected_right
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  full_rank_mod2_exact={full_rank == 4}"
    )

    print(
        f"  core_rank_mod2_exact={core_rank == 3}"
    )

    print(
        f"  F_ret_rank_exact={f_ret_rank == 3}"
    )

    print(
        f"  F_all_rank_exact={f_all_rank == 4}"
    )

    print(
        f"  rank_gap_one={rank_gap_exact}"
    )

    print(
        f"  F_block_gap_one={f_gap_exact}"
    )

    print(
        f"  F_quotient_exact={f_quotient_exact}"
    )

    print(
        f"  full_row_quotient_exact="
        f"{full_row_quotient_exact}"
    )

    print(
        f"  full_column_quotient_exact="
        f"{full_col_quotient_exact}"
    )

    print(
        f"  right_nullity_gap_one="
        f"{right_nullity_gap == 1}"
    )

    print(
        f"  left_nullity_gap_one="
        f"{left_nullity_gap == 1}"
    )

    print(
        f"  projected_right_dimension_exact="
        f"{projected_right_rank == CORE_N - core_rank}"
    )

    print(
        f"  projected_left_dimension_exact="
        f"{projected_left_rank == CORE_N - core_rank}"
    )

    print(
        f"  adj_row_core_membership="
        f"{adj_row_in_core}"
    )

    print(
        f"  adj_col_core_membership="
        f"{adj_col_in_core}"
    )

    print(
        f"  adj_row_projected_full="
        f"{adj_row_in_projected_left}"
    )

    print(
        f"  adj_col_projected_full="
        f"{adj_col_in_projected_right}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 164U COMPLETE")


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