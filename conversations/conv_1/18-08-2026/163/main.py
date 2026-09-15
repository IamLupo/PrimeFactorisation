#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 163 — EXACT FULL-SYSTEM / SCHUR-COMPLEMENT ADJUGATE-LIFT AUDIT
==============================================================================

Corrected and hardened version of Experiment 162.

Main question:

    Is the distinguished one-dimensional mod-2 adjugate direction of the
    12x12 F-core already represented in the original 14x14 system, or is it
    created by the Schur-complement boundary elimination?

We compare:

    full 14x14 system M,
    reduced 12x12 F-core A,
    normalized adjugate AdjR = adj(A) / Delta_11.

All arithmetic is exact.

No SymPy.
No floating point.
No recurrence fitting.
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

FULL_N = 14
CORE_N = 12

BOUNDARY_COLUMNS = (12, 13)
BOUNDARY_ROWS = (5, 10)


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
            f"Expected {FULL_N} full-system rows, got {len(M)}."
        )

    return M, metadata


# ============================================================================
# SCHUR COMPLEMENT / F-CORE
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
            "Boundary block H is singular."
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

    # Boundary coefficients:
    # h = h_coeff * x
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

    core = []
    retained_rows = []

    for idx in range(FULL_N):

        if idx in BOUNDARY_ROWS:
            continue

        row = M[idx]
        out = []

        for j in range(CORE_N):

            value = (
                Fraction(row[j])
                + Fraction(row[12]) * h_coeff[0][j]
                + Fraction(row[13]) * h_coeff[1][j]
            )

            if value.denominator != 1:
                raise ArithmeticError(
                    "Schur complement produced a non-integral entry."
                )

            out.append(value.numerator)

        core.append(out)
        retained_rows.append(idx)

    if len(core) != CORE_N:
        raise ArithmeticError(
            "Expected a 12x12 F-core."
        )

    return core, H, h_coeff, retained_rows


# ============================================================================
# MOD-2 LINEAR ALGEBRA
# ============================================================================

def transpose(A):
    return [
        [A[i][j] for i in range(len(A))]
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
        j
        for j in range(n)
        if j not in pivot_set
    ]

    basis_vectors = []

    for free in free_cols:

        v = [0] * n
        v[free] = 1

        for i, pivot_col in enumerate(pivots):
            if R[i][free]:
                v[pivot_col] = 1

        basis_vectors.append(v)

    return basis_vectors


def in_span_mod2(v, basis):

    if not basis:
        return not any(v)

    return (
        rank_mod2(basis)
        == rank_mod2(basis + [v])
    )


def mat_vec_mod2(A, v):

    out = []

    for row in A:

        s = 0

        for a, b in zip(row, v):
            s ^= (a & 1) & (b & 1)

        out.append(s)

    return out


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

    C = []

    for i in range(len(A)):

        row = []

        for j in range(len(A)):

            value = bareiss_det(
                maximal_minor(A, i, j)
            )

            if (i + j) & 1:
                value = -value

            row.append(value)

        C.append(row)

    return C


def normalize_matrix(C):

    R = []

    for row in C:

        out = []

        for value in row:

            if value % DELTA_11 != 0:
                raise ArithmeticError(
                    "Delta_11 does not divide every cofactor."
                )

            out.append(value // DELTA_11)

        R.append(out)

    return R


def get_adjugate_line(AdjR):

    AdjR2 = [
        [x & 1 for x in row]
        for row in AdjR
    ]

    row_space = row_space_mod2(AdjR2)
    col_space = row_space_mod2(transpose(AdjR2))

    if len(row_space) != 1:
        raise ArithmeticError(
            "Expected rank-one adjugate row space."
        )

    if len(col_space) != 1:
        raise ArithmeticError(
            "Expected rank-one adjugate column space."
        )

    return row_space[0], col_space[0], AdjR2


# ============================================================================
# FULL RIGHT-KERNEL LIFT
# ============================================================================

def lift_core_right_to_full(core_vector, h_coeff):

    x = [Fraction(v) for v in core_vector]

    h0 = sum(
        h_coeff[0][j] * x[j]
        for j in range(CORE_N)
    )

    h1 = sum(
        h_coeff[1][j] * x[j]
        for j in range(CORE_N)
    )

    if h0.denominator != 1:
        raise ArithmeticError(
            "Boundary completion h0 is non-integral."
        )

    if h1.denominator != 1:
        raise ArithmeticError(
            "Boundary completion h1 is non-integral."
        )

    full = [
        int(x[j])
        for j in range(CORE_N)
    ]

    full.extend([
        h0.numerator,
        h1.numerator,
    ])

    if len(full) != FULL_N:
        raise ArithmeticError(
            "Bad full-vector length."
        )

    return full


# ============================================================================
# FIND EXACT FULL COMPLETIONS
# ============================================================================

def matching_full_right_vectors(full_kernel, target_core):

    matches = []

    for v in full_kernel:

        if v[:CORE_N] == target_core:
            matches.append(v)

    return matches


def matching_full_left_vectors(full_left_kernel,
                               retained_rows,
                               target_core):

    matches = []

    retained_set = set(retained_rows)

    for v in full_left_kernel:

        projected = [
            v[idx]
            for idx in retained_rows
        ]

        if projected == target_core:
            matches.append(v)

    return matches


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 163 — EXACT FULL-SYSTEM / "
        "SCHUR-COMPLEMENT ADJUGATE-LIFT AUDIT"
    )
    print("=" * 78)

    M, metadata = build_full_system()

    A, H, h_coeff, retained_rows = build_core(M)

    C = cofactor_matrix(A)
    R = normalize_matrix(C)
    AdjR = transpose(R)

    M2 = [
        [x & 1 for x in row]
        for row in M
    ]

    A2 = [
        [x & 1 for x in row]
        for row in A
    ]

    AdjR2 = [
        [x & 1 for x in row]
        for row in AdjR
    ]

    # ------------------------------------------------------------------
    # 1. FULL / CORE RANKS
    # ------------------------------------------------------------------

    full_rank = rank_mod2(M2)
    core_rank = rank_mod2(A2)
    adj_rank = rank_mod2(AdjR2)

    full_right_kernel = nullspace_mod2(M2)
    full_left_kernel = nullspace_mod2(transpose(M2))

    core_right_kernel = nullspace_mod2(A2)
    core_left_kernel = nullspace_mod2(transpose(A2))

    print()
    print("=" * 78)
    print("1. EXACT SYSTEM RANK / NULLITY")
    print("=" * 78)

    print(
        f"  full_shape={FULL_N}x{FULL_N}"
    )

    print(
        f"  full_rank_mod2={full_rank}"
    )

    print(
        f"  full_right_nullity={len(full_right_kernel)}"
    )

    print(
        f"  full_left_nullity={len(full_left_kernel)}"
    )

    print(
        f"  core_shape={CORE_N}x{CORE_N}"
    )

    print(
        f"  core_rank_mod2={core_rank}"
    )

    print(
        f"  core_right_nullity={len(core_right_kernel)}"
    )

    print(
        f"  core_left_nullity={len(core_left_kernel)}"
    )

    print(
        f"  adjugate_rank_mod2={adj_rank}"
    )

    # ------------------------------------------------------------------
    # 2. SCHUR BLOCK
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. SCHUR COMPLEMENT DATA")
    print("=" * 78)

    print(
        f"  boundary_rows={BOUNDARY_ROWS}"
    )

    print(
        f"  boundary_columns={BOUNDARY_COLUMNS}"
    )

    print(
        f"  H={H}"
    )

    print(
        f"  retained_rows={retained_rows}"
    )

    # ------------------------------------------------------------------
    # 3. ADJUGATE DIRECTIONS
    # ------------------------------------------------------------------

    adj_row_dir, adj_col_dir, _ = (
        get_adjugate_line(AdjR)
    )

    print()
    print("=" * 78)
    print("3. DISTINGUISHED ADJUGATE DIRECTIONS")
    print("=" * 78)

    print(
        f"  adj_row_direction={adj_row_dir}"
    )

    print(
        f"  adj_col_direction={adj_col_dir}"
    )

    # ------------------------------------------------------------------
    # 4. CORE KERNEL MEMBERSHIP
    # ------------------------------------------------------------------

    adj_row_core_member = in_span_mod2(
        adj_row_dir,
        core_left_kernel,
    )

    adj_col_core_member = in_span_mod2(
        adj_col_dir,
        core_right_kernel,
    )

    print()
    print("=" * 78)
    print("4. CORE KERNEL MEMBERSHIP")
    print("=" * 78)

    print(
        f"  adj_row_in_core_left_kernel="
        f"{adj_row_core_member}"
    )

    print(
        f"  adj_col_in_core_right_kernel="
        f"{adj_col_core_member}"
    )

    # ------------------------------------------------------------------
    # 5. PROJECT FULL RIGHT KERNEL
    # ------------------------------------------------------------------

    projected_full_right = [
        v[:CORE_N]
        for v in full_right_kernel
    ]

    projected_right_rank = (
        rank_mod2(projected_full_right)
        if projected_full_right
        else 0
    )

    adj_col_in_projected_right = in_span_mod2(
        adj_col_dir,
        projected_full_right,
    )

    print()
    print("=" * 78)
    print("5. PROJECTED FULL RIGHT KERNEL")
    print("=" * 78)

    print(
        f"  projected_rank={projected_right_rank}"
    )

    print(
        f"  adj_col_in_projected_full_right="
        f"{adj_col_in_projected_right}"
    )

    # ------------------------------------------------------------------
    # 6. EXPLICIT RIGHT COMPLETION
    # ------------------------------------------------------------------

    right_candidate = lift_core_right_to_full(
        adj_col_dir,
        h_coeff,
    )

    right_candidate_mod2 = [
        x & 1
        for x in right_candidate
    ]

    right_residual = mat_vec_mod2(
        M2,
        right_candidate_mod2,
    )

    right_completion_exact = (
        not any(right_residual)
    )

    print()
    print("=" * 78)
    print("6. EXPLICIT RIGHT BOUNDARY COMPLETION")
    print("=" * 78)

    print(
        f"  core_vector={adj_col_dir}"
    )

    print(
        f"  full_vector_mod2={right_candidate_mod2}"
    )

    print(
        f"  M*v_mod2={right_residual}"
    )

    print(
        f"  exact_full_right_kernel="
        f"{right_completion_exact}"
    )

    # ------------------------------------------------------------------
    # 7. PROJECT FULL LEFT KERNEL
    # ------------------------------------------------------------------

    projected_full_left = [
        [
            v[idx]
            for idx in retained_rows
        ]
        for v in full_left_kernel
    ]

    projected_left_rank = (
        rank_mod2(projected_full_left)
        if projected_full_left
        else 0
    )

    adj_row_in_projected_left = in_span_mod2(
        adj_row_dir,
        projected_full_left,
    )

    print()
    print("=" * 78)
    print("7. PROJECTED FULL LEFT KERNEL")
    print("=" * 78)

    print(
        f"  projected_rank={projected_left_rank}"
    )

    print(
        f"  adj_row_in_projected_full_left="
        f"{adj_row_in_projected_left}"
    )

    # ------------------------------------------------------------------
    # 8. EXACT LEFT COMPLETION SEARCH
    # ------------------------------------------------------------------

    left_matches = (
        matching_full_left_vectors(
            full_left_kernel,
            retained_rows,
            adj_row_dir,
        )
    )

    left_completion_exact = (
        len(left_matches) > 0
    )

    print()
    print("=" * 78)
    print("8. LEFT BOUNDARY COMPLETION")
    print("=" * 78)

    print(
        f"  matching_completions={len(left_matches)}"
    )

    for i, v in enumerate(left_matches):

        print(
            f"  completion_{i}={v}"
        )

    # ------------------------------------------------------------------
    # 9. DIMENSION COMPARISON
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. DIMENSION / INTERSECTION AUDIT")
    print("=" * 78)

    print(
        f"  core_left_nullity={len(core_left_kernel)}"
    )

    print(
        f"  projected_full_left_rank={projected_left_rank}"
    )

    print(
        f"  core_right_nullity={len(core_right_kernel)}"
    )

    print(
        f"  projected_full_right_rank={projected_right_rank}"
    )

    print(
        f"  adj_line_rank={adj_rank}"
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
The reduced 12x12 F-core is obtained from the original 14x14 system
by eliminating two boundary variables.

Experiment 163 asks a precise provenance question:

    Is the distinguished one-dimensional adjugate line of the
    reduced core already visible in the kernel of the original
    14-dimensional system?

There are two independent directions:

    right:
        core column direction -> full variable nullspace;

    left:
        core row direction -> full equation nullspace.

For the right direction the Schur complement gives an explicit
boundary completion.

For the left direction the experiment searches the exact full
nullspace for a completion whose retained coordinates equal the
distinguished core adjugate row.

This distinguishes:

    intrinsic full-system structure
        from
    structure created by the boundary elimination.

A positive lift would be materially stronger than the previous
core-only result.

A failure would show that the mod-2 adjugate line is genuinely
specific to the reduced F-core representation.

No causal interpretation is attached to either outcome.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        full_rank == 14
        and core_rank == 3
        and adj_rank == 1
        and adj_col_core_member
        and adj_row_core_member
        and adj_col_in_projected_right
        and right_completion_exact
        and adj_row_in_projected_left
        and left_completion_exact
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  full_system_valid={full_rank == 14}"
    )

    print(
        f"  core_rank_valid={core_rank == 3}"
    )

    print(
        f"  adjugate_rank_one={adj_rank == 1}"
    )

    print(
        f"  adj_column_core_kernel={adj_col_core_member}"
    )

    print(
        f"  adj_row_core_kernel={adj_row_core_member}"
    )

    print(
        f"  adj_column_projected_full={adj_col_in_projected_right}"
    )

    print(
        f"  explicit_right_completion={right_completion_exact}"
    )

    print(
        f"  adj_row_projected_full={adj_row_in_projected_left}"
    )

    print(
        f"  left_completion={left_completion_exact}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 163 COMPLETE")


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