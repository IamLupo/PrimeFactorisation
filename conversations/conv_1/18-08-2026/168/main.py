#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 168R — EXACT SCHUR NUMERATOR / H-ADJUGATE /
MOD-2 RANK-ONE ORIGIN AUDIT
==============================================================================

Corrected replacement for Experiment 168.

Established facts:

    rank_F2(F_ret) = 4
    rank_F2(A)     = 3
    rank_F2(K)     = 1

with

    A = F_ret + K

and

    K = -(1/2) B adj(H) C,

where

    H = [[1,1],[1,3]]
    det(H) = 2.

This experiment asks whether the rank-one K mod-2 structure is already
visible in the factorization

    B adj(H) C.

No SymPy.
No floating point.
No extrapolation.
No recurrence fitting.
No claim about the original (p,q)-kernel.
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

FULL_N = 14
CORE_N = 12

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


def to_mod2(A):

    assert_rectangular(A)

    return [
        [int(x) & 1 for x in row]
        for row in A
    ]


def transpose(A):

    assert_rectangular(A)

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
                out[i][j] += (
                    a * Fraction(B[k][j])
                )

    return out


def matrix_to_int(A, name):

    out = []

    for i, row in enumerate(A):

        new_row = []

        for j, value in enumerate(row):

            value = Fraction(value)

            if value.denominator != 1:
                raise ArithmeticError(
                    f"{name}[{i}][{j}] is non-integral: {value}"
                )

            new_row.append(value.numerator)

        out.append(new_row)

    return out


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


# ============================================================================
# SCHUR BLOCKS
# ============================================================================

def build_schur_blocks(M):

    retained_indices = [
        i for i in range(FULL_N)
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

    # F_ret : 12x12
    F = [
        [
            row[j]
            for j in range(CORE_N)
        ]
        for row in retained_rows
    ]

    # B : 12x2
    B = [
        [
            row[j]
            for j in BOUNDARY_COLS
        ]
        for row in retained_rows
    ]

    # H : 2x2
    H = [
        [
            row[j]
            for j in BOUNDARY_COLS
        ]
        for row in boundary_rows
    ]

    # C : 2x12
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

    adj_H = [
        [
            H[1][1],
            -H[0][1],
        ],
        [
            -H[1][0],
            H[0][0],
        ],
    ]

    return (
        retained_indices,
        F,
        B,
        H,
        C,
        H_inv,
        adj_H,
        det_H,
    )


# ============================================================================
# MOD-2 RANK
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


def outer_product2(u, v):

    return [
        [
            (u[i] & v[j])
            for j in range(len(v))
        ]
        for i in range(len(u))
    ]


def matrices_equal_mod2(A, B):

    if (
        len(A) != len(B)
        or len(A[0]) != len(B[0])
    ):
        return False

    return (
        to_mod2(A)
        == to_mod2(B)
    )


# ============================================================================
# COFACTOR
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
                    "Delta_11 divisibility failed."
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
        "EXPERIMENT 168R — EXACT SCHUR NUMERATOR / H-ADJUGATE / "
        "MOD-2 RANK-ONE ORIGIN AUDIT"
    )
    print("=" * 78)

    M = build_full_system()

    (
        retained_indices,
        F,
        B,
        H,
        C,
        H_inv,
        adj_H,
        det_H,
    ) = build_schur_blocks(M)

    # ------------------------------------------------------------------
    # 1. H / ADJ(H)
    # ------------------------------------------------------------------

    H2 = to_mod2(H)
    adj_H2 = to_mod2(adj_H)

    print()
    print("=" * 78)
    print("1. EXACT H / ADJ(H)")
    print("=" * 78)

    print(
        f"  H={H}"
    )

    print(
        f"  det_H={det_H}"
    )

    print(
        f"  H_mod2={H2}"
    )

    print(
        f"  adj_H={adj_H}"
    )

    print(
        f"  adj_H_mod2={adj_H2}"
    )

    print(
        f"  rank_H_mod2={rank_mod2(H2)}"
    )

    print(
        f"  rank_adj_H_mod2={rank_mod2(adj_H2)}"
    )

    # ------------------------------------------------------------------
    # 2. EXACT SCHUR NUMERATOR
    # ------------------------------------------------------------------

    S_fraction = matmul(
        B,
        matmul(
            adj_H,
            C,
        ),
    )

    S = matrix_to_int(
        S_fraction,
        "S",
    )

    nonzero_S = [
        x
        for row in S
        for x in row
        if x != 0
    ]

    all_even = all(
        x % 2 == 0
        for row in S
        for x in row
    )

    print()
    print("=" * 78)
    print("2. EXACT SCHUR NUMERATOR")
    print("=" * 78)

    print(
        f"  S_shape={len(S)}x{len(S[0])}"
    )

    print(
        f"  nonzero_entries={len(nonzero_S)}"
    )

    print(
        f"  every_entry_even={all_even}"
    )

    # ------------------------------------------------------------------
    # 3. K FROM S
    # ------------------------------------------------------------------

    if not all_even:
        raise ArithmeticError(
            "S is not entrywise divisible by 2."
        )

    K = [
        [
            -(S[i][j] // 2)
            for j in range(CORE_N)
        ]
        for i in range(CORE_N)
    ]

    # Direct rational Schur correction.
    direct = matmul(
        B,
        matmul(
            H_inv,
            C,
        ),
    )

    direct_K = matrix_to_int(
        [
            [
                -direct[i][j]
                for j in range(CORE_N)
            ]
            for i in range(CORE_N)
        ],
        "direct K",
    )

    K_exact = (
        K == direct_K
    )

    print()
    print("=" * 78)
    print("3. EXACT K RECOVERY")
    print("=" * 78)

    print(
        f"  K=-(B*adj_H*C)/2"
    )

    print(
        f"  exact_recovery={K_exact}"
    )

    # ------------------------------------------------------------------
    # 4. MOD-2 FACTORIZATION CHAIN
    # ------------------------------------------------------------------

    B2 = to_mod2(B)
    C2 = to_mod2(C)

    B_adjH = [
        [
            sum(
                B[i][k] * adj_H[k][j]
                for k in range(2)
            )
            for j in range(2)
        ]
        for i in range(CORE_N)
    ]

    adjH_C = [
        [
            sum(
                adj_H[i][k] * C[k][j]
                for k in range(2)
            )
            for j in range(CORE_N)
        ]
        for i in range(2)
    ]

    B_adjH2 = to_mod2(B_adjH)
    adjH_C2 = to_mod2(adjH_C)
    K2 = to_mod2(K)

    rank_B = rank_mod2(B2)
    rank_C = rank_mod2(C2)
    rank_B_adjH = rank_mod2(B_adjH2)
    rank_adjH_C = rank_mod2(adjH_C2)
    rank_K = rank_mod2(K2)

    print()
    print("=" * 78)
    print("4. MOD-2 FACTORIZATION CHAIN")
    print("=" * 78)

    print(
        f"  rank(B_mod2)={rank_B}"
    )

    print(
        f"  rank(C_mod2)={rank_C}"
    )

    print(
        f"  rank(B*adj_H_mod2)="
        f"{rank_B_adjH}"
    )

    print(
        f"  rank(adj_H*C_mod2)="
        f"{rank_adjH_C}"
    )

    print(
        f"  rank(K_mod2)={rank_K}"
    )

    # ------------------------------------------------------------------
    # 5. OUTER PRODUCT TEST
    # ------------------------------------------------------------------

    left_basis = row_basis_mod2(
        B_adjH2
    )

    right_basis = row_basis_mod2(
        adjH_C2
    )

    outer_exact = False

    left_dir = None
    right_dir = None

    if (
        len(left_basis) == 1
        and len(right_basis) == 1
    ):

        left_dir = left_basis[0]
        right_dir = right_basis[0]

        predicted = outer_product2(
            left_dir,
            right_dir,
        )

        outer_exact = (
            predicted == K2
        )

    print()
    print("=" * 78)
    print("5. INTERMEDIATE OUTER-PRODUCT TEST")
    print("=" * 78)

    print(
        f"  left_factor_rank={len(left_basis)}"
    )

    print(
        f"  right_factor_rank={len(right_basis)}"
    )

    if left_dir is not None:
        print(
            f"  left_direction={left_dir}"
        )

    if right_dir is not None:
        print(
            f"  right_direction={right_dir}"
        )

    print(
        f"  K_outer_product_exact="
        f"{outer_exact}"
    )

    # ------------------------------------------------------------------
    # 6. DIRECT K SUPPORT
    # ------------------------------------------------------------------

    support_K = [
        (i, j)
        for i in range(CORE_N)
        for j in range(CORE_N)
        if K2[i][j]
    ]

    support_Badj = [
        (i, j)
        for i in range(CORE_N)
        for j in range(2)
        if B_adjH2[i][j]
    ]

    support_adjC = [
        (i, j)
        for i in range(2)
        for j in range(CORE_N)
        if adjH_C2[i][j]
    ]

    print()
    print("=" * 78)
    print("6. SUPPORT AUDIT")
    print("=" * 78)

    print(
        f"  K_nonzero_positions={support_K}"
    )

    print(
        f"  B_adjH_nonzero_positions="
        f"{support_Badj}"
    )

    print(
        f"  adjH_C_nonzero_positions="
        f"{support_adjC}"
    )

    # ------------------------------------------------------------------
    # 7. RANK OF A / F
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
    print("7. RANK TRANSITION")
    print("=" * 78)

    print(
        f"  rank_F_ret_mod2={rank_F}"
    )

    print(
        f"  rank_K_mod2={rank_K}"
    )

    print(
        f"  rank_A_mod2={rank_A}"
    )

    print(
        f"  rank_drop={rank_F-rank_A}"
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
The exact Schur formula is

    K = -(1/2) B adj(H) C.

The factor 1/2 is forced by

    det(H)=2.

The purpose of this experiment is to determine where the observed
rank-one residue K mod 2 first appears.

There are three possibilities:

    1. B adj(H) is already rank one modulo 2 and
       adj(H) C is rank one modulo 2.
       Then the rank-one K is structurally inherited from both sides.

    2. The intermediate products have larger rank, but their exact
       integer product is divisible by 2 in such a way that K mod 2
       collapses to rank one.
       Then additional arithmetic cancellation is required.

    3. The intermediate products vanish modulo 2 while K mod 2 remains
       nonzero.
       Then the first visible layer appears only after the 2-adic
       division and cannot be read directly from the mod-2 factors.

The experiment distinguishes these cases without assuming a mechanism.
"""
    )

    # ------------------------------------------------------------------
    # 9. FINAL
    # ------------------------------------------------------------------

    rank_F_ok = (
        rank_F == 4
    )

    rank_A_ok = (
        rank_A == 3
    )

    rank_K_ok = (
        rank_K == 1
    )

    final_ok = (
        det_H == 2
        and all_even
        and K_exact
        and rank_F_ok
        and rank_A_ok
        and rank_K_ok
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  det_H_equals_2={det_H == 2}"
    )

    print(
        f"  S_entrywise_even={all_even}"
    )

    print(
        f"  K_reconstruction_exact={K_exact}"
    )

    print(
        f"  rank_F_ret_four={rank_F_ok}"
    )

    print(
        f"  rank_K_one={rank_K_ok}"
    )

    print(
        f"  rank_A_three={rank_A_ok}"
    )

    print(
        f"  intermediate_outer_product_exact="
        f"{outer_exact}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 168R COMPLETE")


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