#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 175 — EXACT SCHUR TERMINAL CONTRIBUTION /
                2-ADIC K-LAYER PROPAGATION AUDIT
==============================================================================

Experiment 174 established:

    C0 = C1 (mod 2),
    C0 != C1 (mod 4),

and

    (C0-C1)/2 mod 2

is supported only in the p-coordinate.

Experiment 175 asks whether that first boundary-row difference actually
propagates into the Schur correction K.

Write

    L1 = (B adj(H))/2

and let

    L1[:,0] = a
    L1[:,1] = b.

Then exactly

    K = -(a C0 + b C1).

This experiment separates the two terminal-boundary contributions:

    K0 = -a C0
    K1 = -b C1

so that

    K = K0 + K1.

The exact questions are:

    * What are the 2-adic valuations of a and b?
    * Does a dominate b at low 2-adic precision?
    * At what modulus does the C0-C1 difference become visible in K?
    * Does K mod 4 equal the contribution from C0 alone?
    * At which power of 2 does the second boundary row first matter?
    * How does the result change under swapping the two boundary rows?

The goal is to distinguish:

    boundary-row difference
        from
    boundary-row difference after Schur weighting.

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
# HELPERS
# ============================================================================

def q_value(p, r):
    if r < 0 or r >= len(Q[p]):
        return 0
    return Q[p][r]


def basis(p, d):
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


def to_mod(A, modulus):
    assert_rectangular(A)

    return [
        [int(x) % modulus for x in row]
        for row in A
    ]


def to_mod2(A):
    return to_mod(A, 2)


def matrix_support(A):
    return [
        (i, j)
        for i in range(len(A))
        for j in range(len(A[0]))
        if A[i][j] != 0
    ]


def vector_support(v):
    return [
        i
        for i, x in enumerate(v)
        if x != 0
    ]


def v2(x):
    if x == 0:
        return None

    x = abs(x)
    e = 0

    while x % 2 == 0:
        x //= 2
        e += 1

    return e


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


def divide_exact(A, divisor, name):
    out = []

    for i, row in enumerate(A):
        new_row = []

        for j, x in enumerate(row):
            if x % divisor != 0:
                raise ArithmeticError(
                    f"{name}[{i}][{j}]={x} is not divisible by {divisor}."
                )

            new_row.append(x // divisor)

        out.append(new_row)

    return out


# ============================================================================
# MATRIX OPERATIONS
# ============================================================================

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


def transpose(A):
    assert_rectangular(A, "transpose input")

    if not A:
        return []

    return [
        [A[i][j] for i in range(len(A))]
        for j in range(len(A[0]))
    ]


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


# ============================================================================
# FULL SYSTEM
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
                    "full_row": len(M) - 1,
                    "p": p,
                    "p_next": p_next,
                    "r": r,
                    "d": d,
                    "q0": q0,
                    "q1": q1,
                }
            )

    if len(M) != FULL_N:
        raise ArithmeticError(
            f"Expected {FULL_N} rows, got {len(M)}."
        )

    return M, metadata


# ============================================================================
# SCHUR BLOCKS
# ============================================================================

def build_schur_data(M, metadata, boundary_order):

    retained_indices = [
        i for i in range(FULL_N)
        if i not in boundary_order
    ]

    retained = [
        M[i] for i in retained_indices
    ]

    boundary = [
        M[i] for i in boundary_order
    ]

    boundary_metadata = [
        metadata[i] for i in boundary_order
    ]

    F = [
        [
            row[j] for j in range(CORE_N)
        ]
        for row in retained
    ]

    B = [
        [
            row[j] for j in BOUNDARY_COLS
        ]
        for row in retained
    ]

    H = [
        [
            row[j] for j in BOUNDARY_COLS
        ]
        for row in boundary
    ]

    C = [
        [
            row[j] for j in range(CORE_N)
        ]
        for row in boundary
    ]

    det_h = (
        H[0][0] * H[1][1]
        - H[0][1] * H[1][0]
    )

    if abs(det_h) != 2:
        raise ArithmeticError(
            f"Expected |det(H)|=2, got {det_h}."
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

    return {
        "F": F,
        "B": B,
        "H": H,
        "C": C,
        "adj_H": adj_h,
        "det_H": det_h,
        "metadata": boundary_metadata,
    }


# ============================================================================
# SCHUR FIRST LIFT
# ============================================================================

def compute_L1(B, adj_H, det_H):

    numerator = matrix_to_int(
        matmul(B, adj_H),
        "B*adj(H)"
    )

    divisor = abs(det_H)

    return divide_exact(
        numerator,
        divisor,
        "B*adj(H)"
    )


def compute_K(B, adj_H, C, det_H):

    numerator = matrix_to_int(
        matmul(
            B,
            matmul(adj_H, C)
        ),
        "Schur numerator"
    )

    if any(
        x % det_H != 0
        for row in numerator
        for x in row
    ):
        raise ArithmeticError(
            "Schur numerator is not divisible by det(H)."
        )

    return [
        [
            -(x // det_H)
            for x in row
        ]
        for row in numerator
    ]


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 175 — EXACT SCHUR TERMINAL CONTRIBUTION / "
        "2-ADIC K-LAYER PROPAGATION AUDIT"
    )
    print("=" * 78)

    M, metadata = build_full_system()

    data = build_schur_data(
        M,
        metadata,
        BOUNDARY_ROWS,
    )

    F = data["F"]
    B = data["B"]
    H = data["H"]
    C = data["C"]
    adj_H = data["adj_H"]
    det_H = data["det_H"]
    boundary_metadata = data["metadata"]

    # ------------------------------------------------------------------
    # 1. BOUNDARY DATA
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. TERMINAL BOUNDARY DATA")
    print("=" * 78)

    for i, meta in enumerate(boundary_metadata):

        print(
            f"  boundary={i} "
            f"full_row={meta['full_row']} "
            f"p={meta['p']} "
            f"r={meta['r']} "
            f"d={meta['d']} "
            f"q0={meta['q0']}"
        )

    print(
        f"  H={H}"
    )

    print(
        f"  det_H={det_H}"
    )

    # ------------------------------------------------------------------
    # 2. FIRST LIFT L1
    # ------------------------------------------------------------------

    L1 = compute_L1(
        B,
        adj_H,
        det_H,
    )

    L1_mod2 = to_mod2(L1)

    print()
    print("=" * 78)
    print("2. EXACT FIRST-LIFT SELECTOR")
    print("=" * 78)

    print(
        f"  L1_shape={len(L1)}x{len(L1[0])}"
    )

    print(
        f"  L1_mod2_support="
        f"{matrix_support(L1_mod2)}"
    )

    for i, row in enumerate(L1_mod2):
        if any(row):
            print(
                f"  L1_mod2_row_{i}={row}"
            )

    # ------------------------------------------------------------------
    # 3. EXTRACT TERMINAL COEFFICIENTS
    # ------------------------------------------------------------------

    a = [
        row[0]
        for row in L1
    ]

    b = [
        row[1]
        for row in L1
    ]

    print()
    print("=" * 78)
    print("3. TWO TERMINAL SCHUR COEFFICIENT VECTORS")
    print("=" * 78)

    print(
        f"  a=L1_column_0={a}"
    )

    print(
        f"  b=L1_column_1={b}"
    )

    print(
        f"  nonzero_a="
        f"{[i for i, x in enumerate(a) if x != 0]}"
    )

    print(
        f"  nonzero_b="
        f"{[i for i, x in enumerate(b) if x != 0]}"
    )

    # ------------------------------------------------------------------
    # 4. 2-ADIC PROFILE OF THE SELECTOR
    # ------------------------------------------------------------------

    a_v2 = [
        v2(x)
        for x in a
    ]

    b_v2 = [
        v2(x)
        for x in b
    ]

    print()
    print("=" * 78)
    print("4. SELECTOR 2-ADIC PROFILE")
    print("=" * 78)

    print(
        f"  a_v2={a_v2}"
    )

    print(
        f"  b_v2={b_v2}"
    )

    # ------------------------------------------------------------------
    # 5. INDIVIDUAL CONTRIBUTIONS
    # ------------------------------------------------------------------

    C0 = C[0]
    C1 = C[1]

    contribution0 = [
        [
            -a[i] * C0[j]
            for j in range(CORE_N)
        ]
        for i in range(CORE_N)
    ]

    contribution1 = [
        [
            -b[i] * C1[j]
            for j in range(CORE_N)
        ]
        for i in range(CORE_N)
    ]

    K = [
        [
            contribution0[i][j]
            + contribution1[i][j]
            for j in range(CORE_N)
        ]
        for i in range(CORE_N)
    ]

    K_mod2 = to_mod2(K)

    print()
    print("=" * 78)
    print("5. EXACT TWO-BOUNDARY DECOMPOSITION")
    print("=" * 78)

    print(
        f"  K0_support_mod2="
        f"{matrix_support(to_mod2(contribution0))}"
    )

    print(
        f"  K1_support_mod2="
        f"{matrix_support(to_mod2(contribution1))}"
    )

    print(
        f"  K_support_mod2="
        f"{matrix_support(K_mod2)}"
    )

    # ------------------------------------------------------------------
    # 6. MODULAR PROPAGATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. MODULAR CONTRIBUTION COMPARISON")
    print("=" * 78)

    contribution_checks = []

    for e in range(1, 9):

        modulus = 2 ** e

        K0_mod = to_mod(
            contribution0,
            modulus
        )

        K1_mod = to_mod(
            contribution1,
            modulus
        )

        K_mod = to_mod(
            K,
            modulus
        )

        K0_only_exact_mod = (
            K_mod == K0_mod
        )

        K1_only_exact_mod = (
            K_mod == K1_mod
        )

        contribution_checks.append(
            (
                modulus,
                K0_only_exact_mod,
                K1_only_exact_mod,
            )
        )

        print(
            f"  modulus={modulus} "
            f"K_equals_K0={K0_only_exact_mod} "
            f"K_equals_K1={K1_only_exact_mod}"
        )

    # ------------------------------------------------------------------
    # 7. FIRST MODULUS WHERE SECOND BOUNDARY MATTERS
    # ------------------------------------------------------------------

    first_second_boundary_modulus = None

    for modulus in [
        2 ** e
        for e in range(1, 16)
    ]:

        if (
            to_mod(contribution0, modulus)
            != to_mod(K, modulus)
        ):

            first_second_boundary_modulus = modulus
            break

    print()
    print("=" * 78)
    print("7. FIRST VISIBLE SECOND-BOUNDARY CONTRIBUTION")
    print("=" * 78)

    print(
        f"  first_modulus_where_K_differs_from_K0="
        f"{first_second_boundary_modulus}"
    )

    # ------------------------------------------------------------------
    # 8. FIRST VISIBLE CONTRIBUTION OF BOUNDARY 1
    # ------------------------------------------------------------------

    K1_entries = [
        x
        for row in contribution1
        for x in row
        if x != 0
    ]

    if K1_entries:

        min_K1_v2 = min(
            v2(x)
            for x in K1_entries
            if x != 0
        )

    else:

        min_K1_v2 = None

    print()
    print("=" * 78)
    print("8. K1 2-ADIC ORDER")
    print("=" * 78)

    print(
        f"  nonzero_K1_entries="
        f"{len(K1_entries)}"
    )

    print(
        f"  minimum_nonzero_v2(K1)="
        f"{min_K1_v2}"
    )

    # ------------------------------------------------------------------
    # 9. SECOND-BOUNDARY CONTRIBUTION TO THE P-COORDINATE
    # ------------------------------------------------------------------

    second_boundary_p_values = [
        contribution1[i][3]
        for i in range(CORE_N)
        if contribution1[i][3] != 0
    ]

    print()
    print("=" * 78)
    print("9. SECOND-BOUNDARY p-COORDINATE CONTRIBUTION")
    print("=" * 78)

    print(
        f"  nonzero_values="
        f"{second_boundary_p_values}"
    )

    if second_boundary_p_values:

        print(
            f"  valuations="
            f"{[v2(x) for x in second_boundary_p_values]}"
        )

    # ------------------------------------------------------------------
    # 10. SWAP CONTROL
    # ------------------------------------------------------------------

    swapped_order = (
        BOUNDARY_ROWS[1],
        BOUNDARY_ROWS[0],
    )

    swapped = build_schur_data(
        M,
        metadata,
        swapped_order,
    )

    L1_swap = compute_L1(
        swapped["B"],
        swapped["adj_H"],
        swapped["det_H"],
    )

    K_swap = compute_K(
        swapped["B"],
        swapped["adj_H"],
        swapped["C"],
        swapped["det_H"],
    )

    swap_K_exact = (
        K_swap == K
    )

    print()
    print("=" * 78)
    print("10. BOUNDARY-ORDER SWAP CONTROL")
    print("=" * 78)

    print(
        f"  swapped_det_H={swapped['det_H']}"
    )

    print(
        f"  swapped_L1_mod2_support="
        f"{matrix_support(to_mod2(L1_swap))}"
    )

    print(
        f"  swapped_K_exactly_same="
        f"{swap_K_exact}"
    )

    # ------------------------------------------------------------------
    # 11. CORE RANK
    # ------------------------------------------------------------------

    A = [
        [
            F[i][j] + K[i][j]
            for j in range(CORE_N)
        ]
        for i in range(CORE_N)
    ]

    rank_F = rank_mod2(F)
    rank_K = rank_mod2(K_mod2)
    rank_A = rank_mod2(A)

    print()
    print("=" * 78)
    print("11. CORE RANK REFERENCE")
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
    # 12. INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("12. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 174 showed that the two terminal C-rows differ first at
2-adic order one, specifically in the p-coordinate.

Experiment 175 asks whether that difference actually contributes to
the leading Schur correction.

The decomposition

    K = K0 + K1

is exact, where

    K0 = -a C0
    K1 = -b C1.

The first-lift data from Experiment 170R suggests that modulo 2

    a = e_11,
    b = 0.

Therefore

    K mod 2 = K0 mod 2.

The new question is the next precision level.

If K remains equal to K0 modulo 4, then the C0-C1 difference is
suppressed by the Schur coefficient b through the first two
2-adic levels.

If K differs from K0 modulo 4, then the mod-4 distinction between
the terminal C rows propagates immediately into K.

The result separates:

    boundary-row arithmetic
        from
    boundary-row arithmetic after Schur weighting.
"""
    )

    # ------------------------------------------------------------------
    # 13. FINAL
    # ------------------------------------------------------------------

    selector_rank_ok = (
        rank_mod2(L1_mod2) == 1
    )

    K_rank_ok = (
        rank_K == 1
    )

    F_rank_ok = (
        rank_F == 4
    )

    A_rank_ok = (
        rank_A == 3
    )

    final_ok = (
        abs(det_H) == 2
        and selector_rank_ok
        and K_rank_ok
        and F_rank_ok
        and A_rank_ok
        and left_right_equivalent(K, contribution0, contribution1)
        and swap_K_exact
    )

    print()
    print("=" * 78)
    print("13. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  abs_det_H_2={abs(det_H) == 2}"
    )

    print(
        f"  selector_rank_one={selector_rank_ok}"
    )

    print(
        f"  exact_K_decomposition="
        f"{left_right_equivalent(K, contribution0, contribution1)}"
    )

    print(
        f"  K_rank_one={K_rank_ok}"
    )

    print(
        f"  F_rank_four={F_rank_ok}"
    )

    print(
        f"  A_rank_three={A_rank_ok}"
    )

    print(
        f"  boundary_swap_preserves_K="
        f"{swap_K_exact}"
    )

    print(
        f"  first_second_boundary_modulus="
        f"{first_second_boundary_modulus}"
    )

    print(
        f"  minimum_K1_v2="
        f"{min_K1_v2}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 175 COMPLETE")


# ============================================================================
# SAFE FINAL HELPER
# ============================================================================

def left_right_equivalent(K, K0, K1):

    for i in range(len(K)):
        for j in range(len(K[0])):

            if K[i][j] != K0[i][j] + K1[i][j]:
                return False

    return True


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
