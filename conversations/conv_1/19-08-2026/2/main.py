#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 171S — EXACT SCHUR BOUNDARY-ROW / C-ROW SELECTION
PROVENANCE + ROW-PERMUTATION AUDIT
==============================================================================

Hardened replacement for Experiment 171R.

Established:

    L1 = (B adj(H))/2

is rank one modulo 2, with a single nonzero entry. Consequently K mod 2
is one selected row of C.

This experiment:

    1. identifies the two deleted boundary equations;
    2. traces the selected C row back to its full-system metadata;
    3. verifies exact reconstruction of K mod 2;
    4. permutes the two boundary equations;
    5. verifies that det(H) changes sign but the Schur correction K does not;
    6. tracks whether the same underlying boundary equation is selected.

All boolean values are computed before printing. No nested expressions
inside f-strings.

No SymPy.
No floating point.
No extrapolation.
No recurrence fitting.
No claim about the original (p,q)-kernel.
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
        if x & 1
    ]


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
                    f"not divisible by {d}."
                )

            new_row.append(x // d)

        out.append(new_row)

    return out


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
# FULL SYSTEM WITH METADATA
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
# SCHUR BLOCK DATA
# ============================================================================

def build_schur_data(M, metadata, boundary_order):

    retained_full_rows = [
        i
        for i in range(FULL_N)
        if i not in boundary_order
    ]

    retained_rows = [
        M[i]
        for i in retained_full_rows
    ]

    boundary_rows = [
        M[i]
        for i in boundary_order
    ]

    boundary_metadata = [
        metadata[i]
        for i in boundary_order
    ]

    # 12x12
    F = [
        [
            row[j]
            for j in range(CORE_N)
        ]
        for row in retained_rows
    ]

    # 12x2
    B = [
        [
            row[j]
            for j in BOUNDARY_COLS
        ]
        for row in retained_rows
    ]

    # 2x2
    H = [
        [
            row[j]
            for j in BOUNDARY_COLS
        ]
        for row in boundary_rows
    ]

    # 2x12
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
        "retained_full_rows": retained_full_rows,
        "boundary_metadata": boundary_metadata,
        "F": F,
        "B": B,
        "H": H,
        "C": C,
        "adj_H": adj_h,
        "det_H": det_h,
    }


def compute_K(B, adj_H, C, det_H):

    numerator = matrix_to_int(
        matmul(
            B,
            matmul(adj_H, C),
        ),
        "B*adj(H)*C",
    )

    if any(
        x % det_H != 0
        for row in numerator
        for x in row
    ):
        raise ArithmeticError(
            f"Schur numerator is not divisible by det(H)={det_H}."
        )

    return [
        [
            -(x // det_H)
            for x in row
        ]
        for row in numerator
    ]


def compute_L1(B, adj_H, det_H):

    numerator = matrix_to_int(
        matmul(B, adj_H),
        "B*adj(H)",
    )

    abs_det = abs(det_H)

    if any(
        x % abs_det != 0
        for row in numerator
        for x in row
    ):
        raise ArithmeticError(
            "B*adj(H) is not divisible by |det(H)|."
        )

    return [
        [
            x // abs_det
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
        "EXPERIMENT 171S — EXACT SCHUR BOUNDARY-ROW / "
        "C-ROW SELECTION PROVENANCE + ROW-PERMUTATION AUDIT"
    )
    print("=" * 78)

    M, metadata = build_full_system()

    # ==================================================================
    # ORIGINAL ORDER
    # ==================================================================

    original = build_schur_data(
        M,
        metadata,
        BOUNDARY_ROWS,
    )

    F = original["F"]
    B = original["B"]
    H = original["H"]
    C = original["C"]
    adj_H = original["adj_H"]
    det_H = original["det_H"]
    boundary_metadata = original["boundary_metadata"]

    L1 = compute_L1(
        B,
        adj_H,
        det_H,
    )

    L1_mod2 = to_mod2(L1)
    L1_support = matrix_support(L1_mod2)

    K = compute_K(
        B,
        adj_H,
        C,
        det_H,
    )

    K_mod2 = to_mod2(K)

    # ------------------------------------------------------------------
    # 1. BOUNDARY METADATA
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT BOUNDARY-EQUATION METADATA")
    print("=" * 78)

    for local_index, data in enumerate(boundary_metadata):

        print(
            f"  C_row={local_index}: "
            f"full_row={data['full_row']} "
            f"p={data['p']} "
            f"p_next={data['p_next']} "
            f"r={data['r']} "
            f"d={data['d']} "
            f"q0={data['q0']} "
            f"q1={data['q1']}"
        )

    # ------------------------------------------------------------------
    # 2. FIRST LIFT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. ORIGINAL FIRST-LIFT")
    print("=" * 78)

    print(f"  H={H}")
    print(f"  det_H={det_H}")
    print(f"  abs_det_H={abs(det_H)}")
    print(f"  L1_support={L1_support}")

    print("  L1_mod2:")

    for row in L1_mod2:
        print(
            "    "
            + " ".join(str(x) for x in row)
        )

    # ------------------------------------------------------------------
    # 3. SELECTED ROW
    # ------------------------------------------------------------------

    selected_core_row = None
    selected_c_row = None
    selected_meta = None

    if len(L1_support) == 1:

        selected_core_row = L1_support[0][0]
        selected_c_row = L1_support[0][1]

        selected_meta = boundary_metadata[
            selected_c_row
        ]

    single_matrix_unit = (
        len(L1_support) == 1
    )

    metadata_available = (
        selected_meta is not None
    )

    print()
    print("=" * 78)
    print("3. EXACT ROW-SELECTION PROVENANCE")
    print("=" * 78)

    print(
        f"  L1_single_matrix_unit="
        f"{single_matrix_unit}"
    )

    if metadata_available:

        print(
            f"  selected_core_row={selected_core_row}"
        )

        print(
            f"  selected_C_row={selected_c_row}"
        )

        print(
            f"  selected_full_boundary_row="
            f"{selected_meta['full_row']}"
        )

        print(
            f"  selected_p={selected_meta['p']}"
        )

        print(
            f"  selected_p_next={selected_meta['p_next']}"
        )

        print(
            f"  selected_r={selected_meta['r']}"
        )

        print(
            f"  selected_d={selected_meta['d']}"
        )

        print(
            f"  selected_q0={selected_meta['q0']}"
        )

        print(
            f"  selected_q1={selected_meta['q1']}"
        )

    # ------------------------------------------------------------------
    # 4. C ROW STRUCTURE
    # ------------------------------------------------------------------

    C_mod2 = to_mod2(C)

    print()
    print("=" * 78)
    print("4. EXACT C-ROW STRUCTURE")
    print("=" * 78)

    for i, row in enumerate(C_mod2):

        print(
            f"  C_row={i} "
            f"support={vector_support(row)} "
            f"row={row}"
        )

    # ------------------------------------------------------------------
    # 5. SELECTED C ROW RECONSTRUCTION
    # ------------------------------------------------------------------

    selected_reconstructs_K = False

    if metadata_available:

        predicted_K = [
            [
                0 for _ in range(CORE_N)
            ]
            for _ in range(CORE_N)
        ]

        for j in range(CORE_N):

            predicted_K[
                selected_core_row
            ][j] = C_mod2[
                selected_c_row
            ][j]

        selected_reconstructs_K = (
            predicted_K == K_mod2
        )

    print()
    print("=" * 78)
    print("5. SELECTED C-ROW RECONSTRUCTION")
    print("=" * 78)

    print(
        f"  K_support={matrix_support(K_mod2)}"
    )

    print(
        f"  selected_C_row_reconstructs_K="
        f"{selected_reconstructs_K}"
    )

    # ------------------------------------------------------------------
    # 6. UNSELECTED C ROW
    # ------------------------------------------------------------------

    other_row_exists = False
    other_row = None
    other_support = None
    other_prediction_differs = False

    if selected_c_row is not None:

        other_row_exists = True
        other_row = 1 - selected_c_row

        other_support = vector_support(
            C_mod2[other_row]
        )

        other_prediction = [
            [
                0 for _ in range(CORE_N)
            ]
            for _ in range(CORE_N)
        ]

        for j in range(CORE_N):

            other_prediction[
                selected_core_row
            ][j] = C_mod2[
                other_row
            ][j]

        other_prediction_differs = (
            other_prediction != K_mod2
        )

    print()
    print("=" * 78)
    print("6. UNSELECTED C-ROW AUDIT")
    print("=" * 78)

    print(
        f"  other_row_exists="
        f"{other_row_exists}"
    )

    if other_row_exists:

        print(
            f"  unselected_C_row={other_row}"
        )

        print(
            f"  unselected_support="
            f"{other_support}"
        )

    print(
        f"  replacing_selected_row_changes_K="
        f"{other_prediction_differs}"
    )

    # ==================================================================
    # SWAP
    # ==================================================================

    swapped_order = (
        BOUNDARY_ROWS[1],
        BOUNDARY_ROWS[0],
    )

    swapped = build_schur_data(
        M,
        metadata,
        swapped_order,
    )

    B_swap = swapped["B"]
    H_swap = swapped["H"]
    C_swap = swapped["C"]
    adj_swap = swapped["adj_H"]
    det_swap = swapped["det_H"]
    meta_swap = swapped["boundary_metadata"]

    L1_swap = compute_L1(
        B_swap,
        adj_swap,
        det_swap,
    )

    L1_swap_mod2 = to_mod2(
        L1_swap
    )

    K_swap = compute_K(
        B_swap,
        adj_swap,
        C_swap,
        det_swap,
    )

    K_swap_mod2 = to_mod2(
        K_swap
    )

    swap_support = matrix_support(
        L1_swap_mod2
    )

    # ------------------------------------------------------------------
    # 7. SWAP DATA
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. BOUNDARY-ROW SWAP")
    print("=" * 78)

    print(
        f"  swapped_order={swapped_order}"
    )

    print(
        f"  swapped_H={H_swap}"
    )

    print(
        f"  swapped_det_H={det_swap}"
    )

    print(
        f"  swapped_abs_det_H={abs(det_swap)}"
    )

    print(
        f"  swapped_L1_support={swap_support}"
    )

    swapped_K_exact = (
        K_swap == K
    )

    swapped_K_mod2 = (
        K_swap_mod2 == K_mod2
    )

    print(
        f"  K_exactly_invariant={swapped_K_exact}"
    )

    print(
        f"  K_mod2_invariant={swapped_K_mod2}"
    )

    # ------------------------------------------------------------------
    # 8. TRACK SELECTION THROUGH SWAP
    # ------------------------------------------------------------------

    swapped_selected_full_row = None

    if len(swap_support) == 1:

        swapped_c_row = swap_support[0][1]

        swapped_selected_full_row = meta_swap[
            swapped_c_row
        ]["full_row"]

    original_selected_full_row = None

    if metadata_available:

        original_selected_full_row = (
            selected_meta["full_row"]
        )

    same_underlying_boundary_equation = (
        original_selected_full_row is not None
        and swapped_selected_full_row is not None
        and (
            original_selected_full_row
            == swapped_selected_full_row
        )
    )

    print()
    print("=" * 78)
    print("8. SELECTION TRACKING UNDER PERMUTATION")
    print("=" * 78)

    print(
        f"  original_selected_full_row="
        f"{original_selected_full_row}"
    )

    print(
        f"  swapped_selected_full_row="
        f"{swapped_selected_full_row}"
    )

    print(
        f"  same_underlying_boundary_equation="
        f"{same_underlying_boundary_equation}"
    )

    # ------------------------------------------------------------------
    # 9. TERMINAL / TRANSITION INFORMATION
    # ------------------------------------------------------------------

    selected_terminal = False
    selected_transition = None

    if metadata_available:

        selected_terminal = (
            selected_meta["d"] == 0
        )

        selected_transition = (
            f"{selected_meta['p']}"
            f"->{selected_meta['p_next']}"
        )

    print()
    print("=" * 78)
    print("9. SELECTED-BOUNDARY LOCATION")
    print("=" * 78)

    print(
        f"  selected_is_terminal="
        f"{selected_terminal}"
    )

    print(
        f"  selected_transition="
        f"{selected_transition}"
    )

    # ------------------------------------------------------------------
    # 10. CORE RANK
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
    print("10. CORE RANK REFERENCE")
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
    # 11. INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("11. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 170R found a single matrix-unit entry in the first lifted
Schur factor L1.

Therefore K mod 2 is exactly one row of C, placed into one row of the
core.

Experiment 171S follows that selection back to the original full
system.

The boundary-row permutation is an essential control:

    swapping the two boundary equations changes det(H) from +2 to -2,

but the exact Schur correction must remain unchanged.

Thus the correct invariant statements are:

    |det(H)| = 2,
    K_swapped = K_original,
    K_swapped mod 2 = K_original mod 2.

The selection itself is then tracked by full-system row identity rather
than local C-row number. This distinguishes genuine boundary-equation
selection from a mere artifact of ordering the two rows.
"""
    )

    # ------------------------------------------------------------------
    # 12. FINAL
    # ------------------------------------------------------------------

    original_det_ok = (
        abs(det_H) == 2
    )

    swapped_det_ok = (
        abs(det_swap) == 2
    )

    final_ok = (
        original_det_ok
        and swapped_det_ok
        and single_matrix_unit
        and metadata_available
        and selected_reconstructs_K
        and swapped_K_exact
        and swapped_K_mod2
        and rank_K == 1
        and rank_F == 4
        and rank_A == 3
    )

    print()
    print("=" * 78)
    print("12. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  original_abs_det_H_2="
        f"{original_det_ok}"
    )

    print(
        f"  swapped_abs_det_H_2="
        f"{swapped_det_ok}"
    )

    print(
        f"  L1_single_matrix_unit="
        f"{single_matrix_unit}"
    )

    print(
        f"  boundary_metadata_available="
        f"{metadata_available}"
    )

    print(
        f"  selected_C_row_reconstructs_K="
        f"{selected_reconstructs_K}"
    )

    print(
        f"  boundary_swap_preserves_K="
        f"{swapped_K_exact}"
    )

    print(
        f"  boundary_swap_preserves_K_mod2="
        f"{swapped_K_mod2}"
    )

    print(
        f"  rank_K_one="
        f"{rank_K == 1}"
    )

    print(
        f"  rank_F_four="
        f"{rank_F == 4}"
    )

    print(
        f"  rank_A_three="
        f"{rank_A == 3}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 171S COMPLETE")


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