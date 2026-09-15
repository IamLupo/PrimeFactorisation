#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 172 — EXACT TERMINAL-BOUNDARY ROW / 2-ADIC DIFFERENCE AUDIT
==============================================================================

Experiment 171S established:

    C_row 0 and C_row 1 are identical modulo 2;

    the original first-lift selects C_row 0;

    after swapping the two boundary equations, the first-lift selects
    C_row 1;

    the exact Schur correction K is unchanged.

Thus the mod-2 selection is ordering-invariant because

    C_0 = C_1  (mod 2).

Experiment 172 studies how far that equality persists 2-adically.

Main objects:

    C0 = first boundary row restricted to the F-core;
    C1 = second boundary row restricted to the F-core;

    DeltaC = C0 - C1.

We compute:

    * exact DeltaC;
    * entrywise v_2(DeltaC);
    * the minimum v_2 across nonzero entries;
    * DeltaC / 2^e for successive exact layers;
    * the first modulus 2^e at which C0 and C1 differ;
    * whether C0 and C1 remain proportional modulo higher powers of 2;
    * the relationship between DeltaC and the first lifted Schur factor.

We deliberately do not assume that equality modulo 2 implies any
deeper symmetry.

No SymPy.
No floating point.
No extrapolation.
No recurrence fitting.
No connection to the original (p,q)-kernel is asserted.
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
        j
        for j, x in enumerate(v)
        if x != 0
    ]


def v2(x: int):
    if x == 0:
        return None

    x = abs(x)
    e = 0

    while x % 2 == 0:
        x //= 2
        e += 1

    return e


def exact_common_v2(values):
    nonzero = [
        abs(x)
        for x in values
        if x != 0
    ]

    if not nonzero:
        return None

    e = None

    for x in nonzero:
        vx = v2(x)

        if e is None:
            e = vx
        else:
            e = min(e, vx)

    return e


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


def divide_matrix_exact(A, d, name):

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
# SCHUR DATA
# ============================================================================

def build_schur_data(
    M,
    metadata,
    boundary_order,
):

    retained_rows = [
        i
        for i in range(FULL_N)
        if i not in boundary_order
    ]

    retained = [
        M[i]
        for i in retained_rows
    ]

    boundary = [
        M[i]
        for i in boundary_order
    ]

    boundary_metadata = [
        metadata[i]
        for i in boundary_order
    ]

    F = [
        [
            row[j]
            for j in range(CORE_N)
        ]
        for row in retained
    ]

    B = [
        [
            row[j]
            for j in BOUNDARY_COLS
        ]
        for row in retained
    ]

    H = [
        [
            row[j]
            for j in BOUNDARY_COLS
        ]
        for row in boundary
    ]

    C = [
        [
            row[j]
            for j in range(CORE_N)
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
        "boundary_metadata": boundary_metadata,
    }


# ============================================================================
# SCHUR OPERATIONS
# ============================================================================

def compute_K(B, adj_H, C, det_H):

    numerator = matrix_to_int(
        matmul(
            B,
            matmul(adj_H, C),
        ),
        "Schur numerator",
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


def compute_L1(B, adj_H, det_H):

    numerator = matrix_to_int(
        matmul(B, adj_H),
        "B adj(H)",
    )

    e = abs(det_H)

    if any(
        x % e != 0
        for row in numerator
        for x in row
    ):
        raise ArithmeticError(
            "B adj(H) is not divisible by |det(H)|."
        )

    return [
        [
            x // e
            for x in row
        ]
        for row in numerator
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
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 172 — EXACT TERMINAL-BOUNDARY ROW / "
        "2-ADIC DIFFERENCE AUDIT"
    )
    print("=" * 78)

    M, metadata = build_full_system()

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

    C0 = C[0]
    C1 = C[1]

    # ------------------------------------------------------------------
    # 1. BOUNDARY ROW METADATA
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. TERMINAL BOUNDARY ROWS")
    print("=" * 78)

    for i, meta in enumerate(boundary_metadata):

        print(
            f"  C_row={i}: "
            f"full_row={meta['full_row']} "
            f"p={meta['p']} "
            f"p_next={meta['p_next']} "
            f"r={meta['r']} "
            f"d={meta['d']} "
            f"q0={meta['q0']} "
            f"q1={meta['q1']}"
        )

    print(
        f"  H={H}"
    )

    print(
        f"  det_H={det_H}"
    )

    # ------------------------------------------------------------------
    # 2. EXACT C ROWS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. EXACT C ROWS")
    print("=" * 78)

    print(
        f"  C0={C0}"
    )

    print(
        f"  C1={C1}"
    )

    print(
        f"  C0_mod2={to_mod2([C0])[0]}"
    )

    print(
        f"  C1_mod2={to_mod2([C1])[0]}"
    )

    print(
        f"  equal_exact={C0 == C1}"
    )

    print(
        f"  equal_mod2="
        f"{to_mod2([C0]) == to_mod2([C1])}"
    )

    # ------------------------------------------------------------------
    # 3. EXACT DIFFERENCE
    # ------------------------------------------------------------------

    delta_C = [
        C0[j] - C1[j]
        for j in range(CORE_N)
    ]

    delta_v2 = [
        v2(x)
        for x in delta_C
    ]

    min_delta_v2 = exact_common_v2(
        delta_C
    )

    nonzero_delta_positions = [
        j
        for j, x in enumerate(delta_C)
        if x != 0
    ]

    print()
    print("=" * 78)
    print("3. EXACT C0-C1 DIFFERENCE")
    print("=" * 78)

    print(
        f"  DeltaC=C0-C1={delta_C}"
    )

    print(
        f"  nonzero_positions="
        f"{nonzero_delta_positions}"
    )

    print(
        f"  entrywise_v2={delta_v2}"
    )

    print(
        f"  minimum_nonzero_v2={min_delta_v2}"
    )

    # ------------------------------------------------------------------
    # 4. SUCCESSIVE 2-ADIC LAYERS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. SUCCESSIVE 2-ADIC DIFFERENCE LAYERS")
    print("=" * 78)

    layer_results = []

    if min_delta_v2 is not None:

        current_e = min_delta_v2

        for e in range(
            min_delta_v2,
            min_delta_v2 + 8,
        ):

            modulus = 2 ** e

            equal_mod_e = all(
                (x - y) % modulus == 0
                for x, y in zip(C0, C1)
            )

            if all(
                x % modulus == 0
                for x in delta_C
            ):

                stripped = [
                    x // modulus
                    for x in delta_C
                ]

                residues = [
                    x % 2
                    for x in stripped
                ]

            else:

                stripped = None
                residues = None

            layer_results.append(
                (
                    e,
                    modulus,
                    equal_mod_e,
                    residues,
                )
            )

            print(
                f"  e={e}: "
                f"modulus={modulus} "
                f"equal={equal_mod_e} "
                f"stripped_mod2={residues}"
            )

    # ------------------------------------------------------------------
    # 5. FIRST MODULUS WHERE ROWS DIFFER
    # ------------------------------------------------------------------

    first_difference_e = None

    for e in range(1, 64):

        modulus = 2 ** e

        if any(
            (x - y) % modulus != 0
            for x, y in zip(C0, C1)
        ):

            first_difference_e = e
            break

    print()
    print("=" * 78)
    print("5. FIRST 2-ADIC SEPARATION")
    print("=" * 78)

    print(
        f"  first_difference_e={first_difference_e}"
    )

    if first_difference_e is not None:

        print(
            f"  first_difference_modulus="
            f"{2 ** first_difference_e}"
        )

    # ------------------------------------------------------------------
    # 6. COMPARE WITH L1 SELECTION
    # ------------------------------------------------------------------

    L1 = compute_L1(
        B,
        adj_H,
        det_H,
    )

    L1_mod2 = to_mod2(L1)

    L1_support = matrix_support(
        L1_mod2
    )

    K = compute_K(
        B,
        adj_H,
        C,
        det_H,
    )

    K_mod2 = to_mod2(K)

    print()
    print("=" * 78)
    print("6. FIRST-LIFT SELECTION")
    print("=" * 78)

    print(
        f"  L1_support={L1_support}"
    )

    print(
        f"  K_support="
        f"{matrix_support(K_mod2)}"
    )

    # ------------------------------------------------------------------
    # 7. DOES C0-C1 MATTER TO K MOD 2?
    # ------------------------------------------------------------------

    K_difference_under_row_swap = [
        [
            0
            for _ in range(CORE_N)
        ]
        for _ in range(CORE_N)
    ]

    if len(L1_support) == 1:

        selected_row, selected_col = (
            L1_support[0]
        )

        for j in range(CORE_N):

            K_difference_under_row_swap[
                selected_row
            ][j] = (
                to_mod2([delta_C])[0][j]
            )

    K_delta_mod2_zero = (
        all(
            x == 0
            for row in K_difference_under_row_swap
            for x in row
        )
    )

    print()
    print("=" * 78)
    print("7. C-ROW DIFFERENCE VS K MOD 2")
    print("=" * 78)

    print(
        f"  DeltaC_mod2="
        f"{to_mod2([delta_C])[0]}"
    )

    print(
        f"  induced_K_difference_mod2="
        f"{K_difference_under_row_swap}"
    )

    print(
        f"  DeltaC_invisible_to_K_mod2="
        f"{K_delta_mod2_zero}"
    )

    # ------------------------------------------------------------------
    # 8. HIGHER-LAYER K COMPARISON
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. HIGHER-LAYER COMPARISON")
    print("=" * 78)

    # Compare the two possible one-row reconstructions using C0 or C1
    # at increasing powers of two. This is a direct test of whether their
    # difference becomes relevant after dividing by higher powers.

    if len(L1_support) == 1:

        selected_row, selected_col = (
            L1_support[0]
        )

        for e in range(1, 7):

            modulus = 2 ** e

            c0_residue = [
                x % modulus
                for x in C0
            ]

            c1_residue = [
                x % modulus
                for x in C1
            ]

            different = (
                c0_residue != c1_residue
            )

            print(
                f"  modulus={modulus}: "
                f"C0_equals_C1={not different}"
            )

    # ------------------------------------------------------------------
    # 9. CORE RANK
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
    print("9. CORE RANK REFERENCE")
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
    # 10. INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 171S showed that the two terminal boundary rows have the
same mod-2 C-row and that swapping them changes only the local row
index selected by the first Schur lift.

Experiment 172 isolates the deeper arithmetic.

If

    C0 - C1

has a large common factor 2^e, then the two terminal boundary
equations remain indistinguishable through e 2-adic layers.

If the first nonzero layer occurs immediately at e=1, then their
mod-2 agreement is shallow.

The quantity

    DeltaC = C0 - C1

is also the exact correction that would change the selected-row
reconstruction of K when switching between C0 and C1.

Thus:

    DeltaC mod 2 = 0

explains why the row-order swap leaves K mod 2 unchanged.

The higher-layer audit determines whether that coincidence survives
modulo 4, 8, 16, and beyond.

A high-order equality would indicate a genuine 2-adic closeness of the
two terminal boundary equations. A low-order split would indicate that
the equivalence is purely first-order.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    data_exact = (
        len(C0) == CORE_N
        and len(C1) == CORE_N
    )

    mod2_equal = (
        to_mod2([C0])
        == to_mod2([C1])
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
        data_exact
        and abs(det_H) == 2
        and mod2_equal
        and K_rank_ok
        and F_rank_ok
        and A_rank_ok
        and L1_support == [(11, 0)]
        and K_delta_mod2_zero
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  boundary_rows_exact={data_exact}"
    )

    print(
        f"  abs_det_H_2={abs(det_H) == 2}"
    )

    print(
        f"  C0_equals_C1_mod2={mod2_equal}"
    )

    print(
        f"  L1_expected_support="
        f"{L1_support == [(11, 0)]}"
    )

    print(
        f"  DeltaC_invisible_to_K_mod2="
        f"{K_delta_mod2_zero}"
    )

    print(
        f"  rank_K_one={K_rank_ok}"
    )

    print(
        f"  rank_F_four={F_rank_ok}"
    )

    print(
        f"  rank_A_three={A_rank_ok}"
    )

    print(
        f"  minimum_nonzero_v2_DeltaC="
        f"{min_delta_v2}"
    )

    print(
        f"  first_difference_e="
        f"{first_difference_e}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 172 COMPLETE")


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

