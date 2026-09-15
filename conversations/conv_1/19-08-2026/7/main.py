#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 176 — EXACT SECOND-BOUNDARY 2-ADIC CORRECTION /
                K1-LAYER VS TERMINAL-DIFFERENCE AUDIT
==============================================================================

Experiment 175 established the exact decomposition

    K = K0 + K1

with

    K0 = -a C0
    K1 = -b C1,

and found

    a = -e_11
    b = 2 e_11.

Therefore, in the original ordering,

    K0 = e_11 C0
    K1 = -2 e_11 C1.

Hence

    K mod 2 = C0 inserted in row 11,

while the first mod-4 correction is carried by

    K1 / 2 mod 2.

Experiment 176 isolates that correction layer.

It compares:

    A = (K - K0)/2 mod 2,

    B = K1/2 mod 2,

    C = (C0 - C1)/2 mod 2,

    D = C1 mod 2.

The exact questions are:

    1. Is A exactly B?
    2. Is B exactly -C1 mod 2?
    3. Is A related to (C0-C1)/2?
    4. Which coordinates survive in the second correction?
    5. Does the second correction have rank one?
    6. Does the same relation survive after swapping the two boundary rows?
    7. Is the entire mod-4 discrepancy from K0 accounted for by the
       second terminal boundary alone?

No SymPy.
No floating point.
No extrapolation.
No recurrence fitting.
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


def divide_exact(A, divisor, name):
    out = []

    for i, row in enumerate(A):
        new_row = []

        for j, x in enumerate(row):
            if x % divisor != 0:
                raise ArithmeticError(
                    f"{name}[{i}][{j}]={x} "
                    f"is not divisible by {divisor}."
                )

            new_row.append(x // divisor)

        out.append(new_row)

    return out


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

    if len(M) != FULL_N:
        raise ArithmeticError(
            f"Expected {FULL_N} rows, got {len(M)}."
        )

    return M


# ============================================================================
# SCHUR DATA
# ============================================================================

def build_schur_data(M, boundary_order):

    retained_indices = [
        i
        for i in range(FULL_N)
        if i not in boundary_order
    ]

    retained = [
        M[i] for i in retained_indices
    ]

    boundary = [
        M[i] for i in boundary_order
    ]

    F = [
        [
            row[j] for j in range(CORE_N)
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

    return F, B, H, adj_h, C, det_h


def compute_L1(B, adj_H, det_H):

    numerator = matrix_to_int(
        matmul(B, adj_H),
        "B adj(H)"
    )

    return divide_exact(
        numerator,
        abs(det_H),
        "B adj(H)"
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
            "Schur numerator not divisible by det(H)."
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
        "EXPERIMENT 176 — EXACT SECOND-BOUNDARY 2-ADIC CORRECTION / "
        "K1-LAYER VS TERMINAL-DIFFERENCE AUDIT"
    )
    print("=" * 78)

    M = build_full_system()

    F, B, H, adj_H, C, det_H = build_schur_data(
        M,
        BOUNDARY_ROWS,
    )

    L1 = compute_L1(
        B,
        adj_H,
        det_H,
    )

    K = compute_K(
        B,
        adj_H,
        C,
        det_H,
    )

    K_mod2 = to_mod2(K)

    C0 = C[0]
    C1 = C[1]

    # ------------------------------------------------------------------
    # 1. EXACT SELECTOR
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
    print("1. EXACT SCHUR SELECTOR")
    print("=" * 78)

    print(
        f"  det_H={det_H}"
    )

    print(
        f"  nonzero_a={[i for i, x in enumerate(a) if x != 0]}"
    )

    print(
        f"  nonzero_b={[i for i, x in enumerate(b) if x != 0]}"
    )

    print(
        f"  a={a}"
    )

    print(
        f"  b={b}"
    )

    # ------------------------------------------------------------------
    # 2. EXACT K0/K1
    # ------------------------------------------------------------------

    K0 = [
        [
            -a[i] * C0[j]
            for j in range(CORE_N)
        ]
        for i in range(CORE_N)
    ]

    K1 = [
        [
            -b[i] * C1[j]
            for j in range(CORE_N)
        ]
        for i in range(CORE_N)
    ]

    decomposition_exact = (
        all(
            K[i][j]
            == K0[i][j] + K1[i][j]
            for i in range(CORE_N)
            for j in range(CORE_N)
        )
    )

    print()
    print("=" * 78)
    print("2. EXACT K0 / K1 DECOMPOSITION")
    print("=" * 78)

    print(
        f"  decomposition_exact={decomposition_exact}"
    )

    print(
        f"  K0_mod2_support="
        f"{matrix_support(to_mod2(K0))}"
    )

    print(
        f"  K1_mod2_support="
        f"{matrix_support(to_mod2(K1))}"
    )

    print(
        f"  K_mod2_support="
        f"{matrix_support(K_mod2)}"
    )

    # ------------------------------------------------------------------
    # 3. SECOND-BOUNDARY LAYER
    # ------------------------------------------------------------------

    K_minus_K0 = [
        [
            K[i][j] - K0[i][j]
            for j in range(CORE_N)
        ]
        for i in range(CORE_N)
    ]

    exact_same_K1 = (
        K_minus_K0 == K1
    )

    K1_v2_values = [
        v2(x)
        for row in K1
        for x in row
        if x != 0
    ]

    minimum_K1_v2 = min(
        K1_v2_values
    ) if K1_v2_values else None

    K1_over_2 = None

    if minimum_K1_v2 is not None:

        if minimum_K1_v2 >= 1:

            K1_over_2 = divide_exact(
                K1,
                2,
                "K1"
            )

    print()
    print("=" * 78)
    print("3. FIRST NONZERO K1 LAYER")
    print("=" * 78)

    print(
        f"  K_minus_K0_equals_K1={exact_same_K1}"
    )

    print(
        f"  minimum_nonzero_v2_K1="
        f"{minimum_K1_v2}"
    )

    if K1_over_2 is not None:

        K1_over_2_mod2 = to_mod2(
            K1_over_2
        )

        print(
            f"  K1/2_mod2_support="
            f"{matrix_support(K1_over_2_mod2)}"
        )

        for i, row in enumerate(
            K1_over_2_mod2
        ):

            if any(row):

                print(
                    f"  K1/2_mod2_row_{i}="
                    f"{row}"
                )

    # ------------------------------------------------------------------
    # 4. COMPARE WITH C1
    # ------------------------------------------------------------------

    C1_mod2 = to_mod2([C1])[0]

    expected_K1_layer = [
        [0 for _ in range(CORE_N)]
        for _ in range(CORE_N)
    ]

    # Since K1/2 = -(b/2) C1 and b_11=2,
    # modulo 2 this is C1 in row 11.
    selected_row = None

    for i, x in enumerate(b):
        if x != 0:
            selected_row = i
            break

    C1_prediction_exact = False

    if (
        K1_over_2 is not None
        and selected_row is not None
    ):

        for j in range(CORE_N):

            expected_K1_layer[
                selected_row
            ][j] = C1_mod2[j]

        C1_prediction_exact = (
            to_mod2(K1_over_2)
            == expected_K1_layer
        )

    print()
    print("=" * 78)
    print("4. K1/2 VS SECOND TERMINAL ROW")
    print("=" * 78)

    print(
        f"  selected_row={selected_row}"
    )

    print(
        f"  C1_mod2={C1_mod2}"
    )

    print(
        f"  C1_prediction_exact="
        f"{C1_prediction_exact}"
    )

    # ------------------------------------------------------------------
    # 5. COMPARE WITH C0-C1 FIRST LAYER
    # ------------------------------------------------------------------

    delta_C = [
        C0[j] - C1[j]
        for j in range(CORE_N)
    ]

    delta_C_mod2 = [
        (x // 2) % 2
        if x % 2 == 0
        else None
        for x in delta_C
    ]

    induced_from_delta = [
        [0 for _ in range(CORE_N)]
        for _ in range(CORE_N)
    ]

    if selected_row is not None:

        for j in range(CORE_N):

            if delta_C_mod2[j] is not None:

                induced_from_delta[
                    selected_row
                ][j] = delta_C_mod2[j]

    K1_layer_mod2 = (
        to_mod2(K1_over_2)
        if K1_over_2 is not None
        else None
    )

    delta_prediction_equal = (
        K1_layer_mod2 is not None
        and induced_from_delta == K1_layer_mod2
    )

    print()
    print("=" * 78)
    print("5. K1 LAYER VS DELTA-C LAYER")
    print("=" * 78)

    print(
        f"  DeltaC/2_mod2={delta_C_mod2}"
    )

    print(
        f"  induced_K1_layer="
        f"{induced_from_delta}"
    )

    print(
        f"  equal_to_K1_layer="
        f"{delta_prediction_equal}"
    )

    # ------------------------------------------------------------------
    # 6. MODULUS PROPAGATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. K VERSUS K0 BY MODULUS")
    print("=" * 78)

    first_difference_modulus = None

    for e in range(1, 9):

        modulus = 2 ** e

        equal = (
            to_mod(K, modulus)
            == to_mod(K0, modulus)
        )

        if (
            not equal
            and first_difference_modulus is None
        ):
            first_difference_modulus = modulus

        print(
            f"  modulus={modulus} "
            f"K_equals_K0={equal}"
        )

    # ------------------------------------------------------------------
    # 7. K1 LAYER SUPPORT
    # ------------------------------------------------------------------

    if K1_over_2 is not None:

        K1_layer_mod2 = to_mod2(
            K1_over_2
        )

        print()
        print("=" * 78)
        print("7. SECOND-BOUNDARY SUPPORT")
        print("=" * 78)

        print(
            f"  support="
            f"{matrix_support(K1_layer_mod2)}"
        )

        print(
            f"  rank_mod2="
            f"{rank_mod2(K1_layer_mod2)}"
        )

    # ------------------------------------------------------------------
    # 8. SWAP CONTROL
    # ------------------------------------------------------------------

    swapped_order = (
        BOUNDARY_ROWS[1],
        BOUNDARY_ROWS[0],
    )

    swap = build_schur_data(
        M,
        swapped_order,
    )

    L1_swap = compute_L1(
        swap[1],
        swap[3],
        swap[5],
    )

    K_swap = compute_K(
        swap[1],
        swap[3],
        swap[4],
        swap[5],
    )

    swap_exact = (
        K_swap == K
    )

    print()
    print("=" * 78)
    print("8. BOUNDARY-ORDER SWAP CONTROL")
    print("=" * 78)

    print(
        f"  swapped_det_H={swap[5]}"
    )

    print(
        f"  swapped_L1_mod2_support="
        f"{matrix_support(to_mod2(L1_swap))}"
    )

    print(
        f"  swapped_K_exactly_same="
        f"{swap_exact}"
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
Experiment 175 established that

    K = K0 + K1,

with

    K0 mod 2 = K mod 2,

but

    K != K0 mod 4.

The current experiment isolates exactly what appears at the first
higher 2-adic layer.

Because the second Schur coefficient is

    b_11 = 2,

we have

    K1/2 = -C1

on the distinguished row.

Modulo 2, the sign disappears, so

    K1/2 mod 2 = C1 mod 2.

Thus the first mod-4 discrepancy of K from K0 can be compared directly
with the second terminal boundary row.

The comparison with

    (C0-C1)/2 mod 2

is a separate diagnostic. Equality would mean that the first terminal
row difference and the second-boundary Schur correction encode the same
residue. Inequality would distinguish the two effects.

The result should therefore tell us whether the mod-4 correction is
simply the second boundary row itself or a more complicated consequence
of the pair of terminal equations.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    decomposition_ok = (
        decomposition_exact
    )

    K1_layer_ok = (
        K1_over_2 is not None
    )

    C1_layer_ok = (
        C1_prediction_exact
    )

    final_ok = (
        abs(det_H) == 2
        and decomposition_ok
        and K1_layer_ok
        and C1_layer_ok
        and swap_exact
        and rank_K == 1
        and rank_F == 4
        and rank_A == 3
        and first_difference_modulus == 4
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  abs_det_H_2="
        f"{abs(det_H) == 2}"
    )

    print(
        f"  exact_K0_plus_K1="
        f"{decomposition_ok}"
    )

    print(
        f"  K1_has_first_2adic_layer="
        f"{K1_layer_ok}"
    )

    print(
        f"  K1_over_2_matches_C1_mod2="
        f"{C1_layer_ok}"
    )

    print(
        f"  K_differs_from_K0_first_at_mod4="
        f"{first_difference_modulus == 4}"
    )

    print(
        f"  boundary_swap_preserves_K="
        f"{swap_exact}"
    )

    print(
        f"  K_rank_one="
        f"{rank_K == 1}"
    )

    print(
        f"  F_rank_four="
        f"{rank_F == 4}"
    )

    print(
        f"  A_rank_three="
        f"{rank_A == 3}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 176 COMPLETE")


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

