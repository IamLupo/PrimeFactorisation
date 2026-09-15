#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 184 — EXACT TERMINAL 17-REMOVAL / SCHUR-CONTENT COUNTERFACTUAL
==============================================================================

Experiment 183 established:

    the only nonzero q-values divisible by 17 are

        q_1(5),
        q_3(4),

    and every 17-divisible full-system entry is generated from one of
    those two terminal q-values.

The previous script incorrectly counted zero basis values as "17-divisible".
Zeros have undefined p-adic valuation and are not treated as factors here.

Experiment 184 performs an exact counterfactual:

    original source:
        q_1(5), q_3(4);

    modified source:
        q_1(5)/17, q_3(4)/17.

All other q-values are unchanged.

The experiment rebuilds:

    M,
    F,
    B,
    H,
    C,
    K,

for the modified source table.

It then compares:

    * full-system 17-content;
    * C 17-content;
    * K 17-content;
    * exact relation between original K and modified K;
    * exact relation between original C and modified C;
    * whether B and H remain unchanged;
    * whether the modified K becomes primitive at 17.

No determinant calculation is required.
No SymPy.
No floating point.
No extrapolation.
"""

from __future__ import annotations

import sys
from math import gcd


# ============================================================================
# ORIGINAL DATA
# ============================================================================

Q_ORIGINAL = {
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

P = 17


# ============================================================================
# HELPERS
# ============================================================================

def v17(x):
    if x == 0:
        return None

    x = abs(x)
    e = 0

    while x % P == 0:
        x //= P
        e += 1

    return e


def gcd_nonzero(A):
    g = 0

    for row in A:
        for x in row:
            if x != 0:
                g = gcd(g, abs(x))

    return g


def q_value(Q, p, r):
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
                f"{name}: row {i} has width {len(row)}, "
                f"expected {width}."
            )


def matmul(A, B):
    assert_rectangular(A, "left matrix")
    assert_rectangular(B, "right matrix")

    if not A or not B:
        return []

    if len(A[0]) != len(B):
        raise ValueError(
            f"Dimension mismatch: "
            f"{len(A)}x{len(A[0])} times "
            f"{len(B)}x{len(B[0])}."
        )

    out = [
        [0 for _ in range(len(B[0]))]
        for _ in range(len(A))
    ]

    for i in range(len(A)):
        for k in range(len(B)):

            if A[i][k] == 0:
                continue

            for j in range(len(B[0])):
                out[i][j] += A[i][k] * B[k][j]

    return out


def divide_exact_matrix(A, d, name):
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


def scalar_multiple_equal(A, B, c):
    if len(A) != len(B):
        return False

    if A and len(A[0]) != len(B[0]):
        return False

    for i in range(len(A)):
        for j in range(len(A[0])):

            if A[i][j] != c * B[i][j]:
                return False

    return True


# ============================================================================
# FULL SYSTEM
# ============================================================================

def build_full_system(Q):
    M = []

    for p, _ in TRANSITIONS:

        for r in range(D[p] + 1):

            d = D[p] - r

            q0 = q_value(
                Q,
                p,
                r
            )

            q1 = q_value(
                Q,
                p,
                r + 1
            )

            b = basis(
                p,
                d
            )

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

def build_blocks(M):
    retained = [
        i
        for i in range(FULL_N)
        if i not in BOUNDARY_ROWS
    ]

    retained_rows = [
        M[i]
        for i in retained
    ]

    boundary_rows = [
        M[i]
        for i in BOUNDARY_ROWS
    ]

    F = [
        [
            row[j]
            for j in range(CORE_N)
        ]
        for row in retained_rows
    ]

    B = [
        [
            row[j]
            for j in BOUNDARY_COLS
        ]
        for row in retained_rows
    ]

    H = [
        [
            row[j]
            for j in BOUNDARY_COLS
        ]
        for row in boundary_rows
    ]

    C = [
        [
            row[j]
            for j in range(CORE_N)
        ]
        for row in boundary_rows
    ]

    return F, B, H, C


def compute_schur(B, H, C):
    det_H = (
        H[0][0] * H[1][1]
        - H[0][1] * H[1][0]
    )

    if abs(det_H) != 2:
        raise ArithmeticError(
            f"Expected |det(H)|=2, got {det_H}."
        )

    adj_H = [
        [H[1][1], -H[0][1]],
        [-H[1][0], H[0][0]],
    ]

    numerator = matmul(
        B,
        matmul(
            adj_H,
            C
        )
    )

    for row in numerator:
        for x in row:
            if x % det_H != 0:
                raise ArithmeticError(
                    "Schur numerator not divisible by det(H)."
                )

    K = [
        [
            -(x // det_H)
            for x in row
        ]
        for row in numerator
    ]

    return det_H, adj_H, K


# ============================================================================
# BLOCK STATISTICS
# ============================================================================

def block_stats(A):
    vals = [
        x
        for row in A
        for x in row
        if x != 0
    ]

    if not vals:
        return {
            "gcd": 0,
            "min_v17": None,
            "max_v17": None,
            "divisible": 0,
            "nonzero": 0,
        }

    valuations = [
        v17(x)
        for x in vals
    ]

    div_count = sum(
        1
        for x in vals
        if x % P == 0
    )

    return {
        "gcd": gcd_nonzero(A),
        "min_v17": min(valuations),
        "max_v17": max(valuations),
        "divisible": div_count,
        "nonzero": len(vals),
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 184 — EXACT TERMINAL 17-REMOVAL / "
        "SCHUR-CONTENT COUNTERFACTUAL"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. CONSTRUCT MODIFIED SOURCE TABLE
    # ------------------------------------------------------------------

    Q_MODIFIED = {
        p: list(values)
        for p, values in Q_ORIGINAL.items()
    }

    q1_original = Q_ORIGINAL[1][5]
    q3_original = Q_ORIGINAL[3][4]

    if q1_original % 17 != 0:
        raise ArithmeticError(
            "q_1(5) is not divisible by 17."
        )

    if q3_original % 17 != 0:
        raise ArithmeticError(
            "q_3(4) is not divisible by 17."
        )

    Q_MODIFIED[1][5] = q1_original // 17
    Q_MODIFIED[3][4] = q3_original // 17

    # ------------------------------------------------------------------
    # 2. ORIGINAL / MODIFIED TERMINAL SOURCES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. TERMINAL SOURCE COUNTERFACTUAL")
    print("=" * 78)

    print(
        f"  original_q1={q1_original}"
    )

    print(
        f"  modified_q1={Q_MODIFIED[1][5]}"
    )

    print(
        f"  original_q3={q3_original}"
    )

    print(
        f"  modified_q3={Q_MODIFIED[3][4]}"
    )

    print(
        f"  q1_divided_exact="
        f"{q1_original == 17 * Q_MODIFIED[1][5]}"
    )

    print(
        f"  q3_divided_exact="
        f"{q3_original == 17 * Q_MODIFIED[3][4]}"
    )

    # ------------------------------------------------------------------
    # 3. BUILD BOTH SYSTEMS
    # ------------------------------------------------------------------

    M_original = build_full_system(
        Q_ORIGINAL
    )

    M_modified = build_full_system(
        Q_MODIFIED
    )

    F0, B0, H0, C0 = build_blocks(
        M_original
    )

    F1, B1, H1, C1 = build_blocks(
        M_modified
    )

    det0, adj0, K0 = compute_schur(
        B0,
        H0,
        C0
    )

    det1, adj1, K1 = compute_schur(
        B1,
        H1,
        C1
    )

    # ------------------------------------------------------------------
    # 4. UNCHANGED BLOCKS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. UNCHANGED / CHANGED SCHUR BLOCKS")
    print("=" * 78)

    print(
        f"  B_unchanged={B0 == B1}"
    )

    print(
        f"  H_unchanged={H0 == H1}"
    )

    print(
        f"  adj_H_unchanged={adj0 == adj1}"
    )

    print(
        f"  det_H_unchanged={det0 == det1}"
    )

    print(
        f"  C_changed={C0 != C1}"
    )

    print(
        f"  F_changed={F0 != F1}"
    )

    # ------------------------------------------------------------------
    # 5. EXACT C SCALING
    # ------------------------------------------------------------------

    c_terminal_scaling = (
        scalar_multiple_equal(
            C0,
            C1,
            17
        )
    )

    print()
    print("=" * 78)
    print("3. EXACT C SCALING")
    print("=" * 78)

    print(
        f"  C_original=17*C_modified="
        f"{c_terminal_scaling}"
    )

    print(
        f"  gcd(C_original)={block_stats(C0)['gcd']}"
    )

    print(
        f"  gcd(C_modified)={block_stats(C1)['gcd']}"
    )

    # ------------------------------------------------------------------
    # 6. EXACT K SCALING
    # ------------------------------------------------------------------

    k_scaling = (
        scalar_multiple_equal(
            K0,
            K1,
            17
        )
    )

    print()
    print("=" * 78)
    print("4. EXACT K SCALING")
    print("=" * 78)

    print(
        f"  K_original=17*K_modified="
        f"{k_scaling}"
    )

    print(
        f"  gcd(K_original)={block_stats(K0)['gcd']}"
    )

    print(
        f"  gcd(K_modified)={block_stats(K1)['gcd']}"
    )

    # ------------------------------------------------------------------
    # 7. 17-ADIC CONTENT COMPARISON
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. 17-ADIC CONTENT COMPARISON")
    print("=" * 78)

    for name, A in [
        ("M_original", M_original),
        ("M_modified", M_modified),
        ("F_original", F0),
        ("F_modified", F1),
        ("C_original", C0),
        ("C_modified", C1),
        ("K_original", K0),
        ("K_modified", K1),
    ]:

        s = block_stats(A)

        print(
            f"  {name}: "
            f"gcd={s['gcd']} "
            f"min_v17={s['min_v17']} "
            f"max_v17={s['max_v17']} "
            f"divisible={s['divisible']}/{s['nonzero']}"
        )

    # ------------------------------------------------------------------
    # 8. MODIFIED K PRIMITIVITY
    # ------------------------------------------------------------------

    modified_K_stats = block_stats(
        K1
    )

    modified_C_stats = block_stats(
        C1
    )

    modified_K_primitive = (
        modified_K_stats["gcd"] == 1
    )

    modified_C_primitive = (
        modified_C_stats["gcd"] == 1
    )

    print()
    print("=" * 78)
    print("6. MODIFIED SCHUR PRIMITIVITY")
    print("=" * 78)

    print(
        f"  modified_C_primitive="
        f"{modified_C_primitive}"
    )

    print(
        f"  modified_K_primitive="
        f"{modified_K_primitive}"
    )

    # ------------------------------------------------------------------
    # 9. SUPPORT
    # ------------------------------------------------------------------

    support_original_K = [
        (i, j)
        for i in range(FULL_N - 2)
        for j in range(CORE_N)
        if K0[i][j] != 0
    ]

    support_modified_K = [
        (i, j)
        for i in range(FULL_N - 2)
        for j in range(CORE_N)
        if K1[i][j] != 0
    ]

    print()
    print("=" * 78)
    print("7. K SUPPORT COMPARISON")
    print("=" * 78)

    print(
        f"  support_original="
        f"{support_original_K}"
    )

    print(
        f"  support_modified="
        f"{support_modified_K}"
    )

    print(
        f"  support_unchanged="
        f"{support_original_K == support_modified_K}"
    )

    # ------------------------------------------------------------------
    # 10. EXACT TERMINAL ROW
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. DISTINGUISHED TERMINAL K ROW")
    print("=" * 78)

    print(
        f"  original_K_row={K0[11]}"
    )

    print(
        f"  modified_K_row={K1[11]}"
    )

    print(
        f"  exact_row_scaling="
        f"{K0[11] == [17*x for x in K1[11]]}"
    )

    # ------------------------------------------------------------------
    # 11. INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The counterfactual removes the factor 17 only from the two terminal
q-values that were identified in Experiment 183.

Everything else is held fixed.

If the output satisfies

    C_original = 17 C_modified
    K_original = 17 K_modified,

while B and H are unchanged, then the 17-factor is exactly localized
to those two terminal q-values and transferred linearly through the
Schur complement.

If the modified C and K become primitive, then there is no hidden
additional 17-factor produced by the Schur algebra itself.

If F and M also become structurally different only in the rows
containing the modified terminal q-values, the counterfactual gives a
clean source-level localization.

This is a direct intervention on the finite exact data, not a claim
about how future q-values behave.
"""
    )

    # ------------------------------------------------------------------
    # 12. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        c_terminal_scaling
        and k_scaling
        and B0 == B1
        and H0 == H1
        and adj0 == adj1
        and det0 == det1
        and modified_C_primitive
        and modified_K_primitive
        and support_original_K == support_modified_K
    )

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  C_scales_by_17={c_terminal_scaling}"
    )

    print(
        f"  K_scales_by_17={k_scaling}"
    )

    print(
        f"  B_unchanged={B0 == B1}"
    )

    print(
        f"  H_unchanged={H0 == H1}"
    )

    print(
        f"  adj_H_unchanged={adj0 == adj1}"
    )

    print(
        f"  det_H_unchanged={det0 == det1}"
    )

    print(
        f"  modified_C_primitive="
        f"{modified_C_primitive}"
    )

    print(
        f"  modified_K_primitive="
        f"{modified_K_primitive}"
    )

    print(
        f"  K_support_unchanged="
        f"{support_original_K == support_modified_K}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 184 COMPLETE")


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

