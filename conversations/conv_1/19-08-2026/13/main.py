#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 182 — EXACT 17-ADIC SCHUR-CONTENT TRANSFER AUDIT
==============================================================================

Experiment 181 established:

    content(M)      = 1
    content(F_ret)  = 1
    content(B)      = 1
    content(H)      = 1

while

    content(C) = 17
    content(K) = 17.

Moreover:

    v17(Delta_11) = 0
    v17(d12)      = 0.

Experiment 182 tests whether the factor 17 is transferred through the
Schur map exactly.

If

    C = 17 C',

then, because

    K = -B H^{-1} C,

we should have exactly

    K = 17 K'

with

    K' = -B H^{-1} C'.

This is stronger than observing gcd(K)=17.

The experiment therefore performs:

    1. exact divisibility C/17;
    2. exact divisibility K/17;
    3. reconstruction K = 17 K';
    4. independent Schur recomputation from C/17;
    5. entrywise 17-adic valuation comparison;
    6. support comparison between K and the recomputed K';
    7. verification that no additional factor 17 is introduced by
       the B,H side.

No SymPy.
No floating point.
No extrapolation.
"""

from __future__ import annotations

import sys
from math import gcd


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

P17 = 17


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
                f"{name}: row {i} width {len(row)}, "
                f"expected {width}."
            )


def v_p(x, p):
    if x == 0:
        return None

    x = abs(x)
    e = 0

    while x % p == 0:
        x //= p
        e += 1

    return e


def flatten(A):
    return [x for row in A for x in row]


def gcd_nonzero(A):
    g = 0

    for x in flatten(A):
        if x != 0:
            g = gcd(g, abs(x))

    return g


def matrix_support(A):
    return [
        (i, j)
        for i in range(len(A))
        for j in range(len(A[0]))
        if A[i][j] != 0
    ]


def divide_exact_matrix(A, d, name):
    out = []

    for i, row in enumerate(A):
        new_row = []

        for j, x in enumerate(row):

            if x % d != 0:
                raise ArithmeticError(
                    f"{name}[{i}][{j}]={x} "
                    f"is not divisible by {d}."
                )

            new_row.append(x // d)

        out.append(new_row)

    return out


def scalar_multiple_equal(A, B, c):
    assert_rectangular(A, "A")
    assert_rectangular(B, "B")

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
# MATRIX MULTIPLICATION
# ============================================================================

def matmul(A, B):
    assert_rectangular(A, "left")
    assert_rectangular(B, "right")

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
                out[i][j] += (
                    A[i][k] * B[k][j]
                )

    return out


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

def build_blocks(M):

    retained = [
        i for i in range(FULL_N)
        if i not in BOUNDARY_ROWS
    ]

    retained_rows = [
        M[i] for i in retained
    ]

    boundary_rows = [
        M[i] for i in BOUNDARY_ROWS
    ]

    F = [
        [row[j] for j in range(CORE_N)]
        for row in retained_rows
    ]

    B = [
        [row[j] for j in BOUNDARY_COLS]
        for row in retained_rows
    ]

    H = [
        [row[j] for j in BOUNDARY_COLS]
        for row in boundary_rows
    ]

    C = [
        [row[j] for j in range(CORE_N)]
        for row in boundary_rows
    ]

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

    return F, B, H, adj_H, C, det_H


def compute_schur_from_C(
    B,
    adj_H,
    C,
    det_H,
):
    numerator = matmul(
        B,
        matmul(
            adj_H,
            C,
        ),
    )

    for row in numerator:
        for x in row:
            if x % det_H != 0:
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
        "EXPERIMENT 182 — EXACT 17-ADIC SCHUR-CONTENT "
        "TRANSFER AUDIT"
    )
    print("=" * 78)

    M = build_full_system()

    F, B, H, adj_H, C, det_H = build_blocks(M)

    K = compute_schur_from_C(
        B,
        adj_H,
        C,
        det_H,
    )

    # ------------------------------------------------------------------
    # 1. RAW CONTENT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. RAW 17-ADIC CONTENT")
    print("=" * 78)

    gcd_C = gcd_nonzero(C)
    gcd_K = gcd_nonzero(K)

    print(
        f"  gcd(C)={gcd_C}"
    )

    print(
        f"  v17(gcd(C))={v_p(gcd_C,17)}"
    )

    print(
        f"  gcd(K)={gcd_K}"
    )

    print(
        f"  v17(gcd(K))={v_p(gcd_K,17)}"
    )

    # ------------------------------------------------------------------
    # 2. EXACT C/17
    # ------------------------------------------------------------------

    C17_exact = (
        gcd_C % 17 == 0
    )

    C17 = (
        divide_exact_matrix(
            C,
            17,
            "C"
        )
        if C17_exact
        else None
    )

    print()
    print("=" * 78)
    print("2. EXACT C/17")
    print("=" * 78)

    print(
        f"  C_entrywise_divisible_by_17={C17_exact}"
    )

    if C17 is not None:

        print(
            f"  gcd(C/17)={gcd_nonzero(C17)}"
        )

        print(
            f"  C17_support={matrix_support(C17)}"
        )

    # ------------------------------------------------------------------
    # 3. EXACT K/17
    # ------------------------------------------------------------------

    K17_exact = (
        gcd_K % 17 == 0
    )

    K17 = (
        divide_exact_matrix(
            K,
            17,
            "K"
        )
        if K17_exact
        else None
    )

    print()
    print("=" * 78)
    print("3. EXACT K/17")
    print("=" * 78)

    print(
        f"  K_entrywise_divisible_by_17={K17_exact}"
    )

    if K17 is not None:

        print(
            f"  gcd(K/17)={gcd_nonzero(K17)}"
        )

        print(
            f"  K17_support={matrix_support(K17)}"
        )

    # ------------------------------------------------------------------
    # 4. INDEPENDENT SCHUR RECOMPUTATION
    # ------------------------------------------------------------------

    if C17 is not None:

        K_from_C17 = compute_schur_from_C(
            B,
            adj_H,
            C17,
            det_H,
        )

        recomputation_exact = (
            K_from_C17 == K17
        )

    else:

        K_from_C17 = None
        recomputation_exact = False

    print()
    print("=" * 78)
    print("4. SCHUR RECOMPUTATION FROM C/17")
    print("=" * 78)

    print(
        f"  recomputation_exact="
        f"{recomputation_exact}"
    )

    if K_from_C17 is not None:

        print(
            f"  K_from_C17_support="
            f"{matrix_support(K_from_C17)}"
        )

    # ------------------------------------------------------------------
    # 5. ENTRYWISE 17-ADIC VALUATIONS
    # ------------------------------------------------------------------

    c_support = matrix_support(C)
    k_support = matrix_support(K)

    valuation_pairs = []

    for i, j in k_support:

        # Find corresponding source entries contributing through
        # the nonzero terminal rows.
        c0_v = v_p(C[0][j], 17)
        c1_v = v_p(C[1][j], 17)
        k_v = v_p(K[i][j], 17)

        valuation_pairs.append(
            (
                i,
                j,
                c0_v,
                c1_v,
                k_v,
            )
        )

    print()
    print("=" * 78)
    print("5. ENTRYWISE 17-ADIC VALUATION COMPARISON")
    print("=" * 78)

    for item in valuation_pairs:

        print(
            f"  (row={item[0]}, col={item[1]}): "
            f"v17(C0)={item[2]} "
            f"v17(C1)={item[3]} "
            f"v17(K)={item[4]}"
        )

    # ------------------------------------------------------------------
    # 6. MINIMUM K VALUATION
    # ------------------------------------------------------------------

    k_nonzero_values = [
        x
        for row in K
        for x in row
        if x != 0
    ]

    c_nonzero_values = [
        x
        for row in C
        for x in row
        if x != 0
    ]

    min_v_C = min(
        v_p(x, 17)
        for x in c_nonzero_values
    )

    min_v_K = min(
        v_p(x, 17)
        for x in k_nonzero_values
    )

    print()
    print("=" * 78)
    print("6. MINIMUM 17-ADIC ORDERS")
    print("=" * 78)

    print(
        f"  min_v17(C)={min_v_C}"
    )

    print(
        f"  min_v17(K)={min_v_K}"
    )

    # ------------------------------------------------------------------
    # 7. NO EXTRA 17 FROM B/H
    # ------------------------------------------------------------------

    gcd_B = gcd_nonzero(B)
    gcd_H = gcd_nonzero(H)

    print()
    print("=" * 78)
    print("7. B / H 17-ADIC CONTROL")
    print("=" * 78)

    print(
        f"  gcd(B)={gcd_B}"
    )

    print(
        f"  v17(gcd(B))={v_p(gcd_B,17)}"
    )

    print(
        f"  gcd(H)={gcd_H}"
    )

    print(
        f"  v17(gcd(H))={v_p(gcd_H,17)}"
    )

    # ------------------------------------------------------------------
    # 8. SCHUR CONTENT TRANSFER IDENTITY
    # ------------------------------------------------------------------

    content_transfer_exact = (
        C17 is not None
        and K17 is not None
        and K17 == compute_schur_from_C(
            B,
            adj_H,
            C17,
            det_H,
        )
    )

    exact_scalar_transfer = (
        K17_exact
        and content_transfer_exact
        and gcd_nonzero(C17) == 1
        and gcd_nonzero(K17) == 1
    )

    print()
    print("=" * 78)
    print("8. CONTENT-TRANSFER IDENTITY")
    print("=" * 78)

    print(
        f"  K=17*K17_exact="
        f"{scalar_multiple_equal(K, K17, 17) if K17 is not None else False}"
    )

    print(
        f"  Schur(C/17)=K/17="
        f"{content_transfer_exact}"
    )

    print(
        f"  C17_primitive="
        f"{gcd_nonzero(C17) == 1 if C17 is not None else False}"
    )

    print(
        f"  K17_primitive="
        f"{gcd_nonzero(K17) == 1 if K17 is not None else False}"
    )

    # ------------------------------------------------------------------
    # 9. SMITH REFERENCES
    # ------------------------------------------------------------------

    Delta_11 = 9512681472

    d12 = (
        22698304379332609168018772364495783974334437093962503157858693669865148581087339234303140035372612096
    )

    print()
    print("=" * 78)
    print("9. SMITH REFERENCE")
    print("=" * 78)

    print(
        f"  v17(Delta_11)={v_p(Delta_11,17)}"
    )

    print(
        f"  v17(d12)={v_p(d12,17)}"
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
Experiment 181 localized the factor 17:

    C has content 17,
    K has content 17,

while B and H have no 17-content.

Experiment 182 asks whether this is exactly functorial under the Schur
map.

Because

    K = -B H^{-1} C,

integer scaling of C by 17 should scale K by 17 exactly, provided the
Schur division remains integral.

The strongest positive outcome is:

    C = 17 C',
    K = 17 K',
    K' = -B H^{-1} C',

with both C' and K' primitive.

If that holds, the 17-factor is not created by cancellation in the
Schur computation. It is transmitted linearly from the terminal C
block.

The fact that

    v17(Delta_11)=0
    v17(d12)=0

would then further separate this boundary content from the Smith
determinantal invariants.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        C17_exact
        and K17_exact
        and content_transfer_exact
        and gcd_nonzero(C17) == 1
        and gcd_nonzero(K17) == 1
        and min_v_C == 1
        and min_v_K == 1
        and v_p(gcd_B,17) == 0
        and v_p(gcd_H,17) == 0
        and v_p(Delta_11,17) == 0
        and v_p(d12,17) == 0
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  C_divisible_by_17={C17_exact}"
    )

    print(
        f"  K_divisible_by_17={K17_exact}"
    )

    print(
        f"  exact_Schur_content_transfer="
        f"{content_transfer_exact}"
    )

    print(
        f"  C17_primitive="
        f"{gcd_nonzero(C17) == 1 if C17 is not None else False}"
    )

    print(
        f"  K17_primitive="
        f"{gcd_nonzero(K17) == 1 if K17 is not None else False}"
    )

    print(
        f"  min_v17_C_is_1={min_v_C == 1}"
    )

    print(
        f"  min_v17_K_is_1={min_v_K == 1}"
    )

    print(
        f"  B_has_no_17_content="
        f"{v_p(gcd_B,17) == 0}"
    )

    print(
        f"  H_has_no_17_content="
        f"{v_p(gcd_H,17) == 0}"
    )

    print(
        f"  Delta11_has_no_17={v_p(Delta_11,17) == 0}"
    )

    print(
        f"  d12_has_no_17={v_p(d12,17) == 0}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 182 COMPLETE")


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

