#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 140 — EXACT TERMINAL-INVARIANT ISOLATION AUDIT
==============================================================================

This experiment avoids full Smith reduction.

Known exact data from Experiment 138:

    core is 12 x 12 and full rank.

    v_2(det) = 33
    v_3(det) = 6
    v_7(det) = 2

    v_2(Delta_11) = 24
    v_3(Delta_11) = 4
    v_7(Delta_11) = 1

Therefore:

    d_12 = det / Delta_11

has

    v_2(d_12) = 9
    v_3(d_12) = 2
    v_7(d_12) = 1.

The goal here is to verify these terminal valuations independently
using selected maximal minors and exact determinant quotients.

We deliberately do NOT attempt the complete Smith form.

No floating point.
No SymPy.
No extrapolation.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd


# ============================================================================
# EXACT CORE FROM THE ESTABLISHED 12x12 SYSTEM
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


def valuation(n: int, p: int) -> int:
    n = abs(n)

    if n == 0:
        return 10**9

    v = 0

    while n % p == 0:
        n //= p
        v += 1

    return v


def monomial_basis(p: int, d: int):
    return [
        1,
        d,
        d * d,
        p,
        p * d,
        p * p,
    ]


def q_value(p: int, r: int) -> int:
    if r < 0:
        return 0
    if r >= len(Q[p]):
        return 0
    return Q[p][r]


def bareiss_det(A):
    M = [
        [int(x) for x in row]
        for row in A
    ]

    n = len(M)

    if n == 0:
        return 1

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
                    value //= previous

                M[i][j] = value

        for i in range(k + 1, n):
            M[i][k] = 0

        previous = pivot

    return sign * M[-1][-1]


def build_full_system():

    M = []

    for p, p_next in TRANSITIONS:

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


def build_core(M):

    r0 = M[5]
    r1 = M[10]

    H = [
        [r0[12], r0[13]],
        [r1[12], r1[13]],
    ]

    det_h = (
        H[0][0] * H[1][1]
        - H[0][1] * H[1][0]
    )

    inv_h = [
        [
            Fraction(H[1][1], det_h),
            Fraction(-H[0][1], det_h),
        ],
        [
            Fraction(-H[1][0], det_h),
            Fraction(H[0][0], det_h),
        ],
    ]

    rhs_h = [
        q_value(3, 5),
        q_value(5, 4),
    ]

    h_const = [
        inv_h[i][0] * rhs_h[0]
        + inv_h[i][1] * rhs_h[1]
        for i in range(2)
    ]

    h_coeff = [
        [Fraction(0) for _ in range(12)]
        for _ in range(2)
    ]

    for i in range(2):
        for j in range(12):
            h_coeff[i][j] = -(
                inv_h[i][0] * r0[j]
                + inv_h[i][1] * r1[j]
            )

    core = []

    for idx in range(14):

        if idx in (5, 10):
            continue

        row = M[idx]

        if idx < 6:
            p = 1
            r = idx
        elif idx < 11:
            p = 3
            r = idx - 6
        else:
            p = 5
            r = idx - 11

        p_next = p + 2

        out = []

        for j in range(12):

            value = (
                Fraction(row[j])
                + Fraction(row[12]) * h_coeff[0][j]
                + Fraction(row[13]) * h_coeff[1][j]
            )

            if value.denominator != 1:
                raise ArithmeticError(
                    "Non-integral core entry."
                )

            out.append(value.numerator)

        core.append(out)

    return core


def minor(A, omit_row, omit_col):

    rows = [
        i for i in range(len(A))
        if i != omit_row
    ]

    cols = [
        j for j in range(len(A[0]))
        if j != omit_col
    ]

    return [
        [A[i][j] for j in cols]
        for i in rows
    ]


def gcd_maximal_minors(A):

    g = 0
    n = len(A)

    for i in range(n):

        for j in range(n):

            value = bareiss_det(
                minor(A, i, j)
            )

            g = gcd(g, abs(value))

            if g == 1:
                return 1

    return g


def main():

    print("=" * 78)
    print("EXPERIMENT 140 — EXACT TERMINAL-INVARIANT ISOLATION AUDIT")
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. BUILD
    # ------------------------------------------------------------------

    M = build_full_system()
    core = build_core(M)

    core_det = abs(bareiss_det(core))
    delta11 = gcd_maximal_minors(core)

    print()
    print("=" * 78)
    print("1. EXACT CORE")
    print("=" * 78)

    print("  shape=(12,12)")
    print(f"  det_digits={len(str(core_det))}")
    print(f"  Delta_11={delta11}")

    # ------------------------------------------------------------------
    # 2. EXACT VALUATIONS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. DETERMINANT / DELTA-11 / TERMINAL VALUATIONS")
    print("=" * 78)

    terminal = core_det // delta11

    for p in (2, 3, 7):

        v_det = valuation(core_det, p)
        v_delta = valuation(delta11, p)
        v_terminal = valuation(terminal, p)

        print()
        print(f"  prime={p}")
        print(f"    v(det)     = {v_det}")
        print(f"    v(Delta11) = {v_delta}")
        print(f"    v(d12)     = {v_terminal}")
        print(
            f"    exact_sum = "
            f"{v_delta + v_terminal == v_det}"
        )

    # ------------------------------------------------------------------
    # 3. TERMINAL FACTOR COPRIMALITY TEST
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. TERMINAL FACTOR COPRIMALITY")
    print("=" * 78)

    g = gcd(delta11, terminal)

    print(f"  gcd(Delta_11, d_12)={g}")

    # Prime-wise overlap is more informative.
    for p in (2, 3, 7):

        overlap = min(
            valuation(delta11, p),
            valuation(terminal, p),
        )

        print(
            f"  prime={p}: "
            f"shared_valuation={overlap}"
        )

    # ------------------------------------------------------------------
    # 4. SELECTED MAXIMAL MINORS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. SELECTED 11x11 MINOR VALUATIONS")
    print("=" * 78)

    sample_positions = [
        (0, 0),
        (0, 11),
        (5, 5),
        (5, 11),
        (11, 0),
        (11, 11),
    ]

    for i, j in sample_positions:

        value = abs(
            bareiss_det(
                minor(core, i, j)
            )
        )

        print(
            f"  omit_row={i}, omit_col={j}: "
            f"digits={len(str(value))}, "
            f"v2={valuation(value,2)}, "
            f"v3={valuation(value,3)}, "
            f"v7={valuation(value,7)}"
        )

    # ------------------------------------------------------------------
    # 5. TERMINAL FACTOR RECONSTRUCTION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. TERMINAL FACTOR RECONSTRUCTION")
    print("=" * 78)

    expected_delta = 9512681472

    expected_d12 = (
        22698304379332609168018772364495783974334437093962503157858693669865148581087339234303140035372612096
    )

    delta_match = (
        delta11 == expected_delta
    )

    terminal_match = (
        terminal == expected_d12
    )

    print(f"  Delta_11_match={delta_match}")
    print(f"  d_12_match={terminal_match}")

    # ------------------------------------------------------------------
    # 6. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
This experiment deliberately stops one level before the complete
Smith spectrum.

The exact questions are:

    * Is the 101-digit terminal factor obtained independently?
    * How is its prime valuation split from Delta_11?
    * Does the terminal factor share substantial prime content
      with Delta_11?
    * Do selected maximal minors exhibit the same prime spectrum?

The calculation does not try to explain the large factor.

It only isolates it arithmetically and verifies that the result is
not an artifact of a previous Smith-form computation.

No claim about the original (p,q)-kernel is made.
"""
    )

    # ------------------------------------------------------------------
    # 7. FINAL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. FINAL EXACTNESS")
    print("=" * 78)

    all_prime_sums = all(
        valuation(delta11, p)
        + valuation(terminal, p)
        == valuation(core_det, p)
        for p in (2, 3, 7)
    )

    final_ok = (
        core_det != 0
        and delta11 != 0
        and terminal != 0
        and delta_match
        and terminal_match
        and all_prime_sums
    )

    print(f"  core_exact=True")
    print(f"  Delta_11_exact={delta_match}")
    print(f"  terminal_d12_exact={terminal_match}")
    print(f"  valuation_split_exact={all_prime_sums}")
    print(f"  failures={0 if final_ok else 1}")
    print(f"  ALL BASIC CHECKS PASS={final_ok}")

    print()
    print("EXPERIMENT 140 COMPLETE")


if __name__ == "__main__":
    main()

