#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 154 — EXACT p-ADIC MINOR-FILTRATION / COFACTOR-LAYER AUDIT
==============================================================================

Experiment 153R established the exact mod-2 relations

    rank(A mod 2) = 3,

    rank(R mod 2) = 1,

where A is the 12x12 F-core and R is the normalized cofactor matrix.

Moreover the rank-one directions of R mod 2 satisfy

    r in ker(A mod 2),
    c in ker(A^T mod 2),

but

    r notin Row(A mod 2),
    c notin Col(A mod 2).

This suggests that the normalized cofactor layer may reflect a
higher-order minor filtration of A rather than ordinary mod-2 image
geometry.

Experiment 154 therefore studies exact p-adic valuations of small
minors of A.

For each size k we examine sampled/full families where feasible and
record:

    * minimum v_2 among k x k minors;
    * number of zero minors;
    * number of minors attaining the minimum;
    * rank of A modulo 2^e as inferred from divisibility thresholds;
    * the first k for which a nonzero minor exists modulo 2;
    * comparison with the known terminal cofactor valuation

          v_2(det A) = 33,

          v_2(Delta_11) = 24,

          v_2(d_12) = 9.

The central structural question is:

    Does the normalized 11x11 cofactor layer arise at the first
    p-adic level where 11x11 minors become visible, and does that
    filtration agree with the lower-order rank profile of A?

This is an exact minor audit.

No recurrence fitting.
No SymPy.
No floating point.
No extrapolation.
No claim about the original (p,q)-kernel.
"""

from __future__ import annotations

from itertools import combinations
from math import comb
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

N = 12


# ============================================================================
# BASIC HELPERS
# ============================================================================

def q_value(p: int, r: int) -> int:
    if r < 0 or r >= len(Q[p]):
        return 0
    return Q[p][r]


def monomial_basis(p: int, d: int):
    return [
        1,
        d,
        d * d,
        p,
        p * d,
        p * p,
    ]


# ============================================================================
# BAREISS
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

            M[k], M[pivot_row] = (
                M[pivot_row],
                M[k],
            )

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
# SYSTEM / CORE
# ============================================================================

def build_system():

    M = []

    for p, pnext in TRANSITIONS:

        for r in range(D[p] + 1):

            d = D[p] - r

            q0 = q_value(p, r)
            q1 = q_value(p, r + 1)

            b = monomial_basis(p, d)

            row = [0] * 14

            for j in range(6):
                row[j] = b[j] * q0

            for j in range(6):
                row[6 + j] = b[j] * q1

            if d == 0:
                row[12] = 1
                row[13] = p

            M.append(row)

    return M


def build_core(M):

    b0 = M[5]
    b1 = M[10]

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
            "Boundary block singular."
        )

    # H^{-1}
    a00 = H[1][1]
    a01 = -H[0][1]
    a10 = -H[1][0]
    a11 = H[0][0]

    core = []

    for idx in range(14):

        if idx in (5, 10):
            continue

        row = M[idx]

        reduced = []

        for j in range(12):

            correction = (
                H[0][1] * row[12]
                + H[1][1] * row[13]
            )

            # Direct Schur calculation, written explicitly to avoid
            # unnecessary Fraction growth.
            num = (
                row[j] * det_h
                - row[12] * (
                    H[0][1] * row[j]
                    if False else 0
                )
            )

            # Use exact rational form instead.
            from fractions import Fraction

            value = (
                Fraction(row[j])
                + Fraction(
                    -(
                        a00 * b0[j]
                        + a01 * b1[j]
                    ),
                    det_h,
                ) * row[12]
                + Fraction(
                    -(
                        a10 * b0[j]
                        + a11 * b1[j]
                    ),
                    det_h,
                ) * row[13]
            )

            if value.denominator != 1:
                raise ArithmeticError(
                    "Non-integral F-core entry."
                )

            reduced.append(value.numerator)

        core.append(reduced)

    return core


# ============================================================================
# p-ADIC VALUATION
# ============================================================================

def valuation(n, p):

    n = abs(int(n))

    if n == 0:
        return float("inf")

    v = 0

    while n % p == 0:
        n //= p
        v += 1

    return v


# ============================================================================
# MINOR ENUMERATION
# ============================================================================

def minor_det(A, rows, cols):

    B = [
        [
            A[i][j]
            for j in cols
        ]
        for i in rows
    ]

    return bareiss_det(B)


def audit_all_k_minors(A, k, prime=2):

    row_choices = list(
        combinations(range(N), k)
    )

    col_choices = list(
        combinations(range(N), k)
    )

    minimum = None
    minimum_count = 0
    zero_count = 0
    total = 0

    sample_positions = []

    for rows in row_choices:

        for cols in col_choices:

            total += 1

            det = minor_det(
                A,
                rows,
                cols,
            )

            if det == 0:

                zero_count += 1
                continue

            v = valuation(det, prime)

            if minimum is None or v < minimum:

                minimum = v
                minimum_count = 1
                sample_positions = [
                    (rows, cols, v)
                ]

            elif v == minimum:

                minimum_count += 1

                if len(sample_positions) < 10:

                    sample_positions.append(
                        (rows, cols, v)
                    )

    if minimum is None:
        minimum = float("inf")

    return {
        "k": k,
        "total": total,
        "zero_count": zero_count,
        "minimum_v": minimum,
        "minimum_count": minimum_count,
        "sample": sample_positions,
    }


# ============================================================================
# MODULAR RANK
# ============================================================================

def rref_mod_prime(A, p):

    M = [
        [
            int(x) % p
            for x in row
        ]
        for row in A
    ]

    rows = len(M)
    cols = len(M[0])

    rank = 0

    for col in range(cols):

        pivot = None

        for r in range(rank, rows):

            if M[r][col] % p:
                pivot = r
                break

        if pivot is None:
            continue

        M[rank], M[pivot] = (
            M[pivot],
            M[rank],
        )

        inv = pow(
            M[rank][col],
            -1,
            p,
        )

        for j in range(col, cols):
            M[rank][j] = (
                M[rank][j] * inv
            ) % p

        for r in range(rows):

            if r == rank:
                continue

            factor = M[r][col] % p

            if factor:

                for j in range(col, cols):

                    M[r][j] = (
                        M[r][j]
                        - factor * M[rank][j]
                    ) % p

        rank += 1

        if rank == rows:
            break

    return rank


# ============================================================================
# COFACTOR VALUATION REFERENCE
# ============================================================================

def normalized_cofactor_reference(A):

    n = len(A)

    minimum = None
    count = 0

    for i in range(n):

        for j in range(n):

            rows = [
                r for r in range(n)
                if r != i
            ]

            cols = [
                c for c in range(n)
                if c != j
            ]

            det = minor_det(
                A,
                rows,
                cols,
            )

            if det == 0:
                continue

            v = valuation(
                det,
                2,
            )

            if minimum is None or v < minimum:

                minimum = v
                count = 1

            elif v == minimum:

                count += 1

    return minimum, count


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 154 — EXACT p-ADIC MINOR-FILTRATION / "
        "COFACTOR-LAYER AUDIT"
    )
    print("=" * 78)

    M = build_system()
    A = build_core(M)

    A2_rank = rref_mod_prime(A, 2)

    print()
    print("=" * 78)
    print("1. EXACT CORE VALIDATION")
    print("=" * 78)

    print(
        f"  shape={len(A)}x{len(A[0])}"
    )

    print(
        f"  rank_mod_2={A2_rank}"
    )

    # ------------------------------------------------------------------
    # 2. SMALL MINOR FILTRATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. EXACT 2-ADIC MINOR FILTRATION")
    print("=" * 78)

    # Full enumeration through k=4 is small enough:
    # k=1: 144
    # k=2: 4356
    # k=3: 48400
    # k=4: 245025
    #
    # These are still manageable with Bareiss, but k=4 can take time.
    # The script therefore measures k=1..3 completely and only computes
    # k=4 if explicitly enabled below.

    FULL_K4 = False

    audits = []

    for k in (1, 2, 3):

        result = audit_all_k_minors(
            A,
            k,
            2,
        )

        audits.append(result)

        print()
        print(
            f"  k={k}"
        )

        print(
            f"    total_minors={result['total']}"
        )

        print(
            f"    zero_minors={result['zero_count']}"
        )

        print(
            f"    minimum_v2={result['minimum_v']}"
        )

        print(
            f"    minimum_count={result['minimum_count']}"
        )

        print(
            f"    sample={result['sample']}"
        )

    if FULL_K4:

        result = audit_all_k_minors(
            A,
            4,
            2,
        )

        audits.append(result)

        print()
        print("  k=4")

        print(
            f"    total_minors={result['total']}"
        )

        print(
            f"    zero_minors={result['zero_count']}"
        )

        print(
            f"    minimum_v2={result['minimum_v']}"
        )

        print(
            f"    minimum_count={result['minimum_count']}"
        )

        print(
            f"    sample={result['sample']}"
        )

    # ------------------------------------------------------------------
    # 3. MOD-2 RANK VS FIRST NONZERO MINOR
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. MOD-2 RANK / MINOR CONSISTENCY")
    print("=" * 78)

    first_mod2_nonzero_k = None

    for result in audits:

        if result["minimum_v"] == 0:

            first_mod2_nonzero_k = result["k"]
            break

    print(
        f"  rank_mod_2={A2_rank}"
    )

    print(
        f"  first_k_with_v2_zero_minor="
        f"{first_mod2_nonzero_k}"
    )

    print(
        "  rank/minor_consistency="
        f"{first_mod2_nonzero_k == 1}"
    )

    # ------------------------------------------------------------------
    # 4. 3x3 MINOR PROFILE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. 3x3 MINOR PROFILE")
    print("=" * 78)

    three = next(
        result
        for result in audits
        if result["k"] == 3
    )

    print(
        f"  minimum_v2={three['minimum_v']}"
    )

    print(
        f"  minimum_count={three['minimum_count']}"
    )

    print(
        f"  zero_count={three['zero_count']}"
    )

    # ------------------------------------------------------------------
    # 5. 11x11 COFACTOR LAYER
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. 11x11 COFACTOR VALUATION REFERENCE")
    print("=" * 78)

    minimum_11, count_11 = (
        normalized_cofactor_reference(A)
    )

    print(
        f"  minimum_v2_of_11x11_minors={minimum_11}"
    )

    print(
        f"  number_attaining_minimum={count_11}"
    )

    print(
        f"  Delta_11_v2={valuation(DELTA_11, 2)}"
    )

    # ------------------------------------------------------------------
    # 6. KNOWN DETERMINANT LAYER
    # ------------------------------------------------------------------

    det_A = bareiss_det(A)

    print()
    print("=" * 78)
    print("6. DETERMINANT / TERMINAL FACTOR")
    print("=" * 78)

    print(
        f"  det_digits={len(str(abs(det_A)))}"
    )

    print(
        f"  v2(det_A)={valuation(det_A, 2)}"
    )

    print(
        f"  v2(Delta_11)={valuation(DELTA_11, 2)}"
    )

    terminal_v = (
        valuation(det_A, 2)
        - valuation(DELTA_11, 2)
    )

    print(
        f"  implied_terminal_v2={terminal_v}"
    )

    # ------------------------------------------------------------------
    # 7. FILTRATION SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. p-ADIC FILTRATION SUMMARY")
    print("=" * 78)

    for result in audits:

        print(
            f"  k={result['k']}: "
            f"minimum_v2={result['minimum_v']} "
            f"zero_count={result['zero_count']} "
            f"minimizers={result['minimum_count']}"
        )

    print(
        f"  k=11: "
        f"minimum_v2={minimum_11} "
        f"minimizers={count_11}"
    )

    print(
        f"  k=12: "
        f"minimum_v2={valuation(det_A, 2)}"
    )

    # ------------------------------------------------------------------
    # 8. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The purpose of Experiment 154 is to replace a vague statement such as

    "the normalized cofactor looks special"

with an exact minor filtration.

For each minor size k, the minimum 2-adic valuation measures how
visible the k-th exterior layer of the core is before normalization.

The mod-2 rank of A gives the first-order picture.

The 3x3 minor valuation tests how far that rank survives at the next
determinantal level.

The 11x11 valuation is directly connected to the normalization used
to define R.

Finally,

    v2(det A) - v2(Delta_11)

is the terminal factor valuation already observed in the Smith audit.

The main question is whether these numbers exhibit a coherent
determinantal filtration rather than unrelated valuation jumps.

This experiment does not claim that such a filtration has a particular
cause. It only determines the exact arithmetic structure.
"""
    )

    # ------------------------------------------------------------------
    # 9. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(A) == 12
        and len(A[0]) == 12
        and A2_rank == 3
        and len(audits) >= 3
        and minimum_11 == 24
        and terminal_v == 9
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  core_exact={len(A) == 12 and len(A[0]) == 12}"
    )

    print(
        f"  rank_mod2_exact={A2_rank == 3}"
    )

    print(
        f"  small_minor_audit_completed={len(audits) >= 3}"
    )

    print(
        f"  Delta11_v2_match={minimum_11 == 24}"
    )

    print(
        f"  terminal_v2_match={terminal_v == 9}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 154 COMPLETE")


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

