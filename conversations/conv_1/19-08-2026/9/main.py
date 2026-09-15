#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 178 — EXACT TERMINAL SCHUR ROW / COMMON CONTENT / RATIO AUDIT
==============================================================================

Experiment 177 established the exact terminal Schur row

    K_row =
      [s_1, 0, 0, s_p, 0, s_p2, 0, ..., 0]

with

    s_1  = q_1(0) - 2 q_3(0)
    s_p  = q_1(0) - 6 q_3(0)
    s_p2 = q_1(0) - 18 q_3(0).

The exact values are

    s_1  = -347577631
    s_p  = -2033635387
    s_p2 = -7091808655.

All three are odd, so the previous "minimum v2 = 1" expectation
was incorrect. This experiment does not impose that expectation.

Experiment 178 studies the newly observed common content

    gcd(|s_1|, |s_p|, |s_p2|) = 17.

Questions:

    1. Verify the exact gcd.
    2. Remove the common factor 17.
    3. Test whether the primitive coefficients satisfy exact simple
       relations.
    4. Test first and second finite differences in the p-coordinate.
    5. Test whether the three values fit an exact quadratic in p at
       p = 1,3 using the known p^2 structure.
    6. Test ratios and cross-multiplication identities exactly.
    7. Examine prime valuations of the primitive coefficients.
    8. Verify that 17 is the complete integer content.

No interpolation beyond the already-present basis structure is used.

No claim is made for other p-values.
No SymPy.
No floating point.
"""

from __future__ import annotations

import sys
from math import gcd


# ============================================================================
# DATA
# ============================================================================

Q1 = 495451247
Q3 = 421514439

P0 = 1
P1 = 3


# ============================================================================
# HELPERS
# ============================================================================

def v2(x):
    if x == 0:
        return None

    x = abs(x)
    e = 0

    while x % 2 == 0:
        x //= 2
        e += 1

    return e


def vp(x, p):
    if x == 0:
        return None

    x = abs(x)
    e = 0

    while x % p == 0:
        x //= p
        e += 1

    return e


def factor_small_primes(x, primes):
    return {
        p: vp(x, p)
        for p in primes
    }


def exact_ratio(a, b):
    """
    Return the reduced rational ratio a/b as (numerator, denominator).
    """
    g = gcd(abs(a), abs(b))

    if b < 0:
        g = -g

    return a // g, b // g


def pairwise_cross_zero(a, b, c, d):
    return a * d == b * c


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 178 — EXACT TERMINAL SCHUR ROW / "
        "COMMON CONTENT / RATIO AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. EXACT TERMINAL SCHUR COEFFICIENTS
    # ------------------------------------------------------------------

    s1 = Q1 - 2 * Q3
    sp = Q1 - 6 * Q3
    sp2 = Q1 - 18 * Q3

    values = [s1, sp, sp2]

    print()
    print("=" * 78)
    print("1. EXACT TERMINAL SCHUR COEFFICIENTS")
    print("=" * 78)

    print(
        f"  s1={s1}"
    )

    print(
        f"  sp={sp}"
    )

    print(
        f"  sp2={sp2}"
    )

    print(
        f"  v2=[{v2(s1)}, {v2(sp)}, {v2(sp2)}]"
    )

    # ------------------------------------------------------------------
    # 2. EXACT COMMON CONTENT
    # ------------------------------------------------------------------

    g = gcd(abs(s1), abs(sp))
    g = gcd(g, abs(sp2))

    primitive = [
        s1 // g,
        sp // g,
        sp2 // g,
    ]

    primitive_gcd = gcd(
        gcd(abs(primitive[0]), abs(primitive[1])),
        abs(primitive[2]),
    )

    print()
    print("=" * 78)
    print("2. EXACT COMMON CONTENT")
    print("=" * 78)

    print(
        f"  gcd={g}"
    )

    print(
        f"  primitive={primitive}"
    )

    print(
        f"  gcd_after_dividing={primitive_gcd}"
    )

    print(
        f"  gcd_prime_valuation_17={vp(g, 17)}"
    )

    print(
        f"  gcd_prime_valuation_2={vp(g, 2)}"
    )

    # ------------------------------------------------------------------
    # 3. FINITE DIFFERENCE AUDIT
    # ------------------------------------------------------------------

    first_difference = sp - s1
    second_difference = sp2 - sp

    second_difference_again = (
        sp2 - 2 * sp + s1
    )

    print()
    print("=" * 78)
    print("3. EXACT FINITE DIFFERENCE AUDIT")
    print("=" * 78)

    print(
        f"  sp-s1={first_difference}"
    )

    print(
        f"  sp2-sp={second_difference}"
    )

    print(
        f"  second_difference="
        f"{second_difference_again}"
    )

    print(
        f"  gcd(first_differences)="
        f"{gcd(abs(first_difference), abs(second_difference))}"
    )

    # ------------------------------------------------------------------
    # 4. ACCOUNT FOR THE ACTUAL p LOCATIONS
    # ------------------------------------------------------------------

    """
    The active coordinates are not arbitrary values at 1,2,3.

    They are evaluations of the terminal Schur expression at

        p = 1
        p = 3

    with the coefficients

        q1 - 2q3
        q1 - 6q3
        q1 - 18q3.

    We test the exact relation between these three coefficients without
    inventing a new interpolation model.
    """

    # The difference between successive polynomial coefficients in the
    # formal p-pattern is exact.
    delta_p_coeff = sp - s1
    delta_p2_coeff = sp2 - sp

    print()
    print("=" * 78)
    print("4. p-COEFFICIENT PROGRESSION")
    print("=" * 78)

    print(
        f"  delta1={delta_p_coeff}"
    )

    print(
        f"  delta2={delta_p2_coeff}"
    )

    print(
        f"  delta2-delta1="
        f"{delta_p2_coeff-delta_p_coeff}"
    )

    # ------------------------------------------------------------------
    # 5. RATIO AUDIT
    # ------------------------------------------------------------------

    r_sp_s1 = exact_ratio(sp, s1)
    r_sp2_s1 = exact_ratio(sp2, s1)
    r_sp2_sp = exact_ratio(sp2, sp)

    print()
    print("=" * 78)
    print("5. EXACT RATIO AUDIT")
    print("=" * 78)

    print(
        f"  sp/s1={r_sp_s1[0]}/{r_sp_s1[1]}"
    )

    print(
        f"  sp2/s1={r_sp2_s1[0]}/{r_sp2_s1[1]}"
    )

    print(
        f"  sp2/sp={r_sp2_sp[0]}/{r_sp2_sp[1]}"
    )

    print(
        f"  |sp/s1|_denominator={abs(r_sp_s1[1])}"
    )

    print(
        f"  |sp2/s1|_denominator={abs(r_sp2_s1[1])}"
    )

    # ------------------------------------------------------------------
    # 6. CROSS-MULTIPLICATION TESTS
    # ------------------------------------------------------------------

    """
    Rather than comparing floating-point ratios, test exact identities
    of the form

        a/b = c/d.

    Several simple integer candidates are tested.
    """

    ratio_tests = {
        "sp = 2*s1": sp == 2 * s1,
        "sp = 3*s1": sp == 3 * s1,
        "sp2 = 2*sp": sp2 == 2 * sp,
        "sp2 = 3*sp": sp2 == 3 * sp,
        "sp2 = 2*sp - s1": (
            sp2 == 2 * sp - s1
        ),
        "sp2 = 3*sp - 2*s1": (
            sp2 == 3 * sp - 2 * s1
        ),
    }

    print()
    print("=" * 78)
    print("6. SIMPLE EXACT RELATION TESTS")
    print("=" * 78)

    for name, result in ratio_tests.items():

        print(
            f"  {name}: {result}"
        )

    # ------------------------------------------------------------------
    # 7. PRIMITIVE PRIME SPECTRUM
    # ------------------------------------------------------------------

    small_primes = [
        2, 3, 5, 7, 11, 13, 17, 19
    ]

    print()
    print("=" * 78)
    print("7. PRIMITIVE PRIME SPECTRUM")
    print("=" * 78)

    for i, value in enumerate(primitive):

        print(
            f"  primitive_{i}={value} "
            f"valuations="
            f"{factor_small_primes(value, small_primes)}"
        )

    # ------------------------------------------------------------------
    # 8. CONTENT AFTER SELECTED PRIME REMOVAL
    # ------------------------------------------------------------------

    primitive_after_17 = [
        value
        for value in primitive
    ]

    gcd_after_17 = gcd(
        gcd(
            abs(primitive_after_17[0]),
            abs(primitive_after_17[1])
        ),
        abs(primitive_after_17[2])
    )

    print()
    print("=" * 78)
    print("8. CONTENT CONSISTENCY")
    print("=" * 78)

    print(
        f"  gcd_original={g}"
    )

    print(
        f"  gcd_after_primitive_normalization="
        f"{primitive_gcd}"
    )

    print(
        f"  gcd_after_17={gcd_after_17}"
    )

    # ------------------------------------------------------------------
    # 9. EXACT TERMINAL ROW
    # ------------------------------------------------------------------

    K_row = [
        s1,
        0,
        0,
        sp,
        0,
        sp2,
        0,
        0,
        0,
        0,
        0,
        0,
    ]

    primitive_row = [
        x // g
        for x in K_row
    ]

    reconstruction_from_primitive = [
        g * x
        for x in primitive_row
    ]

    reconstruction_exact = (
        reconstruction_from_primitive == K_row
    )

    print()
    print("=" * 78)
    print("9. EXACT ROW CONTENT RECONSTRUCTION")
    print("=" * 78)

    print(
        f"  K_row={K_row}"
    )

    print(
        f"  primitive_row={primitive_row}"
    )

    print(
        f"  g*primitive_row="
        f"{reconstruction_from_primitive}"
    )

    print(
        f"  reconstruction_exact="
        f"{reconstruction_exact}"
    )

    # ------------------------------------------------------------------
    # 10. 2-ADIC CONTENT
    # ------------------------------------------------------------------

    row_v2 = [
        v2(x)
        for x in K_row
        if x != 0
    ]

    row_min_v2 = min(row_v2)

    print()
    print("=" * 78)
    print("10. 2-ADIC CONTENT")
    print("=" * 78)

    print(
        f"  nonzero_row_v2={row_v2}"
    )

    print(
        f"  minimum_row_v2={row_min_v2}"
    )

    print(
        f"  gcd_v2={v2(g)}"
    )

    print(
        f"  odd_content={g//(2**v2(g))}"
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
Experiment 177 showed that all three active terminal Schur
coefficients are already odd, while their common integer content is 17.

Experiment 178 therefore separates two different phenomena:

    2-adic content:
        none, because the row is already odd;

    odd arithmetic content:
        exactly 17.

After division by 17 the three coefficients are primitive.

The finite-difference and ratio tests are intentionally conservative.
They ask whether the primitive triple satisfies an additional exact
low-complexity relation beyond the already-known Schur construction.

A negative result means the factor 17 is common arithmetic content but
does not by itself produce an obvious low-degree relation.

A positive exact relation would identify additional structure in the
terminal Schur row.
"""
    )

    # ------------------------------------------------------------------
    # 12. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        g == 17
        and primitive_gcd == 1
        and reconstruction_exact
        and row_min_v2 == 0
        and v2(g) == 0
    )

    print()
    print("=" * 78)
    print("12. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  gcd_is_17={g == 17}"
    )

    print(
        f"  primitive_content_one={primitive_gcd == 1}"
    )

    print(
        f"  row_reconstruction_exact="
        f"{reconstruction_exact}"
    )

    print(
        f"  row_is_2adic_primitive="
        f"{row_min_v2 == 0}"
    )

    print(
        f"  gcd_is_odd="
        f"{v2(g) == 0}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 178 COMPLETE")


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

