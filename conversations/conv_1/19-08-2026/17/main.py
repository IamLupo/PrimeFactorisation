#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 186 — EXACT 3-ADIC TERMINAL SCHUR LAYER /
                CONGRUENCE-STABILITY AUDIT
==============================================================================

Experiment 185 established:

    v3(q1) = 0
    v3(q3) = 1

after removing the localized factor 17, while

    v3(s0) = v3(s1) = v3(s2) = 0.

Thus the q-source valuation asymmetry is not inherited entrywise by the
Schur coefficients.

Experiment 186 studies the row itself modulo increasing powers of 3.

The exact terminal coefficients are

    s_k = q1 - 2*3^k*q3,
    k = 0,1,2.

Since

    v3(q3)=1,

the correction term

    2*3^k*q3

has valuation

    k+1.

Therefore one expects the following exact congruence hierarchy:

    s0 != s1       generally modulo 9,
    s1 == s2       modulo 27,

and more generally differences gain predictable powers of 3.

The experiment verifies directly:

    s1-s0 = -4 q3
    s2-s1 = -12 q3

and checks equality of the three Schur coefficients modulo

    3, 9, 27, 81, 243, 729,...

It also records pairwise difference valuations and identifies the
precision at which each adjacent pair first separates.

No extrapolation beyond the exact k=0,1,2 data is used.
"""

from __future__ import annotations

import sys


# ============================================================================
# SOURCE DATA
# ============================================================================

Q1 = 29144191
Q3 = 24794967


# ============================================================================
# HELPERS
# ============================================================================

def v3(x):
    if x == 0:
        return None

    x = abs(x)
    e = 0

    while x % 3 == 0:
        x //= 3
        e += 1

    return e


def residues(x, prime=3, max_e=8):
    return [
        (
            prime ** e,
            x % (prime ** e),
        )
        for e in range(1, max_e + 1)
    ]


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 186 — EXACT 3-ADIC TERMINAL SCHUR LAYER / "
        "CONGRUENCE-STABILITY AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. EXACT SCHUR ROW
    # ------------------------------------------------------------------

    s0 = Q1 - 2 * Q3
    s1 = Q1 - 6 * Q3
    s2 = Q1 - 18 * Q3

    S = [s0, s1, s2]

    print()
    print("=" * 78)
    print("1. EXACT TERMINAL SCHUR ROW")
    print("=" * 78)

    for i, s in enumerate(S):

        print(
            f"  s{i}={s} "
            f"v3={v3(s)} "
            f"mod3={s % 3}"
        )

    # ------------------------------------------------------------------
    # 2. ADJACENT DIFFERENCES
    # ------------------------------------------------------------------

    d01 = s1 - s0
    d12 = s2 - s1

    print()
    print("=" * 78)
    print("2. EXACT ADJACENT DIFFERENCES")
    print("=" * 78)

    print(
        f"  s1-s0={d01}"
    )

    print(
        f"  v3(s1-s0)={v3(d01)}"
    )

    print(
        f"  s2-s1={d12}"
    )

    print(
        f"  v3(s2-s1)={v3(d12)}"
    )

    print(
        f"  d01/(-4*q3)="
        f"{d01 == -4 * Q3}"
    )

    print(
        f"  d12/(-12*q3)="
        f"{d12 == -12 * Q3}"
    )

    # ------------------------------------------------------------------
    # 3. NON-ADJACENT DIFFERENCE
    # ------------------------------------------------------------------

    d02 = s2 - s0

    print()
    print("=" * 78)
    print("3. NON-ADJACENT DIFFERENCE")
    print("=" * 78)

    print(
        f"  s2-s0={d02}"
    )

    print(
        f"  v3(s2-s0)={v3(d02)}"
    )

    print(
        f"  d02/(-16*q3)="
        f"{d02 == -16 * Q3}"
    )

    # ------------------------------------------------------------------
    # 4. INCREASING 3-ADIC PRECISION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. PAIRWISE CONGRUENCE BY 3-ADIC PRECISION")
    print("=" * 78)

    for e in range(1, 9):

        modulus = 3 ** e

        eq01 = (
            s0 % modulus
            == s1 % modulus
        )

        eq12 = (
            s1 % modulus
            == s2 % modulus
        )

        eq02 = (
            s0 % modulus
            == s2 % modulus
        )

        all_equal = (
            eq01 and eq12 and eq02
        )

        print(
            f"  modulus={modulus}: "
            f"s0=s1 {eq01} "
            f"s1=s2 {eq12} "
            f"s0=s2 {eq02} "
            f"all_equal={all_equal}"
        )

    # ------------------------------------------------------------------
    # 5. FIRST SEPARATION PRECISION
    # ------------------------------------------------------------------

    def first_difference(a, b):

        for e in range(1, 12):

            modulus = 3 ** e

            if a % modulus != b % modulus:
                return e

        return None

    first01 = first_difference(s0, s1)
    first12 = first_difference(s1, s2)
    first02 = first_difference(s0, s2)

    print()
    print("=" * 78)
    print("5. FIRST 3-ADIC SEPARATION")
    print("=" * 78)

    print(
        f"  first_s0_s1_difference_e={first01}"
    )

    print(
        f"  first_s1_s2_difference_e={first12}"
    )

    print(
        f"  first_s0_s2_difference_e={first02}"
    )

    # ------------------------------------------------------------------
    # 6. EXACT RESIDUE PROFILES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. 3-ADIC RESIDUE PROFILES")
    print("=" * 78)

    for i, s in enumerate(S):

        print(
            f"  s{i}: {residues(s, 3, 8)}"
        )

    # ------------------------------------------------------------------
    # 7. SOURCE CORRECTION ORDERS
    # ------------------------------------------------------------------

    correction0 = 2 * Q3
    correction1 = 2 * 3 * Q3
    correction2 = 2 * 9 * Q3

    print()
    print("=" * 78)
    print("7. EXACT SOURCE CORRECTION ORDERS")
    print("=" * 78)

    print(
        f"  correction_k0={correction0} "
        f"v3={v3(correction0)}"
    )

    print(
        f"  correction_k1={correction1} "
        f"v3={v3(correction1)}"
    )

    print(
        f"  correction_k2={correction2} "
        f"v3={v3(correction2)}"
    )

    # ------------------------------------------------------------------
    # 8. SOURCE VS SCHUR LAYER COMPARISON
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. SOURCE / SCHUR LAYER COMPARISON")
    print("=" * 78)

    print(
        f"  v3(q1)={v3(Q1)}"
    )

    print(
        f"  v3(q3)={v3(Q3)}"
    )

    print(
        f"  v3(correction_k0)={v3(correction0)}"
    )

    print(
        f"  v3(correction_k1)={v3(correction1)}"
    )

    print(
        f"  v3(correction_k2)={v3(correction2)}"
    )

    print(
        f"  v3(s0)={v3(s0)}"
    )

    print(
        f"  v3(s1)={v3(s1)}"
    )

    print(
        f"  v3(s2)={v3(s2)}"
    )

    # ------------------------------------------------------------------
    # 9. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The primitive terminal source pair has

    v3(q1)=0,
    v3(q3)=1.

The Schur coefficients are

    s_k=q1-2*3^k*q3.

Therefore the successive correction terms have valuations

    k+1,

for k=0,1,2.

The exact pairwise difference valuations reveal how many 3-adic digits
of the Schur coefficients agree before the geometric correction becomes
visible.

The experiment distinguishes:

    source valuation,
    correction valuation,
    resulting Schur valuation.

In the present dataset the Schur coefficients remain 3-adic units,
even though their differences acquire increasing powers of 3.

That is a stronger statement than simply saying the coefficients are
all nonzero modulo 3.
"""
    )

    # ------------------------------------------------------------------
    # 10. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        v3(Q1) == 0
        and v3(Q3) == 1
        and v3(s0) == 0
        and v3(s1) == 0
        and v3(s2) == 0
        and v3(d01) == 1
        and v3(d12) == 2
        and v3(d02) == 1
        and d01 == -4 * Q3
        and d12 == -12 * Q3
        and d02 == -16 * Q3
        and first12 == 2
    )

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  q1_is_3adic_unit="
        f"{v3(Q1) == 0}"
    )

    print(
        f"  q3_has_exact_v3_1="
        f"{v3(Q3) == 1}"
    )

    print(
        f"  all_schur_coefficients_3adic_units="
        f"{v3(s0) == v3(s1) == v3(s2) == 0}"
    )

    print(
        f"  adjacent_difference_orders="
        f"{[v3(d01), v3(d12)]}"
    )

    print(
        f"  nonadjacent_difference_order="
        f"{v3(d02)}"
    )

    print(
        f"  exact_difference_formulas="
        f"{d01 == -4*Q3 and d12 == -12*Q3 and d02 == -16*Q3}"
    )

    print(
        f"  s1_s2_first_separation_is_mod9="
        f"{first12 == 2}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 186 COMPLETE")


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

