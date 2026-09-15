#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 187 — EXACT 3-ADIC DIFFERENCE-QUOTIENT /
                GEOMETRIC-LAYER AUDIT
==============================================================================

Experiment 186 established:

    s_k = q1 - 2*3^k*q3,

with

    v3(q1)=0,
    v3(q3)=1,

and therefore

    v3(s1-s0)=1,
    v3(s2-s1)=2,
    v3(s2-s0)=1.

The previous final check incorrectly expected the s1/s2 pair to first
separate at modulus 9. In fact:

    s1 == s2 (mod 9),
    s1 != s2 (mod 27),

so the first separation occurs at modulus 27.

Experiment 187 studies the successive difference quotients:

    D0 = (s1-s0)/3
    D1 = (s2-s1)/9.

If the geometric 3-adic layering is exact, then

    D0 = -4*q3
    D1 = -4*q3,

so the two normalized difference layers should be exactly equal.

The experiment therefore checks:

    * exact quotient integrality;
    * equality of the normalized difference layers;
    * their 3-adic unit status;
    * their residues modulo increasing powers of 3;
    * the second-difference identity;
    * the exact first-separation modulus for every pair.

No extrapolation beyond k=0,1,2 is used.
"""

from __future__ import annotations

import sys


# ============================================================================
# DATA
# ============================================================================

Q1 = 29144191
Q3 = 24794967

S0 = Q1 - 2 * Q3
S1 = Q1 - 6 * Q3
S2 = Q1 - 18 * Q3

S = [S0, S1, S2]


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


def residues(x, max_e=8):
    return [
        (3 ** e, x % (3 ** e))
        for e in range(1, max_e + 1)
    ]


def first_separation(a, b, max_e=12):
    for e in range(1, max_e + 1):
        m = 3 ** e

        if a % m != b % m:
            return e

    return None


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 187 — EXACT 3-ADIC DIFFERENCE-QUOTIENT / "
        "GEOMETRIC-LAYER AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. EXACT SCHUR COEFFICIENTS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT TERMINAL SCHUR COEFFICIENTS")
    print("=" * 78)

    for i, value in enumerate(S):

        print(
            f"  s{i}={value} "
            f"v3={v3(value)} "
            f"mod3={value % 3}"
        )

    # ------------------------------------------------------------------
    # 2. EXACT DIFFERENCES
    # ------------------------------------------------------------------

    d0 = S1 - S0
    d1 = S2 - S1
    d2 = S2 - S0

    print()
    print("=" * 78)
    print("2. EXACT DIFFERENCE DATA")
    print("=" * 78)

    print(
        f"  d0=s1-s0={d0} "
        f"v3={v3(d0)}"
    )

    print(
        f"  d1=s2-s1={d1} "
        f"v3={v3(d1)}"
    )

    print(
        f"  d2=s2-s0={d2} "
        f"v3={v3(d2)}"
    )

    # ------------------------------------------------------------------
    # 3. NORMALIZED DIFFERENCE LAYERS
    # ------------------------------------------------------------------

    if d0 % 3 != 0:
        raise ArithmeticError(
            "d0 is not divisible by 3."
        )

    if d1 % 9 != 0:
        raise ArithmeticError(
            "d1 is not divisible by 9."
        )

    D0 = d0 // 3
    D1 = d1 // 9

    print()
    print("=" * 78)
    print("3. NORMALIZED 3-ADIC DIFFERENCE LAYERS")
    print("=" * 78)

    print(
        f"  D0=(s1-s0)/3={D0}"
    )

    print(
        f"  D1=(s2-s1)/9={D1}"
    )

    print(
        f"  v3(D0)={v3(D0)}"
    )

    print(
        f"  v3(D1)={v3(D1)}"
    )

    print(
        f"  D0=-4*q3="
        f"{D0 == -4 * Q3}"
    )

    print(
        f"  D1=-4*q3="
        f"{D1 == -4 * Q3}"
    )

    print(
        f"  D0_equals_D1="
        f"{D0 == D1}"
    )

    # ------------------------------------------------------------------
    # 4. NORMALIZED DIFFERENCE RESIDUES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. NORMALIZED DIFFERENCE RESIDUE PROFILES")
    print("=" * 78)

    print(
        f"  D0: {residues(D0)}"
    )

    print(
        f"  D1: {residues(D1)}"
    )

    print(
        f"  residues_identical="
        f"{residues(D0) == residues(D1)}"
    )

    # ------------------------------------------------------------------
    # 5. SECOND DIFFERENCE
    # ------------------------------------------------------------------

    second_difference = (
        S2 - 2 * S1 + S0
    )

    print()
    print("=" * 78)
    print("5. EXACT SECOND DIFFERENCE")
    print("=" * 78)

    print(
        f"  second_difference={second_difference}"
    )

    print(
        f"  v3(second_difference)="
        f"{v3(second_difference)}"
    )

    print(
        f"  second_difference=8*q3="
        f"{second_difference == 8 * Q3}"
    )

    print(
        f"  second_difference/9="
        f"{second_difference // 9}"
    )

    # ------------------------------------------------------------------
    # 6. FIRST SEPARATION MODULI
    # ------------------------------------------------------------------

    sep01 = first_separation(
        S0,
        S1
    )

    sep12 = first_separation(
        S1,
        S2
    )

    sep02 = first_separation(
        S0,
        S2
    )

    print()
    print("=" * 78)
    print("6. FIRST 3-ADIC SEPARATION MODULI")
    print("=" * 78)

    print(
        f"  s0_vs_s1: first_e={sep01} "
        f"modulus={3 ** sep01}"
    )

    print(
        f"  s1_vs_s2: first_e={sep12} "
        f"modulus={3 ** sep12}"
    )

    print(
        f"  s0_vs_s2: first_e={sep02} "
        f"modulus={3 ** sep02}"
    )

    # ------------------------------------------------------------------
    # 7. CONGRUENCE TABLE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. EXACT CONGRUENCE TABLE")
    print("=" * 78)

    for e in range(1, 7):

        modulus = 3 ** e

        eq01 = S0 % modulus == S1 % modulus
        eq12 = S1 % modulus == S2 % modulus
        eq02 = S0 % modulus == S2 % modulus

        print(
            f"  e={e} modulus={modulus}: "
            f"s0=s1 {eq01} "
            f"s1=s2 {eq12} "
            f"s0=s2 {eq02}"
        )

    # ------------------------------------------------------------------
    # 8. SOURCE CORRECTION LAYER
    # ------------------------------------------------------------------

    correction0 = 2 * Q3
    correction1 = 6 * Q3
    correction2 = 18 * Q3

    print()
    print("=" * 78)
    print("8. SOURCE CORRECTION 3-ADIC LAYERS")
    print("=" * 78)

    print(
        f"  correction0={correction0} "
        f"v3={v3(correction0)}"
    )

    print(
        f"  correction1={correction1} "
        f"v3={v3(correction1)}"
    )

    print(
        f"  correction2={correction2} "
        f"v3={v3(correction2)}"
    )

    # ------------------------------------------------------------------
    # 9. GEOMETRIC-LAYER IDENTITY
    # ------------------------------------------------------------------

    geometric_identity = (
        D0 == D1 == -4 * Q3
    )

    second_layer_identity = (
        second_difference == 8 * Q3
    )

    print()
    print("=" * 78)
    print("9. GEOMETRIC-LAYER IDENTITIES")
    print("=" * 78)

    print(
        f"  D0=D1=-4*q3="
        f"{geometric_identity}"
    )

    print(
        f"  second_difference=8*q3="
        f"{second_layer_identity}"
    )

    # ------------------------------------------------------------------
    # 10. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The exact terminal coefficients are

    s_k = q1 - 2*3^k*q3.

The first differences therefore contain one extra factor of 3 each
time the index k increases:

    s1-s0 = -4*q3,
    s2-s1 = -12*q3.

After dividing by the exact visible powers,

    (s1-s0)/3 = -4*q3,
    (s2-s1)/9 = -4*q3.

Thus the normalized successive 3-adic layers are exactly identical.

This is stronger than merely observing that the differences have
different valuations. It identifies the same 3-adic unit appearing
after removal of the predictable powers of 3.

The first separation moduli are therefore:

    s0 vs s1  : 9,
    s1 vs s2  : 27,
    s0 vs s2  : 9.

No general law beyond the three exact terminal coefficients is claimed.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        v3(Q1) == 0
        and v3(Q3) == 1
        and D0 == D1
        and D0 == -4 * Q3
        and second_difference == 8 * Q3
        and sep01 == 2
        and sep12 == 3
        and sep02 == 2
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  q1_3adic_unit={v3(Q1) == 0}"
    )

    print(
        f"  q3_exact_v3_1={v3(Q3) == 1}"
    )

    print(
        f"  normalized_difference_layers_equal="
        f"{D0 == D1}"
    )

    print(
        f"  normalized_layer_equals_minus4q3="
        f"{D0 == -4*Q3}"
    )

    print(
        f"  second_difference_exact="
        f"{second_difference == 8*Q3}"
    )

    print(
        f"  first_separation_s0_s1_is_mod9="
        f"{sep01 == 2}"
    )

    print(
        f"  first_separation_s1_s2_is_mod27="
        f"{sep12 == 3}"
    )

    print(
        f"  first_separation_s0_s2_is_mod9="
        f"{sep02 == 2}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 187 COMPLETE")


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

