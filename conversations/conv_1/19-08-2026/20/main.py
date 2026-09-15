#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 189 — EXACT 7-ADIC TERMINAL SCHUR LAYER /
                UNIT-CANCELLATION AUDIT
==============================================================================

Experiment 188 established an exact 3-adic layer-unit recovery:

    u3 = q3/3,

and

    (s1-s0)/3 = (s2-s1)/9 = -4*u3.

Experiment 189 changes prime.

For the primitive terminal sources:

    q1 = 29144191
    q3 = 24794967

both are 7-adic units:

    v7(q1)=0
    v7(q3)=0.

The terminal Schur coefficients remain

    s0 = q1 - 2 q3
    s1 = q1 - 6 q3
    s2 = q1 - 18 q3.

Unlike the 3-adic case, the correction terms

    2*3^k*q3

are also 7-adic units.

Therefore the Schur coefficients may exhibit genuine 7-adic
cancellation.

The experiment determines exactly:

    * v7(q1), v7(q3);
    * v7 of the three Schur coefficients;
    * v7 of their pairwise differences;
    * residue profiles modulo powers of 7;
    * first separation precision of each pair;
    * normalized difference units;
    * whether a common 7-adic layer exists;
    * whether cancellation, rather than inherited valuation, controls
      the Schur coefficients.

No extrapolation beyond the three exact coefficients.
No SymPy.
No floating point.
"""

from __future__ import annotations

import sys


# ============================================================================
# PRIMITIVE SOURCE DATA
# ============================================================================

Q1 = 29144191
Q3 = 24794967

S0 = Q1 - 2 * Q3
S1 = Q1 - 6 * Q3
S2 = Q1 - 18 * Q3


# ============================================================================
# HELPERS
# ============================================================================

def v7(x):
    if x == 0:
        return None

    x = abs(x)
    e = 0

    while x % 7 == 0:
        x //= 7
        e += 1

    return e


def residue_profile(x, max_e=8):
    return [
        (
            7 ** e,
            x % (7 ** e),
        )
        for e in range(1, max_e + 1)
    ]


def first_separation(a, b, max_e=12):
    for e in range(1, max_e + 1):

        modulus = 7 ** e

        if a % modulus != b % modulus:
            return e

    return None


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 189 — EXACT 7-ADIC TERMINAL SCHUR LAYER / "
        "UNIT-CANCELLATION AUDIT"
    )
    print("=" * 78)

    S = [S0, S1, S2]

    # ------------------------------------------------------------------
    # 1. SOURCE 7-ADIC PROFILE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. SOURCE 7-ADIC PROFILE")
    print("=" * 78)

    print(
        f"  q1={Q1} "
        f"v7={v7(Q1)} "
        f"mod7={Q1 % 7}"
    )

    print(
        f"  q3={Q3} "
        f"v7={v7(Q3)} "
        f"mod7={Q3 % 7}"
    )

    # ------------------------------------------------------------------
    # 2. SCHUR COEFFICIENT 7-ADIC PROFILE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. SCHUR COEFFICIENT 7-ADIC PROFILE")
    print("=" * 78)

    for i, value in enumerate(S):

        print(
            f"  s{i}={value} "
            f"v7={v7(value)} "
            f"mod7={value % 7}"
        )

    # ------------------------------------------------------------------
    # 3. CORRECTION TERMS
    # ------------------------------------------------------------------

    corrections = [
        2 * Q3,
        6 * Q3,
        18 * Q3,
    ]

    print()
    print("=" * 78)
    print("3. EXACT 7-ADIC CORRECTION TERMS")
    print("=" * 78)

    for k, value in enumerate(corrections):

        print(
            f"  correction_k{k}={value} "
            f"v7={v7(value)} "
            f"mod7={value % 7}"
        )

    # ------------------------------------------------------------------
    # 4. PAIRWISE DIFFERENCES
    # ------------------------------------------------------------------

    d01 = S1 - S0
    d12 = S2 - S1
    d02 = S2 - S0

    print()
    print("=" * 78)
    print("4. PAIRWISE DIFFERENCES")
    print("=" * 78)

    print(
        f"  d01=s1-s0={d01} "
        f"v7={v7(d01)}"
    )

    print(
        f"  d12=s2-s1={d12} "
        f"v7={v7(d12)}"
    )

    print(
        f"  d02=s2-s0={d02} "
        f"v7={v7(d02)}"
    )

    print(
        f"  d01=-4*q3={d01 == -4 * Q3}"
    )

    print(
        f"  d12=-12*q3={d12 == -12 * Q3}"
    )

    print(
        f"  d02=-16*q3={d02 == -16 * Q3}"
    )

    # ------------------------------------------------------------------
    # 5. FIRST 7-ADIC SEPARATION
    # ------------------------------------------------------------------

    sep01 = first_separation(
        S0,
        S1,
    )

    sep12 = first_separation(
        S1,
        S2,
    )

    sep02 = first_separation(
        S0,
        S2,
    )

    print()
    print("=" * 78)
    print("5. FIRST 7-ADIC SEPARATION")
    print("=" * 78)

    print(
        f"  s0_vs_s1: first_e={sep01} "
        f"modulus={7 ** sep01}"
    )

    print(
        f"  s1_vs_s2: first_e={sep12} "
        f"modulus={7 ** sep12}"
    )

    print(
        f"  s0_vs_s2: first_e={sep02} "
        f"modulus={7 ** sep02}"
    )

    # ------------------------------------------------------------------
    # 6. CONGRUENCE TABLE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. 7-ADIC CONGRUENCE TABLE")
    print("=" * 78)

    for e in range(1, 7):

        modulus = 7 ** e

        eq01 = (
            S0 % modulus
            == S1 % modulus
        )

        eq12 = (
            S1 % modulus
            == S2 % modulus
        )

        eq02 = (
            S0 % modulus
            == S2 % modulus
        )

        print(
            f"  e={e} modulus={modulus}: "
            f"s0=s1 {eq01} "
            f"s1=s2 {eq12} "
            f"s0=s2 {eq02}"
        )

    # ------------------------------------------------------------------
    # 7. RESIDUE PROFILES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. 7-ADIC RESIDUE PROFILES")
    print("=" * 78)

    for i, value in enumerate(S):

        print(
            f"  s{i}: "
            f"{residue_profile(value)}"
        )

    # ------------------------------------------------------------------
    # 8. MOD-7 CANCELLATION TEST
    # ------------------------------------------------------------------

    source_residues = [
        Q1 % 7,
        Q3 % 7,
    ]

    schur_residues = [
        S0 % 7,
        S1 % 7,
        S2 % 7,
    ]

    expected_residues = [
        (
            source_residues[0]
            - 2 * source_residues[1]
        ) % 7,

        (
            source_residues[0]
            - 6 * source_residues[1]
        ) % 7,

        (
            source_residues[0]
            - 18 * source_residues[1]
        ) % 7,
    ]

    print()
    print("=" * 78)
    print("8. MOD-7 CANCELLATION TEST")
    print("=" * 78)

    print(
        f"  source_residues={source_residues}"
    )

    print(
        f"  expected_schur_residues="
        f"{expected_residues}"
    )

    print(
        f"  actual_schur_residues="
        f"{schur_residues}"
    )

    print(
        f"  mod7_reconstruction_exact="
        f"{expected_residues == schur_residues}"
    )

    print(
        f"  any_schur_zero_mod7="
        f"{0 in schur_residues}"
    )

    # ------------------------------------------------------------------
    # 9. 7-ADIC DIFFERENCE UNITS
    # ------------------------------------------------------------------

    difference_data = [
        ("d01", d01),
        ("d12", d12),
        ("d02", d02),
    ]

    print()
    print("=" * 78)
    print("9. 7-ADIC DIFFERENCE UNITS")
    print("=" * 78)

    for label, value in difference_data:

        e = v7(value)

        unit = (
            value // (7 ** e)
            if e is not None
            else 0
        )

        print(
            f"  {label}: "
            f"v7={e} "
            f"unit={unit} "
            f"unit_mod7={unit % 7 if unit else None}"
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
Unlike the 3-adic case, both primitive source values are 7-adic units:

    v7(q1)=0,
    v7(q3)=0.

Therefore the correction terms and q1 live at the same 7-adic order.

Any positive v7 of a Schur coefficient must consequently arise from
cancellation in

    q1 - 2*3^k*q3.

The pairwise differences remain controlled exactly by q3:

    s1-s0 = -4*q3,
    s2-s1 = -12*q3,
    s2-s0 = -16*q3.

The experiment therefore separates two mechanisms:

    source inheritance:
        difference valuations;

    cancellation:
        valuations of the coefficients themselves.

If a Schur coefficient becomes divisible by 7, that divisibility is
not inherited from either source input individually. It is an exact
cancellation phenomenon in the terminal Schur combination.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        v7(Q1) == 0
        and v7(Q3) == 0
        and d01 == -4 * Q3
        and d12 == -12 * Q3
        and d02 == -16 * Q3
        and expected_residues == schur_residues
        and sep01 is not None
        and sep12 is not None
        and sep02 is not None
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  q1_7adic_unit={v7(Q1) == 0}"
    )

    print(
        f"  q3_7adic_unit={v7(Q3) == 0}"
    )

    print(
        f"  exact_difference_formulas="
        f"{d01 == -4*Q3 and d12 == -12*Q3 and d02 == -16*Q3}"
    )

    print(
        f"  exact_mod7_schur_reconstruction="
        f"{expected_residues == schur_residues}"
    )

    print(
        f"  separation_data_available="
        f"{sep01 is not None and sep12 is not None and sep02 is not None}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 189 COMPLETE")


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

