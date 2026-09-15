#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 188 — EXACT 3-ADIC LAYER-UNIT RECOVERY / SCHUR SOURCE AUDIT
==============================================================================

Experiment 187 established:

    s0 = q1 - 2 q3
    s1 = q1 - 6 q3
    s2 = q1 - 18 q3,

with

    v3(q1)=0,
    v3(q3)=1.

The corrected identities are:

    D0 = (s1-s0)/3 = -4*(q3/3)
    D1 = (s2-s1)/9 = -4*(q3/3)

and

    s2 - 2*s1 + s0 = -8*q3.

Define the 3-adic unit

    u3 = q3/3.

This experiment recovers u3 directly from the Schur row and checks
whether the two normalized difference layers give exactly the same
unit.

It then reconstructs:

    q3 = 3*u3
    q1 = s0 + 2*q3

and verifies all three Schur coefficients.

No floating point.
No SymPy.
No extrapolation.
"""

from __future__ import annotations

import sys


# ============================================================================
# EXACT SOURCE
# ============================================================================

Q1 = 29144191
Q3 = 24794967

S0 = Q1 - 2 * Q3
S1 = Q1 - 6 * Q3
S2 = Q1 - 18 * Q3


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


def first_layer_unit(x, power):
    if x % (3 ** power) != 0:
        raise ArithmeticError(
            f"{x} is not divisible by 3^{power}"
        )

    return x // (3 ** power)


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 188 — EXACT 3-ADIC LAYER-UNIT RECOVERY / "
        "SCHUR SOURCE AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT SOURCE DATA")
    print("=" * 78)

    print(f"  q1={Q1}")
    print(f"  q3={Q3}")

    print(f"  v3(q1)={v3(Q1)}")
    print(f"  v3(q3)={v3(Q3)}")

    u3 = Q3 // 3

    print(
        f"  u3=q3/3={u3}"
    )

    print(
        f"  v3(u3)={v3(u3)}"
    )

    # ------------------------------------------------------------------
    # 2. SCHUR ROW
    # ------------------------------------------------------------------

    S = [S0, S1, S2]

    print()
    print("=" * 78)
    print("2. EXACT SCHUR ROW")
    print("=" * 78)

    for i, value in enumerate(S):
        print(
            f"  s{i}={value} "
            f"v3={v3(value)}"
        )

    # ------------------------------------------------------------------
    # 3. NORMALIZED DIFFERENCE LAYERS
    # ------------------------------------------------------------------

    d0 = S1 - S0
    d1 = S2 - S1

    D0 = first_layer_unit(
        d0,
        1,
    )

    D1 = first_layer_unit(
        d1,
        2,
    )

    print()
    print("=" * 78)
    print("3. NORMALIZED 3-ADIC DIFFERENCE LAYERS")
    print("=" * 78)

    print(
        f"  d0=s1-s0={d0}"
    )

    print(
        f"  d1=s2-s1={d1}"
    )

    print(
        f"  D0=(s1-s0)/3={D0}"
    )

    print(
        f"  D1=(s2-s1)/9={D1}"
    )

    print(
        f"  D0_v3={v3(D0)}"
    )

    print(
        f"  D1_v3={v3(D1)}"
    )

    print(
        f"  D0_equals_D1={D0 == D1}"
    )

    print(
        f"  D0=-4*u3={D0 == -4*u3}"
    )

    print(
        f"  D1=-4*u3={D1 == -4*u3}"
    )

    # ------------------------------------------------------------------
    # 4. RECOVER q3
    # ------------------------------------------------------------------

    recovered_u3_from_D0 = -D0 // 4
    recovered_u3_from_D1 = -D1 // 4

    recovered_q3_from_D0 = (
        3 * recovered_u3_from_D0
    )

    recovered_q3_from_D1 = (
        3 * recovered_u3_from_D1
    )

    print()
    print("=" * 78)
    print("4. EXACT q3 / UNIT RECOVERY")
    print("=" * 78)

    print(
        f"  u3_from_D0={recovered_u3_from_D0}"
    )

    print(
        f"  u3_from_D1={recovered_u3_from_D1}"
    )

    print(
        f"  q3_from_D0={recovered_q3_from_D0}"
    )

    print(
        f"  q3_from_D1={recovered_q3_from_D1}"
    )

    print(
        f"  both_u3_recoveries_exact="
        f"{recovered_u3_from_D0 == u3 and recovered_u3_from_D1 == u3}"
    )

    print(
        f"  both_q3_recoveries_exact="
        f"{recovered_q3_from_D0 == Q3 and recovered_q3_from_D1 == Q3}"
    )

    # ------------------------------------------------------------------
    # 5. RECOVER q1
    # ------------------------------------------------------------------

    recovered_q1_a = S0 + 2 * Q3
    recovered_q1_b = S1 + 6 * Q3
    recovered_q1_c = S2 + 18 * Q3

    print()
    print("=" * 78)
    print("5. EXACT q1 RECOVERY")
    print("=" * 78)

    print(
        f"  q1_from_s0={recovered_q1_a}"
    )

    print(
        f"  q1_from_s1={recovered_q1_b}"
    )

    print(
        f"  q1_from_s2={recovered_q1_c}"
    )

    print(
        f"  all_q1_recoveries_exact="
        f"{recovered_q1_a == Q1 and recovered_q1_b == Q1 and recovered_q1_c == Q1}"
    )

    # ------------------------------------------------------------------
    # 6. SECOND DIFFERENCE
    # ------------------------------------------------------------------

    second_difference = (
        S2 - 2 * S1 + S0
    )

    print()
    print("=" * 78)
    print("6. SECOND DIFFERENCE")
    print("=" * 78)

    print(
        f"  second_difference={second_difference}"
    )

    print(
        f"  v3(second_difference)="
        f"{v3(second_difference)}"
    )

    print(
        f"  second_difference=-8*q3="
        f"{second_difference == -8 * Q3}"
    )

    # ------------------------------------------------------------------
    # 7. DIRECT SOURCE-LAYER FORM
    # ------------------------------------------------------------------

    reconstructed = [
        Q1 - 2 * (3 ** 0) * 3 * u3,
        Q1 - 2 * (3 ** 1) * 3 * u3,
        Q1 - 2 * (3 ** 2) * 3 * u3,
    ]

    print()
    print("=" * 78)
    print("7. DIRECT 3-ADIC SOURCE-LAYER FORM")
    print("=" * 78)

    print(
        f"  reconstructed={reconstructed}"
    )

    print(
        f"  original={S}"
    )

    print(
        f"  exact="
        f"{reconstructed == S}"
    )

    # ------------------------------------------------------------------
    # 8. 3-ADIC UNIT RESIDUE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. 3-ADIC UNIT PROFILE")
    print("=" * 78)

    print(
        f"  u3_mod3={u3 % 3}"
    )

    print(
        f"  -4*u3_mod3={(-4*u3) % 3}"
    )

    print(
        f"  D0_mod3={D0 % 3}"
    )

    print(
        f"  D1_mod3={D1 % 3}"
    )

    # ------------------------------------------------------------------
    # 9. INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 187 showed that the two normalized difference layers are
equal, but the previous expected identity was off by the factor 3 coming
from v3(q3)=1.

The correct layer unit is

    u3 = q3 / 3.

The exact identities are

    (s1-s0)/3 = -4*u3,

    (s2-s1)/9 = -4*u3.

Thus the same 3-adic unit is recovered independently from both adjacent
difference layers.

This gives an exact source-recovery statement:

    Schur difference layers
        ->
    u3
        ->
    q3
        ->
    q1.

The experiment remains finite and descriptive; it does not assume that
the same structure occurs for other terminal pairs.
"""
    )

    # ------------------------------------------------------------------
    # 10. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        v3(Q1) == 0
        and v3(Q3) == 1
        and v3(u3) == 0
        and D0 == D1
        and D0 == -4 * u3
        and D1 == -4 * u3
        and recovered_u3_from_D0 == u3
        and recovered_u3_from_D1 == u3
        and recovered_q3_from_D0 == Q3
        and recovered_q3_from_D1 == Q3
        and recovered_q1_a == Q1
        and recovered_q1_b == Q1
        and recovered_q1_c == Q1
        and second_difference == -8 * Q3
        and reconstructed == S
    )

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  q1_v3_0={v3(Q1) == 0}"
    )

    print(
        f"  q3_v3_1={v3(Q3) == 1}"
    )

    print(
        f"  u3_is_3adic_unit={v3(u3) == 0}"
    )

    print(
        f"  normalized_layers_equal={D0 == D1}"
    )

    print(
        f"  normalized_layer_is_minus4u3="
        f"{D0 == -4*u3}"
    )

    print(
        f"  exact_u3_recovery="
        f"{recovered_u3_from_D0 == u3 and recovered_u3_from_D1 == u3}"
    )

    print(
        f"  exact_q3_recovery="
        f"{recovered_q3_from_D0 == Q3 and recovered_q3_from_D1 == Q3}"
    )

    print(
        f"  exact_q1_recovery="
        f"{recovered_q1_a == Q1 and recovered_q1_b == Q1 and recovered_q1_c == Q1}"
    )

    print(
        f"  second_difference_exact="
        f"{second_difference == -8*Q3}"
    )

    print(
        f"  direct_source_layer_reconstruction="
        f"{reconstructed == S}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 188 COMPLETE")


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

