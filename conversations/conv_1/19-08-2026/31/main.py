#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 199 — EXACT 7-ADIC OVERLIFT / SECOND RESIDUAL-DIGIT AUDIT
==============================================================================

Experiment 198 established the affine residual law

    R_e(t) = A_e + t  (mod 7)

for the successive exponent lifts, with

    A = [0, 5, 2]
    t = [0, 2, 5].

The corresponding exponents are

    k = [1, 1, 85, 1555].

Experiment 198 also revealed a higher-order event:

    e=3 -> e=4:
        required v7(F) >= 4,
        actual    v7(F) = 5.

Experiment 199 isolates this extra digit.

For each successful lift k_{e+1}, define

    F = 2*3^(k_{e+1})*q3 - q1.

After the required division by

    7^(e+1),

inspect the next residue modulo 7.

If that residue is zero, the lift automatically extends one more
7-adic level.

Thus the experiment distinguishes:

    ordinary one-digit lifting
        from
    accidental higher-order lifting.

No floating point.
No SymPy.
"""

from __future__ import annotations

import sys


# ============================================================================
# DATA
# ============================================================================

Q1 = 29144191
Q3 = 24794967

BASE = 2
MULT = 3
P = 7

# Established exponent sequence.
EXPONENTS = [1, 1, 85, 1555]


# ============================================================================
# HELPERS
# ============================================================================

def v7(x: int):
    if x == 0:
        return None

    x = abs(x)
    e = 0

    while x % 7 == 0:
        x //= 7
        e += 1

    return e


def orbit_value(k: int) -> int:
    return BASE * (MULT ** k)


def inverse_mod(a: int, m: int) -> int:
    if m <= 1:
        return 0

    a %= m

    if a == 0:
        raise ArithmeticError(
            f"Cannot invert 0 modulo {m}."
        )

    old_r, r = a, m
    old_s, s = 1, 0

    while r:
        q = old_r // r

        old_r, r = (
            r,
            old_r - q * r,
        )

        old_s, s = (
            s,
            old_s - q * s,
        )

    if old_r != 1:
        raise ArithmeticError(
            f"{a} is not invertible modulo {m}."
        )

    return old_s % m


def ratio_mod(a: int, b: int, modulus: int) -> int:
    if modulus == 1:
        return 0

    return (
        a * inverse_mod(b, modulus)
    ) % modulus


def multiplicative_order(
    a: int,
    modulus: int,
) -> int:

    if modulus == 1:
        return 1

    a %= modulus

    if a == 0:
        raise ArithmeticError(
            f"{a} is not a unit modulo {modulus}."
        )

    x = a

    for n in range(1, 1_000_000):

        if x == 1:
            return n

        x = (
            x * a
        ) % modulus

    raise ArithmeticError(
        f"Could not determine order of {a} modulo {modulus}."
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 199 — EXACT 7-ADIC OVERLIFT / "
        "SECOND RESIDUAL-DIGIT AUDIT"
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
    print(f"  v7(q1)={v7(Q1)}")
    print(f"  v7(q3)={v7(Q3)}")

    # ------------------------------------------------------------------
    # 2. ESTABLISHED EXPONENT LEVELS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. ESTABLISHED EXPONENT LEVELS")
    print("=" * 78)

    for i, k in enumerate(EXPONENTS, start=1):

        modulus = P ** i
        order = multiplicative_order(
            MULT,
            modulus,
        )

        rho = ratio_mod(
            Q1,
            Q3,
            modulus,
        )

        c_mod = (
            BASE
            * pow(MULT, k, modulus)
        ) % modulus

        print(
            f"  e={i}: "
            f"modulus={modulus} "
            f"order={order} "
            f"k={k} "
            f"rho={rho} "
            f"orbit={c_mod} "
            f"exact={rho == c_mod}"
        )

    # ------------------------------------------------------------------
    # 3. SUCCESSFUL LIFT RESIDUALS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. SUCCESSFUL LIFT RESIDUALS")
    print("=" * 78)

    lift_rows = []

    for e in range(1, 4):

        k_next = EXPONENTS[e]

        c = orbit_value(
            k_next
        )

        F = (
            c * Q3
            - Q1
        )

        valuation = v7(F)

        required_power = e + 1

        if F % (P ** required_power) != 0:
            raise ArithmeticError(
                f"Lift at e={e} is not valid to level "
                f"{required_power}."
            )

        first_quotient = (
            F // (P ** required_power)
        )

        next_digit = (
            first_quotient % P
        )

        next_valuation = (
            valuation
            if valuation is not None
            else 99
        )

        extra_depth = (
            next_valuation
            - required_power
        )

        row = {
            "e": e,
            "k": k_next,
            "F": F,
            "v7": valuation,
            "required_power": required_power,
            "quotient": first_quotient,
            "next_digit": next_digit,
            "extra_depth": extra_depth,
        }

        lift_rows.append(
            row
        )

        print(
            f"  e={e}: "
            f"k={k_next} "
            f"v7(F)={valuation} "
            f"required_power={required_power} "
            f"next_digit={next_digit} "
            f"extra_depth={extra_depth}"
        )

    # ------------------------------------------------------------------
    # 4. OVERLIFT TEST
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. OVERLIFT TEST")
    print("=" * 78)

    overlift_flags = []

    for row in lift_rows:

        is_overlift = (
            row["next_digit"] == 0
        )

        overlift_flags.append(
            is_overlift
        )

        print(
            f"  e={row['e']}: "
            f"next_digit={row['next_digit']} "
            f"overlift={is_overlift}"
        )

    # ------------------------------------------------------------------
    # 5. NEXT-LEVEL CONSEQUENCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. NEXT-LEVEL CONSEQUENCE")
    print("=" * 78)

    for row in lift_rows:

        next_power = (
            row["required_power"] + 1
        )

        automatic_extension = (
            row["F"] % (P ** next_power)
            == 0
        )

        print(
            f"  e={row['e']}: "
            f"extension_to_7^{next_power}="
            f"{automatic_extension}"
        )

    # ------------------------------------------------------------------
    # 6. NORMALIZED NEXT DIGIT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. NORMALIZED NEXT-DIGIT DATA")
    print("=" * 78)

    for row in lift_rows:

        q = row["quotient"]

        print(
            f"  e={row['e']}: "
            f"F/7^{row['required_power']}="
            f"{q} "
            f"mod7={q % 7}"
        )

    # ------------------------------------------------------------------
    # 7. HIGHEST VERIFIED OVERSHOOT
    # ------------------------------------------------------------------

    maximum_extra_depth = max(
        row["extra_depth"]
        for row in lift_rows
    )

    deepest_rows = [
        row["e"]
        for row in lift_rows
        if row["extra_depth"]
        == maximum_extra_depth
    ]

    print()
    print("=" * 78)
    print("7. HIGHEST VERIFIED OVERLIFT")
    print("=" * 78)

    print(
        f"  maximum_extra_depth={maximum_extra_depth}"
    )

    print(
        f"  deepest_e={deepest_rows}"
    )

    # ------------------------------------------------------------------
    # 8. SOURCE-RATIO INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. SOURCE-RATIO INTERPRETATION")
    print("=" * 78)

    rho343 = ratio_mod(
        Q1,
        Q3,
        343,
    )

    rho2401 = ratio_mod(
        Q1,
        Q3,
        2401,
    )

    print(
        f"  rho_mod343={rho343}"
    )

    print(
        f"  rho_mod2401={rho2401}"
    )

    print(
        """
The successful exponent lift determines a point on the orbit
2*3^k whose residual vanishes to the required 7-adic order.

The quotient

    F / 7^(e+1)

contains the next residual digit.

If this digit is nonzero, the lift stops at the expected precision.

If it is zero, the same exponent already solves the next 7-adic
congruence without changing the exponent.

Therefore an overlift is an exact arithmetic coincidence, distinct
from the ordinary Hensel step.

The experiment checks this directly.
"""
    )

    # ------------------------------------------------------------------
    # 9. FINAL
    # ------------------------------------------------------------------

    expected_v7 = [
        2,
        3,
        5,
    ]

    actual_v7 = [
        row["v7"]
        for row in lift_rows
    ]

    expected_digits = [
        0,
        0,
        0,
    ]

    actual_digits = [
        row["next_digit"]
        for row in lift_rows
    ]

    expected_extra_depth = [
        1,
        1,
        1,
    ]

    actual_extra_depth = [
        row["extra_depth"]
        for row in lift_rows
    ]

    final_ok = (
        v7(Q1) == 0
        and v7(Q3) == 0
        and actual_v7 == expected_v7
        and actual_digits == expected_digits
        and actual_extra_depth == expected_extra_depth
        and overlift_flags == [True, True, True]
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  source_q1_7adic_unit={v7(Q1) == 0}"
    )

    print(
        f"  source_q3_7adic_unit={v7(Q3) == 0}"
    )

    print(
        f"  actual_v7={actual_v7}"
    )

    print(
        f"  expected_v7={expected_v7}"
    )

    print(
        f"  actual_next_digits={actual_digits}"
    )

    print(
        f"  all_next_digits_zero="
        f"{actual_digits == expected_digits}"
    )

    print(
        f"  actual_extra_depth={actual_extra_depth}"
    )

    print(
        f"  all_levels_overlifted="
        f"{overlift_flags == [True, True, True]}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 199 COMPLETE")


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:
        print("\nInterrupted.")

    except Exception as exc:

        print(
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )
        raise

