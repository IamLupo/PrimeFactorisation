#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 200 — EXACT 7-ADIC OVERLIFT SOURCE / QUADRATIC RESIDUAL AUDIT
==============================================================================

Experiment 199 established the exact successful-lift residual valuations:

    e=1 -> k=1:    v7(F)=2
    e=2 -> k=85:   v7(F)=3
    e=3 -> k=1555: v7(F)=5

where

    F(k) = 2*3^k*q3 - q1.

Thus only the final tested lift over-lifts:

    expected precision: 7^4
    actual precision:   7^5.

Experiment 200 isolates the extra digit.

For the successful exponent k at level e, write

    F(k) = 7^(e+1) * U.

The first next digit is

    U mod 7.

When this vanishes, the lift automatically gains one extra 7-adic
digit.

The experiment therefore records:

    * the normalized residual U;
    * U mod 7;
    * U / 7 when U is divisible by 7;
    * the next residual digit mod 7;
    * the exact amount of over-lifting.

It then compares the over-lifted level with the next ordinary Hensel
candidate search.

The purpose is to determine whether the extra digit is simply an
accidental zero in the normalized residual or whether it predicts the
next exponent correction.

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

EXPONENTS = [
    1,
    1,
    85,
    1555,
]


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


def ratio_mod(
    a: int,
    b: int,
    modulus: int,
) -> int:

    if modulus == 1:
        return 0

    return (
        a * inverse_mod(
            b,
            modulus,
        )
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
        f"Order search failed for {a} mod {modulus}."
    )


def orbit_value_mod(
    k: int,
    modulus: int,
) -> int:

    return (
        BASE
        * pow(
            MULT,
            k,
            modulus,
        )
    ) % modulus


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 200 — EXACT 7-ADIC OVERLIFT SOURCE / "
        "QUADRATIC RESIDUAL AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT SOURCE DATA")
    print("=" * 78)

    print(
        f"  q1={Q1}"
    )

    print(
        f"  q3={Q3}"
    )

    print(
        f"  v7(q1)={v7(Q1)}"
    )

    print(
        f"  v7(q3)={v7(Q3)}"
    )

    # ------------------------------------------------------------------
    # 2. SUCCESSFUL EXPONENT LEVELS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. SUCCESSFUL EXPONENT LEVELS")
    print("=" * 78)

    for i, k in enumerate(
        EXPONENTS,
        start=1,
    ):

        modulus = P ** i
        rho = ratio_mod(
            Q1,
            Q3,
            modulus,
        )

        c = orbit_value_mod(
            k,
            modulus,
        )

        F = (
            orbit_value(k) * Q3
            - Q1
        )

        print(
            f"  e={i}: "
            f"modulus={modulus} "
            f"k={k} "
            f"rho={rho} "
            f"orbit={c} "
            f"exact={rho == c} "
            f"v7(F)={v7(F)}"
        )

    # ------------------------------------------------------------------
    # 3. NORMALIZED RESIDUAL DIGITS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. NORMALIZED RESIDUAL DIGITS")
    print("=" * 78)

    records = []

    for e in range(1, 4):

        k = EXPONENTS[e]

        F = (
            orbit_value(k) * Q3
            - Q1
        )

        required_power = e + 1

        if F % (P ** required_power) != 0:
            raise ArithmeticError(
                f"Residual at e={e} is not divisible by "
                f"7^{required_power}."
            )

        U = (
            F // (P ** required_power)
        )

        digit0 = (
            U % P
        )

        if U % P == 0:
            U1 = U // P
            digit1 = (
                U1 % P
            )
        else:
            U1 = None
            digit1 = None

        valuation = v7(F)

        extra_depth = (
            valuation
            - required_power
        )

        record = {
            "e": e,
            "k": k,
            "F": F,
            "required_power": required_power,
            "U": U,
            "digit0": digit0,
            "U1": U1,
            "digit1": digit1,
            "v7": valuation,
            "extra_depth": extra_depth,
        }

        records.append(
            record
        )

        print(
            f"  e={e}: "
            f"k={k} "
            f"v7(F)={valuation} "
            f"F/7^{required_power} mod7={digit0} "
            f"extra_depth={extra_depth}"
        )

    # ------------------------------------------------------------------
    # 4. DEEPER LAYER WHEN OVERLIFT OCCURS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. DEEPER OVERLIFT LAYER")
    print("=" * 78)

    for record in records:

        if record["digit0"] == 0:

            print(
                f"  e={record['e']}: "
                f"U={record['U']}"
            )

            print(
                f"    U/7={record['U1']}"
            )

            print(
                f"    next_digit_after_overlift="
                f"{record['digit1']}"
            )

        else:

            print(
                f"  e={record['e']}: "
                f"no overlift; "
                f"first_next_digit={record['digit0']}"
            )

    # ------------------------------------------------------------------
    # 5. EXACT OVERLIFT LOCATION
    # ------------------------------------------------------------------

    overlift_levels = [
        record["e"]
        for record in records
        if record["extra_depth"] > 0
    ]

    print()
    print("=" * 78)
    print("5. EXACT OVERLIFT LOCATION")
    print("=" * 78)

    print(
        f"  overlift_levels={overlift_levels}"
    )

    print(
        f"  overlift_count="
        f"{len(overlift_levels)}"
    )

    # ------------------------------------------------------------------
    # 6. SECOND RESIDUAL DIGIT AT e=3
    # ------------------------------------------------------------------

    e3_record = next(
        record
        for record in records
        if record["e"] == 3
    )

    e3_second_digit = (
        e3_record["digit1"]
    )

    print()
    print("=" * 78)
    print("6. e=3 OVERLIFT DIGITS")
    print("=" * 78)

    print(
        f"  k=1555"
    )

    print(
        f"  F/7^4={e3_record['U']}"
    )

    print(
        f"  first_digit={e3_record['digit0']}"
    )

    print(
        f"  second_digit={e3_second_digit}"
    )

    print(
        f"  v7(F)={e3_record['v7']}"
    )

    # ------------------------------------------------------------------
    # 7. NEXT MODULUS VALIDATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. OVERLIFT MODULUS VALIDATION")
    print("=" * 78)

    next_validity = []

    for record in records:

        power = (
            record["required_power"]
            + 1
        )

        exact = (
            record["F"]
            % (P ** power)
            == 0
        )

        next_validity.append(
            exact
        )

        print(
            f"  e={record['e']}: "
            f"extension_to_7^{power}="
            f"{exact}"
        )

    # ------------------------------------------------------------------
    # 8. COMPARE WITH NEXT HENSEL DIGIT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. OVERLIFT VS NEXT HENSEL DIGIT")
    print("=" * 78)

    # At e=3, the successful k_4 already satisfies the modulus 7^5.
    #
    # To move beyond it, the next exponent search would normally use
    #
    #     k_4 + t * ord_{7^4}(3).
    #
    # The current exponent therefore has t=0 as a valid extension
    # through 7^5.

    order_e4 = multiplicative_order(
        MULT,
        7 ** 4,
    )

    k4 = EXPONENTS[3]

    print(
        f"  k4={k4}"
    )

    print(
        f"  ord_2401(3)={order_e4}"
    )

    print(
        f"  zero_digit_extension="
        f"{e3_second_digit == 0}"
    )

    print(
        f"  candidate_t0_next_exponent="
        f"{k4}"
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
Experiment 199 showed that only the e=3 -> e=4 lift over-lifts.

Experiment 200 identifies the exact reason.

At that level,

    F(k4) / 7^4

is itself divisible by 7.

Therefore the already-selected exponent k4 solves not merely

    F(k4) == 0 (mod 7^4),

but also

    F(k4) == 0 (mod 7^5).

The first residual digit is zero.

The next residual digit determines whether the same exponent extends
further.

This is an exact separation between:

    Hensel lift choice,
    accidental over-lift,
    subsequent residual obstruction.
"""
    )

    # ------------------------------------------------------------------
    # 10. FINAL
    # ------------------------------------------------------------------

    actual_v7 = [
        record["v7"]
        for record in records
    ]

    actual_digits = [
        record["digit0"]
        for record in records
    ]

    actual_extra = [
        record["extra_depth"]
        for record in records
    ]

    final_ok = (
        v7(Q1) == 0
        and v7(Q3) == 0
        and actual_v7 == [2, 3, 5]
        and actual_digits == [5, 2, 0]
        and actual_extra == [0, 0, 1]
        and overlift_levels == [3]
        and e3_second_digit is not None
        and next_validity == [False, False, True]
    )

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  actual_v7={actual_v7}"
    )

    print(
        f"  actual_first_digits={actual_digits}"
    )

    print(
        f"  actual_extra_depth={actual_extra}"
    )

    print(
        f"  overlift_levels={overlift_levels}"
    )

    print(
        f"  only_e3_overlifts="
        f"{overlift_levels == [3]}"
    )

    print(
        f"  extension_to_7^5_at_e3="
        f"{next_validity[2]}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 200 COMPLETE")


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

