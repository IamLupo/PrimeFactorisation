#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 205 — EXACT 7-ADIC SOURCE-DEFECT DIGIT RECURRENCE AUDIT
==============================================================================

Experiment 204 showed that the previous sign convention was wrong.

With

    rho = q1/q3,

and

    F(k) = 2*3^k*q3 - q1,

we have exactly

    F(k) = q3 * (2*3^k - rho).

Therefore, if

    delta_e =
        (2*3^k - rho) / 7^e mod 7,

then

    A_e = F(k)/7^e mod 7
        = q3 * delta_e mod 7.

For the present data,

    q3 mod 7 = 1,

so

    A_e = delta_e mod 7.

The observed sequence is

    delta = [0,5,4,1,2]

through e=8,...,12.

Experiment 205 does NOT assume a recurrence.

It tests exact low-complexity possibilities for the defect digits:

    delta_{e+1}
        versus
    delta_e,

    delta_{e+1}-delta_e,

    affine maps
        delta_{e+1} = a*delta_e + b mod 7,

and short-period behavior.

This is a finite exact recurrence audit.
"""

from __future__ import annotations

import sys


# ============================================================================
# DATA
# ============================================================================

Q1 = 29144191
Q3 = 24794967

P = 7
BASE = 2
MULT = 3

LEVELS = [
    (8, 3401371),
    (9, 3401371),
    (10, 72578983),
    (11, 798943909),
    (12, 10968052873),
]


# ============================================================================
# HELPERS
# ============================================================================

def valuation_7(x: int):
    if x == 0:
        return None

    x = abs(x)
    e = 0

    while x % P == 0:
        x //= P
        e += 1

    return e


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
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s

    if old_r != 1:
        raise ArithmeticError(
            f"{a} is not invertible modulo {m}."
        )

    return old_s % m


def ratio_mod(modulus: int) -> int:
    return (
        Q1 * inverse_mod(Q3, modulus)
    ) % modulus


def orbit_mod(k: int, modulus: int) -> int:
    return (
        BASE * pow(MULT, k, modulus)
    ) % modulus


def order_3_mod_7e(e: int) -> int:
    return 6 * (P ** (e - 1))


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 205 — EXACT 7-ADIC SOURCE-DEFECT "
        "DIGIT RECURRENCE AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. SOURCE DATA")
    print("=" * 78)

    print(f"  q1={Q1}")
    print(f"  q3={Q3}")
    print(f"  v7(q1)={valuation_7(Q1)}")
    print(f"  v7(q3)={valuation_7(Q3)}")
    print(f"  q3_mod7={Q3 % P}")

    # ------------------------------------------------------------------
    # 2. EXACT DEFECT DIGITS
    # ------------------------------------------------------------------

    rows = []

    print()
    print("=" * 78)
    print("2. EXACT SOURCE-RATIO DEFECT DIGITS")
    print("=" * 78)

    for e, k in LEVELS:

        modulus_e = P ** e
        modulus_next = P ** (e + 1)

        rho_next = ratio_mod(
            modulus_next
        )

        orbit_next = orbit_mod(
            k,
            modulus_next
        )

        defect_residue = (
            orbit_next - rho_next
        ) % modulus_next

        if defect_residue % modulus_e != 0:
            raise ArithmeticError(
                f"Defect at e={e} is not divisible by 7^{e}."
            )

        delta = (
            defect_residue // modulus_e
        ) % P

        # Direct normalized F digit.
        F_residue = (
            orbit_next * Q3 - Q1
        ) % modulus_next

        if F_residue % modulus_e != 0:
            raise ArithmeticError(
                f"F(k) at e={e} is not divisible by 7^{e}."
            )

        A = (
            F_residue // modulus_e
        ) % P

        expected_A = (
            (Q3 % P) * delta
        ) % P

        rows.append(
            {
                "e": e,
                "k": k,
                "delta": delta,
                "A": A,
                "expected_A": expected_A,
            }
        )

        print(
            f"  e={e}: "
            f"k={k} "
            f"delta={delta} "
            f"A={A} "
            f"q3*delta={expected_A} "
            f"exact={A == expected_A}"
        )

    delta_seq = [
        row["delta"]
        for row in rows
    ]

    A_seq = [
        row["A"]
        for row in rows
    ]

    # ------------------------------------------------------------------
    # 3. CORRECTED SOURCE-TO-INTERCEPT IDENTITY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. CORRECT SOURCE-TO-INTERCEPT IDENTITY")
    print("=" * 78)

    identity_ok = all(
        row["A"] == row["expected_A"]
        for row in rows
    )

    print(
        f"  A_e = q3*delta_e mod7 exact={identity_ok}"
    )

    print(
        f"  delta_sequence={delta_seq}"
    )

    print(
        f"  A_sequence={A_seq}"
    )

    # ------------------------------------------------------------------
    # 4. FIRST DIFFERENCE AUDIT
    # ------------------------------------------------------------------

    first_diffs = [
        (
            delta_seq[i + 1]
            - delta_seq[i]
        ) % P
        for i in range(
            len(delta_seq) - 1
        )
    ]

    print()
    print("=" * 78)
    print("4. DEFECT FIRST DIFFERENCES")
    print("=" * 78)

    print(
        f"  first_differences_mod7={first_diffs}"
    )

    # ------------------------------------------------------------------
    # 5. CONSTANT DIFFERENCE TEST
    # ------------------------------------------------------------------

    constant_difference = (
        len(set(first_diffs)) == 1
    )

    print(
        f"  constant_difference={constant_difference}"
    )

    if constant_difference:
        print(
            f"  common_difference={first_diffs[0]}"
        )

    # ------------------------------------------------------------------
    # 6. AFFINE RECURRENCE SEARCH
    # ------------------------------------------------------------------

    affine_solutions = []

    for a in range(P):
        for b in range(P):

            ok = True

            for i in range(
                len(delta_seq) - 1
            ):

                predicted = (
                    a * delta_seq[i]
                    + b
                ) % P

                actual = delta_seq[i + 1]

                if predicted != actual:
                    ok = False
                    break

            if ok:
                affine_solutions.append(
                    (a, b)
                )

    print()
    print("=" * 78)
    print("6. AFFINE RECURRENCE SEARCH")
    print("=" * 78)

    print(
        "  test: delta_(e+1) = a*delta_e + b mod7"
    )

    print(
        f"  affine_solutions={affine_solutions}"
    )

    # ------------------------------------------------------------------
    # 7. SPECIAL RECURRENCE TESTS
    # ------------------------------------------------------------------

    negation = all(
        delta_seq[i + 1]
        == (-delta_seq[i]) % P
        for i in range(
            len(delta_seq) - 1
        )
    )

    increment_one = all(
        delta_seq[i + 1]
        == (delta_seq[i] + 1) % P
        for i in range(
            len(delta_seq) - 1
        )
    )

    multiply_three = all(
        delta_seq[i + 1]
        == (3 * delta_seq[i]) % P
        for i in range(
            len(delta_seq) - 1
        )
    )

    multiply_minus_one = all(
        delta_seq[i + 1]
        == (-delta_seq[i]) % P
        for i in range(
            len(delta_seq) - 1
        )
    )

    print()
    print("=" * 78)
    print("7. SIMPLE MAP TESTS")
    print("=" * 78)

    print(
        f"  delta_next=-delta={multiply_minus_one}"
    )

    print(
        f"  delta_next=delta+1={increment_one}"
    )

    print(
        f"  delta_next=3*delta={multiply_three}"
    )

    print(
        f"  delta_next=negation={negation}"
    )

    # ------------------------------------------------------------------
    # 8. SHORT PERIOD TEST
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. SHORT PERIOD TEST")
    print("=" * 78)

    period_results = {}

    for period in range(
        1,
        len(delta_seq) + 1
    ):

        periodic = True

        for i in range(
            len(delta_seq)
        ):

            if (
                delta_seq[i]
                != delta_seq[
                    i % period
                ]
            ):
                periodic = False
                break

        period_results[period] = periodic

    print(
        f"  periodicity={period_results}"
    )

    # ------------------------------------------------------------------
    # 9. RELATION TO LIFT DIGITS
    # ------------------------------------------------------------------

    lift_digits = [
        (-delta) % P
        for delta in delta_seq
    ]

    print()
    print("=" * 78)
    print("9. DEFECT DIGITS VS LIFT DIGITS")
    print("=" * 78)

    print(
        f"  delta={delta_seq}"
    )

    print(
        f"  lift_digit=-delta mod7={lift_digits}"
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
The previous experiment's negative result was caused by the wrong sign
in the source-defect identity, not by the data.

The correct relation is

    F(k)
      = q3 * (2*3^k - rho),

hence

    A_e = q3*delta_e mod 7.

Because

    q3 == 1 mod 7,

the intercept and source-defect digits are identical:

    A_e = delta_e.

Experiment 205 therefore asks the genuinely new question:

    Is the defect sequence itself generated by a simple recurrence?

The tested possibilities are deliberately small:

    constant difference,
    affine map over F_7,
    negation,
    increment by one,
    multiplication by three,
    short period.

A negative result would mean the stable slope-one Hensel law is real,
but the source-ratio defect digits remain level-dependent data.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        identity_ok
        and A_seq == delta_seq
        and lift_digits == [
            (-x) % P
            for x in delta_seq
        ]
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  corrected_source_intercept_identity={identity_ok}"
    )

    print(
        f"  A_equals_delta={A_seq == delta_seq}"
    )

    print(
        f"  delta_sequence={delta_seq}"
    )

    print(
        f"  lift_digits={lift_digits}"
    )

    print(
        f"  constant_difference={constant_difference}"
    )

    print(
        f"  affine_solutions={affine_solutions}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 205 COMPLETE")


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

