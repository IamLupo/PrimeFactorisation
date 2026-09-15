#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 210R — EXACT PRINCIPAL-EXPONENT / BASE-7 DIGIT AUDIT
==============================================================================

Corrected version of Experiment 210.

The previous script failed because a nested generator expression was placed
directly inside an f-string.

This version:

    * starts only from q1/q3;
    * reconstructs k_1 from the mod-7 orbit;
    * recursively computes every Hensel digit;
    * forms k_E;
    * computes m_E=(k_E-1)/6;
    * extracts ordinary integer base-7 digits of m_E;
    * compares them with the Hensel digits;
    * verifies the orbit independently.

No hard-coded exponent sequence.
No huge powers of 3.
No nested f-strings.
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

END_E = 12


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


def ratio_mod(modulus: int) -> int:
    return (
        Q1
        * inverse_mod(Q3, modulus)
    ) % modulus


def orbit_mod(k: int, modulus: int) -> int:
    return (
        BASE
        * pow(
            MULT,
            k,
            modulus,
        )
    ) % modulus


def order_3_mod_7e(e: int) -> int:
    return (
        6
        * (P ** (e - 1))
    )


def base_exponent_mod7() -> int:

    rho = ratio_mod(7)
    matches = []

    for k in range(6):

        if orbit_mod(k, 7) == rho:
            matches.append(k)

    if len(matches) != 1:
        raise ArithmeticError(
            f"Expected unique base exponent, got {matches}."
        )

    return matches[0]


def defect_digit(e: int, k: int) -> int:

    modulus_e = P ** e
    modulus_next = P ** (e + 1)

    rho_next = ratio_mod(
        modulus_next
    )

    orbit_next = orbit_mod(
        k,
        modulus_next
    )

    diff = (
        orbit_next
        - rho_next
    ) % modulus_next

    if diff % modulus_e != 0:
        raise ArithmeticError(
            f"At e={e}, k={k} does not solve modulo 7^{e}."
        )

    return (
        diff // modulus_e
    ) % P


def base7_digits(n: int, count: int | None = None):
    """
    Base-7 digits, low-to-high.
    """
    if n < 0:
        raise ValueError(
            "base7_digits requires n >= 0."
        )

    if n == 0:
        digits = [0]
    else:
        digits = []

        while n:
            digits.append(
                n % P
            )
            n //= P

    if count is not None:
        while len(digits) < count:
            digits.append(0)

    return digits


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 210R — EXACT PRINCIPAL-EXPONENT / "
        "BASE-7 DIGIT AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. SOURCE DATA")
    print("=" * 78)

    q1_v7 = valuation_7(Q1)
    q3_v7 = valuation_7(Q3)

    print(
        f"  q1={Q1}"
    )

    print(
        f"  q3={Q3}"
    )

    print(
        f"  q3_mod7={Q3 % P}"
    )

    print(
        f"  v7(q1)={q1_v7}"
    )

    print(
        f"  v7(q3)={q3_v7}"
    )

    # ------------------------------------------------------------------
    # 2. BASE CLASS
    # ------------------------------------------------------------------

    k1 = base_exponent_mod7()

    print()
    print("=" * 78)
    print("2. BASE EXPONENT CLASS")
    print("=" * 78)

    rho7 = ratio_mod(7)
    orbit7 = orbit_mod(k1, 7)

    print(
        f"  rho_mod7={rho7}"
    )

    print(
        f"  k1={k1}"
    )

    print(
        f"  k1_mod6={k1 % 6}"
    )

    print(
        f"  orbit_mod7={orbit7}"
    )

    print(
        f"  base_alignment={orbit7 == rho7}"
    )

    # ------------------------------------------------------------------
    # 3. RECURSIVE HENSEL LIFT
    # ------------------------------------------------------------------

    current_k = k1

    exponent_sequence = [
        k1
    ]

    lift_digits = []
    transition_rows = []

    print()
    print("=" * 78)
    print("3. RECURSIVE HENSEL RECONSTRUCTION")
    print("=" * 78)

    for e in range(
        1,
        END_E,
    ):

        order = order_3_mod_7e(
            e
        )

        delta = defect_digit(
            e,
            current_k
        )

        t = (
            -delta
        ) % P

        next_k = (
            current_k
            + t * order
        )

        modulus_next = P ** (
            e + 1
        )

        orbit_value = orbit_mod(
            next_k,
            modulus_next
        )

        rho_value = ratio_mod(
            modulus_next
        )

        alignment = (
            orbit_value
            == rho_value
        )

        if not alignment:
            raise ArithmeticError(
                f"Lift failed at e={e}."
            )

        lift_digits.append(t)

        transition_rows.append(
            {
                "e": e,
                "current_k": current_k,
                "delta": delta,
                "t": t,
                "order": order,
                "next_k": next_k,
                "alignment": alignment,
            }
        )

        exponent_sequence.append(
            next_k
        )

        print(
            f"  e={e}: "
            f"k_e={current_k} "
            f"delta={delta} "
            f"t={t} "
            f"order={order} "
            f"k_(e+1)={next_k} "
            f"alignment={alignment}"
        )

        current_k = next_k

    # ------------------------------------------------------------------
    # 4. PRINCIPAL EXPONENT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. PRINCIPAL EXPONENT COORDINATE")
    print("=" * 78)

    m_sequence = []

    for E, k in enumerate(
        exponent_sequence,
        start=1
    ):

        remainder = (
            k - 1
        ) % 6

        if remainder != 0:
            raise ArithmeticError(
                f"k={k} is not 1 mod 6."
            )

        m = (
            k - 1
        ) // 6

        m_sequence.append(
            m
        )

        print(
            f"  E={E}: "
            f"k={k} "
            f"m=(k-1)/6={m}"
        )

    # ------------------------------------------------------------------
    # 5. DIRECT INTEGER BASE-7 DIGITS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. DIRECT BASE-7 DIGITS OF m")
    print("=" * 78)

    direct_digit_sequences = []

    for E, m in enumerate(
        m_sequence,
        start=1
    ):

        required_count = E - 1

        digits = base7_digits(
            m,
            required_count,
        )

        visible = digits[
            :required_count
        ]

        direct_digit_sequences.append(
            visible
        )

        print(
            f"  E={E}: "
            f"m={m} "
            f"digits_low_to_high={visible}"
        )

    # ------------------------------------------------------------------
    # 6. HENSEL VS INTEGER BASE-7 DIGITS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. HENSEL DIGIT / INTEGER BASE-7 DIGIT MATCH")
    print("=" * 78)

    digit_match_flags = []

    for E in range(
        1,
        END_E + 1,
    ):

        hensel_prefix = (
            lift_digits[: E - 1]
        )

        integer_digits = (
            direct_digit_sequences[
                E - 1
            ]
        )

        exact = (
            hensel_prefix
            == integer_digits
        )

        digit_match_flags.append(
            exact
        )

        print(
            f"  E={E}: "
            f"Hensel={hensel_prefix} "
            f"integer_base7={integer_digits} "
            f"exact={exact}"
        )

    # ------------------------------------------------------------------
    # 7. CLOSED-FORM RECONSTRUCTION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. CLOSED-FORM PRINCIPAL-EXPONENT RECONSTRUCTION")
    print("=" * 78)

    closed_form_flags = []

    for E in range(
        1,
        END_E + 1,
    ):

        digits = (
            lift_digits[: E - 1]
        )

        m_closed = 0

        for j, t in enumerate(
            digits
        ):

            m_closed += (
                t
                * (P ** j)
            )

        k_closed = (
            1
            + 6 * m_closed
        )

        k_actual = exponent_sequence[
            E - 1
        ]

        exact = (
            k_closed
            == k_actual
        )

        closed_form_flags.append(
            exact
        )

        print(
            f"  E={E}: "
            f"m_closed={m_closed} "
            f"k_closed={k_closed} "
            f"k_actual={k_actual} "
            f"exact={exact}"
        )

    # ------------------------------------------------------------------
    # 8. ZERO-DIGIT CONTROL
    # ------------------------------------------------------------------

    zero_digit_events = []

    for row in transition_rows:

        if row["t"] == 0:
            zero_digit_events.append(
                {
                    "e": row["e"],
                    "k": row["current_k"],
                }
            )

    repeated_exponent_count = 0

    for i in range(
        len(exponent_sequence) - 1
    ):

        if (
            exponent_sequence[i]
            == exponent_sequence[i + 1]
        ):
            repeated_exponent_count += 1

    print()
    print("=" * 78)
    print("8. ZERO-DIGIT CONTROL")
    print("=" * 78)

    print(
        f"  zero_digit_events={zero_digit_events}"
    )

    print(
        f"  repeated_exponent_count="
        f"{repeated_exponent_count}"
    )

    print(
        f"  zero_digits_equal_repeats="
        f"{len(zero_digit_events) == repeated_exponent_count}"
    )

    # ------------------------------------------------------------------
    # 9. GLOBAL ORBIT VALIDATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. GLOBAL ORBIT VALIDATION")
    print("=" * 78)

    orbit_checks = []

    for E, k in enumerate(
        exponent_sequence,
        start=1
    ):

        modulus = P ** E

        orbit_value = orbit_mod(
            k,
            modulus
        )

        rho_value = ratio_mod(
            modulus
        )

        exact = (
            orbit_value
            == rho_value
        )

        orbit_checks.append(
            exact
        )

        print(
            f"  E={E}: "
            f"k={k} "
            f"orbit={orbit_value} "
            f"rho={rho_value} "
            f"exact={exact}"
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
The mod-7 orbit fixes

    k_1 = 1 mod 6.

Therefore write

    k_E = 1 + 6 m_E.

The Hensel recurrence gives

    k_(e+1)
      = k_e + 6*t_e*7^(e-1),

so necessarily

    m_E
      = t_1 + 7*t_2 + 7^2*t_3 + ...

Experiment 210 compares this formal Hensel expansion with the ordinary
integer base-7 representation of the actual integer m_E.

Equality means the Hensel digits are literally the canonical base-7
digits of the principal exponent coordinate.

A zero digit has an exact consequence:

    t_e = 0

if and only if

    k_(e+1) = k_e.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    all_digit_match = all(
        digit_match_flags
    )

    all_closed = all(
        closed_form_flags
    )

    all_orbit = all(
        orbit_checks
    )

    zero_repeat_exact = (
        len(zero_digit_events)
        == repeated_exponent_count
    )

    final_ok = (
        q1_v7 == 0
        and q3_v7 == 0
        and all_digit_match
        and all_closed
        and all_orbit
        and zero_repeat_exact
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  Hensel_digits={lift_digits}"
    )

    print(
        f"  principal_m_sequence={m_sequence}"
    )

    print(
        f"  all_integer_base7_digit_matches="
        f"{all_digit_match}"
    )

    print(
        f"  all_closed_form_reconstructions="
        f"{all_closed}"
    )

    print(
        f"  all_orbit_checks_exact="
        f"{all_orbit}"
    )

    print(
        f"  zero_digit_repeat_identity="
        f"{zero_repeat_exact}"
    )

    print(
        f"  zero_digit_events={zero_digit_events}"
    )

    print(
        f"  repeated_exponent_count="
        f"{repeated_exponent_count}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 210R COMPLETE")


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