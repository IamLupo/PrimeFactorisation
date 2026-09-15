#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 196 — EXACT 7-ADIC EXPONENT-LIFT / ORBIT HENSEL AUDIT
==============================================================================

Experiment 195 established:

    rho = q1/q3,

and the orbit

    c_k = 2*3^k.

The unique orbit matches are:

    mod 7:
        k = 1

    mod 49:
        k = 1

    mod 343:
        k = 85

The previous Experiment 195 incorrectly expected no orbit point modulo
343. The complete mod-343 orbit actually contains exactly one match,
namely k=85.

Experiment 196 tracks the exponent that solves

    2*3^k == rho (mod 7^e)

through successive powers of 7.

The exact orders are expected to be

    ord_{7^e}(3) = 6*7^(e-1).

For each e we search one complete multiplicative orbit and record:

    * rho mod 7^e;
    * orbit order;
    * unique matching exponent k_e;
    * congruence relation between k_e and k_(e-1);
    * the lift quotient
          (k_e-k_(e-1))/ord_{7^(e-1)}(3);
    * v7(rho-c_k) at the resulting exponent.

This is an exact finite Hensel-style lifting audit.

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

MAX_E = 4


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


def inverse_mod(a, m):
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


def ratio_mod(a, b, modulus):
    if modulus == 1:
        return 0

    return (
        a * inverse_mod(b, modulus)
    ) % modulus


def orbit_value_mod(k, modulus):
    return (
        BASE
        * pow(MULT, k, modulus)
    ) % modulus


def multiplicative_order(a, modulus):
    """
    Exact order by finite search. For the small powers used here this
    is inexpensive and avoids relying on an unproved formula.
    """
    if modulus == 1:
        return 1

    a %= modulus

    if a == 0:
        raise ArithmeticError(
            f"{a} is not a unit modulo {modulus}."
        )

    x = a

    for n in range(1, 1000000):

        if x == 1:
            return n

        x = (
            x * a
        ) % modulus

    raise ArithmeticError(
        f"Could not determine order of {a} modulo {modulus}."
    )


def complete_orbit_matches(
    rho,
    modulus,
    order,
):
    return [
        k
        for k in range(order)
        if orbit_value_mod(k, modulus) == rho
    ]


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 196 — EXACT 7-ADIC EXPONENT-LIFT / "
        "ORBIT HENSEL AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE RATIO
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. SOURCE RATIO")
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
    # 2. EXPONENT LIFT TABLE
    # ------------------------------------------------------------------

    lift_data = []

    previous_k = None
    previous_order = None

    for e in range(1, MAX_E + 1):

        modulus = P ** e

        rho = ratio_mod(
            Q1,
            Q3,
            modulus,
        )

        order = multiplicative_order(
            MULT,
            modulus,
        )

        matches = complete_orbit_matches(
            rho,
            modulus,
            order,
        )

        if len(matches) != 1:
            raise ArithmeticError(
                f"Expected a unique orbit match at "
                f"7^{e}, got {matches}."
            )

        k = matches[0]

        if previous_k is None:
            lift_multiple = None
        else:
            if previous_order is None:
                raise AssertionError(
                    "Missing previous order."
                )

            delta_k = k - previous_k

            if delta_k % previous_order != 0:
                raise ArithmeticError(
                    f"Exponent {k} does not lift "
                    f"previous exponent {previous_k} "
                    f"modulo previous order {previous_order}."
                )

            lift_multiple = (
                delta_k // previous_order
            )

        exact_value = (
            Q1
            - orbit_value_mod(k, modulus)
            * Q3
        )

        # The expression above uses c_k reduced modulo 7^e, so its
        # valuation as an integer is not meaningful. Compute the exact
        # integer c_k for the valuation diagnostic.
        c_exact = (
            BASE
            * (MULT ** k)
        )

        exact_gap = (
            Q1
            - c_exact * Q3
        )

        lift_data.append(
            {
                "e": e,
                "modulus": modulus,
                "rho": rho,
                "order": order,
                "k": k,
                "matches": matches,
                "previous_k": previous_k,
                "previous_order": previous_order,
                "lift_multiple": lift_multiple,
                "v7_gap": v7(exact_gap),
            }
        )

        previous_k = k
        previous_order = order

    print()
    print("=" * 78)
    print("2. EXACT EXPONENT-LIFT TABLE")
    print("=" * 78)

    for item in lift_data:

        print(
            f"  e={item['e']} "
            f"modulus={item['modulus']} "
            f"rho={item['rho']} "
            f"order={item['order']} "
            f"k_e={item['k']} "
            f"lift_multiple={item['lift_multiple']} "
            f"v7(rho-c_k)={item['v7_gap']}"
        )

    # ------------------------------------------------------------------
    # 3. CONGRUENCE RELATIONS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXPONENT CONGRUENCE RELATIONS")
    print("=" * 78)

    for i in range(1, len(lift_data)):

        prev = lift_data[i - 1]
        curr = lift_data[i]

        prev_order = prev["order"]

        print(
            f"  k_{curr['e']} = "
            f"{curr['k']} "
            f"mod {prev_order}: "
            f"{curr['k'] % prev_order == prev['k'] % prev_order}"
        )

        print(
            f"    previous_k={prev['k']} "
            f"current_k={curr['k']} "
            f"difference={curr['k'] - prev['k']} "
            f"quotient={(curr['k'] - prev['k']) // prev_order}"
        )

    # ------------------------------------------------------------------
    # 4. KNOWN EXPLICIT LIFTS
    # ------------------------------------------------------------------

    expected_k = {
        1: 1,
        2: 1,
        3: 85,
        4: 1555,
    }

    expected_order = {
        1: 6,
        2: 42,
        3: 294,
        4: 2058,
    }

    print()
    print("=" * 78)
    print("4. EXPLICIT LIFT VALUES")
    print("=" * 78)

    for item in lift_data:

        e = item["e"]

        print(
            f"  e={e}: "
            f"k={item['k']} "
            f"expected={expected_k[e]} "
            f"order={item['order']} "
            f"expected_order={expected_order[e]}"
        )

    # ------------------------------------------------------------------
    # 5. MODULAR ORBIT REPRESENTATIVES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. LIFTED ORBIT REPRESENTATIVES")
    print("=" * 78)

    for item in lift_data:

        e = item["e"]
        modulus = item["modulus"]
        k = item["k"]

        c = orbit_value_mod(
            k,
            modulus,
        )

        print(
            f"  e={e}: "
            f"k={k} "
            f"c_k mod {modulus}={c} "
            f"rho={item['rho']} "
            f"equal={c == item['rho']}"
        )

    # ------------------------------------------------------------------
    # 6. INTEGER GAPS AT LIFTED EXPONENTS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. INTEGER SOURCE-RATIO GAPS")
    print("=" * 78)

    for item in lift_data:

        k = item["k"]

        c = (
            BASE
            * (MULT ** k)
        )

        gap = (
            Q1
            - c * Q3
        )

        print(
            f"  k={k}: "
            f"c_k={c} "
            f"gap=q1-c_k*q3={gap} "
            f"v7={v7(gap)}"
        )

    # ------------------------------------------------------------------
    # 7. HENSEL-STYLE DIGIT LIFTS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. HENSEL-STYLE EXPONENT DIGITS")
    print("=" * 78)

    for i in range(1, len(lift_data)):

        prev = lift_data[i - 1]
        curr = lift_data[i]

        digit = (
            curr["lift_multiple"]
        )

        print(
            f"  from e={prev['e']} to e={curr['e']}: "
            f"new_exponent_digit={digit}"
        )

    # ------------------------------------------------------------------
    # 8. INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The source ratio rho determines a unique exponent class k_e at each
tested precision through

    2*3^k_e == rho (mod 7^e).

The exponent therefore lifts through the tower

    mod 7
       ->
    mod 49
       ->
    mod 343
       ->
    mod 2401.

For the present data the lift is

    k_1 = 1,
    k_2 = 1,
    k_3 = 85,
    k_4 = 1555.

The successive orbit orders are

    6,
    42,
    294,
    2058.

Thus the exponent itself has a genuine 7-adic lift structure.

The important distinction is:

    k=1

is the first-order representative, while

    k=85

is the unique representative that matches rho through 3-adic powers
of 7 in the orbit parameter.

The experiment does not claim that this exponent lift extends
indefinitely; it only establishes the exact finite levels tested.
"""
    )

    # ------------------------------------------------------------------
    # 9. FINAL
    # ------------------------------------------------------------------

    expected_profile = [
        (1, 1, 6),
        (2, 1, 42),
        (3, 85, 294),
        (4, 1555, 2058),
    ]

    actual_profile = [
        (
            item["e"],
            item["k"],
            item["order"],
        )
        for item in lift_data
    ]

    lift_consistent = (
        actual_profile
        == expected_profile
    )

    unique_at_all_levels = all(
        len(item["matches"]) == 1
        for item in lift_data
    )

    exact_modular_alignment = all(
        orbit_value_mod(
            item["k"],
            item["modulus"],
        ) == item["rho"]
        for item in lift_data
    )

    final_ok = (
        v7(Q1) == 0
        and v7(Q3) == 0
        and lift_consistent
        and unique_at_all_levels
        and exact_modular_alignment
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  unique_orbit_match_at_all_levels="
        f"{unique_at_all_levels}"
    )

    print(
        f"  exact_modular_alignment="
        f"{exact_modular_alignment}"
    )

    print(
        f"  exponent_lift_profile_exact="
        f"{lift_consistent}"
    )

    print(
        f"  expected_k_sequence="
        f"{[1, 1, 85, 1555]}"
    )

    print(
        f"  actual_k_sequence="
        f"{[item['k'] for item in lift_data]}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 196 COMPLETE")


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

