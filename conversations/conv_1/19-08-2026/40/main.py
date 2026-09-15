#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 207 — EXACT DEFECT-DIGIT -> EXPONENT RECONSTRUCTION AUDIT
==============================================================================

Established by Experiment 206:

    t_e = -delta_e mod 7,

where delta_e is the normalized 7-adic source-ratio defect.

Experiment 207 removes the verified exponent sequence as an input.

It starts from the known base exponent

    k_8 = 3401371

and the directly computed source-defect digits, then reconstructs

    k_(e+1)
      = k_e + (-delta_e mod 7) * ord_{7^e}(3).

For every reconstructed exponent it independently checks

    2*3^k == q1/q3 (mod 7^(e+1)).

Thus the experiment tests whether the defect digits alone are sufficient
to generate the observed 7-adic orbit exponent.

No huge powers of 3 are constructed.
All checks use modular pow().
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

START_E = 8
START_K = 3401371

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

        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s

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


def source_defect_digit(e: int, k: int) -> int:
    """
    Compute

        delta_e =
            (2*3^k - rho) / 7^e mod 7

    using only arithmetic modulo 7^(e+1).
    """
    modulus_e = P ** e
    modulus_next = P ** (e + 1)

    rho_next = ratio_mod(
        modulus_next
    )

    orbit_next = orbit_mod(
        k,
        modulus_next
    )

    defect = (
        orbit_next
        - rho_next
    ) % modulus_next

    if defect % modulus_e != 0:
        raise ArithmeticError(
            f"Defect is not divisible by 7^{e}."
        )

    return (
        defect // modulus_e
    ) % P


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 207 — EXACT DEFECT-DIGIT -> "
        "EXPONENT RECONSTRUCTION AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. SOURCE DATA")
    print("=" * 78)

    print(
        f"  q1={Q1}"
    )

    print(
        f"  q3={Q3}"
    )

    print(
        f"  v7(q1)={valuation_7(Q1)}"
    )

    print(
        f"  v7(q3)={valuation_7(Q3)}"
    )

    # ------------------------------------------------------------------
    # 2. STARTING POINT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. STARTING EXPONENT")
    print("=" * 78)

    print(
        f"  start_e={START_E}"
    )

    print(
        f"  start_k={START_K}"
    )

    start_modulus = P ** START_E

    start_rho = ratio_mod(
        start_modulus
    )

    start_orbit = orbit_mod(
        START_K,
        start_modulus
    )

    print(
        f"  rho_mod_7^e={start_rho}"
    )

    print(
        f"  orbit_mod_7^e={start_orbit}"
    )

    print(
        f"  start_alignment="
        f"{start_rho == start_orbit}"
    )

    # ------------------------------------------------------------------
    # 3. DIRECT DEFECT DIGITS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. DIRECTLY COMPUTED DEFECT DIGITS")
    print("=" * 78)

    direct_digits = []

    temp_k = START_K

    for e in range(
        START_E,
        END_E,
    ):

        delta = source_defect_digit(
            e,
            temp_k
        )

        direct_digits.append(
            delta
        )

        print(
            f"  e={e}: "
            f"k={temp_k} "
            f"delta_e={delta} "
            f"lift_digit={(-delta) % P}"
        )

        order = order_3_mod_7e(
            e
        )

        temp_k += (
            ((-delta) % P)
            * order
        )

    # ------------------------------------------------------------------
    # 4. RECONSTRUCT EXPONENTS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. EXPONENT RECONSTRUCTION FROM DEFECT DIGITS")
    print("=" * 78)

    reconstructed = [
        START_K
    ]

    current_k = START_K

    transition_rows = []

    for i, e in enumerate(
        range(
            START_E,
            END_E,
        )
    ):

        delta = direct_digits[i]

        t = (
            -delta
        ) % P

        order = order_3_mod_7e(
            e
        )

        increment = (
            t * order
        )

        next_k = (
            current_k
            + increment
        )

        transition_rows.append(
            {
                "e": e,
                "current_k": current_k,
                "delta": delta,
                "t": t,
                "order": order,
                "increment": increment,
                "next_k": next_k,
            }
        )

        reconstructed.append(
            next_k
        )

        print(
            f"  e={e}: "
            f"current_k={current_k} "
            f"delta={delta} "
            f"t={t} "
            f"order={order} "
            f"increment={increment} "
            f"next_k={next_k}"
        )

        current_k = next_k

    # ------------------------------------------------------------------
    # 5. INDEPENDENT MODULAR VALIDATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. INDEPENDENT MODULAR VALIDATION")
    print("=" * 78)

    validation_flags = []

    for i, e in enumerate(
        range(
            START_E,
            END_E + 1,
        )
    ):

        k = reconstructed[i]

        modulus = P ** e

        rho = ratio_mod(
            modulus
        )

        orbit = orbit_mod(
            k,
            modulus
        )

        exact = (
            orbit
            == rho
        )

        validation_flags.append(
            exact
        )

        print(
            f"  e={e}: "
            f"k={k} "
            f"orbit={orbit} "
            f"rho={rho} "
            f"exact={exact}"
        )

    # ------------------------------------------------------------------
    # 6. NEXT-PRECISION VALIDATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. TRANSITION-LEVEL VALIDATION")
    print("=" * 78)

    transition_flags = []

    for row in transition_rows:

        next_e = row["e"] + 1
        k = row["next_k"]

        modulus = P ** next_e

        rho = ratio_mod(
            modulus
        )

        orbit = orbit_mod(
            k,
            modulus
        )

        exact = (
            orbit
            == rho
        )

        transition_flags.append(
            exact
        )

        print(
            f"  e={row['e']} -> "
            f"{next_e}: "
            f"k={k} "
            f"orbit={orbit} "
            f"rho={rho} "
            f"exact={exact}"
        )

    # ------------------------------------------------------------------
    # 7. RECONSTRUCTION UNIQUENESS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. DEFECT-DIGIT RECONSTRUCTION UNIQUENESS")
    print("=" * 78)

    uniqueness_flags = []

    for row in transition_rows:

        e = row["e"]
        current_k = row["current_k"]
        target_k = row["next_k"]
        order = row["order"]

        matching_digits = []

        for t in range(P):

            candidate = (
                current_k
                + t * order
            )

            modulus = P ** (
                e + 1
            )

            if (
                orbit_mod(
                    candidate,
                    modulus
                )
                ==
                ratio_mod(
                    modulus
                )
            ):
                matching_digits.append(
                    t
                )

        unique = (
            len(matching_digits)
            == 1
            and matching_digits[0]
            == row["t"]
        )

        uniqueness_flags.append(
            unique
        )

        print(
            f"  e={e}: "
            f"matching_digits={matching_digits} "
            f"chosen={row['t']} "
            f"unique={unique}"
        )

    # ------------------------------------------------------------------
    # 8. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The previous experiment established

    t_e = -delta_e mod 7.

Experiment 207 reverses the direction of the calculation.

Instead of taking the exponent and reading off its lift digit, it takes
the source-ratio defect digit and constructs the next exponent:

    delta_e
        ->
    t_e = -delta_e
        ->
    k_(e+1)
        = k_e + t_e ord_{7^e}(3).

The resulting exponent is then checked independently against

    2*3^k == q1/q3 (mod 7^(e+1)).

Therefore the defect digits are sufficient to reconstruct the entire
verified exponent sequence from the single starting exponent k_8.

This is an exact finite inverse-limit reconstruction statement.
"""
    )

    # ------------------------------------------------------------------
    # 9. FINAL
    # ------------------------------------------------------------------

    expected_exponents = [
        3401371,
        3401371,
        72578983,
        798943909,
        10968052873,
    ]

    exponent_exact = (
        reconstructed
        == expected_exponents
    )

    all_validation = all(
        validation_flags
    )

    all_transitions = all(
        transition_flags
    )

    all_unique = all(
        uniqueness_flags
    )

    final_ok = (
        valuation_7(Q1) == 0
        and valuation_7(Q3) == 0
        and exponent_exact
        and all_validation
        and all_transitions
        and all_unique
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  defect_digits="
        f"{direct_digits}"
    )

    print(
        f"  reconstructed_exponents="
        f"{reconstructed}"
    )

    print(
        f"  expected_exponents="
        f"{expected_exponents}"
    )

    print(
        f"  exponent_reconstruction_exact="
        f"{exponent_exact}"
    )

    print(
        f"  all_modular_validations_exact="
        f"{all_validation}"
    )

    print(
        f"  all_transition_checks_exact="
        f"{all_transitions}"
    )

    print(
        f"  unique_digit_at_every_level="
        f"{all_unique}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 207 COMPLETE")


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

