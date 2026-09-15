#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 197R — EXACT 7-ADIC EXPONENT-LIFT RECURRENCE AUDIT
==============================================================================

Experiment 196 established the unique exponent sequence satisfying

    2*3^k == q1/q3 (mod 7^e)

through e=1,...,4:

    k_1 = 1
    k_2 = 1
    k_3 = 85
    k_4 = 1555.

The successive lift digits are

    0, 2, 5.

Experiment 197R verifies these digits constructively.

At level e, with exponent k_e and

    ord_{7^e}(3) = L_e,

all seven possible lifts are

    k_e + t L_e,     t in {0,...,6}.

Exactly one candidate must satisfy the next congruence modulo 7^(e+1).

The script explicitly enumerates all seven candidates, records their
7-adic residual valuations, and verifies that the unique successful
candidate is the next exponent.

No SymPy.
No floating point.
No nested f-strings.
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

def v7(x: int):
    if x == 0:
        return None

    x = abs(x)
    e = 0

    while x % 7 == 0:
        x //= 7
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

    while r != 0:
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


def orbit_value(k: int) -> int:
    return (
        BASE
        * (MULT ** k)
    )


def orbit_value_mod(
    k: int,
    modulus: int,
) -> int:
    return (
        BASE
        * pow(MULT, k, modulus)
    ) % modulus


def unique_orbit_match(
    rho: int,
    modulus: int,
    order: int,
) -> int:

    matches = []

    for k in range(order):
        value = orbit_value_mod(
            k,
            modulus,
        )

        if value == rho:
            matches.append(k)

    if len(matches) != 1:
        raise ArithmeticError(
            f"Expected exactly one match modulo {modulus}, "
            f"got {matches}."
        )

    return matches[0]


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 197R — EXACT 7-ADIC EXPONENT-LIFT "
        "RECURRENCE AUDIT"
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
    # 2. FIND EXPONENTS AT EACH PRECISION
    # ------------------------------------------------------------------

    exponent_data = []

    print()
    print("=" * 78)
    print("2. EXACT EXPONENT SEQUENCE")
    print("=" * 78)

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

        k = unique_orbit_match(
            rho,
            modulus,
            order,
        )

        c = orbit_value(k)

        residual = (
            Q1
            - c * Q3
        )

        record = {
            "e": e,
            "modulus": modulus,
            "rho": rho,
            "order": order,
            "k": k,
            "c": c,
            "residual": residual,
            "v7": v7(residual),
        }

        exponent_data.append(
            record
        )

        print(
            f"  e={e}: "
            f"modulus={modulus} "
            f"order={order} "
            f"k={k} "
            f"rho={rho} "
            f"v7(residual)={v7(residual)}"
        )

    # ------------------------------------------------------------------
    # 3. ENUMERATE ALL SEVEN LIFTS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. COMPLETE SEVEN-CANDIDATE LIFT AUDIT")
    print("=" * 78)

    lift_records = []

    for i in range(
        len(exponent_data) - 1
    ):

        current = exponent_data[i]
        target = exponent_data[i + 1]

        current_e = current["e"]
        current_k = current["k"]
        previous_order = current["order"]

        candidate_rows = []

        for t in range(7):

            candidate_k = (
                current_k
                + t * previous_order
            )

            candidate_c = orbit_value(
                candidate_k
            )

            candidate_residual = (
                Q1
                - candidate_c * Q3
            )

            candidate_v7 = v7(
                candidate_residual
            )

            target_match = (
                orbit_value_mod(
                    candidate_k,
                    target["modulus"],
                )
                == target["rho"]
            )

            row = {
                "t": t,
                "k": candidate_k,
                "v7": candidate_v7,
                "target_match": target_match,
            }

            candidate_rows.append(
                row
            )

            print(
                f"  e={current_e}->e={current_e+1}: "
                f"t={t} "
                f"k={candidate_k} "
                f"v7={candidate_v7} "
                f"target_match={target_match}"
            )

        successful = [
            row
            for row in candidate_rows
            if row["target_match"]
        ]

        if len(successful) != 1:
            raise ArithmeticError(
                f"Expected one successful lift from e={current_e}, "
                f"got {successful}."
            )

        chosen = successful[0]

        lift_records.append(
            {
                "from_e": current_e,
                "to_e": current_e + 1,
                "current_k": current_k,
                "previous_order": previous_order,
                "chosen_t": chosen["t"],
                "chosen_k": chosen["k"],
                "target_k": target["k"],
                "candidate_rows": candidate_rows,
                "successful_count": len(successful),
            }
        )

    # ------------------------------------------------------------------
    # 4. RECURRENCE CHECK
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. EXACT EXPONENT RECURRENCE")
    print("=" * 78)

    recurrence_checks = []

    for record in lift_records:

        predicted = (
            record["current_k"]
            + record["chosen_t"]
            * record["previous_order"]
        )

        exact = (
            predicted
            == record["target_k"]
        )

        recurrence_checks.append(
            exact
        )

        print(
            f"  k_{record['to_e']} = "
            f"k_{record['from_e']} + "
            f"{record['chosen_t']}*"
            f"{record['previous_order']} "
            f"= {predicted}; "
            f"target={record['target_k']} "
            f"exact={exact}"
        )

    # ------------------------------------------------------------------
    # 5. UNIQUE DIGIT CHECK
    # ------------------------------------------------------------------

    unique_lift_checks = [
        record["successful_count"] == 1
        for record in lift_records
    ]

    unique_lift_exact = all(
        unique_lift_checks
    )

    print()
    print("=" * 78)
    print("5. UNIQUE LIFT-DIGIT CHECK")
    print("=" * 78)

    for record in lift_records:

        print(
            f"  e={record['from_e']}->e={record['to_e']}: "
            f"successful_candidates="
            f"{record['successful_count']} "
            f"unique="
            f"{record['successful_count'] == 1}"
        )

    # ------------------------------------------------------------------
    # 6. VALUATION IMPROVEMENT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. RESIDUAL VALUATION IMPROVEMENT")
    print("=" * 78)

    valuation_checks = []

    for i in range(
        len(exponent_data) - 1
    ):

        current = exponent_data[i]
        target = exponent_data[i + 1]

        current_v = current["v7"]
        target_v = target["v7"]

        required = (
            current["e"] + 1
        )

        exact = (
            target_v is not None
            and target_v >= required
        )

        valuation_checks.append(
            exact
        )

        print(
            f"  e={current['e']}->"
            f"{current['e']+1}: "
            f"current_v7={current_v} "
            f"next_v7={target_v} "
            f"required>={required} "
            f"exact={exact}"
        )

    valuation_gain_exact = all(
        valuation_checks
    )

    # ------------------------------------------------------------------
    # 7. MODULAR ALIGNMENT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. MODULAR ALIGNMENT")
    print("=" * 78)

    modular_checks = []

    for record in exponent_data:

        lhs = orbit_value_mod(
            record["k"],
            record["modulus"],
        )

        rhs = record["rho"]

        exact = (
            lhs == rhs
        )

        modular_checks.append(
            exact
        )

        print(
            f"  e={record['e']}: "
            f"2*3^k mod {record['modulus']}="
            f"{lhs}; "
            f"rho={rhs}; "
            f"exact={exact}"
        )

    modular_exact = all(
        modular_checks
    )

    # ------------------------------------------------------------------
    # 8. LIFT DIGITS
    # ------------------------------------------------------------------

    lift_digits = [
        record["chosen_t"]
        for record in lift_records
    ]

    exponent_sequence = [
        record["k"]
        for record in exponent_data
    ]

    print()
    print("=" * 78)
    print("8. LIFT DIGITS")
    print("=" * 78)

    print(
        f"  lift_digits={lift_digits}"
    )

    print(
        f"  exponent_sequence={exponent_sequence}"
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
At each precision level, the exponent k_e has seven possible lifts:

    k_e + t * ord_{7^e}(3),
        t in {0,...,6}.

The experiment explicitly checks all seven possibilities.

Exactly one candidate matches the next modulus.

For the present data:

    k_1 = 1
    k_2 = 1
    k_3 = 85
    k_4 = 1555

with lift digits

    0, 2, 5.

Thus the exponent sequence is not merely obtained by independent
searches. Each higher-precision exponent is the unique one-step lift
of the previous exponent.

The residual valuation simultaneously increases enough to certify the
next 7-adic precision.

This is a finite exact Hensel-style lifting statement for the orbit
equation

    2*3^k = q1/q3.
"""
    )

    # ------------------------------------------------------------------
    # 10. FINAL
    # ------------------------------------------------------------------

    expected_exponents = [
        1,
        1,
        85,
        1555,
    ]

    expected_digits = [
        0,
        2,
        5,
    ]

    exponent_sequence_exact = (
        exponent_sequence
        == expected_exponents
    )

    lift_digits_exact = (
        lift_digits
        == expected_digits
    )

    final_ok = (
        v7(Q1) == 0
        and v7(Q3) == 0
        and recurrence_checks
        and all(recurrence_checks)
        and unique_lift_exact
        and valuation_gain_exact
        and modular_exact
        and exponent_sequence_exact
        and lift_digits_exact
    )

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  recurrence_exact={all(recurrence_checks)}"
    )

    print(
        f"  unique_lift_digit_at_each_level="
        f"{unique_lift_exact}"
    )

    print(
        f"  valuation_gain_exact="
        f"{valuation_gain_exact}"
    )

    print(
        f"  modular_alignment_exact="
        f"{modular_exact}"
    )

    print(
        f"  exponent_sequence_exact="
        f"{exponent_sequence_exact}"
    )

    print(
        f"  lift_digits_exact="
        f"{lift_digits_exact}"
    )

    print(
        f"  lift_digits={lift_digits}"
    )

    print(
        f"  exponent_sequence={exponent_sequence}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 197R COMPLETE")


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