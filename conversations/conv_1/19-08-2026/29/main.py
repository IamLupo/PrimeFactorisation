#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 198 — EXACT 7-ADIC RESIDUAL-DIGIT /
                LINEARIZED LIFT AUDIT
==============================================================================

Experiment 197R established the unique exponent lift

    k1 = 1
    k2 = 1
    k3 = 85
    k4 = 1555

with digits

    t = [0, 2, 5].

At level e, suppose

    F(k) = 2*3^k*q3 - q1.

If k_e solves

    F(k_e) == 0 (mod 7^e),

then the seven lifts

    k_e + t*L_e,

where

    L_e = ord_{7^e}(3),

must be tested.

Experiment 198 goes one step further.

For every candidate lift it computes the normalized residual

    F(k_e + t L_e) / 7^e  (mod 7)

whenever the division is exact at that level.

The seven normalized residues provide the exact finite "digit equation".

The successful lift corresponds to residue zero.

The experiment also verifies that the dependence on t is affine modulo 7:

    normalized_residual(t)
        == A + t*B (mod 7),

and extracts A and B directly from the candidate table.

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


def orbit_value_mod(
    k: int,
    modulus: int,
) -> int:
    return (
        BASE * pow(MULT, k, modulus)
    ) % modulus


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

    x = a

    for n in range(1, 1_000_000):

        if x == 1:
            return n

        x = (
            x * a
        ) % modulus

    raise ArithmeticError(
        f"Order search failed for {a} modulo {modulus}."
    )


def unique_exponent(
    rho: int,
    modulus: int,
    order: int,
) -> int:

    matches = []

    for k in range(order):

        if (
            orbit_value_mod(
                k,
                modulus,
            )
            == rho
        ):
            matches.append(k)

    if len(matches) != 1:
        raise ArithmeticError(
            f"Expected one exponent modulo {modulus}, "
            f"got {matches}."
        )

    return matches[0]


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 198 — EXACT 7-ADIC RESIDUAL-DIGIT / "
        "LINEARIZED LIFT AUDIT"
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
    print(f"  v7(q1)={v7(Q1)}")
    print(f"  v7(q3)={v7(Q3)}")

    # ------------------------------------------------------------------
    # 2. BASE EXPONENTS
    # ------------------------------------------------------------------

    levels = []

    print()
    print("=" * 78)
    print("2. EXACT LEVEL DATA")
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

        k = unique_exponent(
            rho,
            modulus,
            order,
        )

        c = orbit_value(k)

        residual = (
            c * Q3 - Q1
        )

        levels.append(
            {
                "e": e,
                "modulus": modulus,
                "rho": rho,
                "order": order,
                "k": k,
                "c": c,
                "residual": residual,
                "v7": v7(residual),
            }
        )

        print(
            f"  e={e}: "
            f"modulus={modulus} "
            f"order={order} "
            f"k={k} "
            f"v7(F)={v7(residual)}"
        )

    # ------------------------------------------------------------------
    # 3. CANDIDATE RESIDUAL-DIGIT TABLE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. NORMALIZED RESIDUAL DIGITS")
    print("=" * 78)

    lift_records = []

    for i in range(
        len(levels) - 1
    ):

        current = levels[i]
        target = levels[i + 1]

        e = current["e"]
        k = current["k"]
        L = current["order"]

        print()
        print(
            f"  e={e} -> e={e+1}: "
            f"k={k}, L={L}"
        )

        candidate_rows = []

        for t in range(7):

            candidate_k = (
                k + t * L
            )

            candidate_c = orbit_value(
                candidate_k
            )

            F = (
                candidate_c * Q3
                - Q1
            )

            divisible = (
                F % (P ** e) == 0
            )

            if not divisible:
                raise ArithmeticError(
                    f"Candidate residual not divisible by 7^{e}: "
                    f"t={t}, k={candidate_k}."
                )

            normalized = (
                (F // (P ** e))
                % P
            )

            target_match = (
                candidate_c
                % target["modulus"]
                == target["rho"]
            )

            candidate_rows.append(
                {
                    "t": t,
                    "k": candidate_k,
                    "F": F,
                    "normalized": normalized,
                    "target_match": target_match,
                }
            )

            print(
                f"    t={t}: "
                f"k={candidate_k} "
                f"normalized={normalized} "
                f"target_match={target_match}"
            )

        lift_records.append(
            {
                "e": e,
                "current_k": k,
                "order": L,
                "target_k": target["k"],
                "rows": candidate_rows,
            }
        )

    # ------------------------------------------------------------------
    # 4. AFFINE MOD-7 TEST
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. AFFINE RESIDUAL-DIGIT TEST")
    print("=" * 78)

    affine_checks = []
    chosen_digits = []

    for record in lift_records:

        rows = record["rows"]

        r0 = rows[0]["normalized"]
        r1 = rows[1]["normalized"]

        B = (
            r1 - r0
        ) % 7

        A = r0

        affine = True

        for row in rows:

            t = row["t"]

            predicted = (
                A
                + t * B
            ) % 7

            if predicted != row["normalized"]:
                affine = False

        affine_checks.append(
            affine
        )

        zero_digits = [
            row["t"]
            for row in rows
            if row["normalized"] == 0
        ]

        if len(zero_digits) != 1:
            raise ArithmeticError(
                f"Expected exactly one zero residual digit, "
                f"got {zero_digits}."
            )

        chosen_digits.append(
            zero_digits[0]
        )

        print(
            f"  e={record['e']}: "
            f"A={A} "
            f"B={B} "
            f"affine={affine} "
            f"zero_digit={zero_digits[0]}"
        )

    # ------------------------------------------------------------------
    # 5. UNIQUE DIGIT / LINEAR COEFFICIENT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. LINEARIZED DIGIT DATA")
    print("=" * 78)

    for i, record in enumerate(
        lift_records
    ):

        rows = record["rows"]

        A = rows[0]["normalized"]

        B = (
            rows[1]["normalized"]
            - rows[0]["normalized"]
        ) % 7

        chosen = chosen_digits[i]

        print(
            f"  e={record['e']}: "
            f"residual(t)=A+tB mod7 "
            f"with A={A}, B={B}, "
            f"chosen_t={chosen}"
        )

        if B % 7 == 0:
            print(
                "    WARNING: zero linear coefficient."
            )

        else:

            check_zero = (
                A
                + chosen * B
            ) % 7

            print(
                f"    A+tB at chosen_t={check_zero}"
            )

    # ------------------------------------------------------------------
    # 6. RESIDUAL VALUATION AFTER LIFT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. CHOSEN-LIFT VALUATION")
    print("=" * 78)

    valuation_checks = []

    for i, record in enumerate(
        lift_records
    ):

        chosen_t = chosen_digits[i]

        chosen_row = next(
            row
            for row in record["rows"]
            if row["t"] == chosen_t
        )

        chosen_v = v7(
            chosen_row["F"]
        )

        target_level = (
            record["e"] + 1
        )

        exact = (
            chosen_v is not None
            and chosen_v >= target_level
        )

        valuation_checks.append(
            exact
        )

        print(
            f"  e={record['e']}: "
            f"chosen_t={chosen_t} "
            f"v7(F)={chosen_v} "
            f"required>={target_level} "
            f"exact={exact}"
        )

    # ------------------------------------------------------------------
    # 7. LIFT DIGIT RECONSTRUCTION
    # ------------------------------------------------------------------

    predicted_exponents = []

    for i, record in enumerate(
        lift_records
    ):

        predicted = (
            record["current_k"]
            + chosen_digits[i]
            * record["order"]
        )

        predicted_exponents.append(
            predicted
        )

    print()
    print("=" * 78)
    print("7. EXPONENT RECONSTRUCTION")
    print("=" * 78)

    for i, record in enumerate(
        lift_records
    ):

        predicted = predicted_exponents[i]

        print(
            f"  from e={record['e']}: "
            f"predicted_next_k={predicted} "
            f"actual={record['target_k']} "
            f"exact={predicted == record['target_k']}"
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
Experiment 197R identified the successful lift digit by brute-force
enumeration of seven candidates.

Experiment 198 extracts the residual equation itself.

After fixing k_e, every lift has the form

    k_e + t*L_e.

Modulo the next power of 7, the normalized residual is an affine
function of t:

    R_e(t) = A_e + B_e t (mod 7).

The unique zero of this equation is exactly the next Hensel lift digit.

Thus the successful exponent can be reconstructed from the local
residual rather than simply selected by search.

The resulting digits are the exact finite Hensel data for the source
ratio orbit.
"""
    )

    # ------------------------------------------------------------------
    # 9. FINAL
    # ------------------------------------------------------------------

    expected_digits = [0, 2, 5]
    expected_exponents = [1, 1, 85, 1555]

    actual_exponents = [
        levels[0]["k"]
    ] + predicted_exponents

    affine_exact = all(
        affine_checks
    )

    valuation_exact = all(
        valuation_checks
    )

    exponent_exact = (
        actual_exponents
        == expected_exponents
    )

    digits_exact = (
        chosen_digits
        == expected_digits
    )

    final_ok = (
        v7(Q1) == 0
        and v7(Q3) == 0
        and affine_exact
        and valuation_exact
        and exponent_exact
        and digits_exact
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  affine_residual_law_exact="
        f"{affine_exact}"
    )

    print(
        f"  unique_zero_digit_at_each_level="
        f"{all(len([r for r in record['rows'] if r['normalized'] == 0]) == 1 for record in lift_records)}"
    )

    print(
        f"  chosen_lift_valuation_exact="
        f"{valuation_exact}"
    )

    print(
        f"  exponent_reconstruction_exact="
        f"{exponent_exact}"
    )

    print(
        f"  lift_digits="
        f"{chosen_digits}"
    )

    print(
        f"  exponent_sequence="
        f"{actual_exponents}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 198 COMPLETE")


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

