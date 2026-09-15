#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 209 — EXACT CANONICAL 7-ADIC EXPONENT DIGIT EXPANSION AUDIT
==============================================================================

Experiment 208 showed that the source ratio alone reconstructs the
entire finite exponent sequence, but its final PASS flag was spoiled
by a bad hard-coded expected list.

Experiment 209 removes ALL hard-coded exponent expectations.

It reconstructs the exponent sequence only from rho=q1/q3 and records
the canonical Hensel digits

    t_e in {0,...,6}

defined by

    k_(e+1)
      = k_e + t_e * ord_{7^e}(3).

Since

    ord_{7^e}(3) = 6*7^(e-1),

we have

    k_E
      = k_1
        + 6 * sum_{j=1}^{E-1} t_j 7^(j-1).

The experiment verifies this identity independently in two ways:

    1. recursive exponent lifting;
    2. direct base-7 digit reconstruction.

No reference exponent sequence is supplied.
No SymPy.
No giant powers of 3.
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
            f"Expected one base exponent, got {matches}."
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
            f"Exponent k={k} does not solve modulo 7^{e}."
        )

    return (
        diff // modulus_e
    ) % P


def reconstruct_from_digits(
    k1: int,
    digits: list[int],
    target_e: int,
) -> int:

    """
    Direct formula:

        k_E = k1 + 6 * sum_{j=1}^{E-1} t_j 7^(j-1)

    digits[j] corresponds to t_(j+1), with the first transition at e=1.
    """

    value = k1

    for j, t in enumerate(digits):

        power = P ** j

        value += (
            6
            * t
            * power
        )

    return value


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 209 — EXACT CANONICAL 7-ADIC "
        "EXPONENT DIGIT EXPANSION AUDIT"
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
    print(f"  q3_mod7={Q3 % P}")
    print(f"  v7(q1)={valuation_7(Q1)}")
    print(f"  v7(q3)={valuation_7(Q3)}")

    # ------------------------------------------------------------------
    # 2. BASE EXPONENT
    # ------------------------------------------------------------------

    k1 = base_exponent_mod7()

    print()
    print("=" * 78)
    print("2. BASE EXPONENT CLASS")
    print("=" * 78)

    print(
        f"  rho_mod7={ratio_mod(7)}"
    )

    print(
        f"  k1={k1}"
    )

    print(
        f"  orbit_mod7={orbit_mod(k1, 7)}"
    )

    # ------------------------------------------------------------------
    # 3. RECURSIVE CANONICAL LIFT
    # ------------------------------------------------------------------

    current_k = k1

    exponent_sequence = [k1]
    transition_digits = []
    transition_rows = []

    print()
    print("=" * 78)
    print("3. CANONICAL RECURSIVE LIFT")
    print("=" * 78)

    for e in range(1, END_E):

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

        alignment = (
            orbit_mod(
                next_k,
                modulus_next
            )
            ==
            ratio_mod(
                modulus_next
            )
        )

        if not alignment:
            raise ArithmeticError(
                f"Lift failed at e={e}."
            )

        transition_digits.append(t)
        exponent_sequence.append(next_k)

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
    # 4. DIRECT DIGIT RECONSTRUCTION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. DIRECT BASE-7 EXPONENT RECONSTRUCTION")
    print("=" * 78)

    direct_sequence = []

    for E in range(
        1,
        END_E + 1
    ):

        digits_for_E = (
            transition_digits[: E - 1]
        )

        direct_k = reconstruct_from_digits(
            k1,
            digits_for_E,
            E
        )

        direct_sequence.append(
            direct_k
        )

        recursive_k = exponent_sequence[
            E - 1
        ]

        exact = (
            direct_k == recursive_k
        )

        print(
            f"  E={E}: "
            f"direct={direct_k} "
            f"recursive={recursive_k} "
            f"exact={exact}"
        )

    # ------------------------------------------------------------------
    # 5. EXPONENT DIGIT BLOCK
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. CANONICAL HENSEL DIGITS")
    print("=" * 78)

    print(
        f"  t_digits={transition_digits}"
    )

    print(
        f"  defect_digits="
        f"{[(-t) % P for t in transition_digits]}"
    )

    print(
        f"  digit_range_valid="
        f"{all(0 <= t < P for t in transition_digits)}"
    )

    # ------------------------------------------------------------------
    # 6. REPEATED-EXPONENT EVENTS
    # ------------------------------------------------------------------

    repeats = []

    for i in range(
        len(exponent_sequence) - 1
    ):

        if (
            exponent_sequence[i]
            == exponent_sequence[i + 1]
        ):
            repeats.append(
                {
                    "from_e": i + 1,
                    "to_e": i + 2,
                    "k": exponent_sequence[i],
                    "t": transition_digits[i],
                }
            )

    print()
    print("=" * 78)
    print("6. REPEATED EXPONENT EVENTS")
    print("=" * 78)

    print(
        f"  repeats={repeats}"
    )

    # ------------------------------------------------------------------
    # 7. GLOBAL ORBIT VALIDATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. GLOBAL ORBIT VALIDATION")
    print("=" * 78)

    orbit_checks = []

    for E, k in enumerate(
        exponent_sequence,
        start=1
    ):

        modulus = P ** E

        orbit = orbit_mod(
            k,
            modulus
        )

        rho = ratio_mod(
            modulus
        )

        exact = (
            orbit == rho
        )

        orbit_checks.append(
            exact
        )

        print(
            f"  E={E}: "
            f"k={k} "
            f"orbit={orbit} "
            f"rho={rho} "
            f"exact={exact}"
        )

    # ------------------------------------------------------------------
    # 8. UNIQUENESS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. LIFT DIGIT UNIQUENESS")
    print("=" * 78)

    uniqueness_flags = []

    for row in transition_rows:

        e = row["e"]
        current_k = row["current_k"]
        order = row["order"]
        chosen_t = row["t"]

        modulus = P ** (
            e + 1
        )

        matches = []

        for t in range(P):

            candidate = (
                current_k
                + t * order
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
                matches.append(t)

        unique = (
            len(matches) == 1
            and matches[0] == chosen_t
        )

        uniqueness_flags.append(
            unique
        )

        print(
            f"  e={e}: "
            f"matches={matches} "
            f"chosen={chosen_t} "
            f"unique={unique}"
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
The finite exponent sequence is now reconstructed in two completely
different ways.

Recursive form:

    k_(e+1)
      = k_e + t_e * 6*7^(e-1).

Closed digit form:

    k_E
      = k_1
        + 6 * sum_{j=1}^{E-1} t_j 7^(j-1).

The experiment starts only from rho and determines k_1 from the
mod-7 orbit. Every subsequent digit is obtained from the source-ratio
defect.

Therefore the exponent is recovered from the source ratio as a
canonical base-7 digit expansion in the exponent coordinate.

Repeated exponents occur exactly when the corresponding digit t_e is
zero.

No exponent sequence is hard-coded into the experiment.
"""
    )

    # ------------------------------------------------------------------
    # 10. FINAL
    # ------------------------------------------------------------------

    all_direct = all(
        a == b
        for a, b in zip(
            direct_sequence,
            exponent_sequence
        )
    )

    all_orbit = all(
        orbit_checks
    )

    all_unique = all(
        uniqueness_flags
    )

    final_ok = (
        valuation_7(Q1) == 0
        and valuation_7(Q3) == 0
        and all_direct
        and all_orbit
        and all_unique
    )

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  recursive_sequence="
        f"{exponent_sequence}"
    )

    print(
        f"  direct_digit_sequence="
        f"{direct_sequence}"
    )

    print(
        f"  recursive_equals_direct="
        f"{all_direct}"
    )

    print(
        f"  all_orbit_checks_exact="
        f"{all_orbit}"
    )

    print(
        f"  unique_digit_every_level="
        f"{all_unique}"
    )

    print(
        f"  t_digits="
        f"{transition_digits}"
    )

    print(
        f"  defect_digits="
        f"{[(-t) % P for t in transition_digits]}"
    )

    print(
        f"  repeated_exponent_count="
        f"{len(repeats)}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 209 COMPLETE")


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

