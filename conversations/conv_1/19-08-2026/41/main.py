#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 208 — EXACT SOURCE-RATIO -> FULL EXPONENT-LIFT RECONSTRUCTION
==============================================================================

Experiment 207 established that, once k_8 is known,

    delta_e -> t_e=-delta_e mod 7
            -> k_(e+1)

reconstructs the verified exponent sequence exactly.

Experiment 208 removes k_8 as an input.

The reconstruction starts only from the source ratio

    rho = q1/q3

and the orbit

    c_k = 2*3^k.

Step 1:
    Solve the complete orbit modulo 7 and obtain the unique base class k_1.

Step 2:
    Recover the defect digit at level e.

Step 3:
    Convert it to the Hensel digit

        t_e = -delta_e mod 7.

Step 4:
    Lift the exponent

        k_(e+1)
          = k_e + t_e * ord_{7^e}(3).

The exponent is therefore reconstructed from rho alone.

No previously known exponent sequence is supplied.

All calculations remain modular.
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
        * inverse_mod(
            Q3,
            modulus,
        )
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


def find_base_exponent_mod7() -> int:
    """
    Find the unique k in {0,...,5} satisfying

        2*3^k == rho (mod 7).
    """
    rho = ratio_mod(7)

    matches = []

    for k in range(6):

        if (
            orbit_mod(k, 7)
            == rho
        ):
            matches.append(k)

    if len(matches) != 1:
        raise ArithmeticError(
            f"Expected exactly one base exponent modulo 7, "
            f"got {matches}."
        )

    return matches[0]


def defect_digit(e: int, k: int) -> int:
    """
    Compute

        delta_e =
            (2*3^k - rho) / 7^e mod 7

    from residues modulo 7^(e+1).
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

    diff = (
        orbit_next
        - rho_next
    ) % modulus_next

    if diff % modulus_e != 0:
        raise ArithmeticError(
            f"At e={e}, current exponent does not solve "
            f"the congruence modulo 7^{e}."
        )

    return (
        diff // modulus_e
    ) % P


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 208 — EXACT SOURCE-RATIO -> "
        "FULL EXPONENT-LIFT RECONSTRUCTION"
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
        f"  q3_mod7={Q3 % P}"
    )

    print(
        f"  v7(q1)={valuation_7(Q1)}"
    )

    print(
        f"  v7(q3)={valuation_7(Q3)}"
    )

    # ------------------------------------------------------------------
    # 2. BASE EXPONENT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. BASE ORBIT CLASS MODULO 7")
    print("=" * 78)

    rho7 = ratio_mod(7)

    k1 = find_base_exponent_mod7()

    print(
        f"  rho_mod7={rho7}"
    )

    print(
        f"  unique_k1={k1}"
    )

    print(
        f"  orbit_k1_mod7="
        f"{orbit_mod(k1, 7)}"
    )

    print(
        f"  exact_base_alignment="
        f"{orbit_mod(k1, 7) == rho7}"
    )

    # ------------------------------------------------------------------
    # 3. BUILD THE FULL LIFT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. FULL EXPONENT RECONSTRUCTION")
    print("=" * 78)

    reconstructed = [
        k1
    ]

    current_k = k1

    transition_rows = []

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

        increment = (
            t
            * order
        )

        next_k = (
            current_k
            + increment
        )

        modulus_next = P ** (
            e + 1
        )

        alignment = (
            orbit_mod(
                next_k,
                modulus_next,
            )
            ==
            ratio_mod(
                modulus_next
            )
        )

        transition_rows.append(
            {
                "e": e,
                "current_k": current_k,
                "order": order,
                "delta": delta,
                "t": t,
                "increment": increment,
                "next_k": next_k,
                "alignment": alignment,
            }
        )

        reconstructed.append(
            next_k
        )

        print(
            f"  e={e}: "
            f"k_e={current_k} "
            f"delta={delta} "
            f"t={t} "
            f"order={order} "
            f"increment={increment} "
            f"k_(e+1)={next_k} "
            f"alignment={alignment}"
        )

        current_k = next_k

    # ------------------------------------------------------------------
    # 4. RECONSTRUCTED EXPONENT SEQUENCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. RECONSTRUCTED EXPONENT SEQUENCE")
    print("=" * 78)

    print(
        f"  reconstructed={reconstructed}"
    )

    # ------------------------------------------------------------------
    # 5. SOURCE-RATIO CHECK AT ALL LEVELS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. SOURCE-RATIO CHECK")
    print("=" * 78)

    alignment_flags = []

    for e, k in enumerate(
        reconstructed,
        start=1,
    ):

        modulus = P ** e

        orbit = orbit_mod(
            k,
            modulus,
        )

        rho = ratio_mod(
            modulus
        )

        exact = (
            orbit == rho
        )

        alignment_flags.append(
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
    # 6. UNIQUE HENSEL DIGITS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. UNIQUE DIGIT CHECK")
    print("=" * 78)

    unique_flags = []

    for row in transition_rows:

        e = row["e"]
        current_k = row["current_k"]
        order = row["order"]
        chosen_t = row["t"]

        modulus_next = P ** (
            e + 1
        )

        matches = []

        for t in range(7):

            candidate_k = (
                current_k
                + t * order
            )

            if (
                orbit_mod(
                    candidate_k,
                    modulus_next,
                )
                ==
                ratio_mod(
                    modulus_next
                )
            ):
                matches.append(t)

        unique = (
            len(matches) == 1
            and matches[0] == chosen_t
        )

        unique_flags.append(
            unique
        )

        print(
            f"  e={e}: "
            f"matches={matches} "
            f"chosen={chosen_t} "
            f"unique={unique}"
        )

    # ------------------------------------------------------------------
    # 7. DEFECT / DIGIT TABLE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. SOURCE DEFECT / LIFT DIGITS")
    print("=" * 78)

    defects = [
        row["delta"]
        for row in transition_rows
    ]

    digits = [
        row["t"]
        for row in transition_rows
    ]

    print(
        f"  defects={defects}"
    )

    print(
        f"  lift_digits={digits}"
    )

    print(
        f"  t_plus_delta_mod7="
        f"{[(t + d) % P for t, d in zip(digits, defects)]}"
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
This experiment starts only from the source ratio

    rho = q1/q3.

No previously computed k_8, k_9, ..., k_12 is supplied.

The reconstruction proceeds in two stages.

First, solve the order-six orbit equation modulo 7:

    2*3^k = rho (mod 7).

This gives the unique base class k_1.

Then repeatedly:

    compute the source defect digit delta_e;

    set
        t_e = -delta_e mod 7;

    lift
        k_(e+1)
          = k_e + t_e*ord_{7^e}(3).

If every resulting exponent independently reproduces rho at the next
7-adic precision, then the entire finite exponent sequence is recovered
from rho alone.

This is stronger than Experiment 207 because the starting exponent is
also reconstructed rather than supplied.
"""
    )

    # ------------------------------------------------------------------
    # 9. FINAL
    # ------------------------------------------------------------------

    expected = [
        1,
        1,
        85,
        1555,
        73585,
        577795,
        3401371,
        3401371,
        72578983,
        798943909,
        10968052873,
    ]

    # The expected sequence above is the known sequence through e=12.
    # The new reconstruction starts from k_1 and should reproduce it.
    #
    # Note: if END_E changes, this explicit expectation should be updated.
    expected_used = expected[
        :len(reconstructed)
    ]

    exponent_exact = (
        reconstructed
        == expected_used
    )

    all_alignment = all(
        alignment_flags
    )

    all_unique = all(
        unique_flags
    )

    final_ok = (
        valuation_7(Q1) == 0
        and valuation_7(Q3) == 0
        and exponent_exact
        and all_alignment
        and all_unique
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  reconstructed={reconstructed}"
    )

    print(
        f"  expected={expected_used}"
    )

    print(
        f"  exponent_reconstruction_exact="
        f"{exponent_exact}"
    )

    print(
        f"  all_source_ratio_alignments_exact="
        f"{all_alignment}"
    )

    print(
        f"  unique_hensel_digit_every_level="
        f"{all_unique}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 208 COMPLETE")


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

