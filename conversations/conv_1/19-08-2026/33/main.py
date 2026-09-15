#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 201 — EXACT 7-ADIC POST-OVERLIFT / NEXT-DIGIT RECURRENCE AUDIT
==============================================================================

Experiment 200 established:

    k4 = 1555
    v7(F(k4)) = 5

for

    F(k) = 2*3^k*q3 - q1.

Thus k4=1555 solves the orbit equation through modulus

    7^5,

even though it was first obtained as the e=4 solution.

Experiment 201 continues from this over-lifted exponent.

At level e=5, the seven possible exponent lifts are

    k4 + t * ord_{7^5}(3),
    t = 0,...,6.

The experiment computes the normalized residual

    F(k4 + t L5) / 7^5  mod 7

and determines:

    * the next unique lift digit t4;
    * whether the residual law remains affine;
    * whether its slope remains 1;
    * the new exponent k5;
    * the resulting valuation v7(F(k5));
    * whether another over-lift occurs.

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

# Known through Experiment 200.
K1 = 1
K2 = 1
K3 = 85
K4 = 1555


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


def orbit_value_mod(k: int, modulus: int) -> int:
    return (
        BASE
        * pow(MULT, k, modulus)
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


def ratio_mod(a: int, b: int, modulus: int) -> int:
    if modulus == 1:
        return 0

    return (
        a
        * inverse_mod(b, modulus)
    ) % modulus


def multiplicative_order(a: int, modulus: int) -> int:
    if modulus == 1:
        return 1

    a %= modulus

    if a == 0:
        raise ArithmeticError(
            f"{a} is not a unit modulo {modulus}."
        )

    x = a

    for n in range(1, 10_000_000):

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
        "EXPERIMENT 201 — EXACT 7-ADIC POST-OVERLIFT / "
        "NEXT-DIGIT RECURRENCE AUDIT"
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
    # 2. KNOWN OVERLIFTED EXPONENT
    # ------------------------------------------------------------------

    F4 = (
        orbit_value(K4)
        * Q3
        - Q1
    )

    print()
    print("=" * 78)
    print("2. KNOWN OVERLIFTED LEVEL")
    print("=" * 78)

    print(
        f"  k4={K4}"
    )

    print(
        f"  v7(F(k4))={v7(F4)}"
    )

    print(
        f"  F(k4) divisible by 7^5="
        f"{F4 % (P ** 5) == 0}"
    )

    print(
        f"  F(k4) divisible by 7^6="
        f"{F4 % (P ** 6) == 0}"
    )

    # ------------------------------------------------------------------
    # 3. NEXT LEVEL DATA
    # ------------------------------------------------------------------

    e = 5

    modulus = P ** e

    next_modulus = P ** (e + 1)

    order = multiplicative_order(
        MULT,
        modulus,
    )

    rho = ratio_mod(
        Q1,
        Q3,
        next_modulus,
    )

    print()
    print("=" * 78)
    print("3. POST-OVERLIFT LEVEL")
    print("=" * 78)

    print(
        f"  e={e}"
    )

    print(
        f"  modulus={modulus}"
    )

    print(
        f"  next_modulus={next_modulus}"
    )

    print(
        f"  ord_{modulus}(3)={order}"
    )

    print(
        f"  rho_mod_{next_modulus}={rho}"
    )

    # ------------------------------------------------------------------
    # 4. SEVEN CANDIDATE LIFTS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. COMPLETE SEVEN-CANDIDATE POST-OVERLIFT AUDIT")
    print("=" * 78)

    candidates = []

    for t in range(7):

        k = K4 + t * order

        F = (
            orbit_value(k)
            * Q3
            - Q1
        )

        divisible = (
            F % (P ** e) == 0
        )

        if not divisible:
            raise ArithmeticError(
                f"Candidate k={k} is not divisible by 7^{e}."
            )

        normalized = (
            (F // (P ** e))
            % P
        )

        target_match = (
            orbit_value_mod(
                k,
                next_modulus,
            )
            == rho
        )

        candidates.append(
            {
                "t": t,
                "k": k,
                "F": F,
                "v7": v7(F),
                "normalized": normalized,
                "target_match": target_match,
            }
        )

        print(
            f"  t={t}: "
            f"k={k} "
            f"normalized={normalized} "
            f"v7(F)={v7(F)} "
            f"target_match={target_match}"
        )

    successful = [
        row
        for row in candidates
        if row["target_match"]
    ]

    if len(successful) != 1:
        raise ArithmeticError(
            f"Expected one successful post-overlift lift, "
            f"got {successful}."
        )

    chosen = successful[0]

    # ------------------------------------------------------------------
    # 5. AFFINE RESIDUAL LAW
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. POST-OVERLIFT AFFINE RESIDUAL LAW")
    print("=" * 78)

    A = candidates[0]["normalized"]

    B = (
        candidates[1]["normalized"]
        - candidates[0]["normalized"]
    ) % P

    affine_checks = []

    for row in candidates:

        predicted = (
            A
            + row["t"] * B
        ) % P

        exact = (
            predicted
            == row["normalized"]
        )

        affine_checks.append(
            exact
        )

        print(
            f"  t={row['t']}: "
            f"actual={row['normalized']} "
            f"predicted={predicted} "
            f"exact={exact}"
        )

    affine_exact = all(
        affine_checks
    )

    print(
        f"  A={A}"
    )

    print(
        f"  B={B}"
    )

    # ------------------------------------------------------------------
    # 6. NEXT EXPONENT
    # ------------------------------------------------------------------

    k5 = chosen["k"]

    print()
    print("=" * 78)
    print("6. NEXT EXPONENT")
    print("=" * 78)

    print(
        f"  chosen_t={chosen['t']}"
    )

    print(
        f"  k5={k5}"
    )

    print(
        f"  k5 = k4 + t*ord_{modulus}(3)="
        f"{K4} + {chosen['t']}*{order}"
    )

    print(
        f"  recurrence_exact="
        f"{k5 == K4 + chosen['t'] * order}"
    )

    # ------------------------------------------------------------------
    # 7. NEW RESIDUAL VALUATION
    # ------------------------------------------------------------------

    F5 = chosen["F"]
    v5 = v7(F5)

    print()
    print("=" * 78)
    print("7. NEW RESIDUAL VALUATION")
    print("=" * 78)

    print(
        f"  v7(F(k5))={v5}"
    )

    print(
        f"  required>=6="
        f"{v5 is not None and v5 >= 6}"
    )

    print(
        f"  solves_mod_7^6="
        f"{F5 % (P ** 6) == 0}"
    )

    print(
        f"  solves_mod_7^7="
        f"{F5 % (P ** 7) == 0}"
    )

    # ------------------------------------------------------------------
    # 8. NEW OVERLIFT TEST
    # ------------------------------------------------------------------

    next_digit = (
        (F5 // (P ** 6))
        % P
        if F5 % (P ** 6) == 0
        else None
    )

    print()
    print("=" * 78)
    print("8. NEW OVERLIFT TEST")
    print("=" * 78)

    print(
        f"  F(k5)/7^6 mod7={next_digit}"
    )

    if next_digit == 0:
        print(
            "  new_overlift=True"
        )
    else:
        print(
            "  new_overlift=False"
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
Experiment 200 found a single over-lift:

    k4=1555,
    v7(F(k4))=5.

Experiment 201 asks what happens after that exceptional extra digit.

The seven candidates

    k4 + t*ord_{7^5}(3)

form the complete next Hensel lift family.

Their normalized residuals determine the unique new exponent digit.

This separates two questions:

    over-lift:
        why k4 survives one extra power of 7;

    next lift:
        which exponent correction is required after that extra power.

The experiment also checks whether the post-overlift residual law
remains affine with slope 1, as observed at the previous levels.
"""
    )

    # ------------------------------------------------------------------
    # 10. FINAL
    # ------------------------------------------------------------------

    unique = (
        len(successful) == 1
    )

    recurrence_exact = (
        k5
        == K4 + chosen["t"] * order
    )

    expected_no_assumption = (
        unique
        and affine_exact
        and recurrence_exact
        and F5 % (P ** 6) == 0
        and chosen["normalized"] == 0
    )

    final_ok = (
        v7(Q1) == 0
        and v7(Q3) == 0
        and v7(F4) == 5
        and unique
        and affine_exact
        and recurrence_exact
        and expected_no_assumption
    )

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  previous_overlift_v7_is_5="
        f"{v7(F4) == 5}"
    )

    print(
        f"  unique_next_lift="
        f"{unique}"
    )

    print(
        f"  post_overlift_affine_exact="
        f"{affine_exact}"
    )

    print(
        f"  next_exponent_recurrence_exact="
        f"{recurrence_exact}"
    )

    print(
        f"  chosen_lift_solves_mod_7^6="
        f"{F5 % (P ** 6) == 0}"
    )

    print(
        f"  chosen_normalized_digit_zero="
        f"{chosen['normalized'] == 0}"
    )

    print(
        f"  new_overlift="
        f"{next_digit == 0}"
    )

    print(
        f"  chosen_t={chosen['t']}"
    )

    print(
        f"  k5={k5}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 201 COMPLETE")


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

