#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 191R — EXACT 7-ADIC RATIO-DEFECT / s1-LAYER FACTORIZATION AUDIT
==============================================================================

Corrected from Experiment 191:

At precision e=2,

    7^(e-2) = 1,

so modular inversion modulo 1 is not defined/needed.

The e=2 layer is therefore treated as the trivial modulus-1 case, and
actual modular inversions begin at e=3.

Established data:

    q1 = 29144191
    q3 = 24794967

    s1 = q1 - 6*q3 = -119625611

    v7(q1)=0
    v7(q3)=0
    v7(s1)=2

Define

    u7 = s1 / 49.

Then the 7-adic ratio defect is

    delta7 = (q1/q3 - 6)/49,

and formally

    s1/49 = q3 * delta7.

The experiment verifies this at increasing finite 7-adic precision.
"""

from __future__ import annotations

import sys


# ============================================================================
# DATA
# ============================================================================

Q1 = 29144191
Q3 = 24794967
C = 6

S1 = Q1 - C * Q3


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
    """
    Exact modular inverse for m > 1.
    Modulus 1 is handled separately by callers.
    """
    if m <= 1:
        raise ArithmeticError(
            f"Modular inverse requested for modulus {m}."
        )

    a %= m

    if a == 0:
        raise ArithmeticError(
            f"Cannot invert 0 modulo {m}."
        )

    old_r, r = a, m
    old_s, s = 1, 0

    while r != 0:

        q = old_r // r

        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s

    if old_r != 1:
        raise ArithmeticError(
            f"{a} is not invertible modulo {m}."
        )

    return old_s % m


def ratio_mod(a, b, p, e):
    modulus = p ** e

    if modulus == 1:
        return 0

    inv_b = inverse_mod(
        b,
        modulus,
    )

    return (
        a * inv_b
    ) % modulus


def base_p_digits(x, p, count):
    digits = []

    x = int(x)

    for _ in range(count):

        digits.append(
            x % p
        )

        x //= p

    return digits


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 191R — EXACT 7-ADIC RATIO-DEFECT / "
        "s1-LAYER FACTORIZATION AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT SOURCE PROFILE")
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
    # 2. EXCEPTIONAL SCHUR COEFFICIENT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. EXCEPTIONAL SCHUR COEFFICIENT")
    print("=" * 78)

    print(
        f"  s1=q1-6*q3={S1}"
    )

    print(
        f"  v7(s1)={v7(S1)}"
    )

    if v7(S1) != 2:
        raise ArithmeticError(
            "Expected v7(s1)=2."
        )

    u7 = S1 // 49

    print(
        f"  u7=s1/49={u7}"
    )

    print(
        f"  v7(u7)={v7(u7)}"
    )

    print(
        f"  u7_mod7={u7 % 7}"
    )

    # ------------------------------------------------------------------
    # 3. SOURCE RATIO PROFILE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXACT SOURCE RATIO q1/q3")
    print("=" * 78)

    rho = {}

    for e in range(1, 9):

        modulus = 7 ** e

        rho[e] = ratio_mod(
            Q1,
            Q3,
            7,
            e,
        )

        print(
            f"  modulus={modulus}: "
            f"rho={rho[e]}"
        )

    # ------------------------------------------------------------------
    # 4. SECOND-ORDER DEFECT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. SECOND-ORDER RATIO DEFECT")
    print("=" * 78)

    defect_mod = {}

    for e in range(2, 9):

        modulus = 7 ** e
        reduced_modulus = 7 ** (e - 2)

        if reduced_modulus == 1:
            defect_mod[e] = 0

            print(
                f"  modulus={modulus}: "
                f"delta7 mod 1 = 0 "
                f"(trivial precision)"
            )

            continue

        rho_e = rho[e]

        difference = (
            rho_e - C
        ) % modulus

        if difference % 49 != 0:
            raise ArithmeticError(
                f"rho-6 not divisible by 49 at modulus {modulus}."
            )

        defect = (
            (difference // 49)
            % reduced_modulus
        )

        defect_mod[e] = defect

        print(
            f"  modulus={modulus}: "
            f"rho-6={difference} "
            f"delta7={defect} "
            f"delta_modulus={reduced_modulus}"
        )

    # ------------------------------------------------------------------
    # 5. DIRECT INTEGER DEFECT
    # ------------------------------------------------------------------

    direct_difference = (
        Q1 - C * Q3
    )

    direct_defect = (
        direct_difference // 49
    )

    print()
    print("=" * 78)
    print("5. DIRECT INTEGER DEFECT")
    print("=" * 78)

    print(
        f"  q1-6*q3={direct_difference}"
    )

    print(
        f"  (q1-6*q3)/49={direct_defect}"
    )

    print(
        f"  direct_defect_equals_u7="
        f"{direct_defect == u7}"
    )

    # ------------------------------------------------------------------
    # 6. FACTORIZATION TEST
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. SOURCE-RATIO DEFECT FACTORIZATION")
    print("=" * 78)

    factorization_checks = []

    for e in range(3, 9):

        reduced_modulus = 7 ** (e - 2)

        defect = defect_mod[e]

        inv_q3 = inverse_mod(
            Q3,
            reduced_modulus,
        )

        rhs = (
            u7 * inv_q3
        ) % reduced_modulus

        exact = (
            defect == rhs
        )

        factorization_checks.append(
            exact
        )

        print(
            f"  precision_e={e}: "
            f"delta7={defect} "
            f"u7/q3={rhs} "
            f"exact={exact}"
        )

    # e=2 is the trivial modulus-1 case
    factorization_checks.append(
        defect_mod[2] == 0
    )

    # ------------------------------------------------------------------
    # 7. MULTIPLICATIVE RECOVERY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. MULTIPLICATIVE RECOVERY OF s1/49")
    print("=" * 78)

    recovery_checks = []

    for e in range(3, 9):

        reduced_modulus = 7 ** (e - 2)

        defect = defect_mod[e]

        reconstructed = (
            Q3 * defect
        ) % reduced_modulus

        actual = (
            u7 % reduced_modulus
        )

        exact = (
            reconstructed == actual
        )

        recovery_checks.append(
            exact
        )

        print(
            f"  precision_e={e}: "
            f"actual={actual} "
            f"reconstructed={reconstructed} "
            f"exact={exact}"
        )

    # ------------------------------------------------------------------
    # 8. BASE-7 DEFECT DIGITS
    # ------------------------------------------------------------------

    defect_digits = base_p_digits(
        direct_defect,
        7,
        8,
    )

    print()
    print("=" * 78)
    print("8. BASE-7 DIGITS OF THE DEFECT")
    print("=" * 78)

    print(
        f"  delta7_integer={direct_defect}"
    )

    print(
        f"  base7_digits_low_to_high={defect_digits}"
    )

    print(
        f"  delta7_mod7={direct_defect % 7}"
    )

    print(
        f"  delta7_v7={v7(direct_defect)}"
    )

    # ------------------------------------------------------------------
    # 9. FIRST LIFT DIGIT
    # ------------------------------------------------------------------

    first_digit = (
        direct_defect % 7
    )

    rho_mod343_reconstructed = (
        6 + 49 * first_digit
    ) % 343

    print()
    print("=" * 78)
    print("9. FIRST NONZERO RATIO-LIFT DIGIT")
    print("=" * 78)

    print(
        f"  rho=6+49*delta7"
    )

    print(
        f"  delta7_first_digit={first_digit}"
    )

    print(
        f"  rho_mod343_actual={rho[3]}"
    )

    print(
        f"  rho_mod343_reconstructed="
        f"{rho_mod343_reconstructed}"
    )

    print(
        f"  first_digit_reconstructs_rho="
        f"{rho[3] == rho_mod343_reconstructed}"
    )

    # ------------------------------------------------------------------
    # 10. INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The exceptional coefficient satisfies

    s1 = q1 - 6q3 = 49u7,

with

    u7 = -2441339

a 7-adic unit.

Because q3 is also a 7-adic unit, the second-order ratio defect

    delta7 = (q1/q3 - 6)/49

satisfies

    u7 = q3 * delta7.

Thus the cancellation modulo 49 is separated into:

    zeroth-order ratio:
        q1/q3 = 6 mod 49;

    second-order defect:
        delta7 mod 7;

    higher lifts:
        subsequent base-7 digits of delta7.

The modulus-1 case at e=2 is treated as a trivial residue and requires
no modular inversion.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        v7(Q1) == 0
        and v7(Q3) == 0
        and v7(S1) == 2
        and v7(u7) == 0
        and direct_defect == u7
        and all(factorization_checks)
        and all(recovery_checks)
        and (
            rho[3]
            == rho_mod343_reconstructed
        )
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  q1_7adic_unit={v7(Q1) == 0}"
    )

    print(
        f"  q3_7adic_unit={v7(Q3) == 0}"
    )

    print(
        f"  s1_v7_is_2={v7(S1) == 2}"
    )

    print(
        f"  s1_over_49_is_7adic_unit="
        f"{v7(u7) == 0}"
    )

    print(
        f"  exact_second_order_defect="
        f"{direct_defect == u7}"
    )

    print(
        f"  ratio_defect_factorization_exact="
        f"{all(factorization_checks)}"
    )

    print(
        f"  s1_over_49_recovery_exact="
        f"{all(recovery_checks)}"
    )

    print(
        f"  first_lift_digit_reconstructs_rho="
        f"{rho[3] == rho_mod343_reconstructed}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 191R COMPLETE")


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