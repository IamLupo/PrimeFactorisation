#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 190 — EXACT 7-ADIC CANCELLATION / SOURCE-RATIO LIFT AUDIT
==============================================================================

Experiment 189 found:

    v7(q1)=0
    v7(q3)=0

but

    v7(s1)=2,

where

    s1 = q1 - 6 q3.

Therefore the middle terminal Schur coefficient is divisible by 49
because of exact 7-adic cancellation.

Since q3 is a 7-adic unit, we can invert q3 modulo powers of 7 and
form the exact source ratio

    rho = q1 / q3  (mod 7^e).

The cancellation condition becomes

    s_k = q1 - c_k q3,

where

    c0 = 2
    c1 = 6
    c2 = 18,

and therefore

    v7(s_k) >= e
        iff
    rho == c_k (mod 7^e).

Experiment 190 determines exactly:

    * rho modulo 7, 49, 343, ...
    * comparison with 2, 6, 18;
    * the exact cancellation depth of each s_k;
    * the normalized unit s1 / 49;
    * the first modulus where rho separates from 6;
    * whether the exceptional s1 divisibility is fully explained by
      the lifted source ratio.

No floating point.
No SymPy.
No extrapolation.
"""

from __future__ import annotations

import sys


# ============================================================================
# PRIMITIVE SOURCE DATA
# ============================================================================

Q1 = 29144191
Q3 = 24794967

COEFFICIENTS = {
    "s0": 2,
    "s1": 6,
    "s2": 18,
}

S = {
    "s0": Q1 - 2 * Q3,
    "s1": Q1 - 6 * Q3,
    "s2": Q1 - 18 * Q3,
}


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
    Exact modular inverse using extended Euclid.
    """
    a %= m

    if a == 0:
        raise ArithmeticError(
            f"No inverse for 0 modulo {m}."
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


def ratio_mod_prime_power(a, b, p, e):
    modulus = p ** e

    inv_b = inverse_mod(
        b,
        modulus,
    )

    return (
        a * inv_b
    ) % modulus


def first_separation(a, b, p, max_e=12):
    for e in range(1, max_e + 1):

        modulus = p ** e

        if a % modulus != b % modulus:
            return e

    return None


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 190 — EXACT 7-ADIC CANCELLATION / "
        "SOURCE-RATIO LIFT AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE PROFILE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. SOURCE 7-ADIC PROFILE")
    print("=" * 78)

    print(
        f"  q1={Q1} "
        f"v7={v7(Q1)} "
        f"mod7={Q1 % 7}"
    )

    print(
        f"  q3={Q3} "
        f"v7={v7(Q3)} "
        f"mod7={Q3 % 7}"
    )

    if v7(Q3) != 0:
        raise ArithmeticError(
            "q3 must be a 7-adic unit."
        )

    # ------------------------------------------------------------------
    # 2. SCHUR COEFFICIENTS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. EXACT SCHUR COEFFICIENTS")
    print("=" * 78)

    for name in ("s0", "s1", "s2"):

        print(
            f"  {name}={S[name]} "
            f"v7={v7(S[name])} "
            f"mod7={S[name] % 7}"
        )

    # ------------------------------------------------------------------
    # 3. SOURCE RATIO rho = q1/q3
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXACT SOURCE RATIO q1/q3 MOD 7^e")
    print("=" * 78)

    ratio_profiles = {}

    for e in range(1, 8):

        modulus = 7 ** e

        rho = ratio_mod_prime_power(
            Q1,
            Q3,
            7,
            e,
        )

        ratio_profiles[e] = rho

        print(
            f"  modulus={modulus}: "
            f"rho={rho}"
        )

    # ------------------------------------------------------------------
    # 4. COMPARE rho WITH SCHUR COEFFICIENTS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. RATIO VS SCHUR CANCELLATION CONSTANTS")
    print("=" * 78)

    for name, coeff in COEFFICIENTS.items():

        print()
        print(
            f"  {name}: coefficient={coeff}"
        )

        for e in range(1, 7):

            modulus = 7 ** e

            rho = ratio_profiles[e]

            same = (
                rho == (coeff % modulus)
            )

            print(
                f"    mod={modulus}: "
                f"rho=={coeff} -> {same}"
            )

    # ------------------------------------------------------------------
    # 5. CANCELLATION DEPTH
    # ------------------------------------------------------------------

    cancellation_depths = {}

    for name, coeff in COEFFICIENTS.items():

        value = S[name]

        cancellation_depths[name] = (
            v7(value)
        )

    print()
    print("=" * 78)
    print("5. EXACT CANCELLATION DEPTH")
    print("=" * 78)

    for name in ("s0", "s1", "s2"):

        print(
            f"  {name}: "
            f"v7={cancellation_depths[name]}"
        )

    # ------------------------------------------------------------------
    # 6. FIRST RATIO SEPARATION FROM EACH COEFFICIENT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. FIRST SOURCE-RATIO SEPARATION")
    print("=" * 78)

    for name, coeff in COEFFICIENTS.items():

        sep = first_separation(
            Q1 * inverse_mod(Q3, 7),
            coeff,
            7,
            max_e=1,
        )

        # The direct expression above only concerns mod 7.
        # For the actual lift, determine the first e for which
        # q1/q3 != coeff modulo 7^e.
        sep_full = None

        for e in range(1, 10):

            modulus = 7 ** e

            rho = ratio_mod_prime_power(
                Q1,
                Q3,
                7,
                e,
            )

            if rho != coeff % modulus:

                sep_full = e
                break

        print(
            f"  {name}: "
            f"first_rho_separation_e={sep_full} "
            f"modulus={7 ** sep_full if sep_full else None}"
        )

    # ------------------------------------------------------------------
    # 7. s1 / 49 NORMALIZED UNIT
    # ------------------------------------------------------------------

    s1 = S["s1"]
    s1_v7 = v7(s1)

    if s1_v7 is None or s1_v7 < 2:
        raise ArithmeticError(
            "Expected s1 to have v7 >= 2."
        )

    s1_layer = s1 // (7 ** 2)

    print()
    print("=" * 78)
    print("7. s1 NORMALIZED 7-ADIC UNIT")
    print("=" * 78)

    print(
        f"  s1={s1}"
    )

    print(
        f"  s1/49={s1_layer}"
    )

    print(
        f"  v7(s1/49)={v7(s1_layer)}"
    )

    print(
        f"  (s1/49) mod7={s1_layer % 7}"
    )

    # ------------------------------------------------------------------
    # 8. DIRECT CONGRUENCE q1 - 6q3
    # ------------------------------------------------------------------

    c = 6

    print()
    print("=" * 78)
    print("8. DIRECT 7-ADIC CANCELLATION TEST FOR s1")
    print("=" * 78)

    for e in range(1, 7):

        modulus = 7 ** e

        lhs = Q1 % modulus
        rhs = (
            c * Q3
        ) % modulus

        equal = lhs == rhs

        print(
            f"  modulus={modulus}: "
            f"q1={lhs} "
            f"6*q3={rhs} "
            f"equal={equal}"
        )

    # ------------------------------------------------------------------
    # 9. CANCELLATION EXPLANATION
    # ------------------------------------------------------------------

    cancellation_explained = (
        v7(S["s1"]) == 2
        and ratio_profiles[1] == 6 % 7
        and ratio_profiles[2] == 6 % 49
        and ratio_profiles[3] != 6 % 343
    )

    print()
    print("=" * 78)
    print("9. CANCELLATION EXPLANATION")
    print("=" * 78)

    print(
        f"  rho_equals_6_mod7="
        f"{ratio_profiles[1] == 6 % 7}"
    )

    print(
        f"  rho_equals_6_mod49="
        f"{ratio_profiles[2] == 6 % 49}"
    )

    print(
        f"  rho_differs_from_6_mod343="
        f"{ratio_profiles[3] != 6 % 343}"
    )

    print(
        f"  exact_v7_s1_is_2="
        f"{v7(S['s1']) == 2}"
    )

    print(
        f"  cancellation_fully_explained="
        f"{cancellation_explained}"
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
Experiment 189 found the exceptional valuation

    v7(s1)=2,

while q1 and q3 are both 7-adic units.

The exact source-ratio criterion is

    s1 = q1 - 6 q3,

so

    v7(s1) >= e

exactly when

    q1/q3 == 6 (mod 7^e).

Experiment 190 therefore determines the cancellation depth directly
in the 7-adic source ratio.

If

    q1/q3 == 6 (mod 49)

but

    q1/q3 != 6 (mod 343),

then the valuation

    v7(s1)=2

is completely explained by a second-order 7-adic coincidence.

The normalized quotient

    s1/49

then gives the first nonzero 7-adic obstruction beyond the
cancellation layer.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        v7(Q1) == 0
        and v7(Q3) == 0
        and v7(S["s1"]) == 2
        and ratio_profiles[1] == 6 % 7
        and ratio_profiles[2] == 6 % 49
        and ratio_profiles[3] != 6 % 343
        and cancellation_explained
        and v7(s1_layer) == 0
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
        f"  s1_v7_is_2={v7(S['s1']) == 2}"
    )

    print(
        f"  rho_equals_6_mod7="
        f"{ratio_profiles[1] == 6 % 7}"
    )

    print(
        f"  rho_equals_6_mod49="
        f"{ratio_profiles[2] == 6 % 49}"
    )

    print(
        f"  rho_differs_6_mod343="
        f"{ratio_profiles[3] != 6 % 343}"
    )

    print(
        f"  s1_over_49_is_7adic_unit="
        f"{v7(s1_layer) == 0}"
    )

    print(
        f"  cancellation_fully_explained="
        f"{cancellation_explained}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 190 COMPLETE")


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

