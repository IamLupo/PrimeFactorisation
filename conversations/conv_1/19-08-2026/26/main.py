#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 195 — EXACT 49-ADIC ORBIT ALIGNMENT / SOURCE-RATIO RETURN AUDIT
==============================================================================

Experiment 194 established:

    ord_7(3)   = 6
    ord_49(3)  = 42
    ord_343(3) = 294

and

    rho = q1/q3,

with

    rho == 6 (mod 49),
    rho != 6 (mod 343).

The orbit is

    c_k = 2*3^k.

Thus c_1 = 6 is the exceptional orbit point.

Experiment 195 asks a sharper question:

    Is rho exactly aligned with the orbit point c_1 modulo 49
    because 1 is its unique exponent class modulo 42?

We search the complete exponent classes

    k = 0,...,41

modulo 49 and 343.

Then we determine:

    * all k with c_k == rho mod 49;
    * all k with c_k == rho mod 343;
    * whether the mod-49 match is unique;
    * whether the mod-343 match disappears;
    * the first return k=43 to c_1 modulo 49;
    * comparison of rho-c_k valuations for all 42 classes.

This directly distinguishes:

    orbit alignment modulo 49

from

    accidental first-order agreement modulo 7.

No extrapolation.
"""

from __future__ import annotations

import sys


Q1 = 29144191
Q3 = 24794967

BASE = 2
MULT = 3

P = 7


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


def orbit_value(k):
    return BASE * (MULT ** k)


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 195 — EXACT 49-ADIC ORBIT ALIGNMENT / "
        "SOURCE-RATIO RETURN AUDIT"
    )
    print("=" * 78)

    rho7 = ratio_mod(
        Q1,
        Q3,
        7,
    )

    rho49 = ratio_mod(
        Q1,
        Q3,
        49,
    )

    rho343 = ratio_mod(
        Q1,
        Q3,
        343,
    )

    # ------------------------------------------------------------------
    # 1. SOURCE RATIO
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. SOURCE RATIO")
    print("=" * 78)

    print(f"  q1={Q1}")
    print(f"  q3={Q3}")
    print(f"  rho_mod7={rho7}")
    print(f"  rho_mod49={rho49}")
    print(f"  rho_mod343={rho343}")

    # ------------------------------------------------------------------
    # 2. COMPLETE MOD-49 ORBIT
    # ------------------------------------------------------------------

    orbit49 = []

    for k in range(42):

        c = orbit_value(k)

        orbit49.append(
            (
                k,
                c % 49,
            )
        )

    matches49 = [
        k
        for k, residue in orbit49
        if residue == rho49
    ]

    print()
    print("=" * 78)
    print("2. COMPLETE MOD-49 ORBIT")
    print("=" * 78)

    print(
        f"  rho_mod49={rho49}"
    )

    print(
        f"  matching_k_mod49={matches49}"
    )

    print(
        f"  match_count_mod49={len(matches49)}"
    )

    # ------------------------------------------------------------------
    # 3. COMPLETE MOD-343 ORBIT
    # ------------------------------------------------------------------

    orbit343 = []

    for k in range(294):

        c = orbit_value(k)

        orbit343.append(
            (
                k,
                c % 343,
            )
        )

    matches343 = [
        k
        for k, residue in orbit343
        if residue == rho343
    ]

    print()
    print("=" * 78)
    print("3. COMPLETE MOD-343 ORBIT")
    print("=" * 78)

    print(
        f"  rho_mod343={rho343}"
    )

    print(
        f"  matching_k_mod343={matches343}"
    )

    print(
        f"  match_count_mod343={len(matches343)}"
    )

    # ------------------------------------------------------------------
    # 4. FIRST MOD-49 RETURN TO c1
    # ------------------------------------------------------------------

    c1_mod49 = (
        orbit_value(1) % 49
    )

    first_return49 = None

    for k in range(2, 100):

        if (
            orbit_value(k) % 49
            == c1_mod49
        ):

            first_return49 = k
            break

    print()
    print("=" * 78)
    print("4. FIRST RETURN TO c1 MOD 49")
    print("=" * 78)

    print(
        f"  c1_mod49={c1_mod49}"
    )

    print(
        f"  first_return_k={first_return49}"
    )

    # ------------------------------------------------------------------
    # 5. COMPLETE 42-CLASS VALUATION PROFILE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. v7(rho-c_k), k=0,...,41")
    print("=" * 78)

    valuation_profile = []

    for k in range(42):

        c = orbit_value(k)

        value = (
            Q1 - c * Q3
        )

        valuation = v7(value)

        valuation_profile.append(
            (
                k,
                c,
                valuation,
            )
        )

        print(
            f"  k={k}: "
            f"c_k={c} "
            f"v7(rho-c_k)={valuation}"
        )

    # ------------------------------------------------------------------
    # 6. MAXIMUM VALUATION
    # ------------------------------------------------------------------

    finite = [
        item
        for item in valuation_profile
        if item[2] is not None
    ]

    maximum_v = max(
        item[2]
        for item in finite
    )

    maximizers = [
        item[0]
        for item in finite
        if item[2] == maximum_v
    ]

    print()
    print("=" * 78)
    print("6. COMPLETE ORBIT MAXIMUM")
    print("=" * 78)

    print(
        f"  maximum_v7={maximum_v}"
    )

    print(
        f"  maximizer_k={maximizers}"
    )

    # ------------------------------------------------------------------
    # 7. SOURCE-RATIO ALIGNMENT
    # ------------------------------------------------------------------

    mod49_alignment = (
        rho49
        == orbit_value(1) % 49
    )

    mod343_alignment = (
        rho343
        == orbit_value(1) % 343
    )

    print()
    print("=" * 78)
    print("7. c1 SOURCE-RATIO ALIGNMENT")
    print("=" * 78)

    print(
        f"  c1_mod49={orbit_value(1) % 49}"
    )

    print(
        f"  rho_mod49={rho49}"
    )

    print(
        f"  exact_alignment_mod49="
        f"{mod49_alignment}"
    )

    print(
        f"  c1_mod343={orbit_value(1) % 343}"
    )

    print(
        f"  rho_mod343={rho343}"
    )

    print(
        f"  exact_alignment_mod343="
        f"{mod343_alignment}"
    )

    # ------------------------------------------------------------------
    # 8. NEXT RETURN k=43
    # ------------------------------------------------------------------

    k43 = 43

    c43 = orbit_value(
        k43
    )

    print()
    print("=" * 78)
    print("8. k=43 RETURN")
    print("=" * 78)

    print(
        f"  c43_mod49={c43 % 49}"
    )

    print(
        f"  c1_mod49={orbit_value(1) % 49}"
    )

    print(
        f"  same_mod49="
        f"{c43 % 49 == orbit_value(1) % 49}"
    )

    print(
        f"  c43_mod343={c43 % 343}"
    )

    print(
        f"  c1_mod343={orbit_value(1) % 343}"
    )

    print(
        f"  same_mod343="
        f"{c43 % 343 == orbit_value(1) % 343}"
    )

    print(
        f"  v7(c43-c1)="
        f"{v7(c43 - orbit_value(1))}"
    )

    print(
        f"  v7(rho-c43)="
        f"{v7(Q1 - c43 * Q3)}"
    )

    # ------------------------------------------------------------------
    # 9. INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The complete mod-49 orbit has 42 distinct exponents because

    ord_49(3)=42.

Therefore a match between rho and c_k modulo 49 is an exact
orbit-class statement.

If rho matches only c_1 modulo 49, then the exceptional middle
Schur coefficient is attached to a unique orbit point at this
precision.

At modulus 343 the match should disappear, because Experiment 190
established

    rho != 6 (mod 343).

The return at k=43 is different: it returns the orbit itself to c1
modulo 49, but not modulo 343.

Thus k=1 is the unique source-ratio alignment in the 42-class
49-adic orbit, even though the orbit itself returns to the same
residue after 42 steps.
"""
    )

    # ------------------------------------------------------------------
    # 10. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        rho49 == 6
        and rho343 == 104
        and matches49 == [1]
        and matches343 == []
        and first_return49 == 43
        and maximum_v == 2
        and 1 in maximizers
        and c43 % 49 == orbit_value(1) % 49
        and c43 % 343 != orbit_value(1) % 343
        and mod49_alignment
        and not mod343_alignment
        and v7(Q1 - orbit_value(1) * Q3) == 2
        and v7(Q1 - c43 * Q3) == 2
    )

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  rho_matches_unique_k1_mod49="
        f"{matches49 == [1]}"
    )

    print(
        f"  rho_matches_no_orbit_point_mod343="
        f"{matches343 == []}"
    )

    print(
        f"  first_return_to_c1_mod49_is_43="
        f"{first_return49 == 43}"
    )

    print(
        f"  maximum_orbit_v7_is_2="
        f"{maximum_v == 2}"
    )

    print(
        f"  k1_is_maximizer="
        f"{1 in maximizers}"
    )

    print(
        f"  k43_returns_mod49="
        f"{c43 % 49 == orbit_value(1) % 49}"
    )

    print(
        f"  k43_differs_mod343="
        f"{c43 % 343 != orbit_value(1) % 343}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 195 COMPLETE")


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

