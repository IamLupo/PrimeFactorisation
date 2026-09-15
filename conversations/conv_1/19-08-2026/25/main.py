#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 194 — EXACT PERIOD-6 ORBIT / 49-ADIC RETURN-DEPTH AUDIT
==============================================================================

Experiment 193 established:

    ord_7(3)  = 6
    ord_49(3) = 42

and for the orbit

    c_k = 2*3^k

the source ratio rho=q1/q3 satisfies:

    v7(rho-c_1)=2,

while

    v7(rho-c_7)=1.

Since

    7 == 1 (mod 6),

the indices k=1 and k=7 return to the same residue modulo 7.

However, because

    ord_49(3)=42,

the return modulo 49 occurs only after 42 steps.

Experiment 194 studies this hierarchy explicitly.

For indices

    k = 1 + 6m,

we compare

    c_k = 2*3^k

to c_1

modulo 7, 49, and 343.

The exact questions are:

    * when does the orbit return to the same mod-7 class?
    * when does it return to the same mod-49 class?
    * what is the 7-adic valuation of c_k-c_1?
    * how does that compare with rho-c_k?
    * is the k=7 valuation-one event explained entirely by orbit
      recurrence, rather than by a second independent coincidence?

The finite range is deliberately chosen to include a complete
mod-49 orbit period.

No floating point.
No SymPy.
No extrapolation.
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

K_MAX = 42


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


def ratio_mod(a, b, modulus):
    if modulus == 1:
        return 0

    return (
        a * inverse_mod(b, modulus)
    ) % modulus


def orbit_value(k):
    return (
        BASE
        * (MULT ** k)
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 194 — EXACT PERIOD-6 ORBIT / "
        "49-ADIC RETURN-DEPTH AUDIT"
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
        f"  rho_mod7={rho7}"
    )

    print(
        f"  rho_mod49={rho49}"
    )

    print(
        f"  rho_mod343={rho343}"
    )

    # ------------------------------------------------------------------
    # 2. ORDERS
    # ------------------------------------------------------------------

    def multiplicative_order(a, m):

        x = a % m

        for n in range(1, 10000):

            if x == 1:
                return n

            x = (
                x * a
            ) % m

        raise ArithmeticError(
            f"Could not determine order of {a} mod {m}."
        )

    ord7 = multiplicative_order(
        3,
        7,
    )

    ord49 = multiplicative_order(
        3,
        49,
    )

    ord343 = multiplicative_order(
        3,
        343,
    )

    print()
    print("=" * 78)
    print("2. MULTIPLICATIVE ORDERS")
    print("=" * 78)

    print(
        f"  ord_7(3)={ord7}"
    )

    print(
        f"  ord_49(3)={ord49}"
    )

    print(
        f"  ord_343(3)={ord343}"
    )

    # ------------------------------------------------------------------
    # 3. REFERENCE k=1
    # ------------------------------------------------------------------

    c1 = orbit_value(1)

    print()
    print("=" * 78)
    print("3. REFERENCE ORBIT ELEMENT")
    print("=" * 78)

    print(
        f"  k=1"
    )

    print(
        f"  c1={c1}"
    )

    print(
        f"  rho-c1="
        f"{Q1 - c1 * Q3}"
    )

    print(
        f"  v7(rho-c1)="
        f"{v7(Q1 - c1 * Q3)}"
    )

    # ------------------------------------------------------------------
    # 4. PERIOD-6 RETURNS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. k = 1 + 6m RETURN SEQUENCE")
    print("=" * 78)

    returns = []

    for m in range(0, 8):

        k = 1 + 6 * m

        c = orbit_value(k)

        delta_orbit = c - c1

        rho_delta = (
            Q1 - c * Q3
        )

        v_orbit = v7(
            delta_orbit
        )

        v_rho = v7(
            rho_delta
        )

        returns.append(
            (
                m,
                k,
                c,
                delta_orbit,
                v_orbit,
                rho_delta,
                v_rho,
            )
        )

        print(
            f"  m={m} k={k}: "
            f"c_k={c} "
            f"v7(c_k-c1)={v_orbit} "
            f"v7(rho-c_k)={v_rho}"
        )

    # ------------------------------------------------------------------
    # 5. RETURN CONGRUENCE TABLE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. RETURN CONGRUENCES")
    print("=" * 78)

    for m, k, c, delta_orbit, v_orbit, rho_delta, v_rho in returns:

        print(
            f"  k={k}: "
            f"mod7_same={c % 7 == c1 % 7} "
            f"mod49_same={c % 49 == c1 % 49} "
            f"mod343_same={c % 343 == c1 % 343}"
        )

    # ------------------------------------------------------------------
    # 6. COMPLETE MOD-49 PERIOD
    # ------------------------------------------------------------------

    first_return_49 = None

    for k in range(1, K_MAX + 1):

        c = orbit_value(k)

        if (
            k > 1
            and c % 49 == c1 % 49
        ):
            first_return_49 = k
            break

    print()
    print("=" * 78)
    print("6. FIRST RETURN TO c1")
    print("=" * 78)

    print(
        f"  first_k_same_mod7="
        f"{1 + 6}"
    )

    print(
        f"  first_k_same_mod49="
        f"{first_return_49}"
    )

    # ------------------------------------------------------------------
    # 7. ORBIT DIFFERENCE VALUATION PATTERN
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. v7(c_k-c1) FOR PERIOD-6 RETURNS")
    print("=" * 78)

    for m, k, c, delta_orbit, v_orbit, rho_delta, v_rho in returns:

        print(
            f"  k={k}: "
            f"c_k-c1={delta_orbit} "
            f"v7={v_orbit}"
        )

    # ------------------------------------------------------------------
    # 8. RHO-C_k VS ORBIT RETURN
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. SOURCE-RATIO DISTANCE VS ORBIT RETURN")
    print("=" * 78)

    for m, k, c, delta_orbit, v_orbit, rho_delta, v_rho in returns:

        print(
            f"  k={k}: "
            f"v7(c_k-c1)={v_orbit} "
            f"v7(rho-c_k)={v_rho}"
        )

    # ------------------------------------------------------------------
    # 9. k=7 EXPLANATION
    # ------------------------------------------------------------------

    k7 = 7
    c7 = orbit_value(k7)

    orbit_gap_7 = (
        c7 - c1
    )

    rho_gap_7 = (
        Q1 - c7 * Q3
    )

    print()
    print("=" * 78)
    print("9. k=7 FIRST RETURN ANALYSIS")
    print("=" * 78)

    print(
        f"  c1={c1}"
    )

    print(
        f"  c7={c7}"
    )

    print(
        f"  c7-c1={orbit_gap_7}"
    )

    print(
        f"  v7(c7-c1)={v7(orbit_gap_7)}"
    )

    print(
        f"  rho-c7={rho_gap_7}"
    )

    print(
        f"  v7(rho-c7)={v7(rho_gap_7)}"
    )

    print(
        f"  c7_mod7_equals_c1="
        f"{c7 % 7 == c1 % 7}"
    )

    print(
        f"  c7_mod49_equals_c1="
        f"{c7 % 49 == c1 % 49}"
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
The finite-field orbit has period

    ord_7(3)=6.

Therefore

    c_(1+6m) == c_1 (mod 7).

This explains why k=7 returns to the same residue class as k=1
modulo 7.

However,

    ord_49(3)=42,

so the orbit does not return to c1 modulo 49 until a full 42-step
period.

The experiment therefore separates two phenomena:

    finite-field recurrence:
        k -> k+6,

    49-adic recurrence:
        k -> k+42.

The observed value

    v7(rho-c7)=1

is therefore consistent with a mod-7 orbit return without a mod-49
return.

The exceptional

    v7(rho-c1)=2

remains stronger than the generic period-6 return.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    k7_v = v7(
        rho_gap_7
    )

    final_ok = (
        v7(Q1) == 0
        and v7(Q3) == 0
        and ord7 == 6
        and ord49 == 42
        and ord343 == 294
        and c7 % 7 == c1 % 7
        and c7 % 49 != c1 % 49
        and k7_v == 1
        and v7(Q1 - c1 * Q3) == 2
        and first_return_49 == 43
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
        f"  ord_7_3_is_6={ord7 == 6}"
    )

    print(
        f"  ord_49_3_is_42={ord49 == 42}"
    )

    print(
        f"  ord_343_3_is_294={ord343 == 294}"
    )

    print(
        f"  k7_same_mod7={c7 % 7 == c1 % 7}"
    )

    print(
        f"  k7_differs_mod49={c7 % 49 != c1 % 49}"
    )

    print(
        f"  k7_ratio_distance_v7_is_1={k7_v == 1}"
    )

    print(
        f"  k1_ratio_distance_v7_is_2="
        f"{v7(Q1 - c1 * Q3) == 2}"
    )

    print(
        f"  first_return_mod49_is_43="
        f"{first_return_49 == 43}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 194 COMPLETE")


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

