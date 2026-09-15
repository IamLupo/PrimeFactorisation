#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 193 — EXACT 7-ADIC ORBIT / 2*3^k SCHUR-CONSTANT AUDIT
==============================================================================

Established:

    rho = q1/q3  in Z_7,

and the terminal Schur constants are

    c_k = 2 * 3^k,

for k = 0,1,2.

Experiment 192 found the unique cancellation

    v7(rho-6)=2,

while

    v7(rho-2)=0,
    v7(rho-18)=0.

Experiment 193 asks whether this should be viewed as an isolated
coincidence with c=6, or as a relation between rho and the multiplicative
3-orbit

    2, 6, 18, 54, 162, ...

modulo powers of 7.

For each k in a finite exact range, compute

    c_k = 2*3^k,

and

    v7(rho-c_k)

via the equivalent integer quantity

    q1 - c_k*q3.

Then determine:

    * the valuation profile of the orbit;
    * the unique closest orbit element;
    * the first k beyond 0,1,2 producing another 7-adic approximation;
    * the orbit residues modulo 7, 49, 343;
    * the multiplicative order of 3 modulo 7 and 49.

This is an exact finite orbit audit only.
No extrapolation is claimed.
"""

from __future__ import annotations

import sys


# ============================================================================
# DATA
# ============================================================================

Q1 = 29144191
Q3 = 24794967


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


def ratio_mod(a, b, modulus):
    if modulus == 1:
        return 0

    inv = inverse_mod(b, modulus)
    return (a * inv) % modulus


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

        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s

    if old_r != 1:
        raise ArithmeticError(
            f"{a} is not invertible modulo {m}."
        )

    return old_s % m


def multiplicative_order(a, m):
    """
    Small exact order search. Requires gcd(a,m)=1.
    """
    if a % m == 0:
        raise ArithmeticError(
            f"{a} is not a unit modulo {m}."
        )

    x = a % m

    for n in range(1, 10000):

        if x == 1:
            return n

        x = (
            x * a
        ) % m

    raise ArithmeticError(
        f"Order search exceeded bound for {a} mod {m}."
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 193 — EXACT 7-ADIC ORBIT / "
        "2*3^k SCHUR-CONSTANT AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE RATIO
    # ------------------------------------------------------------------

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
        f"  v7(q1)={v7(Q1)}"
    )

    print(
        f"  v7(q3)={v7(Q3)}"
    )

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
    # 2. MULTIPLICATIVE ORDER
    # ------------------------------------------------------------------

    order7 = multiplicative_order(
        3,
        7,
    )

    order49 = multiplicative_order(
        3,
        49,
    )

    print()
    print("=" * 78)
    print("2. MULTIPLICATIVE ORDER OF 3")
    print("=" * 78)

    print(
        f"  ord_7(3)={order7}"
    )

    print(
        f"  ord_49(3)={order49}"
    )

    # ------------------------------------------------------------------
    # 3. ORBIT VALUES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXACT ORBIT c_k = 2*3^k")
    print("=" * 78)

    orbit = []

    for k in range(0, 10):

        c = 2 * (3 ** k)

        epsilon = (
            Q1
            - c * Q3
        )

        valuation = v7(
            epsilon
        )

        orbit.append(
            (
                k,
                c,
                epsilon,
                valuation,
            )
        )

        print(
            f"  k={k}: "
            f"c_k={c} "
            f"epsilon=q1-c_k*q3={epsilon} "
            f"v7={valuation}"
        )

    # ------------------------------------------------------------------
    # 4. VALUATION PROFILE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. ORBIT 7-ADIC VALUATION PROFILE")
    print("=" * 78)

    for k, c, epsilon, valuation in orbit:

        print(
            f"  k={k}: v7(rho-c_k)={valuation}"
        )

    finite_valuations = [
        (
            k,
            valuation,
        )
        for k, _, _, valuation in orbit
        if valuation is not None
    ]

    max_v = max(
        valuation
        for _, valuation in finite_valuations
    )

    closest = [
        k
        for k, valuation in finite_valuations
        if valuation == max_v
    ]

    print(
        f"  maximum_finite_v7={max_v}"
    )

    print(
        f"  closest_k_values={closest}"
    )

    # ------------------------------------------------------------------
    # 5. ORBIT RESIDUES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. ORBIT RESIDUES MOD 7,49,343")
    print("=" * 78)

    for modulus in (
        7,
        49,
        343,
    ):

        residues = [
            (
                k,
                c % modulus,
            )
            for k, c, _, _ in orbit
        ]

        print(
            f"  modulus={modulus}: "
            f"{residues}"
        )

    # ------------------------------------------------------------------
    # 6. MATCHING ORBIT VALUES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. ORBIT MATCHES TO rho")
    print("=" * 78)

    for modulus in (
        7,
        49,
        343,
    ):

        rho = ratio_mod(
            Q1,
            Q3,
            modulus,
        )

        matches = [
            k
            for k, c, _, _ in orbit
            if c % modulus == rho
        ]

        print(
            f"  modulus={modulus}: "
            f"rho={rho} "
            f"matching_k={matches}"
        )

    # ------------------------------------------------------------------
    # 7. EXCEPTIONAL k=1 LIFT
    # ------------------------------------------------------------------

    k1 = orbit[1]
    k1_c = k1[1]
    k1_epsilon = k1[2]
    k1_v = k1[3]

    print()
    print("=" * 78)
    print("7. k=1 EXCEPTIONAL LIFT")
    print("=" * 78)

    print(
        f"  c1={k1_c}"
    )

    print(
        f"  epsilon1={k1_epsilon}"
    )

    print(
        f"  v7(epsilon1)={k1_v}"
    )

    if k1_v is not None:

        normalized = (
            k1_epsilon
            // (7 ** k1_v)
        )

        print(
            f"  normalized_unit={normalized}"
        )

        print(
            f"  normalized_unit_mod7="
            f"{normalized % 7}"
        )

    # ------------------------------------------------------------------
    # 8. ORBIT SEPARATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. ADJACENT ORBIT CONSTANT SEPARATION")
    print("=" * 78)

    for k in range(0, 8):

        c_a = 2 * (3 ** k)
        c_b = 2 * (3 ** (k + 1))

        diff = c_b - c_a

        print(
            f"  k={k}: "
            f"c_(k+1)-c_k={diff} "
            f"v7={v7(diff)}"
        )

    # ------------------------------------------------------------------
    # 9. PROJECTIVE INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. PROJECTIVE ORBIT INTERPRETATION")
    print("=" * 78)

    print(
        """
The terminal Schur constants are the first three elements of the
multiplicative orbit

    c_k = 2*3^k.

Because q3 is a 7-adic unit,

    v7(q1-c_k*q3)
      = v7(rho-c_k).

Thus the entire question is equivalent to locating the source ratio
rho in the 7-adic completion relative to the orbit generated by
multiplication by 3.

The previous experiments showed that

    rho is exceptionally close to c_1 = 6.

Experiment 193 tests whether that proximity is unique over a larger
finite orbit segment or whether later orbit points also produce
nontrivial 7-adic cancellation.

The multiplicative-order data provide the finite-field periodic
background for the orbit.
"""
    )

    # ------------------------------------------------------------------
    # 10. FINAL
    # ------------------------------------------------------------------

    k0_v = orbit[0][3]
    k1_v = orbit[1][3]
    k2_v = orbit[2][3]

    final_ok = (
        v7(Q1) == 0
        and v7(Q3) == 0
        and order7 == 6
        and order49 == 42
        and k0_v == 0
        and k1_v == 2
        and k2_v == 0
        and max_v >= 2
        and 1 in closest
    )

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  q1_7adic_unit={v7(Q1) == 0}"
    )

    print(
        f"  q3_7adic_unit={v7(Q3) == 0}"
    )

    print(
        f"  ord_7_3_is_6={order7 == 6}"
    )

    print(
        f"  ord_49_3_is_42={order49 == 42}"
    )

    print(
        f"  k0_v7_is_0={k0_v == 0}"
    )

    print(
        f"  k1_v7_is_2={k1_v == 2}"
    )

    print(
        f"  k2_v7_is_0={k2_v == 0}"
    )

    print(
        f"  k1_is_among_closest="
        f"{1 in closest}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 193 COMPLETE")


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
