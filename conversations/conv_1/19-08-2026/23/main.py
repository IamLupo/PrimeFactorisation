#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 192 — EXACT 7-ADIC PROJECTIVE SCHUR-ROW AUDIT
==============================================================================

Established:

    q1 = 29144191
    q3 = 24794967

    s0 = q1 - 2*q3
    s1 = q1 - 6*q3
    s2 = q1 - 18*q3

and

    v7(q1)=v7(q3)=0
    v7(s0)=0
    v7(s1)=2
    v7(s2)=0.

Also:

    q1/q3 = 6 + 49*delta7,

with

    delta7 = -2441339,
    delta7 mod 7 = 2.

Experiment 192 compares the three Schur coefficients through the
projective source-ratio quantity

    rho = q1/q3.

For each Schur constant c in {2,6,18}, define

    epsilon_c = rho - c.

Then

    s_c = q3 * epsilon_c.

Because q3 is a 7-adic unit,

    v7(s_c) = v7(epsilon_c).

The experiment therefore removes q3 entirely and studies only the
distance of rho from the three Schur constants in Z_7.

It determines:

    * exact epsilon_c residues;
    * v7(epsilon_c);
    * normalized unit epsilon_c / 7^v7;
    * first separating precision;
    * pairwise differences among the three constants;
    * whether c=6 is uniquely closest to rho 7-adically.

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

SCHUR_CONSTANTS = {
    "c0": 2,
    "c1": 6,
    "c2": 18,
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


def ratio_mod(a, b, p, e):
    modulus = p ** e

    if modulus == 1:
        return 0

    return (
        a * inverse_mod(b, modulus)
    ) % modulus


def first_equal_precision(a, b, p, max_e=12):
    for e in range(1, max_e + 1):
        m = p ** e

        if a % m == b % m:
            return e

    return None


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 192 — EXACT 7-ADIC PROJECTIVE SCHUR-ROW AUDIT"
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

    rho_7 = ratio_mod(
        Q1,
        Q3,
        7,
        1,
    )

    rho_49 = ratio_mod(
        Q1,
        Q3,
        7,
        2,
    )

    rho_343 = ratio_mod(
        Q1,
        Q3,
        7,
        3,
    )

    print(
        f"  rho_mod7={rho_7}"
    )

    print(
        f"  rho_mod49={rho_49}"
    )

    print(
        f"  rho_mod343={rho_343}"
    )

    # ------------------------------------------------------------------
    # 2. EXACT SCHUR COEFFICIENTS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. EXACT SCHUR COEFFICIENTS")
    print("=" * 78)

    s0 = Q1 - 2 * Q3
    s1 = Q1 - 6 * Q3
    s2 = Q1 - 18 * Q3

    S = {
        "c0": s0,
        "c1": s1,
        "c2": s2,
    }

    for name, value in S.items():

        print(
            f"  {name}={value} "
            f"v7={v7(value)} "
            f"mod7={value % 7}"
        )

    # ------------------------------------------------------------------
    # 3. PROJECTIVE DISTANCES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. PROJECTIVE 7-ADIC DISTANCES")
    print("=" * 78)

    distances = {}

    for name, c in SCHUR_CONSTANTS.items():

        epsilon = Q1 - c * Q3
        e = v7(epsilon)

        distances[name] = (
            epsilon,
            e,
        )

        print(
            f"  {name}: c={c} "
            f"epsilon={epsilon} "
            f"v7(epsilon)={e}"
        )

        if e is not None:
            unit = (
                epsilon // (7 ** e)
            )

            print(
                f"       normalized_unit={unit} "
                f"unit_mod7={unit % 7}"
            )

    # ------------------------------------------------------------------
    # 4. SCHUR / PROJECTIVE VALUATION MATCH
    # ------------------------------------------------------------------

    valuation_match = True

    print()
    print("=" * 78)
    print("4. SCHUR / PROJECTIVE VALUATION MATCH")
    print("=" * 78)

    for name in ("c0", "c1", "c2"):

        s_val = v7(
            S[name]
        )

        eps_val = distances[name][1]

        exact = (
            s_val == eps_val
        )

        valuation_match &= exact

        print(
            f"  {name}: "
            f"v7(s)={s_val} "
            f"v7(epsilon)={eps_val} "
            f"exact={exact}"
        )

    # ------------------------------------------------------------------
    # 5. UNIQUE CLOSEST CONSTANT
    # ------------------------------------------------------------------

    finite_vals = {
        name: value[1]
        for name, value in distances.items()
        if value[1] is not None
    }

    max_distance = max(
        finite_vals.values()
    )

    closest = [
        name
        for name, value in finite_vals.items()
        if value == max_distance
    ]

    print()
    print("=" * 78)
    print("5. CLOSEST 7-ADIC SCHUR CONSTANT")
    print("=" * 78)

    print(
        f"  valuations={finite_vals}"
    )

    print(
        f"  maximum_v7_distance={max_distance}"
    )

    print(
        f"  closest_constants={closest}"
    )

    # ------------------------------------------------------------------
    # 6. FIRST PRECISION WHERE rho MATCHES EACH CONSTANT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. FIRST MATCH PRECISION")
    print("=" * 78)

    for name, c in SCHUR_CONSTANTS.items():

        matches = []

        for e in range(1, 6):

            rho = ratio_mod(
                Q1,
                Q3,
                7,
                e,
            )

            if rho == c % (7 ** e):
                matches.append(e)

        print(
            f"  {name}: c={c} "
            f"matching_precisions={matches}"
        )

    # ------------------------------------------------------------------
    # 7. CONSTANT SEPARATION
    # ------------------------------------------------------------------

    c0 = SCHUR_CONSTANTS["c0"]
    c1 = SCHUR_CONSTANTS["c1"]
    c2 = SCHUR_CONSTANTS["c2"]

    dc01 = c1 - c0
    dc12 = c2 - c1
    dc02 = c2 - c0

    print()
    print("=" * 78)
    print("7. 7-ADIC SEPARATION OF SCHUR CONSTANTS")
    print("=" * 78)

    print(
        f"  c1-c0={dc01} "
        f"v7={v7(dc01)}"
    )

    print(
        f"  c2-c1={dc12} "
        f"v7={v7(dc12)}"
    )

    print(
        f"  c2-c0={dc02} "
        f"v7={v7(dc02)}"
    )

    # ------------------------------------------------------------------
    # 8. EXACT rho-6 LIFT
    # ------------------------------------------------------------------

    epsilon6 = Q1 - 6 * Q3
    epsilon6_v = v7(epsilon6)
    epsilon6_unit = (
        epsilon6 // (7 ** epsilon6_v)
    )

    print()
    print("=" * 78)
    print("8. EXCEPTIONAL rho-6 LIFT")
    print("=" * 78)

    print(
        f"  rho-6 numerator={epsilon6}"
    )

    print(
        f"  v7(rho-6 numerator)={epsilon6_v}"
    )

    print(
        f"  normalized_unit={epsilon6_unit}"
    )

    print(
        f"  normalized_unit_mod7="
        f"{epsilon6_unit % 7}"
    )

    print(
        f"  rho_mod343={rho_343}"
    )

    print(
        f"  rho_343=6+49*unit_mod7="
        f"{rho_343 == (6 + 49*(epsilon6_unit % 7)) % 343}"
    )

    # ------------------------------------------------------------------
    # 9. FULL TERMINAL ROW PROJECTIVE VIEW
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. PROJECTIVE TERMINAL ROW")
    print("=" * 78)

    print(
        f"  source_ratio_mod343={rho_343}"
    )

    print(
        f"  distance_to_2_mod343="
        f"{(rho_343 - 2) % 343}"
    )

    print(
        f"  distance_to_6_mod343="
        f"{(rho_343 - 6) % 343}"
    )

    print(
        f"  distance_to_18_mod343="
        f"{(rho_343 - 18) % 343}"
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
Because q3 is a 7-adic unit,

    s_c = q3 * (rho-c)

has exactly the same 7-adic valuation as rho-c.

Thus the three terminal Schur coefficients can be viewed projectively:
they measure the 7-adic distance from the source ratio rho to

    2, 6, 18.

The result identifies which Schur constant is exceptionally close to
the source ratio.

In the observed data, c=6 is the unique constant with

    v7(rho-6)=2,

while the other two distances have valuation zero.

That means the exceptional factor 49 in the middle Schur coefficient
is a genuinely projective 7-adic proximity phenomenon.

This reframes the anomaly without introducing any new matrix mechanism.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        v7(Q1) == 0
        and v7(Q3) == 0
        and distances["c0"][1] == 0
        and distances["c1"][1] == 2
        and distances["c2"][1] == 0
        and valuation_match
        and closest == ["c1"]
        and rho_343 == (
            6
            + 49 * (
                (epsilon6 // 49) % 7
            )
        ) % 343
        and v7(dc01) == 0
        and v7(dc12) == 0
        and v7(dc02) == 0
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
        f"  c2_distance_v7_is_0={distances['c0'][1] == 0}"
    )

    print(
        f"  c6_distance_v7_is_2={distances['c1'][1] == 2}"
    )

    print(
        f"  c18_distance_v7_is_0={distances['c2'][1] == 0}"
    )

    print(
        f"  schur_projective_valuation_match="
        f"{valuation_match}"
    )

    print(
        f"  unique_closest_constant={closest == ['c1']}"
    )

    print(
        f"  rho_mod343_lift_exact="
        f"{rho_343 == (6 + 49*((epsilon6 // 49) % 7)) % 343}"
    )

    print(
        f"  Schur_constants_pairwise_7adic_units="
        f"{v7(dc01) == 0 and v7(dc12) == 0 and v7(dc02) == 0}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 192 COMPLETE")


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

