#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 204 — EXACT SOURCE-RATIO / INTERCEPT CORRESPONDENCE AUDIT
==============================================================================

Experiment 203FR established a stable local lift law:

    R_e(t) = A_e + t  (mod 7),

with

    t_e = -A_e mod 7.

The recent verified intercepts are:

    A = [0, 5, 4, 1]

for the lifts

    e=8 -> 9
    e=9 -> 10
    e=10 -> 11
    e=11 -> 12.

Experiment 204 asks whether A_e can be recovered directly from the
source ratio

    rho = q1/q3

at the corresponding 7-adic precision.

The normalized residual is

    A_e = F(k_e) / 7^e mod 7,

where

    F(k) = 2*3^k*q3 - q1.

Rather than treating A_e as an independent state variable, this
experiment compares it with exact normalized ratio defects.

For each verified exponent k_e, write

    rho - 2*3^k_e = 7^e * delta_e

in Z_7.

Since q3 is a 7-adic unit,

    F(k_e) / 7^e
      = q3 * (2*3^k_e - rho) / 7^e
      = -q3 * delta_e.

Thus A_e should equal the source-ratio defect multiplied by -q3.

The experiment checks this exactly modulo 7 and determines whether
the observed intercept sequence is simply a source-ratio digit sequence
in disguise.

No giant powers of 3 are constructed.
All orbit values are computed modulo powers of 7.
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

# Verified exponent sequence from 203FR.
LEVELS = [
    (8, 3401371),
    (9, 3401371),
    (10, 72578983),
    (11, 798943909),
    (12, 10968052873),
]


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


def ratio_mod(
    modulus: int,
) -> int:
    if modulus == 1:
        return 0

    return (
        Q1
        * inverse_mod(
            Q3,
            modulus,
        )
    ) % modulus


def orbit_mod(
    k: int,
    modulus: int,
) -> int:
    return (
        BASE
        * pow(
            MULT,
            k,
            modulus,
        )
    ) % modulus


def order_3_mod_7e(e: int) -> int:
    return 6 * (
        P ** (e - 1)
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 204 — EXACT SOURCE-RATIO / "
        "INTERCEPT CORRESPONDENCE AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. SOURCE RATIO DATA")
    print("=" * 78)

    print(
        f"  q1={Q1}"
    )

    print(
        f"  q3={Q3}"
    )

    print(
        f"  v7(q1)={valuation_7(Q1)}"
    )

    print(
        f"  v7(q3)={valuation_7(Q3)}"
    )

    # ------------------------------------------------------------------
    # 2. LEVEL-BY-LEVEL SOURCE-RATIO DEFECT
    # ------------------------------------------------------------------

    rows = []

    print()
    print("=" * 78)
    print("2. EXACT SOURCE-RATIO DEFECTS")
    print("=" * 78)

    for e, k in LEVELS:

        modulus_e = P ** e
        modulus_next = P ** (e + 1)

        rho_e = ratio_mod(
            modulus_e,
        )

        rho_next = ratio_mod(
            modulus_next,
        )

        orbit_e = orbit_mod(
            k,
            modulus_e,
        )

        orbit_next = orbit_mod(
            k,
            modulus_next,
        )

        # Since orbit_e == rho_e by construction, the next precision
        # difference is divisible by 7^e.
        residue_next = (
            orbit_next - rho_next
        ) % modulus_next

        if residue_next % modulus_e != 0:
            raise ArithmeticError(
                f"Source-ratio defect not divisible by 7^{e} "
                f"at level e={e}."
            )

        delta_e = (
            residue_next // modulus_e
        ) % P

        # F(k)/7^e mod 7 can be recovered entirely modulo 7^(e+1).
        F_mod = (
            orbit_next
            * Q3
            - Q1
        ) % modulus_next

        if F_mod % modulus_e != 0:
            raise ArithmeticError(
                f"F(k) not divisible by 7^{e}."
            )

        A_e = (
            F_mod // modulus_e
        ) % P

        rows.append(
            {
                "e": e,
                "k": k,
                "modulus": modulus_e,
                "rho_e": rho_e,
                "rho_next": rho_next,
                "orbit_e": orbit_e,
                "orbit_next": orbit_next,
                "delta": delta_e,
                "A": A_e,
            }
        )

        print(
            f"  e={e}: "
            f"k={k} "
            f"rho_mod_7^e={rho_e} "
            f"orbit_mod_7^e={orbit_e} "
            f"rho_defect_digit={delta_e} "
            f"A_e={A_e}"
        )

    # ------------------------------------------------------------------
    # 3. SOURCE-RATIO TO INTERCEPT MAP
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXACT A_e / SOURCE-DEFECT RELATION")
    print("=" * 78)

    relation_checks = []

    q3_mod7 = Q3 % P

    print(
        f"  q3_mod7={q3_mod7}"
    )

    for row in rows:

        delta = row["delta"]
        A = row["A"]

        predicted = (
            (-q3_mod7)
            * delta
        ) % P

        exact = (
            predicted
            == A
        )

        relation_checks.append(
            exact
        )

        print(
            f"  e={row['e']}: "
            f"delta={delta} "
            f"-q3*delta mod7={predicted} "
            f"A={A} "
            f"exact={exact}"
        )

    # ------------------------------------------------------------------
    # 4. LIFT DIGIT RECOVERY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. LIFT DIGIT FROM SOURCE-RATIO DEFECT")
    print("=" * 78)

    digit_checks = []

    for row in rows:

        A = row["A"]

        t = (
            -A
        ) % P

        digit_checks.append(
            t
        )

        print(
            f"  e={row['e']}: "
            f"A={A} "
            f"t=-A mod7={t}"
        )

    # ------------------------------------------------------------------
    # 5. SOURCE-RATIO MOD-7 DIGITS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. SOURCE-RATIO DEFECT SEQUENCE")
    print("=" * 78)

    delta_sequence = [
        row["delta"]
        for row in rows
    ]

    A_sequence = [
        row["A"]
        for row in rows
    ]

    print(
        f"  delta_sequence={delta_sequence}"
    )

    print(
        f"  A_sequence={A_sequence}"
    )

    print(
        f"  q3_mod7={q3_mod7}"
    )

    print(
        f"  -q3_mod7="
        f"{(-q3_mod7) % P}"
    )

    # ------------------------------------------------------------------
    # 6. DIFFERENCE AUDIT
    # ------------------------------------------------------------------

    delta_diffs = [
        (
            delta_sequence[i + 1]
            - delta_sequence[i]
        ) % P
        for i in range(
            len(delta_sequence) - 1
        )
    ]

    A_diffs = [
        (
            A_sequence[i + 1]
            - A_sequence[i]
        ) % P
        for i in range(
            len(A_sequence) - 1
        )
    ]

    print()
    print("=" * 78)
    print("6. DEFECT / INTERCEPT DIFFERENCES")
    print("=" * 78)

    print(
        f"  delta_differences={delta_diffs}"
    )

    print(
        f"  A_differences={A_diffs}"
    )

    # ------------------------------------------------------------------
    # 7. PROJECTIVE INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The intercept

    A_e = F(k_e) / 7^e mod 7

is not an independent object.

Because

    F(k) = q3 * (2*3^k - rho),

the normalized intercept satisfies

    A_e
      = q3 * (2*3^k - rho) / 7^e mod 7.

With the defect convention

    delta_e
      = (2*3^k - rho) / 7^e mod 7,

we obtain

    A_e = q3 * delta_e mod 7.

The sign depends only on the chosen defect orientation.

This experiment uses the explicit orientation above and checks the
identity exactly.

A positive result means the apparently level-dependent intercept
sequence is simply the source-ratio defect sequence viewed through
multiplication by the fixed 7-adic unit q3.

That would explain why the slope stays stable while the intercept
changes.
"""
    )

    # ------------------------------------------------------------------
    # 8. FINAL
    # ------------------------------------------------------------------

    relation_exact = all(
        relation_checks
    )

    digits_exact = (
        [
            (
                -row["A"]
            ) % P
            for row in rows
        ]
        == [
            0,
            2,
            3,
            6,
            3,
        ]
    )

    expected_A = [
        0,
        5,
        4,
        1,
        4,
    ]

    A_exact = (
        A_sequence
        == expected_A
    )

    final_ok = (
        valuation_7(Q1) == 0
        and valuation_7(Q3) == 0
        and relation_exact
        and A_exact
    )

    print()
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  source_defect_to_intercept_exact="
        f"{relation_exact}"
    )

    print(
        f"  A_sequence={A_sequence}"
    )

    print(
        f"  expected_A_sequence={expected_A}"
    )

    print(
        f"  A_sequence_exact="
        f"{A_exact}"
    )

    print(
        f"  recovered_lift_digits="
        f"{digit_checks}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 204 COMPLETE")


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

