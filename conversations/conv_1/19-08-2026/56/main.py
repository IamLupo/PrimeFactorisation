#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 223R — EXACT TERMINAL SOURCE /
                 FIRST PRIME-ADIC OBSTRUCTION AUDIT
==============================================================================

Corrected single-file main.py.

Fixes:
    * imports Fraction explicitly;
    * contains no external Python dependencies;
    * uses only main.py;
    * keeps the exact q-normalized source analysis from Experiment 223.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd
import sys


# ============================================================================
# KNOWN INSTANCE
# ============================================================================

P = 2
Q = 3
N = P * Q

Q1 = 29144191
Q3 = 24794967

SEVEN = 7


# ============================================================================
# EXACT HELPERS
# ============================================================================

def valuation_p(
    x: int,
    p: int,
) -> int | None:

    x = int(x)
    p = int(p)

    if p < 2:
        raise ValueError(
            "p must be >= 2."
        )

    if x == 0:
        return None

    x = abs(x)
    v = 0

    while x % p == 0:
        x //= p
        v += 1

    return v


def gcd_int(
    a: int,
    b: int,
) -> int:

    a = abs(int(a))
    b = abs(int(b))

    while b:
        a, b = (
            b,
            a % b,
        )

    return int(a)


def divisibility_report(
    x: int,
    base: int,
    max_e: int,
):

    rows = []

    for e in range(
        1,
        max_e + 1,
    ):

        modulus = int(
            base ** e
        )

        rows.append(
            {
                "e": e,
                "modulus": modulus,
                "divisible": (
                    x % modulus == 0
                ),
                "residue": (
                    x % modulus
                ),
            }
        )

    return rows


def exact_division(
    x: int,
    d: int,
) -> tuple[bool, int | None]:

    if d == 0:
        return (
            False,
            None,
        )

    if x % d != 0:
        return (
            False,
            None,
        )

    return (
        True,
        int(x // d),
    )


# ============================================================================
# TERMINAL SCHUR ROW
# ============================================================================

def schur_row(
    q1: int,
    q3: int,
) -> list[int]:

    return [
        int(
            q1 - 2 * q3
        ),
        int(
            q1 - 6 * q3
        ),
        int(
            q1 - 18 * q3
        ),
    ]


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 223R — EXACT TERMINAL SOURCE / "
        "FIRST PRIME-ADIC OBSTRUCTION AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE DATA
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. SOURCE DATA")
    print("=" * 78)

    print(
        f"  p={P}"
    )

    print(
        f"  q={Q}"
    )

    print(
        f"  n={N}"
    )

    print(
        f"  q1={Q1}"
    )

    print(
        f"  q3={Q3}"
    )

    # ------------------------------------------------------------------
    # 2. NORMALIZED q-SOURCE
    # ------------------------------------------------------------------

    if Q3 % Q != 0:
        raise ArithmeticError(
            "q3 is not divisible by q."
        )

    U = int(
        Q3 // Q
    )

    print()
    print("=" * 78)
    print("2. NORMALIZED q-SOURCE")
    print("=" * 78)

    print(
        f"  u=q3/q={U}"
    )

    print(
        f"  v_q(q3)="
        f"{valuation_p(Q3,Q)}"
    )

    print(
        f"  v_p(u)="
        f"{valuation_p(U,P)}"
    )

    print(
        f"  v_q(u)="
        f"{valuation_p(U,Q)}"
    )

    # ------------------------------------------------------------------
    # 3. FIRST SOURCE OBSTRUCTIONS
    # ------------------------------------------------------------------

    q1_minus_1 = int(
        Q1 - 1
    )

    u_minus_1 = int(
        U - 1
    )

    print()
    print("=" * 78)
    print("3. FIRST SOURCE OBSTRUCTIONS")
    print("=" * 78)

    print(
        f"  q1-1={q1_minus_1}"
    )

    print(
        f"  u-1={u_minus_1}"
    )

    print(
        f"  v_p(q1-1)="
        f"{valuation_p(q1_minus_1,P)}"
    )

    print(
        f"  v_q(q1-1)="
        f"{valuation_p(q1_minus_1,Q)}"
    )

    print(
        f"  v_p(u-1)="
        f"{valuation_p(u_minus_1,P)}"
    )

    print(
        f"  v_q(u-1)="
        f"{valuation_p(u_minus_1,Q)}"
    )

    # ------------------------------------------------------------------
    # 4. q1-1 PRECISION LADDER
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. q1-1 p-ADIC / q-ADIC PRECISION LADDER")
    print("=" * 78)

    print(
        "  p-adic:"
    )

    for row in divisibility_report(
        q1_minus_1,
        P,
        8,
    ):

        print(
            f"    e={row['e']} "
            f"modulus={row['modulus']} "
            f"divisible={row['divisible']} "
            f"residue={row['residue']}"
        )

    print(
        "  q-adic:"
    )

    for row in divisibility_report(
        q1_minus_1,
        Q,
        8,
    ):

        print(
            f"    e={row['e']} "
            f"modulus={row['modulus']} "
            f"divisible={row['divisible']} "
            f"residue={row['residue']}"
        )

    # ------------------------------------------------------------------
    # 5. u-1 PRECISION LADDER
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. u-1 p-ADIC / q-ADIC PRECISION LADDER")
    print("=" * 78)

    print(
        "  p-adic:"
    )

    for row in divisibility_report(
        u_minus_1,
        P,
        8,
    ):

        print(
            f"    e={row['e']} "
            f"modulus={row['modulus']} "
            f"divisible={row['divisible']} "
            f"residue={row['residue']}"
        )

    print(
        "  q-adic:"
    )

    for row in divisibility_report(
        u_minus_1,
        Q,
        8,
    ):

        print(
            f"    e={row['e']} "
            f"modulus={row['modulus']} "
            f"divisible={row['divisible']} "
            f"residue={row['residue']}"
        )

    # ------------------------------------------------------------------
    # 6. COMBINED pq PRECISION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. COMBINED pq PRECISION")
    print("=" * 78)

    print(
        f"  q1-1 mod pq="
        f"{(Q1 - 1) % N}"
    )

    print(
        f"  u-1 mod pq="
        f"{(U - 1) % N}"
    )

    print(
        f"  q1_equals_1_mod_pq="
        f"{Q1 % N == 1 % N}"
    )

    print(
        f"  u_equals_1_mod_pq="
        f"{U % N == 1 % N}"
    )

    # ------------------------------------------------------------------
    # 7. q3 RECONSTRUCTION
    # ------------------------------------------------------------------

    reconstructed_q3 = int(
        Q * U
    )

    q3_reconstruction_exact = (
        reconstructed_q3 == Q3
    )

    print()
    print("=" * 78)
    print("7. q3 = q*u EXACT RECONSTRUCTION")
    print("=" * 78)

    print(
        f"  q*u={reconstructed_q3}"
    )

    print(
        f"  q3={Q3}"
    )

    print(
        f"  exact="
        f"{q3_reconstruction_exact}"
    )

    # ------------------------------------------------------------------
    # 8. TERMINAL SCHUR ROW
    # ------------------------------------------------------------------

    s0, s1, s2 = schur_row(
        Q1,
        Q3,
    )

    print()
    print("=" * 78)
    print("8. TERMINAL SCHUR ROW")
    print("=" * 78)

    print(
        f"  s0={s0}"
    )

    print(
        f"  s1={s1}"
    )

    print(
        f"  s2={s2}"
    )

    # ------------------------------------------------------------------
    # 9. SCHUR ROW IN q,u COORDINATES
    # ------------------------------------------------------------------

    reconstructed_s = [
        int(
            Q1 - 2 * Q * U
        ),
        int(
            Q1 - 6 * Q * U
        ),
        int(
            Q1 - 18 * Q * U
        ),
    ]

    schur_reconstruction_exact = (
        reconstructed_s
        == [s0, s1, s2]
    )

    print()
    print("=" * 78)
    print("9. SCHUR ROW AFTER q-NORMALIZATION")
    print("=" * 78)

    print(
        f"  reconstructed_from_q_u="
        f"{reconstructed_s}"
    )

    print(
        f"  original="
        f"{[s0,s1,s2]}"
    )

    print(
        f"  exact="
        f"{schur_reconstruction_exact}"
    )

    # ------------------------------------------------------------------
    # 10. q-ADIC SCHUR PROFILE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. q-ADIC SCHUR COEFFICIENT PROFILE")
    print("=" * 78)

    for name, value in (
        ("s0", s0),
        ("s1", s1),
        ("s2", s2),
    ):

        print(
            f"  {name}: "
            f"v_q={valuation_p(value,Q)} "
            f"v_p={valuation_p(value,P)} "
            f"mod_q={value % Q}"
        )

    # ------------------------------------------------------------------
    # 11. DIFFERENCE LAYER
    # ------------------------------------------------------------------

    d01 = int(
        s1 - s0
    )

    d12 = int(
        s2 - s1
    )

    d02 = int(
        s2 - s0
    )

    difference_exact = (
        d01 == -4 * Q * U
        and
        d12 == -12 * Q * U
        and
        d02 == -16 * Q * U
    )

    geometric_difference_law = (
        d12 == 3 * d01
    )

    print()
    print("=" * 78)
    print("11. SCHUR DIFFERENCE LAYERS")
    print("=" * 78)

    print(
        f"  d01={d01}"
    )

    print(
        f"  d12={d12}"
    )

    print(
        f"  d02={d02}"
    )

    print(
        f"  d01/(-4*q*u)="
        f"{Fraction(d01, -4 * Q * U)}"
    )

    print(
        f"  d12/(-12*q*u)="
        f"{Fraction(d12, -12 * Q * U)}"
    )

    print(
        f"  d02/(-16*q*u)="
        f"{Fraction(d02, -16 * Q * U)}"
    )

    print(
        f"  difference_formulas_exact="
        f"{difference_exact}"
    )

    print(
        f"  d12=3*d01="
        f"{geometric_difference_law}"
    )

    # ------------------------------------------------------------------
    # 12. FIRST-OBSTRUCTION SUMMARY
    # ------------------------------------------------------------------

    vp_q1 = valuation_p(
        q1_minus_1,
        P,
    )

    vq_q1 = valuation_p(
        q1_minus_1,
        Q,
    )

    vp_u = valuation_p(
        u_minus_1,
        P,
    )

    vq_u = valuation_p(
        u_minus_1,
        Q,
    )

    q1_breaks_at_q2 = (
        Q1 % (Q * Q)
        != 1 % (Q * Q)
    )

    u_survives_q2 = (
        U % (Q * Q)
        == 1 % (Q * Q)
    )

    print()
    print("=" * 78)
    print("12. FIRST-OBSTRUCTION SUMMARY")
    print("=" * 78)

    print(
        f"  v_p(q1-1)={vp_q1}"
    )

    print(
        f"  v_q(q1-1)={vq_q1}"
    )

    print(
        f"  v_p(u-1)={vp_u}"
    )

    print(
        f"  v_q(u-1)={vq_u}"
    )

    print(
        f"  q1_equals_1_mod_pq="
        f"{Q1 % N == 1 % N}"
    )

    print(
        f"  q1_breaks_at_q2="
        f"{q1_breaks_at_q2}"
    )

    print(
        f"  u_equals_1_mod_q2="
        f"{u_survives_q2}"
    )

    # ------------------------------------------------------------------
    # 13. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("13. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The source layer has an exact normalization

    q3 = q*u,

where u is a p-adic and q-adic unit.

The first congruence layer is

    q1 = 1 (mod p),
    q1 = 1 (mod q),

and

    u = 1 (mod p),
    u = 1 (mod q).

The next question is the depth at which the two normalized source
coordinates leave this unit neighborhood.

The present computation therefore separates:

    q-adic divisibility of q3,

from

    higher q-adic structure of the normalized unit u.

The Schur row can then be written entirely in the coordinates

    (q1,u):

        [q1-2qu,
         q1-6qu,
         q1-18qu].

If this same normalization and congruence pattern appears in further
genuine n=pq cases, it will provide a plausible coordinate system for
the eventual symbolic family theorem.

At this stage it remains a finite observation for n=6.
"""
    )

    # ------------------------------------------------------------------
    # 14. FINAL
    # ------------------------------------------------------------------

    base_congruence = (
        Q1 % N == 1
        and
        U % N == 1
    )

    q_order_one = (
        valuation_p(
            Q3,
            Q,
        )
        == 1
    )

    obstruction_pattern = (
        q1_breaks_at_q2
        and
        u_survives_q2
    )

    final_ok = (
        q3_reconstruction_exact
        and
        schur_reconstruction_exact
        and
        difference_exact
        and
        geometric_difference_law
        and
        base_congruence
        and
        q_order_one
        and
        obstruction_pattern
    )

    print()
    print("=" * 78)
    print("14. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  q3=q*u_exact="
        f"{q3_reconstruction_exact}"
    )

    print(
        f"  Schur_q_u_reconstruction_exact="
        f"{schur_reconstruction_exact}"
    )

    print(
        f"  difference_formulas_exact="
        f"{difference_exact}"
    )

    print(
        f"  geometric_difference_law="
        f"{geometric_difference_law}"
    )

    print(
        f"  q3_exact_q_order_1="
        f"{q_order_one}"
    )

    print(
        f"  q1_and_u_equal_1_mod_pq="
        f"{base_congruence}"
    )

    print(
        f"  q1_breaks_mod_q2_while_u_survives="
        f"{obstruction_pattern}"
    )

    print(
        "  general_n_pq_theorem_established=False"
    )

    print(
        "  reason=single_known_n6_instance"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 223R COMPLETE")


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:

        print(
            "\nInterrupted."
        )

        sys.exit(130)

    except Exception as exc:

        print(
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )

        raise