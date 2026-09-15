#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 221R — EXACT n=pq TERMINAL SOURCE /
                 PRIME-PARAMETER PROVENANCE AUDIT
==============================================================================

Corrected single-file main.py.

Fixes:
    * PRINCIPAL_BASE is explicitly defined.
    * No undefined d12_equals_3d01 variable.
    * No external Python files.
    * All exact arithmetic remains integer/rational.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd
import sys


# ============================================================================
# KNOWN n=pq INSTANCE
# ============================================================================

P = 2
Q = 3
N = P * Q

Q1 = 29144191
Q3 = 24794967

SEVEN = 7
PRINCIPAL_BASE = 729

END_E = 12


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


def gcd_many(
    values: list[int],
) -> int:

    values = [
        abs(int(x))
        for x in values
        if int(x) != 0
    ]

    if not values:
        return 0

    g = values[0]

    for value in values[1:]:
        g = gcd(
            g,
            value,
        )

    return int(g)


def exact_ratio(
    a: int,
    b: int,
) -> Fraction:

    return Fraction(
        int(a),
        int(b),
    )


def inverse_mod(
    a: int,
    m: int,
) -> int:

    a = int(a)
    m = int(m)

    if m <= 1:
        return 0

    a %= m

    if a == 0:
        raise ArithmeticError(
            f"Cannot invert 0 modulo {m}."
        )

    old_r = a
    r = m

    old_s = 1
    s = 0

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

    return int(
        old_s % m
    )


def ratio_mod(
    a: int,
    b: int,
    modulus: int,
) -> int | None:

    modulus = int(modulus)

    bb = int(b) % modulus

    if gcd(
        bb,
        modulus,
    ) != 1:

        return None

    return int(
        (
            (int(a) % modulus)
            * inverse_mod(
                bb,
                modulus,
            )
        )
        % modulus
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


def recover_source_from_schur(
    s0: int,
    s1: int,
    s2: int,
):

    q3_a = Fraction(
        s0 - s1,
        4,
    )

    q3_b = Fraction(
        s1 - s2,
        12,
    )

    q1_a = (
        Fraction(s0)
        + 2 * q3_a
    )

    q1_b = Fraction(
        3 * s1 - s2,
        2,
    )

    return (
        q1_a,
        q1_b,
        q3_a,
        q3_b,
    )


# ============================================================================
# 7-ADIC PROJECTIVE DATA
# ============================================================================

def orbit_mod(
    k: int,
    modulus: int,
) -> int:

    return int(
        (
            2
            * pow(
                3,
                int(k),
                int(modulus),
            )
        )
        % int(modulus)
    )


def find_k1(
    q1: int,
    q3: int,
) -> int:

    rho = ratio_mod(
        q1,
        q3,
        SEVEN,
    )

    if rho is None:
        raise ArithmeticError(
            "q3 is not invertible modulo 7."
        )

    matches = []

    for k in range(6):

        if orbit_mod(
            k,
            SEVEN,
        ) == rho:

            matches.append(
                int(k)
            )

    if len(matches) != 1:
        raise ArithmeticError(
            f"Expected unique k1; rho={rho}, matches={matches}."
        )

    return int(
        matches[0]
    )


# ============================================================================
# 7-ADIC LOGARITHMIC REFERENCE
# ============================================================================

def fraction_valuation_7(
    x: Fraction,
) -> int | None:

    if x == 0:
        return None

    vn = valuation_p(
        x.numerator,
        7,
    )

    vd = valuation_p(
        x.denominator,
        7,
    )

    vn = 0 if vn is None else vn
    vd = 0 if vd is None else vd

    return int(
        vn - vd
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 221R — EXACT n=pq TERMINAL SOURCE / "
        "PRIME-PARAMETER PROVENANCE AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. PRIME PARAMETERS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. n=pq PARAMETERS")
    print("=" * 78)

    print(
        f"  p={P}"
    )

    print(
        f"  q={Q}"
    )

    print(
        f"  n=p*q={N}"
    )

    print(
        f"  q1={Q1}"
    )

    print(
        f"  q3={Q3}"
    )

    # ------------------------------------------------------------------
    # 2. SOURCE CONTENT
    # ------------------------------------------------------------------

    source_gcd = gcd(
        Q1,
        Q3,
    )

    print()
    print("=" * 78)
    print("2. SOURCE CONTENT / PRIME VALUATIONS")
    print("=" * 78)

    print(
        f"  gcd(q1,q3)={source_gcd}"
    )

    for value_name, value in (
        ("q1", Q1),
        ("q3", Q3),
    ):

        print(
            f"  {value_name}: "
            f"v2={valuation_p(value,2)} "
            f"v3={valuation_p(value,3)} "
            f"v5={valuation_p(value,5)} "
            f"v7={valuation_p(value,7)} "
            f"v17={valuation_p(value,17)}"
        )

    print(
        f"  q1_mod_p={Q1 % P}"
    )

    print(
        f"  q3_mod_p={Q3 % P}"
    )

    print(
        f"  q1_mod_q={Q1 % Q}"
    )

    print(
        f"  q3_mod_q={Q3 % Q}"
    )

    print(
        f"  q1_mod_n={Q1 % N}"
    )

    print(
        f"  q3_mod_n={Q3 % N}"
    )

    # ------------------------------------------------------------------
    # 3. LOW-COMPLEXITY SOURCE COMBINATIONS
    # ------------------------------------------------------------------

    combinations = {
        "q1-q3": Q1 - Q3,
        "q1+q3": Q1 + Q3,
        "q1-q*q3": Q1 - Q * Q3,
        "q1-p*q3": Q1 - P * Q3,
        "q1-n*q3": Q1 - N * Q3,
        "p*q1-q3": P * Q1 - Q3,
        "q*q1-q3": Q * Q1 - Q3,
        "n*q1-q3": N * Q1 - Q3,
    }

    print()
    print("=" * 78)
    print("3. LOW-COMPLEXITY SOURCE RELATIONS")
    print("=" * 78)

    for name, value in combinations.items():

        print(
            f"  {name}={value} "
            f"gcd_with_q1={gcd(value,Q1)} "
            f"gcd_with_q3={gcd(value,Q3)}"
        )

    # ------------------------------------------------------------------
    # 4. TERMINAL SCHUR ROW
    # ------------------------------------------------------------------

    s0, s1, s2 = schur_row(
        Q1,
        Q3,
    )

    print()
    print("=" * 78)
    print("4. EXACT TERMINAL SCHUR ROW")
    print("=" * 78)

    print(
        f"  s0=q1-2*q3={s0}"
    )

    print(
        f"  s1=q1-6*q3={s1}"
    )

    print(
        f"  s2=q1-18*q3={s2}"
    )

    print(
        f"  gcd(s0,s1,s2)="
        f"{gcd_many([s0,s1,s2])}"
    )

    # ------------------------------------------------------------------
    # 5. DIFFERENCE GEOMETRY
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

    geometry_exact = (
        d01 == -4 * Q3
        and d12 == -12 * Q3
        and d02 == -16 * Q3
    )

    geometric_difference_law = (
        d12 == 3 * d01
    )

    print()
    print("=" * 78)
    print("5. EXACT SCHUR DIFFERENCE GEOMETRY")
    print("=" * 78)

    print(
        f"  d01=s1-s0={d01}"
    )

    print(
        f"  d12=s2-s1={d12}"
    )

    print(
        f"  d02=s2-s0={d02}"
    )

    print(
        f"  d01/(-4*q3)="
        f"{Fraction(d01,-4*Q3)}"
    )

    print(
        f"  d12/(-12*q3)="
        f"{Fraction(d12,-12*Q3)}"
    )

    print(
        f"  d02/(-16*q3)="
        f"{Fraction(d02,-16*Q3)}"
    )

    print(
        f"  d12=3*d01="
        f"{geometric_difference_law}"
    )

    # ------------------------------------------------------------------
    # 6. SOURCE RECOVERY
    # ------------------------------------------------------------------

    (
        q1_a,
        q1_b,
        q3_a,
        q3_b,
    ) = recover_source_from_schur(
        s0,
        s1,
        s2,
    )

    source_recovery_exact = (
        q1_a == Q1
        and q1_b == Q1
        and q3_a == Q3
        and q3_b == Q3
    )

    print()
    print("=" * 78)
    print("6. EXACT SOURCE RECOVERY FROM SCHUR ROW")
    print("=" * 78)

    print(
        f"  q3_from_s0_s1={q3_a}"
    )

    print(
        f"  q3_from_s1_s2={q3_b}"
    )

    print(
        f"  q1_from_s0_s1={q1_a}"
    )

    print(
        f"  q1_from_s1_s2={q1_b}"
    )

    print(
        f"  source_recovery_exact="
        f"{source_recovery_exact}"
    )

    # ------------------------------------------------------------------
    # 7. PRIME PARAMETER TESTS
    # ------------------------------------------------------------------

    prime_parameter_tests = {
        "q1_divisible_by_p": Q1 % P == 0,
        "q1_divisible_by_q": Q1 % Q == 0,
        "q3_divisible_by_p": Q3 % P == 0,
        "q3_divisible_by_q": Q3 % Q == 0,

        "s0_divisible_by_p": s0 % P == 0,
        "s0_divisible_by_q": s0 % Q == 0,
        "s1_divisible_by_p": s1 % P == 0,
        "s1_divisible_by_q": s1 % Q == 0,
        "s2_divisible_by_p": s2 % P == 0,
        "s2_divisible_by_q": s2 % Q == 0,
    }

    print()
    print("=" * 78)
    print("7. PRIME-PARAMETER DIVISIBILITY AUDIT")
    print("=" * 78)

    for name, value in prime_parameter_tests.items():

        print(
            f"  {name}={value}"
        )

    # ------------------------------------------------------------------
    # 8. EXACT RATIO AUDIT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. EXACT RATIO AUDIT")
    print("=" * 78)

    ratio_tests = {
        "q1/q3": exact_ratio(
            Q1,
            Q3,
        ),
        "q1/(q*q3)": exact_ratio(
            Q1,
            Q * Q3,
        ),
        "q1/(p*q3)": exact_ratio(
            Q1,
            P * Q3,
        ),
        "q1/(n*q3)": exact_ratio(
            Q1,
            N * Q3,
        ),
        "s1/s0": exact_ratio(
            s1,
            s0,
        ),
        "s2/s1": exact_ratio(
            s2,
            s1,
        ),
    }

    for name, value in ratio_tests.items():

        print(
            f"  {name}={value}"
        )

    # ------------------------------------------------------------------
    # 9. 7-ADIC PROJECTIVE REFERENCE
    # ------------------------------------------------------------------

    rho7 = ratio_mod(
        Q1,
        Q3,
        7,
    )

    k1 = find_k1(
        Q1,
        Q3,
    )

    normalized_729_increment = int(
        (
            PRINCIPAL_BASE - 1
        )
        // 7
    )

    predicted_slope = int(
        (
            rho7
            * normalized_729_increment
        )
        % 7
    )

    projective_slope_exact = (
        predicted_slope == 1
    )

    print()
    print("=" * 78)
    print("9. 7-ADIC PROJECTIVE REFERENCE")
    print("=" * 78)

    print(
        f"  rho_mod7={rho7}"
    )

    print(
        f"  k1={k1}"
    )

    print(
        f"  PRINCIPAL_BASE={PRINCIPAL_BASE}"
    )

    print(
        f"  (729-1)/7={normalized_729_increment}"
    )

    print(
        f"  ((729-1)/7)_mod7="
        f"{normalized_729_increment % 7}"
    )

    print(
        f"  predicted_slope="
        f"{predicted_slope}"
    )

    print(
        f"  projective_slope_matches_known_n6="
        f"{projective_slope_exact}"
    )

    # ------------------------------------------------------------------
    # 10. NORMALIZATION BY p,q,n
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. NORMALIZATION BY p, q, AND n")
    print("=" * 78)

    normalization_candidates = {
        "q1/p": Fraction(
            Q1,
            P,
        ),
        "q1/q": Fraction(
            Q1,
            Q,
        ),
        "q1/n": Fraction(
            Q1,
            N,
        ),
        "q3/p": Fraction(
            Q3,
            P,
        ),
        "q3/q": Fraction(
            Q3,
            Q,
        ),
        "q3/n": Fraction(
            Q3,
            N,
        ),
    }

    for name, value in normalization_candidates.items():

        print(
            f"  {name}={value} "
            f"integer={value.denominator == 1}"
        )

    # ------------------------------------------------------------------
    # 11. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("11. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The present experiment asks whether the exact terminal source pair
contains an obvious low-complexity signature of the prime parameters
p=2 and q=3.

The Schur row itself is exactly

    s0 = q1 - 2q3,
    s1 = q1 - 6q3,
    s2 = q1 - 18q3.

Consequently

    s1-s0 = -4q3,
    s2-s1 = -12q3,

and q1,q3 can be recovered exactly from the Schur coefficients.

The unresolved family problem is upstream:

    q1 = q1(p,q),
    q3 = q3(p,q).

This script deliberately does not guess such formulas.

A positive low-complexity identity involving p,q would be a useful clue
for the eventual symbolic n=pq theorem.

A negative result simply means that the source formula has to be
searched at a higher structural level.
"""
    )

    # ------------------------------------------------------------------
    # 12. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        source_recovery_exact
        and geometry_exact
        and geometric_difference_law
        and projective_slope_exact
    )

    print()
    print("=" * 78)
    print("12. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  source_recovery_exact="
        f"{source_recovery_exact}"
    )

    print(
        f"  Schur_difference_geometry_exact="
        f"{geometry_exact}"
    )

    print(
        f"  geometric_difference_law="
        f"{geometric_difference_law}"
    )

    print(
        f"  projective_slope_formula_exact="
        f"{projective_slope_exact}"
    )

    print(
        "  n=pq_family_formula_established=False"
    )

    print(
        "  reason=only the n=6 source pair is currently available"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 221R COMPLETE")


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