#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 224 — EXACT NORMALIZED SOURCE-COORDINATE RELATION AUDIT
==============================================================================

Known case:

    p = 2
    q = 3
    n = 6

    q1 = 29144191
    q3 = 24794967

Normalize:

    q3 = q*u

and then

    q1 = 1 + p*q*a
    u  = 1 + p*q*b.

For the current data:

    a = (q1-1)/(pq)
    b = (u-1)/(pq).

Experiment 224 asks whether a and b have simple exact arithmetic
relations.

It tests:

    a-b
    a+b
    a-q*b
    a-p*b
    a-(p+q)*b
    a-n*b
    q*a-p*b
    p*a-q*b

and exact ratios.

It also examines their p-adic and q-adic valuations and their
congruence profiles.

No family theorem is claimed.
Only main.py is used.
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


# ============================================================================
# EXACT HELPERS
# ============================================================================

def valuation_p(
    x: int,
    p: int,
) -> int | None:

    x = int(x)
    p = int(p)

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
        a, b = b, a % b

    return int(a)


def exact_ratio(
    a: int,
    b: int,
) -> Fraction:

    return Fraction(
        int(a),
        int(b),
    )


def residues(
    x: int,
    p: int,
    levels: int,
):

    result = []

    for e in range(
        1,
        levels + 1,
    ):

        modulus = int(
            p ** e
        )

        result.append(
            (
                modulus,
                int(x % modulus),
            )
        )

    return result


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 224 — EXACT NORMALIZED "
        "SOURCE-COORDINATE RELATION AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE NORMALIZATION
    # ------------------------------------------------------------------

    if Q3 % Q != 0:
        raise ArithmeticError(
            "q3 is not divisible by q."
        )

    U = int(
        Q3 // Q
    )

    pq = int(
        P * Q
    )

    if (Q1 - 1) % pq != 0:
        raise ArithmeticError(
            "q1-1 is not divisible by pq."
        )

    if (U - 1) % pq != 0:
        raise ArithmeticError(
            "u-1 is not divisible by pq."
        )

    A = int(
        (Q1 - 1) // pq
    )

    B = int(
        (U - 1) // pq
    )

    print()
    print("=" * 78)
    print("1. NORMALIZED SOURCE COORDINATES")
    print("=" * 78)

    print(
        f"  p={P}"
    )

    print(
        f"  q={Q}"
    )

    print(
        f"  pq={pq}"
    )

    print(
        f"  q3=q*u: u={U}"
    )

    print(
        f"  q1=1+pq*a: a={A}"
    )

    print(
        f"  u=1+pq*b: b={B}"
    )

    print(
        f"  q1_reconstruction="
        f"{1 + pq*A}"
    )

    print(
        f"  u_reconstruction="
        f"{1 + pq*B}"
    )

    print(
        f"  q1_exact="
        f"{1 + pq*A == Q1}"
    )

    print(
        f"  u_exact="
        f"{1 + pq*B == U}"
    )

    # ------------------------------------------------------------------
    # 2. a,b PRIME-ADIC PROFILE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. NORMALIZED COORDINATE VALUATIONS")
    print("=" * 78)

    for name, value in (
        ("a", A),
        ("b", B),
        ("a-b", A-B),
        ("a+b", A+B),
    ):

        print(
            f"  {name}={value} "
            f"v_p={valuation_p(value,P)} "
            f"v_q={valuation_p(value,Q)}"
        )

    # ------------------------------------------------------------------
    # 3. EXACT LINEAR RELATION SEARCH
    # ------------------------------------------------------------------

    combinations = {
        "a-b": A - B,
        "a+b": A + B,
        "a-q*b": A - Q * B,
        "a-p*b": A - P * B,
        "a-(p+q)*b": A - (P+Q) * B,
        "a-n*b": A - N * B,
        "q*a-p*b": Q * A - P * B,
        "p*a-q*b": P * A - Q * B,
        "q*a+b": Q * A + B,
        "p*a+b": P * A + B,
        "a+q*b": A + Q * B,
        "a+p*b": A + P * B,
    }

    print()
    print("=" * 78)
    print("3. EXACT LOW-COMPLEXITY RELATION SEARCH")
    print("=" * 78)

    for name, value in combinations.items():

        print(
            f"  {name}={value} "
            f"gcd_with_a={gcd_int(value,A)} "
            f"gcd_with_b={gcd_int(value,B)}"
        )

    # ------------------------------------------------------------------
    # 4. EXACT RATIOS
    # ------------------------------------------------------------------

    ratios = {
        "a/b": exact_ratio(A,B),
        "a/(q*b)": exact_ratio(A,Q*B),
        "a/(p*b)": exact_ratio(A,P*B),
        "a/((p+q)*b)": exact_ratio(
            A,
            (P+Q)*B,
        ),
        "a/(n*b)": exact_ratio(
            A,
            N*B,
        ),
    }

    print()
    print("=" * 78)
    print("4. EXACT NORMALIZED RATIOS")
    print("=" * 78)

    for name, value in ratios.items():

        print(
            f"  {name}={value}"
        )

    # ------------------------------------------------------------------
    # 5. FIRST-ORDER MOD p,q PROFILES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. NORMALIZED COORDINATES MOD p AND q")
    print("=" * 78)

    print(
        f"  a mod p={A % P}"
    )

    print(
        f"  a mod q={A % Q}"
    )

    print(
        f"  b mod p={B % P}"
    )

    print(
        f"  b mod q={B % Q}"
    )

    print(
        f"  a/b mod q="
        f"{(A * pow(B, -1, Q)) % Q if gcd(B,Q)==1 else 'undefined'}"
    )

    print(
        f"  a/b mod p="
        f"{(A * pow(B, -1, P)) % P if gcd(B,P)==1 else 'undefined'}"
    )

    # ------------------------------------------------------------------
    # 6. HIGHER p-ADIC PROFILE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. HIGHER p-ADIC PROFILES")
    print("=" * 78)

    print(
        "  a:"
    )

    for modulus, residue in residues(
        A,
        P,
        8,
    ):

        print(
            f"    mod {modulus}: {residue}"
        )

    print(
        "  b:"
    )

    for modulus, residue in residues(
        B,
        P,
        8,
    ):

        print(
            f"    mod {modulus}: {residue}"
        )

    # ------------------------------------------------------------------
    # 7. HIGHER q-ADIC PROFILE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. HIGHER q-ADIC PROFILES")
    print("=" * 78)

    print(
        "  a:"
    )

    for modulus, residue in residues(
        A,
        Q,
        8,
    ):

        print(
            f"    mod {modulus}: {residue}"
        )

    print(
        "  b:"
    )

    for modulus, residue in residues(
        B,
        Q,
        8,
    ):

        print(
            f"    mod {modulus}: {residue}"
        )

    # ------------------------------------------------------------------
    # 8. MOD pq PROFILE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. MOD pq PROFILE")
    print("=" * 78)

    print(
        f"  a mod pq={A % pq}"
    )

    print(
        f"  b mod pq={B % pq}"
    )

    print(
        f"  a-b mod pq={(A-B) % pq}"
    )

    print(
        f"  a+b mod pq={(A+B) % pq}"
    )

    # ------------------------------------------------------------------
    # 9. RECONSTRUCT q1,q3 DIRECTLY FROM a,b
    # ------------------------------------------------------------------

    q1_recovered = int(
        1 + pq * A
    )

    u_recovered = int(
        1 + pq * B
    )

    q3_recovered = int(
        Q * u_recovered
    )

    print()
    print("=" * 78)
    print("9. SOURCE RECONSTRUCTION FROM (a,b)")
    print("=" * 78)

    print(
        f"  q1_recovered={q1_recovered}"
    )

    print(
        f"  q1_original={Q1}"
    )

    print(
        f"  q3_recovered={q3_recovered}"
    )

    print(
        f"  q3_original={Q3}"
    )

    print(
        f"  q1_exact={q1_recovered == Q1}"
    )

    print(
        f"  q3_exact={q3_recovered == Q3}"
    )

    # ------------------------------------------------------------------
    # 10. SCHUR ROW IN NORMALIZED COORDINATES
    # ------------------------------------------------------------------

    s0 = int(
        Q1 - 2 * Q3
    )

    s1 = int(
        Q1 - 6 * Q3
    )

    s2 = int(
        Q1 - 18 * Q3
    )

    reconstructed_s = [
        int(
            1
            + pq*A
            - 2*Q*(1 + pq*B)
        ),
        int(
            1
            + pq*A
            - 6*Q*(1 + pq*B)
        ),
        int(
            1
            + pq*A
            - 18*Q*(1 + pq*B)
        ),
    ]

    schur_exact = (
        reconstructed_s
        == [s0,s1,s2]
    )

    print()
    print("=" * 78)
    print("10. SCHUR ROW FROM NORMALIZED COORDINATES")
    print("=" * 78)

    print(
        f"  reconstructed={reconstructed_s}"
    )

    print(
        f"  original={[s0,s1,s2]}"
    )

    print(
        f"  exact={schur_exact}"
    )

    # ------------------------------------------------------------------
    # 11. CANDIDATE RELATIONS THAT ACTUALLY HIT
    # ------------------------------------------------------------------

    exact_relations = []

    for name, value in combinations.items():

        if value == 0:
            exact_relations.append(
                name
            )

    ratio_one_relations = []

    for name, value in ratios.items():

        if value == 1:
            ratio_one_relations.append(
                name
            )

    print()
    print("=" * 78)
    print("11. EXACT IDENTITIES FOUND")
    print("=" * 78)

    print(
        f"  zero_relations={exact_relations}"
    )

    print(
        f"  unit_ratio_relations={ratio_one_relations}"
    )

    # ------------------------------------------------------------------
    # 12. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("12. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 223R established

    q3 = q*u,

    q1 = 1 + pq*a,

    u  = 1 + pq*b.

Experiment 224 asks whether the two normalized coordinates a and b
are themselves related by a simple expression involving p and q.

This is the next natural source-level question because a symbolic
n=pq theorem would ideally reduce the large terminal integers to a
small number of normalized quantities.

The important possible outcomes are:

    exact linear relation:
        a = c1*b + c2;

    exact multiplicative relation:
        a = c*b;

    simple modular relation:
        a = f(b) mod p or mod q;

    or no low-complexity relation.

A negative result is still useful: it tells us that the source layer
cannot be compressed to an elementary relation between these two
normalized coordinates.

A positive result would provide a strong clue toward the eventual
closed form for q1(p,q) and q3(p,q).

This remains a finite observation for n=6.
"""
    )

    # ------------------------------------------------------------------
    # 13. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        q1_recovered == Q1
        and q3_recovered == Q3
        and schur_exact
    )

    print()
    print("=" * 78)
    print("13. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  normalized_coordinates_exact="
        f"{q1_recovered == Q1 and u_recovered == U if 'U' in globals() else True}"
    )

    print(
        f"  q1_reconstruction_exact="
        f"{q1_recovered == Q1}"
    )

    print(
        f"  q3_reconstruction_exact="
        f"{q3_recovered == Q3}"
    )

    print(
        f"  Schur_normalized_reconstruction_exact="
        f"{schur_exact}"
    )

    print(
        f"  exact_zero_relations="
        f"{exact_relations}"
    )

    print(
        f"  exact_unit_ratios="
        f"{ratio_one_relations}"
    )

    print(
        "  general_n_pq_formula_established=False"
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
    print("EXPERIMENT 224 COMPLETE")


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

