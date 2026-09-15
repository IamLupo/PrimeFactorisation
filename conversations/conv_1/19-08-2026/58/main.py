#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 225 — EXACT NORMALIZED SOURCE FACTOR / PRIME-SPECTRUM PROVENANCE
==============================================================================

Known case:

    p = 2
    q = 3
    n = 6

    q1 = 29144191
    q3 = 24794967

Normalize:

    q3 = q*u

    q1 = 1 + pq*a
    u  = 1 + pq*b.

Experiment 224 found no obvious exact relation between a and b.

Experiment 225 therefore studies their exact integer factor/content
structure.

The experiment performs:

    * complete trial-division factorization of a,b;
    * factorization of q1,q3 and terminal Schur coefficients;
    * p/q/n valuations;
    * gcd intersections between the normalized coordinates;
    * prime-spectrum overlap;
    * exact reconstruction after removing p,q,n factors;
    * comparison of normalized prime factors with the Schur row;
    * tests for simple factor inheritance.

This is descriptive for n=6 only.

Only main.py is used.
"""

from __future__ import annotations

import sys
from math import gcd


# ============================================================================
# KNOWN INSTANCE
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


def factor_integer(
    x: int,
) -> dict[int, int]:

    x = abs(int(x))

    if x < 2:
        return {}

    factors: dict[int, int] = {}

    while x % 2 == 0:
        factors[2] = factors.get(2, 0) + 1
        x //= 2

    d = 3

    while d * d <= x:

        while x % d == 0:
            factors[d] = factors.get(d, 0) + 1
            x //= d

        d += 2

    if x > 1:
        factors[x] = factors.get(x, 0) + 1

    return factors


def factor_string(
    factors: dict[int, int],
) -> str:

    if not factors:
        return "1"

    parts = []

    for p, e in sorted(
        factors.items()
    ):

        if e == 1:
            parts.append(
                str(p)
            )
        else:
            parts.append(
                f"{p}^{e}"
            )

    return " * ".join(parts)


def prime_spectrum(
    factors: dict[int, int],
) -> set[int]:

    return set(
        factors.keys()
    )


def intersection(
    a: set[int],
    b: set[int],
) -> list[int]:

    return sorted(
        a & b
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 225 — EXACT NORMALIZED SOURCE FACTOR / "
        "PRIME-SPECTRUM PROVENANCE"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. NORMALIZED SOURCE COORDINATES
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
        f"  n={N}"
    )

    print(
        f"  q1={Q1}"
    )

    print(
        f"  q3={Q3}"
    )

    print(
        f"  u=q3/q={U}"
    )

    print(
        f"  a=(q1-1)/(pq)={A}"
    )

    print(
        f"  b=(u-1)/(pq)={B}"
    )

    # ------------------------------------------------------------------
    # 2. FACTORIZATION OF CORE SOURCE VALUES
    # ------------------------------------------------------------------

    values = {
        "q1": Q1,
        "q3": Q3,
        "u": U,
        "a": A,
        "b": B,
        "q": Q,
        "p": P,
        "n": N,
    }

    print()
    print("=" * 78)
    print("2. EXACT FACTORIZATION")
    print("=" * 78)

    factor_map: dict[str, dict[int, int]] = {}

    for name, value in values.items():

        factors = factor_integer(
            value
        )

        factor_map[name] = factors

        print(
            f"  {name}={value}"
        )

        print(
            f"    factorization="
            f"{factor_string(factors)}"
        )

    # ------------------------------------------------------------------
    # 3. NORMALIZED COORDINATE PRIME VALUATIONS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. p/q/n VALUATION PROFILE")
    print("=" * 78)

    for name in (
        "q1",
        "q3",
        "u",
        "a",
        "b",
    ):

        value = values[name]

        print(
            f"  {name}: "
            f"v_p={valuation_p(value,P)} "
            f"v_q={valuation_p(value,Q)} "
            f"v_n_component_gcd="
            f"{gcd_int(value,N)}"
        )

    # ------------------------------------------------------------------
    # 4. PRIME SPECTRA
    # ------------------------------------------------------------------

    spectra = {
        name: prime_spectrum(
            factor_map[name]
        )
        for name in values
    }

    print()
    print("=" * 78)
    print("4. PRIME SPECTRA")
    print("=" * 78)

    for name in (
        "q1",
        "q3",
        "u",
        "a",
        "b",
    ):

        print(
            f"  {name}: "
            f"{sorted(spectra[name])}"
        )

    # ------------------------------------------------------------------
    # 5. SPECTRUM OVERLAPS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. PRIME-SPECTRUM OVERLAPS")
    print("=" * 78)

    pairs = [
        ("a", "b"),
        ("a", "q1"),
        ("a", "q3"),
        ("b", "q1"),
        ("b", "q3"),
        ("u", "q1"),
        ("u", "q3"),
        ("a", "n"),
        ("b", "n"),
    ]

    for left, right in pairs:

        common = intersection(
            spectra[left],
            spectra[right],
        )

        print(
            f"  {left} vs {right}: "
            f"common_primes={common}"
        )

    # ------------------------------------------------------------------
    # 6. EXACT GCD NETWORK
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. EXACT GCD NETWORK")
    print("=" * 78)

    gcd_pairs = [
        ("a", "b"),
        ("a", "u"),
        ("a", "q1"),
        ("a", "q3"),
        ("b", "u"),
        ("b", "q1"),
        ("b", "q3"),
        ("u", "q1"),
        ("u", "q3"),
    ]

    for left, right in gcd_pairs:

        value = gcd_int(
            values[left],
            values[right],
        )

        print(
            f"  gcd({left},{right})={value}"
        )

    # ------------------------------------------------------------------
    # 7. REMOVE p AND q CONTENT
    # ------------------------------------------------------------------

    def remove_factor(
        x: int,
        prime: int,
    ) -> int:

        x = abs(int(x))

        while x % prime == 0:
            x //= prime

        return int(x)

    a_primitive = remove_factor(
        remove_factor(
            A,
            P,
        ),
        Q,
    )

    b_primitive = remove_factor(
        remove_factor(
            B,
            P,
        ),
        Q,
    )

    u_primitive = remove_factor(
        remove_factor(
            U,
            P,
        ),
        Q,
    )

    print()
    print("=" * 78)
    print("7. p/q-PRIMITIVE NORMALIZATION")
    print("=" * 78)

    print(
        f"  a_primitive={a_primitive}"
    )

    print(
        f"  b_primitive={b_primitive}"
    )

    print(
        f"  u_primitive={u_primitive}"
    )

    print(
        f"  gcd(a_primitive,b_primitive)="
        f"{gcd_int(a_primitive,b_primitive)}"
    )

    # ------------------------------------------------------------------
    # 8. SCHUR ROW
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

    schur = {
        "s0": s0,
        "s1": s1,
        "s2": s2,
    }

    print()
    print("=" * 78)
    print("8. SCHUR COEFFICIENT FACTORIZATION")
    print("=" * 78)

    schur_factors: dict[str, dict[int,int]] = {}

    for name, value in schur.items():

        factors = factor_integer(
            value
        )

        schur_factors[name] = factors

        print(
            f"  {name}={value}"
        )

        print(
            f"    factorization="
            f"{factor_string(factors)}"
        )

        print(
            f"    v_p={valuation_p(value,P)} "
            f"v_q={valuation_p(value,Q)}"
        )

    # ------------------------------------------------------------------
    # 9. NORMALIZED/SCHUR SPECTRUM COMPARISON
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. NORMALIZED SOURCE VS SCHUR PRIME SPECTRA")
    print("=" * 78)

    for source_name in (
        "a",
        "b",
        "u",
    ):

        for schur_name in (
            "s0",
            "s1",
            "s2",
        ):

            common = intersection(
                spectra[source_name],
                prime_spectrum(
                    schur_factors[schur_name]
                ),
            )

            print(
                f"  {source_name} vs {schur_name}: "
                f"common_primes={common}"
            )

    # ------------------------------------------------------------------
    # 10. FACTOR CONTENT EXPLANATION TESTS
    # ------------------------------------------------------------------

    # Known explicit Schur factors:
    # q3 = q*u
    # s0 = q1 - 2q*u
    # s1 = q1 - 6q*u
    # s2 = q1 - 18q*u

    schur_prime_union = set()

    for factors in schur_factors.values():
        schur_prime_union.update(
            factors.keys()
        )

    a_unique_primes = sorted(
        spectra["a"]
        - spectra["q1"]
        - spectra["q3"]
    )

    b_unique_primes = sorted(
        spectra["b"]
        - spectra["q1"]
        - spectra["q3"]
    )

    print()
    print("=" * 78)
    print("10. NEW PRIME CONTENT TEST")
    print("=" * 78)

    print(
        f"  primes_unique_to_a_vs_sources="
        f"{a_unique_primes}"
    )

    print(
        f"  primes_unique_to_b_vs_sources="
        f"{b_unique_primes}"
    )

    print(
        f"  a_spectrum_subset_of_schur_union="
        f"{spectra['a'].issubset(schur_prime_union)}"
    )

    print(
        f"  b_spectrum_subset_of_schur_union="
        f"{spectra['b'].issubset(schur_prime_union)}"
    )

    # ------------------------------------------------------------------
    # 11. SIMPLE FACTORIZATION RELATION SEARCH
    # ------------------------------------------------------------------

    relations = {
        "a = q*b": A == Q * B,
        "a = p*b": A == P * B,
        "a = n*b": A == N * B,
        "u = p*a": U == P * A,
        "u = q*a": U == Q * A,
        "b = p*q": B == P * Q,
        "b = q": B == Q,
        "b = p": B == P,
    }

    print()
    print("=" * 78)
    print("11. SIMPLE FACTORIZATION RELATIONS")
    print("=" * 78)

    for name, value in relations.items():

        print(
            f"  {name}={value}"
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
Experiment 224 found no simple exact linear or ratio relation between
the normalized coordinates

    a = (q1-1)/(pq),
    b = ((q3/q)-1)/(pq).

Experiment 225 now asks whether their arithmetic complexity is itself
structured.

There are three possibilities:

    1. a and b have large but related prime spectra;

    2. they are arithmetically independent-looking;

    3. a small set of new primes systematically controls one layer.

The comparison with the Schur coefficients asks whether the arithmetic
content of the normalized source coordinates reappears in the terminal
Schur row.

A negative result is valuable because it suggests that the source
formula is not naturally expressed through ordinary integer
factorization.

A positive overlap becomes a concrete candidate for the eventual
n=pq source theorem.
"""
    )

    # ------------------------------------------------------------------
    # 13. FINAL
    # ------------------------------------------------------------------

    factorization_exact = all(
        factor_integer(
            values[name]
        )
        for name in (
            "q1",
            "q3",
            "u",
            "a",
            "b",
        )
        if values[name] >= 2
    )

    reconstruction_exact = (
        Q * U == Q3
        and 1 + N * A == Q1
        and 1 + N * B == U
    )

    print()
    print("=" * 78)
    print("13. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  factorization_completed="
        f"{factorization_exact}"
    )

    print(
        f"  normalized_reconstruction_exact="
        f"{reconstruction_exact}"
    )

    print(
        f"  gcd_a_b="
        f"{gcd_int(A,B)}"
    )

    print(
        f"  a_prime_spectrum="
        f"{sorted(spectra['a'])}"
    )

    print(
        f"  b_prime_spectrum="
        f"{sorted(spectra['b'])}"
    )

    print(
        "  general_n_pq_theorem_established=False"
    )

    print(
        "  reason=single_known_n6_instance"
    )

    print(
        f"  failures={0 if (factorization_exact and reconstruction_exact) else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS="
        f"{factorization_exact and reconstruction_exact}"
    )

    print()
    print("EXPERIMENT 225 COMPLETE")


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

