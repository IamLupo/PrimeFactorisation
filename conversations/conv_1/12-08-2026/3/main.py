#!/usr/bin/env python3

"""
==============================================================================
KAPPA AUXILIARY -> u COLLAPSE EXPERIMENT
==============================================================================

Purpose
-------

Investigate whether auxiliary primes

    r, s, t, ...

can impose enough structure on

    A = F(p)F(q)

to reduce the possible value of

    u = p + q

for a known n = pq.

The central identities are:

    F(x) = x^2 - x + 1

    A = F(p)F(q)

    A = u^2 - (n+1)u + n^2 - n + 1

and therefore

    D = 4A - 3(n-1)^2
      = (2u - (n+1))^2.

Thus A determines u through a perfect-square discriminant.

IMPORTANT:

Exact K values are an oracle:

    1-K_r = A F(r)

so

    A = (1-K_r)/F(r).

That does NOT demonstrate that A can be computed from n.

This experiment therefore separates:

    MODE 1: exact-oracle information
    MODE 2: n-only structural information / collision analysis

No CSV files are written.
"""

from __future__ import annotations

import math
import random
from collections import defaultdict
from itertools import combinations


# ============================================================================
# BASIC FUNCTIONS
# ============================================================================

def F(x: int) -> int:
    return x * x - x + 1


def A_of_pq(p: int, q: int) -> int:
    return F(p) * F(q)


def u_of_pq(p: int, q: int) -> int:
    return p + q


def discriminant_from_A(n: int, A: int) -> int:
    return 4 * A - 3 * (n - 1) ** 2


def u_candidates_from_A(n: int, A: int):
    """
    Recover possible u from exact A.

        D = (2u-(n+1))^2

    Hence:

        u = ((n+1) +/- sqrt(D))/2

    Returns integer candidates only.
    """

    D = discriminant_from_A(n, A)

    if D < 0:
        return []

    w = math.isqrt(D)

    if w * w != D:
        return []

    out = []

    for sign in (-1, 1):
        numerator = n + 1 + sign * w

        if numerator % 2 == 0:
            u = numerator // 2

            if u >= 2:
                out.append(u)

    return sorted(set(out))


def factor_pair_from_u(n: int, u: int):
    """
    Given n and u=p+q, solve

        x^2-u*x+n=0.
    """

    D = u * u - 4 * n

    if D < 0:
        return None

    w = math.isqrt(D)

    if w * w != D:
        return None

    if (u + w) % 2:
        return None

    p = (u - w) // 2
    q = (u + w) // 2

    if p * q != n:
        return None

    return tuple(sorted((p, q)))


def primes_up_to(limit: int):
    sieve = [True] * (limit + 1)

    if limit >= 0:
        sieve[0] = False
    if limit >= 1:
        sieve[1] = False

    for i in range(2, math.isqrt(limit) + 1):
        if sieve[i]:
            for j in range(i * i, limit + 1, i):
                sieve[j] = False

    return [i for i, ok in enumerate(sieve) if ok]


# ============================================================================
# TEST PRIME GENERATION
# ============================================================================

PRIMES = primes_up_to(500)

AUXILIARY_SEQUENCES = {
    "first": PRIMES[:9],

    "odd_first": [
        p for p in PRIMES
        if p != 2
    ][:8],

    "every_other": PRIMES[::2][:7],

    "larger_first": [
        p for p in PRIMES
        if p >= 11
    ][:6],
}

MODULI = [
    3, 5, 7, 11, 13, 17, 19,
    23, 29, 31, 37, 41, 43, 47
]


# ============================================================================
# SEMIPRIME GENERATION
# ============================================================================

def is_prime(n: int) -> bool:
    if n < 2:
        return False

    if n % 2 == 0:
        return n == 2

    d = 3

    while d * d <= n:
        if n % d == 0:
            return False
        d += 2

    return True


def random_prime_with_bits(bits: int, rng: random.Random) -> int:
    """
    Generate a random odd prime with approximately 'bits' bits.
    """

    low = 1 << (bits - 1)
    high = (1 << bits) - 1

    while True:
        x = rng.randint(low, high)
        x |= 1

        if is_prime(x):
            return x


def generate_semiprime(bits: int, rng: random.Random, balanced=True):
    """
    Generate p*q with approximately 'bits' total bits.

    balanced=True:
        p and q have approximately bits/2 bits.

    balanced=False:
        deliberately allow more asymmetric factors.
    """

    if balanced:
        b1 = bits // 2
        b2 = bits - b1
    else:
        # Random split, but avoid extremely tiny factors.
        b1 = rng.randint(max(8, bits // 4), bits * 3 // 4)
        b2 = bits - b1

    p = random_prime_with_bits(b1, rng)
    q = random_prime_with_bits(b2, rng)

    return tuple(sorted((p, q)))


# ============================================================================
# AUXILIARY PRODUCT
# ============================================================================

def auxiliary_product(rs):
    B = 1

    for r in rs:
        B *= F(r)

    return B


def bits(x: int) -> int:
    return x.bit_length()


# ============================================================================
# EXACT ORACLE EXPERIMENT
# ============================================================================

def oracle_experiment(p: int, q: int, rs):
    """
    Pretend we know exact K for every auxiliary r.

        K_r = 1 - A F(r)

    Recover A and then u.
    """

    n = p * q
    A = A_of_pq(p, q)

    print()
    print("EXACT ORACLE")
    print("-" * 78)

    print(f"n = {n}")
    print(f"p = {p}")
    print(f"q = {q}")
    print(f"true u = {p + q}")
    print(f"A = {A}")
    print(f"A bits = {bits(A)}")

    recovered_values = []

    for r in rs:
        K = 1 - A * F(r)

        recovered_A = (1 - K) // F(r)

        recovered_values.append(recovered_A)

        print(
            f"r={r:3d}  "
            f"F(r)={F(r):8d}  "
            f"K bits={bits(abs(K)):3d}  "
            f"A recovered={recovered_A == A}"
        )

    if len(set(recovered_values)) == 1:
        print("PASS: all auxiliary oracle values recover the same A.")

    us = u_candidates_from_A(n, A)

    print()
    print(f"u candidates from exact A: {us}")

    if p + q in us:
        print("PASS: true u recovered.")

    factor_candidates = []

    for u in us:
        pair = factor_pair_from_u(n, u)

        if pair is not None:
            factor_candidates.append(pair)

    print(f"factor candidates: {factor_candidates}")


# ============================================================================
# DISCRIMINANT TEST
# ============================================================================

def discriminant_experiment(p: int, q: int):
    n = p * q
    u = p + q
    A = A_of_pq(p, q)

    D = discriminant_from_A(n, A)
    expected = (2 * u - (n + 1)) ** 2

    print()
    print("DISCRIMINANT STRUCTURE")
    print("-" * 78)

    print(f"n = {n}")
    print(f"u = {u}")
    print(f"A = {A}")
    print(f"D = {D}")
    print(f"(2u-(n+1))^2 = {expected}")
    print(f"identity holds = {D == expected}")


# ============================================================================
# A-RESIDUE -> POSSIBLE u TEST
# ============================================================================

def possible_u_from_A_residue(
    n: int,
    modulus: int,
    A_residue: int,
    search_limit=None,
):
    """
    Find u values satisfying the polynomial relation modulo modulus.

        A = u^2 -(n+1)u+n^2-n+1

    Therefore:

        u^2 -(n+1)u+n^2-n+1-A = 0 mod m.

    This is intentionally a residue-level test.

    It does NOT factor n.

    For demonstration we enumerate u modulo m, not up to sqrt(n).
    """

    if search_limit is None:
        search_limit = modulus

    survivors = []

    for u in range(search_limit):
        lhs = (
            u * u
            - (n + 1) * u
            + n * n
            - n
            + 1
            - A_residue
        )

        if lhs % modulus == 0:
            survivors.append(u)

    return survivors


# ============================================================================
# MODULAR u-CONSTRAINT EXPERIMENT
# ============================================================================

def modular_u_constraint_experiment(p: int, q: int, rs):
    """
    Use the TRUE A only to create a controlled residue experiment.

    The important quantity is not how many A bits we know,
    but how many u residues survive.

    This isolates whether modular information about A has
    strong consequences for u.
    """

    n = p * q
    A = A_of_pq(p, q)
    true_u = p + q

    print()
    print("MODULAR A -> u CONSTRAINT EXPERIMENT")
    print("-" * 78)

    print(f"n bits = {bits(n)}")
    print(f"A bits = {bits(A)}")
    print(f"true u = {true_u}")

    print()
    print(
        f"{'mod':>6}"
        f"{'A mod m':>14}"
        f"{'u residues':>14}"
        f"{'true u mod m':>16}"
        f"{'survivors':>12}"
    )
    print("-" * 78)

    for m in MODULI:
        a_mod = A % m

        survivors = possible_u_from_A_residue(
            n,
            m,
            a_mod,
            m,
        )

        true_residue = true_u % m

        print(
            f"{m:6d}"
            f"{a_mod:14d}"
            f"{','.join(map(str, survivors)):>14}"
            f"{true_residue:16d}"
            f"{len(survivors):12d}"
        )


# ============================================================================
# COMBINED MODULUS EXPERIMENT
# ============================================================================

def combined_modulus_u_experiment(p: int, q: int):
    """
    Combine small moduli using CRT.

    Instead of reporting only A information, determine how many
    u residues modulo M survive.

    Since M grows rapidly, stop before it becomes impractical.
    """

    n = p * q
    A = A_of_pq(p, q)
    true_u = p + q

    M = 1

    print()
    print("COMBINED MODULUS -> u SURVIVAL")
    print("-" * 78)

    print(
        f"{'m':>5}"
        f"{'M bits':>8}"
        f"{'u residues mod M':>18}"
        f"{'true u mod M':>18}"
    )
    print("-" * 78)

    for m in MODULI:
        if math.gcd(M, m) != 1:
            continue

        new_M = M * m

        # We cannot enumerate huge M.
        # Stop once residue enumeration becomes unreasonable.
        if new_M > 2_000_000:
            break

        M = new_M

        survivors = []

        A_mod = A % M

        for u in range(M):
            lhs = (
                u * u
                - (n + 1) * u
                + n * n
                - n
                + 1
                - A_mod
            )

            if lhs % M == 0:
                survivors.append(u)

        true_residue = true_u % M

        print(
            f"{m:5d}"
            f"{bits(M):8d}"
            f"{len(survivors):18d}"
            f"{true_residue:18d}"
        )


# ============================================================================
# FACTOR-PAIR COLLISION SEARCH
# ============================================================================

def factor_pairs(n: int):
    """
    Experimental divisor enumeration.

    Used only for small n in collision experiments.
    """

    result = []

    limit = math.isqrt(n)

    for d in range(2, limit + 1):
        if n % d == 0:
            q = n // d

            if d <= q:
                result.append((d, q))

    return result


def kappa_signature(p: int, q: int, rs):
    """
    Signature consisting of exact K values.

    This demonstrates an important fact:
    for fixed n, different factor pairs generally produce
    different K signatures.
    """

    A = A_of_pq(p, q)

    return tuple(
        1 - A * F(r)
        for r in rs
    )


def collision_experiment(max_n=10000):
    """
    Search small composite n having multiple factor pairs.

    For each n compare:

      A
      A mod m
      K signatures for auxiliary r

    This tells us whether auxiliary observations distinguish
    different factorizations.

    It does NOT claim those observations are computable from n.
    """

    print()
    print("FIXED-n FACTOR-PAIR COLLISION EXPERIMENT")
    print("-" * 78)

    rs = AUXILIARY_SEQUENCES["first"]

    interesting = []

    for n in range(4, max_n + 1):
        pairs = factor_pairs(n)

        if len(pairs) < 2:
            continue

        A_values = {
            A_of_pq(p, q)
            for p, q in pairs
        }

        signatures = {
            kappa_signature(p, q, rs)
            for p, q in pairs
        }

        if len(A_values) < len(pairs):
            interesting.append(
                (
                    n,
                    pairs,
                    A_values,
                    signatures,
                )
            )

    print(f"range checked: 4 .. {max_n}")
    print(f"n with >=2 factor pairs: ", end="")

    count = sum(
        1 for n in range(4, max_n + 1)
        if len(factor_pairs(n)) >= 2
    )

    print(count)

    print(f"A collisions found: {len(interesting)}")

    if interesting:
        print()
        for n, pairs, A_values, signatures in interesting[:20]:
            print(f"n={n}")
            print(f"  pairs      = {pairs}")
            print(f"  A values   = {sorted(A_values)}")
            print(f"  signatures = {len(signatures)}")

    print()
    print(
        "Interpretation: different factor pairs can be distinguished "
        "by exact auxiliary K values, but this remains oracle information."
    )


# ============================================================================
# AUXILIARY PRODUCT INFORMATION
# ============================================================================

def auxiliary_growth():
    print()
    print("AUXILIARY PRODUCT GROWTH")
    print("-" * 78)

    print(
        f"{'sequence':<18}"
        f"{'count':>8}"
        f"{'B bits':>10}"
        f"{'B':>25}"
    )

    print("-" * 78)

    for name, rs in AUXILIARY_SEQUENCES.items():
        B = auxiliary_product(rs)

        print(
            f"{name:<18}"
            f"{len(rs):8d}"
            f"{bits(B):10d}"
            f"{B:25d}"
        )

        print(f"  r = {tuple(rs)}")


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA AUXILIARY -> u COLLAPSE EXPERIMENT")
    print("=" * 78)

    print()
    print("Core identities:")
    print()
    print("    A = F(p)F(q)")
    print()
    print("    A = u^2 -(n+1)u+n^2-n+1")
    print()
    print("    D = 4A - 3(n-1)^2")
    print()
    print("    D = (2u-(n+1))^2")
    print()

    print("Auxiliary sequences:")
    auxiliary_growth()

    # ------------------------------------------------------------------
    # Representative 50-bit target
    # ------------------------------------------------------------------

    rng = random.Random(20260812)

    p, q = generate_semiprime(
        50,
        rng,
        balanced=True,
    )

    n = p * q
    A = A_of_pq(p, q)

    print()
    print("REPRESENTATIVE TARGET")
    print("-" * 78)

    print(f"p = {p}")
    print(f"q = {q}")
    print(f"n = {n}")
    print(f"n bits = {bits(n)}")
    print(f"A = {A}")
    print(f"A bits = {bits(A)}")
    print(f"u = {p + q}")

    # ------------------------------------------------------------------
    # Exact oracle
    # ------------------------------------------------------------------

    oracle_experiment(
        p,
        q,
        AUXILIARY_SEQUENCES["first"],
    )

    # ------------------------------------------------------------------
    # Exact discriminant
    # ------------------------------------------------------------------

    discriminant_experiment(p, q)

    # ------------------------------------------------------------------
    # Individual modular constraints
    # ------------------------------------------------------------------

    modular_u_constraint_experiment(
        p,
        q,
        AUXILIARY_SEQUENCES["first"],
    )

    # ------------------------------------------------------------------
    # Combined modular constraints
    # ------------------------------------------------------------------

    combined_modulus_u_experiment(
        p,
        q,
    )

    # ------------------------------------------------------------------
    # Repeat on several targets
    # ------------------------------------------------------------------

    print()
    print("MULTI-TARGET SUMMARY")
    print("-" * 78)

    print(
        f"{'target':>7}"
        f"{'n bits':>8}"
        f"{'A bits':>8}"
        f"{'u':>18}"
    )

    for i in range(8):

        p, q = generate_semiprime(
            50,
            rng,
            balanced=(i < 4),
        )

        n = p * q
        A = A_of_pq(p, q)

        print(
            f"{i+1:7d}"
            f"{bits(n):8d}"
            f"{bits(A):8d}"
            f"{p+q:18d}"
        )

    # ------------------------------------------------------------------
    # Small-n collision experiment
    # ------------------------------------------------------------------

    collision_experiment(
        max_n=5000,
    )

    # ------------------------------------------------------------------
    # Final interpretation
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL INTERPRETATION")
    print("=" * 78)

    print()
    print("1. EXACT ORACLE")
    print("-" * 78)
    print(
        "Exact K_r plus known r gives A immediately:"
    )
    print()
    print(
        "    A = (1-K_r)/F(r)"
    )

    print()
    print("2. A -> u")
    print("-" * 78)
    print(
        "Once A is known, the discriminant condition forces u=p+q."
    )

    print()
    print("3. MODULAR EXPERIMENT")
    print("-" * 78)
    print(
        "The useful quantity is not merely the number of bits known "
        "about A."
    )

    print()
    print(
        "Instead we count how many u residues survive:"
    )

    print()
    print(
        "    A mod M"
    )
    print(
        "       ->"
    )
    print(
        "    u^2 -(n+1)u+n^2-n+1 = A mod M"
    )

    print()
    print("4. CRITICAL LIMITATION")
    print("-" * 78)
    print(
        "If the K_r values cannot themselves be computed efficiently "
        "from n, then their information cannot be counted as an "
        "n-only factoring method."
    )

    print()
    print("5. NEXT TARGET")
    print("-" * 78)
    print(
        "The strongest future experiment is therefore:"
    )

    print()
    print(
        "    n"
    )
    print(
        "     -> auxiliary structure"
    )
    print(
        "     -> constraints on A"
    )
    print(
        "     -> surviving u"
    )
    print(
        "     -> factor pair"
    )

    print()
    print(
        "without supplying p, q, A, or K as an oracle."
    )

    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)


if __name__ == "__main__":
    main()
