#!/usr/bin/env python3

import math
import random
import time


# ========================================================================
# START EXPERIMENT 129
# CRT information threshold for unique factor identification
# ========================================================================

print("=" * 72)
print("START EXPERIMENT 129")
print("CRT information threshold for unique factor identification")
print("=" * 72)


# ------------------------------------------------------------------------
# RADIX POOL
#
# Pairwise coprime primes.
#
# We gradually add radices and observe how the product M grows.
# ------------------------------------------------------------------------

RADIX_POOL = [
    3,
    5,
    7,
    11,
    13,
    17,
    19,
    23,
    29,
    31,
    37,
    41,
    43,
]

print()
print("RADIX POOL")
print(RADIX_POOL)


# ------------------------------------------------------------------------
# MILLER-RABIN
# ------------------------------------------------------------------------

def is_probable_prime(n):

    if n < 2:
        return False

    small = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    )

    for p in small:

        if n == p:
            return True

        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        d //= 2
        s += 1

    bases = (
        2,
        325,
        9375,
        28178,
        450775,
        9780504,
        1795265022,
    )

    for a in bases:

        if a % n == 0:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):

            x = (x * x) % n

            if x == n - 1:
                break

        else:
            return False

    return True


# ------------------------------------------------------------------------
# RANDOM PRIME
# ------------------------------------------------------------------------

def random_prime(bits):

    while True:

        x = random.getrandbits(bits)

        x |= 1
        x |= 1 << (bits - 1)

        if is_probable_prime(x):
            return x


# ------------------------------------------------------------------------
# SEMIPRIME
# ------------------------------------------------------------------------

def random_semiprime(bits):

    pb = bits // 2
    qb = bits - pb

    while True:

        p = random_prime(pb)
        q = random_prime(qb)

        if p == q:
            continue

        p, q = sorted((p, q))

        n = p * q

        if n.bit_length() == bits:
            return n, p, q


# ------------------------------------------------------------------------
# CRT
# ------------------------------------------------------------------------

def crt(residues, moduli):

    x = 0
    M = 1

    for a, m in zip(residues, moduli):

        inv = pow(M, -1, m)

        t = ((a - x) * inv) % m

        x += M * t
        M *= m
        x %= M

    return x, M


# ------------------------------------------------------------------------
# NUMBER OF VALUES IN A RESIDUE CLASS
#
# Count:
#
#     p = A + kM
#
# satisfying
#
#     2 <= p <= sqrt(n).
# ------------------------------------------------------------------------

def residue_class_count(n, A, M):

    limit = math.isqrt(n)

    if A == 0:

        first = M

    else:

        first = A

        if first < 2:

            first += (
                (2 - first + M - 1)
                // M
            ) * M

    if first > limit:

        return 0

    return (
        (limit - first) // M
    ) + 1


# ------------------------------------------------------------------------
# VERIFY THE TRUE FACTOR IS IN THE CLASS
# ------------------------------------------------------------------------

def true_class(n, p, radices):

    residues = [
        p % r
        for r in radices
    ]

    A, M = crt(
        residues,
        radices
    )

    count = residue_class_count(
        n,
        A,
        M
    )

    return {
        "residues": residues,
        "A": A,
        "M": M,
        "count": count,
    }


# ------------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------------

random.seed(129)

TEST_BITS = [
    30,
    36,
    42,
    48,
    54,
]

for bits in TEST_BITS:

    print()
    print("=" * 72)
    print("GENERATING", bits, "-BIT SEMIPRIME")
    print("=" * 72)

    n, p, q = random_semiprime(bits)

    sqrt_n = math.isqrt(n)

    print()
    print("-" * 72)

    print("n bits =", n.bit_length())
    print("n      =", n)
    print("p      =", p)
    print("q      =", q)

    print()
    print("sqrt(n) =", sqrt_n)

    print()
    print(
        "THEORETICAL UNIQUENESS THRESHOLD:"
    )

    print(
        "M > sqrt(n) =",
        sqrt_n
    )

    print()
    print("RADIX GROWTH")

    M = 1

    for index, r in enumerate(RADIX_POOL, start=1):

        M *= r

        info = true_class(
            n,
            p,
            RADIX_POOL[:index]
        )

        count = info["count"]

        unique = (
            count <= 1
        )

        larger_than_sqrt = (
            M > sqrt_n
        )

        print()
        print(
            f"{index:2d} radices"
        )

        print(
            "radix added =",
            r
        )

        print(
            "M =",
            M
        )

        print(
            "M / sqrt(n) =",
            f"{M / sqrt_n:.6f}"
        )

        print(
            "possible p values in true residue class =",
            count
        )

        print(
            "M > sqrt(n) =",
            larger_than_sqrt
        )

        print(
            "unique p =",
            unique
        )

        if unique:

            print(
                ">>> UNIQUE FACTOR RESIDUE CLASS <<<"
            )

            break

    # ------------------------------------------------------------
    # Benchmark reconstruction.
    # ------------------------------------------------------------

    print()
    print("DIRECT RECONSTRUCTION TEST")

    # Find the minimum radix prefix producing a unique p.
    M = 1
    threshold_index = None
    threshold_info = None

    for index, r in enumerate(
        RADIX_POOL,
        start=1
    ):

        M *= r

        info = true_class(
            n,
            p,
            RADIX_POOL[:index]
        )

        if info["count"] == 1:

            threshold_index = index
            threshold_info = info

            break

    if threshold_index is None:

        print(
            "No unique residue class within radix pool."
        )

    else:

        print(
            "minimum number of radices =",
            threshold_index
        )

        print(
            "radices =",
            RADIX_POOL[:threshold_index]
        )

        print(
            "M =",
            threshold_info["M"]
        )

        print(
            "A =",
            threshold_info["A"]
        )

        print(
            "true p =",
            p
        )

        recovered_candidates = []

        A = threshold_info["A"]
        M = threshold_info["M"]

        limit = math.isqrt(n)

        x = 0

        while True:

            candidate = A + x * M

            if candidate > limit:
                break

            if candidate >= 2:

                recovered_candidates.append(
                    candidate
                )

            x += 1

        print(
            "candidate count =",
            len(recovered_candidates)
        )

        print(
            "candidates =",
            recovered_candidates
        )

        print(
            "TRUE p recovered =",
            p in recovered_candidates
        )

    # ------------------------------------------------------------
    # INFORMATION MEASURE
    #
    # Number of bits represented by M:
    #
    #     log2(M)
    #
    # Compare this with:
    #
    #     log2(sqrt(n))
    #       = log2(n)/2.
    # ------------------------------------------------------------

    print()
    print("INFORMATION")

    M = 1

    for index, r in enumerate(
        RADIX_POOL,
        start=1
    ):

        M *= r

        bits_M = math.log2(M)

        bits_needed = (
            math.log2(n) / 2
        )

        print(
            f"{index:2d} radices: "
            f"log2(M)={bits_M:8.3f} "
            f"needed={bits_needed:8.3f} "
            f"ratio={M / sqrt_n:12.6f}"
        )


# ========================================================================
# FINISHED EXPERIMENT 129
# ========================================================================

print()
print("=" * 72)
print("FINISHED EXPERIMENT 129")
print("=" * 72)
