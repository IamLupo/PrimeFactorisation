#!/usr/bin/env python3

"""
================================================================================
FAST MULTI-MODULUS FACTOR-BRANCH DISCRIMINATION EXPERIMENT
================================================================================

Original hypothesis:

For odd primes p,q,

    n = p*q.

Modulo 4:

    n == 1:
        (p,q) = (1,1) or (3,3)

    n == 3:
        (p,q) = (1,3) or (3,1)

Question:

Can additional independent prime moduli distinguish the possible
factor branches?

Important:

We DO NOT enumerate Cartesian products of all residue combinations.

Instead we work directly with the actual prime pool and progressively
filter prime pairs using their residue fingerprints.

This keeps memory/time manageable.
================================================================================
"""

from __future__ import annotations

import math
import random
from collections import defaultdict, Counter


# =============================================================================
# CONFIGURATION
# =============================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

ACTUAL_TRIALS = 300
RANDOM_CONTROLS = 300

SEED = 1_511_464_998

# IMPORTANT:
# modulo 4 is the original branch modulus.
#
# The rest are independent prime moduli.
MODULI = [
    4,
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
]


# =============================================================================
# PRIME GENERATION
# =============================================================================

def sieve_primes(lo: int, hi: int) -> list[int]:
    if hi < 2:
        return []

    sieve = bytearray(b"\x01") * (hi + 1)
    sieve[0:2] = b"\x00\x00"

    limit = math.isqrt(hi)

    for p in range(2, limit + 1):
        if sieve[p]:
            start = p * p
            sieve[start:hi + 1:p] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [
        x
        for x in range(max(2, lo), hi + 1)
        if sieve[x]
    ]


# =============================================================================
# DATA
# =============================================================================

class Anchor:
    __slots__ = ("p", "q", "n")

    def __init__(self, p: int, q: int):
        if p > q:
            p, q = q, p

        self.p = p
        self.q = q
        self.n = p * q


# =============================================================================
# RANDOM UNIQUE SEMIPRIMES
# =============================================================================

def generate_anchors(
    count: int,
    primes: list[int],
    rng: random.Random,
) -> list[Anchor]:

    result = []
    seen = set()

    while len(result) < count:

        p = rng.choice(primes)
        q = rng.choice(primes)

        if p > q:
            p, q = q, p

        key = (p, q)

        if key in seen:
            continue

        seen.add(key)
        result.append(Anchor(p, q))

    return result


# =============================================================================
# FINGERPRINTS
# =============================================================================

def residue(p: int, moduli: tuple[int, ...]):
    return tuple(p % r for r in moduli)


def product_fingerprint(
    n: int,
    moduli: tuple[int, ...],
):
    return tuple(n % r for r in moduli)


def factor_fingerprint_unordered(
    p: int,
    q: int,
    moduli: tuple[int, ...],
):
    result = []

    for r in moduli:

        a = p % r
        b = q % r

        if a <= b:
            result.append((a, b))
        else:
            result.append((b, a))

    return tuple(result)


# =============================================================================
# BRANCH REPRESENTATION
# =============================================================================

def branch_mod4(anchor: Anchor):
    """
    Unordered modulo-4 branch.

        (1,1)
        (1,3)

    are the only possibilities for unordered odd factors.
    """

    a = anchor.p % 4
    b = anchor.q % 4

    if a > b:
        a, b = b, a

    return (a, b)


def possible_mod4_branches(n: int):
    r = n % 4

    if r == 1:
        return {(1, 1), (3, 3)}

    if r == 3:
        return {(1, 3)}

    raise RuntimeError("Odd semiprime has invalid mod-4 result.")


# =============================================================================
# BUILD PRIME INDEX
# =============================================================================

print("=" * 100)
print("FAST MULTI-MODULUS FACTOR-BRANCH DISCRIMINATION EXPERIMENT")
print("=" * 100)

print(f"M                    = {M:,}")
print(f"factor range         = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
print(f"actual trials        = {ACTUAL_TRIALS}")
print(f"random controls      = {RANDOM_CONTROLS}")
print(f"moduli               = {MODULI}")
print(f"seed                 = {SEED:,}")

rng = random.Random(SEED)

print()
print("=" * 100)
print("BUILDING PRIME POOL")
print("=" * 100)

primes = sieve_primes(
    FACTOR_MIN,
    FACTOR_MAX,
)

print(f"factor primes        = {len(primes):,}")


# =============================================================================
# PRIME RESIDUE INDEX
# =============================================================================

print()
print("=" * 100)
print("BUILDING RESIDUE INDEX")
print("=" * 100)

# index[r][residue] -> list of primes

index = {}

for r in MODULI:

    groups = defaultdict(list)

    for p in primes:
        groups[p % r].append(p)

    index[r] = groups

    print(
        f"r={r:>2} "
        f"residue classes={len(groups):>3}"
    )


# =============================================================================
# DATASETS
# =============================================================================

actual = generate_anchors(
    ACTUAL_TRIALS,
    primes,
    rng,
)

random_controls = generate_anchors(
    RANDOM_CONTROLS,
    primes,
    rng,
)

print(f"actual anchors       = {len(actual)}")
print(f"random controls      = {len(random_controls)}")


# =============================================================================
# DIRECT FILTERING
# =============================================================================

def count_matching_prime_pairs(
    n: int,
    moduli: tuple[int, ...],
    wanted_branch: tuple[int, int] | None = None,
    limit: int = 10_000,
):
    """
    Count actual prime pairs (a,b) in the factor range satisfying:

        a*b == n (mod r)

    for every r in moduli.

    Optionally restrict modulo-4 branch.

    We avoid enumerating all prime-pair combinations.

    Candidate pairs are generated by choosing the residue of a and deriving
    the required residue of b modulo each modulus.

    For the final exact test we verify the actual integers.

    'limit' prevents pathological explosion.
    """

    # Start from primes grouped by their complete residue vector.

    vectors = defaultdict(list)

    for p in primes:

        vec = tuple(p % r for r in moduli)

        if wanted_branch is not None:

            a4 = p % 4

            if a4 not in wanted_branch:
                continue

        vectors[vec].append(p)

    # The above still creates only one entry per prime, not a giant
    # Cartesian product.

    # We now construct allowed residue vectors by using the actual
    # prime vectors.

    # For this experiment we mostly need whether a branch survives.
    # A direct anchor-preserving lookup is much cheaper.

    actual_count = 0

    # Make a set of prime residue vectors.
    vector_set = set(vectors.keys())

    # Iterate actual primes as candidate first factors.
    #
    # We only need to ask whether a matching second prime exists.

    for a in primes:

        if wanted_branch is not None:

            a4 = a % 4

            # unordered branch
            if a4 not in wanted_branch:
                continue

        required = []

        valid = True

        for r in moduli:

            ar = a % r
            nr = n % r

            # Need:

            # a*b == n (mod r)

            # Therefore if gcd(a,r)=1:

            # b == n * a^-1 (mod r)

            g = math.gcd(ar, r)

            if g != 1:

                valid = False
                break

            b_r = (nr * pow(ar, -1, r)) % r

            required.append(b_r)

        if not valid:
            continue

        # Search only residue-compatible primes.

        candidates = None

        for r, b_r in zip(moduli, required):

            group = index[r].get(b_r)

            if not group:
                candidates = []
                break

            if candidates is None:
                candidates = group
            else:
                # Intersect by actual integer.
                s = set(group)
                candidates = [
                    x for x in candidates
                    if x in s
                ]

            if not candidates:
                break

        if not candidates:
            continue

        for b in candidates:

            if a > b:
                continue

            if a == b:
                continue

            if a * b % r != n % r:
                continue

            if wanted_branch is not None:

                actual_branch = tuple(
                    sorted((a % 4, b % 4))
                )

                if actual_branch != wanted_branch:
                    continue

            actual_count += 1

            if actual_count >= limit:
                return actual_count

    return actual_count


# =============================================================================
# FAST BRANCH SURVIVAL
# =============================================================================

print()
print("=" * 100)
print("BRANCH SURVIVAL")
print("=" * 100)

print(
    f"{'K':>2} "
    f"{'MODULI':<38} "
    f"{'11-BRANCH':>12} "
    f"{'13-BRANCH':>12} "
    f"{'AMBIG':>8} "
    f"{'RESOLVED':>9}"
)

print("-" * 100)


for k in range(1, len(MODULI) + 1):

    current = tuple(MODULI[:k])

    branch_11 = 0
    branch_13 = 0

    ambiguous = 0
    resolved = 0

    for anchor in actual:

        branches = possible_mod4_branches(anchor.n)

        survivors = 0

        if (1, 1) in branches:

            c = count_matching_prime_pairs(
                anchor.n,
                current,
                wanted_branch=(1, 1),
                limit=1,
            )

            if c > 0:
                branch_11 += 1
                survivors += 1

        if (3, 3) in branches:

            c = count_matching_prime_pairs(
                anchor.n,
                current,
                wanted_branch=(3, 3),
                limit=1,
            )

            if c > 0:
                branch_13 += 1
                survivors += 1

        if survivors > 1:
            ambiguous += 1

        elif survivors == 1:
            resolved += 1

    print(
        f"{k:>2} "
        f"{str(current):<38} "
        f"{branch_11:>12} "
        f"{branch_13:>12} "
        f"{ambiguous:>8} "
        f"{resolved:>9}"
    )


# =============================================================================
# ACTUAL FACTOR BRANCH
# =============================================================================

print()
print("=" * 100)
print("ACTUAL BRANCH COUNTS")
print("=" * 100)

branch_counts = Counter()

for anchor in actual:
    branch_counts[branch_mod4(anchor)] += 1

for branch, count in sorted(branch_counts.items()):
    print(
        f"branch={branch} count={count}"
    )


# =============================================================================
# PRODUCT FINGERPRINT COLLISION GROWTH
# =============================================================================

print()
print("=" * 100)
print("PRODUCT FINGERPRINT COLLISION GROWTH")
print("=" * 100)

print(
    f"{'K':>2} "
    f"{'MODULI':<38} "
    f"{'ACT UNIQUE':>12} "
    f"{'MAX COLL':>10} "
    f"{'RND UNIQUE':>12} "
    f"{'RND MAX':>10}"
)

print("-" * 100)

for k in range(1, len(MODULI) + 1):

    current = tuple(MODULI[:k])

    actual_counter = Counter(
        product_fingerprint(a.n, current)
        for a in actual
    )

    random_counter = Counter(
        product_fingerprint(a.n, current)
        for a in random_controls
    )

    print(
        f"{k:>2} "
        f"{str(current):<38} "
        f"{len(actual_counter):>12} "
        f"{max(actual_counter.values()):>10} "
        f"{len(random_counter):>12} "
        f"{max(random_counter.values()):>10}"
    )


# =============================================================================
# FACTOR FINGERPRINT COLLISION GROWTH
# =============================================================================

print()
print("=" * 100)
print("FACTOR FINGERPRINT COLLISION GROWTH")
print("=" * 100)

print(
    f"{'K':>2} "
    f"{'MODULI':<38} "
    f"{'ACT UNIQUE':>12} "
    f"{'ACT MAX':>10} "
    f"{'RND UNIQUE':>12} "
    f"{'RND MAX':>10}"
)

print("-" * 100)

for k in range(1, len(MODULI) + 1):

    current = tuple(MODULI[:k])

    actual_counter = Counter(
        factor_fingerprint_unordered(
            a.p,
            a.q,
            current,
        )
        for a in actual
    )

    random_counter = Counter(
        factor_fingerprint_unordered(
            a.p,
            a.q,
            current,
        )
        for a in random_controls
    )

    print(
        f"{k:>2} "
        f"{str(current):<38} "
        f"{len(actual_counter):>12} "
        f"{max(actual_counter.values()):>10} "
        f"{len(random_counter):>12} "
        f"{max(random_counter.values()):>10}"
    )


# =============================================================================
# EXACT EXAMPLES
# =============================================================================

print()
print("=" * 100)
print("EXAMPLE ANCHORS")
print("=" * 100)

examples = actual[:20]

for anchor in examples:

    print()
    print(
        f"n={anchor.n:,} "
        f"p={anchor.p:,} "
        f"q={anchor.q:,} "
        f"n mod4={anchor.n % 4} "
        f"actual branch={branch_mod4(anchor)}"
    )

    for k in (
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    ):

        current = tuple(MODULI[:k])

        possible = possible_mod4_branches(
            anchor.n
        )

        survivors = []

        for branch in sorted(possible):

            count = count_matching_prime_pairs(
                anchor.n,
                current,
                wanted_branch=branch,
                limit=1,
            )

            if count > 0:
                survivors.append(branch)

        print(
            f"    K={k:>2} "
            f"surviving branches={survivors}"
        )


# =============================================================================
# SANITY
# =============================================================================

print()
print("=" * 100)
print("SANITY CHECK")
print("=" * 100)

failures = 0

for anchor in actual:

    if anchor.p * anchor.q != anchor.n:
        failures += 1

    if anchor.n % 4 != (
        anchor.p * anchor.q
    ) % 4:
        failures += 1

print(
    f"identity failures = {failures}"
)


# =============================================================================
# FINAL
# =============================================================================

print()
print("=" * 100)
print("INTERPRETATION")
print("=" * 100)

print(
"""
The key distinction is:

    PRODUCT fingerprint
        (n mod r1, n mod r2, ...)

versus

    FACTOR fingerprint
        ((p mod r1,q mod r1), ...)

Modulo 4 alone gives the classical branches:

    n == 1 mod 4:
        (1,1) or (3,3)

    n == 3 mod 4:
        (1,3)

The important measurement is not the size of the residue space.

It is:

    how many candidate factor branches remain possible
    after progressively adding independent moduli?

If additional moduli eliminate one of the modulo-4 branches,
then they provide genuine factor-level discrimination.

If both branches remain populated no matter how many moduli are
added, then the extra modular information is only refining the
product/factor fingerprint without resolving that branch.

This experiment deliberately avoids constructing the enormous
Cartesian product of residue possibilities.
"""
)

print()
print("=" * 100)
print("EXPERIMENT COMPLETE")
print("=" * 100)
