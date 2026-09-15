#!/usr/bin/env python3

"""
================================================================================
MULTI-MODULUS FACTOR-BRANCH TEST — HASHED FINGERPRINT VERSION
================================================================================

Question:

For n = p*q and

    n == 1 (mod 4)

the two modulo-4 factor branches are

    B11 = (p mod 4, q mod 4) = (1,1)
    B33 = (p mod 4, q mod 4) = (3,3)

Can additional independent moduli distinguish those branches?

We progressively add:

    4
    3
    5
    7
    11
    13
    17
    19
    23
    29
    31
    37
    41

IMPORTANT:

We do NOT enumerate Cartesian products.

For each prime p we precompute its complete residue fingerprint.
For an anchor n and a possible p, the required residue of q modulo r is

    q == n * inverse(p) (mod r)

when p is invertible modulo r.

The existence question is then reduced to dictionary lookups.

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
    sieve = bytearray(b"\x01") * (hi + 1)

    sieve[0:2] = b"\x00\x00"

    for p in range(2, math.isqrt(hi) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:hi + 1:p] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [
        x
        for x in range(lo, hi + 1)
        if sieve[x]
    ]


# =============================================================================
# ANCHORS
# =============================================================================

class Anchor:
    __slots__ = ("p", "q", "n")

    def __init__(self, p: int, q: int):
        if p > q:
            p, q = q, p

        self.p = p
        self.q = q
        self.n = p * q


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
# MOD-4 BRANCH
# =============================================================================

def actual_branch(anchor: Anchor):
    a = anchor.p % 4
    b = anchor.q % 4

    if a > b:
        a, b = b, a

    return (a, b)


def possible_branches(n: int):
    r = n % 4

    if r == 1:
        return [(1, 1), (3, 3)]

    if r == 3:
        return [(1, 3)]

    raise ValueError("Invalid odd-product residue modulo 4.")


# =============================================================================
# PRIME FINGERPRINTS
# =============================================================================

def make_fingerprint_map(
    primes: list[int],
    moduli: tuple[int, ...],
):
    """
    fingerprint -> set/list of primes

    A fingerprint contains all prime residues.
    """

    result = defaultdict(list)

    for p in primes:

        fp = tuple(
            p % r
            for r in moduli
        )

        result[fp].append(p)

    return result


# =============================================================================
# BRANCH-SPECIFIC PRIME MAPS
# =============================================================================

def make_branch_maps(
    primes: list[int],
    moduli: tuple[int, ...],
):
    """
    Separate primes according to p mod 4.

        branch_maps[1]
        branch_maps[3]

    Each maps:

        fingerprint_without_mod4 -> primes
    """

    # Mod-4 is handled separately.

    non4 = tuple(
        r
        for r in moduli
        if r != 4
    )

    maps = {
        1: defaultdict(list),
        3: defaultdict(list),
    }

    for p in primes:

        b = p % 4

        if b not in maps:
            continue

        fp = tuple(
            p % r
            for r in non4
        )

        maps[b][fp].append(p)

    return maps


# =============================================================================
# COMPLETE REQUIRED Q FINGERPRINT
# =============================================================================

def required_q_fingerprint(
    n: int,
    p: int,
    moduli: tuple[int, ...],
):
    """
    For every modulus r:

        p*q == n (mod r)

    requires

        q == n * p^(-1) (mod r)

    provided gcd(p,r)=1.
    """

    out = []

    for r in moduli:

        pr = p % r

        if math.gcd(pr, r) != 1:
            return None

        q_r = (
            (n % r)
            * pow(pr, -1, r)
        ) % r

        out.append(q_r)

    return tuple(out)


# =============================================================================
# TEST BRANCH
# =============================================================================

def branch_has_candidate(
    n: int,
    branch: tuple[int, int],
    moduli: tuple[int, ...],
    branch_maps,
):
    """
    Determine whether there exists ANY prime pair (a,b) from the factor
    range satisfying:

        a mod 4 == branch[0]
        b mod 4 == branch[1]

    and

        a*b == n mod r

    for every modulus r != 4.

    We only need existence, so stop immediately when a candidate exists.
    """

    a_class, b_class = branch

    non4 = tuple(
        r
        for r in moduli
        if r != 4
    )

    # Only iterate over primes belonging to the first branch.
    #
    # This is <= 8363 candidates.

    prime_map = branch_maps[a_class]

    for fp_a, candidates_a in prime_map.items():

        for a in candidates_a:

            # Required residues for b.

            req = required_q_fingerprint(
                n,
                a,
                non4,
            )

            if req is None:
                continue

            if req not in branch_maps[b_class]:
                continue

            candidates_b = branch_maps[b_class][req]

            # Need actual pair ordering only for uniqueness;
            # existence itself is enough here.

            if not candidates_b:
                continue

            return True

    return False


# =============================================================================
# EVEN FASTER VERSION:
# BUILD COMPLETE PRODUCT FINGERPRINT MAP
# =============================================================================

def build_product_map(
    primes: list[int],
    moduli: tuple[int, ...],
):
    """
    Map product fingerprints to whether each modulo-4 branch is realizable.

    Key:

        ((n mod r1), ...)

    Value:

        bit 1 -> (1,1) branch exists
        bit 2 -> (3,3) branch exists
        bit 3 -> (1,3) branch exists

    This lets us answer branch-survival for an n using ONE lookup.
    """

    non4 = tuple(
        r
        for r in moduli
        if r != 4
    )

    result = {}

    # Build residue lookup for each prime branch.

    by_branch = {
        1: defaultdict(list),
        3: defaultdict(list),
    }

    for p in primes:

        b = p % 4

        fp = tuple(
            p % r
            for r in non4
        )

        by_branch[b][fp].append(p)

    # Candidate factor pairs are generated by residue fingerprints,
    # not raw Cartesian products.

    #
    # We do not need to retain all actual pairs.
    #
    # We instead generate the product fingerprint for each residue
    # fingerprint combination encountered.
    #

    for branch_a, branch_b in (
        (1, 1),
        (3, 3),
        (1, 3),
    ):

        seen_products = set()

        for fp_a in by_branch[branch_a].keys():

            # For every possible b fingerprint.
            #
            # The number of fingerprints is only the number of primes,
            # not the huge residue Cartesian space.

            for fp_b in by_branch[branch_b].keys():

                prod_fp = tuple(
                    (a * b) % r
                    for a, b, r
                    in zip(fp_a, fp_b, non4)
                )

                seen_products.add(prod_fp)

        bit = {
            (1, 1): 1,
            (3, 3): 2,
            (1, 3): 4,
        }[(branch_a, branch_b)]

        for fp in seen_products:
            result[fp] = result.get(fp, 0) | bit

    return result


# =============================================================================
# DIRECT PRIME-INDEX METHOD
# =============================================================================

def build_branch_product_sets(
    primes: list[int],
    moduli: tuple[int, ...],
):
    """
    Create sets of product fingerprints for each modulo-4 branch.

    This is still finite and manageable for this factor range, but we
    progressively build it one modulus at a time to avoid a giant
    intermediate Cartesian space.
    """

    non4 = tuple(
        r
        for r in moduli
        if r != 4
    )

    # Prime residue fingerprints grouped by branch.

    branch_fps = {
        1: set(),
        3: set(),
    }

    for p in primes:

        b = p % 4

        fp = tuple(
            p % r
            for r in non4
        )

        branch_fps[b].add(fp)

    return branch_fps


# =============================================================================
# PROGRESSIVE BRANCH SURVIVAL VIA ACTUAL PRIME PAIRS
# =============================================================================

def branch_survival_for_anchor(
    anchor: Anchor,
    primes: list[int],
    selected_moduli: tuple[int, ...],
):
    """
    Exact test.

    Searches all possible first factors, but performs only O(number of
    primes) dictionary membership checks.

    No pair Cartesian product.
    """

    n = anchor.n

    # Relevant branch set.
    branches = possible_branches(n)

    non4 = tuple(
        r
        for r in selected_moduli
        if r != 4
    )

    result = []

    # Index primes by complete residue vector for the selected moduli.

    maps = {
        1: defaultdict(bool),
        3: defaultdict(bool),
    }

    for p in primes:

        b = p % 4

        maps[b][
            tuple(p % r for r in non4)
        ] = True

    for branch in branches:

        a_class, b_class = branch

        survives = False

        for a in primes:

            if a % 4 != a_class:
                continue

            req = []

            valid = True

            for r in non4:

                ar = a % r

                if math.gcd(ar, r) != 1:
                    valid = False
                    break

                qr = (
                    (n % r)
                    * pow(ar, -1, r)
                ) % r

                req.append(qr)

            if not valid:
                continue

            if tuple(req) in maps[b_class]:
                survives = True
                break

        if survives:
            result.append(branch)

    return result


# =============================================================================
# MAIN
# =============================================================================

print("=" * 100)
print("FAST MULTI-MODULUS FACTOR-BRANCH SURVIVAL EXPERIMENT")
print("=" * 100)

print(f"M                    = {M:,}")
print(
    f"factor range         = "
    f"{FACTOR_MIN:,} - {FACTOR_MAX:,}"
)
print(f"actual trials        = {ACTUAL_TRIALS}")
print(f"random controls      = {RANDOM_CONTROLS}")
print(f"moduli               = {MODULI}")
print(f"seed                 = {SEED:,}")


# =============================================================================
# PRIME POOL
# =============================================================================

print()
print("=" * 100)
print("BUILDING PRIME POOL")
print("=" * 100)

primes = sieve_primes(
    FACTOR_MIN,
    FACTOR_MAX,
)

print(
    f"factor primes        = {len(primes):,}"
)


# =============================================================================
# DATA
# =============================================================================

rng = random.Random(SEED)

actual = generate_anchors(
    ACTUAL_TRIALS,
    primes,
    rng,
)

controls = generate_anchors(
    RANDOM_CONTROLS,
    primes,
    rng,
)

print(
    f"actual anchors       = {len(actual)}"
)

print(
    f"random controls      = {len(controls)}"
)


# =============================================================================
# PROGRESSIVE TEST
# =============================================================================

print()
print("=" * 100)
print("PROGRESSIVE MODULUS BRANCH SURVIVAL")
print("=" * 100)

print(
    f"{'K':>2} "
    f"{'MODULI':<55} "
    f"{'BRANCH11':>10} "
    f"{'BRANCH33':>10} "
    f"{'AMBIG':>8} "
    f"{'RESOLVED':>10}"
)

print("-" * 100)


# We intentionally use a small number of progressive stages first.
#
# The last stages use all moduli but still do NOT construct a Cartesian
# product.

for k in range(1, len(MODULI) + 1):

    selected = tuple(
        MODULI[:k]
    )

    branch11 = 0
    branch33 = 0

    ambiguous = 0
    resolved = 0

    for anchor in actual:

        surviving = branch_survival_for_anchor(
            anchor,
            primes,
            selected,
        )

        count = len(surviving)

        if (1, 1) in surviving:
            branch11 += 1

        if (3, 3) in surviving:
            branch33 += 1

        if count >= 2:
            ambiguous += 1

        elif count == 1:
            resolved += 1

    print(
        f"{k:>2} "
        f"{str(selected):<55} "
        f"{branch11:>10} "
        f"{branch33:>10} "
        f"{ambiguous:>8} "
        f"{resolved:>10}"
    )


# =============================================================================
# ACTUAL BRANCH
# =============================================================================

print()
print("=" * 100)
print("ACTUAL MOD-4 BRANCH DISTRIBUTION")
print("=" * 100)

counter = Counter(
    actual_branch(a)
    for a in actual
)

for branch, count in sorted(counter.items()):
    print(
        f"{branch}: {count}"
    )


# =============================================================================
# FULL FINGERPRINT UNIQUENESS
# =============================================================================

print()
print("=" * 100)
print("FACTOR FINGERPRINT UNIQUENESS")
print("=" * 100)

print(
    f"{'K':>2} "
    f"{'MODULI':<55} "
    f"{'ACT UNIQUE':>12} "
    f"{'ACT MAX':>10} "
    f"{'RND UNIQUE':>12} "
    f"{'RND MAX':>10}"
)

print("-" * 100)

for k in range(1, len(MODULI) + 1):

    selected = tuple(MODULI[:k])

    actual_fp = Counter(
        (
            tuple(a.p % r for r in selected),
            tuple(a.q % r for r in selected),
        )
        for a in actual
    )

    random_fp = Counter(
        (
            tuple(a.p % r for r in selected),
            tuple(a.q % r for r in selected),
        )
        for a in controls
    )

    print(
        f"{k:>2} "
        f"{str(selected):<55} "
        f"{len(actual_fp):>12} "
        f"{max(actual_fp.values()):>10} "
        f"{len(random_fp):>12} "
        f"{max(random_fp.values()):>10}"
    )


# =============================================================================
# EXAMPLES
# =============================================================================

print()
print("=" * 100)
print("EXAMPLE ANCHORS")
print("=" * 100)

shown = 0

for anchor in actual:

    if anchor.n % 4 != 1:
        continue

    print()
    print(
        f"n={anchor.n:,} "
        f"p={anchor.p:,} "
        f"q={anchor.q:,}"
    )

    print(
        f"actual branch = {actual_branch(anchor)}"
    )

    for k in range(1, len(MODULI) + 1):

        selected = tuple(
            MODULI[:k]
        )

        surviving = branch_survival_for_anchor(
            anchor,
            primes,
            selected,
        )

        print(
            f"  K={k:>2} "
            f"surviving={surviving}"
        )

    shown += 1

    if shown >= 10:
        break


# =============================================================================
# SANITY
# =============================================================================

print()
print("=" * 100)
print("SANITY")
print("=" * 100)

failures = 0

for a in actual:

    if a.p * a.q != a.n:
        failures += 1

    if (
        a.n % 4
        !=
        (a.p * a.q) % 4
    ):
        failures += 1

print(
    f"identity failures = {failures}"
)


# =============================================================================
# INTERPRETATION
# =============================================================================

print()
print("=" * 100)
print("INTERPRETATION")
print("=" * 100)

print(
"""
This experiment tests your original idea directly.

For n == 1 (mod 4):

    branch 11:
        p == 1 (mod 4)
        q == 1 (mod 4)

    branch 33:
        p == 3 (mod 4)
        q == 3 (mod 4)

For each additional modulus r we require:

    p*q == n (mod r)

so once p is selected,

    q == n * p^(-1) (mod r).

The question is whether one branch stops having ANY prime realization.

IMPORTANT:

A branch surviving does NOT mean the actual factorization is ambiguous
in a cryptographic sense.

It means only that another prime pair in the tested factor interval
can reproduce the same modular observations.

That distinction is central.

If both branches survive all supplied independent moduli, then adding
more of these moduli is not resolving the original modulo-4 branch.

If one branch disappears at some K, then that K provides genuine
branch discrimination over the chosen factor range.

The final full-modulus case should also be compared against the product
modulus itself. Since the moduli are coprime,

    3*5*7*11*13*17*19*23*29*31*37*41

contains vastly more residue information than M alone.
"""
)

print()
print("=" * 100)
print("EXPERIMENT COMPLETE")
print("=" * 100)
