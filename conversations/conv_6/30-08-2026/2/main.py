#!/usr/bin/env python3
"""
================================================================================
MULTI-MODULUS FACTOR SIGNATURE / BRANCH DISCRIMINATION EXPERIMENT
================================================================================

Question:

For an odd semiprime n = p*q, modulo 4 gives

    n == 1 (mod 4):
        (p,q) = (1,1) or (3,3)

    n == 3 (mod 4):
        (p,q) = (1,3) or (3,1)

Does adding further prime moduli progressively distinguish the factor classes?

We compare:

    PRODUCT fingerprint:
        (n mod r)

against

    FACTOR fingerprint:
        ((p mod r), (q mod r))

The experiment progressively adds moduli.

No huge NxN product map is constructed.
================================================================================
"""

from __future__ import annotations

import math
import random
from collections import defaultdict, Counter
from dataclasses import dataclass


# =============================================================================
# CONFIGURATION
# =============================================================================

M = 111_546_435
TWO_M = 2 * M

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

ACTUAL_TRIALS = 300
RANDOM_CONTROLS = 300

SEED = 1_511_464_998

# Start with modulo 4 because that is the original hypothesis.
#
# Then progressively add independent prime moduli.
#
# None of these need to divide M. That is intentional.
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

# Number of random factor pairs used for the baseline.
GLOBAL_RANDOM_PAIRS = 50_000


# =============================================================================
# PRIME TESTING
# =============================================================================

def sieve_primes(lo: int, hi: int) -> list[int]:
    """Return all primes in [lo, hi]."""
    if hi < 2:
        return []

    sieve = bytearray(b"\x01") * (hi + 1)
    sieve[0:2] = b"\x00\x00"

    limit = int(math.isqrt(hi))

    for p in range(2, limit + 1):
        if sieve[p]:
            start = p * p
            sieve[start:hi + 1:p] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [n for n in range(max(2, lo), hi + 1) if sieve[n]]


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass(frozen=True)
class Anchor:
    p: int
    q: int

    @property
    def n(self) -> int:
        return self.p * self.q


@dataclass
class Stats:
    product_unique: int = 0
    product_collision: int = 0

    factor_unique: int = 0
    factor_collision: int = 0

    actual_signature_seen: int = 0

    mod4_branch_count: int = 0
    mod4_branch_ambiguous: int = 0
    mod4_branch_resolved: int = 0

    possible_factor_signatures_total: int = 0
    possible_factor_signatures_max: int = 0


# =============================================================================
# PRIME POOL
# =============================================================================

print("=" * 100)
print("MULTI-MODULUS FACTOR SIGNATURE / BRANCH DISCRIMINATION EXPERIMENT")
print("=" * 100)

print(f"M                    = {M:,}")
print(f"2M                   = {TWO_M:,}")
print(f"factor range         = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
print(f"actual trials        = {ACTUAL_TRIALS}")
print(f"random controls      = {RANDOM_CONTROLS}")
print(f"global random pairs  = {GLOBAL_RANDOM_PAIRS:,}")
print(f"moduli               = {MODULI}")
print(f"seed                 = {SEED:,}")

rng = random.Random(SEED)

print()
print("=" * 100)
print("BUILDING PRIME POOL")
print("=" * 100)

primes = sieve_primes(FACTOR_MIN, FACTOR_MAX)

print(f"factor primes        = {len(primes):,}")

if len(primes) < 10:
    raise RuntimeError("Prime pool is unexpectedly small.")


# =============================================================================
# RANDOM ACTUAL ANCHORS
# =============================================================================

def generate_unique_semiprimes(
    count: int,
    prime_pool: list[int],
    rng_: random.Random,
) -> list[Anchor]:

    seen: set[tuple[int, int]] = set()
    result: list[Anchor] = []

    while len(result) < count:
        p = rng_.choice(prime_pool)
        q = rng_.choice(prime_pool)

        if p > q:
            p, q = q, p

        key = (p, q)

        if key in seen:
            continue

        seen.add(key)
        result.append(Anchor(p, q))

    return result


actual = generate_unique_semiprimes(
    ACTUAL_TRIALS,
    primes,
    rng,
)

random_controls = generate_unique_semiprimes(
    RANDOM_CONTROLS,
    primes,
    rng,
)

print(f"actual anchors       = {len(actual)}")
print(f"random controls      = {len(random_controls)}")


# =============================================================================
# PRIME RESIDUE POOLS
# =============================================================================

print()
print("=" * 100)
print("BUILDING PRIME RESIDUE CLASS TABLES")
print("=" * 100)

# For each modulus:
#
# residue -> number of primes
#
# This lets us determine how many residue classes are actually realizable
# by primes in the selected factor interval.

prime_residue_counts: dict[int, Counter[int]] = {}

for r in MODULI:
    c = Counter(p % r for p in primes)
    prime_residue_counts[r] = c

    print(
        f"r={r:>2} "
        f"realizable residues={len(c):>3} "
        f"total primes={sum(c.values()):>5}"
    )


# =============================================================================
# SIGNATURE FUNCTIONS
# =============================================================================

def product_signature(n: int, moduli: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(n % r for r in moduli)


def factor_signature(
    p: int,
    q: int,
    moduli: tuple[int, ...],
) -> tuple[tuple[int, int], ...]:
    return tuple((p % r, q % r) for r in moduli)


def unordered_factor_signature(
    p: int,
    q: int,
    moduli: tuple[int, ...],
) -> tuple[tuple[int, int], ...]:
    """
    Canonicalize (p,q) so that factor ordering does not matter.
    """
    out = []

    for r in moduli:
        a = p % r
        b = q % r

        if a <= b:
            out.append((a, b))
        else:
            out.append((b, a))

    return tuple(out)


# =============================================================================
# MOD-4 BRANCH LOGIC
# =============================================================================

def mod4_branches(n: int) -> tuple[tuple[int, int], ...]:
    """
    Return the possible factor residue signatures modulo 4.

    n == 1:
        (1,1), (3,3)

    n == 3:
        (1,3), (3,1)

    Since we use unordered signatures later, (1,3) and (3,1)
    are equivalent for branch counting.
    """
    r = n % 4

    if r == 1:
        return ((1, 1), (3, 3))

    if r == 3:
        return ((1, 3),)

    raise ValueError("n must be odd.")


# =============================================================================
# COMPATIBLE FACTOR SIGNATURES
# =============================================================================

def compatible_factor_signatures(
    n: int,
    moduli: tuple[int, ...],
) -> set[tuple[tuple[int, int], ...]]:
    """
    Enumerate all residue-pair signatures compatible with

        a*b == n (mod r)

    for every r in moduli,

    where a,b must be residues actually attained by primes
    in the factor range.

    We do NOT enumerate actual prime pairs.

    This operates only in residue space.
    """

    partial: list[list[tuple[int, int]]] = [[]]

    for r in moduli:
        nr = n % r
        residues = sorted(prime_residue_counts[r].keys())

        local: list[tuple[int, int]] = []

        for a in residues:
            for b in residues:
                if (a * b) % r == nr:

                    # Use unordered factor signature.
                    if a <= b:
                        local.append((a, b))
                    else:
                        local.append((b, a))

        local = sorted(set(local))

        next_partial: list[list[tuple[int, int]]] = []

        for prefix in partial:
            for pair in local:
                next_partial.append(prefix + [pair])

        partial = next_partial

        if not partial:
            break

    return {tuple(x) for x in partial}


# =============================================================================
# BRANCH DISCRIMINATION
# =============================================================================

def actual_mod4_branch(anchor: Anchor) -> tuple[int, int]:
    a = anchor.p % 4
    b = anchor.q % 4

    return tuple(sorted((a, b)))  # type: ignore


def branch_is_resolved(
    n: int,
    moduli: tuple[int, ...],
) -> bool:
    """
    For the n-compatible residue signatures, determine whether
    modulo-4 factor branches have been reduced to one possibility.

    We look only at the modulo-4 pair contained in the compatible
    signatures.
    """

    possible = compatible_factor_signatures(n, moduli)

    if not possible:
        return False

    branches = set(sig[0] for sig in possible)

    return len(branches) == 1


# =============================================================================
# CORE ANALYSIS
# =============================================================================

def analyze_dataset(
    anchors: list[Anchor],
    moduli: tuple[int, ...],
) -> tuple[dict, set]:
    """
    Analyze a dataset under one modulus set.
    """

    product_counts: Counter[tuple[int, ...]] = Counter()
    factor_counts: Counter[tuple[tuple[int, int], ...]] = Counter()

    possible_counts: list[int] = []
    branch_ambiguous = 0
    branch_resolved = 0

    for anchor in anchors:

        psig = product_signature(anchor.n, moduli)

        fsig = unordered_factor_signature(
            anchor.p,
            anchor.q,
            moduli,
        )

        product_counts[psig] += 1
        factor_counts[fsig] += 1

        possible = compatible_factor_signatures(
            anchor.n,
            moduli,
        )

        possible_count = len(possible)
        possible_counts.append(possible_count)

        branches = set(sig[0] for sig in possible)

        if len(branches) > 1:
            branch_ambiguous += 1
        elif len(branches) == 1:
            branch_resolved += 1

    stats = {
        "product_unique": len(product_counts),
        "product_max_collision": (
            max(product_counts.values())
            if product_counts
            else 0
        ),
        "factor_unique": len(factor_counts),
        "factor_max_collision": (
            max(factor_counts.values())
            if factor_counts
            else 0
        ),
        "mean_possible_factor_signatures": (
            sum(possible_counts) / len(possible_counts)
            if possible_counts
            else 0.0
        ),
        "max_possible_factor_signatures": (
            max(possible_counts)
            if possible_counts
            else 0
        ),
        "branch_ambiguous": branch_ambiguous,
        "branch_resolved": branch_resolved,
    }

    return stats, factor_counts


# =============================================================================
# GLOBAL RANDOM PAIR BASELINE
# =============================================================================

print()
print("=" * 100)
print("BUILDING GLOBAL RANDOM PRIME-PAIR BASELINE")
print("=" * 100)

random_pairs: list[Anchor] = []

for _ in range(GLOBAL_RANDOM_PAIRS):
    p = rng.choice(primes)
    q = rng.choice(primes)

    if p > q:
        p, q = q, p

    random_pairs.append(Anchor(p, q))

print(f"random pairs         = {len(random_pairs):,}")


# =============================================================================
# MAIN PROGRESSIVE EXPERIMENT
# =============================================================================

print()
print("=" * 100)
print("PROGRESSIVE MODULUS ANALYSIS")
print("=" * 100)

print(
    "The factor-signature column counts distinct unordered residue signatures"
)
print(
    "((p mod r, q mod r)) across the supplied modulus set."
)
print()

print(
    f"{'K':>2} "
    f"{'MODULI':<38} "
    f"{'ACT P-UNIQ':>11} "
    f"{'ACT F-UNIQ':>11} "
    f"{'RND P-UNIQ':>11} "
    f"{'RND F-UNIQ':>11} "
    f"{'MEAN POSS':>10} "
    f"{'AMBIG':>8} "
    f"{'RESOLVED':>9}"
)

print("-" * 125)

for k in range(1, len(MODULI) + 1):

    current = tuple(MODULI[:k])

    actual_stats, _ = analyze_dataset(actual, current)
    random_stats, _ = analyze_dataset(random_pairs[:ACTUAL_TRIALS], current)

    print(
        f"{k:>2} "
        f"{str(current):<38} "
        f"{actual_stats['product_unique']:>11} "
        f"{actual_stats['factor_unique']:>11} "
        f"{random_stats['product_unique']:>11} "
        f"{random_stats['factor_unique']:>11} "
        f"{actual_stats['mean_possible_factor_signatures']:>10.3f} "
        f"{actual_stats['branch_ambiguous']:>8} "
        f"{actual_stats['branch_resolved']:>9}"
    )


# =============================================================================
# MOD-4 SPECIFIC ANALYSIS
# =============================================================================

print()
print("=" * 100)
print("MOD-4 BRANCH ANALYSIS")
print("=" * 100)

mod4 = (4,)

for label, dataset in (
    ("ACTUAL", actual),
    ("RANDOM", random_controls),
):
    print()
    print(label)

    c = Counter()

    for anchor in dataset:
        sig = actual_mod4_branch(anchor)
        c[sig] += 1

    print(f"  n == 1 mod 4       = {sum(1 for a in dataset if a.n % 4 == 1)}")
    print(f"  n == 3 mod 4       = {sum(1 for a in dataset if a.n % 4 == 3)}")

    print("  factor signatures:")

    for sig in sorted(c):
        print(
            f"      {sig} : {c[sig]}"
        )


# =============================================================================
# EXACT MOD-4 BRANCH VS EXTENDED FINGERPRINT
# =============================================================================

print()
print("=" * 100)
print("BRANCH SURVIVAL BY MODULUS SET")
print("=" * 100)

for k in range(1, len(MODULI) + 1):

    current = tuple(MODULI[:k])

    ambiguous = 0
    resolved = 0
    impossible = 0

    for anchor in actual:

        possible = compatible_factor_signatures(
            anchor.n,
            current,
        )

        if not possible:
            impossible += 1
            continue

        branches = set(sig[0] for sig in possible)

        if len(branches) == 1:
            resolved += 1
        else:
            ambiguous += 1

    print(
        f"K={k:>2} "
        f"moduli={current} "
        f"resolved={resolved:>3} "
        f"ambiguous={ambiguous:>3} "
        f"impossible={impossible:>3}"
    )


# =============================================================================
# FIND INTERESTING ANCHORS
# =============================================================================

print()
print("=" * 100)
print("MOST DIFFICULT ANCHORS TO DISCRIMINATE")
print("=" * 100)

final_moduli = tuple(MODULI)

ranked = []

for anchor in actual:

    possible = compatible_factor_signatures(
        anchor.n,
        final_moduli,
    )

    ranked.append(
        (
            len(possible),
            anchor,
            possible,
        )
    )

ranked.sort(reverse=True, key=lambda x: x[0])


for possible_count, anchor, possible in ranked[:20]:

    print()
    print(
        f"n={anchor.n:,} "
        f"p={anchor.p:,} "
        f"q={anchor.q:,} "
        f"n mod 4={anchor.n % 4}"
    )

    print(
        f"actual factor signature = "
        f"{unordered_factor_signature(anchor.p, anchor.q, final_moduli)}"
    )

    print(
        f"compatible signatures   = {possible_count}"
    )

    if possible_count <= 20:
        for sig in sorted(possible):
            print(f"    {sig}")


# =============================================================================
# MOST EASILY DISCRIMINATED ANCHORS
# =============================================================================

print()
print("=" * 100)
print("MOST EASILY DISCRIMINATED ANCHORS")
print("=" * 100)

for possible_count, anchor, possible in sorted(
    ranked,
    key=lambda x: x[0],
)[:20]:

    print(
        f"n={anchor.n:,} "
        f"p={anchor.p:,} "
        f"q={anchor.q:,} "
        f"n mod4={anchor.n % 4} "
        f"compatible signatures={possible_count}"
    )


# =============================================================================
# RANDOM CONTROL COMPARISON
# =============================================================================

print()
print("=" * 100)
print("FINAL MODULUS SET: ACTUAL vs RANDOM")
print("=" * 100)

actual_final, actual_factor_counts = analyze_dataset(
    actual,
    final_moduli,
)

random_final, random_factor_counts = analyze_dataset(
    random_controls,
    final_moduli,
)

for name, stats in (
    ("ACTUAL", actual_final),
    ("RANDOM", random_final),
):

    print()
    print(name)

    print(
        f"  product unique                     = "
        f"{stats['product_unique']}"
    )

    print(
        f"  product max collision               = "
        f"{stats['product_max_collision']}"
    )

    print(
        f"  factor unique                      = "
        f"{stats['factor_unique']}"
    )

    print(
        f"  factor max collision               = "
        f"{stats['factor_max_collision']}"
    )

    print(
        f"  mean compatible factor signatures  = "
        f"{stats['mean_possible_factor_signatures']:.4f}"
    )

    print(
        f"  max compatible factor signatures   = "
        f"{stats['max_possible_factor_signatures']}"
    )

    print(
        f"  branch ambiguous                    = "
        f"{stats['branch_ambiguous']}"
    )

    print(
        f"  branch resolved                     = "
        f"{stats['branch_resolved']}"
    )


# =============================================================================
# FINGERPRINT GROWTH
# =============================================================================

print()
print("=" * 100)
print("FINGERPRINT GROWTH")
print("=" * 100)

print(
    f"{'K':>2} "
    f"{'MODULI':<38} "
    f"{'ACT FACTOR UNIQUE':>19} "
    f"{'RND FACTOR UNIQUE':>19} "
    f"{'ACT/RND':>10}"
)

print("-" * 100)

for k in range(1, len(MODULI) + 1):

    current = tuple(MODULI[:k])

    a_stats, _ = analyze_dataset(actual, current)
    r_stats, _ = analyze_dataset(random_controls, current)

    ratio = (
        a_stats["factor_unique"] / r_stats["factor_unique"]
        if r_stats["factor_unique"]
        else 0.0
    )

    print(
        f"{k:>2} "
        f"{str(current):<38} "
        f"{a_stats['factor_unique']:>19} "
        f"{r_stats['factor_unique']:>19} "
        f"{ratio:>10.4f}"
    )


# =============================================================================
# SANITY CHECKS
# =============================================================================

print()
print("=" * 100)
print("SANITY CHECKS")
print("=" * 100)

failures = 0

for anchor in actual:

    for r in final_moduli:

        if (anchor.p * anchor.q) % r != anchor.n % r:
            print(
                f"FAIL product identity: "
                f"n={anchor.n} r={r}"
            )
            failures += 1

print(
    f"product residue identity failures = {failures}"
)


# Verify the mod-4 rule itself.

rule_failures = 0

for anchor in actual:

    n4 = anchor.n % 4
    p4 = anchor.p % 4
    q4 = anchor.q % 4

    valid = False

    if n4 == 1 and {p4, q4} in ({1}, {3}):
        valid = True

    if n4 == 3 and {p4, q4} == {1, 3}:
        valid = True

    if not valid:
        rule_failures += 1
        print(
            f"FAIL mod4 factor rule: "
            f"p={anchor.p} q={anchor.q} "
            f"p4={p4} q4={q4} n4={n4}"
        )

print(
    f"mod-4 branch rule failures = {rule_failures}"
)


# =============================================================================
# FINAL INTERPRETATION
# =============================================================================

print()
print("=" * 100)
print("INTERPRETATION")
print("=" * 100)

print(
"""
This experiment separates three levels of information.

LEVEL 1
-------

Product fingerprint:

    n -> (n mod r1, n mod r2, ...)

This tells us the residue of the PRODUCT.

LEVEL 2
-------

Factor fingerprint:

    (p,q) -> ((p mod r1,q mod r1),
              (p mod r2,q mod r2), ...)

This retains information about the individual factors.

LEVEL 3
-------

Factor-branch discrimination.

For modulo 4:

    n == 1 mod 4
        -> (1,1) or (3,3)

    n == 3 mod 4
        -> (1,3)

after ignoring factor ordering.

The question is whether additional moduli reduce the possible
factor signatures enough to determine the actual factor branch.

A particularly interesting outcome would be:

    modulo 4:
        2 possible branches

    + one additional modulus:
        often still 2 branches

    + several additional moduli:
        branch collapses to 1

That would demonstrate progressive factor-level discrimination.

However, there is an important distinction:

    More moduli always make the PRODUCT fingerprint finer.

That is expected.

The nontrivial question is whether the same additional moduli
also make the INDIVIDUAL FACTOR fingerprint sufficiently restrictive
to distinguish factor branches.

The experiment therefore reports both quantities separately.
"""
)

print()
print("=" * 100)
print("EXPERIMENT COMPLETE")
print("=" * 100)
