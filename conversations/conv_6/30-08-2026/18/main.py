#!/usr/bin/env python3
"""
================================================================================
THREE-MODULUS DISCRIMINANT / SUM-OF-FACTORS RECOVERY EXPERIMENT
================================================================================

Goal
----
Instead of searching directly for p, search for

    s = p + q

using

    D = s^2 - 4*n = (p-q)^2.

For the true factorization:

    D must be a perfect square.

For every modulus r:

    D mod r must be a quadratic residue modulo r.

Therefore the three selected moduli provide three independent
quadratic-residue filters.

We additionally exploit the mod-4 branch immediately:

    n == 1 (mod 4)  -> p,q have equal mod-4 classes
                       -> s == 2 (mod 4)

    n == 3 (mod 4)  -> p,q have opposite mod-4 classes
                       -> s == 0 (mod 4)

The experiment compares:

    BASELINE
        every admissible s in the factor-range sum interval

    MOD-4 FILTER
        only s with the correct parity class modulo 4

    r1 FILTER
        quadratic-residue condition modulo r1

    r1+r2 FILTER

    r1+r2+r3 FILTER

    EXACT
        D is an actual perfect square and produces the factors

Important:
------------
The recovery stage does NOT use the precomputed factor-prime population.
The only use of primes in the recovery stage is for choosing r1,r2,r3.

This is therefore a genuine test of whether the modular structure itself
reduces the factor search.
"""

from __future__ import annotations

import math
import random
from bisect import bisect_right
from dataclasses import dataclass


# =============================================================================
# CONFIGURATION
# =============================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MODULUS_MIN = 300
MODULUS_MAX = 3_000

ACTUAL_ANCHORS = 300

CLOSE_RATIO = 0.20

SEED = 1_511_464_998

PROGRESS_EVERY = 25


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


@dataclass(frozen=True)
class Triple:
    r1: int
    r2: int
    r3: int

    @property
    def product(self) -> int:
        return self.r1 * self.r2 * self.r3


@dataclass
class RecoveryResult:
    n: int
    p: int
    q: int

    triple: Triple
    R: int

    s_low: int
    s_high: int

    baseline_count: int
    mod4_count: int
    qr1_count: int
    qr2_count: int
    qr3_count: int

    exact_square_count: int

    recovered: bool
    recovered_p: int | None
    recovered_q: int | None

    checked_after_r3: int
    exact_checks: int


# =============================================================================
# PRIME GENERATION
# =============================================================================

def sieve(limit: int) -> list[int]:
    """Return all primes <= limit."""
    if limit < 2:
        return []

    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"

    root = int(math.isqrt(limit))

    for p in range(2, root + 1):
        if is_prime[p]:
            start = p * p
            is_prime[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i in range(2, limit + 1) if is_prime[i]]


# =============================================================================
# ANCHOR GENERATION
# =============================================================================

def generate_anchors(
    factor_primes: list[int],
    count: int,
    rng: random.Random,
) -> list[Anchor]:
    """
    Deterministically generate distinct unordered prime pairs.
    """
    anchors: list[Anchor] = []
    used: set[tuple[int, int]] = set()

    while len(anchors) < count:
        p, q = rng.sample(factor_primes, 2)

        if p > q:
            p, q = q, p

        key = (p, q)

        if key in used:
            continue

        used.add(key)
        anchors.append(Anchor(p, q))

    return anchors


# =============================================================================
# MAXIMUM-PRODUCT CLOSE PRIME TRIPLE
# =============================================================================

def choose_max_product_triple(
    n: int,
    modulus_primes: list[int],
    close_ratio: float,
) -> Triple:
    """
    Find r1 < r2 < r3 such that:

        r1*r2*r3 < n

    and

        r3/r1 <= 1 + close_ratio

    while maximizing r1*r2*r3.

    The search is optimized using sorted primes and binary search.
    """
    best: Triple | None = None
    best_product = -1

    primes = modulus_primes
    L = len(primes)

    for i in range(L - 2):
        r1 = primes[i]

        # Anything with r1 too large could still be possible,
        # but if even the largest possible product cannot beat best,
        # we can stop.
        if i + 2 >= L:
            break

        for j in range(i + 1, L - 1):
            r2 = primes[j]

            # Close-ratio restriction:
            #
            # r3 <= r1 * (1 + close_ratio)
            upper_ratio = int(math.floor(r1 * (1.0 + close_ratio)))

            k_limit_ratio = bisect_right(primes, upper_ratio) - 1

            if k_limit_ratio <= j:
                continue

            k_max = k_limit_ratio

            # Product must remain < n.
            max_r3_by_n = (n - 1) // (r1 * r2)

            k_limit_product = bisect_right(primes, max_r3_by_n) - 1

            if k_limit_product <= j:
                continue

            k = min(k_max, k_limit_product)

            if k <= j:
                continue

            r3 = primes[k]
            product = r1 * r2 * r3

            if product >= n:
                continue

            if product > best_product:
                best_product = product
                best = Triple(r1, r2, r3)

    if best is None:
        raise RuntimeError(
            f"No valid modulus triple found for n={n}"
        )

    return best


# =============================================================================
# MOD-4 SUM RESTRICTION
# =============================================================================

def sum_mod4_target(n: int) -> int:
    """
    For odd p,q:

        n == 1 mod 4:
            p,q are same mod-4 class
            p+q == 2 mod 4

        n == 3 mod 4:
            p,q are opposite mod-4 classes
            p+q == 0 mod 4
    """
    nr = n & 3

    if nr == 1:
        return 2

    if nr == 3:
        return 0

    raise ValueError(
        f"Expected odd semiprime n, got n mod 4 = {nr}"
    )


# =============================================================================
# QUADRATIC RESIDUE TABLES
# =============================================================================

def build_qr_table(r: int) -> bytearray:
    """
    qr[x] == 1 iff x is a quadratic residue modulo prime r.
    Includes 0.
    """
    qr = bytearray(r)

    for x in range(r):
        qr[(x * x) % r] = 1

    return qr


# =============================================================================
# SUM RANGE
# =============================================================================

def admissible_sum_range(p: int, q: int) -> tuple[int, int]:
    """
    Since both factors are in [FACTOR_MIN, FACTOR_MAX],

        s = p+q

    is bounded above by 2*FACTOR_MAX.

    The lower bound can be strengthened using p*q=n:

        s >= 2*sqrt(n)

    by AM-GM.

    We use ceil(2*sqrt(n)).
    """
    n = p * q

    low = math.isqrt(4 * n)

    # ceil(sqrt(4n))
    if low * low < 4 * n:
        low += 1

    high = 2 * FACTOR_MAX

    return low, high


# =============================================================================
# EXACT FACTOR RECOVERY FROM SUM
# =============================================================================

def recover_from_sum(
    n: int,
    s: int,
) -> tuple[int, int] | None:
    """
    Given s=p+q:

        D = s^2 - 4n

    If D=d^2, then

        p=(s-d)/2
        q=(s+d)/2.
    """
    D = s * s - 4 * n

    if D < 0:
        return None

    d = math.isqrt(D)

    if d * d != D:
        return None

    if (s - d) & 1:
        return None

    p = (s - d) // 2
    q = (s + d) // 2

    if p < FACTOR_MIN or q > FACTOR_MAX:
        return None

    if p * q != n:
        return None

    return p, q


# =============================================================================
# SINGLE ANCHOR RECOVERY
# =============================================================================

def analyze_anchor(
    anchor: Anchor,
    triple: Triple,
    qr_tables: dict[int, bytearray],
) -> RecoveryResult:

    n = anchor.n

    r1 = triple.r1
    r2 = triple.r2
    r3 = triple.r3

    qr1 = qr_tables[r1]
    qr2 = qr_tables[r2]
    qr3 = qr_tables[r3]

    s_low, s_high = admissible_sum_range(anchor.p, anchor.q)

    baseline_count = 0
    mod4_count = 0
    qr1_count = 0
    qr2_count = 0
    qr3_count = 0
    exact_square_count = 0

    recovered_p = None
    recovered_q = None

    exact_checks = 0

    target_mod4 = sum_mod4_target(n)

    # Since the factors are odd, their sum is even.
    start = s_low

    if start & 1:
        start += 1

    for s in range(start, s_high + 1, 2):

        baseline_count += 1

        # ------------------------------------------------------------
        # FIRST FILTER: MOD 4 BRANCH
        # ------------------------------------------------------------
        if (s & 3) != target_mod4:
            continue

        mod4_count += 1

        D = s * s - 4 * n

        if D < 0:
            continue

        # ------------------------------------------------------------
        # SECOND FILTER: r1
        # ------------------------------------------------------------
        if not qr1[D % r1]:
            continue

        qr1_count += 1

        # ------------------------------------------------------------
        # THIRD FILTER: r2
        # ------------------------------------------------------------
        if not qr2[D % r2]:
            continue

        qr2_count += 1

        # ------------------------------------------------------------
        # FOURTH FILTER: r3
        # ------------------------------------------------------------
        if not qr3[D % r3]:
            continue

        qr3_count += 1
        exact_checks += 1

        # ------------------------------------------------------------
        # EXACT TEST
        # ------------------------------------------------------------
        d = math.isqrt(D)

        if d * d != D:
            continue

        exact_square_count += 1

        result = recover_from_sum(n, s)

        if result is None:
            continue

        rp, rq = result

        # Avoid accidentally accepting the factors in reversed order.
        if rp > rq:
            rp, rq = rq, rp

        recovered_p = rp
        recovered_q = rq

    recovered = (
        recovered_p == anchor.p
        and recovered_q == anchor.q
    )

    return RecoveryResult(
        n=n,
        p=anchor.p,
        q=anchor.q,
        triple=triple,
        R=triple.product,
        s_low=s_low,
        s_high=s_high,
        baseline_count=baseline_count,
        mod4_count=mod4_count,
        qr1_count=qr1_count,
        qr2_count=qr2_count,
        qr3_count=qr3_count,
        exact_square_count=exact_square_count,
        recovered=recovered,
        recovered_p=recovered_p,
        recovered_q=recovered_q,
        checked_after_r3=qr3_count,
        exact_checks=exact_checks,
    )


# =============================================================================
# REPORT HELPERS
# =============================================================================

def percentage(a: int, b: int) -> float:
    if b == 0:
        return 0.0
    return 100.0 * a / b


def reduction(before: int, after: int) -> float:
    if before == 0:
        return 0.0
    return 100.0 * (1.0 - after / before)


# =============================================================================
# MAIN
# =============================================================================

def run() -> None:

    rng = random.Random(SEED)

    print("=" * 100)
    print("THREE-MODULUS DISCRIMINANT / SUM-OF-FACTORS RECOVERY EXPERIMENT")
    print("=" * 100)
    print(f"M                         = {M:,}")
    print(f"factor range             = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"actual anchors            = {ACTUAL_ANCHORS}")
    print(
        f"modulus prime range      = "
        f"{MODULUS_MIN:,} - {MODULUS_MAX:,}"
    )
    print(f"close ratio              = {CLOSE_RATIO:.0%}")
    print(f"seed                      = {SEED:,}")

    # -------------------------------------------------------------------------
    # PRIME POOLS
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    factor_primes_all = sieve(FACTOR_MAX)
    factor_primes = [
        p for p in factor_primes_all
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    modulus_primes_all = sieve(MODULUS_MAX)
    modulus_primes = [
        p for p in modulus_primes_all
        if MODULUS_MIN <= p <= MODULUS_MAX
    ]

    print(f"factor primes             = {len(factor_primes):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")

    # -------------------------------------------------------------------------
    # ACTUAL ANCHORS
    # -------------------------------------------------------------------------

    anchors = generate_anchors(
        factor_primes,
        ACTUAL_ANCHORS,
        rng,
    )

    print(f"actual anchors            = {len(anchors):,}")

    # -------------------------------------------------------------------------
    # QUADRATIC RESIDUE TABLES
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("BUILDING QUADRATIC-RESIDUE TABLES")
    print("=" * 100)

    qr_tables: dict[int, bytearray] = {}

    # Tables are tiny because r <= 3000.
    for r in modulus_primes:
        qr_tables[r] = build_qr_table(r)

    print("done")

    # -------------------------------------------------------------------------
    # SELECT THREE CLOSE MAXIMUM-PRODUCT MODULI
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("SELECTING MAXIMUM-PRODUCT CLOSE PRIME TRIPLES")
    print("=" * 100)

    triples: list[Triple] = []

    for idx, anchor in enumerate(anchors, 1):
        triple = choose_max_product_triple(
            anchor.n,
            modulus_primes,
            CLOSE_RATIO,
        )

        triples.append(triple)

        if idx % PROGRESS_EVERY == 0 or idx == len(anchors):
            print(f"anchor {idx:3d}/{len(anchors)}")

    # -------------------------------------------------------------------------
    # RECOVERY
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("DIRECT RECOVERY")
    print("=" * 100)

    results: list[RecoveryResult] = []

    for idx, (anchor, triple) in enumerate(
        zip(anchors, triples),
        1,
    ):
        result = analyze_anchor(
            anchor,
            triple,
            qr_tables,
        )

        results.append(result)

        if idx % PROGRESS_EVERY == 0 or idx == len(anchors):
            print(f"anchor {idx:3d}/{len(anchors)}")

    # -------------------------------------------------------------------------
    # SANITY
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("SANITY CHECK")
    print("=" * 100)

    failures = 0

    for result in results:
        if not result.recovered:
            failures += 1

    print(f"recovery failures         = {failures}")
    print(
        f"recovered correctly       = "
        f"{len(results) - failures}/{len(results)}"
    )

    # -------------------------------------------------------------------------
    # AGGREGATE COUNTS
    # -------------------------------------------------------------------------

    sum_baseline = sum(r.baseline_count for r in results)
    sum_mod4 = sum(r.mod4_count for r in results)
    sum_qr1 = sum(r.qr1_count for r in results)
    sum_qr2 = sum(r.qr2_count for r in results)
    sum_qr3 = sum(r.qr3_count for r in results)
    sum_exact = sum(r.exact_square_count for r in results)

    print()
    print("=" * 100)
    print("GLOBAL SEARCH REDUCTION")
    print("=" * 100)

    print(
        f"baseline s candidates     = "
        f"{sum_baseline:,}"
    )

    print(
        f"after mod-4               = "
        f"{sum_mod4:,}"
    )

    print(
        f"after r1 QR filter        = "
        f"{sum_qr1:,}"
    )

    print(
        f"after r1+r2 QR filters    = "
        f"{sum_qr2:,}"
    )

    print(
        f"after r1+r2+r3 filters    = "
        f"{sum_qr3:,}"
    )

    print(
        f"exact square solutions     = "
        f"{sum_exact:,}"
    )

    print()
    print(
        f"mod-4 reduction           = "
        f"{reduction(sum_baseline, sum_mod4):.4f}%"
    )

    print(
        f"r1 reduction              = "
        f"{reduction(sum_mod4, sum_qr1):.4f}%"
    )

    print(
        f"r2 additional reduction   = "
        f"{reduction(sum_qr1, sum_qr2):.4f}%"
    )

    print(
        f"r3 additional reduction   = "
        f"{reduction(sum_qr2, sum_qr3):.4f}%"
    )

    print(
        f"TOTAL modular reduction   = "
        f"{reduction(sum_baseline, sum_qr3):.4f}%"
    )

    # -------------------------------------------------------------------------
    # PER-ANCHOR AVERAGES
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("MEAN CANDIDATE COUNTS PER ANCHOR")
    print("=" * 100)

    count = len(results)

    print(
        f"baseline                   = "
        f"{sum_baseline / count:.3f}"
    )

    print(
        f"after mod-4                = "
        f"{sum_mod4 / count:.3f}"
    )

    print(
        f"after r1                   = "
        f"{sum_qr1 / count:.3f}"
    )

    print(
        f"after r1+r2               = "
        f"{sum_qr2 / count:.3f}"
    )

    print(
        f"after r1+r2+r3           = "
        f"{sum_qr3 / count:.3f}"
    )

    # -------------------------------------------------------------------------
    # MODULUS INFORMATION
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("SELECTED MODULUS TRIPLES")
    print("=" * 100)

    ratios = [
        r.R / r.n
        for r in results
    ]

    print(
        f"mean R/n                  = "
        f"{sum(ratios) / count:.10f}"
    )

    print(
        f"minimum R/n              = "
        f"{min(ratios):.10f}"
    )

    print(
        f"maximum R/n              = "
        f"{max(ratios):.10f}"
    )

    # -------------------------------------------------------------------------
    # EXACT RECOVERY DETAILS
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("RECOVERY QUALITY")
    print("=" * 100)

    successful = sum(1 for r in results if r.recovered)

    unique_after_r3 = 0

    for r in results:
        if r.exact_square_count == 1:
            unique_after_r3 += 1

    print(
        f"anchors correctly recovered = "
        f"{successful}/{count}"
    )

    print(
        f"anchors with exactly one "
        f"perfect-square candidate = "
        f"{unique_after_r3}/{count}"
    )

    print(
        f"anchors needing <=1 exact "
        f"square check after filters = "
        f"{sum(1 for r in results if r.exact_checks <= 1)}/{count}"
    )

    # -------------------------------------------------------------------------
    # DISTRIBUTION OF FINAL FILTER LOAD
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("FINAL MODULAR FILTER LOAD")
    print("=" * 100)

    distribution: dict[int, int] = {}

    for r in results:
        distribution[r.qr3_count] = (
            distribution.get(r.qr3_count, 0) + 1
        )

    for candidates in sorted(distribution):
        print(
            f"r1+r2+r3 candidates = "
            f"{candidates:4d} "
            f"anchors = {distribution[candidates]:4d}"
        )

    # -------------------------------------------------------------------------
    # ANCHOR DETAILS
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("ANCHOR RESULTS")
    print("=" * 100)

    print(
        " ID        p        q "
        "       r1      r2      r3 "
        "        R/n "
        " BASE   MOD4   QR1   QR2   QR3   EXACT"
    )
    print("-" * 100)

    for idx, r in enumerate(results, 1):
        print(
            f"{idx:3d} "
            f"{r.p:8,d} "
            f"{r.q:8,d} "
            f"{r.triple.r1:8,d} "
            f"{r.triple.r2:8,d} "
            f"{r.triple.r3:8,d} "
            f"{r.R / r.n:12.9f} "
            f"{r.baseline_count:6d} "
            f"{r.mod4_count:6d} "
            f"{r.qr1_count:5d} "
            f"{r.qr2_count:5d} "
            f"{r.qr3_count:5d} "
            f"{r.exact_square_count:6d}"
        )

    # -------------------------------------------------------------------------
    # BEST EXAMPLES
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("MOST REDUCED EXAMPLES")
    print("=" * 100)

    ranked = sorted(
        results,
        key=lambda x: (
            x.qr3_count / max(1, x.baseline_count),
            x.qr3_count,
        ),
    )

    for r in ranked[:20]:
        print(
            f"n={r.n:,} "
            f"p={r.p:,} q={r.q:,} "
            f"mods=({r.triple.r1},{r.triple.r2},{r.triple.r3}) "
            f"R={r.R:,} "
            f"R/n={r.R / r.n:.10f} "
            f"base={r.baseline_count} "
            f"mod4={r.mod4_count} "
            f"qr3={r.qr3_count} "
            f"exact={r.exact_square_count} "
            f"recovered={r.recovered}"
        )

    # -------------------------------------------------------------------------
    # WORST EXAMPLES
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("LEAST REDUCED EXAMPLES")
    print("=" * 100)

    ranked_worst = sorted(
        results,
        key=lambda x: (
            -(x.qr3_count / max(1, x.baseline_count)),
            -x.qr3_count,
        ),
    )

    for r in ranked_worst[:20]:
        print(
            f"n={r.n:,} "
            f"p={r.p:,} q={r.q:,} "
            f"mods=({r.triple.r1},{r.triple.r2},{r.triple.r3}) "
            f"R/n={r.R / r.n:.10f} "
            f"base={r.baseline_count} "
            f"mod4={r.mod4_count} "
            f"qr3={r.qr3_count} "
            f"exact={r.exact_square_count}"
        )

    # -------------------------------------------------------------------------
    # KEY TEST
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("KEY TEST")
    print("=" * 100)

    print(
        """
The central question is whether the three close primes can be used
as a practical filtering mechanism.

For each candidate sum s:

    D = s^2 - 4n

must be an exact square.

Before calculating the exact square root, we test:

    D is a quadratic residue mod r1
    D is a quadratic residue mod r2
    D is a quadratic residue mod r3.

Because the true factors satisfy

    D = (p-q)^2,

the true solution must pass all three filters.

The experiment therefore measures:

    candidate s values
        ->
    modularly admissible s values
        ->
    exact factorization candidates.

A genuinely useful modular exploitation would show that the
three close moduli reduce the search space dramatically while
retaining the true factor pair.

The particularly interesting quantity is:

    average final QR candidate count

versus

    average original admissible s count.

The experiment does NOT assume that uniqueness itself constitutes
a factoring algorithm. The relevant question is whether the
modular filters reduce the actual computational search enough to
be useful without already enumerating the factor population.
"""
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
