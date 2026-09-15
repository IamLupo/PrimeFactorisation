#!/usr/bin/env python3

"""
====================================================================================================
THREE-CLOSE-PRIME HYPERBOLA / BATCH-DIVISOR ISOLATION EXPERIMENT
====================================================================================================

FOLLOW-UP TO:
    THREE-CLOSE-PRIME QUOTIENT-CELL DEFECT / ONE-REMAINDER ELIMINATION

IMPORTANT OBSERVATION
---------------------

The previous identity was:

    a = (D-k*r1*c)/(l*r2+c)

with:

    D = n-k*l*r1*r2.

But:

    D-k*r1*c
        = n-k*r1(l*r2+c)

and therefore:

    a = n/q-k*r1

where:

    q = l*r2+c.

Thus:

    exact integer division
        <=> 
    q divides n.

So the ~4.28 "division hits" per anchor were not a new
factorization condition. They were simply divisors of n.

This experiment therefore removes the individual divisor test.

For each hyperbola strip we have an integer interval:

    q_lo <= q <= q_hi

Any true factor q in this interval satisfies:

    q | n.

Instead of testing every q independently, we construct:

    P = product(q_lo ... q_hi) mod n

and compute:

    gcd(P,n).

If the gcd is 1:

    no candidate in the interval divides n.

If the gcd is nontrivial:

    at least one candidate in the interval divides n.

If the gcd is n:

    both factors may be contained in that interval,
    so the interval is recursively split.

The recursion isolates actual divisors without requiring one
individual modulo operation per candidate.

The experiment measures:

    integer q candidates
    prime q candidates
    q values actually multiplied
    batch gcd calls
    leaf intervals
    exact factor recoveries

This is an algorithmic experiment rather than a claim of a new
factoring algorithm.

The main question is:

    can batch isolation reduce the arithmetic work associated
    with testing ~100,000 q values to roughly O(number of useful
    intervals) gcd operations?

We also compare against the ordinary direct divisor test.

====================================================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from itertools import combinations


# ================================================================================================
# CONFIGURATION
# ================================================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

ANCHORS = 300

MODULUS_MIN = 300
MODULUS_MAX = 3_000

CLOSE_RATIO = 0.20

SEED = 1_511_464_998

PROGRESS_STEP = 25

# Maximum number of values in one product batch.
BATCH_LIMIT = 512


# ================================================================================================
# DATA
# ================================================================================================

@dataclass(frozen=True)
class Anchor:
    p: int
    q: int
    n: int


@dataclass
class Result:
    anchor_id: int

    p: int
    q: int
    n: int

    r1: int
    r2: int
    r3: int

    R: int
    gap: int

    q_candidates: int
    prime_q_candidates: int

    direct_divisor_tests: int

    batch_values: int
    batch_gcd_calls: int
    batch_leaf_intervals: int

    exact_solutions: int
    recovered: bool


# ================================================================================================
# PRIME SIEVE
# ================================================================================================

def sieve(limit: int) -> list[int]:

    if limit < 2:
        return []

    flags = bytearray(
        b"\x01"
    ) * (limit + 1)

    flags[0:2] = b"\x00\x00"

    root = math.isqrt(limit)

    for p in range(
        2,
        root + 1,
    ):

        if not flags[p]:
            continue

        start = p*p

        count = (
            (limit-start)//p
            + 1
        )

        flags[
            start:
            limit+1:
            p
        ] = b"\x00" * count

    return [
        i
        for i,flag in enumerate(flags)
        if flag
    ]


# ================================================================================================
# ANCHORS
# ================================================================================================

def build_anchors(
    primes: list[int],
    count: int,
    seed: int,
) -> list[Anchor]:

    rng = random.Random(seed)

    result = []

    while len(result) < count:

        p,q = rng.sample(
            primes,
            2,
        )

        if p > q:
            p,q = q,p

        result.append(
            Anchor(
                p=p,
                q=q,
                n=p*q,
            )
        )

    return result


# ================================================================================================
# CLOSE MODULUS TRIPLE
# ================================================================================================

def choose_triple(
    n: int,
    modulus_primes: list[int],
) -> tuple[int,int,int]:

    target = max(
        MODULUS_MIN,
        int(n ** (1.0/3.0)),
    )

    radius = max(
        20,
        int(
            target*CLOSE_RATIO
        ),
    )

    candidates = [
        r
        for r in modulus_primes
        if abs(r-target) <= radius
    ]

    if len(candidates) < 3:
        raise RuntimeError(
            f"not enough close primes for n={n}"
        )

    best = None
    best_score = None

    for r1,r2,r3 in combinations(
        candidates,
        3,
    ):

        R = r1*r2*r3

        if R >= n:
            continue

        score = (
            R/n,
            -(r3-r1),
        )

        if (
            best_score is None
            or score > best_score
        ):

            best_score = score
            best = (
                r1,
                r2,
                r3,
            )

    if best is None:
        raise RuntimeError(
            f"could not find triple for n={n}"
        )

    return best


# ================================================================================================
# HYPERBOLA Q INTERVAL
# ================================================================================================

def q_interval_for_k(
    n: int,
    r1: int,
    k: int,
) -> tuple[int,int] | None:

    """
    p lies in:

        k*r1 <= p <= (k+1)*r1-1

    so:

        floor(n/p_hi) <= q <= floor(n/p_lo)

    We additionally clamp to the factor interval.
    """

    p_lo = max(
        FACTOR_MIN,
        k*r1,
    )

    p_hi = min(
        FACTOR_MAX,
        (k+1)*r1-1,
    )

    if p_lo > p_hi:
        return None

    # Since q = n/p decreases as p increases:

    q_lo = (
        n//p_hi
    )

    q_hi = (
        n//p_lo
    )

    q_lo = max(
        FACTOR_MIN,
        q_lo,
    )

    q_hi = min(
        FACTOR_MAX,
        q_hi,
    )

    if q_lo > q_hi:
        return None

    return q_lo,q_hi


# ================================================================================================
# BATCH PRODUCT MOD N
# ================================================================================================

def product_mod(
    n: int,
    lo: int,
    hi: int,
) -> int:

    """
    Return:

        product(lo ... hi) mod n

    This deliberately uses a sequential modular product so that the
    number of multiplicative steps can be measured directly.
    """

    value = 1 % n

    for x in range(
        lo,
        hi+1,
    ):

        value = (
            value*x
        ) % n

    return value


# ================================================================================================
# BATCH GCD RECURSION
# ================================================================================================

@dataclass
class BatchStats:
    values: int = 0
    gcd_calls: int = 0
    leaves: int = 0


def isolate_interval(
    n: int,
    lo: int,
    hi: int,
    stats: BatchStats,
    solutions: set[int],
) -> None:

    """
    Recursively isolate divisors of n contained in [lo,hi].

    Important cases:

        gcd == 1:
            no divisor in interval.

        gcd == n:
            more than one relevant divisor may exist, so split.

        1 < gcd < n:
            a nontrivial divisor was found directly.

    A singleton interval is checked exactly.
    """

    if lo > hi:
        return

    size = hi-lo+1

    # Singleton.
    if size == 1:

        x = lo

        stats.values += 1
        stats.gcd_calls += 1

        if x != 0 and n % x == 0:

            if (
                x > 1
                and x < n
            ):
                solutions.add(x)

        stats.leaves += 1

        return

    # For very large intervals, split into bounded batches first.
    if size > BATCH_LIMIT:

        mid = (
            lo+hi
        )//2

        isolate_interval(
            n,
            lo,
            mid,
            stats,
            solutions,
        )

        isolate_interval(
            n,
            mid+1,
            hi,
            stats,
            solutions,
        )

        return

    # Batch product.
    stats.values += size

    prod = product_mod(
        n,
        lo,
        hi,
    )

    stats.gcd_calls += 1

    g = math.gcd(
        prod,
        n,
    )

    if g == 1:
        return

    if g != n:

        if (
            g > 1
            and g < n
        ):
            solutions.add(g)

        # There may still be another factor in the interval.
        # Continue recursively so that both factors are visible.
        mid = (
            lo+hi
        )//2

        isolate_interval(
            n,
            lo,
            mid,
            stats,
            solutions,
        )

        isolate_interval(
            n,
            mid+1,
            hi,
            stats,
            solutions,
        )

        return

    # g == n.
    #
    # Both factors may be inside the interval, or the product contains
    # multiple copies/factor combinations. Split.

    mid = (
        lo+hi
    )//2

    isolate_interval(
        n,
        lo,
        mid,
        stats,
        solutions,
    )

    isolate_interval(
        n,
        mid+1,
        hi,
        stats,
        solutions,
    )


# ================================================================================================
# DIRECT BASELINE
# ================================================================================================

def direct_factor_scan(
    n: int,
    intervals: list[tuple[int,int]],
) -> tuple[int,set[int]]:

    tests = 0
    factors = set()

    for lo,hi in intervals:

        for q in range(
            lo,
            hi+1,
        ):

            tests += 1

            if q > 1 and n % q == 0:

                if q < n:
                    factors.add(q)

    return tests,factors


# ================================================================================================
# ANCHOR SEARCH
# ================================================================================================

def run_anchor(
    anchor_id: int,
    anchor: Anchor,
    triple: tuple[int,int,int],
) -> Result:

    n = anchor.n

    r1,r2,r3 = triple

    R = (
        r1*r2*r3
    )

    gap = n-R

    # Determine valid k range.

    k_min = (
        FACTOR_MIN//r1
    )

    if k_min*r1 < FACTOR_MIN:
        k_min += 1

    k_max = (
        FACTOR_MAX//r1
    )

    if k_max*r1 > FACTOR_MAX:
        k_max -= 1

    intervals = []

    for k in range(
        max(1,k_min),
        k_max+1,
    ):

        qr = q_interval_for_k(
            n,
            r1,
            k,
        )

        if qr is None:
            continue

        intervals.append(qr)

    # Merge overlapping q intervals.
    #
    # This is deliberately measured because the hyperbola mapping can
    # produce adjacent/overlapping candidate intervals for neighboring k.

    intervals.sort()

    merged = []

    for lo,hi in intervals:

        if not merged:
            merged.append(
                [lo,hi]
            )
            continue

        prev_lo,prev_hi = merged[-1]

        if lo <= prev_hi+1:

            merged[-1][1] = max(
                prev_hi,
                hi,
            )

        else:

            merged.append(
                [lo,hi]
            )

    merged_intervals = [
        (lo,hi)
        for lo,hi in merged
    ]

    q_candidates = sum(
        hi-lo+1
        for lo,hi
        in merged_intervals
    )

    # Prime q candidates for comparison.
    #
    # This is only a diagnostic baseline, not required for recovery.
    local_prime_count = 0

    # Fast enough at this experimental scale.
    for lo,hi in merged_intervals:

        for q in range(
            lo,
            hi+1,
        ):

            if is_prime_small(q):
                local_prime_count += 1

    # Direct baseline.
    direct_tests,direct_factors = (
        direct_factor_scan(
            n,
            merged_intervals,
        )
    )

    # Batch isolation.
    stats = BatchStats()
    batch_factors = set()

    for lo,hi in merged_intervals:

        isolate_interval(
            n,
            lo,
            hi,
            stats,
            batch_factors,
        )

    # Combine.
    all_factors = (
        direct_factors
        |
        batch_factors
    )

    exact_pairs = set()

    for f in all_factors:

        if f <= 1 or f >= n:
            continue

        if n % f != 0:
            continue

        other = n//f

        if (
            FACTOR_MIN <= f <= FACTOR_MAX
            and
            FACTOR_MIN <= other <= FACTOR_MAX
        ):

            exact_pairs.add(
                (
                    min(f,other),
                    max(f,other),
                )
            )

    # Make the anchor pair itself count only when exact.
    recovered = (
        (
            min(anchor.p,anchor.q),
            max(anchor.p,anchor.q),
        )
        in exact_pairs
    )

    return Result(
        anchor_id=anchor_id,
        p=anchor.p,
        q=anchor.q,
        n=n,
        r1=r1,
        r2=r2,
        r3=r3,
        R=R,
        gap=gap,
        q_candidates=q_candidates,
        prime_q_candidates=local_prime_count,
        direct_divisor_tests=direct_tests,
        batch_values=stats.values,
        batch_gcd_calls=stats.gcd_calls,
        batch_leaf_intervals=stats.leaves,
        exact_solutions=len(exact_pairs),
        recovered=recovered,
    )


# ================================================================================================
# SMALL PRIME TEST
# ================================================================================================

def is_prime_small(n: int) -> bool:

    if n < 2:
        return False

    if n % 2 == 0:
        return n == 2

    d = 3

    while d*d <= n:

        if n % d == 0:
            return False

        d += 2

    return True


# ================================================================================================
# CORRELATION
# ================================================================================================

def correlation(
    xs: list[float],
    ys: list[float],
) -> float:

    if len(xs) < 2:
        return float("nan")

    mx = statistics.mean(xs)
    my = statistics.mean(ys)

    dx = [
        x-mx
        for x in xs
    ]

    dy = [
        y-my
        for y in ys
    ]

    den_x = math.sqrt(
        sum(
            x*x
            for x in dx
        )
    )

    den_y = math.sqrt(
        sum(
            y*y
            for y in dy
        )
    )

    if den_x == 0 or den_y == 0:
        return float("nan")

    return (
        sum(
            x*y
            for x,y in zip(
                dx,
                dy,
            )
        )
        /
        (
            den_x*den_y
        )
    )


# ================================================================================================
# MAIN
# ================================================================================================

def run() -> None:

    total_start = time.perf_counter()

    print("="*100)
    print(
        "THREE-CLOSE-PRIME HYPERBOLA / "
        "BATCH-DIVISOR ISOLATION EXPERIMENT"
    )
    print("="*100)

    print(
        f"M                         = "
        f"{M:,}"
    )

    print(
        f"factor range              = "
        f"{FACTOR_MIN:,} - {FACTOR_MAX:,}"
    )

    print(
        f"anchors                   = "
        f"{ANCHORS}"
    )

    print(
        f"modulus prime range       = "
        f"{MODULUS_MIN:,} - {MODULUS_MAX:,}"
    )

    print(
        f"close ratio               = "
        f"{CLOSE_RATIO:.0%}"
    )

    print(
        f"seed                      = "
        f"{SEED:,}"
    )

    # ============================================================================================
    # PRIME POOLS
    # ============================================================================================

    print()
    print("="*100)
    print("BUILDING PRIME POOLS")
    print("="*100)

    all_primes = sieve(
        FACTOR_MAX
    )

    factor_primes = [
        p
        for p in all_primes
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    modulus_primes = [
        p
        for p in sieve(MODULUS_MAX)
        if MODULUS_MIN <= p <= MODULUS_MAX
    ]

    print(
        f"factor primes             = "
        f"{len(factor_primes):,}"
    )

    print(
        f"modulus primes            = "
        f"{len(modulus_primes):,}"
    )

    # ============================================================================================
    # ANCHORS
    # ============================================================================================

    print()
    print("="*100)
    print("BUILDING ANCHORS")
    print("="*100)

    anchors = build_anchors(
        factor_primes,
        ANCHORS,
        SEED,
    )

    print(
        f"actual anchors            = "
        f"{len(anchors):,}"
    )

    # ============================================================================================
    # MODULI
    # ============================================================================================

    print()
    print("="*100)
    print(
        "SELECTING ANCHOR-SPECIFIC CLOSE TRIPLES"
    )
    print("="*100)

    triples = []

    for i,anchor in enumerate(
        anchors,
        1,
    ):

        triple = choose_triple(
            anchor.n,
            modulus_primes,
        )

        triples.append(
            triple
        )

        if (
            i % PROGRESS_STEP == 0
            or i == len(anchors)
        ):

            print(
                f"anchor {i:3d}/{len(anchors)}"
            )

    # ============================================================================================
    # VALIDATION
    # ============================================================================================

    print()
    print("="*100)
    print("VALIDATING MODULUS TRIPLES")
    print("="*100)

    failures = 0

    for anchor,triple in zip(
        anchors,
        triples,
    ):

        r1,r2,r3 = triple

        if not (
            r1 < r2 < r3
        ):
            failures += 1

        if (
            r1*r2*r3
            >= anchor.n
        ):
            failures += 1

    print(
        f"triple validation failures = "
        f"{failures}"
    )

    if failures:
        raise RuntimeError(
            "modulus validation failed"
        )

    print(
        f"unique modulus triples      = "
        f"{len(set(triples))}"
    )

    # ============================================================================================
    # SEARCH
    # ============================================================================================

    print()
    print("="*100)
    print(
        "RUNNING BATCH-DIVISOR ISOLATION SEARCH"
    )
    print("="*100)

    results = []

    search_start = time.perf_counter()

    for i,(anchor,triple) in enumerate(
        zip(
            anchors,
            triples,
        ),
        1,
    ):
        result = run_anchor(
            i,
            anchor,
            triple,
        )

        results.append(
            result
        )

        if (
            i % PROGRESS_STEP == 0
            or i == len(anchors)
        ):

            print(
                f"anchor {i:3d}/{len(anchors)}"
            )

    search_time = (
        time.perf_counter()
        -
        search_start
    )

    # ============================================================================================
    # AVERAGES
    # ============================================================================================

    avg_q = statistics.mean([
        r.q_candidates
        for r in results
    ])

    avg_prime_q = statistics.mean([
        r.prime_q_candidates
        for r in results
    ])

    avg_direct = statistics.mean([
        r.direct_divisor_tests
        for r in results
    ])

    avg_batch_values = statistics.mean([
        r.batch_values
        for r in results
    ])

    avg_gcd_calls = statistics.mean([
        r.batch_gcd_calls
        for r in results
    ])

    avg_leaves = statistics.mean([
        r.batch_leaf_intervals
        for r in results
    ])

    avg_exact = statistics.mean([
        r.exact_solutions
        for r in results
    ])

    # ============================================================================================
    # SUMMARY
    # ============================================================================================

    print()
    print("="*100)
    print("SUMMARY")
    print("="*100)

    print(
        f"anchors analyzed            = "
        f"{len(results)}"
    )

    print(
        f"average q candidates        = "
        f"{avg_q:,.3f}"
    )

    print(
        f"average prime q candidates  = "
        f"{avg_prime_q:,.3f}"
    )

    print(
        f"average direct divisor tests= "
        f"{avg_direct:,.3f}"
    )

    print(
        f"average batch values        = "
        f"{avg_batch_values:,.3f}"
    )

    print(
        f"average batch gcd calls     = "
        f"{avg_gcd_calls:,.3f}"
    )

    print(
        f"average leaf intervals      = "
        f"{avg_leaves:,.3f}"
    )

    print(
        f"average exact solutions     = "
        f"{avg_exact:,.3f}"
    )

    # ============================================================================================
    # BATCH REDUCTION
    # ============================================================================================

    print()
    print("="*100)
    print("BATCH REDUCTION")
    print("="*100)

    print(
        f"batch gcd / direct tests    = "
        f"{avg_gcd_calls/avg_direct:.9f}"
    )

    print(
        f"batch values / direct tests = "
        f"{avg_batch_values/avg_direct:.9f}"
    )

    print(
        f"q candidates / 8,363        = "
        f"{avg_q/8363:.9f}"
    )

    print(
        f"prime q / 8,363            = "
        f"{avg_prime_q/8363:.9f}"
    )

    # ============================================================================================
    # RECOVERY
    # ============================================================================================

    recovered = sum(
        r.recovered
        for r in results
    )

    print()
    print("="*100)
    print("RECOVERY")
    print("="*100)

    print(
        f"correctly recovered         = "
        f"{recovered}/{len(results)}"
    )

    print(
        f"recovery rate               = "
        f"{100*recovered/len(results):.4f}%"
    )

    # ============================================================================================
    # GAP CORRELATIONS
    # ============================================================================================

    gaps = [
        r.gap
        for r in results
    ]

    print()
    print("="*100)
    print("GAP CORRELATION")
    print("="*100)

    print(
        f"mean gap                    = "
        f"{statistics.mean(gaps):,.3f}"
    )

    print(
        f"minimum gap                 = "
        f"{min(gaps):,}"
    )

    print(
        f"maximum gap                 = "
        f"{max(gaps):,}"
    )

    print(
        f"corr(gap,q candidates)      = "
        f"{correlation(
            gaps,
            [
                r.q_candidates
                for r in results
            ],
        ): .6f}"
    )

    print(
        f"corr(gap,batch gcd calls)   = "
        f"{correlation(
            gaps,
            [
                r.batch_gcd_calls
                for r in results
            ],
        ): .6f}"
    )

    # ============================================================================================
    # STRONGEST COLLAPSES
    # ============================================================================================

    print()
    print("="*100)
    print("STRONGEST BATCH COLLAPSES")
    print("="*100)

    strongest = sorted(
        results,
        key=lambda r: (
            r.batch_gcd_calls,
            r.q_candidates,
        ),
    )[:20]

    for r in strongest:

        print(
            f"n={r.n:,} "
            f"p={r.p:,} "
            f"q={r.q:,} "
            f"mods=({r.r1},{r.r2},{r.r3}) "
            f"R/n={r.R/r.n:.10f} "
            f"gap={r.gap:,} "
            f"q={r.q_candidates:,} "
            f"primeQ={r.prime_q_candidates:,} "
            f"gcd={r.batch_gcd_calls:,} "
            f"leaves={r.batch_leaf_intervals:,} "
            f"exact={r.exact_solutions}"
        )

    # ============================================================================================
    # WEAKEST
    # ============================================================================================

    print()
    print("="*100)
    print("WEAKEST BATCH COLLAPSES")
    print("="*100)

    weakest = sorted(
        results,
        key=lambda r: (
            r.batch_gcd_calls,
            r.q_candidates,
        ),
        reverse=True,
    )[:20]

    for r in weakest:

        print(
            f"n={r.n:,} "
            f"p={r.p:,} "
            f"q={r.q:,} "
            f"mods=({r.r1},{r.r2},{r.r3}) "
            f"R/n={r.R/r.n:.10f} "
            f"gap={r.gap:,} "
            f"q={r.q_candidates:,} "
            f"primeQ={r.prime_q_candidates:,} "
            f"gcd={r.batch_gcd_calls:,} "
            f"leaves={r.batch_leaf_intervals:,} "
            f"exact={r.exact_solutions}"
        )

    # ============================================================================================
    # ANCHOR TABLE
    # ============================================================================================

    print()
    print("="*100)
    print("ANCHOR RESULTS")
    print("="*100)

    print(
        f"{'ID':>3} "
        f"{'p':>7} "
        f"{'q':>7} "
        f"{'r1':>5} "
        f"{'r2':>5} "
        f"{'r3':>5} "
        f"{'Q':>8} "
        f"{'PRIMEQ':>8} "
        f"{'DIRECT':>8} "
        f"{'GCD':>6} "
        f"{'LEAF':>6} "
        f"{'EXACT':>5}"
    )

    print("-"*100)

    for r in results:

        print(
            f"{r.anchor_id:3d} "
            f"{r.p:7,d} "
            f"{r.q:7,d} "
            f"{r.r1:5d} "
            f"{r.r2:5d} "
            f"{r.r3:5d} "
            f"{r.q_candidates:8,d} "
            f"{r.prime_q_candidates:8,d} "
            f"{r.direct_divisor_tests:8,d} "
            f"{r.batch_gcd_calls:6,d} "
            f"{r.batch_leaf_intervals:6,d} "
            f"{r.exact_solutions:5,d}"
        )

    # ============================================================================================
    # TIMING
    # ============================================================================================

    total_time = (
        time.perf_counter()
        -
        total_start
    )

    print()
    print("="*100)
    print("TIMING")
    print("="*100)

    print(
        f"search runtime              = "
        f"{search_time:.3f} seconds"
    )

    print(
        f"total runtime               = "
        f"{total_time:.3f} seconds"
    )

    # ============================================================================================
    # INTERPRETATION
    # ============================================================================================

    print()
    print("="*100)
    print("MATHEMATICAL INTERPRETATION")
    print("="*100)

    print(
        """
The previous defect experiment appeared to produce a powerful
condition:

    a = (D-k*r1*c)/(l*r2+c)

with only ~4.3 exact division hits per anchor.

But the identity simplifies:

    D = n-k*l*r1*r2

therefore:

    D-k*r1*c
      = n-k*r1(l*r2+c)

Set:

    q = l*r2+c.

Then:

    a = n/q-k*r1.

Consequently:

    exact division
        <=>
    q | n.

So the previous "division collapse" was simply rediscovering
the divisors of n.

The present experiment asks whether the hyperbola intervals can
make that divisor search computationally sparse.

For each quotient k:

    k*r1 <= p < (k+1)*r1

and therefore:

    n/((k+1)*r1) < q <= n/(k*r1).

This produces an explicit q interval.

Instead of testing:

    n % q

for every q independently, we form:

    P_interval =
        product(q_lo ... q_hi) mod n

and compute:

    gcd(P_interval,n).

If the gcd is 1, the complete interval contains no divisor.

If the gcd is nontrivial, the interval contains at least one divisor
and is recursively split.

Thus:

    q interval
         |
         v
    batch product
         |
         v
       gcd
      /   \
     1   nontrivial
         |
         v
      recursive
        split
         |
         v
      exact factor

The important comparison is:

    direct divisor tests
        versus
    batch gcd calls.

However, there is an equally important second comparison:

    number of q values multiplied.

A reduction in gcd calls alone does NOT mean the underlying
factorization complexity disappeared, because forming the product
still processes all candidate q values.

Therefore the experiment deliberately reports both.

A genuinely interesting result would require:

    batch gcd calls << direct tests

AND ideally:

    batch values << direct tests.

If only the gcd count collapses while the number of multiplied q
values remains essentially unchanged, then this is a batching
optimization, not a new mathematical factorization mechanism.

The experiment also gives a clean negative test:

    if every anchor requires nearly all q values to be incorporated,
    then the hyperbola interval contains no hidden sparse divisor
    structure beyond the fact that exactly one q divides n.

Every recovered candidate is finally checked with:

    p*q == n.

====================================================================================================
EXPERIMENT COMPLETE
====================================================================================================
"""
    )


if __name__ == "__main__":
    run()

