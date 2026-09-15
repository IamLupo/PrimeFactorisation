#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 94
BATCH-GCD SIGMA_1 ACCESS / FACTOR-DISCOVERY ACCELERATION TEST

GOAL
----
Experiment 93 showed that naive divisor accumulation only learns sigma_1(N)
when it reaches a real factor.

This experiment asks a more interesting computational question:

    Can we detect the hidden factor p of N in large batches,
    without testing every d individually?

If p is found, then:

    sigma_1(N) = (p+1)(q+1)

and:

    s = sigma_1(N) - N - 1

so the factorization is immediately recovered.

METHODS
-------
1. TRIAL:
       test d=2,3,... until a factor is found.

2. BATCH-GCD:
       group consecutive integers into blocks;
       build a product tree;
       compute gcd(N, block_product);
       recursively descend only into blocks that contain a divisor.

3. BLOCK MODULO:
       optional lightweight comparison using product modulo N.

The experiment compares:

    candidate operations
    gcd calls
    modular multiplications
    wall time

The key question is whether batch detection reduces the practical cost
of finding the first factor on balanced semiprimes.

IMPORTANT
---------
This is an engineering experiment, not a claim that batch GCD asymptotically
beats factoring.

For the current ~10^13 targets, the factor is around 2-4 million, so this
should run quickly enough to benchmark.

NO CSV
NO SKLEARN
==============================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from typing import Optional, Tuple, List


# ============================================================================
# CONFIG
# ============================================================================

NUM_TARGETS = 20

P_MIN = 2_000_000
P_MAX = 4_200_000

RNG_SEED = 94094

# Block sizes to compare.
BLOCK_SIZES = [
    32,
    128,
    512,
    2048,
    8192,
    32768,
]

PRINT_TARGETS = 12


# ============================================================================
# DATA
# ============================================================================

@dataclass(frozen=True)
class Target:
    p: int
    q: int
    n: int
    s: int


@dataclass
class Stats:
    gcd_calls: int = 0
    modular_ops: int = 0
    leaves_tested: int = 0
    blocks_tested: int = 0
    product_nodes: int = 0


# ============================================================================
# PRIME SIEVE
# ============================================================================

def sieve_primes(
    lo: int,
    hi: int,
) -> List[int]:

    lo = max(
        2,
        lo,
    )

    sieve = bytearray(
        b"\x01" * (hi + 1)
    )

    sieve[0:2] = b"\x00\x00"

    for p in range(
        2,
        math.isqrt(hi) + 1,
    ):
        if sieve[p]:

            start = p * p

            sieve[
                start:hi + 1:p
            ] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [
        x
        for x in range(
            lo,
            hi + 1,
        )
        if sieve[x]
    ]


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets(
    primes: List[int],
    count: int,
) -> List[Target]:

    rng = random.Random(
        RNG_SEED
    )

    out = []
    seen = set()

    while len(out) < count:

        p = primes[
            rng.randrange(
                len(primes)
            )
        ]

        q = primes[
            rng.randrange(
                len(primes)
            )
        ]

        if p == q:
            continue

        if p > q:
            p, q = q, p

        if (
            p,
            q,
        ) in seen:
            continue

        seen.add(
            (
                p,
                q,
            )
        )

        out.append(
            Target(
                p=p,
                q=q,
                n=p * q,
                s=p + q,
            )
        )

    return out


# ============================================================================
# BASELINE TRIAL DIVISION
# ============================================================================

def trial_factor(
    n: int,
    stats: Stats,
) -> Optional[int]:

    limit = math.isqrt(n)

    for d in range(
        2,
        limit + 1,
    ):

        stats.leaves_tested += 1

        if n % d == 0:
            return d

    return None


# ============================================================================
# PRODUCT TREE
# ============================================================================

def product_tree(
    lo: int,
    hi: int,
    n: int,
    stats: Stats,
) -> int:
    """
    Product of integers [lo, hi], reduced modulo n.

    Recursive product tree.

    This is the core batch operation.
    """

    if lo > hi:
        return 1

    if lo == hi:

        stats.product_nodes += 1

        return lo % n

    mid = (
        lo + hi
    ) // 2

    left = product_tree(
        lo,
        mid,
        n,
        stats,
    )

    right = product_tree(
        mid + 1,
        hi,
        n,
        stats,
    )

    stats.modular_ops += 1

    return (
        left * right
    ) % n


# ============================================================================
# FACTOR SEARCH WITH PRODUCT TREE
# ============================================================================

def batch_factor_recursive(
    n: int,
    lo: int,
    hi: int,
    stats: Stats,
) -> Optional[int]:

    """
    Search [lo,hi] for a nontrivial divisor using a product tree.

    At each interval:

        g = gcd(N, product(interval))

    If g=1:
        no divisor in interval.

    If interval is a singleton:
        return that divisor.

    If g=N:
        there may be several divisors in the interval, so descend.

    This is deliberately conservative and handles the possible gcd=N case.
    """

    if lo > hi:
        return None

    stats.blocks_tested += 1

    if lo == hi:

        stats.gcd_calls += 1

        if (
            lo > 1
            and lo < n
            and n % lo == 0
        ):
            return lo

        stats.leaves_tested += 1

        return None

    product = product_tree(
        lo,
        hi,
        n,
        stats,
    )

    stats.gcd_calls += 1

    g = math.gcd(
        n,
        product,
    )

    if g == 1:
        return None

    if (
        g != n
        and g > 1
    ):

        # We already have a nontrivial gcd.
        return g

    mid = (
        lo + hi
    ) // 2

    left = batch_factor_recursive(
        n,
        lo,
        mid,
        stats,
    )

    if left is not None:
        return left

    return batch_factor_recursive(
        n,
        mid + 1,
        hi,
        stats,
    )


def batch_factor(
    n: int,
    block_size: int,
    stats: Stats,
) -> Optional[int]:

    """
    Top-level blocked batch factor search.

    Rather than constructing one giant tree, process chunks.
    """

    limit = math.isqrt(n)

    lo = 2

    while lo <= limit:

        hi = min(
            limit,
            lo + block_size - 1,
        )

        stats.blocks_tested += 1

        # Compute the product modulo N directly.
        product = 1

        for x in range(
            lo,
            hi + 1,
        ):
            product *= x

            # Frequent reduction keeps Python integers bounded.
            product %= n

            stats.modular_ops += 1

        stats.gcd_calls += 1

        g = math.gcd(
            n,
            product,
        )

        if (
            g > 1
            and g < n
        ):
            return g

        if g == n:

            # There is at least one divisor in this block.
            factor = batch_factor_recursive(
                n,
                lo,
                hi,
                stats,
            )

            if factor is not None:
                return factor

        lo = (
            hi + 1
        )

    return None


# ============================================================================
# SIGMA / S RECOVERY
# ============================================================================

def recover_from_factor(
    n: int,
    p: int,
) -> Tuple[int, int, Tuple[int, int]]:

    if p <= 1:
        raise ValueError(
            "invalid factor"
        )

    if n % p != 0:
        raise ValueError(
            "factor does not divide n"
        )

    q = n // p

    sigma1 = (
        p + 1
    ) * (
        q + 1
    )

    s = (
        sigma1
        - n
        - 1
    )

    return (
        sigma1,
        s,
        (
            min(p, q),
            max(p, q),
        ),
    )


# ============================================================================
# BENCHMARK ONE METHOD
# ============================================================================

def benchmark_target(
    target: Target,
    block_size: int,
) -> dict:

    stats = Stats()

    t0 = time.perf_counter()

    factor = batch_factor(
        target.n,
        block_size,
        stats,
    )

    elapsed = (
        time.perf_counter()
        - t0
    )

    if factor is None:
        recovered = False

        sigma1 = None
        s = None
        factors = None

    else:

        (
            sigma1,
            s,
            factors,
        ) = recover_from_factor(
            target.n,
            factor,
        )

        recovered = (
            factors
            == (
                target.p,
                target.q,
            )
        )

    return {
        "factor": factor,
        "elapsed": elapsed,
        "gcd_calls": stats.gcd_calls,
        "modular_ops": stats.modular_ops,
        "leaves_tested": stats.leaves_tested,
        "blocks_tested": stats.blocks_tested,
        "product_nodes": stats.product_nodes,
        "recovered": recovered,
        "sigma1": sigma1,
        "s": s,
        "factors": factors,
    }


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 94")
    print("BATCH-GCD SIGMA_1 ACCESS / FACTOR-DISCOVERY TEST")
    print("TRIAL DIVISION VS BATCHED DIVISIBILITY")
    print("STRICT FACTOR / SIGMA1 / S VERIFICATION")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    # ------------------------------------------------------------------
    # Prime population
    # ------------------------------------------------------------------

    t0 = time.perf_counter()

    primes = sieve_primes(
        P_MIN,
        P_MAX,
    )

    print("\n1. PRIME POPULATION")
    print("-" * 78)

    print(
        f"prime population = "
        f"{len(primes)}"
    )

    print(
        f"generation time = "
        f"{time.perf_counter()-t0:.6f}s"
    )

    # ------------------------------------------------------------------
    # Targets
    # ------------------------------------------------------------------

    targets = generate_targets(
        primes,
        NUM_TARGETS,
    )

    print(
        f"total targets = "
        f"{len(targets)}"
    )

    for i, t in enumerate(
        targets[:PRINT_TARGETS],
        1,
    ):

        print(
            f"target {i:3d}: "
            f"p={t.p} "
            f"q={t.q} "
            f"n={t.n} "
            f"s={t.s}"
        )

    # ------------------------------------------------------------------
    # Baseline
    # ------------------------------------------------------------------

    print("\n2. TRIAL DIVISION BASELINE")
    print("-" * 78)

    baseline_rows = []

    for i, t in enumerate(
        targets,
        1,
    ):

        stats = Stats()

        t0 = time.perf_counter()

        factor = trial_factor(
            t.n,
            stats,
        )

        elapsed = (
            time.perf_counter()
            - t0
        )

        ok = (
            factor is not None
            and (
                factor
                == t.p
            )
        )

        baseline_rows.append(
            {
                "target": t,
                "factor": factor,
                "elapsed": elapsed,
                "tests": stats.leaves_tested,
                "ok": ok,
            }
        )

        print(
            f"target {i:2d}: "
            f"factor={factor} "
            f"tests={stats.leaves_tested} "
            f"time={elapsed:.6f}s"
        )

    # ------------------------------------------------------------------
    # Batch methods
    # ------------------------------------------------------------------

    print("\n3. BATCH-GCD RESULTS")
    print("-" * 78)

    all_results = {}

    for block_size in BLOCK_SIZES:

        print(
            f"\nBLOCK SIZE = "
            f"{block_size}"
        )

        rows = []

        for i, t in enumerate(
            targets,
            1,
        ):

            result = benchmark_target(
                t,
                block_size,
            )

            rows.append(
                result
            )

            print(
                f"target {i:2d}: "
                f"factor={result['factor']} "
                f"blocks={result['blocks_tested']} "
                f"gcds={result['gcd_calls']} "
                f"ops={result['modular_ops']} "
                f"time={result['elapsed']:.6f}s "
                f"ok={result['recovered']}"
            )

        all_results[
            block_size
        ] = rows

    # ------------------------------------------------------------------
    # Aggregate comparison
    # ------------------------------------------------------------------

    print("\n4. AGGREGATE COMPARISON")
    print("-" * 78)

    baseline_time = statistics.fmean(
        row["elapsed"]
        for row in baseline_rows
    )

    baseline_tests = statistics.fmean(
        row["tests"]
        for row in baseline_rows
    )

    print(
        f"TRIAL "
        f"mean_time={baseline_time:.6f}s "
        f"mean_candidates={baseline_tests:.1f}"
    )

    print(
        "BLOCK | mean_time | speedup | "
        "mean_gcds | mean_ops | success"
    )

    for block_size, rows in (
        all_results.items()
    ):

        mean_time = statistics.fmean(
            row["elapsed"]
            for row in rows
        )

        mean_gcds = statistics.fmean(
            row["gcd_calls"]
            for row in rows
        )

        mean_ops = statistics.fmean(
            row["modular_ops"]
            for row in rows
        )

        successes = sum(
            row["recovered"]
            for row in rows
        )

        speedup = (
            baseline_time
            / mean_time
            if mean_time > 0
            else float("inf")
        )

        print(
            f"{block_size:5d} | "
            f"{mean_time:.6f} | "
            f"{speedup:.3f}x | "
            f"{mean_gcds:.1f} | "
            f"{mean_ops:.1f} | "
            f"{successes}/{len(rows)}"
        )

    # ------------------------------------------------------------------
    # Sigma1 and s reconstruction
    # ------------------------------------------------------------------

    print("\n5. SIGMA1 / S RECONSTRUCTION")
    print("-" * 78)

    best_block = min(
        BLOCK_SIZES,
        key=lambda b:
            statistics.fmean(
                r["elapsed"]
                for r in all_results[b]
            ),
    )

    print(
        f"fastest tested block size = "
        f"{best_block}"
    )

    rows = all_results[
        best_block
    ]

    sigma_success = 0
    factor_success = 0

    for i, (
        t,
        row,
    ) in enumerate(
        zip(
            targets,
            rows,
        ),
        1,
    ):

        if not row["recovered"]:
            continue

        sigma_success += (
            row["sigma1"]
            == (
                (t.p + 1)
                * (t.q + 1)
            )
        )

        factor_success += (
            row["factors"]
            == (
                t.p,
                t.q,
            )
        )

        if i <= PRINT_TARGETS:

            print(
                f"target {i:2d}: "
                f"sigma1={row['sigma1']} "
                f"s={row['s']} "
                f"true_s={t.s} "
                f"factors={row['factors']}"
            )

    print(
        f"sigma1 exact checks = "
        f"{sigma_success}/{NUM_TARGETS}"
    )

    print(
        f"factor exact checks = "
        f"{factor_success}/{NUM_TARGETS}"
    )

    # ------------------------------------------------------------------
    # Scaling diagnostic
    # ------------------------------------------------------------------

    print("\n6. SCALING DIAGNOSTIC")
    print("-" * 78)

    print(
        "For each method we compare:"
    )

    print(
        "  wall time"
    )

    print(
        "  number of modular multiplications"
    )

    print(
        "  number of GCD calls"
    )

    print(
        "  number of individual divisor leaves"
    )

    print()
    print(
        "The interesting outcome is NOT merely fewer GCD calls."
    )

    print(
        "It is a reduction in total work needed before the first "
        "factor is discovered."
    )

    # ------------------------------------------------------------------
    # Final diagnostic
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("7. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "Experiment 93 showed that naive partial sigma accumulation "
        "waits for the first factor."
    )

    print()
    print(
        "Experiment 94 asks whether batch arithmetic can bypass the "
        "individual divisor-by-divisor search."
    )

    print()
    print(
        "POSITIVE:"
    )

    print(
        "  batch-GCD finds p substantially faster than trial division."
    )

    print()
    print(
        "STRONG POSITIVE:"
    )

    print(
        "  the speedup survives across block sizes and targets, "
        "while exact sigma1 and factor recovery remain 100%."
    )

    print()
    print(
        "NEGATIVE:"
    )

    print(
        "  batch-GCD merely moves the same ~sqrt(N) work around."
    )

    print()
    print(
        "If negative, we should stop treating sigma1 as an independent "
        "computational route and focus on genuinely sub-sqrt factor "
        "discovery or relation generation."
    )

    print(
        f"\ntotal runtime = "
        f"{time.perf_counter()-start:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 94 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

