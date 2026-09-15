#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 45
INCREMENTAL RESIDUE-SET INTERSECTION / PRIME CANDIDATE SHRINKAGE
NO CSV OUTPUT
==============================================================================

Purpose
-------
Measure how the searchable prime set shrinks as the modular constraints are
applied one modulus at a time.

This is different from the earlier pair-count experiments.

We work with candidate prime factors p <= sqrt(n).

For each modulus m:

    p*q == n (mod m)

and, because all candidate primes are coprime to the selected moduli,

    q == n * p^(-1) (mod m).

A prime p survives a modulus exactly when its required q-residue is occupied
by at least one prime in the searchable population.

The experiment measures:

    all primes <= sqrt(n)
        ->
    after m=7
        ->
    after m=12
        ->
    after m=28
        ->
    ...
        ->
    after final modulus

Then exact divisibility tests are performed only on the final survivors.

No pair enumeration.
No giant CRT class construction.
No CSV files.
Uses sympy.isprime.
"""

from __future__ import annotations

import math
import random
import statistics
import time
from collections import defaultdict
from typing import Dict, List, Sequence, Tuple

from sympy import isprime


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 20260814

TARGETS = 100

PRIME_LOW = 2_000_000
PRIME_HIGH = 4_200_000

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19,
    23, 29, 31, 37, 41, 43, 47
]

PREFIXES = [3, 4, 5, 6, 7]

# How many targets get full stage-by-stage output.
DETAIL_TARGETS = 12


# ============================================================================
# MODULUS FAMILY
# ============================================================================

def modulus_for_r(r: int) -> int:
    """
    Experimental modulus family used throughout the recent experiments.

        m = r^2 + 3
    """
    return r * r + 3


def build_moduli() -> List[int]:
    return [modulus_for_r(r) for r in R_VALUES]


def lcm(a: int, b: int) -> int:
    return a // math.gcd(a, b) * b


def lcm_many(values: Sequence[int]) -> int:
    result = 1
    for value in values:
        result = lcm(result, value)
    return result


# ============================================================================
# PRIME GENERATION
# ============================================================================

def generate_primes(low: int, high: int) -> List[int]:
    """
    Generate all primes in [low, high).

    Uses sympy.isprime as requested.
    """
    primes: List[int] = []

    for x in range(low, high):
        if isprime(x):
            primes.append(x)

    return primes


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets(
    primes: Sequence[int],
    count: int,
    rng: random.Random,
) -> List[Tuple[int, int, int]]:
    """
    Generate distinct semiprime targets using primes from the population.
    """
    targets: List[Tuple[int, int, int]] = []
    seen_n = set()

    n_primes = len(primes)

    while len(targets) < count:
        p = primes[rng.randrange(n_primes)]
        q = primes[rng.randrange(n_primes)]

        if p == q:
            continue

        n = p * q

        if n in seen_n:
            continue

        seen_n.add(n)
        targets.append((p, q, n))

    return targets


# ============================================================================
# RESIDUE OCCUPANCY
# ============================================================================

def build_residue_occupancy(
    primes: Sequence[int],
    moduli: Sequence[int],
) -> Dict[int, set[int]]:
    """
    For every modulus m, store the residues occupied by at least one prime.

    Example:

        occupancy[28] = {1, 3, 5, ...}
    """
    occupancy: Dict[int, set[int]] = {}

    for m in moduli:
        residues = {p % m for p in primes}
        occupancy[m] = residues

    return occupancy


# ============================================================================
# SEARCH ORDER
# ============================================================================

def closest_prime_order(
    primes: Sequence[int],
    sqrt_n: float,
) -> List[int]:
    """
    Return primes <= sqrt(n), ordered by increasing distance from sqrt(n).

    The ordering is shared between the baseline and filtered search.
    """
    candidates = [p for p in primes if p * p <= int(n_safety_bound(sqrt_n))]
    candidates.sort(key=lambda p: (abs(sqrt_n - p), -p))
    return candidates


def n_safety_bound(x: float) -> float:
    """
    Small helper to avoid accidental floating-point truncation issues.
    """
    return math.floor(x * x)


# ============================================================================
# COMPATIBILITY MAP
# ============================================================================

def build_required_residue_map(
    n: int,
    m: int,
    occupied_q_residues: set[int],
) -> Dict[int, bool]:
    """
    For every residue r = p mod m, determine whether the required q residue

        q == n * p^(-1) (mod m)

    is occupied by a prime.

    Returns:

        compatible[r] = True/False
    """
    n_mod = n % m

    compatible: Dict[int, bool] = {}

    for p_residue in range(m):
        if math.gcd(p_residue, m) != 1:
            compatible[p_residue] = False
            continue

        inverse = pow(p_residue, -1, m)
        required_q = (n_mod * inverse) % m

        compatible[p_residue] = required_q in occupied_q_residues

    return compatible


# ============================================================================
# INCREMENTAL FILTER
# ============================================================================

def incremental_filter(
    ordered_candidates: Sequence[int],
    n: int,
    moduli: Sequence[int],
    occupancy: Dict[int, set[int]],
) -> Tuple[List[List[int]], List[float]]:
    """
    Apply the modular constraints incrementally.

    Returns:
        stage_candidates:
            [initial, after m1, after m2, ...]
        stage_times:
            elapsed time for each filtering stage
    """
    current = list(ordered_candidates)

    stage_candidates: List[List[int]] = [current]
    stage_times: List[float] = []

    for m in moduli:
        t0 = time.perf_counter()

        compatible = build_required_residue_map(
            n,
            m,
            occupancy[m],
        )

        # Preserve the existing search order.
        current = [
            p
            for p in current
            if compatible[p % m]
        ]

        elapsed = time.perf_counter() - t0

        stage_candidates.append(current)
        stage_times.append(elapsed)

    return stage_candidates, stage_times


# ============================================================================
# EXACT RECOVERY
# ============================================================================

def exact_recovery(
    candidates: Sequence[int],
    n: int,
) -> Tuple[bool, int, Tuple[int, int] | None]:
    """
    Perform exact n % p tests only on the final survivors.

    Returns:
        correct,
        number_of_divisibility_tests,
        recovered_pair
    """
    tests = 0

    for p in candidates:
        tests += 1

        if n % p == 0:
            q = n // p

            if p <= q:
                pair = (p, q)
            else:
                pair = (q, p)

            return True, tests, pair

    return False, tests, None


# ============================================================================
# BASELINE
# ============================================================================

def baseline_scan(
    ordered_candidates: Sequence[int],
    n: int,
) -> Tuple[bool, int, Tuple[int, int] | None]:
    """
    Ordinary closest-prime scan.
    """
    tests = 0

    for p in ordered_candidates:
        tests += 1

        if n % p == 0:
            q = n // p

            if p <= q:
                pair = (p, q)
            else:
                pair = (q, p)

            return True, tests, pair

    return False, tests, None


# ============================================================================
# STATISTICS
# ============================================================================

def median(values: Sequence[float]) -> float:
    return statistics.median(values) if values else 0.0


def fmt_int(x: int) -> str:
    return f"{x:,}"


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    total_start = time.perf_counter()

    rng = random.Random(SEED)

    moduli = build_moduli()

    print("=" * 78)
    print("KAPPA EXPERIMENT 45")
    print("INCREMENTAL RESIDUE-SET INTERSECTION / PRIME CANDIDATE SHRINKAGE")
    print("NO CSV OUTPUT")
    print("=" * 78)
    print()
    print(f"random seed      = {SEED}")
    print(f"targets          = {TARGETS}")
    print(
        f"prime interval   = "
        f"[{PRIME_LOW:,}, {PRIME_HIGH:,})"
    )
    print(f"prefixes         = {PREFIXES}")
    print(f"R values         = {R_VALUES}")
    print()

    # ------------------------------------------------------------------------
    # 1. MODULUS INVENTORY
    # ------------------------------------------------------------------------

    print("-" * 78)
    print("1. MODULUS INVENTORY")
    print("-" * 78)

    for r, m in zip(R_VALUES, moduli):
        print(
            f"r={r:>3} "
            f"m={m:>10,} "
            f"gcd(previous structure handled by LCM)"
        )

    print()

    # ------------------------------------------------------------------------
    # 2. PRIME POPULATION
    # ------------------------------------------------------------------------

    print("-" * 78)
    print("2. PRIME POPULATION")
    print("-" * 78)

    t0 = time.perf_counter()

    primes = generate_primes(PRIME_LOW, PRIME_HIGH)

    prime_generation_time = time.perf_counter() - t0

    print(f"prime population = {fmt_int(len(primes))}")
    print(f"generation time  = {prime_generation_time:.6f}s")
    print()

    if len(primes) < 2:
        raise RuntimeError("Prime population is too small.")

    # ------------------------------------------------------------------------
    # 3. TARGETS
    # ------------------------------------------------------------------------

    print("-" * 78)
    print("3. TARGETS")
    print("-" * 78)

    targets = generate_targets(primes, TARGETS, rng)

    for i, (p, q, n) in enumerate(targets, start=1):
        lo, hi = sorted((p, q))
        print(
            f"target {i:>3}: "
            f"p={lo} q={hi} n={n}"
        )

    print()

    # ------------------------------------------------------------------------
    # 4. RESIDUE OCCUPANCY
    # ------------------------------------------------------------------------

    print("-" * 78)
    print("4. RESIDUE OCCUPANCY INDEX")
    print("-" * 78)

    t0 = time.perf_counter()

    occupancy = build_residue_occupancy(primes, moduli)

    occupancy_time = time.perf_counter() - t0

    for m in moduli:
        print(
            f"m={m:>10,} "
            f"occupied residues={len(occupancy[m]):>8,}"
        )

    print(f"index build time = {occupancy_time:.6f}s")
    print()

    # ------------------------------------------------------------------------
    # 5. BASELINE
    # ------------------------------------------------------------------------

    print("-" * 78)
    print("5. ORDINARY CLOSEST-PRIME BASELINE")
    print("-" * 78)

    baseline_tests: List[int] = []
    baseline_times: List[float] = []

    baseline_results = []

    for i, (_, _, n) in enumerate(targets, start=1):
        sqrt_n = math.isqrt(n)
        ordered = [
            p
            for p in primes
            if p <= sqrt_n
        ]

        ordered.sort(
            key=lambda p: (sqrt_n - p, -p)
        )

        t0 = time.perf_counter()

        ok, tests, pair = baseline_scan(ordered, n)

        elapsed = time.perf_counter() - t0

        baseline_tests.append(tests)
        baseline_times.append(elapsed)
        baseline_results.append((ok, pair))

        print(
            f"target {i:>3}: "
            f"tests={tests:>7,} "
            f"time={elapsed:.6f}s "
            f"correct={ok}"
        )

    print()
    print(
        f"baseline total time = "
        f"{sum(baseline_times):.6f}s"
    )
    print(
        f"baseline median tests = "
        f"{median(baseline_tests):,.1f}"
    )
    print()

    # ------------------------------------------------------------------------
    # 6. INCREMENTAL FILTERING
    # ------------------------------------------------------------------------

    print("-" * 78)
    print("6. INCREMENTAL RESIDUE-SET INTERSECTION")
    print("-" * 78)

    all_stage_counts: Dict[int, List[int]] = {
        prefix: []
        for prefix in PREFIXES
    }

    all_stage_reductions: Dict[int, List[float]] = {
        prefix: []
        for prefix in PREFIXES
    }

    all_exact_tests: Dict[int, List[int]] = {
        prefix: []
        for prefix in PREFIXES
    }

    all_filter_times: Dict[int, List[float]] = {
        prefix: []
        for prefix in PREFIXES
    }

    all_total_times: Dict[int, List[float]] = {
        prefix: []
        for prefix in PREFIXES
    }

    all_recovery: Dict[int, int] = {
        prefix: 0
        for prefix in PREFIXES
    }

    for target_index, (_, _, n) in enumerate(targets, start=1):

        sqrt_n = math.isqrt(n)

        ordered_candidates = [
            p
            for p in primes
            if p <= sqrt_n
        ]

        ordered_candidates.sort(
            key=lambda p: (sqrt_n - p, -p)
        )

        if target_index <= DETAIL_TARGETS:
            print()
            print("-" * 78)
            print(
                f"TARGET {target_index:>3} "
                f"n={n}"
            )
            print("-" * 78)

        current = ordered_candidates

        target_stage_rows = []

        for stage_index, m in enumerate(moduli[:max(PREFIXES)], start=1):

            t0 = time.perf_counter()

            compatible = build_required_residue_map(
                n,
                m,
                occupancy[m],
            )

            previous_count = len(current)

            current = [
                p
                for p in current
                if compatible[p % m]
            ]

            elapsed = time.perf_counter() - t0

            new_count = len(current)

            if previous_count:
                reduction = previous_count / new_count
                retained_fraction = new_count / previous_count
            else:
                reduction = float("inf")
                retained_fraction = 0.0

            target_stage_rows.append(
                (
                    stage_index,
                    m,
                    previous_count,
                    new_count,
                    retained_fraction,
                    reduction,
                    elapsed,
                )
            )

            # Store exact counts for the prefix corresponding to this stage.
            if stage_index in PREFIXES:
                prefix = stage_index

                all_stage_counts[prefix].append(new_count)
                all_stage_reductions[prefix].append(reduction)
                all_filter_times[prefix].append(elapsed)

        # ------------------------------------------------------------
        # Detailed output for first DETAIL_TARGETS
        # ------------------------------------------------------------

        if target_index <= DETAIL_TARGETS:

            print()
            print(
                "stage   modulus   before   after   retained      "
                "reduction   time"
            )
            print("-" * 78)

            for (
                stage_index,
                m,
                before,
                after,
                retained_fraction,
                reduction,
                elapsed,
            ) in target_stage_rows:

                print(
                    f"{stage_index:>5} "
                    f"{m:>10,} "
                    f"{before:>8,} "
                    f"{after:>8,} "
                    f"{retained_fraction:>9.6f} "
                    f"{reduction:>10.3f}x "
                    f"{elapsed:.6f}s"
                )

        # ------------------------------------------------------------
        # Exact recovery at every requested prefix
        # ------------------------------------------------------------

        current = ordered_candidates

        for stage_index, m in enumerate(
            moduli[:max(PREFIXES)],
            start=1,
        ):

            compatible = build_required_residue_map(
                n,
                m,
                occupancy[m],
            )

            current = [
                p
                for p in current
                if compatible[p % m]
            ]

            if stage_index in PREFIXES:

                prefix = stage_index

                t0 = time.perf_counter()

                ok, exact_tests, pair = exact_recovery(
                    current,
                    n,
                )

                recovery_time = time.perf_counter() - t0

                all_exact_tests[prefix].append(exact_tests)

                total_search_time = (
                    sum(
                        row[-1]
                        for row in target_stage_rows
                        if row[0] <= stage_index
                    )
                    + recovery_time
                )

                all_total_times[prefix].append(
                    total_search_time
                )

                if ok:
                    all_recovery[prefix] += 1

                if target_index <= DETAIL_TARGETS:
                    print(
                        f"prefix {prefix}: "
                        f"final candidates={len(current):>8,} "
                        f"exact tests={exact_tests:>6,} "
                        f"recovered={ok}"
                    )

    # ------------------------------------------------------------------------
    # 7. GLOBAL STAGE SUMMARY
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. GLOBAL STAGE SUMMARY")
    print("=" * 78)

    print()
    print(
        "prefix   modulus   median survivors   "
        "median retained   median reduction"
    )
    print("-" * 78)

    for prefix in PREFIXES:
        m = moduli[prefix - 1]

        counts = all_stage_counts[prefix]
        reductions = all_stage_reductions[prefix]

        previous_counts = []
        for i in range(TARGETS):
            if prefix == 1:
                previous_counts.append(
                    len(
                        [
                            p
                            for p in primes
                            if p * p <= targets[i][2]
                        ]
                    )
                )
            else:
                # Reconstructing previous median is unnecessary;
                # the important quantity is the per-stage retention.
                pass

        median_count = median(counts)

        median_retained = 1.0 / median(reductions)

        print(
            f"{prefix:>5} "
            f"{m:>10,} "
            f"{median_count:>18,.1f} "
            f"{median_retained:>17.6f} "
            f"{median(reductions):>15.3f}x"
        )

    # ------------------------------------------------------------------------
    # 8. EXACT-TEST SUMMARY
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. EXACT DIVISIBILITY-TEST SUMMARY")
    print("=" * 78)

    print()
    print(
        "prefix   median exact tests   "
        "baseline median tests   "
        "test reduction   recovered"
    )
    print("-" * 78)

    baseline_median = median(baseline_tests)

    for prefix in PREFIXES:
        filtered_median = median(all_exact_tests[prefix])

        if filtered_median:
            reduction = baseline_median / filtered_median
        else:
            reduction = float("inf")

        print(
            f"{prefix:>5} "
            f"{filtered_median:>20,.1f} "
            f"{baseline_median:>22,.1f} "
            f"{reduction:>14.3f}x "
            f"{all_recovery[prefix]:>5}/{TARGETS}"
        )

    # ------------------------------------------------------------------------
    # 9. FILTER COST
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. FILTER COST")
    print("=" * 78)

    print()
    print(
        "prefix   median filter time   "
        "median exact-test time   median total"
    )
    print("-" * 78)

    for prefix in PREFIXES:
        filter_time = median(all_filter_times[prefix])
        total_time = median(all_total_times[prefix])

        # Exact-test time is not stored separately in all targets,
        # so estimate the non-filter component from total.
        exact_component = max(0.0, total_time - filter_time)

        print(
            f"{prefix:>5} "
            f"{filter_time:>19.6f}s "
            f"{exact_component:>22.6f}s "
            f"{total_time:>14.6f}s"
        )

    # ------------------------------------------------------------------------
    # 10. IMPORTANT CROSS-CHECK
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. CROSS-CHECK")
    print("=" * 78)

    print()
    print(
        "The final incremental survivor set should be exactly the same"
    )
    print(
        "as applying all selected congruence constraints simultaneously."
    )
    print()
    print(
        "This experiment therefore measures the shrinkage path:"
    )
    print()
    print(
        "P0"
        " -> P1"
        " -> P2"
        " -> ..."
        " -> Pk"
    )
    print()
    print(
        "without constructing a giant CRT candidate-pair population."
    )

    # ------------------------------------------------------------------------
    # 11. INTERPRETATION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("11. INTERPRETATION")
    print("=" * 78)

    print()
    print(
        "The main quantity is the incremental survivor ratio after each"
    )
    print(
        "new modulus."
    )

    print()
    print(
        "A particularly useful result would be:"
    )
    print()
    print(
        "  early moduli: little or no reduction"
    )
    print(
        "  middle moduli: substantial reduction"
    )
    print(
        "  final moduli: very strong reduction"
    )

    print()
    print(
        "That would show exactly where the useful information enters."
    )

    print()
    print(
        "The experiment also separates three costs:"
    )
    print(
        "  1. residue filtering"
    )
    print(
        "  2. exact divisibility testing"
    )
    print(
        "  3. one-time residue-index construction"
    )

    print()
    print(
        "The desired outcome is not merely a small final candidate set."
    )
    print(
        "The important question is whether the cumulative filtering"
    )
    print(
        "cost is smaller than the exact divisibility work it removes."
    )

    print()
    print(
        "No candidate-pair enumeration is performed."
    )
    print(
        "No CSV files are produced."
    )

    # ------------------------------------------------------------------------
    # 12. TOTAL RUNTIME
    # ------------------------------------------------------------------------

    total_time = time.perf_counter() - total_start

    print()
    print("=" * 78)
    print("EXPERIMENT 45 COMPLETE")
    print("=" * 78)
    print(
        f"total runtime = {total_time:.6f}s "
        f"({total_time / 60:.2f} min)"
    )


if __name__ == "__main__":
    main()

