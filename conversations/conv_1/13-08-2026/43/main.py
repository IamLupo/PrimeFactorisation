#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 43
ZERO-COPY CRT FILTER / NO-SORT SEARCH BENCHMARK
NO CSV OUTPUT
==============================================================================

Goal
----
Experiment 42 showed:

    prefix 7
        ~178x fewer divisibility tests
        but still slower wall-clock time

This experiment isolates the implementation overhead.

It deliberately avoids:

    * candidate-pair construction
    * candidate-prime lists
    * sorting candidate primes
    * repeated CRT reconstruction
    * repeated modular inverse computation

Instead:

    1. Build one global prime ordering by distance from sqrt(n).
    2. Precompute each prime's residue class and multiplicative inverse
       modulo the generalized CRT modulus.
    3. For each target, stream through primes exactly once.
    4. Test the CRT compatibility condition using precomputed values.
    5. Perform n % p only on compatible primes.
    6. Stop immediately when the factor is found.

This is intended to answer:

    "Can the large reduction in divisibility tests become a real
     runtime advantage if CRT filtering is implemented without
     constructing or sorting candidate lists?"
==============================================================================
"""

import math
import random
import statistics
import time

import sympy


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 20260814
random.seed(SEED)

PRIME_LOW = 2_000_000
PRIME_HIGH = 4_200_000

TARGETS = 100

R_VALUES = [
    2, 3, 5, 7, 11,
    13, 17, 19, 23, 29,
    31, 37, 41, 43, 47
]

# The useful prefixes from the previous experiment.
PREFIXES = [6, 7]

# Keep the same target-generation model used throughout the experiments.
# We retain the first 12 control targets exactly, then generate additional
# deterministic targets from the same prime population.
FIXED_TARGETS = [
    (3815363, 2494883),
    (2647429, 2840951),
    (2561969, 2902357),
    (2563679, 2283637),
    (3914461, 2161871),
    (3915257, 2498137),
    (2796503, 4003861),
    (2542229, 2958773),
    (3459413, 2601979),
    (2199217, 2834411),
    (3974689, 3730781),
    (2335969, 3463513),
]


# ============================================================================
# MODULI
# ============================================================================

MODULI = [
    r * r + 3
    for r in R_VALUES
]


# ============================================================================
# BASIC HELPERS
# ============================================================================

def generate_primes():
    """
    Generate the complete prime population in the configured interval.
    Uses sympy.isprime as requested.
    """
    primes = []

    for n in range(PRIME_LOW | 1, PRIME_HIGH, 2):
        if sympy.isprime(n):
            primes.append(n)

    return primes


def lcm_list(values):
    result = 1

    for value in values:
        result = math.lcm(result, value)

    return result


def modular_inverse(a, m):
    """
    Return a^{-1} mod m.

    All primes in this experiment are coprime to the selected M.
    """
    return pow(a, -1, m)


def target_pair_is_valid(p, q):
    return (
        PRIME_LOW <= p < PRIME_HIGH
        and PRIME_LOW <= q < PRIME_HIGH
        and sympy.isprime(p)
        and sympy.isprime(q)
        and p != q
    )


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets(primes, count):
    """
    Generate deterministic target factor pairs.

    The first 12 reproduce the established control set.

    Remaining targets are generated from the prime population without
    duplicates.
    """

    result = []

    seen = set()

    for p, q in FIXED_TARGETS:
        pair = tuple(sorted((p, q)))

        if pair not in seen:
            result.append(pair)
            seen.add(pair)

        if len(result) >= count:
            return result[:count]

    while len(result) < count:

        p = random.choice(primes)
        q = random.choice(primes)

        if p == q:
            continue

        pair = tuple(sorted((p, q)))

        if pair in seen:
            continue

        seen.add(pair)
        result.append(pair)

    return result


# ============================================================================
# GENERALIZED CRT PREFIX INFORMATION
# ============================================================================

def build_prefix_data():
    """
    Build generalized-CRT data.

    The selected moduli are NOT pairwise coprime, therefore:

        M = lcm(moduli)

    is the correct combined modulus.
    """

    result = {}

    for prefix in PREFIXES:

        selected = MODULI[:prefix]
        M = lcm_list(selected)

        result[prefix] = {
            "moduli": selected,
            "M": M,
        }

    return result


# ============================================================================
# PRIME ORDERING
# ============================================================================

def build_closest_order(primes):
    """
    One global ordering.

    For targets in this experiment:

        n = p*q

    and the relevant factor is generally near sqrt(n).

    We cannot use one fixed ordering for every n because sqrt(n) changes,
    but we can precompute sorted orders lazily per distinct target sqrt.

    The important part of this experiment is that NO candidate list is
    constructed after CRT filtering.
    """

    return primes


def order_for_target(primes, n):
    """
    Return prime population ordered by distance to sqrt(n).

    This is intentionally built once per target.

    It is NOT a CRT candidate list.

    The same ordering is used by both baseline and CRT-filtered searches.
    """

    root = math.isqrt(n)

    # Python's sort returns a new list. That is acceptable here because
    # this is the common search order shared by both methods.
    return sorted(
        primes,
        key=lambda p: abs(p - root)
    )


# ============================================================================
# PRIME RESIDUE / INVERSE INDEX
# ============================================================================

def build_prefix_index(primes, M):
    """
    Build a compact residue index for a generalized CRT modulus.

    For every occupied residue class r = p mod M we store:

        inverse[r] = r^{-1} mod M

    Since p is prime in our interval and every tested p is coprime to M,
    the inverse exists.

    We also store the set of occupied prime residue classes.
    """

    inverse_by_residue = {}
    occupied = set()

    for p in primes:
        r = p % M

        if r not in inverse_by_residue:
            inverse_by_residue[r] = pow(r, -1, M)

        occupied.add(r)

    return {
        "inverse": inverse_by_residue,
        "occupied": occupied,
        "classes": len(occupied),
    }


# ============================================================================
# TARGET CRT VALUE
# ============================================================================

def compatible_q_residue(n_mod_M, p_residue, inverse_by_residue, M):
    """
    N-only relation:

        p*q = n (mod M)

    for invertible p gives:

        q = n*p^{-1} (mod M)

    The inverse is already cached.
    """

    inv = inverse_by_residue[p_residue]

    return (n_mod_M * inv) % M


# ============================================================================
# BASELINE SEARCH
# ============================================================================

def ordinary_scan(ordered_primes, n, control_pair):
    """
    Ordinary prime scan.

    Every tested prime performs:

        n % p

    until the correct factor is found.
    """

    start = time.perf_counter()

    tests = 0
    recovered = None

    for p in ordered_primes:

        if p * p > n:
            continue

        tests += 1

        if n % p == 0:
            q = n // p

            if (
                target_pair_is_valid(p, q)
                and p * q == n
            ):
                recovered = tuple(sorted((p, q)))
                break

    elapsed = time.perf_counter() - start

    return {
        "tests": tests,
        "time": elapsed,
        "correct": recovered == tuple(sorted(control_pair)),
        "pair": recovered,
    }


# ============================================================================
# ZERO-COPY CRT SEARCH
# ============================================================================

def crt_stream_scan(
    ordered_primes,
    n,
    control_pair,
    prefix_data,
    prefix_index,
):
    """
    Zero-copy CRT-filtered search.

    IMPORTANT:

    There is no candidate list.

    There is no sorting after filtering.

    There is no repeated modular inverse.

    The prime stream remains in the same global closest-first order.

    For each p:

        1. compute p mod M
        2. derive q mod M using cached inverse
        3. test whether that residue belongs to the prime population
        4. only then perform n % p

    This directly measures the computational idea we want to test.
    """

    M = prefix_data["M"]

    inverse_by_residue = prefix_index["inverse"]
    occupied = prefix_index["occupied"]

    n_mod_M = n % M

    tests = 0
    residue_checks = 0
    residue_hits = 0

    recovered = None

    start = time.perf_counter()

    for p in ordered_primes:

        if p * p > n:
            continue

        residue_checks += 1

        rp = p % M

        inv = inverse_by_residue.get(rp)

        if inv is None:
            continue

        rq = (n_mod_M * inv) % M

        if rq not in occupied:
            continue

        residue_hits += 1

        tests += 1

        if n % p == 0:

            q = n // p

            if (
                target_pair_is_valid(p, q)
                and p * q == n
            ):
                recovered = tuple(sorted((p, q)))
                break

    elapsed = time.perf_counter() - start

    return {
        "tests": tests,
        "residue_checks": residue_checks,
        "residue_hits": residue_hits,
        "time": elapsed,
        "correct": recovered == tuple(sorted(control_pair)),
        "pair": recovered,
    }


# ============================================================================
# NO-SORT BASELINE VARIANT
# ============================================================================

def unsorted_prime_scan(primes, n, control_pair):
    """
    Additional benchmark.

    This does NOT sort by sqrt distance.

    It scans the original prime population order.

    This gives a lower-level comparison against the cost of ordering.
    """

    start = time.perf_counter()

    tests = 0
    recovered = None

    for p in primes:

        if p * p > n:
            continue

        tests += 1

        if n % p == 0:

            q = n // p

            if (
                target_pair_is_valid(p, q)
                and p * q == n
            ):
                recovered = tuple(sorted((p, q)))
                break

    elapsed = time.perf_counter() - start

    return {
        "tests": tests,
        "time": elapsed,
        "correct": recovered == tuple(sorted(control_pair)),
        "pair": recovered,
    }


# ============================================================================
# FORMATTING
# ============================================================================

def median(values):
    if not values:
        return 0.0

    return statistics.median(values)


def format_ratio(x):
    if x == 0:
        return "inf"

    return f"{x:.3f}x"


# ============================================================================
# MAIN
# ============================================================================

def main():

    experiment_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 43")
    print("ZERO-COPY CRT FILTER / NO-SORT SEARCH BENCHMARK")
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
    # MODULI
    # ------------------------------------------------------------------------

    print("-" * 78)
    print("1. MODULUS INVENTORY")
    print("-" * 78)

    for r, m in zip(R_VALUES, MODULI):

        print(
            f"r={r:3d} "
            f"m={m:9,d}"
        )

    # ------------------------------------------------------------------------
    # PRIME POPULATION
    # ------------------------------------------------------------------------

    print()
    print("-" * 78)
    print("2. PRIME POPULATION")
    print("-" * 78)

    start = time.perf_counter()

    primes = generate_primes()

    population_time = time.perf_counter() - start

    print(f"prime population = {len(primes):,}")
    print(f"generation time  = {population_time:.6f}s")

    # ------------------------------------------------------------------------
    # TARGETS
    # ------------------------------------------------------------------------

    print()
    print("-" * 78)
    print("3. TARGETS")
    print("-" * 78)

    targets = generate_targets(primes, TARGETS)

    for i, (p, q) in enumerate(targets, start=1):

        n = p * q

        print(
            f"target {i:3d}: "
            f"p={p} "
            f"q={q} "
            f"n={n}"
        )

    # ------------------------------------------------------------------------
    # PREFIX DATA
    # ------------------------------------------------------------------------

    prefix_data = build_prefix_data()

    print()
    print("-" * 78)
    print("4. GENERALIZED CRT PREFIXES")
    print("-" * 78)

    for prefix in PREFIXES:

        data = prefix_data[prefix]

        print(
            f"prefix={prefix:2d} "
            f"M={data['M']:12,d} "
            f"bits={data['M'].bit_length():2d} "
            f"moduli={data['moduli']}"
        )

    # ------------------------------------------------------------------------
    # BASELINE PRIME ORDER PREPARATION
    # ------------------------------------------------------------------------

    print()
    print("-" * 78)
    print("5. BASELINE ORDER PREPARATION")
    print("-" * 78)

    print(
        "The same closest-to-sqrt prime order is shared by the "
        "ordinary and CRT searches."
    )

    # ------------------------------------------------------------------------
    # INDEX CONSTRUCTION
    # ------------------------------------------------------------------------

    indexes = {}

    print()
    print("-" * 78)
    print("6. ZERO-COPY CRT INDEX CONSTRUCTION")
    print("-" * 78)

    total_index_time = 0.0

    for prefix in PREFIXES:

        data = prefix_data[prefix]

        start = time.perf_counter()

        index = build_prefix_index(
            primes,
            data["M"],
        )

        elapsed = time.perf_counter() - start

        indexes[prefix] = index
        total_index_time += elapsed

        print(
            f"prefix={prefix:2d} "
            f"M={data['M']:12,d} "
            f"occupied={index['classes']:8,d} "
            f"index_time={elapsed:.6f}s"
        )

    print()
    print(
        f"total index construction = "
        f"{total_index_time:.6f}s"
    )

    # ------------------------------------------------------------------------
    # ORDINARY BASELINE
    # ------------------------------------------------------------------------

    print()
    print("-" * 78)
    print("7. ORDINARY CLOSEST-PRIME SCAN")
    print("-" * 78)

    baseline_results = []

    for i, (p, q) in enumerate(targets, start=1):

        n = p * q

        ordered = order_for_target(
            primes,
            n,
        )

        result = ordinary_scan(
            ordered,
            n,
            (p, q),
        )

        baseline_results.append(result)

        print(
            f"target {i:3d}: "
            f"tests={result['tests']:8,d} "
            f"time={result['time']:.6f}s "
            f"correct={result['correct']}"
        )

    baseline_total = sum(
        x["time"]
        for x in baseline_results
    )

    baseline_median_tests = median(
        [x["tests"] for x in baseline_results]
    )

    baseline_median_time = median(
        [x["time"] for x in baseline_results]
    )

    print()
    print(
        f"baseline batch time = "
        f"{baseline_total:.6f}s"
    )

    # ------------------------------------------------------------------------
    # PREFIX BENCHMARKS
    # ------------------------------------------------------------------------

    all_summary = {}

    for prefix in PREFIXES:

        print()
        print("=" * 78)
        print(f"8. PREFIX {prefix} ZERO-COPY CRT SEARCH")
        print("=" * 78)

        results = []

        batch_start = time.perf_counter()

        for i, (p, q) in enumerate(targets, start=1):

            n = p * q

            ordered = order_for_target(
                primes,
                n,
            )

            result = crt_stream_scan(
                ordered,
                n,
                (p, q),
                prefix_data[prefix],
                indexes[prefix],
            )

            results.append(result)

            print(
                f"target {i:3d}: "
                f"residue_checks={result['residue_checks']:8,d} "
                f"residue_hits={result['residue_hits']:8,d} "
                f"tests={result['tests']:8,d} "
                f"time={result['time']:.6f}s "
                f"correct={result['correct']}"
            )

        batch_elapsed = time.perf_counter() - batch_start

        total_time = sum(
            x["time"]
            for x in results
        )

        med_tests = median(
            [x["tests"] for x in results]
        )

        med_residue_hits = median(
            [x["residue_hits"] for x in results]
        )

        med_checks = median(
            [x["residue_checks"] for x in results]
        )

        med_time = median(
            [x["time"] for x in results]
        )

        recovered = sum(
            1
            for x in results
            if x["correct"]
        )

        all_summary[prefix] = {
            "results": results,
            "batch_time": batch_elapsed,
            "total_time": total_time,
            "median_tests": med_tests,
            "median_hits": med_residue_hits,
            "median_checks": med_checks,
            "median_time": med_time,
            "recovered": recovered,
        }

        print()
        print(
            f"prefix {prefix}: "
            f"batch={batch_elapsed:.6f}s "
            f"median_tests={med_tests:,.1f} "
            f"median_time={med_time:.6f}s "
            f"recovered={recovered}/{TARGETS}"
        )

    # ------------------------------------------------------------------------
    # UNSORTED COMPARISON
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. UNSORTED PRIME-ORDER CONTROL")
    print("=" * 78)

    unsorted_results = []

    for i, (p, q) in enumerate(targets, start=1):

        n = p * q

        result = unsorted_prime_scan(
            primes,
            n,
            (p, q),
        )

        unsorted_results.append(result)

        if i <= 12:

            print(
                f"target {i:3d}: "
                f"tests={result['tests']:8,d} "
                f"time={result['time']:.6f}s "
                f"correct={result['correct']}"
            )

    unsorted_total = sum(
        x["time"]
        for x in unsorted_results
    )

    unsorted_median_tests = median(
        [x["tests"] for x in unsorted_results]
    )

    unsorted_median_time = median(
        [x["time"] for x in unsorted_results]
    )

    print()
    print(
        f"unsorted batch time = "
        f"{unsorted_total:.6f}s"
    )

    # ------------------------------------------------------------------------
    # FINAL SUMMARY
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. FINAL COMPARISON")
    print("=" * 78)

    print()
    print(
        f"{'method':30s}"
        f"{'median tests':>16s}"
        f"{'median time':>16s}"
        f"{'batch time':>16s}"
        f"{'speed vs baseline':>20s}"
    )

    print("-" * 78)

    print(
        f"{'ordinary closest-prime scan':30s}"
        f"{baseline_median_tests:16,.1f}"
        f"{baseline_median_time:16.6f}"
        f"{baseline_total:16.6f}"
        f"{1.000:20.3f}x"
    )

    print(
        f"{'ordinary unsorted scan':30s}"
        f"{unsorted_median_tests:16,.1f}"
        f"{unsorted_median_time:16.6f}"
        f"{unsorted_total:16.6f}"
        f"{baseline_total / unsorted_total:20.3f}x"
    )

    for prefix in PREFIXES:

        summary = all_summary[prefix]

        amortized_total = (
            total_index_time
            + summary["batch_time"]
        )

        speed = (
            baseline_total / amortized_total
            if amortized_total > 0
            else float("inf")
        )

        print(
            f"{'CRT zero-copy prefix ' + str(prefix):30s}"
            f"{summary['median_tests']:16,.1f}"
            f"{summary['median_time']:16.6f}"
            f"{amortized_total:16.6f}"
            f"{speed:20.3f}x"
        )

    # ------------------------------------------------------------------------
    # TEST REDUCTION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("11. DIVISIBILITY-TEST REDUCTION")
    print("=" * 78)

    for prefix in PREFIXES:

        summary = all_summary[prefix]

        if summary["median_tests"] == 0:
            ratio = float("inf")
        else:
            ratio = (
                baseline_median_tests
                / summary["median_tests"]
            )

        print(
            f"prefix {prefix}: "
            f"{baseline_median_tests:,.1f} -> "
            f"{summary['median_tests']:,.1f} "
            f"tests "
            f"reduction={ratio:.3f}x"
        )

    # ------------------------------------------------------------------------
    # INDEX AMORTIZATION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("12. INDEX AMORTIZATION")
    print("=" * 78)

    for prefix in PREFIXES:

        summary = all_summary[prefix]

        per_target_baseline = (
            baseline_total / TARGETS
        )

        per_target_crt = (
            summary["batch_time"] / TARGETS
        )

        if per_target_crt < per_target_baseline:

            saving = (
                per_target_baseline
                - per_target_crt
            )

            if saving > 0:
                breakeven = (
                    total_index_time / saving
                )
            else:
                breakeven = float("inf")

            print(
                f"prefix {prefix}: "
                f"potential break-even "
                f"after ~{breakeven:.1f} targets"
            )

        else:

            print(
                f"prefix {prefix}: "
                f"no runtime break-even at current "
                f"implementation"
            )

    # ------------------------------------------------------------------------
    # CORRECTNESS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("13. CORRECTNESS")
    print("=" * 78)

    baseline_correct = sum(
        1
        for x in baseline_results
        if x["correct"]
    )

    print(
        f"ordinary closest scan: "
        f"{baseline_correct}/{TARGETS}"
    )

    print(
        f"ordinary unsorted scan: "
        f"{sum(x['correct'] for x in unsorted_results)}/{TARGETS}"
    )

    for prefix in PREFIXES:

        print(
            f"CRT prefix {prefix}: "
            f"{all_summary[prefix]['recovered']}/{TARGETS}"
        )

    # ------------------------------------------------------------------------
    # IMPORTANT DIAGNOSTIC
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("14. DIAGNOSTIC")
    print("=" * 78)

    print(
        """
This experiment is intentionally different from Experiment 42.

Experiment 42 already demonstrated a large reduction in
divisibility tests.

Experiment 43 removes the main implementation overheads:

    * no candidate-prime list
    * no candidate sorting
    * no repeated inverse calculation
    * no pair enumeration
    * one shared closest-prime ordering
    * one reusable residue/inverse index

The key comparison is therefore:

    fewer n % p operations
        versus
    cost of the residue-class membership calculation.

If the CRT method is still slower, then the remaining gap is
mostly intrinsic to the Python-level filtering operation.

If prefix 7 becomes faster, then the earlier slowdown was primarily
an implementation artifact rather than a mathematical limitation.

A correct recovery for every target remains a necessary control.

No CSV files are produced.
"""
    )

    # ------------------------------------------------------------------------
    # DONE
    # ------------------------------------------------------------------------

    total_runtime = (
        time.perf_counter()
        - experiment_start
    )

    print()
    print("=" * 78)
    print("EXPERIMENT 43 COMPLETE")
    print("=" * 78)
    print(
        f"total runtime = "
        f"{total_runtime:.6f}s"
    )


if __name__ == "__main__":
    main()

