#!/usr/bin/env python3

import math
import random
import time
from bisect import bisect_left
from collections import defaultdict
from statistics import median, fmean

from sympy import isprime


# ============================================================================
# KAPPA EXPERIMENT 43
# FAST CRT-COMPATIBILITY LOOKUP
#
# GOAL
# ----
# Experiment 41/42 showed that a large generalized-CRT modulus can reduce
# exact divisibility tests by ~176x.
#
# Experiment 42 also showed that the current compatibility computation is
# so expensive that wall-clock time gets WORSE.
#
# This experiment asks:
#
#   Can the N-dependent compatibility query itself be made cheap enough
#   that the exact-test reduction becomes an actual runtime improvement?
#
# METHODS
# -------
#   1. BASELINE
#        ordinary prime scan
#
#   2. CRT_SCAN
#        for each prime p:
#            r = p mod M
#            qres = n * inverse(r) mod M
#        retain p if qres exists
#
#   3. CRT_INVERSE_INDEX
#        precompute inverse residue -> prime
#        so q compatibility becomes dictionary membership
#
#   4. CRT_DIRECT_LOOKUP
#        because M > prime interval, residue classes are expected to be
#        singleton. Use direct residue -> prime mapping.
#
#   5. CRT_SORTED
#        sorted residue representation.
#
# IMPORTANT
# ---------
# This is a search-engineering experiment, not a new mathematical identity.
#
# NO CSV OUTPUT
# NO SKLEARN
# ============================================================================


SEED = 20260815
random.seed(SEED)

PRIME_LOW = 2_000_000
PRIME_HIGH = 4_200_000

TARGET_COUNT = 1000

# The Experiment 41 modulus.
R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47
]

MODULI = [
    r * r + 3
    for r in R_VALUES
]

PREFIX = 7


# ============================================================================
# FIXED LEGACY TARGETS
# ============================================================================

FIXED_TARGETS = [
    (3318013, 4042603),
    (2129167, 3402323),
    (3978749, 2224517),
    (3685051, 4020281),
    (2452649, 2399627),
    (2593039, 2996527),
    (2149859, 2772097),
    (2514401, 2060543),
    (2883973, 2675423),
    (3960137, 2828887),
    (3793241, 3497381),
    (4011353, 2193509),
]


# ============================================================================
# UTILITIES
# ============================================================================

def lcm(a, b):
    return abs(a * b) // math.gcd(a, b)


def lcm_list(values):
    result = 1
    for value in values:
        result = lcm(result, value)
    return result


def combined_modulus(prefix):
    return lcm_list(MODULI[:prefix])


def generate_prime_population():
    return [
        x
        for x in range(PRIME_LOW, PRIME_HIGH)
        if isprime(x)
    ]


def generate_targets(primes, count):
    """
    Mix legacy targets with deterministic random targets.
    All factors lie in the prime population.
    """
    targets = []

    seen = set()

    for p, q in FIXED_TARGETS:
        if len(targets) >= count:
            break

        key = tuple(sorted((p, q)))

        if key in seen:
            continue

        seen.add(key)
        targets.append((p, q))

    while len(targets) < count:
        p, q = random.sample(primes, 2)

        if p == q:
            continue

        key = tuple(sorted((p, q)))

        if key in seen:
            continue

        seen.add(key)
        targets.append((p, q))

    return targets


def order_candidates(candidates, n):
    sqrt_n = math.isqrt(n)
    return sorted(
        candidates,
        key=lambda p: abs(p - sqrt_n)
    )


# ============================================================================
# BASELINE
# ============================================================================

def baseline_search(primes, n):
    start = time.perf_counter()

    sqrt_n = math.isqrt(n)

    ordered = sorted(
        primes,
        key=lambda p: abs(p - sqrt_n)
    )

    tests = 0

    for p in ordered:
        tests += 1

        if n % p != 0:
            continue

        q = n // p

        if q == p:
            continue

        if PRIME_LOW <= q < PRIME_HIGH and isprime(q):
            return {
                "p": min(p, q),
                "q": max(p, q),
                "tests": tests,
                "time": time.perf_counter() - start,
            }

    return {
        "p": None,
        "q": None,
        "tests": tests,
        "time": time.perf_counter() - start,
    }


# ============================================================================
# METHOD 2: CRT_SCAN
# ============================================================================

def build_residue_inverse_scan(primes, M):
    """
    Store:
        p -> p^-1 mod M

    Only units are usable.
    """
    table = []

    start = time.perf_counter()

    for p in primes:
        r = p % M

        if math.gcd(r, M) != 1:
            continue

        inv = pow(r, -1, M)

        table.append((p, inv))

    return table, time.perf_counter() - start


def crt_scan_search(prime_inverse_table, prime_set, n, M):
    """
    Scan all primes and determine whether the compatible q residue exists.
    This is intentionally the "simple but expensive" compatibility method.
    """
    start = time.perf_counter()

    candidates = []

    for p, inv in prime_inverse_table:
        q_res = (n * inv) % M

        if q_res in prime_set:
            candidates.append(p)

    compatibility_time = time.perf_counter() - start

    ordered = order_candidates(candidates, n)

    exact_start = time.perf_counter()
    tests = 0

    for p in ordered:
        tests += 1

        if n % p != 0:
            continue

        q = n // p

        if (
            p != q
            and PRIME_LOW <= q < PRIME_HIGH
            and isprime(q)
            and p * q == n
        ):
            return {
                "p": min(p, q),
                "q": max(p, q),
                "candidate_count": len(candidates),
                "compatibility_time": compatibility_time,
                "exact_time": time.perf_counter() - exact_start,
                "tests": tests,
                "total_time": time.perf_counter() - start,
            }

    return {
        "p": None,
        "q": None,
        "candidate_count": len(candidates),
        "compatibility_time": compatibility_time,
        "exact_time": time.perf_counter() - exact_start,
        "tests": tests,
        "total_time": time.perf_counter() - start,
    }


# ============================================================================
# METHOD 3: DIRECT RESIDUE -> PRIME INDEX
# ============================================================================

def build_residue_to_prime(primes, M):
    """
    Because M > the entire prime search interval in this experiment,
    distinct primes should have distinct residues.

    We explicitly verify that.

    Returns:
        residue_to_prime
        collision count
        preparation time
    """
    start = time.perf_counter()

    residue_to_prime = {}
    collisions = 0

    for p in primes:
        r = p % M

        old = residue_to_prime.get(r)

        if old is not None and old != p:
            collisions += 1

        residue_to_prime[r] = p

    return (
        residue_to_prime,
        collisions,
        time.perf_counter() - start,
    )


def direct_lookup_search(residue_to_prime, inverse_map, n, M):
    """
    For each prime residue r:

        qres = n * r^{-1} mod M

    Then directly look up q prime.

    This still iterates through residues, but eliminates expensive repeated
    primality/population checks.
    """
    start = time.perf_counter()

    candidates = []

    for r, p in residue_to_prime.items():
        inv = inverse_map.get(r)

        if inv is None:
            continue

        q_res = (n * inv) % M

        if q_res in residue_to_prime:
            candidates.append(p)

    compatibility_time = time.perf_counter() - start

    ordered = order_candidates(candidates, n)

    exact_start = time.perf_counter()
    tests = 0

    for p in ordered:
        tests += 1

        if n % p != 0:
            continue

        q = n // p

        if (
            p != q
            and PRIME_LOW <= q < PRIME_HIGH
            and isprime(q)
            and p * q == n
        ):
            return {
                "p": min(p, q),
                "q": max(p, q),
                "candidate_count": len(candidates),
                "compatibility_time": compatibility_time,
                "exact_time": time.perf_counter() - exact_start,
                "tests": tests,
                "total_time": time.perf_counter() - start,
            }

    return {
        "p": None,
        "q": None,
        "candidate_count": len(candidates),
        "compatibility_time": compatibility_time,
        "exact_time": time.perf_counter() - exact_start,
        "tests": tests,
        "total_time": time.perf_counter() - start,
    }


# ============================================================================
# METHOD 4: PRECOMPUTED N-INDEPENDENT TRANSFORM INDEX
# ============================================================================

def build_inverse_residue_index(primes, M):
    """
    Construct:

        inverse_residue -> prime

    and

        residue -> inverse_residue

    The important point is that all p-dependent arithmetic is performed once.
    """
    start = time.perf_counter()

    inverse_to_prime = {}
    residue_to_inverse = {}
    usable = 0

    for p in primes:
        r = p % M

        if math.gcd(r, M) != 1:
            continue

        inv = pow(r, -1, M)

        inverse_to_prime[inv] = p
        residue_to_inverse[r] = inv
        usable += 1

    return (
        inverse_to_prime,
        residue_to_inverse,
        usable,
        time.perf_counter() - start,
    )


def inverse_index_search(
    inverse_to_prime,
    residue_to_inverse,
    prime_by_residue,
    n,
    M,
):
    """
    Core trick:

        q = n * p^{-1} mod M

    Let a = p^{-1}.

    Then:

        q_res = n*a mod M.

    We need q_res to correspond to an actual prime.

    This avoids calculating p^{-1} inside the target loop.
    """
    start = time.perf_counter()

    candidates = []

    for inv, p in inverse_to_prime.items():
        q_res = (n * inv) % M

        q = prime_by_residue.get(q_res)

        if q is not None:
            candidates.append(p)

    compatibility_time = time.perf_counter() - start

    ordered = order_candidates(candidates, n)

    exact_start = time.perf_counter()
    tests = 0

    for p in ordered:
        tests += 1

        if n % p != 0:
            continue

        q = n // p

        if (
            p != q
            and PRIME_LOW <= q < PRIME_HIGH
            and isprime(q)
            and p * q == n
        ):
            return {
                "p": min(p, q),
                "q": max(p, q),
                "candidate_count": len(candidates),
                "compatibility_time": compatibility_time,
                "exact_time": time.perf_counter() - exact_start,
                "tests": tests,
                "total_time": time.perf_counter() - start,
            }

    return {
        "p": None,
        "q": None,
        "candidate_count": len(candidates),
        "compatibility_time": compatibility_time,
        "exact_time": time.perf_counter() - exact_start,
        "tests": tests,
        "total_time": time.perf_counter() - start,
    }


# ============================================================================
# METHOD 5: SORTED RESIDUE REPRESENTATION
# ============================================================================

def build_sorted_residue_index(primes, M):
    """
    Sorted representation of usable residue/inverse pairs.
    """
    start = time.perf_counter()

    rows = []

    for p in primes:
        r = p % M

        if math.gcd(r, M) != 1:
            continue

        inv = pow(r, -1, M)
        rows.append((inv, p))

    rows.sort()

    inverses = [x[0] for x in rows]
    primes_for_inverse = [x[1] for x in rows]

    return (
        inverses,
        primes_for_inverse,
        time.perf_counter() - start,
    )


def sorted_index_search(
    inverses,
    primes_for_inverse,
    residue_to_prime,
    n,
    M,
):
    """
    Sorted inverse index + binary search for compatible q residue.
    """
    start = time.perf_counter()

    candidates = []

    for inv, p in zip(inverses, primes_for_inverse):
        q_res = (n * inv) % M

        # Binary search the residue dictionary through a sorted-key list
        # constructed externally.
        q = residue_to_prime.get(q_res)

        if q is not None:
            candidates.append(p)

    compatibility_time = time.perf_counter() - start

    ordered = order_candidates(candidates, n)

    exact_start = time.perf_counter()
    tests = 0

    for p in ordered:
        tests += 1

        if n % p != 0:
            continue

        q = n // p

        if (
            p != q
            and PRIME_LOW <= q < PRIME_HIGH
            and isprime(q)
            and p * q == n
        ):
            return {
                "p": min(p, q),
                "q": max(p, q),
                "candidate_count": len(candidates),
                "compatibility_time": compatibility_time,
                "exact_time": time.perf_counter() - exact_start,
                "tests": tests,
                "total_time": time.perf_counter() - start,
            }

    return {
        "p": None,
        "q": None,
        "candidate_count": len(candidates),
        "compatibility_time": compatibility_time,
        "exact_time": time.perf_counter() - exact_start,
        "tests": tests,
        "total_time": time.perf_counter() - start,
    }


# ============================================================================
# BENCHMARK HELPERS
# ============================================================================

def median_or_zero(values):
    return median(values) if values else 0.0


def summarize(name, results, baseline_results):
    total = len(results)

    correct = 0

    candidates = []
    tests = []
    compat = []
    exact = []
    total_times = []

    for result, base in zip(results, baseline_results):
        if result["p"] is not None and result["q"] is not None:
            if (
                result["p"] == base["p"]
                and result["q"] == base["q"]
            ):
                correct += 1

        candidates.append(result["candidate_count"])
        tests.append(result["tests"])
        compat.append(result["compatibility_time"])
        exact.append(result["exact_time"])
        total_times.append(result["total_time"])

    baseline_time = sum(
        x["time"]
        for x in baseline_results
    )

    method_time = sum(total_times)

    return {
        "name": name,
        "correct": correct,
        "total": total,
        "candidate_median": median_or_zero(candidates),
        "tests_median": median_or_zero(tests),
        "compat_mean": fmean(compat),
        "exact_mean": fmean(exact),
        "total_mean": fmean(total_times),
        "batch_time": method_time,
        "speed_ratio": (
            baseline_time / method_time
            if method_time > 0
            else float("inf")
        ),
    }


# ============================================================================
# BATCH BENCHMARK
# ============================================================================

def benchmark_batch(
    targets,
    primes,
    baseline_results,
    indexes,
    label,
):
    print()
    print("=" * 78)
    print(f"BATCH BENCHMARK: {label}")
    print("=" * 78)

    (
        residue_to_prime,
        residue_inverse_map,
        inverse_to_prime,
        sorted_inverses,
        sorted_primes,
        prime_set,
        inverse_table,
        M,
    ) = indexes

    results_scan = []
    results_inverse = []
    results_direct = []
    results_sorted = []

    for idx, (p_true, q_true) in enumerate(targets):

        n = p_true * q_true

        # ------------------------------------------------------------
        # CRT_SCAN
        # ------------------------------------------------------------

        scan_result = crt_scan_search(
            inverse_table,
            prime_set,
            n,
            M,
        )

        results_scan.append(scan_result)

        # ------------------------------------------------------------
        # CRT_INVERSE_INDEX
        # ------------------------------------------------------------

        inverse_result = inverse_index_search(
            inverse_to_prime,
            residue_inverse_map,
            residue_to_prime,
            n,
            M,
        )

        results_inverse.append(inverse_result)

        # ------------------------------------------------------------
        # CRT_DIRECT_LOOKUP
        # ------------------------------------------------------------

        direct_result = direct_lookup_search(
            residue_to_prime,
            residue_inverse_map,
            n,
            M,
        )

        results_direct.append(direct_result)

        # ------------------------------------------------------------
        # CRT_SORTED
        # ------------------------------------------------------------

        sorted_result = sorted_index_search(
            sorted_inverses,
            sorted_primes,
            residue_to_prime,
            n,
            M,
        )

        results_sorted.append(sorted_result)

        if (idx + 1) % max(1, len(targets) // 10) == 0:
            print(
                f"processed {idx + 1:5d}/{len(targets)}"
            )

    summaries = [
        summarize(
            "CRT_SCAN",
            results_scan,
            baseline_results,
        ),
        summarize(
            "CRT_INVERSE_INDEX",
            results_inverse,
            baseline_results,
        ),
        summarize(
            "CRT_DIRECT_LOOKUP",
            results_direct,
            baseline_results,
        ),
        summarize(
            "CRT_SORTED",
            results_sorted,
            baseline_results,
        ),
    ]

    print()
    print(
        f"{'method':<22s}"
        f"{'correct':>10s}"
        f"{'cand med':>14s}"
        f"{'tests med':>14s}"
        f"{'compat':>14s}"
        f"{'exact':>14s}"
        f"{'total':>14s}"
        f"{'speed':>12s}"
    )

    print("-" * 110)

    for s in summaries:
        print(
            f"{s['name']:<22s}"
            f"{s['correct']:>6d}/{s['total']:<3d}"
            f"{s['candidate_median']:>14,.1f}"
            f"{s['tests_median']:>14,.1f}"
            f"{s['compat_mean']:>14.6f}"
            f"{s['exact_mean']:>14.6f}"
            f"{s['total_mean']:>14.6f}"
            f"{s['speed_ratio']:>11.3f}x"
        )

    return summaries


# ============================================================================
# MAIN
# ============================================================================

def main():

    experiment_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 43")
    print("FAST CRT-COMPATIBILITY LOOKUP")
    print("PRECOMPUTED RESIDUE / INVERSE INDEX")
    print("AMORTIZED MANY-TARGET BENCHMARK")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    # ------------------------------------------------------------------------
    # MODULUS
    # ------------------------------------------------------------------------

    M = combined_modulus(PREFIX)

    print()
    print("1. MODULUS")
    print("-" * 78)
    print(f"R values           = {R_VALUES}")
    print(f"moduli             = {MODULI}")
    print(f"prefix             = {PREFIX}")
    print(f"M                  = {M:,}")
    print(f"M bits             = {M.bit_length():d}")
    print(
        f"prime interval width = "
        f"{PRIME_HIGH - PRIME_LOW:,}"
    )

    print()

    # ------------------------------------------------------------------------
    # PRIME POPULATION
    # ------------------------------------------------------------------------

    print("2. PRIME POPULATION")
    print("-" * 78)

    start = time.perf_counter()

    primes = generate_prime_population()

    prime_generation_time = time.perf_counter() - start

    print(
        f"population         = {len(primes):,}"
    )
    print(
        f"generation time    = {prime_generation_time:.6f}s"
    )

    # ------------------------------------------------------------------------
    # TARGETS
    # ------------------------------------------------------------------------

    print()
    print("3. TARGETS")
    print("-" * 78)

    targets = generate_targets(
        primes,
        TARGET_COUNT,
    )

    print(
        f"targets            = {len(targets):,}"
    )

    for i, (p, q) in enumerate(
        targets[:12],
        1,
    ):
        print(
            f"target {i:3d}: "
            f"p={p} "
            f"q={q} "
            f"n={p*q}"
        )

    # ------------------------------------------------------------------------
    # BASELINE
    # ------------------------------------------------------------------------

    print()
    print("4. BASELINE")
    print("-" * 78)

    baseline_results = []

    start = time.perf_counter()

    for i, (p, q) in enumerate(
        targets,
        1,
    ):
        result = baseline_search(
            primes,
            p * q,
        )

        baseline_results.append(result)

        if i % max(1, TARGET_COUNT // 10) == 0:
            print(
                f"baseline {i:5d}/{TARGET_COUNT}"
            )

    baseline_total = (
        time.perf_counter() - start
    )

    baseline_tests = [
        x["tests"]
        for x in baseline_results
    ]

    baseline_times = [
        x["time"]
        for x in baseline_results
    ]

    print()
    print(
        f"baseline batch time      = "
        f"{baseline_total:.6f}s"
    )
    print(
        f"baseline mean/target     = "
        f"{fmean(baseline_times):.6f}s"
    )
    print(
        f"baseline median/target   = "
        f"{median(baseline_times):.6f}s"
    )
    print(
        f"baseline median tests    = "
        f"{median(baseline_tests):,.1f}"
    )

    # ------------------------------------------------------------------------
    # INDEX BUILD
    # ------------------------------------------------------------------------

    print()
    print("5. INDEX CONSTRUCTION")
    print("-" * 78)

    # residue -> prime
    residue_to_prime, collisions, residue_time = (
        build_residue_to_prime(
            primes,
            M,
        )
    )

    # prime -> inverse
    inverse_table, inverse_table_time = (
        build_residue_inverse_scan(
            primes,
            M,
        )
    )

    inverse_to_prime, residue_inverse_map, usable, inverse_time = (
        build_inverse_residue_index(
            primes,
            M,
        )
    )

    (
        sorted_inverses,
        sorted_primes,
        sorted_time,
    ) = build_sorted_residue_index(
        primes,
        M,
    )

    prime_set = set(primes)

    print(
        f"residue->prime preparation = "
        f"{residue_time:.6f}s"
    )
    print(
        f"p->inverse preparation      = "
        f"{inverse_table_time:.6f}s"
    )
    print(
        f"inverse index preparation   = "
        f"{inverse_time:.6f}s"
    )
    print(
        f"sorted index preparation    = "
        f"{sorted_time:.6f}s"
    )

    print(
        f"residue classes             = "
        f"{len(residue_to_prime):,}"
    )
    print(
        f"residue collisions          = "
        f"{collisions:,}"
    )
    print(
        f"usable inverse entries      = "
        f"{usable:,}"
    )

    if collisions == 0:
        print(
            "singleton residue property = PASS"
        )
    else:
        print(
            "singleton residue property = FAIL"
        )

    # ------------------------------------------------------------------------
    # VERIFY INDEX CORRECTNESS
    # ------------------------------------------------------------------------

    print()
    print("6. INDEX CORRECTNESS")
    print("-" * 78)

    index_failures = 0

    sample = random.sample(
        primes,
        min(5000, len(primes)),
    )

    for p in sample:
        r = p % M

        if residue_to_prime.get(r) != p:
            index_failures += 1
            continue

        if math.gcd(r, M) == 1:
            inv = pow(r, -1, M)

            if inverse_to_prime.get(inv) != p:
                index_failures += 1

    print(
        f"sampled primes           = {len(sample):,}"
    )
    print(
        f"index failures            = {index_failures}"
    )
    print(
        f"status                    = "
        f"{'PASS' if index_failures == 0 else 'FAIL'}"
    )

    if index_failures:
        raise RuntimeError(
            "Index validation failed."
        )

    indexes = (
        residue_to_prime,
        residue_inverse_map,
        inverse_to_prime,
        sorted_inverses,
        sorted_primes,
        prime_set,
        inverse_table,
        M,
    )

    # ------------------------------------------------------------------------
    # BATCH BENCHMARKS
    # ------------------------------------------------------------------------

    print()
    print("7. MANY-TARGET BENCHMARK")
    print("-" * 78)

    summaries = benchmark_batch(
        targets,
        primes,
        baseline_results,
        indexes,
        f"{TARGET_COUNT} targets",
    )

    # ------------------------------------------------------------------------
    # BREAKDOWN
    # ------------------------------------------------------------------------

    print()
    print("8. COST BREAKDOWN")
    print("-" * 78)

    total_index_prep = (
        residue_time
        + inverse_table_time
        + inverse_time
        + sorted_time
    )

    print(
        f"one-time index construction = "
        f"{total_index_prep:.6f}s"
    )

    for s in summaries:
        print()
        print(s["name"])

        print(
            f"  mean compatibility = "
            f"{s['compat_mean']:.6f}s"
        )

        print(
            f"  mean exact search  = "
            f"{s['exact_mean']:.6f}s"
        )

        print(
            f"  mean total         = "
            f"{s['total_mean']:.6f}s"
        )

        print(
            f"  speed vs baseline  = "
            f"{s['speed_ratio']:.6f}x"
        )

    # ------------------------------------------------------------------------
    # AMORTIZATION
    # ------------------------------------------------------------------------

    print()
    print("9. AMORTIZATION")
    print("-" * 78)

    baseline_mean = fmean(
        x["time"]
        for x in baseline_results
    )

    for s in summaries:

        per_target = s["total_mean"]

        if per_target < baseline_mean:
            break_even = (
                total_index_prep
                / (baseline_mean - per_target)
            )
        else:
            break_even = math.inf

        if math.isinf(break_even):
            text = "NEVER"
        else:
            text = f"{break_even:.2f} targets"

        print(
            f"{s['name']:<22s} "
            f"break-even = {text}"
        )

    # ------------------------------------------------------------------------
    # TEST-COUNT REDUCTION
    # ------------------------------------------------------------------------

    print()
    print("10. STRUCTURAL RESULT")
    print("-" * 78)

    base_median = median(
        x["tests"]
        for x in baseline_results
    )

    print(
        f"baseline median exact tests = "
        f"{base_median:,.1f}"
    )

    for s in summaries:

        if s["tests_median"] > 0:
            reduction = (
                base_median / s["tests_median"]
            )
        else:
            reduction = float("inf")

        print(
            f"{s['name']:<22s} "
            f"median exact tests="
            f"{s['tests_median']:,.1f} "
            f"reduction={reduction:.3f}x"
        )

    # ------------------------------------------------------------------------
    # INTERPRETATION
    # ------------------------------------------------------------------------

    print()
    print("11. FINAL INTERPRETATION")
    print("-" * 78)

    print(
        """
The experiment distinguishes three layers:

    A. structural filtering
       How many exact n % p operations disappear?

    B. compatibility computation
       How expensive is it to discover those surviving p values?

    C. total factor-search time
       Does the complete method actually beat ordinary scanning?

Experiment 41/42 already established A.

Experiment 43 is primarily testing B.

A successful result therefore requires:

    very large exact-test reduction
    +
    cheap N-dependent compatibility
    +
    lower total runtime

The most important case is CRT_INVERSE_INDEX.

For prefix 7:

    M = 106,261,428

which is larger than the complete prime-search interval.

Therefore a successful implementation may exploit the near-singleton
residue property:

    p -> p mod M

without constructing candidate p,q pairs.

If compatibility remains expensive despite all precomputation,
then the research result is still mathematically useful:

    CRT filtering gives a huge reduction in exact factor tests,

but the missing ingredient is a faster way to enumerate the
compatible residue classes.

No claim of practical factorization speedup is made unless the
measured total wall time actually beats the baseline.
"""
    )

    # ------------------------------------------------------------------------
    # FINAL
    # ------------------------------------------------------------------------

    total_time = (
        time.perf_counter()
        - experiment_start
    )

    print()
    print("=" * 78)
    print("EXPERIMENT 43 COMPLETE")
    print("=" * 78)
    print(
        f"total runtime = "
        f"{total_time:.6f}s "
        f"({total_time/60:.2f} min)"
    )


if __name__ == "__main__":
    main()

