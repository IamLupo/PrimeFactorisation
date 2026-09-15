#!/usr/bin/env python3

import math
import random
import time
from statistics import median

from sympy import isprime


# ============================================================================
# KAPPA EXPERIMENT 42
# AMORTIZED CRT FILTER VS ORDINARY PRIME SCAN
# NO CSV OUTPUT
# ============================================================================

SEED = 20260814
random.seed(SEED)

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47
]

PRIME_LOW = 2_000_000
PRIME_HIGH = 4_200_000

# Keep this moderate so the experiment finishes comfortably.
TARGETS = 100

# Test both the useful threshold from Experiment 41
# and the strongest prefix.
PREFIXES = [6, 7]

# Prime generation uses sympy.isprime as requested.
# ============================================================================


# ============================================================================
# BASIC NUMBER THEORY
# ============================================================================

def random_prime():
    while True:
        x = random.randrange(PRIME_LOW, PRIME_HIGH)

        if x % 2 == 0:
            x += 1

        if x >= PRIME_HIGH:
            continue

        if isprime(x):
            return x


def make_prime_population():
    print("Generating complete prime population...")

    primes = []

    for x in range(PRIME_LOW, PRIME_HIGH):
        if isprime(x):
            primes.append(x)

    return primes


def make_targets(count):
    targets = []

    while len(targets) < count:
        p = random_prime()
        q = random_prime()

        if p == q:
            continue

        n = p * q
        targets.append((p, q, n))

    return targets


# ============================================================================
# MODULUS / GENERALIZED CRT
# ============================================================================

MODULI = [r * r + 3 for r in R_VALUES]


def lcm(a, b):
    return a // math.gcd(a, b) * b


def prefix_lcm(prefix):
    M = 1

    for m in MODULI[:prefix]:
        M = lcm(M, m)

    return M


PREFIX_INFO = {
    prefix: {
        "moduli": MODULI[:prefix],
        "M": prefix_lcm(prefix),
    }
    for prefix in PREFIXES
}


# ============================================================================
# PRIME RESIDUE INDEX
# ============================================================================

def build_prefix_index(primes, prefix):
    """
    Build:

        residue -> tuple(primes)

    modulo the generalized CRT modulus M.

    Also store an inverse for every prime residue.

    Because all primes are coprime to every selected modulus here,
    every prime is a unit modulo M.
    """

    M = PREFIX_INFO[prefix]["M"]

    groups = {}

    for p in primes:
        a = p % M

        bucket = groups.get(a)

        if bucket is None:
            groups[a] = [p]
        else:
            bucket.append(p)

    # Freeze buckets for cheaper repeated lookups.
    for a in list(groups):
        groups[a] = tuple(groups[a])

    inverse = {}

    for a in groups:
        inverse[a] = pow(a, -1, M)

    return groups, inverse


# ============================================================================
# BASELINE SEARCH
# ============================================================================

def baseline_scan(primes, n):
    """
    Ordinary scan, ordered by closeness to sqrt(n).

    Returns:
        (p, q, tests, elapsed)
    """

    root = math.isqrt(n)

    # Find insertion point near sqrt(n).
    lo = 0
    hi = len(primes)

    while lo < hi:
        mid = (lo + hi) // 2

        if primes[mid] < root:
            lo = mid + 1
        else:
            hi = mid

    center = lo

    order = []

    left = center - 1
    right = center

    while left >= 0 or right < len(primes):

        if left >= 0:
            order.append(primes[left])
            left -= 1

        if right < len(primes):
            order.append(primes[right])
            right += 1

    start = time.perf_counter()

    tests = 0

    for p in order:
        tests += 1

        if n % p == 0:
            q = n // p

            if isprime(q):
                elapsed = time.perf_counter() - start
                return p, q, tests, elapsed

    elapsed = time.perf_counter() - start

    return None, None, tests, elapsed


# ============================================================================
# CRT FILTERED SEARCH
# ============================================================================

def crt_filtered_search(primes, n, groups, inverse, M):
    """
    Search using:

        q == n * p^(-1) mod M

    for each occupied prime residue class.

    No pair list is built.

    A prime p is retained only when the required q-residue class
    exists in the prime index.

    Exact divisibility is then tested on the surviving p's.
    """

    start = time.perf_counter()

    n_mod = n % M

    compatible_primes = []
    residue_hits = 0

    # Iterate residue classes rather than all prime pairs.
    for a, bucket in groups.items():

        required = (n_mod * inverse[a]) % M

        q_bucket = groups.get(required)

        if q_bucket is None:
            continue

        residue_hits += 1
        compatible_primes.extend(bucket)

    # Search compatible primes by closeness to sqrt(n).
    root = math.isqrt(n)

    compatible_primes.sort(
        key=lambda p: abs(p - root)
    )

    filter_time = time.perf_counter() - start

    divisibility_tests = 0

    for p in compatible_primes:
        divisibility_tests += 1

        if n % p == 0:
            q = n // p

            if isprime(q):
                elapsed = time.perf_counter() - start

                return {
                    "p": p,
                    "q": q,
                    "candidate_primes": len(compatible_primes),
                    "residue_hits": residue_hits,
                    "tests": divisibility_tests,
                    "filter_time": filter_time,
                    "elapsed": elapsed,
                }

    elapsed = time.perf_counter() - start

    return {
        "p": None,
        "q": None,
        "candidate_primes": len(compatible_primes),
        "residue_hits": residue_hits,
        "tests": divisibility_tests,
        "filter_time": filter_time,
        "elapsed": elapsed,
    }


# ============================================================================
# FORMATTING
# ============================================================================

def fmt_int(x):
    return f"{x:,}"


def fmt_float(x):
    return f"{x:.6f}"


# ============================================================================
# MAIN
# ============================================================================

def main():

    experiment_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 42")
    print("AMORTIZED CRT FILTER VS ORDINARY PRIME SCAN")
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
    # MODULUS INVENTORY
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. MODULUS INVENTORY")
    print("=" * 78)

    for r, m in zip(R_VALUES, MODULI):
        print(
            f"r={r:3d} "
            f"m={m:9,d}"
        )

    print()

    # ------------------------------------------------------------------------
    # PRIME POPULATION
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("2. PRIME POPULATION")
    print("=" * 78)

    start = time.perf_counter()

    primes = make_prime_population()

    prime_generation_time = time.perf_counter() - start

    print(
        f"prime population = {fmt_int(len(primes))}"
    )
    print(
        f"generation time  = "
        f"{fmt_float(prime_generation_time)}s"
    )
    print()

    # ------------------------------------------------------------------------
    # TARGETS
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("3. TARGETS")
    print("=" * 78)

    target_start = time.perf_counter()

    targets = make_targets(TARGETS)

    target_time = time.perf_counter() - target_start

    for i, (p, q, n) in enumerate(targets, 1):
        print(
            f"target {i:3d}: "
            f"p={p} "
            f"q={q} "
            f"n={n}"
        )

    print(
        f"target generation time = "
        f"{fmt_float(target_time)}s"
    )
    print()

    # ------------------------------------------------------------------------
    # BASELINE BENCHMARK
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("4. ORDINARY PRIME SCAN BENCHMARK")
    print("=" * 78)

    baseline_results = []

    baseline_start = time.perf_counter()

    for i, (p_true, q_true, n) in enumerate(targets, 1):

        p, q, tests, elapsed = baseline_scan(
            primes,
            n
        )

        correct = (
            {p, q} == {p_true, q_true}
        )

        baseline_results.append({
            "tests": tests,
            "elapsed": elapsed,
            "correct": correct,
        })

        print(
            f"target {i:3d}: "
            f"tests={fmt_int(tests):>9s} "
            f"time={elapsed:.6f}s "
            f"correct={correct}"
        )

    baseline_total = time.perf_counter() - baseline_start

    print()
    print(
        f"baseline batch time = "
        f"{baseline_total:.6f}s"
    )
    print()

    # ------------------------------------------------------------------------
    # BUILD REUSABLE CRT INDEXES
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("5. ONE-TIME CRT INDEX CONSTRUCTION")
    print("=" * 78)

    indexes = {}

    total_index_time = 0.0

    for prefix in PREFIXES:

        M = PREFIX_INFO[prefix]["M"]

        start = time.perf_counter()

        groups, inverse = build_prefix_index(
            primes,
            prefix
        )

        elapsed = time.perf_counter() - start

        indexes[prefix] = {
            "groups": groups,
            "inverse": inverse,
        }

        total_index_time += elapsed

        singleton = sum(
            1
            for bucket in groups.values()
            if len(bucket) == 1
        )

        max_class = max(
            len(bucket)
            for bucket in groups.values()
        )

        print(
            f"prefix={prefix:2d} "
            f"M={M:12,d} "
            f"classes={len(groups):9,d} "
            f"singleton={singleton:9,d} "
            f"max_class={max_class:3d} "
            f"time={elapsed:.4f}s"
        )

    print()
    print(
        f"total one-time index build = "
        f"{total_index_time:.6f}s"
    )
    print()

    # ------------------------------------------------------------------------
    # AMORTIZED BENCHMARK
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("6. AMORTIZED CRT BENCHMARK")
    print("=" * 78)

    all_prefix_results = {}

    for prefix in PREFIXES:

        info = PREFIX_INFO[prefix]
        M = info["M"]

        groups = indexes[prefix]["groups"]
        inverse = indexes[prefix]["inverse"]

        print()
        print("-" * 78)
        print(
            f"PREFIX {prefix} "
            f"(M={fmt_int(M)})"
        )
        print("-" * 78)

        results = []

        batch_start = time.perf_counter()

        for i, (p_true, q_true, n) in enumerate(
            targets,
            1
        ):

            result = crt_filtered_search(
                primes,
                n,
                groups,
                inverse,
                M
            )

            correct = (
                {result["p"], result["q"]}
                == {p_true, q_true}
            )

            result["correct"] = correct
            results.append(result)

            print(
                f"target {i:3d}: "
                f"candidates="
                f"{fmt_int(result['candidate_primes']):>9s} "
                f"residue-hits="
                f"{fmt_int(result['residue_hits']):>7s} "
                f"tests="
                f"{fmt_int(result['tests']):>8s} "
                f"time="
                f"{result['elapsed']:.6f}s "
                f"correct={correct}"
            )

        batch_time = time.perf_counter() - batch_start

        all_prefix_results[prefix] = {
            "results": results,
            "batch_time": batch_time,
        }

        print()
        print(
            f"batch time excluding index build = "
            f"{batch_time:.6f}s"
        )

    # ------------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. AMORTIZED SUMMARY")
    print("=" * 78)
    print()

    print(
        f"{'prefix':>6s} "
        f"{'M':>12s} "
        f"{'median tests':>15s} "
        f"{'median time':>15s} "
        f"{'batch time':>15s} "
        f"{'vs baseline':>15s} "
        f"{'recovery':>10s}"
    )

    print("-" * 78)

    baseline_tests = [
        x["tests"]
        for x in baseline_results
    ]

    baseline_per_target = median(
        x["elapsed"]
        for x in baseline_results
    )

    for prefix in PREFIXES:

        results = all_prefix_results[prefix]["results"]

        median_tests = median(
            x["tests"]
            for x in results
        )

        median_time = median(
            x["elapsed"]
            for x in results
        )

        batch_time = (
            total_index_time
            if False
            else all_prefix_results[prefix]["batch_time"]
        )

        total_with_index = (
            total_index_time
            + batch_time
        )

        reduction = (
            median(baseline_tests)
            / median_tests
            if median_tests
            else float("inf")
        )

        recovered = sum(
            1
            for x in results
            if x["correct"]
        )

        print(
            f"{prefix:6d} "
            f"{M:12,d} "
            f"{median_tests:15,.1f} "
            f"{median_time:15.6f} "
            f"{total_with_index:15.6f} "
            f"{reduction:14.2f}x "
            f"{recovered:5d}/{TARGETS}"
        )

    # ------------------------------------------------------------------------
    # PURE SEARCH COST, EXCLUDING ONE-TIME INDEX BUILD
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. SEARCH-ONLY COMPARISON")
    print("=" * 78)
    print()

    baseline_batch = baseline_total

    print(
        f"ordinary scan batch                  = "
        f"{baseline_batch:.6f}s"
    )

    for prefix in PREFIXES:

        crt_batch = all_prefix_results[prefix]["batch_time"]

        ratio = (
            baseline_batch / crt_batch
            if crt_batch > 0
            else float("inf")
        )

        print(
            f"prefix {prefix}: "
            f"CRT search batch                  = "
            f"{crt_batch:.6f}s "
            f"speed ratio={ratio:.3f}x"
        )

    # ------------------------------------------------------------------------
    # AMORTIZATION THRESHOLD
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. INDEX AMORTIZATION")
    print("=" * 78)
    print()

    print(
        "For each prefix, estimate the number of target moduli"
    )
    print(
        "needed before the one-time CRT index cost is paid back."
    )
    print()

    for prefix in PREFIXES:

        crt_batch = all_prefix_results[prefix]["batch_time"]

        crt_per_target = crt_batch / TARGETS

        if crt_per_target < baseline_per_target:
            savings = (
                baseline_per_target
                - crt_per_target
            )

            break_even = (
                total_index_time / savings
                if savings > 0
                else float("inf")
            )

            status = (
                f"break-even ≈ {break_even:.1f} targets"
            )

        else:
            status = "no speed break-even at current implementation"

        print(
            f"prefix {prefix}: "
            f"CRT/target={crt_per_target:.6f}s "
            f"baseline/target={baseline_per_target:.6f}s "
            f"{status}"
        )

    # ------------------------------------------------------------------------
    # CORRECTNESS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. CORRECTNESS")
    print("=" * 78)
    print()

    for prefix in PREFIXES:

        recovered = sum(
            1
            for x in all_prefix_results[prefix]["results"]
            if x["correct"]
        )

        print(
            f"prefix {prefix}: "
            f"{recovered}/{TARGETS} correctly recovered"
        )

    # ------------------------------------------------------------------------
    # INTERPRETATION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("11. INTERPRETATION")
    print("=" * 78)
    print()

    print(
        "The prime residue index is constructed exactly once."
    )

    print(
        "All subsequent targets reuse that index."
    )

    print()

    print(
        "Therefore this experiment separates:"
    )

    print(
        "  ONE-TIME COST:"
    )
    print(
        "      building the residue / inverse index"
    )

    print(
        "  REPEATED COST:"
    )
    print(
        "      filtering and testing each new n"
    )

    print()

    print(
        "The important comparison is not only the reduction"
    )
    print(
        "in divisibility tests, but whether that reduction"
    )
    print(
        "survives after the index-construction cost is amortized."
    )

    print()

    print(
        "No candidate-pair list is constructed."
    )

    print(
        "No CSV files are produced."
    )

    total_runtime = time.perf_counter() - experiment_start

    print()
    print("=" * 78)
    print("EXPERIMENT 42 COMPLETE")
    print("=" * 78)

    print(
        f"total runtime = "
        f"{total_runtime:.3f}s "
        f"({total_runtime / 60:.2f} min)"
    )


if __name__ == "__main__":
    main()

