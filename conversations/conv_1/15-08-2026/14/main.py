#!/usr/bin/env python3

import math
import random
import time
from collections import defaultdict
from statistics import median, fmean

from sympy import isprime


# ============================================================================
# KAPPA EXPERIMENT 42
#
# AMORTIZED CRT FILTER BENCHMARK
#
# Question:
#   Does the dramatic divisibility-test reduction from Experiment 41
#   survive when the CRT / residue infrastructure is precomputed once
#   and reused across many independent semiprimes?
#
# Comparisons:
#   1. ORDINARY PRIME SCAN
#   2. CRT PREFIX-6
#   3. CRT PREFIX-7
#
# For each CRT family we separately report:
#
#   A. preprocessing time
#   B. N-dependent compatibility-pass time
#   C. exact divisibility-test time
#   D. total amortized time/target
#
# This distinction is critical:
#
#   Experiment 41:
#       "177x fewer exact divisibility tests"
#
#   Experiment 42 asks:
#       "Does the COMPLETE ALGORITHM become faster after including
#        the cost of finding those compatible residue classes?"
#
# NO CSV OUTPUT
# NO SKLEARN
# ============================================================================


SEED = 20260815

PRIME_LOW = 2_000_000
PRIME_HIGH = 4_200_000

# Enough targets to make one-time preprocessing amortization visible.
NUM_TARGETS = 200

PREFIXES = [5, 6, 7]

SEARCH_ORDER = "closest"

# Fixed legacy targets from Experiment 41.
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

# Same family as Experiment 41.
R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47
]

MODULI = [r * r + 3 for r in R_VALUES]


# ============================================================================
# BASIC UTILITIES
# ============================================================================

def lcm(a: int, b: int) -> int:
    return abs(a * b) // math.gcd(a, b)


def lcm_list(values) -> int:
    out = 1
    for x in values:
        out = lcm(out, x)
    return out


def combined_modulus(prefix: int) -> int:
    return lcm_list(MODULI[:prefix])


def generate_prime_population():
    return [
        p
        for p in range(PRIME_LOW, PRIME_HIGH)
        if isprime(p)
    ]


def order_candidates(candidates, n: int):
    if SEARCH_ORDER == "ascending":
        return sorted(candidates)

    if SEARCH_ORDER == "descending":
        return sorted(candidates, reverse=True)

    if SEARCH_ORDER == "closest":
        sqrt_n = math.isqrt(n)
        return sorted(
            candidates,
            key=lambda p: abs(p - sqrt_n),
        )

    raise ValueError(f"Unknown SEARCH_ORDER={SEARCH_ORDER!r}")


# ============================================================================
# TARGET GENERATION
# ============================================================================

def build_targets(primes, count: int):
    """
    Build deterministic semiprimes.

    We retain the 12 Experiment-41 targets and then add reproducible
    random semiprimes from the same prime population.
    """

    rng = random.Random(SEED)

    targets = []

    seen_n = set()

    # ------------------------------------------------------------
    # Legacy targets first
    # ------------------------------------------------------------

    for p, q in FIXED_TARGETS:
        n = p * q

        if n not in seen_n:
            targets.append((p, q, n))
            seen_n.add(n)

    # ------------------------------------------------------------
    # Additional targets
    # ------------------------------------------------------------

    prime_count = len(primes)

    while len(targets) < count:
        i = rng.randrange(prime_count)
        j = rng.randrange(prime_count)

        if i == j:
            continue

        p = primes[i]
        q = primes[j]

        n = p * q

        if n in seen_n:
            continue

        targets.append((p, q, n))
        seen_n.add(n)

    return targets


# ============================================================================
# ORDINARY BASELINE
# ============================================================================

def build_baseline_order(primes, targets):
    """
    Sort the entire prime population once for each target.

    This matches the operational structure of Experiment 41.
    """

    ordered = []

    for _, _, n in targets:
        ordered.append(order_candidates(primes, n))

    return ordered


def baseline_search(ordered_primes, n: int):
    """
    Exact ordinary scan.

    Returns:
        divisibility_tests
        elapsed
        recovered
    """

    start = time.perf_counter()

    tests = 0

    for p in ordered_primes:
        tests += 1

        if n % p != 0:
            continue

        q = n // p

        if q == p:
            continue

        if not (PRIME_LOW <= q < PRIME_HIGH):
            continue

        if not isprime(q):
            continue

        return {
            "tests": tests,
            "elapsed": time.perf_counter() - start,
            "p": min(p, q),
            "q": max(p, q),
        }

    return {
        "tests": tests,
        "elapsed": time.perf_counter() - start,
        "p": None,
        "q": None,
    }


# ============================================================================
# CRT INFRASTRUCTURE
# ============================================================================

def build_residue_index(primes, M: int):
    """
    Map residue -> tuple/list of primes having that residue.
    """

    index = defaultdict(list)

    for p in primes:
        index[p % M].append(p)

    return index


def build_inverse_cache(primes, M: int):
    """
    Cache p^{-1} mod M for every prime p coprime to M.

    This separates the one-time inverse cost from the per-target
    compatibility pass.
    """

    inverses = {}

    for p in primes:
        r = p % M

        if math.gcd(r, M) != 1:
            continue

        inverses[p] = pow(r, -1, M)

    return inverses


def prepare_crt_family(primes, prefix: int):
    """
    One-time preprocessing for a prefix.
    """

    M = combined_modulus(prefix)

    start = time.perf_counter()

    residue_index = build_residue_index(
        primes,
        M,
    )

    inverse_cache = build_inverse_cache(
        primes,
        M,
    )

    ordered_primes = sorted(primes)

    prep_time = time.perf_counter() - start

    # Some diagnostics about residue uniqueness.
    class_sizes = [
        len(v)
        for v in residue_index.values()
    ]

    singleton_classes = sum(
        1 for x in class_sizes if x == 1
    )

    max_class_size = (
        max(class_sizes)
        if class_sizes
        else 0
    )

    return {
        "prefix": prefix,
        "M": M,
        "residue_index": residue_index,
        "inverse_cache": inverse_cache,
        "ordered_primes": ordered_primes,
        "prep_time": prep_time,
        "num_classes": len(residue_index),
        "singleton_classes": singleton_classes,
        "max_class_size": max_class_size,
        "class_sizes": class_sizes,
    }


# ============================================================================
# CRT FILTERED SEARCH
# ============================================================================

def crt_filtered_search(family, n: int):
    """
    N-dependent compatibility pass:

        q == n * p^{-1} (mod M)

    For every possible p residue class, determine whether the implied
    q residue occurs among actual primes.

    Then perform exact divisibility testing only on surviving p values.

    This deliberately charges the compatibility-pass cost separately.
    """

    M = family["M"]
    residue_index = family["residue_index"]
    inverse_cache = family["inverse_cache"]

    # ------------------------------------------------------------
    # A. N-dependent compatibility pass
    # ------------------------------------------------------------

    compatibility_start = time.perf_counter()

    n_mod = n % M

    compatible_primes = []

    residue_classes_checked = 0
    compatible_classes = set()

    for p_class, p_values in residue_index.items():

        residue_classes_checked += 1

        if math.gcd(p_class, M) != 1:
            continue

        inv = pow(p_class, -1, M)
        q_class = (n_mod * inv) % M

        if q_class not in residue_index:
            continue

        compatible_classes.add(p_class)
        compatible_primes.extend(p_values)

    compatibility_time = (
        time.perf_counter() - compatibility_start
    )

    # ------------------------------------------------------------
    # B. Exact divisibility checks
    # ------------------------------------------------------------

    search_start = time.perf_counter()

    ordered = order_candidates(
        compatible_primes,
        n,
    )

    divisibility_tests = 0

    for p in ordered:

        divisibility_tests += 1

        if n % p != 0:
            continue

        q = n // p

        if q == p:
            continue

        if not (PRIME_LOW <= q < PRIME_HIGH):
            continue

        if not isprime(q):
            continue

        return {
            "M": M,
            "candidate_primes": len(compatible_primes),
            "compatible_classes": len(compatible_classes),
            "classes_checked": residue_classes_checked,
            "compatibility_time": compatibility_time,
            "search_time": time.perf_counter() - search_start,
            "divisibility_tests": divisibility_tests,
            "p": min(p, q),
            "q": max(p, q),
        }

    return {
        "M": M,
        "candidate_primes": len(compatible_primes),
        "compatible_classes": len(compatible_classes),
        "classes_checked": residue_classes_checked,
        "compatibility_time": compatibility_time,
        "search_time": time.perf_counter() - search_start,
        "divisibility_tests": divisibility_tests,
        "p": None,
        "q": None,
    }


# ============================================================================
# ULTRA-OPTIMIZED VARIANT
# ============================================================================

def build_prime_to_residue_records(primes, M):
    """
    Optional optimization.

    Instead of recomputing p % M and gcd information during each target,
    cache:

        (p, p_mod_M, inverse_mod_M)

    once.

    This does NOT remove the fundamental O(#primes) compatibility pass,
    but removes repeated modular setup work.
    """

    start = time.perf_counter()

    records = []

    for p in primes:

        r = p % M

        if math.gcd(r, M) != 1:
            continue

        inv = pow(r, -1, M)

        records.append(
            (p, r, inv)
        )

    prep_time = time.perf_counter() - start

    return records, prep_time


def cached_record_search(records, residue_index, n):
    """
    Same compatibility logic using completely cached prime records.

    This is the most favorable implementation we test before considering
    more sophisticated data structures.
    """

    M = None

    if not residue_index:
        return None

    # We store M externally in the caller.
    raise RuntimeError(
        "M must be passed explicitly to cached_record_search"
    )


def cached_record_search_with_M(records, residue_index, M, n):
    start_compat = time.perf_counter()

    n_mod = n % M

    compatible_primes = []
    compatible_classes = set()

    for p, _, inv in records:

        q_class = (n_mod * inv) % M

        if q_class not in residue_index:
            continue

        compatible_classes.add(p % M)
        compatible_primes.append(p)

    compat_time = time.perf_counter() - start_compat

    start_exact = time.perf_counter()

    ordered = order_candidates(
        compatible_primes,
        n,
    )

    tests = 0

    for p in ordered:

        tests += 1

        if n % p != 0:
            continue

        q = n // p

        if q == p:
            continue

        if not (PRIME_LOW <= q < PRIME_HIGH):
            continue

        if not isprime(q):
            continue

        return {
            "compatibility_time": compat_time,
            "search_time": time.perf_counter() - start_exact,
            "candidate_primes": len(compatible_primes),
            "compatible_classes": len(compatible_classes),
            "tests": tests,
            "p": min(p, q),
            "q": max(p, q),
        }

    return {
        "compatibility_time": compat_time,
        "search_time": time.perf_counter() - start_exact,
        "candidate_primes": len(compatible_primes),
        "compatible_classes": len(compatible_classes),
        "tests": tests,
        "p": None,
        "q": None,
    }


# ============================================================================
# SUMMARY
# ============================================================================

def summarize(values):
    if not values:
        return {
            "mean": 0.0,
            "median": 0.0,
        }

    return {
        "mean": fmean(values),
        "median": median(values),
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    experiment_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 42")
    print("AMORTIZED CRT FILTER BENCHMARK")
    print("PRECOMPUTE ONCE / MANY TARGETS")
    print("EXPERIMENT 41 FOLLOW-UP")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)
    print()

    print(f"seed              = {SEED}")
    print(f"prime interval    = [{PRIME_LOW:,}, {PRIME_HIGH:,})")
    print(f"targets            = {NUM_TARGETS}")
    print(f"prefixes           = {PREFIXES}")
    print()

    # ========================================================================
    # 1. MODULUS INVENTORY
    # ========================================================================

    print("=" * 78)
    print("1. MODULUS INVENTORY")
    print("=" * 78)

    for prefix in PREFIXES:
        mods = MODULI[:prefix]
        M = combined_modulus(prefix)

        print(
            f"prefix={prefix:2d} "
            f"M={M:,} "
            f"mods={mods}"
        )

    print()

    # ========================================================================
    # 2. PRIME POPULATION
    # ========================================================================

    print("=" * 78)
    print("2. PRIME POPULATION")
    print("=" * 78)

    t0 = time.perf_counter()
    primes = generate_prime_population()
    prime_generation = time.perf_counter() - t0

    print(f"prime population = {len(primes):,}")
    print(f"generation time  = {prime_generation:.6f}s")
    print()

    # ========================================================================
    # 3. TARGETS
    # ========================================================================

    print("=" * 78)
    print("3. TARGETS")
    print("=" * 78)

    targets = build_targets(
        primes,
        NUM_TARGETS,
    )

    print(f"generated targets = {len(targets)}")

    print()
    print("first 12 targets:")
    for i, (p, q, n) in enumerate(targets[:12], 1):
        print(
            f"target {i:3d}: "
            f"p={p} q={q} n={n}"
        )

    print()

    # ========================================================================
    # 4. BASELINE
    # ========================================================================

    print("=" * 78)
    print("4. ORDINARY PRIME SCAN")
    print("=" * 78)

    baseline_results = []

    baseline_start = time.perf_counter()

    for idx, (_, _, n) in enumerate(targets, 1):

        ordered = order_candidates(
            primes,
            n,
        )

        result = baseline_search(
            ordered,
            n,
        )

        baseline_results.append(result)

        if idx % 25 == 0:
            print(f"baseline processed {idx:4d}/{NUM_TARGETS}")

    baseline_total = time.perf_counter() - baseline_start

    base_tests = [
        x["tests"]
        for x in baseline_results
    ]

    base_times = [
        x["elapsed"]
        for x in baseline_results
    ]

    base_summary = summarize(base_times)

    print()
    print(
        f"baseline total time   = {baseline_total:.6f}s"
    )
    print(
        f"baseline mean/target  = {base_summary['mean']:.6f}s"
    )
    print(
        f"baseline median/target= {base_summary['median']:.6f}s"
    )
    print(
        f"baseline median tests = {median(base_tests):,.1f}"
    )
    print()

    # ========================================================================
    # 5. CRT FAMILY PREPARATION
    # ========================================================================

    families = {}

    print("=" * 78)
    print("5. ONE-TIME CRT PREPROCESSING")
    print("=" * 78)

    for prefix in PREFIXES:

        family_start = time.perf_counter()

        family = prepare_crt_family(
            primes,
            prefix,
        )

        family_total = time.perf_counter() - family_start

        family["prep_time_total"] = family_total

        families[prefix] = family

        print(
            f"prefix={prefix:2d} "
            f"M={family['M']:,}"
        )
        print(
            f"  preparation       = {family_total:.6f}s"
        )
        print(
            f"  residue classes   = {family['num_classes']:,}"
        )
        print(
            f"  singleton classes = "
            f"{family['singleton_classes']:,}"
        )
        print(
            f"  max class size    = "
            f"{family['max_class_size']:,}"
        )
        print()

    # ========================================================================
    # 6. PER-TARGET CRT BENCHMARK
    # ========================================================================

    all_results = {}

    for prefix in PREFIXES:

        family = families[prefix]

        print("=" * 78)
        print(
            f"6. CRT PREFIX {prefix} "
            f"M={family['M']:,}"
        )
        print("=" * 78)

        results = []

        total_target_start = time.perf_counter()

        for idx, (_, _, n) in enumerate(targets, 1):

            result = crt_filtered_search(
                family,
                n,
            )

            baseline = baseline_results[idx - 1]

            result["total_per_target"] = (
                result["compatibility_time"]
                + result["search_time"]
            )

            result["with_prep_amortized"] = (
                result["total_per_target"]
                + family["prep_time_total"] / NUM_TARGETS
            )

            result["test_reduction"] = (
                baseline["tests"]
                / max(1, result["divisibility_tests"])
            )

            result["correct"] = (
                result["p"] == min(targets[idx - 1][0],
                                    targets[idx - 1][1])
                and
                result["q"] == max(targets[idx - 1][0],
                                    targets[idx - 1][1])
            )

            results.append(result)

            if idx % 25 == 0:
                print(
                    f"processed {idx:4d}/{NUM_TARGETS}"
                )

        family_total_targets = (
            time.perf_counter()
            - total_target_start
        )

        all_results[prefix] = results

        compat_times = [
            x["compatibility_time"]
            for x in results
        ]

        exact_times = [
            x["search_time"]
            for x in results
        ]

        total_times = [
            x["total_per_target"]
            for x in results
        ]

        amortized_times = [
            x["with_prep_amortized"]
            for x in results
        ]

        candidate_counts = [
            x["candidate_primes"]
            for x in results
        ]

        exact_tests = [
            x["divisibility_tests"]
            for x in results
        ]

        reductions = [
            x["test_reduction"]
            for x in results
        ]

        correct_count = sum(
            x["correct"]
            for x in results
        )

        print()
        print(
            f"targets total search wall time = "
            f"{family_total_targets:.6f}s"
        )
        print(
            f"mean compatibility pass       = "
            f"{fmean(compat_times):.6f}s"
        )
        print(
            f"mean exact-search time        = "
            f"{fmean(exact_times):.6f}s"
        )
        print(
            f"mean CRT total/target         = "
            f"{fmean(total_times):.6f}s"
        )
        print(
            f"mean amortized/target         = "
            f"{fmean(amortized_times):.6f}s"
        )
        print(
            f"median candidate primes       = "
            f"{median(candidate_counts):,.1f}"
        )
        print(
            f"median exact divisibility     = "
            f"{median(exact_tests):,.1f}"
        )
        print(
            f"median test reduction         = "
            f"{median(reductions):,.3f}x"
        )
        print(
            f"exact recovery                = "
            f"{correct_count}/{NUM_TARGETS}"
        )

    # ========================================================================
    # 7. BREAK-EVEN ANALYSIS
    # ========================================================================

    print()
    print("=" * 78)
    print("7. AMORTIZATION / BREAK-EVEN ANALYSIS")
    print("=" * 78)

    baseline_mean = fmean(base_times)

    print(
        f"baseline mean time/target = "
        f"{baseline_mean:.6f}s"
    )
    print()

    for prefix in PREFIXES:

        family = families[prefix]
        results = all_results[prefix]

        mean_crt_without_prep = fmean(
            x["total_per_target"]
            for x in results
        )

        prep = family["prep_time_total"]

        savings_per_target = (
            baseline_mean
            - mean_crt_without_prep
        )

        print(
            f"PREFIX {prefix}"
        )
        print(
            f"  preprocessing = {prep:.6f}s"
        )
        print(
            f"  mean CRT target cost excl. prep = "
            f"{mean_crt_without_prep:.6f}s"
        )

        if savings_per_target <= 0:
            print(
                "  break-even = NEVER "
                "(CRT target cost is not below baseline)"
            )
        else:
            break_even = prep / savings_per_target

            print(
                f"  break-even targets ~= "
                f"{break_even:.2f}"
            )

        print()

    # ========================================================================
    # 8. ACTUAL TOTAL COST AT DIFFERENT BATCH SIZES
    # ========================================================================

    print("=" * 78)
    print("8. BATCH-SIZE SCALING")
    print("=" * 78)

    batch_sizes = [1, 10, 50, 100, NUM_TARGETS]

    for batch in batch_sizes:

        if batch > NUM_TARGETS:
            continue

        baseline_batch = (
            sum(base_times[:batch])
        )

        print()
        print(f"BATCH = {batch}")

        print(
            f"  baseline = "
            f"{baseline_batch:.6f}s"
        )

        for prefix in PREFIXES:

            family = families[prefix]
            results = all_results[prefix]

            crt_batch = (
                family["prep_time_total"]
                +
                sum(
                    x["total_per_target"]
                    for x in results[:batch]
                )
            )

            ratio = (
                baseline_batch / crt_batch
                if crt_batch > 0
                else float("inf")
            )

            print(
                f"  prefix {prefix}: "
                f"{crt_batch:.6f}s "
                f"speed ratio={ratio:.3f}x"
            )

    # ========================================================================
    # 9. TEST-COUNT VS WALL-CLOCK
    # ========================================================================

    print()
    print("=" * 78)
    print("9. TEST-COUNT VS REAL TIME")
    print("=" * 78)

    for prefix in PREFIXES:

        results = all_results[prefix]

        med_tests = median(
            x["divisibility_tests"]
            for x in results
        )

        med_reduction = median(
            x["test_reduction"]
            for x in results
        )

        mean_time = fmean(
            x["total_per_target"]
            for x in results
        )

        print(
            f"prefix {prefix}: "
            f"median exact tests={med_tests:,.1f} "
            f"reduction={med_reduction:.3f}x "
            f"CRT target cost={mean_time:.6f}s"
        )

    # ========================================================================
    # 10. STRUCTURAL INTERPRETATION
    # ========================================================================

    print()
    print("=" * 78)
    print("10. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 42 separates three different resources:

    1. ONE-TIME PREPROCESSING
       Build the residue index and modular inverses.

    2. N-DEPENDENT COMPATIBILITY WORK
       For a new n, determine which p residues imply a q residue
       that actually occurs in the prime population.

    3. EXACT FACTOR CHECKING
       Perform n % p only on the surviving p candidates.

The key question is NOT:

    "How many exact divisibility tests disappeared?"

That was already demonstrated in Experiment 41.

The new question is:

    "After charging the compatibility computation and amortizing
     the fixed residue infrastructure across many targets, is the
     entire factoring procedure faster than an ordinary prime scan?"

A strong result would have the following pattern:

    low preprocessing cost
    + low per-target compatibility cost
    + dramatic reduction in exact tests
    + lower total wall-clock time

If the exact-test reduction remains enormous but total runtime
does not improve, the method is a useful structural sieve but
not yet a practical factoring speedup.

The prefix-7 regime is especially important because its modulus
is much larger than the entire factor-search interval. In that
regime the residue classes are expected to become nearly singleton,
so we explicitly report class sizes rather than hiding this effect.
"""
    )

    # ========================================================================
    # 11. FINAL SUMMARY
    # ========================================================================

    print()
    print("=" * 78)
    print("11. FINAL SUMMARY")
    print("=" * 78)

    for prefix in PREFIXES:

        family = families[prefix]
        results = all_results[prefix]

        mean_total = fmean(
            x["total_per_target"]
            for x in results
        )

        mean_amortized = fmean(
            x["with_prep_amortized"]
            for x in results
        )

        mean_tests = fmean(
            x["divisibility_tests"]
            for x in results
        )

        exact_speedup = (
            fmean(base_tests)
            / mean_tests
        )

        wall_speedup = (
            baseline_mean
            / mean_amortized
        )

        print(
            f"prefix {prefix}: "
            f"M={family['M']:,}"
        )
        print(
            f"  exact-test reduction = "
            f"{exact_speedup:.3f}x"
        )
        print(
            f"  mean total/target    = "
            f"{mean_total:.6f}s"
        )
        print(
            f"  amortized/target     = "
            f"{mean_amortized:.6f}s"
        )
        print(
            f"  wall-clock ratio     = "
            f"{wall_speedup:.3f}x"
        )
        print()

    # ========================================================================
    # 12. FINAL
    # ========================================================================

    total_time = time.perf_counter() - experiment_start

    print("=" * 78)
    print("EXPERIMENT 42 COMPLETE")
    print("=" * 78)
    print(
        f"total runtime = "
        f"{total_time:.6f}s "
        f"({total_time / 60:.2f} min)"
    )


if __name__ == "__main__":
    main()

