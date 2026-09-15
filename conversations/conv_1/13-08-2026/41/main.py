#!/usr/bin/env python3

import math
import random
import time
from collections import defaultdict
from statistics import median

from sympy import isprime


# ============================================================================
# KAPPA EXPERIMENT 41
# CRT-FILTERED SEARCH VS ORDINARY PRIME SCAN
#
# FIXED:
#   The modulus family is NOT pairwise coprime.
#   Therefore the combined modulus is LCM(moduli), not product(moduli).
#
# NO CSV OUTPUT
# ============================================================================


SEED = 20260814
random.seed(SEED)

TARGETS = 12

PRIME_LOW = 2_000_000
PRIME_HIGH = 4_200_000

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47
]

PREFIXES = [3, 4, 5, 6, 7]

SEARCH_ORDER = "closest"


# ============================================================================
# FIXED TARGETS
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
# BASIC UTILITIES
# ============================================================================

def lcm(a, b):
    return abs(a * b) // math.gcd(a, b)


def lcm_list(values):
    result = 1

    for x in values:
        result = lcm(result, x)

    return result


def generate_prime_population():
    """
    Generate all primes in [PRIME_LOW, PRIME_HIGH).

    Uses sympy.isprime.
    """
    return [
        n
        for n in range(PRIME_LOW, PRIME_HIGH)
        if isprime(n)
    ]


# ============================================================================
# MODULI
# ============================================================================

MODULI = [r * r + 3 for r in R_VALUES]


def prefix_moduli(prefix):
    return MODULI[:prefix]


def combined_modulus(prefix):
    """
    IMPORTANT:
    The moduli overlap heavily, so their product is NOT the CRT modulus.

    We use the least common multiple.
    """
    return lcm_list(prefix_moduli(prefix))


# ============================================================================
# PRIME INDEX
# ============================================================================

def build_prime_index(primes, M):
    """
    Index primes by p mod M.

    Because all primes in the search interval are larger than the
    small moduli, this is sufficient for the N-only CRT filter.
    """

    index = defaultdict(list)

    for p in primes:
        index[p % M].append(p)

    return index


# ============================================================================
# N-ONLY COMPATIBILITY
# ============================================================================

def compatible_q_class(n, p_class, M):
    """
    Given

        p*q = n (mod M)

    and gcd(p, M)=1,

    q = n * p^{-1} (mod M).
    """

    if math.gcd(p_class, M) != 1:
        return None

    return (n * pow(p_class, -1, M)) % M


def build_compatible_prime_classes(prime_index, n, M):
    """
    Determine which p-residue classes can participate in

        p*q == n (mod M)

    while requiring both p and q to occur in the prime population.

    Returns:
        compatible p classes
        compatible q classes
    """

    compatible_p = []
    compatible_q = []

    for p_class in prime_index.keys():

        if math.gcd(p_class, M) != 1:
            continue

        q_class = compatible_q_class(
            n,
            p_class,
            M
        )

        if q_class is None:
            continue

        if q_class not in prime_index:
            continue

        compatible_p.append(p_class)
        compatible_q.append(q_class)

    return compatible_p, compatible_q


# ============================================================================
# SEARCH ORDER
# ============================================================================

def order_candidates(candidates, n, mode):
    if mode == "ascending":
        return sorted(candidates)

    if mode == "descending":
        return sorted(candidates, reverse=True)

    if mode == "closest":
        sqrt_n = math.isqrt(n)

        return sorted(
            candidates,
            key=lambda p: abs(p - sqrt_n)
        )

    raise ValueError(
        f"Unknown SEARCH_ORDER={mode!r}"
    )


# ============================================================================
# BASELINE SEARCH
# ============================================================================

def baseline_search(primes, n):
    """
    Ordinary prime scan.

    No CRT filtering.
    """

    ordered = order_candidates(
        primes,
        n,
        SEARCH_ORDER
    )

    start = time.perf_counter()

    divisibility_tests = 0

    for p in ordered:

        divisibility_tests += 1

        if n % p != 0:
            continue

        q = n // p

        if q == p:
            continue

        if PRIME_LOW <= q < PRIME_HIGH and isprime(q):

            elapsed = time.perf_counter() - start

            return {
                "p": min(p, q),
                "q": max(p, q),
                "tests": divisibility_tests,
                "elapsed": elapsed,
            }

    elapsed = time.perf_counter() - start

    return {
        "p": None,
        "q": None,
        "tests": divisibility_tests,
        "elapsed": elapsed,
    }


# ============================================================================
# CRT-FILTERED SEARCH
# ============================================================================

def crt_filtered_search(primes, n, prefix):
    """
    N-only CRT filtering.

    No candidate-pair list is constructed.

    Each prime p is examined only through its residue class modulo M.
    """

    total_start = time.perf_counter()

    M = combined_modulus(prefix)

    # ------------------------------------------------------------
    # Build prime residue index
    # ------------------------------------------------------------

    prep_start = time.perf_counter()

    prime_index = build_prime_index(
        primes,
        M
    )

    compatible_p_classes, compatible_q_classes = (
        build_compatible_prime_classes(
            prime_index,
            n,
            M
        )
    )

    prep_time = time.perf_counter() - prep_start

    compatible_p_class_set = set(
        compatible_p_classes
    )

    # ------------------------------------------------------------
    # Streaming candidate generation
    #
    # We intentionally DO NOT construct a candidate-pair list.
    # ------------------------------------------------------------

    search_start = time.perf_counter()

    candidate_prime_count = sum(
        len(prime_index[c])
        for c in compatible_p_class_set
    )

    # Only construct the p candidate list.
    candidate_primes = []

    for c in compatible_p_class_set:
        candidate_primes.extend(
            prime_index[c]
        )

    ordered = order_candidates(
        candidate_primes,
        n,
        SEARCH_ORDER
    )

    divisibility_tests = 0

    for p in ordered:

        divisibility_tests += 1

        if n % p != 0:
            continue

        q = n // p

        if q == p:
            continue

        if not (
            PRIME_LOW <= q < PRIME_HIGH
        ):
            continue

        if not isprime(q):
            continue

        # --------------------------------------------------------
        # Exact verification
        # --------------------------------------------------------

        if p * q != n:
            continue

        if (
            p % M in compatible_p_class_set
            and q in prime_index.get(q % M, [])
        ):
            search_time = (
                time.perf_counter()
                - search_start
            )

            total_time = (
                time.perf_counter()
                - total_start
            )

            return {
                "M": M,
                "candidate_primes": candidate_prime_count,
                "compatible_p_classes": len(
                    compatible_p_class_set
                ),
                "divisibility_tests": divisibility_tests,
                "p": min(p, q),
                "q": max(p, q),
                "prep_time": prep_time,
                "search_time": search_time,
                "total_time": total_time,
            }

    search_time = (
        time.perf_counter()
        - search_start
    )

    total_time = (
        time.perf_counter()
        - total_start
    )

    return {
        "M": M,
        "candidate_primes": candidate_prime_count,
        "compatible_p_classes": len(
            compatible_p_class_set
        ),
        "divisibility_tests": divisibility_tests,
        "p": None,
        "q": None,
        "prep_time": prep_time,
        "search_time": search_time,
        "total_time": total_time,
    }


# ============================================================================
# SUMMARY HELPERS
# ============================================================================

def median_int(values):
    if not values:
        return 0

    return median(values)


def speedup(old, new):
    if new == 0:
        return float("inf")

    return old / new


# ============================================================================
# MAIN
# ============================================================================

def main():

    experiment_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 41")
    print("CRT-FILTERED SEARCH VS ORDINARY PRIME SCAN")
    print("FIXED GENERALIZED-CRT MODULUS")
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
    print(f"search order     = {SEARCH_ORDER}")
    print()

    # ========================================================================
    # 1. MODULUS INVENTORY
    # ========================================================================

    print("=" * 78)
    print("1. MODULUS INVENTORY")
    print("=" * 78)

    for r, m in zip(R_VALUES, MODULI):
        print(
            f"r={r:3d} "
            f"m={m:10,d}"
        )

    print()

    print("GENERALIZED CRT PREFIX MODULI")
    print("-" * 78)

    for prefix in PREFIXES:

        mods = prefix_moduli(prefix)
        M = combined_modulus(prefix)

        print(
            f"prefix={prefix:2d} "
            f"moduli={mods} "
            f"LCM={M:,}"
        )

    print()

    # ========================================================================
    # 2. PRIME POPULATION
    # ========================================================================

    print("=" * 78)
    print("2. PRIME POPULATION")
    print("=" * 78)

    start = time.perf_counter()

    primes = generate_prime_population()

    generation_time = (
        time.perf_counter() - start
    )

    print(
        f"prime population = "
        f"{len(primes):,}"
    )

    print(
        f"generation time  = "
        f"{generation_time:.3f}s"
    )

    print()

    # ========================================================================
    # 3. TARGETS
    # ========================================================================

    print("=" * 78)
    print("3. TARGETS")
    print("=" * 78)

    targets = []

    for target_id, (p, q) in enumerate(
        FIXED_TARGETS,
        1
    ):

        n = p * q

        targets.append(
            (p, q, n)
        )

        print(
            f"target {target_id:2d}: "
            f"p={p} "
            f"q={q} "
            f"n={n}"
        )

    print()

    # ========================================================================
    # 4. BASELINE
    # ========================================================================

    print("=" * 78)
    print("4. BASELINE PRIME SCAN")
    print("=" * 78)

    print(
        "Ordinary prime scan:"
    )

    print(
        "  test n % p for primes ordered by closeness to sqrt(n)"
    )

    print()

    baseline_results = []

    for target_id, (p, q, n) in enumerate(
        targets,
        1
    ):

        result = baseline_search(
            primes,
            n
        )

        correct = (
            result["p"] == min(p, q)
            and result["q"] == max(p, q)
        )

        baseline_results.append(result)

        print(
            f"target {target_id:2d}: "
            f"tests={result['tests']:8,d} "
            f"time={result['elapsed']:.6f}s "
            f"correct={correct}"
        )

    print()

    # ========================================================================
    # 5. CRT FILTERED SEARCH
    # ========================================================================

    print("=" * 78)
    print("5. GENERALIZED-CRT FILTERED SEARCH")
    print("=" * 78)

    all_results = {
        prefix: []
        for prefix in PREFIXES
    }

    for prefix in PREFIXES:

        M = combined_modulus(prefix)

        print()
        print("-" * 78)
        print(
            f"PREFIX {prefix} "
            f"LCM(M)={M:,}"
        )
        print("-" * 78)

        for target_id, (p, q, n) in enumerate(
            targets,
            1
        ):

            result = crt_filtered_search(
                primes,
                n,
                prefix
            )

            all_results[prefix].append(
                result
            )

            correct = (
                result["p"] == min(p, q)
                and result["q"] == max(p, q)
            )

            baseline = baseline_results[
                target_id - 1
            ]

            reduction = speedup(
                baseline["tests"],
                result["divisibility_tests"]
            )

            print(
                f"target {target_id:2d}: "
                f"classes={result['compatible_p_classes']:7,d} "
                f"candidate_primes="
                f"{result['candidate_primes']:8,d} "
                f"CRT_tests="
                f"{result['divisibility_tests']:8,d} "
                f"reduction="
                f"{reduction:9.3f}x "
                f"correct={correct}"
            )

            print(
                f"             "
                f"prep={result['prep_time']:.4f}s "
                f"search={result['search_time']:.4f}s "
                f"total={result['total_time']:.4f}s"
            )

    # ========================================================================
    # 6. GLOBAL SUMMARY
    # ========================================================================

    print()
    print("=" * 78)
    print("6. GLOBAL SUMMARY")
    print("=" * 78)
    print()

    print(
        f"{'prefix':>6s} "
        f"{'M':>14s} "
        f"{'median candidates':>20s} "
        f"{'median tests':>16s} "
        f"{'median baseline':>18s} "
        f"{'median reduction':>18s}"
    )

    print("-" * 78)

    for prefix in PREFIXES:

        results = all_results[prefix]

        candidates = [
            x["candidate_primes"]
            for x in results
        ]

        tests = [
            x["divisibility_tests"]
            for x in results
        ]

        base_tests = [
            x["tests"]
            for x in baseline_results
        ]

        reductions = [
            speedup(
                base_tests[i],
                tests[i]
            )
            for i in range(TARGETS)
        ]

        print(
            f"{prefix:6d} "
            f"{combined_modulus(prefix):14,d} "
            f"{median_int(candidates):20,.1f} "
            f"{median_int(tests):16,.1f} "
            f"{median_int(base_tests):18,.1f} "
            f"{median(reductions):18,.3f}x"
        )

    # ========================================================================
    # 7. RECOVERY
    # ========================================================================

    print()
    print("=" * 78)
    print("7. RECOVERY CHECK")
    print("=" * 78)

    for prefix in PREFIXES:

        recovered = 0

        for result, target in zip(
            all_results[prefix],
            targets
        ):

            p, q, _ = target

            if (
                result["p"] == min(p, q)
                and result["q"] == max(p, q)
            ):
                recovered += 1

        print(
            f"prefix {prefix:2d}: "
            f"{recovered}/{TARGETS} "
            f"correctly recovered"
        )

    # ========================================================================
    # 8. IMPORTANT STRUCTURAL CHECK
    # ========================================================================

    print()
    print("=" * 78)
    print("8. STRUCTURAL CHECK")
    print("=" * 78)

    print("""
The earlier implementation failed because the modulus family is
not pairwise coprime.

For example:

    lcm(7, 12, 28) = 84

not

    7 * 12 * 28 = 2352.

This experiment therefore works directly modulo the generalized
CRT modulus:

    M = lcm(selected moduli).

The N-only condition is:

    p*q == n (mod M)

and, for a unit p,

    q == n * p^(-1) (mod M).

This is mathematically equivalent to imposing all selected congruences
simultaneously, while correctly handling their shared factors.
""")

    # ========================================================================
    # 9. INTERPRETATION
    # ========================================================================

    print()
    print("=" * 78)
    print("9. INTERPRETATION")
    print("=" * 78)

    print("""
This experiment measures actual search work rather than pair counts.

BASELINE:
    scan primes and test n % p.

CRT FILTER:
    first use n modulo the generalized CRT modulus to reject primes
    whose residue class cannot have a compatible prime q.

Then test n % p only on the surviving p candidates.

The important quantities are:

    candidate-prime reduction
    divisibility-test reduction
    preprocessing time
    total runtime
    exact recovery correctness

A useful result requires more than a smaller residue-compatible
pair population.

The CRT method should reduce the number of actual divisibility
tests enough to compensate for its residue-index construction cost.

The prefix at which this becomes beneficial is the practical
threshold to investigate next.
""")

    # ========================================================================
    # FINAL
    # ========================================================================

    total_time = (
        time.perf_counter()
        - experiment_start
    )

    print()
    print("=" * 78)
    print("EXPERIMENT 41 COMPLETE")
    print("=" * 78)

    print(
        f"total runtime = "
        f"{total_time:.3f}s "
        f"({total_time / 60:.2f} min)"
    )


if __name__ == "__main__":
    main()