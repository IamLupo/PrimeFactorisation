#!/usr/bin/env python3

import math
import random
import time

from collections import Counter, defaultdict
from sympy import isprime


# ==============================================================================
# KAPPA EXPERIMENT 40
# STREAMING N-ONLY CRT RECOVERY
# NO CSV OUTPUT
#
# MEMORY-SAFE VERSION OF EXPERIMENT 39
# ==============================================================================

SEED = 20260814
random.seed(SEED)

TARGETS = 12

PRIME_LOW = 2_000_000
PRIME_HIGH = 4_200_000

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47
]

PREFIXES = [3, 4, 5, 6, 7]

SEARCH_ORDER = "smallest-prime"

# Only print examples when a candidate class is tiny.
PRINT_EXAMPLES_AT_MOST = 12


# ==============================================================================
# MODULI
# ==============================================================================

MODULI = [r * r + 3 for r in R_VALUES]


def lcm_many(values):
    result = 1

    for x in values:
        result = math.lcm(result, x)

    return result


PREFIX_MODULI = {
    k: MODULI[:k]
    for k in PREFIXES
}

PREFIX_M = {
    k: lcm_many(PREFIX_MODULI[k])
    for k in PREFIXES
}


# ==============================================================================
# PRIME GENERATION
# ==============================================================================

def generate_primes(low, high):
    primes = []

    if low <= 2 < high:
        primes.append(2)

    start = max(3, low | 1)

    for n in range(start, high, 2):
        if isprime(n):
            primes.append(n)

    return primes


# ==============================================================================
# TARGET GENERATION
# ==============================================================================

def random_target(primes):
    while True:
        p = random.choice(primes)
        q = random.choice(primes)

        if p != q:
            return p, q, p * q


def generate_targets(primes, count):
    targets = []
    seen = set()

    while len(targets) < count:

        p, q, n = random_target(primes)

        pair = tuple(sorted((p, q)))

        if pair in seen:
            continue

        seen.add(pair)
        targets.append((p, q, n))

    return targets


# ==============================================================================
# RESIDUE KEY
# ==============================================================================

def residue_key(x, mods):
    return tuple(x % m for m in mods)


# ==============================================================================
# PRIME RESIDUE INDEX
# ==============================================================================
#
# We store:
#
#     residue vector -> count
#
# rather than:
#
#     residue vector -> list of all primes
#
# for the candidate-count experiment.
#
# This saves substantial memory.
# ==============================================================================

def build_residue_counts(primes, mods):

    counts = Counter()

    for p in primes:
        counts[residue_key(p, mods)] += 1

    return counts


# ==============================================================================
# N-ONLY COMPATIBLE SECOND RESIDUE
# ==============================================================================

def compatible_q_key(n, p_key, mods):
    """
    For each modulus m:

        p*q = n (mod m)

    If p is invertible modulo every m, then

        q = n * p^{-1} (mod m).

    Returns q's full residue vector or None.
    """

    q_key = []

    for a, m in zip(p_key, mods):

        if math.gcd(a, m) != 1:
            return None

        q = (n % m) * pow(a, -1, m) % m
        q_key.append(q)

    return tuple(q_key)


# ==============================================================================
# STREAMING CANDIDATE COUNT
# ==============================================================================

def streaming_candidate_count(
    n,
    primes,
    residue_counts,
    mods,
):
    """
    Count N-only compatible unordered prime pairs without constructing them.

    For each prime p:

        q_residue = n / p mod M

    and residue_counts[q_residue] tells us how many possible q primes
    occupy that class.

    Ordered count is accumulated first.

    We then divide by two because every unordered pair appears once
    with p=a and once with p=b.
    """

    ordered_count = 0
    compatible_prime_count = 0

    for p in primes:

        p_key = residue_key(p, mods)

        q_key = compatible_q_key(
            n,
            p_key,
            mods
        )

        if q_key is None:
            continue

        q_count = residue_counts.get(
            q_key,
            0
        )

        if q_count == 0:
            continue

        ordered_count += q_count
        compatible_prime_count += 1

    # Because p != q for the semiprime targets in this experiment,
    # every unordered pair is represented twice in the ordered sum.
    unordered_count = ordered_count // 2

    return (
        ordered_count,
        unordered_count,
        compatible_prime_count,
    )


# ==============================================================================
# STREAMING EXACT FACTOR RECOVERY
# ==============================================================================

def streaming_exact_recovery(
    n,
    primes,
    prime_set,
    residue_counts,
    mods,
):
    """
    Recover the factorization without enumerating candidate pairs.

    For every prime p whose residue class admits at least one q,
    directly test:

        q = n // p

    when divisible.

    Thus exact recovery costs O(number of primes), not O(number of
    compatible pairs).
    """

    residue_hits = 0
    divisibility_tests = 0

    for p in primes:

        p_key = residue_key(p, mods)

        q_key = compatible_q_key(
            n,
            p_key,
            mods
        )

        if q_key is None:
            continue

        if residue_counts.get(q_key, 0) == 0:
            continue

        residue_hits += 1

        if n % p != 0:
            continue

        divisibility_tests += 1

        q = n // p

        if q == p:
            continue

        if q in prime_set:

            pair = tuple(sorted((p, q)))

            return (
                pair,
                residue_hits,
                divisibility_tests,
            )

    return (
        None,
        residue_hits,
        divisibility_tests,
    )


# ==============================================================================
# TRUE RESIDUE CONTROL
# ==============================================================================

def target_residue_valid(p, q, n, mods):
    pk = residue_key(p, mods)
    qk = residue_key(q, mods)

    wanted = compatible_q_key(
        n,
        pk,
        mods
    )

    return wanted == qk


# ==============================================================================
# CLASS STATISTICS
# ==============================================================================

def class_statistics(residue_counts):

    occupied = len(residue_counts)

    singleton = 0
    max_class = 0
    total_items = 0

    histogram = Counter()

    for count in residue_counts.values():

        total_items += count
        histogram[count] += 1

        if count == 1:
            singleton += 1

        if count > max_class:
            max_class = count

    return (
        occupied,
        singleton,
        max_class,
        histogram,
        total_items,
    )


# ==============================================================================
# HELPER
# ==============================================================================

def median(values):

    values = sorted(values)

    n = len(values)

    if n == 0:
        return 0.0

    if n % 2:
        return float(values[n // 2])

    return (
        values[n // 2 - 1]
        + values[n // 2]
    ) / 2.0


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    total_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 40")
    print("STREAMING N-ONLY CRT RECOVERY")
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

    # ------------------------------------------------------------------
    # MODULUS INVENTORY
    # ------------------------------------------------------------------

    print("=" * 78)
    print("1. MODULUS INVENTORY")
    print("=" * 78)

    for r, m in zip(R_VALUES, MODULI):

        print(
            f"r={r:3d} "
            f"m={m:9,d}"
        )

    print()

    # ------------------------------------------------------------------
    # PRIME POPULATION
    # ------------------------------------------------------------------

    print("=" * 78)
    print("2. PRIME POPULATION")
    print("=" * 78)

    start = time.perf_counter()

    primes = generate_primes(
        PRIME_LOW,
        PRIME_HIGH
    )

    prime_set = set(primes)

    elapsed = time.perf_counter() - start

    pair_count = (
        len(primes) * (len(primes) - 1) // 2
    )

    print(
        f"prime population = "
        f"{len(primes):,}"
    )

    print(
        f"unordered prime pairs = "
        f"{pair_count:,}"
    )

    print(
        f"generation time = "
        f"{elapsed:.3f}s"
    )

    print()

    # ------------------------------------------------------------------
    # TARGETS
    # ------------------------------------------------------------------

    print("=" * 78)
    print("3. TARGETS")
    print("=" * 78)

    targets = generate_targets(
        primes,
        TARGETS
    )

    for i, (p, q, n) in enumerate(
        targets,
        start=1
    ):

        print(
            f"target {i:2d}: "
            f"p={p} "
            f"q={q} "
            f"n={n}"
        )

    print()

    # ------------------------------------------------------------------
    # BUILD RESIDUE COUNTS
    # ------------------------------------------------------------------

    print("=" * 78)
    print("4. RESIDUE-CLASS COUNT TABLES")
    print("=" * 78)

    residue_counts = {}

    for prefix in PREFIXES:

        mods = PREFIX_MODULI[prefix]
        M = PREFIX_M[prefix]

        start = time.perf_counter()

        counts = build_residue_counts(
            primes,
            mods
        )

        elapsed = time.perf_counter() - start

        (
            occupied,
            singleton,
            max_class,
            histogram,
            total,
        ) = class_statistics(counts)

        residue_counts[prefix] = counts

        print(
            f"prefix={prefix:2d} "
            f"M={M:13,d} "
            f"classes={occupied:8,d} "
            f"singleton={singleton:8,d} "
            f"max_class={max_class:4d} "
            f"time={elapsed:.3f}s"
        )

    print()

    # ------------------------------------------------------------------
    # EXPERIMENT
    # ------------------------------------------------------------------

    print("=" * 78)
    print("5. STREAMING N-ONLY RECOVERY")
    print("=" * 78)

    print()
    print(
        "No candidate-pair list is constructed."
    )
    print(
        "Candidate counts are computed by residue-class multiplication."
    )
    print(
        "Exact recovery tests n % p and q = n // p only when the"
    )
    print(
        "CRT residue class contains at least one prime."
    )
    print()

    global_stats = {
        prefix: []
        for prefix in PREFIXES
    }

    for target_index, (p, q, n) in enumerate(
        targets,
        start=1
    ):

        true_pair = tuple(sorted((p, q)))

        print("-" * 78)
        print(
            f"TARGET {target_index:2d} "
            f"(control factors "
            f"{true_pair[0]}, {true_pair[1]})"
        )
        print(f"n = {n}")
        print("-" * 78)

        for prefix in PREFIXES:

            mods = PREFIX_MODULI[prefix]
            M = PREFIX_M[prefix]

            counts = residue_counts[prefix]

            print()
            print(
                f"PREFIX {prefix:2d} "
                f"(M={M:,})"
            )

            # ----------------------------------------------------------
            # Candidate counting
            # ----------------------------------------------------------

            start = time.perf_counter()

            (
                ordered_count,
                unordered_count,
                compatible_prime_count,
            ) = streaming_candidate_count(
                n,
                primes,
                counts,
                mods
            )

            count_time = (
                time.perf_counter() - start
            )

            fraction = (
                unordered_count
                / pair_count
            )

            reference = 1.0 / M

            ratio = (
                fraction / reference
            )

            # ----------------------------------------------------------
            # Exact recovery
            # ----------------------------------------------------------

            start = time.perf_counter()

            (
                recovered,
                residue_hits,
                divisibility_tests,
            ) = streaming_exact_recovery(
                n,
                primes,
                prime_set,
                counts,
                mods
            )

            recovery_time = (
                time.perf_counter() - start
            )

            recovered_ok = (
                recovered == true_pair
            )

            control_valid = target_residue_valid(
                p,
                q,
                n,
                mods
            )

            print(
                f"  ordered candidate count      = "
                f"{ordered_count:,}"
            )

            print(
                f"  unordered candidate count    = "
                f"{unordered_count:,}"
            )

            print(
                f"  candidate fraction            = "
                f"{fraction:.12g}"
            )

            print(
                f"  1/M                           = "
                f"{reference:.12g}"
            )

            print(
                f"  ratio to 1/M                  = "
                f"{ratio:.8g}"
            )

            print(
                f"  primes with compatible q-class= "
                f"{compatible_prime_count:,}"
            )

            print(
                f"  residue-hit primes             = "
                f"{residue_hits:,}"
            )

            print(
                f"  divisibility tests             = "
                f"{divisibility_tests:,}"
            )

            print(
                f"  exact recovery                 = "
                f"{recovered_ok}"
            )

            print(
                f"  true residue control            = "
                f"{control_valid}"
            )

            if recovered is not None:

                print(
                    f"  recovered pair                 = "
                    f"({recovered[0]}, {recovered[1]})"
                )

            print(
                f"  count time                     = "
                f"{count_time:.4f}s"
            )

            print(
                f"  recovery time                  = "
                f"{recovery_time:.4f}s"
            )

            # ----------------------------------------------------------
            # Small-population diagnostics
            # ----------------------------------------------------------

            if unordered_count <= PRINT_EXAMPLES_AT_MOST:

                print(
                    "  NOTE: candidate population is tiny."
                )

            global_stats[prefix].append({
                "candidate_count": unordered_count,
                "fraction": fraction,
                "ratio": ratio,
                "compatible_prime_count":
                    compatible_prime_count,
                "residue_hits":
                    residue_hits,
                "divisibility_tests":
                    divisibility_tests,
                "recovered":
                    recovered_ok,
            })

    # ------------------------------------------------------------------
    # GLOBAL SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. GLOBAL SUMMARY")
    print("=" * 78)
    print()

    print(
        f"{'prefix':>6s}"
        f"{'M':>14s}"
        f"{'median candidates':>22s}"
        f"{'median fraction':>20s}"
        f"{'median ratio':>16s}"
        f"{'recovered':>12s}"
    )

    print("-" * 78)

    for prefix in PREFIXES:

        rows = global_stats[prefix]

        candidates = [
            x["candidate_count"]
            for x in rows
        ]

        fractions = [
            x["fraction"]
            for x in rows
        ]

        ratios = [
            x["ratio"]
            for x in rows
        ]

        recovered = sum(
            x["recovered"]
            for x in rows
        )

        print(
            f"{prefix:6d}"
            f"{PREFIX_M[prefix]:14,d}"
            f"{median(candidates):22,.1f}"
            f"{median(fractions):20.12g}"
            f"{median(ratios):16.8g}"
            f"{recovered:8d}/{TARGETS:<3d}"
        )

    # ------------------------------------------------------------------
    # STREAMING COST
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. STREAMING RECOVERY COST")
    print("=" * 78)
    print()

    for prefix in PREFIXES:

        rows = global_stats[prefix]

        compatible = [
            x["compatible_prime_count"]
            for x in rows
        ]

        residue_hits = [
            x["residue_hits"]
            for x in rows
        ]

        divisions = [
            x["divisibility_tests"]
            for x in rows
        ]

        print(
            f"prefix {prefix:2d}: "
            f"median compatible-primes = "
            f"{median(compatible):,.1f}, "
            f"median residue-hits = "
            f"{median(residue_hits):,.1f}, "
            f"median divisibility-tests = "
            f"{median(divisions):,.1f}"
        )

    # ------------------------------------------------------------------
    # IMPORTANT COMPARISON
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. KEY COMPARISON")
    print("=" * 78)
    print()

    print(
        "The candidate population may still be large, but exact recovery"
    )
    print(
        "does not require enumerating that population."
    )
    print()

    print(
        "The computational distinction is:"
    )

    print()
    print(
        "PAIR ENUMERATION:"
    )
    print(
        "    O(number of compatible prime pairs)"
    )

    print()
    print(
        "STREAMING RECOVERY:"
    )
    print(
        "    O(number of primes × number of moduli)"
    )

    print()
    print(
        "This is the central point of Experiment 40."
    )

    # ------------------------------------------------------------------
    # FINAL
    # ------------------------------------------------------------------

    total_time = (
        time.perf_counter()
        - total_start
    )

    print()
    print("=" * 78)
    print("EXPERIMENT 40 COMPLETE")
    print("=" * 78)

    print(
        f"total runtime = "
        f"{total_time:.3f}s "
        f"({total_time / 60:.2f} min)"
    )


if __name__ == "__main__":
    main()

