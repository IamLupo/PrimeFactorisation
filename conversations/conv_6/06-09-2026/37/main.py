#!/usr/bin/env python3

import math
import time
from itertools import combinations


# ============================================================================
# START EXPERIMENT 193
#
# EMPIRICAL RUNTIME-OPTIMAL QR PRIME SELECTION
#
# Experiment 192 optimized:
#
#     modular density
#
# Experiment 193 optimizes:
#
#     predicted wall-clock runtime
#
# For a prime set S:
#
#     rho(S) = product rho_p
#
# and approximately:
#
#     T(S) =
#         sieve_cost(S)
#         +
#         survivor_count(S) * Fermat_test_cost
#
# We first empirically measure the marking cost of every individual prime.
#
# Then all 2^16 subsets are evaluated.
#
# The best predicted subsets are benchmarked for real.
#
# ============================================================================


FILTER_PRIMES = [
    3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47, 53, 59
]

BLOCK_SIZE = 1 << 20

# Number of y values used to measure sieve costs.
COST_RANGE = 1 << 22

# Actual Fermat benchmark range.
MAX_SEARCH_Y = 25_000_000

# Maximum number of benchmark subsets.
TOP_BENCHMARKS = 12

# Existing instances.
INSTANCES = {
    48: (8390069, 33547589),
    54: (124517461, 144517463),
    60: (1058841403, 1088841421),
    66: (8569934017, 8609934041),
}

# Standard comparison wheels.
STANDARD_WHEELS = [
    [3],
    [3, 5],
    [3, 5, 7],
    [3, 5, 7, 11],
    [3, 5, 7, 11, 13],
    [3, 5, 7, 11, 13, 17],
    [3, 5, 7, 11, 13, 17, 19],
    [3, 5, 7, 11, 13, 17, 19, 23],
]


# ============================================================================
# MILLER-RABIN
# ============================================================================

MR_BASES = [
    2,
    325,
    9375,
    28178,
    450775,
    9780504,
    1795265022,
]


def is_prime(n: int) -> bool:

    if n < 2:
        return False

    small_primes = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37, 41, 43, 47,
        53, 59, 61, 67
    )

    for p in small_primes:
        if n == p:
            return True

        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while (d & 1) == 0:
        d >>= 1
        s += 1

    for a in MR_BASES:

        if a % n == 0:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):

            x = (x * x) % n

            if x == n - 1:
                break

        else:
            return False

    return True


# ============================================================================
# LEGENDRE / DENSITY
# ============================================================================

def legendre_symbol(a: int, p: int) -> int:

    a %= p

    if a == 0:
        return 0

    x = pow(a, (p - 1) // 2, p)

    if x == 1:
        return 1

    if x == p - 1:
        return -1

    raise RuntimeError("invalid Legendre symbol")


def prime_density(n: int, p: int):

    chi = legendre_symbol(-n, p)

    if n % p == 0:
        allowed = (p + 1) // 2
    else:
        allowed = (p + chi) // 2

    return allowed / p


def subset_density(n: int, subset):

    density = 1.0

    for p in subset:
        density *= prime_density(n, p)

    return density


# ============================================================================
# BAD RESIDUES
# ============================================================================

def build_bad_residues(n: int, p: int):

    qr = {
        (x * x) % p
        for x in range(p)
    }

    return [
        r
        for r in range(p)
        if (n + r * r) % p not in qr
    ]


# ============================================================================
# SINGLE-PRIME MARKING COST
# ============================================================================

def measure_prime_cost(n: int, p: int):

    """
    Measure the actual cost of applying ONE prime to a complete segment
    range.

    The Fermat square-root test is NOT performed here.
    """

    bad = build_bad_residues(n, p)

    start = time.perf_counter()

    blocks = 0
    total_zero = 0

    block_start = 0

    while block_start < COST_RANGE:

        block_end = min(
            block_start + BLOCK_SIZE,
            COST_RANGE
        )

        block_len = block_end - block_start

        alive = bytearray(
            b"\x01"
        ) * block_len

        for residue in bad:

            first = (
                block_start
                + ((residue - block_start) % p)
            )

            if first >= block_end:
                continue

            idx = first - block_start

            count = (
                (block_len - 1 - idx) // p
            ) + 1

            alive[idx:block_len:p] = (
                b"\x00"
            ) * count

        total_zero += block_len - alive.count(1)

        blocks += 1
        block_start = block_end

    elapsed = time.perf_counter() - start

    return elapsed, total_zero, bad


# ============================================================================
# FERMAT TEST COST
# ============================================================================

def measure_fermat_test_cost(n: int, count: int):

    """
    Measure average cost of the actual candidate test.

    We deliberately use nonsolutions.
    """

    start = time.perf_counter()

    hits = 0

    for y in range(count):

        value = n + y * y

        x = math.isqrt(value)

        if x * x == value:
            hits += 1

    elapsed = time.perf_counter() - start

    return elapsed / count, hits


# ============================================================================
# EXHAUSTIVE SUBSET ENUMERATION
# ============================================================================

def enumerate_subsets(
    n: int,
    y_limit: int,
    prime_costs,
    test_cost,
):
    """
    Enumerate all 65536 subsets.

    Cost model:

        predicted =
            sum(single-prime marking costs)
            +
            (y_limit + 1) * rho(S) * test_cost
    """

    items = FILTER_PRIMES
    count = len(items)

    results = []

    for mask in range(1 << count):

        density = 1.0
        sieve_cost = 0.0
        subset = []

        for i, p in enumerate(items):

            if mask & (1 << i):

                subset.append(p)

                density *= prime_density(n, p)

                sieve_cost += prime_costs[p]

        predicted_survivors = (
            (y_limit + 1) * density
        )

        predicted_total = (
            sieve_cost
            +
            predicted_survivors * test_cost
        )

        results.append(
            (
                predicted_total,
                density,
                sieve_cost,
                subset,
                predicted_survivors,
            )
        )

    results.sort(key=lambda x: x[0])

    return results


# ============================================================================
# ACTUAL SEGMENTED FERMAT
# ============================================================================

def segmented_fermat(
    n: int,
    primes,
    max_y: int,
):

    bad_by_prime = {
        p: build_bad_residues(n, p)
        for p in primes
    }

    survivors = 0

    start = time.perf_counter()

    block_start = 0

    while block_start <= max_y:

        block_end = min(
            block_start + BLOCK_SIZE,
            max_y + 1
        )

        block_len = block_end - block_start

        alive = bytearray(
            b"\x01"
        ) * block_len

        for p in primes:

            for residue in bad_by_prime[p]:

                first = (
                    block_start
                    + ((residue - block_start) % p)
                )

                if first >= block_end:
                    continue

                idx = first - block_start

                count = (
                    (block_len - 1 - idx) // p
                ) + 1

                alive[idx:block_len:p] = (
                    b"\x00"
                ) * count

        pos = 0

        while True:

            pos = alive.find(1, pos)

            if pos < 0:
                break

            y = block_start + pos

            survivors += 1

            value = n + y * y
            x = math.isqrt(value)

            if x * x == value:

                p = x - y
                q = x + y

                if p > 1 and q > 1 and p * q == n:

                    elapsed = time.perf_counter() - start

                    return (
                        y,
                        (p, q),
                        survivors,
                        elapsed,
                    )

            pos += 1

        block_start = block_end

    elapsed = time.perf_counter() - start

    return (
        None,
        None,
        survivors,
        elapsed,
    )


# ============================================================================
# STANDARD BENCHMARK
# ============================================================================

def benchmark_wheel(
    n: int,
    true_y: int,
    p: int,
    q: int,
    wheel,
):

    found_y, result, survivors, elapsed = (
        segmented_fermat(
            n,
            wheel,
            MAX_SEARCH_Y,
        )
    )

    correct = (
        result == (p, q)
        or result == (q, p)
    )

    return {
        "wheel": wheel,
        "density": subset_density(n, wheel),
        "survivors": survivors,
        "found_y": found_y,
        "result": result,
        "correct": correct,
        "time": elapsed,
    }


# ============================================================================
# MAIN INSTANCE
# ============================================================================

def run_instance(bits: int, p: int, q: int):

    n = p * q
    true_y = (q - p) // 2

    print()
    print("=" * 72)
    print(f"START INSTANCE {bits}-BIT")
    print("=" * 72)

    print(f"p = {p}")
    print(f"q = {q}")
    print(f"true y = {true_y}")
    print(f"MAX_SEARCH_Y = {MAX_SEARCH_Y}")

    # ------------------------------------------------------------------------
    # Skip impossible recovery benchmarks.
    # ------------------------------------------------------------------------

    if true_y > MAX_SEARCH_Y:

        print()
        print(
            "    RECOVERY SKIPPED: true y exceeds MAX_SEARCH_Y"
        )

        print()
        print(f"FINISHED INSTANCE {bits}-BIT")

        return

    # ------------------------------------------------------------------------
    # Measure Fermat-test cost.
    # ------------------------------------------------------------------------

    test_samples = min(
        true_y + 1,
        1_000_000
    )

    test_cost, hits = measure_fermat_test_cost(
        n,
        test_samples,
    )

    print()
    print("FERMAT TEST COST")
    print(f"    samples = {test_samples}")
    print(f"    average = {test_cost:.12e}s")
    print(f"    accidental squares = {hits}")

    # ------------------------------------------------------------------------
    # Measure every prime.
    # ------------------------------------------------------------------------

    print()
    print("INDIVIDUAL PRIME COSTS")

    prime_costs = {}

    for prime in FILTER_PRIMES:

        elapsed, zero_count, bad = (
            measure_prime_cost(
                n,
                prime,
            )
        )

        prime_costs[prime] = elapsed

        rho = prime_density(
            n,
            prime,
        )

        print(
            f"prime={prime:2d} "
            f"density={rho:.12f} "
            f"bad_classes={len(bad):2d} "
            f"cost={elapsed:.8f}s"
        )

    # ------------------------------------------------------------------------
    # Exhaustive subset optimization.
    # ------------------------------------------------------------------------

    print()
    print("EXHAUSTIVE SUBSET OPTIMIZATION")

    results = enumerate_subsets(
        n,
        true_y,
        prime_costs,
        test_cost,
    )

    print(f"    subsets evaluated = {len(results)}")

    print()
    print("TOP PREDICTED SETS")

    for rank, item in enumerate(
        results[:TOP_BENCHMARKS],
        start=1,
    ):

        predicted_total, density, sieve_cost, subset, predicted_survivors = item

        print()
        print(f"RANK {rank}")
        print(f"    primes = {subset}")
        print(f"    density = {density:.12f}")
        print(f"    predicted survivors = {predicted_survivors:.2f}")
        print(f"    predicted sieve cost = {sieve_cost:.8f}s")
        print(f"    predicted total = {predicted_total:.8f}s")
        print(
            f"    reduction factor = "
            f"{1.0 / density:.3f}x"
        )

    # ------------------------------------------------------------------------
    # Benchmark selected subsets.
    # ------------------------------------------------------------------------

    print()
    print("ACTUAL BENCHMARKS")

    candidate_wheels = []

    for wheel in STANDARD_WHEELS:

        if wheel not in candidate_wheels:

            candidate_wheels.append(wheel)

    for item in results[:TOP_BENCHMARKS]:

        wheel = item[3]

        if wheel not in candidate_wheels:

            candidate_wheels.append(wheel)

    actual_results = []

    for wheel in candidate_wheels:

        result = benchmark_wheel(
            n,
            true_y,
            p,
            q,
            wheel,
        )

        actual_results.append(result)

        print()
        print(f"WHEEL {wheel}")
        print(
            f"    density = "
            f"{result['density']:.12f}"
        )
        print(
            f"    survivors = "
            f"{result['survivors']}"
        )
        print(
            f"    found y = "
            f"{result['found_y']}"
        )
        print(
            f"    result = "
            f"{result['result']}"
        )
        print(
            f"    correct = "
            f"{result['correct']}"
        )
        print(
            f"    actual time = "
            f"{result['time']:.8f}s"
        )

    # ------------------------------------------------------------------------
    # Best actual wheel.
    # ------------------------------------------------------------------------

    successful = [
        r
        for r in actual_results
        if r["correct"]
    ]

    if successful:

        successful.sort(
            key=lambda r: r["time"]
        )

        best = successful[0]

        print()
        print("BEST ACTUAL WHEEL")
        print(
            f"    primes = "
            f"{best['wheel']}"
        )
        print(
            f"    time = "
            f"{best['time']:.8f}s"
        )
        print(
            f"    survivors = "
            f"{best['survivors']}"
        )
        print(
            f"    density = "
            f"{best['density']:.12f}"
        )

    print()
    print(f"FINISHED INSTANCE {bits}-BIT")


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("START EXPERIMENT 193")
    print()
    print("Empirical runtime-optimal QR prime selection")
    print()
    print(f"FILTER_PRIMES = {FILTER_PRIMES}")
    print(f"BLOCK_SIZE = {BLOCK_SIZE}")
    print(f"COST_RANGE = {COST_RANGE}")
    print(f"MAX_SEARCH_Y = {MAX_SEARCH_Y}")
    print(f"TOP_BENCHMARKS = {TOP_BENCHMARKS}")

    # ------------------------------------------------------------------------
    # 48-bit reference first.
    # ------------------------------------------------------------------------

    p, q = INSTANCES[48]

    run_instance(
        48,
        p,
        q,
    )

    # ------------------------------------------------------------------------
    # Additional recovery instances.
    # ------------------------------------------------------------------------

    for bits in [54, 60, 66]:

        p, q = INSTANCES[bits]

        run_instance(
            bits,
            p,
            q,
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 193")
    print("=" * 72)


if __name__ == "__main__":
    main()
