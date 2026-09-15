#!/usr/bin/env python3


# ============================================================================
# START EXPERIMENT 196
#
# EMPIRICAL INCREMENTAL QR PRIME SELECTION
#
# Experiment 195 showed that:
#
#     density
#
# and
#
#     simple workload = sum(bad/p)
#
# are not sufficient to predict real runtime.
#
# Experiment 196 directly measures:
#
#     cost(S + {p}) - cost(S)
#
# for every possible next prime p.
#
# At each stage:
#
#     current set S
#
# is fixed.
#
# Every unused prime is tested by actually running the segmented marking
# operation on the SAME measurement interval.
#
# We record:
#
#     - actual sieve time
#     - survivor count
#     - density
#     - incremental time
#     - incremental survivors removed
#     - survivors removed per microsecond
#
# We then construct an empirical greedy chain:
#
#     S0 = {}
#     S1 = S0 + best prime
#     S2 = S1 + best prime
#     ...
#
# where "best" is determined by actual measured runtime on the measurement
# interval, NOT by analytical density.
#
# Finally, the best members of the greedy chain and the best alternatives are
# benchmarked over the full Fermat search interval.
#
# ============================================================================

import math
import time


# ============================================================================
# CONFIGURATION
# ============================================================================

PRIMES = [
    3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47, 53, 59
]

BLOCK_SIZE = 1 << 20

# Keep calibration small enough that the experiment remains practical.
#
# 4 million y-values = 4 blocks.
MEASURE_RANGE = 1 << 22

# Full Fermat benchmark.
MAX_Y = 25_000_000

# Only build the greedy chain to this depth.
GREEDY_DEPTH = 12

# At each greedy step, benchmark the best N choices from that step in the
# final/full-range test.
TOP_ALTERNATIVES_PER_STEP = 2

# Existing test instances.
INSTANCES = {
    48: (8390069, 33547589),
    54: (124517461, 144517463),
    60: (1058841403, 1088841421),
    66: (8569934017, 8609934041),
}


# ============================================================================
# LEGENDRE / ANALYTICAL DENSITY
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

    raise RuntimeError(
        f"Invalid Legendre symbol for p={p}"
    )


def prime_density(n: int, p: int) -> float:

    chi = legendre_symbol(-n, p)

    if n % p == 0:
        allowed = (p + 1) // 2
    else:
        allowed = (p + chi) // 2

    return allowed / p


def subset_density(n: int, subset) -> float:

    rho = 1.0

    for p in subset:
        rho *= prime_density(n, p)

    return rho


# ============================================================================
# BAD RESIDUE TABLES
# ============================================================================

def build_bad_table(n: int):

    bad = {}

    for p in PRIMES:

        qr = {
            (x * x) % p
            for x in range(p)
        }

        bad[p] = [
            r
            for r in range(p)
            if (
                n + r * r
            ) % p not in qr
        ]

    return bad


# ============================================================================
# APPLY A SIEVE TO ONE MEASUREMENT RANGE
# ============================================================================

def measure_subset(
    subset,
    bad_table,
    limit,
):
    """
    Apply exactly the same bytearray marking implementation used in the
    actual Fermat search.

    This is deliberately NOT an abstract cost model.
    """

    start = time.perf_counter()

    survivors = 0

    block_start = 0

    while block_start < limit:

        block_end = min(
            block_start + BLOCK_SIZE,
            limit
        )

        block_len = block_end - block_start

        alive = bytearray(
            b"\x01"
        ) * block_len

        for p in subset:

            for residue in bad_table[p]:

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

        survivors += alive.count(1)

        block_start = block_end

    elapsed = time.perf_counter() - start

    return elapsed, survivors


# ============================================================================
# ACTUAL FERMAT TEST
# ============================================================================

def fermat_test_y(n: int, y: int):

    value = n + y * y

    x = math.isqrt(value)

    if x * x != value:
        return None

    p = x - y
    q = x + y

    if p > 1 and q > 1 and p * q == n:
        return p, q

    return None


# ============================================================================
# FULL SEGMENTED FERMAT
# ============================================================================

def full_fermat_search(
    n: int,
    subset,
    bad_table,
    max_y,
):

    start = time.perf_counter()

    survivors = 0

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

        for p in subset:

            for residue in bad_table[p]:

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

            pos = alive.find(
                1,
                pos
            )

            if pos < 0:
                break

            y = block_start + pos

            survivors += 1

            result = fermat_test_y(
                n,
                y
            )

            if result is not None:

                elapsed = (
                    time.perf_counter()
                    - start
                )

                return {
                    "time": elapsed,
                    "survivors": survivors,
                    "found_y": y,
                    "result": result,
                }

            pos += 1

        block_start = block_end

    elapsed = time.perf_counter() - start

    return {
        "time": elapsed,
        "survivors": survivors,
        "found_y": None,
        "result": None,
    }


# ============================================================================
# MEASURE INDIVIDUAL PRIME
# ============================================================================

def measure_single_primes(
    bad_table,
):

    print()
    print("=" * 72)
    print("INDIVIDUAL PRIME CALIBRATION")
    print("=" * 72)

    results = {}

    for p in PRIMES:

        elapsed, survivors = measure_subset(
            (p,),
            bad_table,
            MEASURE_RANGE
        )

        results[p] = {
            "time": elapsed,
            "survivors": survivors,
        }

        density = survivors / MEASURE_RANGE

        print(
            f"prime={p:2d} "
            f"density={density:.12f} "
            f"survivors={survivors} "
            f"time={elapsed:.8f}s"
        )

    return results


# ============================================================================
# MEASURE ALL POSSIBLE EXTENSIONS
# ============================================================================

def measure_extensions(
    current,
    bad_table,
    baseline_time,
    baseline_survivors,
):

    unused = [
        p
        for p in PRIMES
        if p not in current
    ]

    results = []

    for p in unused:

        candidate = tuple(
            sorted(
                current + (p,)
            )
        )

        elapsed, survivors = measure_subset(
            candidate,
            bad_table,
            MEASURE_RANGE
        )

        incremental_time = (
            elapsed - baseline_time
        )

        survivors_removed = (
            baseline_survivors - survivors
        )

        if incremental_time > 0:
            removal_per_microsecond = (
                survivors_removed
                / (incremental_time * 1_000_000)
            )
        else:
            removal_per_microsecond = float("inf")

        results.append(
            {
                "prime": p,
                "subset": candidate,
                "time": elapsed,
                "survivors": survivors,
                "incremental_time": incremental_time,
                "survivors_removed": survivors_removed,
                "removal_per_microsecond":
                    removal_per_microsecond,
            }
        )

    return results


# ============================================================================
# CHOOSE BEST EXTENSION
# ============================================================================

def choose_best(results):

    """
    Primary criterion:
        minimum actual measured time.

    Secondary:
        minimum survivor count.
    """

    return min(
        results,
        key=lambda r: (
            r["time"],
            r["survivors"]
        )
    )


# ============================================================================
# GREEDY CHAIN
# ============================================================================

def build_greedy_chain(
    n: int,
    bad_table,
):

    print()
    print("=" * 72)
    print("EMPIRICAL GREEDY CHAIN")
    print("=" * 72)

    chain = []

    current = tuple()

    baseline_time, baseline_survivors = (
        measure_subset(
            current,
            bad_table,
            MEASURE_RANGE
        )
    )

    chain.append(
        {
            "depth": 0,
            "subset": current,
            "time": baseline_time,
            "survivors": baseline_survivors,
        }
    )

    print()
    print(
        f"DEPTH 0 "
        f"primes=[] "
        f"time={baseline_time:.8f}s "
        f"survivors={baseline_survivors}"
    )

    for depth in range(
        1,
        min(GREEDY_DEPTH, len(PRIMES)) + 1
    ):

        extension_results = measure_extensions(
            current,
            bad_table,
            baseline_time,
            baseline_survivors,
        )

        print()
        print(
            f"STEP {depth}: "
            f"CURRENT = {list(current)}"
        )

        ranked = sorted(
            extension_results,
            key=lambda r: (
                r["time"],
                r["survivors"]
            )
        )

        for rank, result in enumerate(
            ranked[:5],
            start=1
        ):

            print(
                f"    OPTION {rank}: "
                f"+{result['prime']:2d} "
                f"time={result['time']:.8f}s "
                f"survivors={result['survivors']} "
                f"delta_time="
                f"{result['incremental_time']:.8f}s "
                f"removed="
                f"{result['survivors_removed']} "
                f"removed/us="
                f"{result['removal_per_microsecond']:.2f}"
            )

        best = choose_best(
            extension_results
        )

        current = best["subset"]

        baseline_time = best["time"]
        baseline_survivors = best["survivors"]

        chain.append(
            {
                "depth": depth,
                "subset": current,
                "time": baseline_time,
                "survivors": baseline_survivors,
            }
        )

        print(
            f"    SELECTED +{best['prime']}"
        )
        print(
            f"    NEW SET = {list(current)}"
        )
        print(
            f"    time = {baseline_time:.8f}s"
        )
        print(
            f"    survivors = {baseline_survivors}"
        )

    return chain


# ============================================================================
# ANALYTICAL BEST-SURVIVOR CHAIN
# ============================================================================

def analytical_chain(
    n: int,
):

    current = []

    result = [
        (
            0,
            tuple(),
            1.0,
        )
    ]

    for depth in range(
        1,
        min(GREEDY_DEPTH, len(PRIMES)) + 1
    ):

        best = None

        for p in PRIMES:

            if p in current:
                continue

            candidate = tuple(
                sorted(
                    current + [p]
                )
            )

            density = subset_density(
                n,
                candidate
            )

            if best is None or density < best[0]:
                best = (
                    density,
                    p,
                    candidate,
                )

        density, p, candidate = best

        current = list(candidate)

        result.append(
            (
                depth,
                candidate,
                density,
            )
        )

    return result


# ============================================================================
# FULL BENCHMARK
# ============================================================================

def benchmark_subset(
    n: int,
    p: int,
    q: int,
    subset,
    bad_table,
):

    result = full_fermat_search(
        n,
        subset,
        bad_table,
        MAX_Y,
    )

    correct = (
        result["result"] == (p, q)
        or
        result["result"] == (q, p)
    )

    return {
        "subset": subset,
        "time": result["time"],
        "survivors": result["survivors"],
        "found_y": result["found_y"],
        "result": result["result"],
        "correct": correct,
    }


# ============================================================================
# BENCHMARK GREEDY CHAIN
# ============================================================================

def benchmark_chain(
    n: int,
    p: int,
    q: int,
    chain,
    bad_table,
):

    print()
    print("=" * 72)
    print("FULL-RANGE GREEDY CHAIN BENCHMARK")
    print("=" * 72)

    results = []

    for item in chain:

        subset = item["subset"]

        result = benchmark_subset(
            n,
            p,
            q,
            subset,
            bad_table,
        )

        results.append(result)

        print()
        print(
            f"DEPTH {item['depth']}"
        )
        print(
            f"    primes = {list(subset)}"
        )
        print(
            f"    calibration time = "
            f"{item['time']:.8f}s"
        )
        print(
            f"    calibration survivors = "
            f"{item['survivors']}"
        )
        print(
            f"    full survivors = "
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
            f"    full time = "
            f"{result['time']:.8f}s"
        )

    return results


# ============================================================================
# BENCHMARK LOCAL ALTERNATIVES
# ============================================================================

def benchmark_step_alternatives(
    n: int,
    p: int,
    q: int,
    chain,
    bad_table,
):

    print()
    print("=" * 72)
    print("LOCAL ALTERNATIVE BENCHMARKS")
    print("=" * 72)

    all_results = []

    current = tuple()

    baseline_time, baseline_survivors = (
        measure_subset(
            current,
            bad_table,
            MEASURE_RANGE
        )
    )

    for depth in range(
        1,
        min(
            GREEDY_DEPTH,
            len(PRIMES)
        ) + 1
    ):

        alternatives = measure_extensions(
            current,
            bad_table,
            baseline_time,
            baseline_survivors,
        )

        alternatives.sort(
            key=lambda x: x["time"]
        )

        selected = alternatives[
            :TOP_ALTERNATIVES_PER_STEP
        ]

        print()
        print(
            f"STEP {depth} "
            f"FROM {list(current)}"
        )

        actual = []

        for candidate in selected:

            result = benchmark_subset(
                n,
                p,
                q,
                candidate["subset"],
                bad_table,
            )

            actual.append(result)

            print(
                f"    candidate = "
                f"{list(candidate['subset'])}"
            )
            print(
                f"        calibration time = "
                f"{candidate['time']:.8f}s"
            )
            print(
                f"        survivors = "
                f"{result['survivors']}"
            )
            print(
                f"        full time = "
                f"{result['time']:.8f}s"
            )
            print(
                f"        correct = "
                f"{result['correct']}"
            )

        all_results.extend(
            actual
        )

        # Follow the empirically best calibration choice.
        best = choose_best(
            alternatives
        )

        current = best["subset"]

        baseline_time = best["time"]
        baseline_survivors = best["survivors"]

    return all_results


# ============================================================================
# INSTANCE
# ============================================================================

def run_instance(
    bits: int,
    p: int,
    q: int,
):

    n = p * q

    true_y = (q - p) // 2

    print()
    print("#" * 72)
    print(f"START INSTANCE {bits}-BIT")
    print("#" * 72)

    print(f"p = {p}")
    print(f"q = {q}")
    print(f"true y = {true_y}")
    print(f"MAX_Y = {MAX_Y}")

    if true_y > MAX_Y:

        print()
        print(
            "SKIPPED: true y exceeds MAX_Y"
        )

        print(
            f"FINISHED INSTANCE {bits}-BIT"
        )

        return

    bad_table = build_bad_table(n)

    # ------------------------------------------------------------------------
    # Analytical information.
    # ------------------------------------------------------------------------

    print()
    print("PRIME ANALYTICAL DATA")

    for prime in PRIMES:

        density = prime_density(
            n,
            prime
        )

        print(
            f"prime={prime:2d} "
            f"density={density:.12f}"
        )

    # ------------------------------------------------------------------------
    # Individual calibration.
    # ------------------------------------------------------------------------

    measure_single_primes(
        bad_table
    )

    # ------------------------------------------------------------------------
    # Greedy chain.
    # ------------------------------------------------------------------------

    chain = build_greedy_chain(
        n,
        bad_table,
    )

    # ------------------------------------------------------------------------
    # Analytical density chain for comparison.
    # ------------------------------------------------------------------------

    density_chain = analytical_chain(
        n
    )

    print()
    print("=" * 72)
    print("ANALYTICAL DENSITY GREEDY CHAIN")
    print("=" * 72)

    for depth, subset, density in density_chain:

        print(
            f"DEPTH {depth:2d}: "
            f"primes={list(subset)} "
            f"density={density:.12f}"
        )

    # ------------------------------------------------------------------------
    # Full benchmark of empirical chain.
    # ------------------------------------------------------------------------

    benchmark_chain(
        n,
        p,
        q,
        chain,
        bad_table,
    )

    # ------------------------------------------------------------------------
    # Benchmark local alternatives.
    # ------------------------------------------------------------------------

    benchmark_step_alternatives(
        n,
        p,
        q,
        chain,
        bad_table,
    )

    print()
    print(f"FINISHED INSTANCE {bits}-BIT")


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("START EXPERIMENT 196")
    print()
    print(
        "Empirical incremental QR prime selection"
    )
    print()

    print(
        f"PRIMES = {PRIMES}"
    )
    print(
        f"BLOCK_SIZE = {BLOCK_SIZE}"
    )
    print(
        f"MEASURE_RANGE = {MEASURE_RANGE}"
    )
    print(
        f"MAX_Y = {MAX_Y}"
    )
    print(
        f"GREEDY_DEPTH = {GREEDY_DEPTH}"
    )

    # Start with 48-bit. This gives enough room for a meaningful full search
    # while remaining fast enough to benchmark many candidate sets.
    run_instance(
        48,
        *INSTANCES[48],
    )

    # Then test 54-bit using the controlled-gap instance.
    run_instance(
        54,
        *INSTANCES[54],
    )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 196")
    print("=" * 72)


if __name__ == "__main__":
    main()
