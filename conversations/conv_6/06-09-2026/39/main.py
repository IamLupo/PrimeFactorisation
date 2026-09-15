#!/usr/bin/env python3

import math
import time


# ============================================================================
# START EXPERIMENT 195
#
# QR SIEVE PARETO OPTIMIZATION
#
# Experiment 194:
#
#   pairwise regression was underdetermined:
#
#       137 coefficients
#       100 samples
#
# resulting in meaningless negative runtime predictions.
#
# Experiment 195 removes regression completely.
#
# For every prime p:
#
#   bad_classes[p] = number of y residues eliminated modulo p
#
# For a segmented block of length B, each bad residue marks approximately
#
#       B / p
#
# positions.
#
# Therefore a physically motivated marking-work estimate is:
#
#       W(S) = sum( bad_classes[p] / p )
#
# over p in S.
#
# We combine:
#
#       1. exact analytical survivor density
#       2. exact estimated number of Fermat tests
#       3. measurable marking workload
#
# and identify Pareto-optimal subsets.
#
# Then we benchmark the Pareto frontier directly.
#
# The goal is to determine whether:
#
#       density + actual marking workload
#
# predicts real runtime better than the previous additive model.
#
# ============================================================================


PRIMES = [
    3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47, 53, 59
]

BLOCK_SIZE = 1 << 20

MAX_Y = 25_000_000

# Number of top Pareto candidates to benchmark.
MAX_BENCHMARKS = 30

# Existing benchmark instances.
INSTANCES = {
    48: (8390069, 33547589),
    54: (124517461, 144517463),
    60: (1058841403, 1088841421),
    66: (8569934017, 8609934041),
}


# ============================================================================
# LEGENDRE SYMBOL
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
        f"Invalid Legendre result {x}"
    )


# ============================================================================
# ANALYTICAL DENSITY
# ============================================================================

def prime_density(n: int, p: int):

    chi = legendre_symbol(-n, p)

    if n % p == 0:
        allowed = (p + 1) // 2
    else:
        allowed = (p + chi) // 2

    bad = p - allowed

    density = allowed / p

    return density, allowed, bad


def subset_density(n: int, subset):

    rho = 1.0

    for p in subset:
        rho *= prime_density(n, p)[0]

    return rho


# ============================================================================
# PRIME DATA
# ============================================================================

def build_prime_data(n: int):

    data = {}

    for p in PRIMES:

        density, allowed, bad = prime_density(
            n,
            p
        )

        # Important:
        #
        # Each bad residue marks approximately B/p positions.
        #
        # Thus bad/p is a first-order estimate of the marking work.
        workload = bad / p

        data[p] = {
            "density": density,
            "allowed": allowed,
            "bad": bad,
            "workload": workload,
        }

    return data


# ============================================================================
# FERMAT TEST COST
# ============================================================================

def measure_fermat_cost(n: int):

    samples = min(
        MAX_Y + 1,
        1_000_000
    )

    start = time.perf_counter()

    accidental = 0

    for y in range(samples):

        value = n + y * y

        x = math.isqrt(value)

        if x * x == value:
            accidental += 1

    elapsed = time.perf_counter() - start

    return elapsed / samples, accidental


# ============================================================================
# MARKING-WORK MODEL
# ============================================================================

def workload_for_subset(
    subset,
    prime_data,
):

    return sum(
        prime_data[p]["workload"]
        for p in subset
    )


def expected_survivors(
    y_limit: int,
    density: float,
):

    return (y_limit + 1) * density


# ============================================================================
# SCORE
# ============================================================================

def analytical_record(
    n: int,
    subset,
    prime_data,
    fermat_test_cost,
):

    density = 1.0
    workload = 0.0

    bad_classes = 0

    for p in subset:

        density *= prime_data[p]["density"]

        workload += prime_data[p]["workload"]

        bad_classes += prime_data[p]["bad"]

    survivors = expected_survivors(
        MAX_Y,
        density
    )

    # Estimated search/test component.
    fermat_work = (
        survivors
        * fermat_test_cost
    )

    # We intentionally do NOT add workload and Fermat time directly here.
    # Their physical scaling differs.
    #
    # Instead, this record is used for Pareto analysis.

    return {
        "subset": tuple(subset),
        "density": density,
        "survivors": survivors,
        "workload": workload,
        "bad_classes": bad_classes,
        "fermat_time": fermat_work,
    }


# ============================================================================
# PARETO FRONTIER
# ============================================================================

def dominates(a, b):
    """
    A dominates B if A is no worse in both:
        workload
        survivor count

    and strictly better in at least one.
    """

    no_worse = (
        a["workload"] <= b["workload"]
        and
        a["survivors"] <= b["survivors"]
    )

    strictly_better = (
        a["workload"] < b["workload"]
        or
        a["survivors"] < b["survivors"]
    )

    return no_worse and strictly_better


def pareto_frontier(records):

    """
    Exact O(N^2) would be unnecessary for 65536 records.

    Sort by workload.

    For increasing workload, a point is Pareto-optimal if its survivor
    count is lower than everything seen before.
    """

    ordered = sorted(
        records,
        key=lambda x: (
            x["workload"],
            x["survivors"]
        )
    )

    frontier = []

    best_survivors = float("inf")

    for record in ordered:

        if record["survivors"] < best_survivors:

            frontier.append(record)

            best_survivors = record["survivors"]

    return frontier


# ============================================================================
# ALL SUBSETS
# ============================================================================

def enumerate_all_subsets(
    n: int,
    prime_data,
    fermat_test_cost,
):

    records = []

    for mask in range(
        1 << len(PRIMES)
    ):

        subset = []

        density = 1.0
        workload = 0.0
        bad_classes = 0

        for i, p in enumerate(PRIMES):

            if mask & (1 << i):

                subset.append(p)

                density *= (
                    prime_data[p]["density"]
                )

                workload += (
                    prime_data[p]["workload"]
                )

                bad_classes += (
                    prime_data[p]["bad"]
                )

        survivors = (
            MAX_Y + 1
        ) * density

        fermat_component = (
            survivors
            * fermat_test_cost
        )

        records.append(
            {
                "subset": tuple(subset),
                "density": density,
                "survivors": survivors,
                "workload": workload,
                "bad_classes": bad_classes,
                "fermat_time": fermat_component,
            }
        )

    return records


# ============================================================================
# DIFFERENT PRACTICAL RANKINGS
# ============================================================================

def rank_by_density(records):

    return sorted(
        records,
        key=lambda x: x["density"]
    )


def rank_by_survivors(records):

    return sorted(
        records,
        key=lambda x: x["survivors"]
    )


def rank_by_workload(records):

    return sorted(
        records,
        key=lambda x: x["workload"]
    )


# ============================================================================
# NORMALIZED WORK SCORE
# ============================================================================

def combined_score(
    record,
    density_weight: float,
):

    """
    A dimensionless score.

    survivor count is normalized relative to MAX_Y.
    workload is already O(number of modular marking operations).

    This is only a ranking device; the actual result is determined by
    direct benchmarks.
    """

    survivor_fraction = (
        record["survivors"]
        / (MAX_Y + 1)
    )

    return (
        record["workload"]
        +
        density_weight * survivor_fraction
    )


# ============================================================================
# BUILD BAD RESIDUE TABLES
# ============================================================================

def build_bad_residue_table(
    n: int,
    subset,
):

    table = {}

    for p in subset:

        qr = {
            (x * x) % p
            for x in range(p)
        }

        table[p] = [
            r
            for r in range(p)
            if (
                n + r * r
            ) % p not in qr
        ]

    return table


# ============================================================================
# ACTUAL SEGMENTED FERMAT
# ============================================================================

def segmented_fermat(
    n: int,
    subset,
    bad_table,
):

    start = time.perf_counter()

    survivors = 0

    block_start = 0

    while block_start <= MAX_Y:

        block_end = min(
            block_start + BLOCK_SIZE,
            MAX_Y + 1
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
                    (block_len - 1 - idx)
                    // p
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

            value = n + y * y

            x = math.isqrt(value)

            if x * x == value:

                p = x - y
                q = x + y

                if (
                    p > 1
                    and q > 1
                    and p * q == n
                ):

                    elapsed = (
                        time.perf_counter()
                        - start
                    )

                    return {
                        "found_y": y,
                        "result": (p, q),
                        "survivors": survivors,
                        "time": elapsed,
                    }

            pos += 1

        block_start = block_end

    elapsed = (
        time.perf_counter()
        - start
    )

    return {
        "found_y": None,
        "result": None,
        "survivors": survivors,
        "time": elapsed,
    }


# ============================================================================
# BENCHMARK
# ============================================================================

def benchmark_subset(
    n: int,
    p: int,
    q: int,
    subset,
):

    bad_table = build_bad_residue_table(
        n,
        subset
    )

    result = segmented_fermat(
        n,
        subset,
        bad_table
    )

    correct = (
        result["result"] == (p, q)
        or
        result["result"] == (q, p)
    )

    return {
        "subset": tuple(subset),
        "time": result["time"],
        "survivors": result["survivors"],
        "found_y": result["found_y"],
        "result": result["result"],
        "correct": correct,
    }


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
    print("=" * 72)
    print(f"START INSTANCE {bits}-BIT")
    print("=" * 72)

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

    # ------------------------------------------------------------------------
    # Fermat test cost.
    # ------------------------------------------------------------------------

    test_cost, accidental = (
        measure_fermat_cost(n)
    )

    print()
    print("FERMAT TEST COST")
    print(
        f"    cost = "
        f"{test_cost:.12e}s"
    )
    print(
        f"    accidental squares = "
        f"{accidental}"
    )

    # ------------------------------------------------------------------------
    # Prime data.
    # ------------------------------------------------------------------------

    prime_data = build_prime_data(n)

    print()
    print("PRIME DATA")

    for p0 in PRIMES:

        d = prime_data[p0]

        print(
            f"prime={p0:2d} "
            f"density={d['density']:.12f} "
            f"bad={d['bad']:2d} "
            f"workload={d['workload']:.8f}"
        )

    # ------------------------------------------------------------------------
    # Enumerate all 65536 subsets.
    # ------------------------------------------------------------------------

    print()
    print("ENUMERATING ALL SUBSETS")

    records = enumerate_all_subsets(
        n,
        prime_data,
        test_cost,
    )

    print(
        f"    total subsets = "
        f"{len(records)}"
    )

    # ------------------------------------------------------------------------
    # Pareto frontier.
    # ------------------------------------------------------------------------

    frontier = pareto_frontier(
        records
    )

    print()
    print("PARETO FRONTIER")
    print(
        f"    frontier size = "
        f"{len(frontier)}"
    )

    print()

    # Print only a manageable portion.
    for rank, record in enumerate(
        frontier[:40],
        start=1,
    ):

        print(
            f"FRONTIER {rank:2d}: "
            f"primes={list(record['subset'])} "
            f"density={record['density']:.10f} "
            f"survivors={record['survivors']:.1f} "
            f"workload={record['workload']:.6f}"
        )

    # ------------------------------------------------------------------------
    # Alternative rankings.
    # ------------------------------------------------------------------------

    density_rank = rank_by_density(
        records
    )

    workload_rank = rank_by_workload(
        records
    )

    # ------------------------------------------------------------------------
    # Combine analytically promising candidates.
    # ------------------------------------------------------------------------

    candidates = []

    def add_candidate(record):

        subset = record["subset"]

        if subset not in {
            x["subset"]
            for x in candidates
        }:

            candidates.append(record)

    # Pareto frontier.
    for record in frontier[:MAX_BENCHMARKS]:
        add_candidate(record)

    # Density-optimal.
    for record in density_rank[:10]:
        add_candidate(record)

    # Workload-optimal among non-empty subsets.
    for record in workload_rank[:10]:
        add_candidate(record)

    # Standard wheels.
    standard = [
        [3],
        [3, 5],
        [3, 5, 7],
        [3, 5, 7, 11],
        [3, 5, 7, 11, 13],
        [3, 5, 7, 11, 13, 17],
        [3, 5, 7, 11, 13, 17, 19],
        [3, 5, 7, 11, 13, 17, 19, 23],
    ]

    record_by_subset = {
        record["subset"]: record
        for record in records
    }

    for wheel in standard:

        record = record_by_subset[
            tuple(wheel)
        ]

        add_candidate(record)

    # ------------------------------------------------------------------------
    # Additional combined-score candidates.
    # ------------------------------------------------------------------------

    for weight in [
        0.25,
        0.5,
        1.0,
        2.0,
        4.0,
        8.0,
        16.0,
    ]:

        ranked = sorted(
            records,
            key=lambda r:
                combined_score(
                    r,
                    weight
                )
        )

        for record in ranked[:5]:

            add_candidate(record)

    # ------------------------------------------------------------------------
    # Sort final candidate list by analytical workload first, then density.
    # ------------------------------------------------------------------------

    candidates.sort(
        key=lambda r: (
            r["workload"],
            r["density"]
        )
    )

    print()
    print("BENCHMARK CANDIDATES")
    print(
        f"    count = "
        f"{len(candidates)}"
    )

    # ------------------------------------------------------------------------
    # Actual benchmarks.
    # ------------------------------------------------------------------------

    actual = []

    for rank, record in enumerate(
        candidates[:MAX_BENCHMARKS],
        start=1,
    ):

        subset = record["subset"]

        result = benchmark_subset(
            n,
            p,
            q,
            subset,
        )

        actual.append(
            {
                **record,
                **result,
            }
        )

        print()
        print(f"BENCHMARK {rank}")
        print(
            f"    primes = "
            f"{list(subset)}"
        )
        print(
            f"    analytical density = "
            f"{record['density']:.12f}"
        )
        print(
            f"    predicted survivors = "
            f"{record['survivors']:.2f}"
        )
        print(
            f"    workload = "
            f"{record['workload']:.8f}"
        )
        print(
            f"    actual survivors = "
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
    # Best actual result.
    # ------------------------------------------------------------------------

    successful = [
        x
        for x in actual
        if x["correct"]
    ]

    if successful:

        successful.sort(
            key=lambda x: x["time"]
        )

        best = successful[0]

        print()
        print("BEST ACTUAL WHEEL")
        print(
            f"    primes = "
            f"{list(best['subset'])}"
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
        print(
            f"    workload = "
            f"{best['workload']:.8f}"
        )

    print()
    print(f"FINISHED INSTANCE {bits}-BIT")


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("START EXPERIMENT 195")
    print()
    print(
        "QR sieve Pareto optimization"
    )
    print()

    print(
        f"PRIMES = {PRIMES}"
    )
    print(
        f"BLOCK_SIZE = {BLOCK_SIZE}"
    )
    print(
        f"MAX_Y = {MAX_Y}"
    )
    print(
        f"MAX_BENCHMARKS = {MAX_BENCHMARKS}"
    )

    # 48 and 54 are sufficient to establish the model without making the
    # benchmark excessively long.
    #
    # 60 and 66 are then included to test whether the Pareto relationship
    # remains useful at larger n.

    for bits in [
        48,
        54,
        60,
        66,
    ]:

        p, q = INSTANCES[bits]

        run_instance(
            bits,
            p,
            q,
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 195")
    print("=" * 72)


if __name__ == "__main__":
    main()
