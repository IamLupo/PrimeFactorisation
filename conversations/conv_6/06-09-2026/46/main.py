#!/usr/bin/env python3

import math
import time
from itertools import combinations
from typing import Dict, Tuple, List

import numpy as np


# ============================================================
# START EXPERIMENT 202
#
# EXACT SURVIVOR SEARCH + SELECTIVE EMPIRICAL VALIDATION
#
# We exhaustively enumerate QR prime subsets, but DO NOT run
# the expensive Fermat square test for every subset.
#
# For each subset we calculate the exact packed survivor mask.
#
# We then empirically benchmark only:
#
#   1. lowest-survivor subsets
#   2. distributionally diverse subsets
#   3. best subset from every previous depth
#
# The objective remains:
#
#       actual end-to-end Fermat runtime
#
# ============================================================


# ------------------------------------------------------------
# TEST INSTANCES
# ------------------------------------------------------------

INSTANCES = {
    48: {
        "p": 8_390_069,
        "q": 33_547_589,
        "n": 281_466_586_493_641,
    },

    54: {
        "p": 124_517_461,
        "q": 144_517_463,
        "n": 17_994_947_562_921_443,
    },

    60: {
        "p": 1_058_841_403,
        "q": 1_088_841_421,
        "n": 1_152_910_377_856_153_663,
    },

    66: {
        "p": 8_569_934_017,
        "q": 8_609_934_041,
        "n": 73_786_566_622_092_172_697,
    },
}


# ------------------------------------------------------------
# PRIME POOL
# ------------------------------------------------------------

PRIMES = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
]


# ------------------------------------------------------------
# SEARCH PARAMETERS
# ------------------------------------------------------------

MAX_DEPTH = 8

# How many lowest-survivor candidates to benchmark.
TOP_BY_SURVIVORS = 12

# Additional candidates sampled across the survivor ranking.
DISTRIBUTION_SAMPLES = 12

# End-to-end timing repetitions.
FULL_REPEATS = 2

# How many historical depth winners to retain.
KEEP_DEPTH_WINNERS = 2

# Print counts.
PRINT_TOP = 15


# ------------------------------------------------------------
# POPCOUNT
# ------------------------------------------------------------

POPCOUNT8 = np.array(
    [bin(i).count("1") for i in range(256)],
    dtype=np.uint8,
)


def popcount(bits: np.ndarray) -> int:
    return int(
        POPCOUNT8[bits].sum(
            dtype=np.int64
        )
    )


# ------------------------------------------------------------
# FERMAT HELPERS
# ------------------------------------------------------------

def fermat_x_start(n: int) -> int:
    return math.isqrt(n - 1) + 1


def planted_x(p: int, q: int) -> int:
    if (p + q) & 1:
        raise ValueError(
            "p + q must be even"
        )

    return (p + q) // 2


def planted_y(p: int, q: int) -> int:
    if (q - p) & 1:
        raise ValueError(
            "q - p must be even"
        )

    return (q - p) // 2


# ------------------------------------------------------------
# QR MASK
# ------------------------------------------------------------

def build_qr_table(
    prime: int,
    n: int,
) -> np.ndarray:

    squares = np.zeros(
        prime,
        dtype=np.uint8,
    )

    r = np.arange(
        prime,
        dtype=np.int64,
    )

    squares[
        (r * r) % prime
    ] = 1

    residues = np.arange(
        prime,
        dtype=np.int64,
    )

    values = (
        residues * residues
        - (n % prime)
    ) % prime

    return squares[values]


def build_mask(
    start_x: int,
    length: int,
    prime: int,
    qr_table: np.ndarray,
) -> np.ndarray:

    x = np.arange(
        start_x,
        start_x + length,
        dtype=np.int64,
    )

    allowed = qr_table[
        x % prime
    ]

    return np.packbits(
        allowed,
        bitorder="little",
    )


# ------------------------------------------------------------
# ONES MASK
# ------------------------------------------------------------

def make_ones(length: int) -> np.ndarray:

    size = (
        length + 7
    ) // 8

    bits = np.full(
        size,
        0xFF,
        dtype=np.uint8,
    )

    padding = (
        size * 8 - length
    )

    if padding:
        bits[-1] &= (
            0xFF >> padding
        )

    return bits


# ------------------------------------------------------------
# EXACT SUBSET ENUMERATION
# ------------------------------------------------------------

def enumerate_depth(
    depth: int,
    masks: Dict[int, np.ndarray],
    initial: np.ndarray,
) -> List[
    Tuple[
        int,
        Tuple[int, ...],
        np.ndarray,
    ]
]:

    # --------------------------------------------------------
    # Depth 1
    # --------------------------------------------------------

    current = []

    for prime in PRIMES:

        bits = (
            initial.copy()
        )

        bits &= masks[prime]

        survivors = popcount(bits)

        current.append(
            (
                survivors,
                (prime,),
                bits,
            )
        )

    if depth == 1:
        return current

    # --------------------------------------------------------
    # Incrementally extend.
    # --------------------------------------------------------

    for current_depth in range(
        2,
        depth + 1,
    ):

        next_level = []

        for (
            parent_survivors,
            subset,
            parent_bits,
        ) in current:

            last = subset[-1]

            for prime in PRIMES:

                if prime <= last:
                    continue

                child_bits = (
                    parent_bits.copy()
                )

                child_bits &= masks[prime]

                survivors = popcount(
                    child_bits
                )

                child_subset = (
                    subset + (prime,)
                )

                next_level.append(
                    (
                        survivors,
                        child_subset,
                        child_bits,
                    )
                )

        current = next_level

    return current


# ------------------------------------------------------------
# END-TO-END ATTACK USING AN ALREADY COMPUTED MASK
# ------------------------------------------------------------

def attack_mask(
    n: int,
    start_x: int,
    true_x: int,
    bits: np.ndarray,
    length: int,
) -> Tuple[
    float,
    int,
    int,
    int,
    bool,
]:

    t0 = time.perf_counter()

    survivors = popcount(bits)

    unpacked = np.unpackbits(
        bits,
        bitorder="little",
    )[:length]

    offsets = np.flatnonzero(
        unpacked
    )

    square_tests = 0
    found_x = -1

    for offset_np in offsets:

        offset = int(offset_np)

        x = (
            start_x
            + offset
        )

        square_tests += 1

        d = (
            x * x
            - n
        )

        root = math.isqrt(d)

        if root * root == d:

            found_x = x
            break

    elapsed = (
        time.perf_counter()
        - t0
    )

    return (
        elapsed,
        survivors,
        square_tests,
        found_x,
        found_x == true_x,
    )


# ------------------------------------------------------------
# REPEATED BENCHMARK
# ------------------------------------------------------------

def benchmark_mask(
    n: int,
    start_x: int,
    true_x: int,
    bits: np.ndarray,
    length: int,
) -> Tuple[
    float,
    int,
    int,
    int,
    bool,
]:

    best = None

    for _ in range(
        FULL_REPEATS
    ):

        result = attack_mask(
            n,
            start_x,
            true_x,
            bits,
            length,
        )

        if (
            best is None
            or result[0] < best[0]
        ):
            best = result

    assert best is not None

    return best


# ------------------------------------------------------------
# SELECT VALIDATION SET
# ------------------------------------------------------------

def choose_validation_set(
    results: List[
        Tuple[
            int,
            Tuple[int, ...],
            np.ndarray,
        ]
    ],
) -> List[
    Tuple[
        int,
        Tuple[int, ...],
        np.ndarray,
    ]
]:

    if not results:
        return []

    ordered = sorted(
        results,
        key=lambda x: x[0],
    )

    selected = {}

    # --------------------------------------------------------
    # Best survivor-count candidates.
    # --------------------------------------------------------

    for item in ordered[
        :TOP_BY_SURVIVORS
    ]:

        survivors, subset, bits = item

        selected[subset] = item

    # --------------------------------------------------------
    # Sample the rest of the survivor distribution.
    #
    # This protects against the assumption:
    #
    #     fewer survivors == faster runtime
    #
    # --------------------------------------------------------

    count = len(ordered)

    if count > 1:

        for i in range(
            DISTRIBUTION_SAMPLES
        ):

            fraction = (
                i
                / (
                    DISTRIBUTION_SAMPLES - 1
                )
            )

            index = int(
                fraction * (count - 1)
            )

            item = ordered[index]

            selected[item[1]] = item

    return list(
        selected.values()
    )


# ------------------------------------------------------------
# BASELINE
# ------------------------------------------------------------

def measure_baseline(
    n: int,
    start_x: int,
    true_x: int,
) -> Tuple[
    float,
    int,
    int,
]:

    best_time = float("inf")
    best_tests = 0
    best_found = -1

    for _ in range(
        FULL_REPEATS
    ):

        t0 = time.perf_counter()

        tests = 0
        found_x = -1

        for x in range(
            start_x,
            true_x + 1,
        ):

            tests += 1

            d = (
                x * x
                - n
            )

            root = math.isqrt(d)

            if root * root == d:

                found_x = x
                break

        elapsed = (
            time.perf_counter()
            - t0
        )

        if elapsed < best_time:

            best_time = elapsed
            best_tests = tests
            best_found = found_x

    return (
        best_time,
        best_tests,
        best_found,
    )


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print("START EXPERIMENT 202")
    print(
        "Exact survivor search + selective empirical validation"
    )
    print("=" * 78)

    for bitsize in sorted(
        INSTANCES
    ):

        inst = INSTANCES[bitsize]

        p = inst["p"]
        q = inst["q"]
        n = inst["n"]

        start_x = fermat_x_start(n)
        true_x = planted_x(p, q)
        true_y = planted_y(p, q)

        length = (
            true_x
            - start_x
            + 1
        )

        print()
        print("#" * 78)
        print(
            f"INSTANCE {bitsize}-BIT"
        )
        print("#" * 78)

        print(
            "p       =",
            p,
        )

        print(
            "q       =",
            q,
        )

        print(
            "n       =",
            n,
        )

        print(
            "start_x =",
            start_x,
        )

        print(
            "true_x  =",
            true_x,
        )

        print(
            "true_y  =",
            true_y,
        )

        print(
            "gap     =",
            true_x - start_x,
        )

        print(
            "range   =",
            length,
        )

        # ----------------------------------------------------
        # Build masks.
        # ----------------------------------------------------

        print()
        print(
            "Building QR masks..."
        )

        masks = {}

        for prime in PRIMES:

            qr = build_qr_table(
                prime,
                n,
            )

            masks[prime] = build_mask(
                start_x,
                length,
                prime,
                qr,
            )

        # ----------------------------------------------------
        # Verify planted x.
        # ----------------------------------------------------

        offset = (
            true_x - start_x
        )

        byte_index = (
            offset >> 3
        )

        bit = (
            1 << (offset & 7)
        )

        rejected = []

        for prime in PRIMES:

            if not (
                masks[prime][byte_index]
                & bit
            ):

                rejected.append(
                    prime
                )

        if rejected:

            raise RuntimeError(
                "true_x rejected by "
                f"{rejected}"
            )

        print(
            "True Fermat x survives every QR test."
        )

        initial = make_ones(
            length
        )

        # ----------------------------------------------------
        # Baseline.
        # ----------------------------------------------------

        print()
        print(
            "Measuring no-sieve baseline..."
        )

        (
            baseline_time,
            baseline_tests,
            baseline_found,
        ) = measure_baseline(
            n,
            start_x,
            true_x,
        )

        print(
            f"  runtime      = "
            f"{baseline_time:.9f}s"
        )

        print(
            f"  square tests = "
            f"{baseline_tests:,}"
        )

        print(
            f"  found_x      = "
            f"{baseline_found}"
        )

        # ----------------------------------------------------
        # Keep global best.
        # ----------------------------------------------------

        global_best = None

        historical_winners = []

        # ----------------------------------------------------
        # Search depths.
        # ----------------------------------------------------

        for depth in range(
            1,
            MAX_DEPTH + 1,
        ):

            expected = math.comb(
                len(PRIMES),
                depth,
            )

            print()
            print("=" * 78)
            print(
                f"DEPTH {depth}"
            )
            print("=" * 78)

            print(
                "expected subsets =",
                expected,
            )

            # ------------------------------------------------
            # Exact enumeration.
            # ------------------------------------------------

            enum_t0 = (
                time.perf_counter()
            )

            results = (
                enumerate_depth(
                    depth,
                    masks,
                    initial,
                )
            )

            enum_elapsed = (
                time.perf_counter()
                - enum_t0
            )

            if len(results) != expected:

                raise RuntimeError(
                    "Enumeration count mismatch: "
                    f"{len(results)} != {expected}"
                )

            results.sort(
                key=lambda x: x[0]
            )

            print(
                "enumeration time =",
                f"{enum_elapsed:.3f}s",
            )

            # ------------------------------------------------
            # Print best survivor subsets.
            # ------------------------------------------------

            print()
            print(
                "TOP BY SURVIVOR COUNT"
            )

            for rank, (
                survivors,
                subset,
                _,
            ) in enumerate(
                results[:PRINT_TOP],
                1,
            ):

                print(
                    f"{rank:2d}. "
                    f"{str(list(subset)):<38} "
                    f"survivors={survivors:,}"
                )

            # ------------------------------------------------
            # Select validation set.
            # ------------------------------------------------

            validation = (
                choose_validation_set(
                    results
                )
            )

            # Add previous depth winners so they remain in
            # consideration even if their survivor count is not
            # exceptional at a later stage.
            for old in historical_winners:

                validation_subset = old[1]

                for item in results:

                    if (
                        item[1]
                        == validation_subset
                    ):

                        validation.append(
                            item
                        )

                        break

            # Deduplicate.
            validation_map = {}

            for item in validation:
                validation_map[
                    item[1]
                ] = item

            validation = list(
                validation_map.values()
            )

            print()
            print(
                "EMPIRICAL VALIDATION COUNT =",
                len(validation),
            )

            # ------------------------------------------------
            # Benchmark.
            # ------------------------------------------------

            actual = []

            benchmark_t0 = (
                time.perf_counter()
            )

            for (
                survivors_hint,
                subset,
                bits,
            ) in validation:

                (
                    elapsed,
                    survivors,
                    tests,
                    found_x,
                    correct,
                ) = benchmark_mask(
                    n,
                    start_x,
                    true_x,
                    bits,
                    length,
                )

                actual.append(
                    (
                        elapsed,
                        subset,
                        survivors,
                        tests,
                        found_x,
                        correct,
                    )
                )

                if correct:

                    if (
                        global_best is None
                        or elapsed
                        < global_best[0]
                    ):

                        global_best = (
                            elapsed,
                            subset,
                            survivors,
                            tests,
                            found_x,
                            correct,
                        )

            benchmark_elapsed = (
                time.perf_counter()
                - benchmark_t0
            )

            actual.sort(
                key=lambda x: x[0]
            )

            print(
                "benchmark time =",
                f"{benchmark_elapsed:.3f}s",
            )

            # ------------------------------------------------
            # Results.
            # ------------------------------------------------

            print()
            print(
                "BEST EMPIRICAL RESULTS"
            )

            for rank, result in enumerate(
                actual[:PRINT_TOP],
                1,
            ):

                (
                    elapsed,
                    subset,
                    survivors,
                    tests,
                    found_x,
                    correct,
                ) = result

                speedup = (
                    baseline_time
                    / elapsed
                )

                print(
                    f"{rank:2d}. "
                    f"{str(list(subset)):<38} "
                    f"time={elapsed:.9f}s "
                    f"speedup={speedup:.3f}x "
                    f"survivors={survivors:,} "
                    f"tests={tests:,} "
                    f"correct={correct}"
                )

            # ------------------------------------------------
            # Keep depth winner.
            # ------------------------------------------------

            if actual:

                best_depth = actual[0]

                historical_winners.append(
                    best_depth
                )

                if len(
                    historical_winners
                ) > KEEP_DEPTH_WINNERS * MAX_DEPTH:

                    historical_winners = (
                        historical_winners[
                            -KEEP_DEPTH_WINNERS * MAX_DEPTH:
                        ]
                    )

                print()
                print(
                    "DEPTH BEST"
                )

                print(
                    "  subset    =",
                    list(best_depth[1]),
                )

                print(
                    "  runtime   =",
                    f"{best_depth[0]:.9f}s",
                )

                print(
                    "  survivors =",
                    best_depth[2],
                )

                print(
                    "  tests     =",
                    best_depth[3],
                )

        # ----------------------------------------------------
        # Global best.
        # ----------------------------------------------------

        print()
        print("=" * 78)
        print(
            f"GLOBAL RESULT {bitsize}-BIT"
        )
        print("=" * 78)

        if global_best is not None:

            (
                elapsed,
                subset,
                survivors,
                tests,
                found_x,
                correct,
            ) = global_best

            print(
                "BEST SUBSET =",
                list(subset),
            )

            print(
                "RUNTIME     =",
                f"{elapsed:.9f}s",
            )

            print(
                "SPEEDUP     =",
                f"{baseline_time / elapsed:.3f}x",
            )

            print(
                "SURVIVORS   =",
                survivors,
            )

            print(
                "SQUARE TESTS=",
                tests,
            )

            print(
                "FOUND X     =",
                found_x,
            )

            print(
                "EXPECTED X  =",
                true_x,
            )

            print(
                "CORRECT     =",
                correct,
            )

        # ----------------------------------------------------
        # Fermat sanity.
        # ----------------------------------------------------

        lhs = (
            true_x * true_x
            - n
        )

        rhs = (
            true_y * true_y
        )

        print()
        print(
            "FERMAT SANITY CHECK"
        )

        print(
            "  x^2 - n =",
            lhs,
        )

        print(
            "  y^2     =",
            rhs,
        )

        print(
            "  equal   =",
            lhs == rhs,
        )

    print()
    print("=" * 78)
    print("FINISHED EXPERIMENT 202")
    print("=" * 78)


if __name__ == "__main__":
    main()
