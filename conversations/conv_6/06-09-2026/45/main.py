#!/usr/bin/env python3

import math
import time
from itertools import combinations
from typing import Dict, Tuple, List

import numpy as np


# ============================================================
# START EXPERIMENT 201
#
# INCREMENTAL EXHAUSTIVE QR SUBSET SEARCH
#
# The previous experiment showed that the analytical timing
# model predicts the right region, but not the exact winner.
#
# Therefore this experiment:
#
#   1. Enumerates subsets exhaustively.
#   2. Builds child masks incrementally from parent masks.
#   3. Uses EXACT survivor counts for ranking.
#   4. End-to-end benchmarks the top candidates.
#   5. Tests both:
#        - lowest survivor candidates
#        - diverse candidates from the ranking
#
# No regression.
# No guessed runtime model for final selection.
#
# ============================================================


# ------------------------------------------------------------
# INSTANCES
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
# SEARCH SETTINGS
# ------------------------------------------------------------

# Explore subsets through this size.
MAX_DEPTH = 7

# Number of candidates benchmarked per depth.
TOP_BENCHMARK = 100

# Additional candidates chosen from different survivor-count
# quantiles. This reduces the chance that the true runtime
# optimum is missed because it has slightly more survivors.
QUANTILE_CANDIDATES = 25

# Repeat full benchmark this many times and retain the fastest.
FULL_REPEATS = 2

# Print this many candidates.
PRINT_TOP = 20

# Popcount cache.
POPCOUNT8 = np.array(
    [bin(i).count("1") for i in range(256)],
    dtype=np.uint8,
)


# ------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------

def popcount(bits: np.ndarray) -> int:
    return int(
        POPCOUNT8[bits].sum(
            dtype=np.int64
        )
    )


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


def build_qr_table(
    prime: int,
    n: int,
) -> np.ndarray:

    p = prime

    squares = np.zeros(
        p,
        dtype=np.uint8,
    )

    r = np.arange(
        p,
        dtype=np.int64,
    )

    squares[
        (r * r) % p
    ] = 1

    residues = np.arange(
        p,
        dtype=np.int64,
    )

    values = (
        residues * residues
        - (n % p)
    ) % p

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


def make_ones(length: int) -> np.ndarray:

    size = (
        (length + 7) // 8
    )

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
# END-TO-END ATTACK
# ------------------------------------------------------------

def attack_subset(
    n: int,
    start_x: int,
    true_x: int,
    subset: Tuple[int, ...],
    masks: Dict[int, np.ndarray],
    length: int,
) -> Tuple[
    float,
    int,
    int,
    int,
    bool,
]:

    t0 = time.perf_counter()

    bits = make_ones(length)

    for prime in subset:
        bits &= masks[prime]

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

    for off_np in offsets:

        off = int(off_np)
        x = start_x + off

        d = (
            x * x
            - n
        )

        square_tests += 1

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
# BENCHMARK
# ------------------------------------------------------------

def benchmark(
    n: int,
    start_x: int,
    true_x: int,
    subset: Tuple[int, ...],
    masks: Dict[int, np.ndarray],
    length: int,
) -> Tuple[
    float,
    int,
    int,
    int,
    bool,
]:

    best = None

    for _ in range(FULL_REPEATS):

        result = attack_subset(
            n,
            start_x,
            true_x,
            subset,
            masks,
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
# INCREMENTAL SUBSET ENUMERATION
# ------------------------------------------------------------

def enumerate_depth(
    depth: int,
    prime_masks: Dict[int, np.ndarray],
    initial: np.ndarray,
) -> List[Tuple[int, Tuple[int, ...], np.ndarray]]:
    """
    Enumerate all subsets of a given depth.

    Returns:

        survivor_count,
        subset,
        resulting_bitset

    Instead of reconstructing every subset from scratch, this
    recursively combines the parent's bitset with one new prime.
    """

    results = []

    # --------------------------------------------------------
    # Depth 1.
    # --------------------------------------------------------

    if depth == 1:

        for p in PRIMES:

            bits = initial.copy()
            bits &= prime_masks[p]

            survivors = popcount(bits)

            results.append(
                (
                    survivors,
                    (p,),
                    bits,
                )
            )

        return results

    # --------------------------------------------------------
    # Recursive incremental enumeration.
    # --------------------------------------------------------

    previous = enumerate_depth(
        depth - 1,
        prime_masks,
        initial,
    )

    for survivors_parent, subset_parent, bits_parent in previous:

        last = subset_parent[-1]

        start_index = 0

        for i, p in enumerate(PRIMES):

            if p > last:
                start_index = i
                break
        else:
            continue

        for p in PRIMES[start_index:]:

            if p <= last:
                continue

            child_subset = (
                subset_parent + (p,)
            )

            child_bits = (
                bits_parent.copy()
            )

            child_bits &= prime_masks[p]

            child_survivors = popcount(
                child_bits
            )

            results.append(
                (
                    child_survivors,
                    child_subset,
                    child_bits,
                )
            )

    return results


# ------------------------------------------------------------
# BETTER ITERATIVE ENUMERATOR
# ------------------------------------------------------------

def enumerate_depth_iterative(
    depth: int,
    prime_masks: Dict[int, np.ndarray],
    initial: np.ndarray,
) -> List[Tuple[int, Tuple[int, ...], np.ndarray]]:

    current = []

    # Depth 1.
    for p in PRIMES:

        bits = initial.copy()
        bits &= prime_masks[p]

        current.append(
            (
                popcount(bits),
                (p,),
                bits,
            )
        )

    if depth == 1:
        return current

    for current_depth in range(
        2,
        depth + 1,
    ):

        next_level = []

        for parent_survivors, subset, parent_bits in current:

            last = subset[-1]

            for p in PRIMES:

                if p <= last:
                    continue

                bits = parent_bits.copy()
                bits &= prime_masks[p]

                survivors = popcount(bits)

                next_level.append(
                    (
                        survivors,
                        subset + (p,),
                        bits,
                    )
                )

        current = next_level

    return current


# ------------------------------------------------------------
# SELECT DIVERSE CANDIDATES
# ------------------------------------------------------------

def choose_candidates(
    results: List[
        Tuple[int, Tuple[int, ...], np.ndarray]
    ],
) -> List[
    Tuple[int, Tuple[int, ...]]
]:

    if not results:
        return []

    # Sort by exact survivor count.
    ordered = sorted(
        results,
        key=lambda x: x[0],
    )

    selected = {}

    # --------------------------------------------------------
    # Best TOP_BENCHMARK.
    # --------------------------------------------------------

    for entry in ordered[
        :TOP_BENCHMARK
    ]:

        survivors, subset, bits = entry

        selected[subset] = (
            survivors,
            subset,
        )

    # --------------------------------------------------------
    # Quantile sampling.
    # --------------------------------------------------------

    count = len(ordered)

    if count > 1:

        for i in range(
            QUANTILE_CANDIDATES
        ):

            fraction = (
                i
                / (
                    QUANTILE_CANDIDATES - 1
                )
            )

            index = int(
                fraction * (count - 1)
            )

            survivors, subset, bits = (
                ordered[index]
            )

            selected[subset] = (
                survivors,
                subset,
            )

    return list(
        selected.values()
    )


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print("START EXPERIMENT 201")
    print("Incremental exhaustive QR subset search")
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
            "Building packed QR masks..."
        )

        masks = {}

        for prime in PRIMES:

            qr_table = build_qr_table(
                prime,
                n,
            )

            masks[prime] = build_mask(
                start_x,
                length,
                prime,
                qr_table,
            )

        # ----------------------------------------------------
        # Verify target.
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
                f"true_x rejected by {rejected}"
            )

        print(
            "True Fermat x survives every QR test."
        )

        initial = make_ones(
            length
        )

        global_best = None

        # ----------------------------------------------------
        # Depth loop.
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

            t0 = time.perf_counter()

            results = (
                enumerate_depth_iterative(
                    depth,
                    masks,
                    initial,
                )
            )

            enumeration_time = (
                time.perf_counter()
                - t0
            )

            if len(results) != expected:

                raise RuntimeError(
                    "Enumeration count mismatch: "
                    f"got {len(results)}, "
                    f"expected {expected}"
                )

            # ------------------------------------------------
            # Sort by survivor count.
            # ------------------------------------------------

            results.sort(
                key=lambda x: x[0]
            )

            print(
                "enumeration time =",
                f"{enumeration_time:.3f}s",
            )

            print()
            print(
                f"TOP {PRINT_TOP} BY SURVIVOR COUNT"
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
            # Diverse empirical benchmark set.
            # ------------------------------------------------

            candidates = choose_candidates(
                results
            )

            print()
            print(
                "Benchmarking",
                len(candidates),
                "candidates..."
            )

            actual = []

            benchmark_t0 = (
                time.perf_counter()
            )

            for survivors_hint, subset in candidates:

                (
                    elapsed,
                    survivors,
                    square_tests,
                    found_x,
                    correct,
                ) = benchmark(
                    n,
                    start_x,
                    true_x,
                    subset,
                    masks,
                    length,
                )

                actual.append(
                    (
                        elapsed,
                        subset,
                        survivors,
                        square_tests,
                        found_x,
                        correct,
                    )
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

            print()
            print(
                "BEST END-TO-END RESULTS"
            )

            for rank, result in enumerate(
                actual[:PRINT_TOP],
                1,
            ):

                (
                    elapsed,
                    subset,
                    survivors,
                    square_tests,
                    found_x,
                    correct,
                ) = result

                print(
                    f"{rank:2d}. "
                    f"{str(list(subset)):<38} "
                    f"time={elapsed:.9f}s "
                    f"survivors={survivors:,} "
                    f"tests={square_tests:,} "
                    f"correct={correct}"
                )

            if actual:

                best_depth = actual[0]

                if (
                    global_best is None
                    or best_depth[0]
                    < global_best[0]
                ):
                    global_best = best_depth

                print()
                print(
                    "DEPTH BEST"
                )
                print(
                    "  subset     =",
                    list(best_depth[1]),
                )
                print(
                    "  runtime    =",
                    f"{best_depth[0]:.9f}s",
                )
                print(
                    "  survivors  =",
                    best_depth[2],
                )
                print(
                    "  square_tests =",
                    best_depth[3],
                )

        # ----------------------------------------------------
        # Final summary.
        # ----------------------------------------------------

        print()
        print("=" * 78)
        print(
            f"FINAL SUMMARY {bitsize}-BIT"
        )
        print("=" * 78)

        if global_best is not None:

            (
                elapsed,
                subset,
                survivors,
                square_tests,
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
                "SURVIVORS   =",
                survivors,
            )
            print(
                "SQUARE TESTS=",
                square_tests,
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
        # Fermat identity.
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
    print("FINISHED EXPERIMENT 201")
    print("=" * 78)


if __name__ == "__main__":
    main()
