#!/usr/bin/env python3

import math
import time
from itertools import combinations
from typing import Dict, Tuple, List

import numpy as np


# ============================================================
# START EXPERIMENT 204
#
# STREAMING EXHAUSTIVE QR OPTIMIZATION
#
# This experiment avoids storing every subset's bitset.
#
# For every subset:
#
#   1. Build its packed QR mask.
#   2. Count exact survivors.
#   3. Keep only a tiny top-K list.
#   4. Discard the mask immediately.
#
# Then only those finalists are reconstructed and subjected
# to the expensive end-to-end Fermat test.
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
# SETTINGS
# ------------------------------------------------------------

MAX_DEPTH = 10

# Only this many lowest-survivor subsets are retained.
TOP_SURVIVORS = 20

# Number of retained candidates that receive the expensive
# full Fermat benchmark.
FULL_BENCHMARK_COUNT = 12

FULL_REPEATS = 2

PRINT_TOP = 15


# ------------------------------------------------------------
# POPCOUNT TABLE
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
            "p + q is odd"
        )

    return (p + q) // 2


def planted_y(p: int, q: int) -> int:
    if (q - p) & 1:
        raise ValueError(
            "q - p is odd"
        )

    return (q - p) // 2


# ------------------------------------------------------------
# QR TABLE
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

    values = (
        r * r
        - (n % prime)
    ) % prime

    return squares[values]


# ------------------------------------------------------------
# PACKED PRIME MASK
# ------------------------------------------------------------

def build_prime_mask(
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
# INITIAL ALL-ONES MASK
# ------------------------------------------------------------

def make_ones(
    length: int,
) -> np.ndarray:

    size = (
        length + 7
    ) // 8

    bits = np.full(
        size,
        0xFF,
        dtype=np.uint8,
    )

    padding = (
        size * 8
        - length
    )

    if padding:
        bits[-1] &= (
            0xFF >> padding
        )

    return bits


# ------------------------------------------------------------
# BUILD SUBSET MASK
# ------------------------------------------------------------

def make_subset_mask(
    subset: Tuple[int, ...],
    prime_masks: Dict[int, np.ndarray],
    initial: np.ndarray,
) -> np.ndarray:
    """
    Construct the packed survivor mask for one subset.

    IMPORTANT:
    This function is only called for finalists, not every
    subset in the exhaustive enumeration.
    """

    bits = initial.copy()

    for prime in subset:
        bits &= prime_masks[prime]

    return bits


# ------------------------------------------------------------
# STREAMING SUBSET ENUMERATION
# ------------------------------------------------------------

def stream_subsets(
    depth: int,
    prime_masks: Dict[int, np.ndarray],
    initial: np.ndarray,
):
    """
    Generate subsets one at a time.

    Yields:

        subset
        packed mask
        survivor count

    No full level is retained.
    """

    def recurse(
        next_index: int,
        remaining: int,
        subset: Tuple[int, ...],
        bits: np.ndarray,
    ):

        if remaining == 0:

            yield (
                subset,
                bits,
                popcount(bits),
            )

            return

        # Ensure enough primes remain to complete the subset.
        max_index = (
            len(PRIMES)
            - remaining
            + 1
        )

        for index in range(
            next_index,
            max_index,
        ):

            prime = PRIMES[index]

            child_bits = bits.copy()

            child_bits &= prime_masks[
                prime
            ]

            yield from recurse(
                index + 1,
                remaining - 1,
                subset + (prime,),
                child_bits,
            )

    yield from recurse(
        0,
        depth,
        tuple(),
        initial,
    )


# ------------------------------------------------------------
# TOP-K INSERTION
# ------------------------------------------------------------

def insert_top_k(
    top: List[
        Tuple[
            int,
            Tuple[int, ...],
        ]
    ],
    survivors: int,
    subset: Tuple[int, ...],
):
    """
    Keep only the K smallest survivor counts.

    The packed mask is deliberately NOT stored.
    """

    top.append(
        (
            survivors,
            subset,
        )
    )

    top.sort(
        key=lambda x: x[0]
    )

    if len(top) > TOP_SURVIVORS:
        del top[TOP_SURVIVORS:]


# ------------------------------------------------------------
# END-TO-END FERMAT ATTACK
# ------------------------------------------------------------

def attack_bits(
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

        offset = int(
            offset_np
        )

        x = start_x + offset

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
# FULL BENCHMARK
# ------------------------------------------------------------

def benchmark_bits(
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

        result = attack_bits(
            n,
            start_x,
            true_x,
            bits,
            length,
        )

        if (
            best is None
            or result[0]
            < best[0]
        ):
            best = result

    assert best is not None

    return best


# ------------------------------------------------------------
# BASELINE
# ------------------------------------------------------------

def benchmark_baseline(
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
    best_x = -1

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
            best_x = found_x

    return (
        best_time,
        best_tests,
        best_x,
    )


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print(
        "START EXPERIMENT 204"
    )
    print(
        "Streaming exhaustive QR optimization"
    )
    print("=" * 78)

    for bitsize in sorted(
        INSTANCES
    ):

        inst = INSTANCES[
            bitsize
        ]

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

        if true_x < start_x:
            raise RuntimeError(
                "Fermat target is before search start."
            )

        # ----------------------------------------------------
        # BUILD MASKS
        # ----------------------------------------------------

        print()
        print(
            "Building prime masks..."
        )

        prime_masks = {}

        for prime in PRIMES:

            qr = build_qr_table(
                prime,
                n,
            )

            prime_masks[
                prime
            ] = build_prime_mask(
                start_x,
                length,
                prime,
                qr,
            )

        # ----------------------------------------------------
        # VERIFY TRUE X
        # ----------------------------------------------------

        offset = (
            true_x
            - start_x
        )

        byte_index = (
            offset >> 3
        )

        bit = (
            1 << (
                offset & 7
            )
        )

        rejected = []

        for prime in PRIMES:

            if not (
                prime_masks[
                    prime
                ][byte_index]
                & bit
            ):

                rejected.append(
                    prime
                )

        if rejected:

            raise RuntimeError(
                "True x rejected by "
                f"{rejected}"
            )

        print(
            "True Fermat x survives every QR mask."
        )

        initial = make_ones(
            length
        )

        # ----------------------------------------------------
        # BASELINE
        # ----------------------------------------------------

        print()
        print(
            "Measuring baseline..."
        )

        (
            baseline_time,
            baseline_tests,
            baseline_x,
        ) = benchmark_baseline(
            n,
            start_x,
            true_x,
        )

        print(
            "  runtime =",
            f"{baseline_time:.9f}s",
        )

        print(
            "  tests   =",
            baseline_tests,
        )

        print(
            "  found_x =",
            baseline_x,
        )

        # ----------------------------------------------------
        # GLOBAL BEST
        # ----------------------------------------------------

        global_best = None

        # ----------------------------------------------------
        # DEPTH SEARCH
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

            top = []

            # ------------------------------------------------
            # Stream all subsets.
            # ------------------------------------------------

            t0 = time.perf_counter()

            count = 0

            for (
                subset,
                bits,
                survivors,
            ) in stream_subsets(
                depth,
                prime_masks,
                initial,
            ):

                count += 1

                insert_top_k(
                    top,
                    survivors,
                    subset,
                )

            enum_time = (
                time.perf_counter()
                - t0
            )

            if count != expected:

                raise RuntimeError(
                    "Subset count mismatch: "
                    f"{count} != {expected}"
                )

            print(
                "enumeration time =",
                f"{enum_time:.3f}s",
            )

            # ------------------------------------------------
            # Top survivor subsets.
            # ------------------------------------------------

            print()
            print(
                "TOP SURVIVOR COUNTS"
            )

            for rank, (
                survivors,
                subset,
            ) in enumerate(
                top[:PRINT_TOP],
                1,
            ):

                print(
                    f"{rank:2d}. "
                    f"{str(list(subset)):<40} "
                    f"survivors={survivors:,}"
                )

            # ------------------------------------------------
            # Rebuild masks only for finalists.
            # ------------------------------------------------

            print()
            print(
                "BENCHMARKING",
                min(
                    FULL_BENCHMARK_COUNT,
                    len(top),
                ),
                "FINALISTS"
            )

            actual = []

            for survivors_hint, subset in top[
                :FULL_BENCHMARK_COUNT
            ]:

                bits = make_subset_mask(
                    subset,
                    prime_masks,
                    initial,
                )

                (
                    elapsed,
                    survivors,
                    tests,
                    found_x,
                    correct,
                ) = benchmark_bits(
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

                print(
                    f"{str(list(subset)):<40} "
                    f"time={elapsed:.9f}s "
                    f"survivors={survivors:,} "
                    f"tests={tests:,} "
                    f"correct={correct}"
                )

                if not correct:

                    raise RuntimeError(
                        "A finalist failed to "
                        "recover the planted Fermat x."
                    )

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

            # ------------------------------------------------
            # Depth winner.
            # ------------------------------------------------

            if actual:

                actual.sort(
                    key=lambda x: x[0]
                )

                winner = actual[0]

                print()
                print(
                    "DEPTH BEST"
                )

                print(
                    "  subset    =",
                    list(winner[1]),
                )

                print(
                    "  runtime   =",
                    f"{winner[0]:.9f}s",
                )

                print(
                    "  speedup   =",
                    f"{baseline_time / winner[0]:.3f}x",
                )

                print(
                    "  survivors =",
                    winner[2],
                )

                print(
                    "  tests     =",
                    winner[3],
                )

        # ----------------------------------------------------
        # INSTANCE GLOBAL WINNER
        # ----------------------------------------------------

        print()
        print("=" * 78)
        print(
            f"GLOBAL BEST {bitsize}-BIT"
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
                "subset     =",
                list(subset),
            )

            print(
                "runtime    =",
                f"{elapsed:.9f}s",
            )

            print(
                "speedup    =",
                f"{baseline_time / elapsed:.3f}x",
            )

            print(
                "survivors  =",
                survivors,
            )

            print(
                "square_tests =",
                tests,
            )

            print(
                "found_x    =",
                found_x,
            )

            print(
                "expected_x =",
                true_x,
            )

            print(
                "correct    =",
                correct,
            )

        # ----------------------------------------------------
        # FERMAT SANITY CHECK
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
    print(
        "FINISHED EXPERIMENT 204"
    )
    print("=" * 78)


if __name__ == "__main__":
    main()