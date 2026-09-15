#!/usr/bin/env python3

import math
import time
from itertools import combinations
from typing import Dict, Tuple, List

import numpy as np


# ============================================================
# START EXPERIMENT 199
#
# END-TO-END FERMAT QR SIEVE
#
# Unlike Experiment 198, this experiment does NOT optimize
# sieve time alone.
#
# For every candidate subset:
#
#   1. Apply the packed QR sieve to x.
#   2. Extract surviving x values.
#   3. Test x^2 - n for being a perfect square.
#   4. Stop when the planted Fermat solution is reached.
#
# Objective:
#
#       minimize actual wall-clock attack time
#
# This measures the real workload:
#
#       sieve + candidate testing
#
# rather than only the cost of bitwise AND operations.
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
# Candidate primes
# ------------------------------------------------------------

PRIMES = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
]


# ------------------------------------------------------------
# Search parameters
# ------------------------------------------------------------

# Exhaustively test all subsets up to this depth.
#
# Number of subsets:
#
#   depth <= 1 :       20
#   depth <= 2 :      210
#   depth <= 3 :    1,350
#   depth <= 4 :    6,195
#   depth <= 5 :   21,699
#
# This is intentionally limited to keep the experiment useful.
MAX_DEPTH = 4

# Number of repetitions for timing.
REPEATS = 3

# Print best candidates after each depth.
PRINT_TOP = 15


# ------------------------------------------------------------
# Datatype
# ------------------------------------------------------------

class PrimeMask:
    __slots__ = (
        "prime",
        "allowed",
    )

    def __init__(
        self,
        prime: int,
        allowed: np.ndarray,
    ):
        self.prime = prime
        self.allowed = allowed


# ------------------------------------------------------------
# Fermat helpers
# ------------------------------------------------------------

def fermat_x_start(n: int) -> int:
    return math.isqrt(n - 1) + 1


def planted_x(p: int, q: int) -> int:
    return (p + q) // 2


def planted_y(p: int, q: int) -> int:
    return (q - p) // 2


# ------------------------------------------------------------
# QR table
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# Packed mask
# ------------------------------------------------------------

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
# Popcount
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
# Convert packed result to candidate x values
# ------------------------------------------------------------

def extract_candidates(
    bits: np.ndarray,
    start_x: int,
    length: int,
) -> np.ndarray:
    """
    Return all x values whose bit is set.

    np.unpackbits is used only AFTER the sieve, so candidate
    extraction cost is part of the actual attack workload.
    """

    values = np.unpackbits(
        bits,
        bitorder="little",
    )

    values = values[:length]

    indices = np.flatnonzero(
        values
    )

    return (
        indices.astype(np.int64)
        + start_x
    )


# ------------------------------------------------------------
# Perfect-square test
# ------------------------------------------------------------

def is_square(value: int) -> Tuple[bool, int]:
    """
    Return:
        (is_square, sqrt)
    """

    if value < 0:
        return False, 0

    root = math.isqrt(value)

    return (
        root * root == value,
        root,
    )


# ------------------------------------------------------------
# End-to-end attack
# ------------------------------------------------------------

def attack_subset(
    n: int,
    start_x: int,
    true_x: int,
    masks: Dict[int, np.ndarray],
    primes: Tuple[int, ...],
    full_length: int,
) -> Tuple[
    float,
    int,
    int,
    bool,
    int,
]:
    """
    Complete Fermat attack using a QR subset.

    Returns:

        elapsed
        QR survivor count
        square tests performed
        found factor
        found x
    """

    t0 = time.perf_counter()

    # --------------------------------------------------------
    # Start with all candidates alive.
    # --------------------------------------------------------

    size = (
        (full_length + 7) // 8
    )

    bits = np.full(
        size,
        0xFF,
        dtype=np.uint8,
    )

    extra = size * 8 - full_length

    if extra:
        bits[-1] &= (
            0xFF >> extra
        )

    # --------------------------------------------------------
    # Apply QR masks.
    # --------------------------------------------------------

    for prime in primes:
        bits &= masks[prime]

    survivor_count = popcount(bits)

    # --------------------------------------------------------
    # Extract surviving x values.
    # --------------------------------------------------------

    candidates = extract_candidates(
        bits,
        start_x,
        full_length,
    )

    # --------------------------------------------------------
    # Test candidates in ascending Fermat order.
    # --------------------------------------------------------

    square_tests = 0

    found = False
    found_x = -1

    for x_np in candidates:

        x = int(x_np)

        value = (
            x * x
            - n
        )

        square_tests += 1

        ok, root = is_square(
            value
        )

        if ok:

            found = True
            found_x = x

            # The true Fermat solution should be the first
            # square encountered.
            break

    elapsed = (
        time.perf_counter()
        - t0
    )

    return (
        elapsed,
        survivor_count,
        square_tests,
        found,
        found_x,
    )


# ------------------------------------------------------------
# Benchmark repeated subset
# ------------------------------------------------------------

def benchmark_subset(
    n: int,
    start_x: int,
    true_x: int,
    masks: Dict[int, np.ndarray],
    primes: Tuple[int, ...],
    full_length: int,
) -> Tuple[
    float,
    int,
    int,
    bool,
    int,
]:

    best = None

    for _ in range(REPEATS):

        result = attack_subset(
            n,
            start_x,
            true_x,
            masks,
            primes,
            full_length,
        )

        if (
            best is None
            or result[0] < best[0]
        ):
            best = result

    assert best is not None

    return best


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print("START EXPERIMENT 199")
    print("End-to-end Fermat QR sieve benchmark")
    print("=" * 78)

    for bits in sorted(INSTANCES):

        inst = INSTANCES[bits]

        p = inst["p"]
        q = inst["q"]
        n = inst["n"]

        start_x = fermat_x_start(n)
        true_x = planted_x(p, q)
        true_y = planted_y(p, q)

        full_length = (
            true_x
            - start_x
            + 1
        )

        print()
        print("#" * 78)
        print(f"INSTANCE {bits}-BIT")
        print("#" * 78)

        print("p       =", p)
        print("q       =", q)
        print("n       =", n)
        print("start_x =", start_x)
        print("true_x  =", true_x)
        print("true_y  =", true_y)
        print("gap     =", true_x - start_x)
        print("range   =", full_length)

        if true_x < start_x:
            raise RuntimeError(
                "Invalid Fermat interval."
            )

        # ----------------------------------------------------
        # Build masks.
        # ----------------------------------------------------

        print()
        print("Building masks...")

        masks: Dict[int, np.ndarray] = {}

        for prime in PRIMES:

            table = build_qr_table(
                prime,
                n,
            )

            masks[prime] = build_mask(
                start_x,
                full_length,
                prime,
                table,
            )

        # ----------------------------------------------------
        # Verify target survives.
        # ----------------------------------------------------

        true_offset = (
            true_x - start_x
        )

        true_byte = (
            true_offset >> 3
        )

        true_bit = (
            1 << (true_offset & 7)
        )

        rejected = []

        for prime in PRIMES:

            if not (
                masks[prime][true_byte]
                & true_bit
            ):
                rejected.append(prime)

        if rejected:

            raise RuntimeError(
                f"true_x rejected by {rejected}"
            )

        print(
            "True Fermat x survives all QR masks."
        )

        # ----------------------------------------------------
        # Baseline: no QR filtering.
        # ----------------------------------------------------

        print()
        print("BASELINE: no QR sieve")

        baseline_times = []

        for _ in range(REPEATS):

            t0 = time.perf_counter()

            found = False
            found_x = -1
            tests = 0

            for x in range(
                start_x,
                true_x + 1,
            ):

                tests += 1

                value = (
                    x * x
                    - n
                )

                root = math.isqrt(
                    value
                )

                if (
                    root * root
                    == value
                ):

                    found = True
                    found_x = x
                    break

            elapsed = (
                time.perf_counter()
                - t0
            )

            baseline_times.append(
                (
                    elapsed,
                    tests,
                    found,
                    found_x,
                )
            )

        baseline = min(
            baseline_times,
            key=lambda x: x[0]
        )

        print(
            "  runtime      =",
            f"{baseline[0]:.9f}s",
        )

        print(
            "  square tests =",
            baseline[1],
        )

        print(
            "  found_x      =",
            baseline[3],
        )

        # ----------------------------------------------------
        # Exhaustive search by subset depth.
        # ----------------------------------------------------

        global_best = None
        depth_best = {}

        for depth in range(
            1,
            MAX_DEPTH + 1,
        ):

            print()
            print("=" * 78)
            print(
                f"EXHAUSTIVE DEPTH {depth}"
            )
            print("=" * 78)

            results = []

            total = 0

            for subset in combinations(
                PRIMES,
                depth,
            ):

                total += 1

                (
                    elapsed,
                    survivors,
                    square_tests,
                    found,
                    found_x,
                ) = benchmark_subset(
                    n,
                    start_x,
                    true_x,
                    masks,
                    subset,
                    full_length,
                )

                if not found:
                    continue

                result = (
                    elapsed,
                    subset,
                    survivors,
                    square_tests,
                    found_x,
                )

                results.append(result)

                if (
                    global_best is None
                    or elapsed
                    < global_best[0]
                ):
                    global_best = result

            results.sort(
                key=lambda r: r[0]
            )

            depth_best[depth] = (
                results[0]
                if results
                else None
            )

            print(
                "tested subsets =",
                total,
            )

            print()

            if not results:

                print(
                    "No successful subsets."
                )

                continue

            print(
                f"BEST {PRINT_TOP} "
                "END-TO-END SUBSETS"
            )

            for rank, result in enumerate(
                results[:PRINT_TOP],
                1,
            ):

                (
                    elapsed,
                    subset,
                    survivors,
                    square_tests,
                    found_x,
                ) = result

                print(
                    f"{rank:2d}. "
                    f"{str(list(subset)):<35} "
                    f"time={elapsed:.9f}s "
                    f"survivors={survivors:>9,d} "
                    f"square_tests={square_tests:>9,d} "
                    f"x={found_x}"
                )

        # ----------------------------------------------------
        # Final summary.
        # ----------------------------------------------------

        print()
        print("=" * 78)
        print(
            "INSTANCE SUMMARY"
        )
        print("=" * 78)

        print(
            "No-sieve baseline:",
            f"{baseline[0]:.9f}s",
        )

        for depth in sorted(
            depth_best
        ):

            result = depth_best[depth]

            if result is None:
                continue

            (
                elapsed,
                subset,
                survivors,
                square_tests,
                found_x,
            ) = result

            speedup = (
                baseline[0]
                / elapsed
            )

            print(
                f"depth={depth:2d} "
                f"subset={list(subset)} "
                f"time={elapsed:.9f}s "
                f"speedup={speedup:.3f}x "
                f"survivors={survivors:,} "
                f"square_tests={square_tests:,}"
            )

        if global_best is not None:

            (
                elapsed,
                subset,
                survivors,
                square_tests,
                found_x,
            ) = global_best

            speedup = (
                baseline[0]
                / elapsed
            )

            print()
            print(
                "GLOBAL BEST THROUGH "
                f"DEPTH {MAX_DEPTH}"
            )

            print(
                "  primes       =",
                list(subset),
            )

            print(
                "  runtime      =",
                f"{elapsed:.9f}s",
            )

            print(
                "  speedup      =",
                f"{speedup:.3f}x",
            )

            print(
                "  QR survivors =",
                survivors,
            )

            print(
                "  square tests =",
                square_tests,
            )

            print(
                "  found_x      =",
                found_x,
            )

            print(
                "  expected_x   =",
                true_x,
            )

            print(
                "  correct      =",
                found_x == true_x,
            )

    print()
    print("=" * 78)
    print("FINISHED EXPERIMENT 199")
    print("=" * 78)


if __name__ == "__main__":
    main()
