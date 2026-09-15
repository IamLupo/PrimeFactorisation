#!/usr/bin/env python3

import math
import time
from itertools import combinations
from typing import Dict, Tuple, List

import numpy as np


# ============================================================
# START EXPERIMENT 200
#
# COST-AWARE EXHAUSTIVE QR SUBSET SEARCH
#
# We search for the subset S of QR primes minimizing:
#
#     total attack time
#
# instead of minimizing:
#
#     survivor count
#
# or:
#
#     sieve time
#
# For every subset we first calculate the EXACT number of QR
# survivors using packed bitsets.
#
# Then we estimate:
#
#     T ~= k * T_sieve_prime
#          + survivors * T_isqrt
#
# where k = number of primes.
#
# Only the best predicted subsets at each depth are subjected
# to the expensive end-to-end Fermat attack.
#
# This lets us examine tens of thousands of subsets without
# performing millions of math.isqrt() calls for each one.
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

# Exhaustively enumerate every subset up to this size.
#
# Through depth 6 with 20 primes:
#
#   C(20,1) + ... + C(20,6) = 72,999
#
MAX_DEPTH = 6

# Number of predicted candidates actually benchmarked per depth.
TOP_PREDICTED = 25

# Number of timing repetitions for end-to-end candidates.
FULL_REPEATS = 2

# Number of repetitions when measuring primitive bitset
# AND cost.
AND_REPEATS = 500

# Number of repetitions for the baseline isqrt cost estimate.
BASELINE_REPEATS = 1


# ------------------------------------------------------------
# POPCOUNT LOOKUP
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
            "p+q must be even"
        )

    return (p + q) // 2


def planted_y(p: int, q: int) -> int:
    if (q - p) & 1:
        raise ValueError(
            "q-p must be even"
        )

    return (q - p) // 2


# ------------------------------------------------------------
# QR TABLE
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

    residues = np.arange(
        p,
        dtype=np.int64,
    )

    squares[
        (residues * residues) % p
    ] = 1

    values = (
        residues * residues
        - (n % p)
    ) % p

    return squares[values]


# ------------------------------------------------------------
# PACKED MASK
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
# ALL-ONES BITSET
# ------------------------------------------------------------

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
    length: int,
    masks: Dict[int, np.ndarray],
    subset: Tuple[int, ...],
) -> Tuple[
    float,
    int,
    int,
    int,
    bool,
]:
    """
    Full Fermat attack:

        QR sieve
        -> extract survivors
        -> math.isqrt()

    Returns:

        elapsed
        survivors
        square_tests
        found_x
        correct
    """

    t0 = time.perf_counter()

    bits = make_ones(length)

    for prime in subset:
        bits &= masks[prime]

    survivors = popcount(bits)

    unpacked = np.unpackbits(
        bits,
        bitorder="little",
    )[:length]

    candidate_offsets = np.flatnonzero(
        unpacked
    )

    square_tests = 0
    found_x = -1

    for offset_np in candidate_offsets:

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

    correct = (
        found_x == true_x
    )

    return (
        elapsed,
        survivors,
        square_tests,
        found_x,
        correct,
    )


# ------------------------------------------------------------
# MEASURE PURE SIEVE COST
# ------------------------------------------------------------

def measure_prime_and_cost(
    masks: Dict[int, np.ndarray],
) -> float:
    """
    Measure average cost of one packed NumPy AND.
    """

    sample = next(
        iter(masks.values())
    )

    accumulator = sample.copy()

    best = float("inf")

    for _ in range(AND_REPEATS):

        t0 = time.perf_counter()

        accumulator &= sample

        elapsed = (
            time.perf_counter()
            - t0
        )

        if elapsed < best:
            best = elapsed

    # Prevent optimizer/elimination surprises.
    if accumulator.size == 0:
        raise RuntimeError(
            "Unexpected empty mask."
        )

    return best


# ------------------------------------------------------------
# BASELINE ISQRT COST
# ------------------------------------------------------------

def measure_isqrt_cost(
    n: int,
    start_x: int,
    true_x: int,
) -> float:
    """
    Measure average time per Fermat square test.

    We intentionally test a representative range of values,
    but do not stop at the actual factor.
    """

    count = (
        true_x
        - start_x
        + 1
    )

    if count <= 0:
        raise RuntimeError(
            "Invalid Fermat interval."
        )

    # Use at most one million samples for calibration.
    sample_count = min(
        count,
        1_000_000,
    )

    if sample_count == count:

        xs = range(
            start_x,
            true_x + 1,
        )

    else:

        step = (
            count
            // sample_count
        )

        xs = (
            start_x + i * step
            for i in range(sample_count)
        )

    total = 0
    tested = 0

    t0 = time.perf_counter()

    for x in xs:

        d = (
            x * x
            - n
        )

        root = math.isqrt(d)

        # Prevent the calculation from being optimized away.
        total ^= root

        tested += 1

    elapsed = (
        time.perf_counter()
        - t0
    )

    if total < 0:
        print("unused:", total)

    return (
        elapsed
        / tested
    )


# ------------------------------------------------------------
# EXACT SURVIVOR COUNT FOR A SUBSET
# ------------------------------------------------------------

def subset_survivors(
    subset: Tuple[int, ...],
    masks: Dict[int, np.ndarray],
    initial: np.ndarray,
) -> int:

    bits = initial.copy()

    for prime in subset:
        bits &= masks[prime]

    return popcount(bits)


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print("START EXPERIMENT 200")
    print("Cost-aware exhaustive QR subset search")
    print("=" * 78)

    for bit_size in sorted(INSTANCES):

        inst = INSTANCES[bit_size]

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
        print(f"INSTANCE {bit_size}-BIT")
        print("#" * 78)

        print("p       =", p)
        print("q       =", q)
        print("n       =", n)
        print("start_x =", start_x)
        print("true_x  =", true_x)
        print("true_y  =", true_y)
        print("gap     =", true_x - start_x)
        print("range   =", length)

        if true_x < start_x:
            raise RuntimeError(
                "Invalid Fermat interval."
            )

        # ----------------------------------------------------
        # Build masks.
        # ----------------------------------------------------

        print()
        print("Building packed QR masks...")

        masks: Dict[
            int,
            np.ndarray
        ] = {}

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
        # Verify target.
        # ----------------------------------------------------

        offset = (
            true_x
            - start_x
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
                rejected.append(prime)

        if rejected:

            raise RuntimeError(
                f"true_x rejected by {rejected}"
            )

        print(
            "True Fermat x survives every QR mask."
        )

        # ----------------------------------------------------
        # Initial mask.
        # ----------------------------------------------------

        initial = make_ones(
            length
        )

        # ----------------------------------------------------
        # Measure primitive costs.
        # ----------------------------------------------------

        print()
        print("Measuring primitive costs...")

        prime_and_cost = (
            measure_prime_and_cost(
                masks
            )
        )

        isqrt_cost = (
            measure_isqrt_cost(
                n,
                start_x,
                true_x,
            )
        )

        print(
            "  one packed AND =",
            f"{prime_and_cost:.12e}s",
        )

        print(
            "  one isqrt      =",
            f"{isqrt_cost:.12e}s",
        )

        # ----------------------------------------------------
        # Baseline full Fermat.
        # ----------------------------------------------------

        print()
        print("Measuring no-sieve baseline...")

        baseline_runs = []

        for _ in range(
            BASELINE_REPEATS
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

            baseline_runs.append(
                (
                    elapsed,
                    tests,
                    found_x,
                )
            )

        baseline = min(
            baseline_runs,
            key=lambda r: r[0],
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
            baseline[2],
        )

        # ----------------------------------------------------
        # Exhaustive subsets.
        # ----------------------------------------------------

        global_best_actual = None

        depth_summaries = []

        for depth in range(
            1,
            MAX_DEPTH + 1,
        ):

            total_subsets = math.comb(
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
                "subsets =",
                total_subsets,
            )

            predicted = []

            # ------------------------------------------------
            # Exact survivor enumeration.
            # ------------------------------------------------

            enum_t0 = time.perf_counter()

            for subset in combinations(
                PRIMES,
                depth,
            ):

                survivors = subset_survivors(
                    subset,
                    masks,
                    initial,
                )

                # Cost model.
                estimated = (
                    depth
                    * prime_and_cost
                    +
                    survivors
                    * isqrt_cost
                )

                predicted.append(
                    (
                        estimated,
                        subset,
                        survivors,
                    )
                )

            enum_elapsed = (
                time.perf_counter()
                - enum_t0
            )

            predicted.sort(
                key=lambda r: r[0]
            )

            print(
                "enumeration time =",
                f"{enum_elapsed:.3f}s",
            )

            print()
            print(
                f"TOP {min(TOP_PREDICTED, len(predicted))} "
                "PREDICTED"
            )

            for rank, (
                estimated,
                subset,
                survivors,
            ) in enumerate(
                predicted[:TOP_PREDICTED],
                1,
            ):

                print(
                    f"{rank:2d}. "
                    f"{str(list(subset)):<38} "
                    f"pred={estimated:.9f}s "
                    f"survivors={survivors:,}"
                )

            # ------------------------------------------------
            # Benchmark only top predicted subsets.
            # ------------------------------------------------

            actual_results = []

            for (
                estimated,
                subset,
                survivors,
            ) in predicted[
                :TOP_PREDICTED
            ]:

                (
                    elapsed,
                    actual_survivors,
                    square_tests,
                    found_x,
                    correct,
                ) = attack_subset(
                    n,
                    start_x,
                    true_x,
                    length,
                    masks,
                    subset,
                )

                if not correct:
                    print()
                    print(
                        "WARNING: candidate did not "
                        "find the planted factor."
                    )
                    print(
                        "subset =",
                        subset,
                    )
                    print(
                        "found_x =",
                        found_x,
                    )

                actual_results.append(
                    (
                        elapsed,
                        estimated,
                        subset,
                        actual_survivors,
                        square_tests,
                        found_x,
                        correct,
                    )
                )

                if correct:

                    if (
                        global_best_actual is None
                        or elapsed
                        < global_best_actual[0]
                    ):
                        global_best_actual = (
                            elapsed,
                            estimated,
                            subset,
                            actual_survivors,
                            square_tests,
                            found_x,
                            correct,
                        )

            actual_results.sort(
                key=lambda r: r[0]
            )

            print()
            print(
                "ACTUAL BEST AMONG PREDICTED"
            )

            for rank, result in enumerate(
                actual_results[:10],
                1,
            ):

                (
                    elapsed,
                    estimated,
                    subset,
                    survivors,
                    square_tests,
                    found_x,
                    correct,
                ) = result

                print(
                    f"{rank:2d}. "
                    f"{str(list(subset)):<38} "
                    f"actual={elapsed:.9f}s "
                    f"pred={estimated:.9f}s "
                    f"survivors={survivors:,} "
                    f"tests={square_tests:,} "
                    f"correct={correct}"
                )

            if actual_results:

                best_actual = actual_results[0]

                depth_summaries.append(
                    best_actual
                )

                (
                    elapsed,
                    estimated,
                    subset,
                    survivors,
                    square_tests,
                    found_x,
                    correct,
                ) = best_actual

                print()
                print(
                    "DEPTH BEST"
                )
                print(
                    "  subset     =",
                    list(subset),
                )
                print(
                    "  actual     =",
                    f"{elapsed:.9f}s",
                )
                print(
                    "  predicted  =",
                    f"{estimated:.9f}s",
                )
                print(
                    "  survivors  =",
                    survivors,
                )
                print(
                    "  square_test=",
                    square_tests,
                )

        # ----------------------------------------------------
        # Final instance summary.
        # ----------------------------------------------------

        print()
        print("=" * 78)
        print(
            "INSTANCE SUMMARY"
        )
        print("=" * 78)

        print(
            "baseline =",
            f"{baseline[0]:.9f}s",
        )

        for result in depth_summaries:

            (
                elapsed,
                estimated,
                subset,
                survivors,
                square_tests,
                found_x,
                correct,
            ) = result

            speedup = (
                baseline[0]
                / elapsed
            )

            print(
                f"depth={len(subset):2d} "
                f"subset={list(subset)} "
                f"actual={elapsed:.9f}s "
                f"speedup={speedup:.3f}x "
                f"survivors={survivors:,} "
                f"tests={square_tests:,}"
            )

        if global_best_actual is not None:

            (
                elapsed,
                estimated,
                subset,
                survivors,
                square_tests,
                found_x,
                correct,
            ) = global_best_actual

            speedup = (
                baseline[0]
                / elapsed
            )

            print()
            print(
                "GLOBAL BEST ACTUAL"
            )
            print(
                "  subset     =",
                list(subset),
            )
            print(
                "  runtime    =",
                f"{elapsed:.9f}s",
            )
            print(
                "  speedup    =",
                f"{speedup:.3f}x",
            )
            print(
                "  survivors  =",
                survivors,
            )
            print(
                "  square_test=",
                square_tests,
            )
            print(
                "  found_x    =",
                found_x,
            )
            print(
                "  expected_x =",
                true_x,
            )

        # ----------------------------------------------------
        # Sanity check.
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
            "FERMAT SANITY CHECK:"
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
    print("FINISHED EXPERIMENT 200")
    print("=" * 78)


if __name__ == "__main__":
    main()

