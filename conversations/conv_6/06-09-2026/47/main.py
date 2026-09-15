#!/usr/bin/env python3

import math
import time
import random
from itertools import combinations
from typing import Dict, Tuple, List

import numpy as np


# ============================================================
# START EXPERIMENT 203
#
# EMPIRICAL COST CALIBRATION + QR SUBSET OPTIMIZATION
#
# Goal:
#
#   Find the practical QR-prime subset that minimizes the
#   complete Fermat factorization runtime.
#
# Instead of fully testing every subset, we:
#
#   1. Calculate exact QR survivor counts for every subset.
#   2. Empirically measure square-test cost using samples.
#   3. Build a simple measured cost estimate.
#   4. Benchmark only a small elite set end-to-end.
#
# This avoids the huge benchmarking cost seen in Exp 201/202.
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

MAX_DEPTH = 8

# Number of candidates fully benchmarked after ranking.
FULL_BENCHMARKS = 8

# Number of subsets selected from the best survivor counts.
SURVIVOR_ELITE = 5

# Number selected from the predicted-runtime ranking.
MODEL_ELITE = 5

# Number of quantile samples across the entire candidate space.
QUANTILE_SAMPLES = 5

# Number of survivor x-values actually tested when calibrating
# the expensive math.isqrt component.
SQUARE_SAMPLE = 20_000

# Full benchmark repetitions.
FULL_REPEATS = 2

# Random seed for deterministic experiment.
RANDOM_SEED = 203


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
# FERMAT
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
# ALL-ONES
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
# APPLY SUBSET
# ------------------------------------------------------------

def make_subset_mask(
    subset: Tuple[int, ...],
    prime_masks: Dict[int, np.ndarray],
    initial: np.ndarray,
) -> np.ndarray:

    result = initial.copy()

    for prime in subset:
        result &= prime_masks[prime]

    return result


# ------------------------------------------------------------
# EXTRACT SURVIVOR OFFSETS
# ------------------------------------------------------------

def extract_offsets(
    bits: np.ndarray,
    length: int,
) -> np.ndarray:

    values = np.unpackbits(
        bits,
        bitorder="little",
    )[:length]

    return np.flatnonzero(
        values
    )


# ------------------------------------------------------------
# EMPIRICAL SQUARE-TEST CALIBRATION
# ------------------------------------------------------------

def calibrate_isqrt(
    n: int,
    start_x: int,
    true_x: int,
    sample_count: int,
) -> float:
    """
    Measure the actual Python cost of the Fermat square test.

    We deliberately do not test the candidates in their actual
    order here. The purpose is to measure the primitive cost.
    """

    length = (
        true_x
        - start_x
        + 1
    )

    sample_count = min(
        sample_count,
        length,
    )

    rng = random.Random(
        RANDOM_SEED
    )

    offsets = [
        rng.randrange(length)
        for _ in range(sample_count)
    ]

    # Sort so memory behavior is relatively stable.
    offsets.sort()

    checksum = 0

    t0 = time.perf_counter()

    for offset in offsets:

        x = (
            start_x
            + offset
        )

        d = (
            x * x
            - n
        )

        root = math.isqrt(d)

        checksum ^= root

    elapsed = (
        time.perf_counter()
        - t0
    )

    # Keep the computation observable.
    if checksum == -1:
        print("checksum:", checksum)

    return (
        elapsed
        / sample_count
    )


# ------------------------------------------------------------
# CALIBRATE PACKED AND COST
# ------------------------------------------------------------

def calibrate_and_cost(
    prime_masks: Dict[int, np.ndarray],
) -> float:

    first = next(
        iter(prime_masks.values())
    )

    accumulator = first.copy()

    # Use a reasonably large number of operations.
    rounds = 2000

    t0 = time.perf_counter()

    for _ in range(rounds):
        accumulator &= first

    elapsed = (
        time.perf_counter()
        - t0
    )

    if accumulator.size == 0:
        raise RuntimeError(
            "Unexpected empty packed mask."
        )

    return (
        elapsed
        / rounds
    )


# ------------------------------------------------------------
# EXACT SUBSET ENUMERATION
# ------------------------------------------------------------

def enumerate_subsets_of_depth(
    depth: int,
    prime_masks: Dict[int, np.ndarray],
    initial: np.ndarray,
) -> List[
    Tuple[
        int,
        Tuple[int, ...],
        np.ndarray,
    ]
]:
    """
    Exact exhaustive enumeration.

    Every child is formed from the parent's packed mask.
    """

    current = []

    for prime in PRIMES:

        bits = initial.copy()
        bits &= prime_masks[prime]

        current.append(
            (
                popcount(bits),
                (prime,),
                bits,
            )
        )

    if depth == 1:
        return current

    for _ in range(
        2,
        depth + 1,
    ):

        nxt = []

        for (
            parent_survivors,
            subset,
            parent_bits,
        ) in current:

            last = subset[-1]

            for prime in PRIMES:

                if prime <= last:
                    continue

                bits = parent_bits.copy()

                bits &= prime_masks[
                    prime
                ]

                survivors = popcount(
                    bits
                )

                nxt.append(
                    (
                        survivors,
                        subset + (prime,),
                        bits,
                    )
                )

        current = nxt

    return current


# ------------------------------------------------------------
# EMPIRICAL EXTRACTION COST
# ------------------------------------------------------------

def measure_extraction_cost(
    bits: np.ndarray,
    length: int,
    repeats: int = 3,
) -> float:
    """
    Measure unpackbits + flatnonzero.
    """

    best = float("inf")

    for _ in range(repeats):

        t0 = time.perf_counter()

        values = np.unpackbits(
            bits,
            bitorder="little",
        )[:length]

        offsets = np.flatnonzero(
            values
        )

        elapsed = (
            time.perf_counter()
            - t0
        )

        if offsets.size < 0:
            raise RuntimeError(
                "Impossible offset count."
            )

        best = min(
            best,
            elapsed,
        )

    return best


# ------------------------------------------------------------
# FULL ATTACK
# ------------------------------------------------------------

def full_attack(
    n: int,
    start_x: int,
    true_x: int,
    subset: Tuple[int, ...],
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

    offsets = extract_offsets(
        bits,
        length,
    )

    square_tests = 0
    found_x = -1

    for offset_np in offsets:

        offset = int(
            offset_np
        )

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
# SELECT EMPIRICAL CANDIDATES
# ------------------------------------------------------------

def choose_candidates(
    results,
    predicted,
) -> Dict[
    Tuple[int, ...],
    Tuple[int, Tuple[int, ...], np.ndarray],
]:

    selected = {}

    # Best survivor counts.
    survivor_sorted = sorted(
        results,
        key=lambda x: x[0],
    )

    for item in survivor_sorted[
        :SURVIVOR_ELITE
    ]:

        selected[item[1]] = item

    # Best predicted.
    predicted_sorted = sorted(
        predicted,
        key=lambda x: x[0],
    )

    for item in predicted_sorted[
        :MODEL_ELITE
    ]:

        subset = item[1]

        for result in results:

            if result[1] == subset:
                selected[subset] = result
                break

    # Quantile sampling by survivor count.
    count = len(
        survivor_sorted
    )

    if count > 1:

        for i in range(
            QUANTILE_SAMPLES
        ):

            fraction = (
                i
                / (
                    QUANTILE_SAMPLES - 1
                )
            )

            idx = int(
                fraction
                * (count - 1)
            )

            item = survivor_sorted[
                idx
            ]

            selected[
                item[1]
            ] = item

    return selected


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print("START EXPERIMENT 203")
    print(
        "Empirical cost calibration + QR subset optimization"
    )
    print("=" * 78)

    for bit_size in sorted(
        INSTANCES
    ):

        inst = INSTANCES[
            bit_size
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
            f"INSTANCE {bit_size}-BIT"
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

        prime_masks = {}

        for prime in PRIMES:

            qr = build_qr_table(
                prime,
                n,
            )

            prime_masks[
                prime
            ] = build_mask(
                start_x,
                length,
                prime,
                qr,
            )

        # ----------------------------------------------------
        # Target verification.
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
                prime_masks[prime][
                    byte_index
                ]
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
            "True Fermat x survives every QR mask."
        )

        initial = make_ones(
            length
        )

        # ----------------------------------------------------
        # Primitive cost calibration.
        # ----------------------------------------------------

        print()
        print(
            "Calibrating primitive costs..."
        )

        and_cost = calibrate_and_cost(
            prime_masks
        )

        isqrt_cost = calibrate_isqrt(
            n,
            start_x,
            true_x,
            SQUARE_SAMPLE,
        )

        print(
            "  packed AND =",
            f"{and_cost:.12e}s",
        )

        print(
            "  isqrt      =",
            f"{isqrt_cost:.12e}s",
        )

        # ----------------------------------------------------
        # Baseline.
        # ----------------------------------------------------

        baseline_t0 = time.perf_counter()

        baseline_tests = 0

        for x in range(
            start_x,
            true_x + 1,
        ):

            baseline_tests += 1

            d = (
                x * x
                - n
            )

            root = math.isqrt(
                d
            )

            if root * root == d:
                break

        baseline_time = (
            time.perf_counter()
            - baseline_t0
        )

        print()
        print(
            "BASELINE"
        )

        print(
            "  runtime =",
            f"{baseline_time:.9f}s",
        )

        print(
            "  tests   =",
            baseline_tests,
        )

        # ----------------------------------------------------
        # Global winner.
        # ----------------------------------------------------

        global_best = None

        # ----------------------------------------------------
        # Search each depth.
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
            # Enumerate exact survivor masks.
            # ------------------------------------------------

            enum_t0 = time.perf_counter()

            results = (
                enumerate_subsets_of_depth(
                    depth,
                    prime_masks,
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

            print(
                "enumeration time =",
                f"{enum_elapsed:.3f}s",
            )

            # ------------------------------------------------
            # Measure extraction and use empirical timing.
            #
            # Only a small sample of masks gets directly timed.
            # ------------------------------------------------

            extraction_samples = []

            # Lowest survivors.
            survivor_sorted = sorted(
                results,
                key=lambda x: x[0],
            )

            extraction_samples.extend(
                survivor_sorted[
                    :min(
                        5,
                        len(survivor_sorted),
                    )
                ]
            )

            # Some larger survivor counts.
            count = len(
                survivor_sorted
            )

            for i in range(5):

                if count <= 1:
                    break

                fraction = (
                    i / 4
                )

                idx = int(
                    fraction
                    * (count - 1)
                )

                extraction_samples.append(
                    survivor_sorted[idx]
                )

            seen = set()
            extraction_samples_unique = []

            for item in extraction_samples:

                subset = item[1]

                if subset in seen:
                    continue

                seen.add(subset)

                extraction_samples_unique.append(
                    item
                )

            extraction_per_byte = []

            for (
                survivors,
                subset,
                bits,
            ) in extraction_samples_unique:

                extraction_time = (
                    measure_extraction_cost(
                        bits,
                        length,
                    )
                )

                extraction_per_byte.append(
                    extraction_time
                )

            if extraction_per_byte:

                extraction_cost = (
                    sum(
                        extraction_per_byte
                    )
                    / len(
                        extraction_per_byte
                    )
                )

            else:

                extraction_cost = 0.0

            # ------------------------------------------------
            # Predicted total runtime.
            # ------------------------------------------------

            predicted = []

            for (
                survivors,
                subset,
                bits,
            ) in results:

                sieve_cost = (
                    len(subset)
                    * and_cost
                )

                square_cost = (
                    survivors
                    * isqrt_cost
                )

                estimate = (
                    sieve_cost
                    + extraction_cost
                    + square_cost
                )

                predicted.append(
                    (
                        estimate,
                        subset,
                        survivors,
                    )
                )

            predicted.sort(
                key=lambda x: x[0]
            )

            print()
            print(
                "TOP PREDICTED SUBSETS"
            )

            for rank, (
                estimate,
                subset,
                survivors,
            ) in enumerate(
                predicted[:15],
                1,
            ):

                print(
                    f"{rank:2d}. "
                    f"{str(list(subset)):<38} "
                    f"pred={estimate:.9f}s "
                    f"survivors={survivors:,}"
                )

            # ------------------------------------------------
            # Select small empirical validation set.
            # ------------------------------------------------

            selected = choose_candidates(
                results,
                predicted,
            )

            validation_items = list(
                selected.values()
            )

            # Limit total validation.
            validation_items = (
                validation_items[
                    :FULL_BENCHMARKS
                ]
            )

            print()
            print(
                "FULL END-TO-END VALIDATION =",
                len(validation_items),
            )

            actual = []

            for (
                survivor_hint,
                subset,
                bits,
            ) in validation_items:

                # Verify target is retained.
                target_offset = (
                    true_x - start_x
                )

                target_byte = (
                    target_offset >> 3
                )

                target_bit = (
                    1 << (
                        target_offset & 7
                    )
                )

                if not (
                    bits[target_byte]
                    & target_bit
                ):

                    raise RuntimeError(
                        "Selected subset "
                        "eliminated true_x."
                    )

                (
                    elapsed,
                    survivors,
                    tests,
                    found_x,
                    correct,
                ) = (None, None, None, None, None)

                best = None

                for _ in range(
                    FULL_REPEATS
                ):

                    result = full_attack(
                        n,
                        start_x,
                        true_x,
                        subset,
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

                (
                    elapsed,
                    survivors,
                    tests,
                    found_x,
                    correct,
                ) = best

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
                    f"{str(list(subset)):<38} "
                    f"time={elapsed:.9f}s "
                    f"survivors={survivors:,} "
                    f"tests={tests:,} "
                    f"correct={correct}"
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

            # ------------------------------------------------
            # Depth winner.
            # ------------------------------------------------

            actual.sort(
                key=lambda x: x[0]
            )

            if actual:

                best_depth = actual[0]

                speedup = (
                    baseline_time
                    / best_depth[0]
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
                    "  speedup   =",
                    f"{speedup:.3f}x",
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
        # Final result.
        # ----------------------------------------------------

        print()
        print("=" * 78)
        print(
            f"GLOBAL BEST {bit_size}-BIT"
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
    print("FINISHED EXPERIMENT 203")
    print("=" * 78)


if __name__ == "__main__":
    main()

