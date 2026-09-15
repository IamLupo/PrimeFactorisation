#!/usr/bin/env python3

import math
import time
from typing import Dict, Tuple, List

import numpy as np


# ============================================================
# START EXPERIMENT 205
#
# ADAPTIVE MARGINAL-GAIN QR SEARCH
#
# Exhaustive subset enumeration has become too expensive.
#
# Instead we build one subset incrementally.
#
# At each depth:
#
#     S -> S + {p}
#
# For every unused prime p we calculate:
#
#     survivors(S)
#     survivors(S + {p})
#
# and therefore the exact marginal reduction:
#
#     gain = survivors(S) - survivors(S+p)
#
# We rank candidates using:
#
#     gain / incremental sieve cost
#
# but then benchmark the best few candidates with the actual
# end-to-end Fermat attack.
#
# This gives us:
#
#     adaptive search
#     + empirical validation
#     + bounded memory
#     + dramatically less work than exhaustive enumeration.
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

MAX_DEPTH = 12

# How many candidate extensions are actually tested
# end-to-end at each depth.
EXTENSION_BENCHMARKS = 4

# Number of repetitions for end-to-end timing.
FULL_REPEATS = 2

# Stop if adding another prime no longer improves the actual
# measured runtime by at least this fraction.
MIN_IMPROVEMENT = 0.01


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

def build_prime_mask(
    start_x: int,
    length: int,
    prime: int,
    table: np.ndarray,
) -> np.ndarray:

    x = np.arange(
        start_x,
        start_x + length,
        dtype=np.int64,
    )

    allowed = table[
        x % prime
    ]

    return np.packbits(
        allowed,
        bitorder="little",
    )


# ------------------------------------------------------------
# ALL ONES
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
# APPLY ONE PRIME
# ------------------------------------------------------------

def apply_prime(
    bits: np.ndarray,
    prime: int,
    masks: Dict[int, np.ndarray],
) -> np.ndarray:

    child = bits.copy()

    child &= masks[prime]

    return child


# ------------------------------------------------------------
# ATTACK USING PRECOMPUTED MASK
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

    tests = 0
    found_x = -1

    for off_np in offsets:

        off = int(off_np)

        x = (
            start_x + off
        )

        d = (
            x * x
            - n
        )

        tests += 1

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
        tests,
        found_x,
        found_x == true_x,
    )


# ------------------------------------------------------------
# END-TO-END BENCHMARK
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
            or result[0] < best[0]
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
        "START EXPERIMENT 205"
    )
    print(
        "Adaptive marginal-gain QR search"
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

        # ----------------------------------------------------
        # Build prime masks.
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

            masks[prime] = (
                build_prime_mask(
                    start_x,
                    length,
                    prime,
                    qr,
                )
            )

        # ----------------------------------------------------
        # Verify target survives.
        # ----------------------------------------------------

        offset = (
            true_x
            - start_x
        )

        target_byte = (
            offset >> 3
        )

        target_bit = (
            1 << (
                offset & 7
            )
        )

        rejected = []

        for prime in PRIMES:

            if not (
                masks[prime][
                    target_byte
                ]
                & target_bit
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
            "True Fermat x survives "
            "all individual QR tests."
        )

        # ----------------------------------------------------
        # Baseline.
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
        # Current subset.
        # ----------------------------------------------------

        current_subset: Tuple[
            int, ...
        ] = tuple()

        current_bits = make_ones(
            length
        )

        current_survivors = length

        global_best = (
            baseline_time,
            tuple(),
            current_survivors,
            baseline_tests,
            baseline_x,
            True,
        )

        # ----------------------------------------------------
        # Build search depth by depth.
        # ----------------------------------------------------

        for depth in range(
            1,
            MAX_DEPTH + 1,
        ):

            print()
            print("=" * 78)
            print(
                f"DEPTH {depth}"
            )
            print("=" * 78)

            candidates = []

            # ------------------------------------------------
            # Evaluate marginal gains for every unused prime.
            # ------------------------------------------------

            for prime in PRIMES:

                if prime in current_subset:
                    continue

                child_bits = apply_prime(
                    current_bits,
                    prime,
                    masks,
                )

                child_survivors = (
                    popcount(
                        child_bits
                    )
                )

                gain = (
                    current_survivors
                    - child_survivors
                )

                # Score = candidates removed per added prime.
                #
                # Every candidate costs one additional packed
                # AND, so this is a simple marginal efficiency.
                score = (
                    gain
                    / max(
                        1,
                        child_survivors,
                    )
                )

                candidates.append(
                    (
                        score,
                        gain,
                        child_survivors,
                        prime,
                        child_bits,
                    )
                )

            candidates.sort(
                key=lambda x: (
                    -x[0],
                    x[2],
                )
            )

            print()
            print(
                "MARGINAL RANKING"
            )

            for rank, (
                score,
                gain,
                survivors,
                prime,
                _,
            ) in enumerate(
                candidates[:10],
                1,
            ):

                print(
                    f"{rank:2d}. "
                    f"+{prime:2d} "
                    f"gain={gain:>10,d} "
                    f"survivors={survivors:>10,d} "
                    f"score={score:.6f}"
                )

            # ------------------------------------------------
            # Empirically benchmark the best few extensions.
            # ------------------------------------------------

            actual = []

            print()
            print(
                "END-TO-END EXTENSION TESTS"
            )

            for (
                score,
                gain,
                predicted_survivors,
                prime,
                child_bits,
            ) in candidates[
                :EXTENSION_BENCHMARKS
            ]:

                (
                    runtime,
                    survivors,
                    tests,
                    found_x,
                    correct,
                ) = benchmark_bits(
                    n,
                    start_x,
                    true_x,
                    child_bits,
                    length,
                )

                subset = (
                    current_subset
                    + (prime,)
                )

                actual.append(
                    (
                        runtime,
                        subset,
                        survivors,
                        tests,
                        found_x,
                        correct,
                        gain,
                    )
                )

                print(
                    f"{str(list(subset)):<38} "
                    f"time={runtime:.9f}s "
                    f"survivors={survivors:,} "
                    f"tests={tests:,} "
                    f"correct={correct}"
                )

            actual.sort(
                key=lambda x: x[0]
            )

            if not actual:

                print(
                    "No candidate extensions."
                )

                break

            winner = actual[0]

            # ------------------------------------------------
            # Accept empirically best extension.
            # ------------------------------------------------

            (
                winner_runtime,
                winner_subset,
                winner_survivors,
                winner_tests,
                winner_x,
                winner_correct,
                winner_gain,
            ) = winner

            if not winner_correct:

                raise RuntimeError(
                    "Selected subset failed "
                    "to recover true factor."
                )

            improvement = (
                (
                    global_best[0]
                    - winner_runtime
                )
                / global_best[0]
            )

            # Find winner's packed mask.
            winner_bits = None

            for (
                score,
                gain,
                survivors,
                prime,
                child_bits,
            ) in candidates:

                if (
                    current_subset
                    + (prime,)
                    == winner_subset
                ):

                    winner_bits = child_bits
                    break

            if winner_bits is None:

                raise RuntimeError(
                    "Internal winner-mask mismatch."
                )

            # ------------------------------------------------
            # Accept.
            # ------------------------------------------------

            current_subset = (
                winner_subset
            )

            current_bits = winner_bits
            current_survivors = (
                winner_survivors
            )

            print()
            print(
                "ACCEPTED EXTENSION"
            )

            print(
                "  subset    =",
                list(current_subset),
            )

            print(
                "  runtime   =",
                f"{winner_runtime:.9f}s",
            )

            print(
                "  speedup   =",
                f"{baseline_time / winner_runtime:.3f}x",
            )

            print(
                "  survivors =",
                current_survivors,
            )

            print(
                "  square tests =",
                winner_tests,
            )

            # ------------------------------------------------
            # Update global winner.
            # ------------------------------------------------

            if (
                winner_runtime
                < global_best[0]
            ):

                global_best = (
                    winner_runtime,
                    current_subset,
                    current_survivors,
                    winner_tests,
                    winner_x,
                    winner_correct,
                )

            else:

                # We have reached a local minimum.
                print()
                print(
                    "No runtime improvement "
                    "at this depth."
                )

                print(
                    "Stopping adaptive search."
                )

                break

            # ------------------------------------------------
            # Don't continue once improvement becomes tiny.
            # ------------------------------------------------

            if improvement < MIN_IMPROVEMENT:

                print()
                print(
                    "Improvement below threshold:"
                )

                print(
                    f"  improvement = "
                    f"{improvement * 100:.3f}%"
                )

                print(
                    "Stopping."
                )

                break

        # ----------------------------------------------------
        # Final result.
        # ----------------------------------------------------

        print()
        print("=" * 78)
        print(
            f"FINAL RESULT {bitsize}-BIT"
        )
        print("=" * 78)

        (
            best_runtime,
            best_subset,
            best_survivors,
            best_tests,
            best_x,
            best_correct,
        ) = global_best

        print(
            "best subset    =",
            list(best_subset),
        )

        print(
            "runtime        =",
            f"{best_runtime:.9f}s",
        )

        print(
            "speedup        =",
            f"{baseline_time / best_runtime:.3f}x",
        )

        print(
            "survivors      =",
            best_survivors,
        )

        print(
            "square tests   =",
            best_tests,
        )

        print(
            "found_x        =",
            best_x,
        )

        print(
            "expected_x     =",
            true_x,
        )

        print(
            "correct        =",
            best_correct,
        )

        # ----------------------------------------------------
        # Fermat sanity check.
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
        "FINISHED EXPERIMENT 205"
    )
    print("=" * 78)


if __name__ == "__main__":
    main()
