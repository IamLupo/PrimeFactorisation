#!/usr/bin/env python3

"""
==============================================================================
START EXPERIMENT 212

PACKED REPRESENTATION BENCHMARK

Goal
----
Exp 211 showed that Python-int masks are excellent for subset search,
but poor for end-to-end attack timing because converting the final huge
integer back to NumPy bits is expensive.

Exp 212 isolates the representation problem.

Three attack implementations are compared:

    A) BOOL
       NumPy bool mask per x.

    B) PACKED_NUMPY
       NumPy uint8 packed mask.
       All QR intersections happen directly on packed bytes.
       Candidate extraction uses np.unpackbits() only once.

    C) PYTHON_INT
       Python integer bitset.
       QR intersections use integer AND.
       Final candidate extraction converts the integer to bytes.

Each representation is tested on the same Exp 207 subset.

This experiment does NOT search for a new prime subset.

The objective is to identify the fastest correct representation before
continuing subset optimization.

==============================================================================
FINISHED EXPERIMENT 212
==============================================================================
"""

import math
import statistics
import time

import numpy as np


# ============================================================================
# CONFIGURATION
# ============================================================================

REPEATS = 30


# ============================================================================
# BENCHMARK INSTANCES
# ============================================================================

INSTANCES = {
    48: {
        "p": 8_390_069,
        "q": 33_547_589,
        "n": 281_466_586_493_641,
        "start_x": 16_776_966,
        "true_x": 20_968_829,

        "subset": (
            7, 11, 17, 29, 31, 37,
            41, 47, 53, 59, 71, 73
        ),
    },

    54: {
        "p": 124_517_461,
        "q": 144_517_463,
        "n": 17_994_947_562_921_443,
        "start_x": 134_145_249,
        "true_x": 134_517_462,

        "subset": (
            3, 5, 7, 17, 19, 23,
            29, 37, 41, 43, 47, 61, 67
        ),
    },

    60: {
        "p": 1_058_841_403,
        "q": 1_088_841_421,
        "n": 1_152_910_377_856_153_663,
        "start_x": 1_073_736_643,
        "true_x": 1_073_841_412,

        "subset": (
            5, 7, 11, 13, 17, 19,
            23, 41, 47, 53, 59, 73
        ),
    },

    66: {
        "p": 8_569_934_017,
        "q": 8_609_934_041,
        "n": 73_786_566_622_092_172_697,
        "start_x": 8_589_910_746,
        "true_x": 8_589_934_029,

        "subset": (
            3, 5, 11, 17, 19, 31,
            41, 43, 53, 59, 61
        ),
    },
}


# ============================================================================
# PRIME POOL
# ============================================================================

PRIME_POOL = (
    3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47, 53,
    59, 61, 67, 71, 73
)


# ============================================================================
# QR MASK GENERATION
# ============================================================================

def build_bool_mask(
    n: int,
    start_x: int,
    end_x: int,
    prime: int,
):
    n_mod = n % prime

    qr = np.zeros(
        prime,
        dtype=np.bool_,
    )

    for r in range(prime):
        qr[(r * r) % prime] = True

    allowed = np.zeros(
        prime,
        dtype=np.bool_,
    )

    for r in range(prime):

        if qr[
            (r * r - n_mod) % prime
        ]:
            allowed[r] = True

    length = (
        end_x
        - start_x
        + 1
    )

    residues = (
        np.arange(
            length,
            dtype=np.int64,
        )
        + start_x
    ) % prime

    return allowed[residues]


def build_packed_numpy_mask(
    boolean_mask,
):
    return np.packbits(
        boolean_mask,
        bitorder="little",
    )


def packed_numpy_to_int(
    packed,
):
    return int.from_bytes(
        packed.tobytes(),
        byteorder="little",
        signed=False,
    )


# ============================================================================
# MASK CACHE
# ============================================================================

class Cache:

    def __init__(
        self,
        n,
        start_x,
        end_x,
    ):
        self.bool_masks = {}
        self.packed_masks = {}
        self.int_masks = {}

        self.length = (
            end_x
            - start_x
            + 1
        )

        self.byte_length = (
            self.length + 7
        ) // 8

        print()
        print(
            "Building all three mask representations..."
        )

        for p in PRIME_POOL:

            boolean_mask = build_bool_mask(
                n,
                start_x,
                end_x,
                p,
            )

            packed_mask = (
                build_packed_numpy_mask(
                    boolean_mask
                )
            )

            integer_mask = (
                packed_numpy_to_int(
                    packed_mask
                )
            )

            self.bool_masks[p] = (
                boolean_mask
            )

            self.packed_masks[p] = (
                packed_mask
            )

            self.int_masks[p] = (
                integer_mask
            )

        bool_bytes = self.length

        packed_bytes = self.byte_length

        print()
        print(
            "Representation sizes:"
        )

        print(
            "  BOOL          = {:.3f} MB".format(
                bool_bytes
                / (1024.0 * 1024.0)
            )
        )

        print(
            "  PACKED_NUMPY  = {:.3f} MB".format(
                packed_bytes
                / (1024.0 * 1024.0)
            )
        )

        print(
            "  Python-int    = {:.3f} MB".format(
                packed_bytes
                / (1024.0 * 1024.0)
            )
        )


# ============================================================================
# EXACT FERMAT
# ============================================================================

def exact_fermat(
    n: int,
    x: int,
):
    y2 = (
        x * x
        - n
    )

    if y2 < 0:
        return None

    y = math.isqrt(y2)

    if y * y != y2:
        return None

    p = x - y
    q = x + y

    if p > 1 and q > 1 and p * q == n:
        return (
            p,
            q,
        )

    return None


# ============================================================================
# BOOL ATTACK
# ============================================================================

def attack_bool(
    instance,
    cache,
    subset,
):
    start_x = instance["start_x"]
    n = instance["n"]

    mask = np.ones(
        cache.length,
        dtype=np.bool_,
    )

    for p in subset:
        mask &= cache.bool_masks[p]

    survivors = int(
        np.count_nonzero(mask)
    )

    indices = np.flatnonzero(mask)

    for idx in indices:

        x = (
            start_x
            + int(idx)
        )

        factor = exact_fermat(
            n,
            x,
        )

        if factor is not None:
            return (
                factor,
                x,
                survivors,
            )

    return (
        None,
        None,
        survivors,
    )


# ============================================================================
# PACKED NUMPY ATTACK
# ============================================================================

def attack_packed_numpy(
    instance,
    cache,
    subset,
):
    start_x = instance["start_x"]
    n = instance["n"]

    # Important:
    # Work directly on packed uint8 bytes.
    mask = cache.packed_masks[
        subset[0]
    ].copy()

    for p in subset[1:]:

        np.bitwise_and(
            mask,
            cache.packed_masks[p],
            out=mask,
        )

    # Count set bits without unpacking the entire mask.
    #
    # np.unpackbits is only used once at the end.
    survivors = int(
        np.unpackbits(
            mask,
            bitorder="little",
        )[:cache.length].sum()
    )

    bits = np.unpackbits(
        mask,
        bitorder="little",
    )[:cache.length]

    indices = np.flatnonzero(bits)

    for idx in indices:

        x = (
            start_x
            + int(idx)
        )

        factor = exact_fermat(
            n,
            x,
        )

        if factor is not None:
            return (
                factor,
                x,
                survivors,
            )

    return (
        None,
        None,
        survivors,
    )


# ============================================================================
# OPTIMIZED PACKED NUMPY ATTACK
# ============================================================================

def attack_packed_numpy_fast(
    instance,
    cache,
    subset,
):
    """
    Same as PACKED_NUMPY but unpack only once.

    This avoids doing unpackbits twice.
    """

    start_x = instance["start_x"]
    n = instance["n"]

    mask = cache.packed_masks[
        subset[0]
    ].copy()

    for p in subset[1:]:

        np.bitwise_and(
            mask,
            cache.packed_masks[p],
            out=mask,
        )

    bits = np.unpackbits(
        mask,
        bitorder="little",
    )

    bits = bits[
        :cache.length
    ]

    survivors = int(
        bits.sum()
    )

    indices = np.flatnonzero(
        bits
    )

    for idx in indices:

        x = (
            start_x
            + int(idx)
        )

        factor = exact_fermat(
            n,
            x,
        )

        if factor is not None:
            return (
                factor,
                x,
                survivors,
            )

    return (
        None,
        None,
        survivors,
    )


# ============================================================================
# PYTHON INT ATTACK
# ============================================================================

def iter_integer_bits(
    value,
):
    while value:

        low = (
            value
            & -value
        )

        index = (
            low.bit_length()
            - 1
        )

        yield index

        value ^= low


def attack_python_int(
    instance,
    cache,
    subset,
):
    start_x = instance["start_x"]
    n = instance["n"]

    mask = (
        cache.int_masks[
            subset[0]
        ]
    )

    for p in subset[1:]:

        mask &= (
            cache.int_masks[p]
        )

    survivors = (
        mask.bit_count()
    )

    for idx in iter_integer_bits(
        mask
    ):

        if idx >= cache.length:
            break

        x = (
            start_x
            + idx
        )

        factor = exact_fermat(
            n,
            x,
        )

        if factor is not None:
            return (
                factor,
                x,
                survivors,
            )

    return (
        None,
        None,
        survivors,
    )


# ============================================================================
# TIMING
# ============================================================================

def benchmark_attack(
    name,
    fn,
    instance,
    cache,
    subset,
):
    samples = []

    final_factor = None
    final_x = None
    final_survivors = None

    for _ in range(REPEATS):

        t0 = time.perf_counter()

        factor, found_x, survivors = fn(
            instance,
            cache,
            subset,
        )

        t1 = time.perf_counter()

        samples.append(
            t1 - t0
        )

        final_factor = factor
        final_x = found_x
        final_survivors = survivors

    median = statistics.median(
        samples
    )

    mean = statistics.mean(
        samples
    )

    minimum = min(samples)

    maximum = max(samples)

    print()
    print(
        "{}".format(name)
    )

    print(
        "  median    = {:.9f}s".format(
            median
        )
    )

    print(
        "  mean      = {:.9f}s".format(
            mean
        )
    )

    print(
        "  min       = {:.9f}s".format(
            minimum
        )
    )

    print(
        "  max       = {:.9f}s".format(
            maximum
        )
    )

    print(
        "  survivors = {}".format(
            final_survivors
        )
    )

    print(
        "  factor    = {}".format(
            final_factor
        )
    )

    print(
        "  x         = {}".format(
            final_x
        )
    )

    return {
        "median": median,
        "mean": mean,
        "min": minimum,
        "max": maximum,
        "survivors": final_survivors,
        "factor": final_factor,
        "x": final_x,
        "samples": samples,
    }


# ============================================================================
# ONE INSTANCE
# ============================================================================

def run_instance(
    bits,
    instance,
):
    print()
    print("#" * 78)
    print(
        "# {}-BIT INSTANCE".format(bits)
    )
    print("#" * 78)

    cache = Cache(
        instance["n"],
        instance["start_x"],
        instance["true_x"],
    )

    subset = instance["subset"]

    print()
    print(
        "Test subset:"
    )

    print(
        subset
    )

    # --------------------------------------------------------
    # Verify survivor count agrees across representations.
    # --------------------------------------------------------

    bool_survivors = int(
        np.count_nonzero(
            np.logical_and.reduce(
                [
                    cache.bool_masks[p]
                    for p in subset
                ]
            )
        )
    )

    packed_temp = (
        cache.packed_masks[
            subset[0]
        ].copy()
    )

    for p in subset[1:]:

        np.bitwise_and(
            packed_temp,
            cache.packed_masks[p],
            out=packed_temp,
        )

    packed_bits = np.unpackbits(
        packed_temp,
        bitorder="little",
    )[:cache.length]

    packed_survivors = int(
        packed_bits.sum()
    )

    int_mask = (
        cache.int_masks[
            subset[0]
        ]
    )

    for p in subset[1:]:
        int_mask &= (
            cache.int_masks[p]
        )

    int_survivors = (
        int_mask.bit_count()
    )

    print()
    print(
        "Representation survivor check:"
    )

    print(
        "  BOOL         = {}".format(
            bool_survivors
        )
    )

    print(
        "  PACKED_NUMPY = {}".format(
            packed_survivors
        )
    )

    print(
        "  PYTHON_INT   = {}".format(
            int_survivors
        )
    )

    if not (
        bool_survivors
        == packed_survivors
        == int_survivors
    ):
        print()
        print(
            "ERROR: REPRESENTATIONS DISAGREE"
        )

    # --------------------------------------------------------
    # Benchmark BOOL
    # --------------------------------------------------------

    bool_result = benchmark_attack(
        "A) BOOL",
        attack_bool,
        instance,
        cache,
        subset,
    )

    # --------------------------------------------------------
    # Benchmark packed NumPy
    # --------------------------------------------------------

    packed_result = benchmark_attack(
        "B) PACKED_NUMPY",
        attack_packed_numpy,
        instance,
        cache,
        subset,
    )

    # --------------------------------------------------------
    # Benchmark optimized packed NumPy
    # --------------------------------------------------------

    packed_fast_result = benchmark_attack(
        "C) PACKED_NUMPY_FAST",
        attack_packed_numpy_fast,
        instance,
        cache,
        subset,
    )

    # --------------------------------------------------------
    # Benchmark Python int
    # --------------------------------------------------------

    int_result = benchmark_attack(
        "D) PYTHON_INT",
        attack_python_int,
        instance,
        cache,
        subset,
    )

    # --------------------------------------------------------
    # Comparison
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("REPRESENTATION COMPARISON")
    print("=" * 78)

    results = [
        (
            "BOOL",
            bool_result,
        ),
        (
            "PACKED_NUMPY",
            packed_result,
        ),
        (
            "PACKED_NUMPY_FAST",
            packed_fast_result,
        ),
        (
            "PYTHON_INT",
            int_result,
        ),
    ]

    fastest = min(
        result["median"]
        for _, result in results
    )

    for name, result in results:

        print(
            "{:<20} "
            "median={:.9f}s "
            "relative={:.3f}x "
            "survivors={}".format(
                name,
                result["median"],
                result["median"]
                / fastest,
                result["survivors"],
            )
        )

    return {
        "bool": bool_result,
        "packed": packed_result,
        "packed_fast": packed_fast_result,
        "int": int_result,
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "START EXPERIMENT 212"
    )
    print(
        "PACKED REPRESENTATION BENCHMARK"
    )
    print("=" * 78)

    results = {}

    for bits in (
        48,
        54,
        60,
        66,
    ):

        results[bits] = run_instance(
            bits,
            INSTANCES[bits],
        )

    # ========================================================
    # SUMMARY
    # ========================================================

    print()
    print("=" * 78)
    print(
        "EXP 212 SUMMARY"
    )
    print("=" * 78)

    for bits in (
        48,
        54,
        60,
        66,
    ):

        r = results[bits]

        print()
        print(
            "{}-BIT".format(bits)
        )

        for key in (
            "bool",
            "packed",
            "packed_fast",
            "int",
        ):

            result = r[key]

            print(
                "  {:<15} {:.9f}s".format(
                    key,
                    result["median"],
                )
            )

        print(
            "  survivors = {}".format(
                r["packed_fast"]["survivors"]
            )
        )

        print(
            "  factor    = {}".format(
                r["packed_fast"]["factor"]
            )
        )

        print(
            "  x         = {}".format(
                r["packed_fast"]["x"]
            )
        )

    print()
    print("=" * 78)
    print(
        "FINISHED EXPERIMENT 212"
    )
    print("=" * 78)


if __name__ == "__main__":
    main()
