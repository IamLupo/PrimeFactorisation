#!/usr/bin/env python3

"""
==============================================================================
START EXPERIMENT 213

EXPANDED QR-PRIME POOL + ADAPTIVE REPRESENTATION

Exp 212 established that the best mask representation depends on
the size of the search interval.

Exp 213 expands the prime pool from:

    primes <= 73

to:

    primes <= 199

and performs greedy QR-prime selection.

The search evaluates every unused prime at each step and chooses
the prime producing the smallest number of surviving Fermat x values.

The attack representation is selected automatically:

    <= 54 bits  -> NumPy bool
    60 bits     -> benchmark bool vs Python-int
    >= 66 bits  -> Python-int

The true factors are used ONLY for final verification.

==============================================================================
FINISHED EXPERIMENT 213
==============================================================================
"""

import math
import statistics
import time

import numpy as np


# ============================================================================
# CONFIGURATION
# ============================================================================

SEARCH_REPEATS = 3
FINAL_REPEATS = 30

# Maximum number of greedy additions.
#
# This is deliberately separate from the prime-pool size.
MAX_STEPS = 24


# ============================================================================
# PRIME GENERATION
# ============================================================================

def primes_up_to(limit):
    primes = []

    for n in range(2, limit + 1):

        if n < 2:
            continue

        if n == 2:
            primes.append(n)
            continue

        if n % 2 == 0:
            continue

        prime = True

        r = math.isqrt(n)

        d = 3

        while d <= r:

            if n % d == 0:
                prime = False
                break

            d += 2

        if prime:
            primes.append(n)

    return tuple(primes)


PRIME_POOL = primes_up_to(199)


# ============================================================================
# INSTANCES
# ============================================================================

INSTANCES = {
    48: {
        "p": 8_390_069,
        "q": 33_547_589,
        "n": 281_466_586_493_641,
        "start_x": 16_776_966,
        "true_x": 20_968_829,

        "seed": (
            7, 11, 17, 29, 31, 37,
            41, 47, 53, 59, 71, 73
        ),

        "baseline": 1.233306658,
    },

    54: {
        "p": 124_517_461,
        "q": 144_517_463,
        "n": 17_994_947_562_921_443,
        "start_x": 134_145_249,
        "true_x": 134_517_462,

        "seed": (
            3, 5, 7, 17, 19, 23,
            29, 37, 41, 43, 47, 61, 67
        ),

        "baseline": 0.106013845,
    },

    60: {
        "p": 1_058_841_403,
        "q": 1_088_841_421,
        "n": 1_152_910_377_856_153_663,
        "start_x": 1_073_736_643,
        "true_x": 1_073_841_412,

        "seed": (
            5, 7, 11, 13, 17, 19,
            23, 41, 47, 53, 59, 73
        ),

        "baseline": 0.033881383,
    },

    66: {
        "p": 8_569_934_017,
        "q": 8_609_934_041,
        "n": 73_786_566_622_092_172_697,
        "start_x": 8_589_910_746,
        "true_x": 8_589_934_029,

        "seed": (
            3, 5, 11, 17, 19, 31,
            41, 43, 53, 59, 61
        ),

        "baseline": 0.007534830,
    },
}


# ============================================================================
# QR MASK GENERATION
# ============================================================================

def build_bool_qr_mask(
    n,
    start_x,
    end_x,
    prime,
):
    """
    Build:

        allowed[x] = True

    iff:

        x^2 - n

    is a quadratic residue modulo prime.
    """

    n_mod = n % prime

    qr = np.zeros(
        prime,
        dtype=np.bool_,
    )

    for r in range(prime):

        qr[
            (r * r) % prime
        ] = True

    allowed = np.zeros(
        prime,
        dtype=np.bool_,
    )

    for r in range(prime):

        value = (
            r * r
            - n_mod
        ) % prime

        if qr[value]:
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


# ============================================================================
# BOOL -> PYTHON INT
# ============================================================================

def bool_to_int(mask):

    packed = np.packbits(
        mask,
        bitorder="little",
    )

    return int.from_bytes(
        packed.tobytes(),
        byteorder="little",
        signed=False,
    )


# ============================================================================
# CACHE
# ============================================================================

class Cache:

    def __init__(
        self,
        n,
        start_x,
        end_x,
    ):

        self.n = n
        self.start_x = start_x
        self.end_x = end_x

        self.length = (
            end_x
            - start_x
            + 1
        )

        self.bool_masks = {}
        self.int_masks = {}

        print()
        print(
            "Building QR masks for primes <= 199..."
        )

        for prime in PRIME_POOL:

            t0 = time.perf_counter()

            boolean_mask = build_bool_qr_mask(
                n,
                start_x,
                end_x,
                prime,
            )

            integer_mask = bool_to_int(
                boolean_mask
            )

            self.bool_masks[prime] = (
                boolean_mask
            )

            self.int_masks[prime] = (
                integer_mask
            )

            t1 = time.perf_counter()

            print(
                "  prime={:<3d} "
                "build={:.6f}s".format(
                    prime,
                    t1 - t0,
                )
            )


# ============================================================================
# REPRESENTATION SELECTION
# ============================================================================

def preferred_representation(bits):
    if bits <= 54:
        return "bool"

    if bits >= 66:
        return "int"

    return "auto"


def benchmark_representation(
    instance,
    subset,
    cache,
    representation,
    repeats,
):
    """
    Lightweight benchmark used only for the 60-bit crossover.
    """

    samples = []

    for _ in range(repeats):

        t0 = time.perf_counter()

        if representation == "bool":

            mask = cache.bool_masks[
                subset[0]
            ].copy()

            for p in subset[1:]:
                mask &= cache.bool_masks[p]

            survivors = int(
                np.count_nonzero(mask)
            )

        else:

            mask = cache.int_masks[
                subset[0]
            ]

            for p in subset[1:]:
                mask &= cache.int_masks[p]

            survivors = (
                mask.bit_count()
            )

        t1 = time.perf_counter()

        samples.append(
            t1 - t0
        )

    return statistics.median(
        samples
    )


def choose_representation(
    bits,
    instance,
    cache,
):
    selected = preferred_representation(
        bits
    )

    if selected != "auto":

        print(
            "Representation selected: {}".format(
                selected
            )
        )

        return selected

    subset = instance["seed"]

    print()
    print(
        "60-bit representation crossover"
    )

    bool_time = benchmark_representation(
        instance,
        subset,
        cache,
        "bool",
        10,
    )

    int_time = benchmark_representation(
        instance,
        subset,
        cache,
        "int",
        10,
    )

    print(
        "  BOOL = {:.9f}s".format(
            bool_time
        )
    )

    print(
        "  INT  = {:.9f}s".format(
            int_time
        )
    )

    if bool_time <= int_time:
        selected = "bool"
    else:
        selected = "int"

    print(
        "  selected = {}".format(
            selected
        )
    )

    return selected


# ============================================================================
# SUBSET MASK
# ============================================================================

def make_bool_mask(
    subset,
    cache,
):
    mask = cache.bool_masks[
        subset[0]
    ].copy()

    for p in subset[1:]:

        mask &= cache.bool_masks[p]

    return mask


def make_int_mask(
    subset,
    cache,
):
    mask = cache.int_masks[
        subset[0]
    ]

    for p in subset[1:]:

        mask &= cache.int_masks[p]

    return mask


# ============================================================================
# EXACT FERMAT
# ============================================================================

def exact_fermat(
    n,
    x,
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
    subset,
    cache,
):
    n = instance["n"]
    start_x = instance["start_x"]

    mask = make_bool_mask(
        subset,
        cache,
    )

    survivors = int(
        np.count_nonzero(mask)
    )

    indices = np.flatnonzero(
        mask
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

def attack_int(
    instance,
    subset,
    cache,
):
    n = instance["n"]
    start_x = instance["start_x"]

    mask = make_int_mask(
        subset,
        cache,
    )

    survivors = (
        mask.bit_count()
    )

    byte_length = (
        cache.length + 7
    ) // 8

    raw = mask.to_bytes(
        byte_length,
        byteorder="little",
        signed=False,
    )

    bits = np.unpackbits(
        np.frombuffer(
            raw,
            dtype=np.uint8,
        ),
        bitorder="little",
    )

    bits = bits[
        :cache.length
    ]

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
# GENERIC ATTACK
# ============================================================================

def attack(
    instance,
    subset,
    cache,
    representation,
):

    if representation == "bool":

        return attack_bool(
            instance,
            subset,
            cache,
        )

    return attack_int(
        instance,
        subset,
        cache,
    )


# ============================================================================
# BENCHMARK
# ============================================================================

def benchmark(
    instance,
    subset,
    cache,
    representation,
    repeats,
):
    samples = []

    factor = None
    found_x = None
    survivors = None

    for _ in range(repeats):

        t0 = time.perf_counter()

        factor, found_x, survivors = attack(
            instance,
            subset,
            cache,
            representation,
        )

        t1 = time.perf_counter()

        samples.append(
            t1 - t0
        )

    return {
        "median": statistics.median(
            samples
        ),
        "mean": statistics.mean(
            samples
        ),
        "min": min(samples),
        "max": max(samples),
        "survivors": survivors,
        "factor": factor,
        "x": found_x,
    }


# ============================================================================
# CHEAP SURVIVOR COUNT
# ============================================================================

def count_survivors(
    subset,
    cache,
    representation,
):
    if representation == "bool":

        mask = make_bool_mask(
            subset,
            cache,
        )

        return int(
            np.count_nonzero(mask)
        )

    mask = make_int_mask(
        subset,
        cache,
    )

    return mask.bit_count()


# ============================================================================
# GREEDY SEARCH
# ============================================================================

def greedy_search(
    bits,
    instance,
    cache,
    representation,
):
    print()
    print("=" * 78)
    print(
        "EXPANDED POOL GREEDY SEARCH - {} BIT".format(
            bits
        )
    )
    print("=" * 78)

    current = tuple(
        sorted(
            instance["seed"]
        )
    )

    used = set(current)

    # --------------------------------------------------------
    # Seed benchmark
    # --------------------------------------------------------

    current_result = benchmark(
        instance,
        current,
        cache,
        representation,
        FINAL_REPEATS,
    )

    current_time = (
        current_result["median"]
    )

    current_survivors = (
        current_result["survivors"]
    )

    print()
    print(
        "SEED"
    )

    print(
        "  subset    = {}".format(
            current
        )
    )

    print(
        "  size      = {}".format(
            len(current)
        )
    )

    print(
        "  survivors = {}".format(
            current_survivors
        )
    )

    print(
        "  runtime   = {:.9f}s".format(
            current_time
        )
    )

    # --------------------------------------------------------
    # Greedy additions
    # --------------------------------------------------------

    for step in range(
        1,
        MAX_STEPS + 1,
    ):

        # ----------------------------------------------------
        # Stop if every prime in the pool is already selected.
        # ----------------------------------------------------

        if len(used) >= len(PRIME_POOL):

            print()
            print(
                "All primes in pool have been selected."
            )

            break

        print()
        print(
            "STEP {}".format(
                step
            )
        )

        best_prime = None
        best_survivors = None
        best_subset = None

        # ----------------------------------------------------
        # Evaluate every unused prime.
        # ----------------------------------------------------

        for prime in PRIME_POOL:

            if prime in used:
                continue

            candidate = tuple(
                sorted(
                    current
                    + (prime,)
                )
            )

            survivors = count_survivors(
                candidate,
                cache,
                representation,
            )

            if (
                best_survivors is None
                or survivors < best_survivors
            ):

                best_survivors = survivors
                best_prime = prime
                best_subset = candidate

        if best_prime is None:
            break

        # ----------------------------------------------------
        # Benchmark candidate.
        # ----------------------------------------------------

        result = benchmark(
            instance,
            best_subset,
            cache,
            representation,
            SEARCH_REPEATS,
        )

        candidate_time = (
            result["median"]
        )

        if candidate_time > 0:

            speedup = (
                current_time
                / candidate_time
            )

        else:

            speedup = float("inf")

        if current_survivors > 0:

            survivor_ratio = (
                best_survivors
                / current_survivors
            )

        else:

            survivor_ratio = 0.0

        print(
            "  ADD {:>3d} -> "
            "survivors={:<8d} "
            "runtime={:.9f}s "
            "speedup={:.3f}x "
            "survivor_ratio={:.6f}".format(
                best_prime,
                best_survivors,
                candidate_time,
                speedup,
                survivor_ratio,
            )
        )

        # ----------------------------------------------------
        # Accept the greedy survivor optimum.
        # ----------------------------------------------------

        current = best_subset

        used.add(
            best_prime
        )

        current_survivors = (
            best_survivors
        )

        current_time = (
            candidate_time
        )

        # ----------------------------------------------------
        # Stop when the filter leaves one candidate.
        #
        # Further primes cannot remove the true x, so once
        # exactly one x survives there is little reason to
        # continue this greedy survivor search.
        # ----------------------------------------------------

        if current_survivors <= 1:

            print()
            print(
                "ONE SURVIVING x REMAINS."
            )

            print(
                "Stopping greedy expansion."
            )

            break

    # --------------------------------------------------------
    # Final confirmation
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL CONFIRMATION")
    print("=" * 78)

    final = benchmark(
        instance,
        current,
        cache,
        representation,
        FINAL_REPEATS,
    )

    print(
        "final subset = {}".format(
            current
        )
    )

    print(
        "size         = {}".format(
            len(current)
        )
    )

    print(
        "representation = {}".format(
            representation
        )
    )

    print(
        "survivors    = {}".format(
            final["survivors"]
        )
    )

    print(
        "median       = {:.9f}s".format(
            final["median"]
        )
    )

    print(
        "mean         = {:.9f}s".format(
            final["mean"]
        )
    )

    print(
        "min          = {:.9f}s".format(
            final["min"]
        )
    )

    print(
        "max          = {:.9f}s".format(
            final["max"]
        )
    )

    print(
        "factor       = {}".format(
            final["factor"]
        )
    )

    print(
        "found x      = {}".format(
            final["x"]
        )
    )

    print(
        "true x       = {}".format(
            instance["true_x"]
        )
    )

    print(
        "vs Exp205    = {:.3f}x".format(
            instance["baseline"]
            / final["median"]
        )
    )

    if final["x"] != instance["true_x"]:

        print()
        print(
            "WARNING: TRUE X MISMATCH!"
        )

    return {
        "subset": current,
        "final": final,
        "representation": representation,
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "START EXPERIMENT 213"
    )
    print(
        "EXPANDED QR-PRIME POOL <= 199"
    )
    print("=" * 78)

    print()
    print(
        "Prime pool:"
    )

    print(
        PRIME_POOL
    )

    print(
        "Prime count = {}".format(
            len(PRIME_POOL)
        )
    )

    results = {}

    for bits in (
        48,
        54,
        60,
        66,
    ):

        instance = INSTANCES[bits]

        print()
        print("#" * 78)
        print(
            "# {}-BIT INSTANCE".format(
                bits
            )
        )
        print("#" * 78)

        cache = Cache(
            instance["n"],
            instance["start_x"],
            instance["true_x"],
        )

        representation = (
            choose_representation(
                bits,
                instance,
                cache,
            )
        )

        results[bits] = greedy_search(
            bits,
            instance,
            cache,
            representation,
        )

    # ========================================================
    # SUMMARY
    # ========================================================

    print()
    print("=" * 78)
    print(
        "EXP 213 SUMMARY"
    )
    print("=" * 78)

    for bits in (
        48,
        54,
        60,
        66,
    ):

        instance = INSTANCES[bits]

        result = results[bits]

        final = result["final"]

        print()
        print(
            "{}-BIT".format(
                bits
            )
        )

        print(
            "  representation = {}".format(
                result["representation"]
            )
        )

        print(
            "  subset         = {}".format(
                result["subset"]
            )
        )

        print(
            "  size            = {}".format(
                len(result["subset"])
            )
        )

        print(
            "  survivors       = {}".format(
                final["survivors"]
            )
        )

        print(
            "  median          = {:.9f}s".format(
                final["median"]
            )
        )

        print(
            "  vs Exp205       = {:.3f}x".format(
                instance["baseline"]
                / final["median"]
            )
        )

        print(
            "  factor          = {}".format(
                final["factor"]
            )
        )

        print(
            "  x               = {}".format(
                final["x"]
            )
        )

    print()
    print("=" * 78)
    print(
        "FINISHED EXPERIMENT 213"
    )
    print("=" * 78)


if __name__ == "__main__":
    main()