#!/usr/bin/env python3

"""
==============================================================================
START EXPERIMENT 214

EXHAUSTIVE 1/2/3-PRIME EXTENSION SEARCH

Motivation
----------

Exp 213 expanded the QR prime pool to primes <= 199.

It showed that large primes can strongly reduce the survivor count:

    48-bit:
        694 -> 316 -> 141 -> 62 -> 24 -> 9 -> 3 -> 1

but runtime became WORSE.

This means:

    survivor minimization != runtime minimization

Exp 213 was also greedy:

    choose the best single next prime

That can miss interactions where:

    A alone  is not especially useful
    B alone  is not especially useful
    A+B      is very useful

Exp 214 therefore exhaustively evaluates combinations of:

    +1 prime
    +2 primes
    +3 primes

starting from the strong Exp 207 subset.

For every candidate:

    survivor count is computed cheaply

Then the strongest candidates are benchmarked end-to-end.

This gives us a direct test of short-range combinatorial interactions.

Prime pool:
    all primes <= 199

Base subsets:
    Exp 207 final subsets

==============================================================================
FINISHED EXPERIMENT 214
==============================================================================
"""

import itertools
import math
import statistics
import time

import numpy as np


# ============================================================================
# CONFIGURATION
# ============================================================================

POOL_LIMIT = 199

SEARCH_ORDERS = (
    1,
    2,
    3,
)

# Number of lowest-survivor candidates retained for benchmarking
TOP_SURVIVOR_CANDIDATES = 40

SCREEN_REPEATS = 3

FINAL_REPEATS = 30


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

        "base": (
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

        "base": (
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

        "base": (
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

        "base": (
            3, 5, 11, 17, 19, 31,
            41, 43, 53, 59, 61
        ),

        "baseline": 0.007534830,
    },
}


# ============================================================================
# PRIME GENERATION
# ============================================================================

def primes_up_to(limit):

    primes = []

    for n in range(
        2,
        limit + 1,
    ):

        if n == 2:

            primes.append(n)
            continue

        if n % 2 == 0:
            continue

        prime = True

        d = 3

        while d * d <= n:

            if n % d == 0:

                prime = False
                break

            d += 2

        if prime:
            primes.append(n)

    return tuple(primes)


PRIME_POOL = primes_up_to(
    POOL_LIMIT
)


# ============================================================================
# QR MASK
# ============================================================================

def build_bool_mask(
    n,
    start_x,
    end_x,
    prime,
):
    """
    Build boolean QR mask for:

        x^2 - n

    being a quadratic residue modulo prime.
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
            "Building QR masks for primes <= {}".format(
                POOL_LIMIT
            )
        )

        for p in PRIME_POOL:

            t0 = time.perf_counter()

            mask = build_bool_mask(
                n,
                start_x,
                end_x,
                p,
            )

            self.bool_masks[p] = mask

            packed = np.packbits(
                mask,
                bitorder="little",
            )

            self.int_masks[p] = (
                int.from_bytes(
                    packed.tobytes(),
                    "little",
                    signed=False,
                )
            )

            t1 = time.perf_counter()

            print(
                "  prime={:<3d} "
                "build={:.6f}s".format(
                    p,
                    t1 - t0,
                )
            )


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

    if p > 1 and q > 1:

        if p * q == n:

            return (
                p,
                q,
            )

    return None


# ============================================================================
# REPRESENTATION
# ============================================================================

def choose_representation(bits):

    if bits <= 54:

        return "bool"

    return "int"


# ============================================================================
# BASE MASK
# ============================================================================

def make_base_mask(
    base,
    cache,
    representation,
):

    if representation == "bool":

        mask = cache.bool_masks[
            base[0]
        ].copy()

        for p in base[1:]:

            mask &= cache.bool_masks[p]

        return mask

    else:

        mask = cache.int_masks[
            base[0]
        ]

        for p in base[1:]:

            mask &= cache.int_masks[p]

        return mask


# ============================================================================
# COUNT CHILD SURVIVORS
# ============================================================================

def count_child(
    base_mask,
    added,
    cache,
    representation,
):

    if representation == "bool":

        mask = base_mask.copy()

        for p in added:

            mask &= cache.bool_masks[p]

        return (
            int(
                np.count_nonzero(mask)
            ),
            mask,
        )

    else:

        mask = base_mask

        for p in added:

            mask &= cache.int_masks[p]

        return (
            mask.bit_count(),
            mask,
        )


# ============================================================================
# END-TO-END BOOL ATTACK
# ============================================================================

def attack_bool(
    instance,
    subset,
    cache,
):

    n = instance["n"]
    start_x = instance["start_x"]

    mask = make_base_mask(
        subset,
        cache,
        "bool",
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
# END-TO-END INT ATTACK
# ============================================================================

def attack_int(
    instance,
    subset,
    cache,
):

    n = instance["n"]
    start_x = instance["start_x"]

    mask = make_base_mask(
        subset,
        cache,
        "int",
    )

    survivors = (
        mask.bit_count()
    )

    byte_length = (
        cache.length + 7
    ) // 8

    raw = mask.to_bytes(
        byte_length,
        "little",
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
# EXHAUSTIVE EXTENSION SEARCH
# ============================================================================

def extension_search(
    bits,
    instance,
    cache,
    representation,
):

    base = tuple(
        sorted(
            instance["base"]
        )
    )

    base_set = set(base)

    available = tuple(
        p
        for p in PRIME_POOL
        if p not in base_set
    )

    print()
    print("=" * 78)
    print(
        "BASE INSTANCE - {} BIT".format(
            bits
        )
    )
    print("=" * 78)

    print(
        "Base subset:"
    )

    print(
        base
    )

    base_survivors = (
        count_child(
            make_base_mask(
                base,
                cache,
                representation,
            ),
            (),
            cache,
            representation,
        )[0]
    )

    print(
        "Base survivors = {}".format(
            base_survivors
        )
    )

    print()
    print(
        "Available additional primes:"
    )

    print(
        available
    )

    all_candidates = []

    base_mask = make_base_mask(
        base,
        cache,
        representation,
    )

    # --------------------------------------------------------
    # 1, 2, 3 additions
    # --------------------------------------------------------

    for order in SEARCH_ORDERS:

        total = math.comb(
            len(available),
            order,
        )

        print()
        print(
            "SEARCHING {}-PRIME EXTENSIONS "
            "({} candidates)".format(
                order,
                total,
            )
        )

        t0 = time.perf_counter()

        counter = 0

        for added in itertools.combinations(
            available,
            order,
        ):

            survivors, _ = count_child(
                base_mask,
                added,
                cache,
                representation,
            )

            candidate = tuple(
                sorted(
                    base + added
                )
            )

            all_candidates.append(
                (
                    survivors,
                    len(candidate),
                    added,
                    candidate,
                )
            )

            counter += 1

        t1 = time.perf_counter()

        print(
            "  evaluated = {}".format(
                counter
            )
        )

        print(
            "  time      = {:.6f}s".format(
                t1 - t0
            )
        )

    # --------------------------------------------------------
    # Survivor ranking
    # --------------------------------------------------------

    all_candidates.sort(
        key=lambda item: (
            item[0],
            item[1],
            item[2],
        )
    )

    # Deduplicate.
    unique = {}

    for item in all_candidates:

        candidate = item[3]

        if candidate not in unique:

            unique[candidate] = item

    candidates = list(
        unique.values()
    )

    candidates.sort(
        key=lambda item: (
            item[0],
            item[1],
            item[2],
        )
    )

    print()
    print("=" * 78)
    print("LOWEST-SURVIVOR EXTENSIONS")
    print("=" * 78)

    for i, item in enumerate(
        candidates[:30],
        1,
    ):

        survivors = item[0]
        size = item[1]
        added = item[2]
        candidate = item[3]

        print(
            "{:2d}. survivors={:<7d} "
            "added={} "
            "k={:2d} "
            "subset={}".format(
                i,
                survivors,
                added,
                size,
                candidate,
            )
        )

    # --------------------------------------------------------
    # Benchmark candidates
    # --------------------------------------------------------

    benchmark_candidates = (
        candidates[
            :TOP_SURVIVOR_CANDIDATES
        ]
    )

    print()
    print("=" * 78)
    print(
        "END-TO-END BENCHMARK TOP {}".format(
            len(benchmark_candidates)
        )
    )
    print("=" * 78)

    actual = []

    for item in benchmark_candidates:

        candidate = item[3]

        result = benchmark(
            instance,
            candidate,
            cache,
            representation,
            SCREEN_REPEATS,
        )

        actual.append(
            (
                result["median"],
                candidate,
                result,
                item[0],
                item[2],
            )
        )

        print(
            "runtime={:.9f}s "
            "survivors={:<7d} "
            "added={} "
            "k={:2d} "
            "subset={}".format(
                result["median"],
                result["survivors"],
                item[2],
                len(candidate),
                candidate,
            )
        )

    actual.sort(
        key=lambda item: item[0]
    )

    # --------------------------------------------------------
    # Runtime ranking
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("ACTUAL RUNTIME RANKING")
    print("=" * 78)

    for i, item in enumerate(
        actual[:20],
        1,
    ):

        runtime = item[0]
        candidate = item[1]
        result = item[2]
        survivors = item[3]
        added = item[4]

        print(
            "{:2d}. runtime={:.9f}s "
            "survivors={:<7d} "
            "added={} "
            "k={:2d} "
            "subset={}".format(
                i,
                runtime,
                survivors,
                added,
                len(candidate),
                candidate,
            )
        )

    # --------------------------------------------------------
    # Final winner
    # --------------------------------------------------------

    winner = actual[0][1]

    print()
    print("=" * 78)
    print("FINAL WINNER CONFIRMATION")
    print("=" * 78)

    print(
        "Winner:"
    )

    print(
        winner
    )

    final = benchmark(
        instance,
        winner,
        cache,
        representation,
        FINAL_REPEATS,
    )

    print(
        "size      = {}".format(
            len(winner)
        )
    )

    print(
        "survivors = {}".format(
            final["survivors"]
        )
    )

    print(
        "median    = {:.9f}s".format(
            final["median"]
        )
    )

    print(
        "mean      = {:.9f}s".format(
            final["mean"]
        )
    )

    print(
        "min       = {:.9f}s".format(
            final["min"]
        )
    )

    print(
        "max       = {:.9f}s".format(
            final["max"]
        )
    )

    print(
        "factor    = {}".format(
            final["factor"]
        )
    )

    print(
        "found x   = {}".format(
            final["x"]
        )
    )

    print(
        "true x    = {}".format(
            instance["true_x"]
        )
    )

    print(
        "vs Exp205 = {:.3f}x".format(
            instance["baseline"]
            / final["median"]
        )
    )

    return {
        "winner": winner,
        "final": final,
        "base_survivors": base_survivors,
        "candidate_count": len(candidates),
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "START EXPERIMENT 214"
    )
    print(
        "EXHAUSTIVE 1/2/3-PRIME EXTENSION SEARCH"
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
        print(
            "#" * 78
        )

        print(
            "# {}-BIT INSTANCE".format(
                bits
            )
        )

        print(
            "#" * 78
        )

        cache = Cache(
            instance["n"],
            instance["start_x"],
            instance["true_x"],
        )

        representation = (
            choose_representation(
                bits
            )
        )

        print()
        print(
            "Representation = {}".format(
                representation
            )
        )

        results[bits] = extension_search(
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
        "EXP 214 SUMMARY"
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
            "  winner    = {}".format(
                result["winner"]
            )
        )

        print(
            "  size      = {}".format(
                len(result["winner"])
            )
        )

        print(
            "  survivors = {}".format(
                final["survivors"]
            )
        )

        print(
            "  median    = {:.9f}s".format(
                final["median"]
            )
        )

        print(
            "  vs Exp205 = {:.3f}x".format(
                instance["baseline"]
                / final["median"]
            )
        )

        print(
            "  factor    = {}".format(
                final["factor"]
            )
        )

        print(
            "  x         = {}".format(
                final["x"]
            )
        )

    print()
    print("=" * 78)
    print(
        "FINISHED EXPERIMENT 214"
    )
    print("=" * 78)


if __name__ == "__main__":
    main()
