#!/usr/bin/env python3

"""
==============================================================================
START EXPERIMENT 210
DIRECT COST BEAM SEARCH
==============================================================================

Exp 205:
    Greedy end-to-end runtime optimization.

Exp 206:
    1-for-1 neighborhood search.

Exp 207:
    2-for-2 neighborhood search.

Exp 208:
    Beam search by survivor count.

Exp 209:
    Runtime-model beam search.

Problem with Exp 209:
    Calibration requires many full end-to-end attacks and becomes
    unnecessarily expensive before the actual search starts.

Exp 210:
    Remove calibration completely.

Use a direct empirical cost proxy:

        COST = survivors + lambda * number_of_primes

Sweep several lambda values.

This searches several hypotheses about the runtime tradeoff:

    lambda = 0
        pure survivor minimization

    lambda = 1
        one additional prime is equivalent to one survivor

    lambda = 5
        prime count matters more

    lambda = 10
        strongly penalize large subsets

    lambda = 25
        aggressively prefer compact subsets

    lambda = 50
        strong overhead penalty

The actual finalists are then benchmarked end-to-end.

Packed Python-integer masks are used throughout the search.

No expensive calibration stage.

==============================================================================
FINISHED EXPERIMENT 210
==============================================================================
"""

import math
import statistics
import time

import numpy as np


# ============================================================================
# CONFIGURATION
# ============================================================================

BEAM_WIDTH = 300

MAX_SIZE = 16

SCREEN_REPEATS = 3

FINAL_REPEATS = 20

BENCHMARK_TOP = 20

# Multiple runtime hypotheses.
LAMBDAS = (
    0.0,
    0.5,
    1.0,
    2.0,
    5.0,
    10.0,
    25.0,
    50.0,
)


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

        "exp207": (
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

        "exp207": (
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

        "exp207": (
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

        "exp207": (
            3, 5, 11, 17, 19, 31,
            41, 43, 53, 59, 61
        ),

        "baseline": 0.007534830,
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

def build_qr_mask(
    n: int,
    start_x: int,
    end_x: int,
    prime: int,
):
    """
    Build boolean QR mask.
    """

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


def to_int_mask(mask):
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
# MASK CACHE
# ============================================================================

class MaskCache:

    def __init__(
        self,
        n,
        start_x,
        end_x,
    ):
        self.masks = {}

        print()
        print("Building packed QR masks...")

        for p in PRIME_POOL:

            t0 = time.perf_counter()

            boolean_mask = build_qr_mask(
                n,
                start_x,
                end_x,
                p,
            )

            packed = to_int_mask(
                boolean_mask
            )

            t1 = time.perf_counter()

            size_mb = (
                max(
                    1,
                    (packed.bit_length() + 7) // 8,
                )
                / (1024.0 * 1024.0)
            )

            print(
                "  prime={:<3d} "
                "build={:.6f}s "
                "size={:.3f}MB".format(
                    p,
                    t1 - t0,
                    size_mb,
                )
            )

            self.masks[p] = packed

    def get(self, p):
        return self.masks[p]


# ============================================================================
# SUBSET OPERATIONS
# ============================================================================

def subset_mask(
    subset,
    cache,
):
    result = cache.get(
        subset[0]
    )

    for p in subset[1:]:
        result &= cache.get(p)

    return result


def survivors_for_subset(
    subset,
    cache,
):
    return subset_mask(
        subset,
        cache,
    ).bit_count()


# ============================================================================
# EXACT FERMAT
# ============================================================================

def exact_fermat(
    n: int,
    x: int,
):
    y2 = x * x - n

    if y2 < 0:
        return None

    y = math.isqrt(y2)

    if y * y != y2:
        return None

    p = x - y
    q = x + y

    if p > 1 and q > 1 and p * q == n:
        return p, q

    return None


# ============================================================================
# PACKED SURVIVOR ITERATION
# ============================================================================

def iter_set_bits(value):

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


# ============================================================================
# END-TO-END ATTACK
# ============================================================================

def attack(
    instance,
    subset,
    cache,
):
    n = instance["n"]
    start_x = instance["start_x"]

    mask = subset_mask(
        subset,
        cache,
    )

    survivors = mask.bit_count()

    for idx in iter_set_bits(mask):

        x = start_x + idx

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
# BENCHMARK
# ============================================================================

def benchmark(
    instance,
    subset,
    cache,
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
        "samples": samples,
    }


# ============================================================================
# CANDIDATE
# ============================================================================

class Candidate:

    __slots__ = (
        "subset",
        "mask",
        "survivors",
        "cost",
        "lam",
    )

    def __init__(
        self,
        subset,
        mask,
        survivors,
        cost,
        lam,
    ):
        self.subset = subset
        self.mask = mask
        self.survivors = survivors
        self.cost = cost
        self.lam = lam


# ============================================================================
# COST
# ============================================================================

def cost_function(
    subset_size,
    survivors,
    lam,
):
    return (
        survivors
        + lam * subset_size
    )


# ============================================================================
# BEAM EXPANSION
# ============================================================================

def expand_beam(
    beam,
    cache,
    lam,
):
    """
    Expand and rank candidates by:

        survivors + lambda * k
    """

    candidates = []

    for parent in beam:

        subset = parent.subset

        last = subset[-1]

        for p in PRIME_POOL:

            if p <= last:
                continue

            child_subset = (
                subset + (p,)
            )

            child_mask = (
                parent.mask
                & cache.get(p)
            )

            survivors = (
                child_mask.bit_count()
            )

            cost = cost_function(
                len(child_subset),
                survivors,
                lam,
            )

            candidates.append(
                Candidate(
                    child_subset,
                    child_mask,
                    survivors,
                    cost,
                    lam,
                )
            )

    candidates.sort(
        key=lambda c: (
            c.cost,
            c.survivors,
            len(c.subset),
            c.subset,
        )
    )

    return candidates[
        :BEAM_WIDTH
    ]


# ============================================================================
# BEAM SEARCH
# ============================================================================

def run_beam(
    bits,
    instance,
    cache,
    lam,
):
    print()
    print("-" * 78)
    print(
        "LAMBDA = {}".format(
            lam
        )
    )
    print("-" * 78)

    beam = []

    # --------------------------------------------------------
    # k = 1
    # --------------------------------------------------------

    for p in PRIME_POOL:

        mask = cache.get(p)

        survivors = mask.bit_count()

        cost = cost_function(
            1,
            survivors,
            lam,
        )

        beam.append(
            Candidate(
                (p,),
                mask,
                survivors,
                cost,
                lam,
            )
        )

    beam.sort(
        key=lambda c: (
            c.cost,
            c.survivors,
            c.subset,
        )
    )

    beam = beam[
        :BEAM_WIDTH
    ]

    print(
        "k=1 "
        "cost={:.2f} "
        "survivors={} "
        "subset={}".format(
            beam[0].cost,
            beam[0].survivors,
            beam[0].subset,
        )
    )

    best_by_size = [
        beam[0]
    ]

    # --------------------------------------------------------
    # Grow beam.
    # --------------------------------------------------------

    for k in range(
        2,
        MAX_SIZE + 1,
    ):

        t0 = time.perf_counter()

        beam = expand_beam(
            beam,
            cache,
            lam,
        )

        t1 = time.perf_counter()

        if not beam:
            break

        best = beam[0]

        best_by_size.append(
            best
        )

        print(
            "k={:<2d} "
            "cost={:<12.2f} "
            "survivors={:<8d} "
            "time={:.4f}s "
            "beam={}".format(
                k,
                best.cost,
                best.survivors,
                t1 - t0,
                len(beam),
            )
        )

    # --------------------------------------------------------
    # Add Exp 207 reference.
    # --------------------------------------------------------

    exp207 = tuple(
        sorted(instance["exp207"])
    )

    exp207_mask = subset_mask(
        exp207,
        cache,
    )

    exp207_survivors = (
        exp207_mask.bit_count()
    )

    exp207_cost = cost_function(
        len(exp207),
        exp207_survivors,
        lam,
    )

    best_by_size.append(
        Candidate(
            exp207,
            exp207_mask,
            exp207_survivors,
            exp207_cost,
            lam,
        )
    )

    # --------------------------------------------------------
    # Final shortlist.
    # --------------------------------------------------------

    unique = {}

    for c in best_by_size:
        unique[c.subset] = c

    for c in beam:
        unique[c.subset] = c

    shortlist = list(
        unique.values()
    )

    shortlist.sort(
        key=lambda c: (
            c.cost,
            c.survivors,
            len(c.subset),
        )
    )

    return shortlist


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

    cache = MaskCache(
        instance["n"],
        instance["start_x"],
        instance["true_x"],
    )

    # --------------------------------------------------------
    # Run all lambda beams.
    # --------------------------------------------------------

    all_candidates = {}

    for lam in LAMBDAS:

        shortlist = run_beam(
            bits,
            instance,
            cache,
            lam,
        )

        for candidate in shortlist:

            all_candidates[
                candidate.subset
            ] = candidate

    # --------------------------------------------------------
    # Survivor-based reference.
    # --------------------------------------------------------

    survivor_sorted = sorted(
        all_candidates.values(),
        key=lambda c: (
            c.survivors,
            len(c.subset),
        ),
    )

    print()
    print("=" * 78)
    print("LOWEST SURVIVOR CANDIDATES")
    print("=" * 78)

    for i, c in enumerate(
        survivor_sorted[:15],
        1,
    ):

        print(
            "{:2d}. survivors={:<8d} "
            "k={:2d} "
            "lambda={:<5g} "
            "subset={}".format(
                i,
                c.survivors,
                len(c.subset),
                c.lam,
                c.subset,
            )
        )

    # --------------------------------------------------------
    # Benchmark all candidates with low survivor counts.
    # --------------------------------------------------------

    benchmark_pool = (
        survivor_sorted[
            :BENCHMARK_TOP
        ]
    )

    print()
    print("=" * 78)
    print(
        "END-TO-END BENCHMARK TOP {}".format(
            len(benchmark_pool)
        )
    )
    print("=" * 78)

    actual = []

    for c in benchmark_pool:

        result = benchmark(
            instance,
            c.subset,
            cache,
            SCREEN_REPEATS,
        )

        actual.append(
            (
                result["median"],
                c,
                result,
            )
        )

        print(
            "median={:.9f}s "
            "survivors={:<8d} "
            "k={:2d} "
            "lambda={:<5g} "
            "subset={}".format(
                result["median"],
                result["survivors"],
                len(c.subset),
                c.lam,
                c.subset,
            )
        )

    actual.sort(
        key=lambda x: x[0]
    )

    # --------------------------------------------------------
    # Runtime ranking.
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("ACTUAL RUNTIME RANKING")
    print("=" * 78)

    for i, (
        median,
        c,
        result,
    ) in enumerate(
        actual[:15],
        1,
    ):

        print(
            "{:2d}. {:.9f}s "
            "survivors={:<8d} "
            "k={:2d} "
            "lambda={:<5g} "
            "subset={}".format(
                i,
                median,
                result["survivors"],
                len(c.subset),
                c.lam,
                c.subset,
            )
        )

    # --------------------------------------------------------
    # Final winner.
    # --------------------------------------------------------

    winner = actual[0][1]

    print()
    print("=" * 78)
    print("FINAL WINNER CONFIRMATION")
    print("=" * 78)

    final = benchmark(
        instance,
        winner.subset,
        cache,
        FINAL_REPEATS,
    )

    print(
        "Winner subset:"
    )

    print(
        winner.subset
    )

    print(
        "size      = {}".format(
            len(winner.subset)
        )
    )

    print(
        "lambda    = {}".format(
            winner.lam
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
        "survivors = {}".format(
            final["survivors"]
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

    # --------------------------------------------------------
    # Exp 207 comparison.
    # --------------------------------------------------------

    exp207_result = benchmark(
        instance,
        tuple(
            sorted(instance["exp207"])
        ),
        cache,
        FINAL_REPEATS,
    )

    print()
    print(
        "DIRECT EXP 207 COMPARISON"
    )

    print(
        "Exp207 median = {:.9f}s".format(
            exp207_result["median"]
        )
    )

    print(
        "Exp210 median = {:.9f}s".format(
            final["median"]
        )
    )

    print(
        "Exp210 / Exp207 = {:.3f}x".format(
            exp207_result["median"]
            / final["median"]
        )
    )

    print(
        "Exp207 survivors = {}".format(
            exp207_result["survivors"]
        )
    )

    print(
        "Exp210 survivors = {}".format(
            final["survivors"]
        )
    )

    if final["x"] != instance["true_x"]:

        print()
        print(
            "WARNING: TRUE X MISMATCH!"
        )

    return {
        "winner": winner.subset,
        "final": final,
        "exp207": exp207_result,
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print("START EXPERIMENT 210")
    print("DIRECT COST BEAM SEARCH")
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

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("EXP 210 SUMMARY")
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
        exp207 = result["exp207"]

        print()
        print(
            "{}-BIT".format(bits)
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
            "  vs Exp207 = {:.3f}x".format(
                exp207["median"]
                / final["median"]
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
    print("FINISHED EXPERIMENT 210")
    print("=" * 78)


if __name__ == "__main__":
    main()
