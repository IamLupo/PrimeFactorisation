#!/usr/bin/env python3

"""
START EXPERIMENT 208

MEMORY-SAFE BEAM SEARCH FOR QR-PRIME SUBSET OPTIMIZATION

Exp 205:
    Greedy QR-prime selection.

Exp 206:
    1-for-1 local search.

Exp 207:
    2-for-2 local search.

Exp 208:
    Beam search.

IMPORTANT MEMORY FIX
--------------------
The original implementation stored a full NumPy boolean mask inside
every candidate generated during beam expansion.

For the 48-bit instance:

    search interval ~= 4.2 million x values
    one bool mask   ~= 4.2 MB

Thousands of temporary candidate masks therefore caused a memory
explosion.

This version:

    * stores masks only for the CURRENT beam;
    * creates each child mask temporarily;
    * counts survivors immediately;
    * only keeps the masks of candidates that enter the next beam;
    * never stores the entire expanded candidate population.

The search space is therefore much more memory controlled.

FINISHED EXPERIMENT 208
"""

import math
import statistics
import time

import numpy as np


# ============================================================
# CONFIGURATION
# ============================================================

# Keep this moderate because each retained beam element owns
# one full mask.
BEAM_WIDTH = 150

# Number of candidates benchmarked after beam search.
BENCHMARK_TOP = 30

# Runtime screening repetitions.
SCREEN_REPEATS = 5

# Final confirmation repetitions.
FINAL_REPEATS = 30

# Maximum subset size.
MAX_SIZE = 16


# ============================================================
# INSTANCES
# ============================================================

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


# ============================================================
# PRIME POOL
# ============================================================

PRIME_POOL = (
    3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47, 53,
    59, 61, 67, 71, 73
)


# ============================================================
# QR MASK GENERATION
# ============================================================

def build_qr_mask(
    n: int,
    start_x: int,
    end_x: int,
    prime: int,
):
    """
    Return a boolean mask indexed by:

        x = start_x + index

    True iff:

        x^2 - n

    is a quadratic residue modulo prime.
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
        if qr[(r * r - n_mod) % prime]:
            allowed[r] = True

    length = end_x - start_x + 1

    residues = (
        np.arange(
            length,
            dtype=np.int64,
        )
        + start_x
    ) % prime

    return allowed[residues]


# ============================================================
# MASK CACHE
# ============================================================

class MaskCache:

    def __init__(
        self,
        n,
        start_x,
        end_x,
        primes,
    ):

        self.masks = {}

        print()
        print("Building QR masks...")

        for p in primes:

            t0 = time.perf_counter()

            self.masks[p] = build_qr_mask(
                n,
                start_x,
                end_x,
                p,
            )

            t1 = time.perf_counter()

            print(
                "  prime={:<3d} build={:.6f}s".format(
                    p,
                    t1 - t0,
                )
            )

    def get(self, p):
        return self.masks[p]


# ============================================================
# CANDIDATE
# ============================================================

class Candidate:

    __slots__ = (
        "subset",
        "mask",
        "survivors",
    )

    def __init__(
        self,
        subset,
        mask,
        survivors,
    ):
        self.subset = subset
        self.mask = mask
        self.survivors = survivors


# ============================================================
# EXACT FERMAT
# ============================================================

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


# ============================================================
# END-TO-END ATTACK
# ============================================================

def attack(
    instance,
    subset,
    cache,
):
    n = instance["n"]

    start_x = instance["start_x"]
    end_x = instance["true_x"]

    length = end_x - start_x + 1

    mask = np.ones(
        length,
        dtype=np.bool_,
    )

    for p in subset:
        mask &= cache.get(p)

    indices = np.flatnonzero(mask)

    for idx in indices:

        x = start_x + int(idx)

        factor = exact_fermat(
            n,
            x,
        )

        if factor is not None:
            return (
                factor,
                x,
                len(indices),
            )

    return (
        None,
        None,
        len(indices),
    )


# ============================================================
# BENCHMARK
# ============================================================

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
        "median": statistics.median(samples),
        "mean": statistics.mean(samples),
        "min": min(samples),
        "max": max(samples),
        "survivors": survivors,
        "factor": factor,
        "x": found_x,
        "samples": samples,
    }


# ============================================================
# MEMORY-SAFE BEAM EXPANSION
# ============================================================

def expand_beam(
    beam,
    cache,
    prime_pool,
):
    """
    Expand the beam without retaining every child mask.

    OLD:

        generate thousands of children
        keep thousands of full masks
        sort all candidates

    NEW:

        generate child
        create temporary mask
        count survivors
        keep only the best BEAM_WIDTH children

    A small Python list is used as the bounded beam.
    """

    next_beam = []

    for parent in beam:

        subset = parent.subset
        last_prime = subset[-1]

        for p in prime_pool:

            if p <= last_prime:
                continue

            new_subset = subset + (p,)

            # Temporary child mask.
            new_mask = parent.mask & cache.get(p)

            survivors = int(
                np.count_nonzero(new_mask)
            )

            child = Candidate(
                new_subset,
                new_mask,
                survivors,
            )

            # ------------------------------------------------
            # Insert into bounded beam.
            # ------------------------------------------------

            next_beam.append(child)

            if len(next_beam) > BEAM_WIDTH * 2:

                next_beam.sort(
                    key=lambda c: (
                        c.survivors,
                        c.subset,
                    )
                )

                # Discard the worse half immediately.
                next_beam = next_beam[
                    :BEAM_WIDTH
                ]

    # Final trim.
    next_beam.sort(
        key=lambda c: (
            c.survivors,
            c.subset,
        )
    )

    next_beam = next_beam[
        :BEAM_WIDTH
    ]

    return next_beam


# ============================================================
# BEAM SEARCH
# ============================================================

def beam_search(
    bits,
    instance,
    cache,
):
    print()
    print("=" * 78)
    print(
        "BEAM SEARCH - {} BIT".format(bits)
    )
    print("=" * 78)

    # --------------------------------------------------------
    # Initial beam
    # --------------------------------------------------------

    beam = []

    for p in PRIME_POOL:

        mask = cache.get(p).copy()

        survivors = int(
            np.count_nonzero(mask)
        )

        beam.append(
            Candidate(
                (p,),
                mask,
                survivors,
            )
        )

    beam.sort(
        key=lambda c: (
            c.survivors,
            c.subset,
        )
    )

    beam = beam[
        :BEAM_WIDTH
    ]

    best_by_size = [
        beam[0]
    ]

    print(
        "k=1 best_survivors={} subset={}".format(
            beam[0].survivors,
            beam[0].subset,
        )
    )

    # --------------------------------------------------------
    # Progressive search
    # --------------------------------------------------------

    for k in range(
        2,
        MAX_SIZE + 1,
    ):

        t0 = time.perf_counter()

        beam = expand_beam(
            beam,
            cache,
            PRIME_POOL,
        )

        t1 = time.perf_counter()

        if not beam:
            break

        best = beam[0]

        best_by_size.append(
            best
        )

        print(
            "k={} best_survivors={} "
            "subset={} time={:.6f}s "
            "beam={}".format(
                k,
                best.survivors,
                best.subset,
                t1 - t0,
                len(beam),
            )
        )

        # ----------------------------------------------------
        # Progress / memory information.
        # ----------------------------------------------------

        mask_bytes = 0

        for candidate in beam:
            mask_bytes += candidate.mask.nbytes

        mask_mb = (
            mask_bytes
            / (1024.0 * 1024.0)
        )

        print(
            "    beam mask memory ~= {:.1f} MB".format(
                mask_mb
            )
        )

    # --------------------------------------------------------
    # Exp 207 reference
    # --------------------------------------------------------

    exp207 = tuple(
        sorted(instance["exp207"])
    )

    exp207_mask = None

    for p in exp207:

        if exp207_mask is None:
            exp207_mask = cache.get(p).copy()
        else:
            exp207_mask &= cache.get(p)

    exp207_survivors = int(
        np.count_nonzero(exp207_mask)
    )

    print()
    print("Exp 207 reference:")
    print(
        "  subset    = {}".format(
            exp207
        )
    )
    print(
        "  survivors = {}".format(
            exp207_survivors
        )
    )

    # --------------------------------------------------------
    # Build shortlist WITHOUT copying masks unnecessarily.
    # --------------------------------------------------------

    shortlist = []

    shortlist.extend(
        best_by_size
    )

    shortlist.extend(
        beam
    )

    shortlist.sort(
        key=lambda c: (
            c.survivors,
            len(c.subset),
            c.subset,
        )
    )

    # --------------------------------------------------------
    # Explicitly put Exp 207 reference into shortlist.
    # --------------------------------------------------------

    reference = Candidate(
        exp207,
        exp207_mask,
        exp207_survivors,
    )

    shortlist.append(reference)

    # --------------------------------------------------------
    # Deduplicate subsets.
    # --------------------------------------------------------

    unique = {}

    for candidate in shortlist:
        unique[candidate.subset] = candidate

    shortlist = list(
        unique.values()
    )

    shortlist.sort(
        key=lambda c: (
            c.survivors,
            len(c.subset),
            c.subset,
        )
    )

    return shortlist


# ============================================================
# RUN INSTANCE
# ============================================================

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
        PRIME_POOL,
    )

    shortlist = beam_search(
        bits,
        instance,
        cache,
    )

    # --------------------------------------------------------
    # Survivor ranking
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("GLOBAL SURVIVOR SHORTLIST")
    print("=" * 78)

    for i, candidate in enumerate(
        shortlist[:20],
        1,
    ):

        print(
            "{:2d}. k={:2d} "
            "survivors={:<6d} "
            "subset={}".format(
                i,
                len(candidate.subset),
                candidate.survivors,
                candidate.subset,
            )
        )

    # --------------------------------------------------------
    # Benchmark top survivor candidates.
    # --------------------------------------------------------

    candidates = shortlist[
        :BENCHMARK_TOP
    ]

    print()
    print("=" * 78)
    print(
        "END-TO-END BENCHMARK TOP {}".format(
            len(candidates)
        )
    )
    print("=" * 78)

    benchmark_results = []

    for candidate in candidates:

        result = benchmark(
            instance,
            candidate.subset,
            cache,
            SCREEN_REPEATS,
        )

        benchmark_results.append(
            (
                result["median"],
                candidate,
                result,
            )
        )

        print(
            "k={:2d} survivors={:<6d} "
            "median={:.9f}s "
            "subset={}".format(
                len(candidate.subset),
                result["survivors"],
                result["median"],
                candidate.subset,
            )
        )

    benchmark_results.sort(
        key=lambda x: x[0]
    )

    # --------------------------------------------------------
    # Runtime shortlist
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("RUNTIME SHORTLIST")
    print("=" * 78)

    for i, (
        median,
        candidate,
        result,
    ) in enumerate(
        benchmark_results[:15],
        1,
    ):

        print(
            "{:2d}. {:.9f}s "
            "survivors={:<6d} "
            "k={:2d} subset={}".format(
                i,
                median,
                candidate.survivors,
                len(candidate.subset),
                candidate.subset,
            )
        )

    # --------------------------------------------------------
    # Final winner
    # --------------------------------------------------------

    winner = benchmark_results[0][1]

    print()
    print("=" * 78)
    print("FINAL WINNER CONFIRMATION")
    print("=" * 78)

    print(
        "Winner subset:"
    )

    print(
        winner.subset
    )

    final = benchmark(
        instance,
        winner.subset,
        cache,
        FINAL_REPEATS,
    )

    print(
        "size      = {}".format(
            len(winner.subset)
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
    # Exp 207 direct comparison.
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
        "Exp208 median = {:.9f}s".format(
            final["median"]
        )
    )

    print(
        "Exp208 / Exp207 = {:.3f}x".format(
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
        "Exp208 survivors = {}".format(
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


# ============================================================
# MAIN
# ============================================================

def main():

    print("=" * 78)
    print("START EXPERIMENT 208")
    print("MEMORY-SAFE BEAM SEARCH")
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
    # FINAL SUMMARY
    # ========================================================

    print()
    print("=" * 78)
    print("EXP 208 SUMMARY")
    print("=" * 78)

    for bits in (
        48,
        54,
        60,
        66,
    ):

        instance = INSTANCES[bits]
        result = results[bits]

        winner = result["winner"]
        final = result["final"]
        exp207 = result["exp207"]

        print()
        print(
            "{}-BIT".format(bits)
        )

        print(
            "  winner    = {}".format(
                winner
            )
        )

        print(
            "  size      = {}".format(
                len(winner)
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
            "  Exp207    = {:.9f}s".format(
                exp207["median"]
            )
        )

        print(
            "  improvement_vs_207 = {:.3f}x".format(
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
    print("FINISHED EXPERIMENT 208")
    print("=" * 78)


if __name__ == "__main__":
    main()
