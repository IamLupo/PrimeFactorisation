#!/usr/bin/env python3

"""
==============================================================================
START EXPERIMENT 209
RUNTIME-MODEL BEAM SEARCH WITH PACKED BITSETS
==============================================================================

Exp 205:
    Greedy QR-prime selection.

Exp 206:
    1-for-1 local search.

Exp 207:
    2-for-2 local search.

Exp 208:
    Survivor-count beam search.

Problem observed in Exp 208:
    minimizing survivors is not identical to minimizing runtime.

Exp 209:
    fit an empirical runtime model

        T ~= alpha + beta*k + gamma*survivors

    and use that model as the beam-search objective.

IMPORTANT IMPLEMENTATION FIX
-----------------------------
QR masks are stored as Python integers.

The actual attack ALSO operates entirely on packed Python integers:

    mask &= prime_mask
    survivors = mask.bit_count()

Only at the very end are survivor positions extracted.

This avoids the NumPy:

    bool &= uint8

casting error from the previous implementation.

==============================================================================
FINISHED EXPERIMENT 209
==============================================================================
"""

import math
import statistics
import time

import numpy as np


# ============================================================================
# CONFIGURATION
# ============================================================================

BEAM_WIDTH = 500

MAX_SIZE = 16

CALIBRATION_REPEATS = 7

SCREEN_REPEATS = 5

FINAL_REPEATS = 30

BENCHMARK_TOP = 40


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
# QR MASK
# ============================================================================

def build_qr_mask_bool(
    n: int,
    start_x: int,
    end_x: int,
    prime: int,
):
    """
    Construct boolean QR mask for x in:

        [start_x, end_x]
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


def bool_mask_to_int(mask):
    """
    Convert boolean vector into a Python integer.

    Bit i == mask[i].
    """

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
# PACKED MASK CACHE
# ============================================================================

class PackedMaskCache:

    def __init__(
        self,
        n,
        start_x,
        end_x,
        primes,
    ):

        self.masks = {}

        print()
        print(
            "Building PACKED QR masks..."
        )

        for p in primes:

            t0 = time.perf_counter()

            boolean_mask = build_qr_mask_bool(
                n,
                start_x,
                end_x,
                p,
            )

            packed = bool_mask_to_int(
                boolean_mask
            )

            t1 = time.perf_counter()

            packed_bytes = (
                max(
                    1,
                    (packed.bit_length() + 7) // 8
                )
            )

            print(
                "  prime={:<3d} "
                "build={:.6f}s "
                "packed_size={:.3f} MB".format(
                    p,
                    t1 - t0,
                    packed_bytes / (
                        1024.0 * 1024.0
                    ),
                )
            )

            self.masks[p] = packed

    def get(self, p):
        return self.masks[p]


# ============================================================================
# SUBSET MASK
# ============================================================================

def subset_mask(
    subset,
    cache,
):
    """
    Packed intersection of all prime masks.
    """

    result = cache.get(
        subset[0]
    )

    for p in subset[1:]:
        result &= cache.get(p)

    return result


def subset_survivors(
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
# PACKED MASK -> SURVIVOR POSITIONS
# ============================================================================

def survivor_indices(mask):
    """
    Yield set bit indices from a packed Python integer.

    This avoids unpacking the entire search interval into NumPy.
    """

    while mask:

        lowest = (
            mask
            & -mask
        )

        bit_index = (
            lowest.bit_length()
            - 1
        )

        yield bit_index

        mask ^= lowest


# ============================================================================
# END-TO-END ATTACK
# ============================================================================

def attack(
    instance,
    subset,
    cache,
):
    """
    Full attack using packed integers.

    1. Intersect QR masks.
    2. Enumerate surviving x positions.
    3. Perform exact Fermat verification.
    """

    start_x = instance["start_x"]
    n = instance["n"]

    mask = subset_mask(
        subset,
        cache,
    )

    survivors = mask.bit_count()

    for idx in survivor_indices(mask):

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
        "median": statistics.median(samples),
        "mean": statistics.mean(samples),
        "min": min(samples),
        "max": max(samples),
        "survivors": survivors,
        "factor": factor,
        "x": found_x,
        "samples": samples,
    }


# ============================================================================
# CALIBRATION SETS
# ============================================================================

def make_prefixes():

    result = []

    for k in range(
        1,
        MAX_SIZE + 1,
    ):

        result.append(
            tuple(
                PRIME_POOL[:k]
            )
        )

    return result


def make_calibration_sets(
    instance,
):
    """
    Use prefixes and nearby one-for-one variants of Exp 207.
    """

    candidates = []

    candidates.extend(
        make_prefixes()
    )

    exp207 = tuple(
        sorted(instance["exp207"])
    )

    candidates.append(
        exp207
    )

    base = set(exp207)

    outside = [
        p
        for p in PRIME_POOL
        if p not in base
    ]

    # Add nearby 1-for-1 variants.
    for a in sorted(base):

        for b in outside:

            candidate = tuple(
                sorted(
                    (base - {a}) | {b}
                )
            )

            candidates.append(
                candidate
            )

            if len(candidates) >= 50:
                break

        if len(candidates) >= 50:
            break

    # Deduplicate.
    unique = {}

    for subset in candidates:

        if len(subset) == 0:
            continue

        if len(subset) > MAX_SIZE:
            continue

        unique[subset] = True

    return list(
        unique.keys()
    )


# ============================================================================
# FIT RUNTIME MODEL
# ============================================================================

def fit_runtime_model(
    calibration_data,
):
    """
    Fit:

        T = alpha + beta*k + gamma*survivors

    using least squares.
    """

    X = []
    y = []

    for row in calibration_data:

        X.append([
            1.0,
            float(row["k"]),
            float(row["survivors"]),
        ])

        y.append(
            row["median"]
        )

    X = np.asarray(
        X,
        dtype=np.float64,
    )

    y = np.asarray(
        y,
        dtype=np.float64,
    )

    coefficients, _, _, _ = (
        np.linalg.lstsq(
            X,
            y,
            rcond=None,
        )
    )

    alpha = max(
        0.0,
        float(coefficients[0]),
    )

    beta = max(
        0.0,
        float(coefficients[1]),
    )

    gamma = max(
        0.0,
        float(coefficients[2]),
    )

    return (
        alpha,
        beta,
        gamma,
    )


def predict_runtime(
    k,
    survivors,
    model,
):
    alpha, beta, gamma = model

    return (
        alpha
        + beta * k
        + gamma * survivors
    )


# ============================================================================
# CANDIDATE
# ============================================================================

class Candidate:

    __slots__ = (
        "subset",
        "mask",
        "survivors",
        "predicted",
    )

    def __init__(
        self,
        subset,
        mask,
        survivors,
        predicted,
    ):
        self.subset = subset
        self.mask = mask
        self.survivors = survivors
        self.predicted = predicted


# ============================================================================
# RUNTIME BEAM EXPANSION
# ============================================================================

def expand_beam(
    beam,
    cache,
    model,
):
    """
    Expand the current beam.

    Ranking is based on predicted end-to-end runtime.

    The child mask is retained only when the child survives
    the beam cutoff.
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

            predicted = predict_runtime(
                len(child_subset),
                survivors,
                model,
            )

            candidates.append(
                Candidate(
                    child_subset,
                    child_mask,
                    survivors,
                    predicted,
                )
            )

    candidates.sort(
        key=lambda c: (
            c.predicted,
            c.survivors,
            len(c.subset),
            c.subset,
        )
    )

    # Keep only the strongest beam.
    return candidates[
        :BEAM_WIDTH
    ]


# ============================================================================
# RUNTIME-MODEL BEAM SEARCH
# ============================================================================

def run_beam_search(
    bits,
    instance,
    cache,
    model,
):
    print()
    print("=" * 78)
    print(
        "RUNTIME-MODEL BEAM SEARCH - {} BIT".format(
            bits
        )
    )
    print("=" * 78)

    # --------------------------------------------------------
    # k = 1
    # --------------------------------------------------------

    beam = []

    for p in PRIME_POOL:

        mask = cache.get(p)

        survivors = mask.bit_count()

        predicted = predict_runtime(
            1,
            survivors,
            model,
        )

        beam.append(
            Candidate(
                (p,),
                mask,
                survivors,
                predicted,
            )
        )

    beam.sort(
        key=lambda c: (
            c.predicted,
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
        "k=1 predicted={:.9f}s "
        "survivors={} subset={}".format(
            beam[0].predicted,
            beam[0].survivors,
            beam[0].subset,
        )
    )

    # --------------------------------------------------------
    # k = 2 ... MAX_SIZE
    # --------------------------------------------------------

    for k in range(
        2,
        MAX_SIZE + 1,
    ):

        t0 = time.perf_counter()

        beam = expand_beam(
            beam,
            cache,
            model,
        )

        t1 = time.perf_counter()

        if not beam:
            break

        best = beam[0]

        best_by_size.append(
            best
        )

        mask_memory = sum(
            max(
                1,
                (
                    c.mask.bit_length()
                    + 7
                ) // 8,
            )
            for c in beam
        )

        print(
            "k={} predicted={:.9f}s "
            "survivors={} subset={} "
            "time={:.6f}s beam={} "
            "memory~{:.2f}MB".format(
                k,
                best.predicted,
                best.survivors,
                best.subset,
                t1 - t0,
                len(beam),
                mask_memory
                / (
                    1024.0 * 1024.0
                ),
            )
        )

    # --------------------------------------------------------
    # Explicit Exp 207 reference
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

    shortlist = []

    shortlist.extend(
        best_by_size
    )

    shortlist.extend(
        beam
    )

    shortlist.append(
        Candidate(
            exp207,
            exp207_mask,
            exp207_survivors,
            predict_runtime(
                len(exp207),
                exp207_survivors,
                model,
            ),
        )
    )

    # Deduplicate.
    unique = {}

    for candidate in shortlist:
        unique[candidate.subset] = candidate

    shortlist = list(
        unique.values()
    )

    shortlist.sort(
        key=lambda c: (
            c.predicted,
            c.survivors,
            len(c.subset),
            c.subset,
        )
    )

    print()
    print(
        "Exp 207 reference:"
    )

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

    # --------------------------------------------------------
    # Packed cache
    # --------------------------------------------------------

    cache = PackedMaskCache(
        instance["n"],
        instance["start_x"],
        instance["true_x"],
        PRIME_POOL,
    )

    # --------------------------------------------------------
    # Calibration
    # --------------------------------------------------------

    calibration_sets = (
        make_calibration_sets(
            instance
        )
    )

    print()
    print("=" * 78)
    print("RUNTIME MODEL CALIBRATION")
    print("=" * 78)

    calibration_data = []

    for subset in calibration_sets:

        survivors = subset_survivors(
            subset,
            cache,
        )

        result = benchmark(
            instance,
            subset,
            cache,
            CALIBRATION_REPEATS,
        )

        row = {
            "subset": subset,
            "k": len(subset),
            "survivors": survivors,
            "median": result["median"],
        }

        calibration_data.append(
            row
        )

        print(
            "k={:2d} "
            "survivors={:<8d} "
            "median={:.9f}s "
            "subset={}".format(
                row["k"],
                row["survivors"],
                row["median"],
                row["subset"],
            )
        )

    model = fit_runtime_model(
        calibration_data
    )

    alpha, beta, gamma = model

    print()
    print("=" * 78)
    print("FITTED MODEL")
    print("=" * 78)

    print(
        "alpha = {:.12e}".format(
            alpha
        )
    )

    print(
        "beta  = {:.12e}".format(
            beta
        )
    )

    print(
        "gamma = {:.12e}".format(
            gamma
        )
    )

    # --------------------------------------------------------
    # Beam search
    # --------------------------------------------------------

    shortlist = run_beam_search(
        bits,
        instance,
        cache,
        model,
    )

    # --------------------------------------------------------
    # Predicted ranking
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("PREDICTED RUNTIME SHORTLIST")
    print("=" * 78)

    for i, candidate in enumerate(
        shortlist[:20],
        1,
    ):

        print(
            "{:2d}. predicted={:.9f}s "
            "k={:2d} survivors={:<8d} "
            "subset={}".format(
                i,
                candidate.predicted,
                len(candidate.subset),
                candidate.survivors,
                candidate.subset,
            )
        )

    # --------------------------------------------------------
    # Actual benchmark
    # --------------------------------------------------------

    candidates = shortlist[
        :BENCHMARK_TOP
    ]

    print()
    print("=" * 78)
    print(
        "ACTUAL END-TO-END BENCHMARK TOP {}".format(
            len(candidates)
        )
    )
    print("=" * 78)

    actual = []

    for candidate in candidates:

        result = benchmark(
            instance,
            candidate.subset,
            cache,
            SCREEN_REPEATS,
        )

        actual.append(
            (
                result["median"],
                candidate,
                result,
            )
        )

        print(
            "actual={:.9f}s "
            "predicted={:.9f}s "
            "k={:2d} "
            "survivors={:<8d} "
            "subset={}".format(
                result["median"],
                candidate.predicted,
                len(candidate.subset),
                result["survivors"],
                candidate.subset,
            )
        )

    actual.sort(
        key=lambda item: item[0]
    )

    # --------------------------------------------------------
    # Actual ranking
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("ACTUAL RUNTIME RANKING")
    print("=" * 78)

    for i, (
        median,
        candidate,
        result,
    ) in enumerate(
        actual[:20],
        1,
    ):

        print(
            "{:2d}. {:.9f}s "
            "pred={:.9f}s "
            "k={:2d} "
            "survivors={:<8d} "
            "subset={}".format(
                i,
                median,
                candidate.predicted,
                len(candidate.subset),
                candidate.survivors,
                candidate.subset,
            )
        )

    # --------------------------------------------------------
    # Final confirmation
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
    # Compare Exp 207
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
        "Exp209 median = {:.9f}s".format(
            final["median"]
        )
    )

    print(
        "Exp209 / Exp207 = {:.3f}x".format(
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
        "Exp209 survivors = {}".format(
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
        "model": model,
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print("START EXPERIMENT 209")
    print("RUNTIME-MODEL BEAM SEARCH")
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
    print("EXP 209 SUMMARY")
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
    print("FINISHED EXPERIMENT 209")
    print("=" * 78)


if __name__ == "__main__":
    main()
