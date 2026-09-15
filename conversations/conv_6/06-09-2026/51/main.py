#!/usr/bin/env python3

"""
START EXPERIMENT 207

QR-FILTER PAIRWISE LOCAL SEARCH

Purpose
-------
Exp 206 showed that the greedy Exp 205 solutions can sometimes
be improved by changing one prime.

Exp 207 asks a stronger question:

    Can TWO primes be replaced simultaneously in a way that
    produces an improvement which no individual replacement
    can discover?

Neighborhoods tested:

    1. Single removal:
         S - {a}

    2. Single addition:
         S + {b}

    3. One-for-one:
         S - {a} + {b}

    4. Two-for-two:
         S - {a,c} + {b,d}

The 2-for-2 neighborhood is the primary target.

Timing methodology
------------------
Because the end-to-end runtime is now extremely small for the
54/60/66-bit cases, candidates are screened with a small number
of repetitions and the best candidates are then re-measured with
a larger number of repetitions.

An improvement is accepted only after confirmation.

FINISHED EXPERIMENT 207
"""

import math
import statistics
import time
from itertools import combinations

import numpy as np


# ============================================================
# Configuration
# ============================================================

SCREEN_REPEATS = 3
CONFIRM_REPEATS = 15

# Number of candidates retained from the screening stage.
TOP_SCREEN_CANDIDATES = 12

# Require at least this fractional improvement before accepting.
# 0.02 = 2%.
MIN_IMPROVEMENT = 0.02

MAX_ROUNDS = 6


# ============================================================
# Benchmark instances
# ============================================================

INSTANCES = {
    48: {
        "p": 8_390_069,
        "q": 33_547_589,
        "n": 281_466_586_493_641,
        "start_x": 16_776_966,
        "true_x": 20_968_829,

        # Exp 206 final subset
        "selected": [
            7, 11, 17, 31, 37, 41,
            43, 47, 53, 59, 67, 71
        ],

        "baseline": 1.233306658,
    },

    54: {
        "p": 124_517_461,
        "q": 144_517_463,
        "n": 17_994_947_562_921_443,
        "start_x": 134_145_249,
        "true_x": 134_517_462,

        # Exp 206 final subset
        "selected": [
            3, 5, 7, 17, 29, 37, 41,
            43, 47, 53, 59, 61, 67
        ],

        "baseline": 0.106013845,
    },

    60: {
        "p": 1_058_841_403,
        "q": 1_088_841_421,
        "n": 1_152_910_377_856_153_663,
        "start_x": 1_073_736_643,
        "true_x": 1_073_841_412,

        # Exp 206 final subset
        "selected": [
            5, 7, 11, 13, 17, 23,
            29, 41, 53, 59, 67, 73
        ],

        "baseline": 0.033881383,
    },

    66: {
        "p": 8_569_934_017,
        "q": 8_609_934_041,
        "n": 73_786_566_622_092_172_697,
        "start_x": 8_589_910_746,
        "true_x": 8_589_934_029,

        # Exp 206 final subset
        "selected": [
            3, 5, 11, 17, 19, 31,
            43, 53, 59, 61, 67
        ],

        "baseline": 0.007534830,
    },
}


# ============================================================
# Fixed prime search pool
# ============================================================

PRIME_POOL = [
    3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47, 53,
    59, 61, 67, 71, 73
]


# ============================================================
# QR sieve
# ============================================================

def build_qr_mask(n: int, start_x: int, end_x: int, prime: int):
    """
    Boolean array indexed by

        x = start_x + index

    where True means

        x^2 - n

    is a quadratic residue modulo prime.
    """

    n_mod = n % prime

    qr = np.zeros(prime, dtype=np.bool_)

    for a in range(prime):
        qr[(a * a) % prime] = True

    allowed = np.zeros(prime, dtype=np.bool_)

    for r in range(prime):
        value = (r * r - n_mod) % prime

        if qr[value]:
            allowed[r] = True

    length = end_x - start_x + 1

    residues = (
        np.arange(length, dtype=np.int64) + start_x
    ) % prime

    return allowed[residues]


class MaskCache:
    def __init__(self, n, start_x, end_x, primes):
        self.masks = {}

        for p in primes:
            self.masks[p] = build_qr_mask(
                n,
                start_x,
                end_x,
                p,
            )

    def get(self, p):
        return self.masks[p]


# ============================================================
# Exact Fermat check
# ============================================================

def exact_fermat(n: int, x: int):
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
# Attack
# ============================================================

def attack(instance, subset, cache):
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

    candidate_indices = np.flatnonzero(mask)

    for idx in candidate_indices:
        x = start_x + int(idx)

        factor = exact_fermat(n, x)

        if factor is not None:
            return (
                factor,
                x,
                len(candidate_indices),
            )

    return (
        None,
        None,
        len(candidate_indices),
    )


# ============================================================
# Timing
# ============================================================

def benchmark(instance, subset, cache, repeats):
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

        samples.append(t1 - t0)

    return {
        "median": statistics.median(samples),
        "mean": statistics.mean(samples),
        "min": min(samples),
        "max": max(samples),
        "samples": samples,
        "factor": factor,
        "x": found_x,
        "survivors": survivors,
    }


# ============================================================
# Neighborhood generation
# ============================================================

def generate_neighborhood(subset):
    """
    Generate all local moves.

    Includes:

      remove
      add
      1-for-1
      2-for-2
    """

    current = set(subset)

    inside = sorted(current)
    outside = sorted(
        p for p in PRIME_POOL
        if p not in current
    )

    # --------------------------------------------------------
    # Single removals
    # --------------------------------------------------------

    for a in inside:

        candidate = tuple(
            sorted(current - {a})
        )

        yield (
            "REMOVE",
            (a,),
            (),
            candidate,
        )

    # --------------------------------------------------------
    # Single additions
    # --------------------------------------------------------

    for b in outside:

        candidate = tuple(
            sorted(current | {b})
        )

        yield (
            "ADD",
            (),
            (b,),
            candidate,
        )

    # --------------------------------------------------------
    # One-for-one
    # --------------------------------------------------------

    for a in inside:
        for b in outside:

            candidate = tuple(
                sorted((current - {a}) | {b})
            )

            yield (
                "1FOR1",
                (a,),
                (b,),
                candidate,
            )

    # --------------------------------------------------------
    # Two-for-two
    # --------------------------------------------------------

    for removed in combinations(inside, 2):

        base = current - set(removed)

        for added in combinations(outside, 2):

            candidate = tuple(
                sorted(base | set(added))
            )

            yield (
                "2FOR2",
                removed,
                added,
                candidate,
            )


# ============================================================
# Move formatting
# ============================================================

def describe_move(action, removed, added):

    if action == "REMOVE":
        return f"REMOVE {removed[0]}"

    if action == "ADD":
        return f"ADD {added[0]}"

    if action == "1FOR1":
        return (
            f"REPLACE "
            f"{removed[0]} -> {added[0]}"
        )

    if action == "2FOR2":
        return (
            f"2FOR2 "
            f"{removed[0]},{removed[1]}"
            f" -> "
            f"{added[0]},{added[1]}"
        )

    return action


# ============================================================
# Search one instance
# ============================================================

def optimize_instance(bits, instance):

    print()
    print("=" * 78)
    print(f"{bits}-BIT INSTANCE")
    print("=" * 78)

    # --------------------------------------------------------
    # Build cache
    # --------------------------------------------------------

    cache = MaskCache(
        instance["n"],
        instance["start_x"],
        instance["true_x"],
        PRIME_POOL,
    )

    current = tuple(
        sorted(instance["selected"])
    )

    # --------------------------------------------------------
    # Initial benchmark
    # --------------------------------------------------------

    current_result = benchmark(
        instance,
        current,
        cache,
        CONFIRM_REPEATS,
    )

    current_time = current_result["median"]

    print("Initial subset:")
    print(current)

    print(
        f"INITIAL"
        f" median={current_result['median']:.9f}s"
        f" mean={current_result['mean']:.9f}s"
        f" min={current_result['min']:.9f}s"
        f" survivors={current_result['survivors']}"
    )

    # --------------------------------------------------------
    # Local search
    # --------------------------------------------------------

    for round_no in range(1, MAX_ROUNDS + 1):

        print()
        print(f"LOCAL SEARCH ROUND {round_no}")
        print("-" * 78)

        screened = []

        # ----------------------------------------------------
        # Screen entire neighborhood
        # ----------------------------------------------------

        count = 0

        for (
            action,
            removed,
            added,
            candidate,
        ) in generate_neighborhood(current):

            result = benchmark(
                instance,
                candidate,
                cache,
                SCREEN_REPEATS,
            )

            screened.append(
                (
                    result["median"],
                    action,
                    removed,
                    added,
                    candidate,
                    result,
                )
            )

            count += 1

        screened.sort(
            key=lambda item: item[0]
        )

        print(
            f"Screened {count} candidates."
        )

        print()
        print("Best screened candidates:")

        for entry in screened[:8]:

            (
                candidate_time,
                action,
                removed,
                added,
                candidate,
                result,
            ) = entry

            move = describe_move(
                action,
                removed,
                added,
            )

            print(
                f"  {move:<28}"
                f" {candidate_time:.9f}s"
                f"  {current_time / candidate_time:.3f}x"
                f"  survivors={result['survivors']}"
            )

        # ----------------------------------------------------
        # Strong retest of best candidates
        # ----------------------------------------------------

        print()
        print(
            f"CONFIRMING TOP {TOP_SCREEN_CANDIDATES} "
            f"CANDIDATES"
        )

        confirmations = []

        for entry in screened[
            :TOP_SCREEN_CANDIDATES
        ]:

            (
                _screen_time,
                action,
                removed,
                added,
                candidate,
                _screen_result,
            ) = entry

            result = benchmark(
                instance,
                candidate,
                cache,
                CONFIRM_REPEATS,
            )

            confirmations.append(
                (
                    result["median"],
                    action,
                    removed,
                    added,
                    candidate,
                    result,
                )
            )

        confirmations.sort(
            key=lambda item: item[0]
        )

        # ----------------------------------------------------
        # Print confirmation table
        # ----------------------------------------------------

        print()

        for entry in confirmations[:8]:

            (
                candidate_time,
                action,
                removed,
                added,
                candidate,
                result,
            ) = entry

            move = describe_move(
                action,
                removed,
                added,
            )

            print(
                f"  {move:<28}"
                f" median={candidate_time:.9f}s"
                f" speedup={current_time / candidate_time:.3f}x"
                f" survivors={result['survivors']}"
            )

        # ----------------------------------------------------
        # Select winner
        # ----------------------------------------------------

        (
            best_time,
            best_action,
            best_removed,
            best_added,
            best_candidate,
            best_result,
        ) = confirmations[0]

        improvement = (
            (current_time - best_time)
            / current_time
        )

        # ----------------------------------------------------
        # Accept only meaningful improvement
        # ----------------------------------------------------

        if improvement > MIN_IMPROVEMENT:

            print()
            print("IMPROVEMENT CONFIRMED")

            print(
                "  operation:",
                describe_move(
                    best_action,
                    best_removed,
                    best_added,
                ),
            )

            print(
                "  old subset:",
                current,
            )

            print(
                "  new subset:",
                best_candidate,
            )

            print(
                f"  old median: "
                f"{current_time:.9f}s"
            )

            print(
                f"  new median: "
                f"{best_time:.9f}s"
            )

            print(
                f"  local speedup: "
                f"{current_time / best_time:.3f}x"
            )

            print(
                f"  improvement: "
                f"{improvement * 100:.2f}%"
            )

            print(
                f"  survivors: "
                f"{best_result['survivors']}"
            )

            current = best_candidate
            current_time = best_time

        else:

            print()
            print(
                "NO SIGNIFICANT LOCAL IMPROVEMENT."
            )

            print(
                f"Best confirmed improvement: "
                f"{improvement * 100:.2f}%"
            )

            break

    # ========================================================
    # Final validation
    # ========================================================

    print()
    print("FINAL CONFIRMATION")
    print("-" * 78)

    final = benchmark(
        instance,
        current,
        cache,
        30,
    )

    print("Final subset:")
    print(current)

    print(f"size      = {len(current)}")
    print(f"median    = {final['median']:.9f}s")
    print(f"mean      = {final['mean']:.9f}s")
    print(f"min       = {final['min']:.9f}s")
    print(f"max       = {final['max']:.9f}s")
    print(f"survivors = {final['survivors']}")
    print(f"factor    = {final['factor']}")
    print(f"found x   = {final['x']}")
    print(f"true x    = {instance['true_x']}")

    baseline = instance["baseline"]

    print(
        f"vs Exp205 baseline: "
        f"{baseline / final['median']:.3f}x"
    )

    if final["x"] != instance["true_x"]:
        print()
        print("WARNING: TRUE X MISMATCH!")

    return {
        "subset": current,
        "result": final,
    }


# ============================================================
# Main
# ============================================================

def main():

    print("=" * 78)
    print("START EXPERIMENT 207")
    print("QR-FILTER PAIRWISE LOCAL SEARCH")
    print("=" * 78)

    results = {}

    for bits in [48, 54, 60, 66]:

        results[bits] = optimize_instance(
            bits,
            INSTANCES[bits],
        )

    # ========================================================
    # Summary
    # ========================================================

    print()
    print("=" * 78)
    print("EXP 207 SUMMARY")
    print("=" * 78)

    for bits in [48, 54, 60, 66]:

        instance = INSTANCES[bits]
        result = results[bits]

        final_subset = result["subset"]
        final_result = result["result"]

        print()
        print(f"{bits}-bit")

        print(
            "  subset    =",
            final_subset,
        )

        print(
            f"  size      = "
            f"{len(final_subset)}"
        )

        print(
            f"  median    = "
            f"{final_result['median']:.9f}s"
        )

        print(
            f"  survivors = "
            f"{final_result['survivors']}"
        )

        print(
            f"  speedup   = "
            f"{instance['baseline'] / final_result['median']:.3f}x"
        )

        print(
            f"  factor    = "
            f"{final_result['factor']}"
        )

        print(
            f"  x         = "
            f"{final_result['x']}"
        )

    print()
    print("=" * 78)
    print("FINISHED EXPERIMENT 207")
    print("=" * 78)


if __name__ == "__main__":
    main()
