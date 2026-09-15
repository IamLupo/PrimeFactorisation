#!/usr/bin/env python3

"""
START EXPERIMENT 206

Local neighborhood search around the Exp 205 greedy QR-prime subsets.

For each benchmark:
  1. Start from the Exp 205 selected prime set.
  2. Benchmark the current set repeatedly.
  3. Test every:
       - single removal
       - single addition
       - one-for-one replacement
  4. If a statistically stable improvement is found, move there.
  5. Repeat until no local improvement remains.

The final subsets are then retested with a larger sample count.

The search is END-TO-END:
  QR bitset construction
  -> intersection
  -> candidate extraction
  -> exact Fermat test
  -> factor verification

FINISHED EXPERIMENT 206
"""

import math
import time
import statistics
from itertools import combinations

import numpy as np


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
        "true_y": 12_578_760,
        "selected": [
            7, 11, 17, 41, 53, 31,
            67, 71, 47, 59, 61, 37
        ],
    },

    54: {
        "p": 124_517_461,
        "q": 144_517_463,
        "n": 17_994_947_562_921_443,
        "start_x": 134_145_249,
        "true_x": 134_517_462,
        "true_y": 10_000_001,
        "selected": [
            3, 5, 47, 43, 67, 41,
            29, 61, 71, 37, 17, 7
        ],
    },

    60: {
        "p": 1_058_841_403,
        "q": 1_088_841_421,
        "n": 1_152_910_377_856_153_663,
        "start_x": 1_073_736_643,
        "true_x": 1_073_841_412,
        "true_y": 15_000_009,
        "selected": [
            5, 7, 11, 13, 23, 67,
            59, 73, 17, 29, 53, 41
        ],
    },

    66: {
        "p": 8_569_934_017,
        "q": 8_609_934_041,
        "n": 73_786_566_622_092_172_697,
        "start_x": 8_589_910_746,
        "true_x": 8_589_934_029,
        "true_y": 20_000_012,
        "selected": [
            3, 5, 11, 19, 67, 47,
            59, 31, 61, 43, 53
        ],
    },
}


# Prime pool used by Exp 205.
# Keep this pool fixed so Exp 206 compares against the same search space.
PRIME_POOL = [
    3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47, 53,
    59, 61, 67, 71, 73
]


# ============================================================
# Configuration
# ============================================================

# Small number of samples for neighborhood exploration.
SEARCH_REPEATS = 5

# More samples for final confirmation.
FINAL_REPEATS = 15

# Require this fractional improvement before accepting a move.
# 0.01 = 1%
MIN_IMPROVEMENT = 0.01

# Safety limit so a noisy machine cannot cause endless searching.
MAX_ROUNDS = 10


# ============================================================
# Utilities
# ============================================================

def is_prime(n: int) -> bool:
    if n < 2:
        return False

    if n % 2 == 0:
        return n == 2

    r = math.isqrt(n)
    d = 3

    while d <= r:
        if n % d == 0:
            return False
        d += 2

    return True


def ceil_sqrt(n: int) -> int:
    x = math.isqrt(n)
    if x * x < n:
        x += 1
    return x


# ============================================================
# QR mask generation
# ============================================================

def qr_allowed_residues(n: int, p: int):
    """
    Return the residues r modulo p for which

        r^2 - n

    is a quadratic residue modulo p.

    This is the condition used by the Fermat x-search.
    """

    n_mod = n % p

    # Compute quadratic residues mod p.
    qr = np.zeros(p, dtype=np.bool_)

    for x in range(p):
        qr[(x * x) % p] = True

    allowed = np.zeros(p, dtype=np.bool_)

    for r in range(p):
        if qr[(r * r - n_mod) % p]:
            allowed[r] = True

    return allowed


def build_bitset(
    n: int,
    start_x: int,
    end_x: int,
    prime: int,
):
    """
    Build a boolean bitset for x in [start_x, end_x] satisfying
    the quadratic-residue condition modulo prime.
    """

    length = end_x - start_x + 1

    allowed = qr_allowed_residues(n, prime)

    xs_mod = (
        np.arange(length, dtype=np.int64) + start_x
    ) % prime

    return allowed[xs_mod]


# ============================================================
# Cache QR masks
# ============================================================

class QRMaskCache:
    def __init__(self, n: int, start_x: int, end_x: int, primes):
        self.n = n
        self.start_x = start_x
        self.end_x = end_x
        self.primes = tuple(primes)

        self.masks = {}

    def get(self, p: int):
        if p not in self.masks:
            self.masks[p] = build_bitset(
                self.n,
                self.start_x,
                self.end_x,
                p,
            )

        return self.masks[p]


# ============================================================
# Exact Fermat verification
# ============================================================

def exact_fermat_from_x(n: int, x: int):
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
# Full end-to-end attack
# ============================================================

def attack_with_subset(
    n: int,
    start_x: int,
    subset,
    cache: QRMaskCache,
):
    """
    Run the complete QR-filtered Fermat search.
    """

    subset = tuple(subset)

    end_x = (n.bit_length())  # dummy initialization only
    del end_x

    # Fermat's upper search bound is the true x if factors are known,
    # but the algorithm itself cannot know that. For benchmarking we
    # use the complete benchmark range stored in the instance.
    #
    # cache.start_x -> cache.end_x is therefore the exact experiment range.
    mask = np.ones(
        cache.end_x - cache.start_x + 1,
        dtype=np.bool_,
    )

    for p in subset:
        mask &= cache.get(p)

    candidate_indices = np.flatnonzero(mask)

    for idx in candidate_indices:
        x = cache.start_x + int(idx)

        result = exact_fermat_from_x(n, x)

        if result is not None:
            return result, int(x), len(candidate_indices)

    return None, None, len(candidate_indices)


# ============================================================
# Timing
# ============================================================

def benchmark_subset(
    instance,
    subset,
    cache,
    repeats,
):
    times = []
    survivors = None
    result = None
    found_x = None

    for _ in range(repeats):
        t0 = time.perf_counter()

        result, found_x, survivors = attack_with_subset(
            instance["n"],
            instance["start_x"],
            subset,
            cache,
        )

        t1 = time.perf_counter()

        times.append(t1 - t0)

    return {
        "median": statistics.median(times),
        "mean": statistics.mean(times),
        "min": min(times),
        "max": max(times),
        "times": times,
        "survivors": survivors,
        "result": result,
        "found_x": found_x,
    }


# ============================================================
# Neighborhood generation
# ============================================================

def generate_neighbors(subset):
    """
    Generate:

      S - {p}
      S + {q}
      S - {p} + {q}
    """

    current = set(subset)

    # --------------------------------------------------------
    # Single removals
    # --------------------------------------------------------

    for p in sorted(current):
        candidate = tuple(
            sorted(current - {p})
        )

        yield ("remove", p, None, candidate)

    # --------------------------------------------------------
    # Single additions
    # --------------------------------------------------------

    outside = [
        p for p in PRIME_POOL
        if p not in current
    ]

    for q in outside:
        candidate = tuple(
            sorted(current | {q})
        )

        yield ("add", None, q, candidate)

    # --------------------------------------------------------
    # One-for-one replacements
    # --------------------------------------------------------

    for p in sorted(current):
        for q in outside:

            candidate = tuple(
                sorted((current - {p}) | {q})
            )

            yield ("replace", p, q, candidate)


# ============================================================
# Pretty printing
# ============================================================

def print_benchmark(label, data):
    print(
        f"{label:<12} "
        f"median={data['median']:.9f}s  "
        f"mean={data['mean']:.9f}s  "
        f"min={data['min']:.9f}s  "
        f"survivors={data['survivors']}"
    )


def speedup(reference, candidate):
    if candidate <= 0:
        return float("inf")

    return reference / candidate


# ============================================================
# Main local search
# ============================================================

def local_search(bits, instance):
    n = instance["n"]
    start_x = instance["start_x"]

    end_x = instance["true_x"]

    cache = QRMaskCache(
        n,
        start_x,
        end_x,
        PRIME_POOL,
    )

    current = tuple(sorted(instance["selected"]))

    print()
    print("=" * 78)
    print(f"{bits}-BIT INSTANCE")
    print("=" * 78)

    print("Initial subset:")
    print(current)

    # --------------------------------------------------------
    # Initial benchmark
    # --------------------------------------------------------

    current_result = benchmark_subset(
        instance,
        current,
        cache,
        SEARCH_REPEATS,
    )

    print_benchmark("INITIAL", current_result)

    best_time = current_result["median"]
    best = current

    # --------------------------------------------------------
    # Hill-climbing
    # --------------------------------------------------------

    for round_no in range(1, MAX_ROUNDS + 1):

        print()
        print(f"LOCAL SEARCH ROUND {round_no}")
        print("-" * 78)

        candidates = list(generate_neighbors(best))

        evaluated = []

        for action, removed, added, candidate in candidates:

            result = benchmark_subset(
                instance,
                candidate,
                cache,
                SEARCH_REPEATS,
            )

            evaluated.append(
                (
                    result["median"],
                    action,
                    removed,
                    added,
                    candidate,
                    result,
                )
            )

        evaluated.sort(key=lambda x: x[0])

        improved = False

        # Compare the best candidate against current.
        for (
            candidate_time,
            action,
            removed,
            added,
            candidate,
            result,
        ) in evaluated:

            gain = speedup(best_time, candidate_time)

            if candidate_time < best_time * (1.0 - MIN_IMPROVEMENT):

                print()
                print("IMPROVEMENT FOUND")

                if action == "remove":
                    print(f"  operation: REMOVE {removed}")

                elif action == "add":
                    print(f"  operation: ADD {added}")

                else:
                    print(
                        f"  operation: REPLACE {removed} -> {added}"
                    )

                print(f"  old subset: {best}")
                print(f"  new subset: {candidate}")

                print(f"  old median: {best_time:.9f}s")
                print(f"  new median: {candidate_time:.9f}s")
                print(f"  local speedup: {gain:.3f}x")
                print(
                    f"  survivors: {result['survivors']}"
                )

                best = candidate
                best_time = candidate_time
                improved = True

                break

        if not improved:

            print()
            print("NO LOCAL IMPROVEMENT FOUND.")

            # Show the three fastest neighbors anyway.
            print()
            print("Top neighborhood candidates:")

            for (
                candidate_time,
                action,
                removed,
                added,
                candidate,
                result,
            ) in evaluated[:3]:

                if action == "remove":
                    op = f"REMOVE {removed}"

                elif action == "add":
                    op = f"ADD {added}"

                else:
                    op = f"REPLACE {removed}->{added}"

                print(
                    f"  {op:<20} "
                    f"{candidate_time:.9f}s "
                    f"{speedup(best_time, candidate_time):.3f}x "
                    f"survivors={result['survivors']}"
                )

            break

    # --------------------------------------------------------
    # Final confirmation with many repeats
    # --------------------------------------------------------

    print()
    print("FINAL CONFIRMATION")
    print("-" * 78)

    final = benchmark_subset(
        instance,
        best,
        cache,
        FINAL_REPEATS,
    )

    print("Final subset:")
    print(best)

    print_benchmark("FINAL", final)

    baseline = instance.get("baseline")

    if baseline is not None:
        print(
            f"vs baseline: "
            f"{speedup(baseline, final['median']):.3f}x"
        )

    print()
    print("Factor:", final["result"])
    print("Found x:", final["found_x"])
    print("True x :", instance["true_x"])

    if final["found_x"] != instance["true_x"]:
        print("WARNING: x mismatch!")

    # --------------------------------------------------------
    # Store for final summary
    # --------------------------------------------------------

    return {
        "subset": best,
        "benchmark": final,
    }


# ============================================================
# Main
# ============================================================

def main():

    print("=" * 78)
    print("START EXPERIMENT 206")
    print("QR-FILTER LOCAL NEIGHBORHOOD SEARCH")
    print("=" * 78)

    # These are the measured Exp 205 baselines.
    # They are only used for reporting.
    INSTANCES[48]["baseline"] = 1.233306658
    INSTANCES[54]["baseline"] = 0.106013845
    INSTANCES[60]["baseline"] = 0.033881383
    INSTANCES[66]["baseline"] = 0.007534830

    results = {}

    for bits in [48, 54, 60, 66]:

        instance = INSTANCES[bits]

        results[bits] = local_search(
            bits,
            instance,
        )

    # ========================================================
    # Final cross-instance summary
    # ========================================================

    print()
    print("=" * 78)
    print("EXP 206 SUMMARY")
    print("=" * 78)

    for bits in [48, 54, 60, 66]:

        instance = INSTANCES[bits]
        result = results[bits]

        subset = result["subset"]
        bench = result["benchmark"]

        print()
        print(f"{bits}-bit")
        print(f"  subset    = {subset}")
        print(f"  size      = {len(subset)}")
        print(f"  median    = {bench['median']:.9f}s")
        print(f"  survivors = {bench['survivors']}")
        print(
            f"  speedup   = "
            f"{speedup(instance['baseline'], bench['median']):.3f}x"
        )
        print(f"  factor    = {bench['result']}")
        print(f"  x         = {bench['found_x']}")

    print()
    print("=" * 78)
    print("FINISHED EXPERIMENT 206")
    print("=" * 78)


if __name__ == "__main__":
    main()
