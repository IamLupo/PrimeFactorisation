#!/usr/bin/env python3

import itertools
import math
import statistics
import time

import numpy as np


# ============================================================================
# EXPERIMENT 215
# ============================================================================
#
# Goal:
#   Search for fast quadratic-residue prime subsets for Fermat factorization.
#
# Key idea:
#   Do NOT optimize survivor count alone.
#
#   Instead:
#       1. Exhaustively test 1/2/3-prime extensions.
#       2. Construct the Pareto frontier in:
#              (number of primes, survivors)
#       3. Benchmark the frontier and several nearby controls.
#       4. Rank by REAL end-to-end runtime.
#
# Prime 2 is deliberately excluded because for odd n:
#
#       x^2 - n == 0 (mod 2)
#
# for every integer x, so it supplies no filtering information.
#
# ============================================================================


# ============================================================================
# INSTANCES
# ============================================================================

INSTANCES = [
    {
        "name": "48-BIT",
        "n": 281466586493641,
        "p": 8390069,
        "q": 33547589,
        "start_x": 16776966,
        "true_x": 20968829,

        "representation": "bool",

        "base": (
            7, 11, 17, 29, 31, 37,
            41, 47, 53, 59, 71, 73
        ),

        "exp207": (
            7, 11, 17, 29, 31, 37,
            41, 47, 53, 59, 71, 73
        ),
    },

    {
        "name": "54-BIT",
        "n": 17994947562921443,
        "p": 124517461,
        "q": 144517463,
        "start_x": 134145249,
        "true_x": 134517462,

        "representation": "bool",

        "base": (
            3, 5, 7, 17, 19, 23, 29,
            37, 41, 43, 47, 61, 67
        ),

        "exp207": (
            3, 5, 7, 17, 19, 23, 29,
            37, 41, 43, 47, 61, 67
        ),
    },

    {
        "name": "60-BIT",
        "n": 1152910377856153663,
        "p": 1058841403,
        "q": 1088841421,
        "start_x": 1073736643,
        "true_x": 1073841412,

        "representation": "int",

        "base": (
            5, 7, 11, 13, 17, 19,
            23, 41, 47, 53, 59, 73
        ),

        "exp207": (
            5, 7, 11, 13, 17, 19,
            23, 41, 47, 53, 59, 73
        ),
    },

    {
        "name": "66-BIT",
        "n": 73786566622092172697,
        "p": 8569934017,
        "q": 8609934041,
        "start_x": 8589910746,
        "true_x": 8589934029,

        "representation": "int",

        "base": (
            3, 5, 11, 17, 19, 31,
            41, 43, 53, 59, 61
        ),

        "exp207": (
            3, 5, 11, 17, 19, 31,
            41, 43, 53, 59, 61
        ),
    },
]


# ============================================================================
# PRIME POOL
# ============================================================================

PRIME_POOL = (
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113,
    127, 131, 137, 139, 149, 151, 157, 163,
    167, 173, 179, 181, 191, 193, 197, 199,
)


# ============================================================================
# BASIC NUMBER THEORY
# ============================================================================

def is_square(v: int):
    """Return (True, sqrt(v)) if v is a perfect square."""
    if v < 0:
        return False, 0

    r = math.isqrt(v)

    if r * r == v:
        return True, r

    return False, r


def verify_factor(n, p, q):
    return p > 1 and q > 1 and p * q == n


# ============================================================================
# QR MASK GENERATION
# ============================================================================

def build_qr_mask_bool(n: int, start_x: int, length: int, p: int):
    """
    Build a boolean mask:

        mask[i] = True

    iff

        ((start_x + i)^2 - n) mod p

    is a quadratic residue modulo p.

    We use a residue lookup table containing:
        0^2, 1^2, ..., (p-1)^2 mod p
    """

    # Every residue modulo p that is a square.
    qr = np.zeros(p, dtype=np.bool_)

    residues = np.arange(p, dtype=np.int64)

    squares = (residues * residues) % p

    qr[squares] = True

    # Chunking prevents very large temporary allocations.
    mask = np.empty(length, dtype=np.bool_)

    chunk = 1_000_000

    for lo in range(0, length, chunk):
        hi = min(lo + chunk, length)

        xs = np.arange(
            start_x + lo,
            start_x + hi,
            dtype=np.int64,
        )

        # n can exceed signed 64-bit only for much larger experiments.
        # All current benchmark instances fit here.
        values = (xs * xs - n) % p

        mask[lo:hi] = qr[values]

    return mask


def bool_to_python_int(mask: np.ndarray):
    """
    Convert bool mask into a Python integer bitset.

    Bit i corresponds to mask[i].
    """

    packed = np.packbits(
        mask,
        bitorder="little",
    )

    value = int.from_bytes(
        packed.tobytes(),
        byteorder="little",
        signed=False,
    )

    return value


def build_masks(instance):
    """
    Build both representations for every prime.

    Returned:

        {
            "length": ...,
            "bool": { p: np.ndarray },
            "int":  { p: Python int },
        }
    """

    n = instance["n"]
    start_x = instance["start_x"]
    true_x = instance["true_x"]

    length = true_x - start_x + 1

    bool_masks = {}
    int_masks = {}

    print()
    print("Building QR masks for primes <= 199")

    for p in PRIME_POOL:
        t0 = time.perf_counter()

        mask = build_qr_mask_bool(
            n=n,
            start_x=start_x,
            length=length,
            p=p,
        )

        bool_masks[p] = mask

        # Build integer representation once.
        int_masks[p] = bool_to_python_int(mask)

        dt = time.perf_counter() - t0

        print(
            f"  prime={p:<3d} build={dt:.6f}s"
        )

    return {
        "length": length,
        "bool": bool_masks,
        "int": int_masks,
    }


# ============================================================================
# SUBSET MASKS
# ============================================================================

def make_subset_mask(
    masks,
    subset,
    representation,
):
    """
    Intersect all QR masks belonging to a prime subset.
    """

    if representation == "bool":

        # Start with all candidates alive.
        result = np.ones(
            masks["length"],
            dtype=np.bool_,
        )

        for p in subset:
            result &= masks["bool"][p]

            # Once empty, nothing can bring candidates back.
            if not result.any():
                break

        return result

    elif representation == "int":

        # Start with all bits alive.
        result = (1 << masks["length"]) - 1

        for p in subset:
            result &= masks["int"][p]

            if result == 0:
                break

        return result

    else:
        raise ValueError(
            f"Unknown representation: {representation}"
        )


def count_mask_survivors(mask, representation):
    if representation == "bool":
        return int(np.count_nonzero(mask))

    if representation == "int":
        return mask.bit_count()

    raise ValueError(
        f"Unknown representation: {representation}"
    )


def count_survivors(
    masks,
    subset,
    representation,
):
    mask = make_subset_mask(
        masks,
        subset,
        representation,
    )

    return count_mask_survivors(
        mask,
        representation,
    )


# ============================================================================
# FACTOR RECOVERY
# ============================================================================

def recover_from_bool_mask(
    n,
    start_x,
    mask,
):
    """
    Scan surviving x positions and perform the exact Fermat test.
    """

    indices = np.flatnonzero(mask)

    for idx in indices:
        x = start_x + int(idx)

        y2 = x * x - n

        if y2 < 0:
            continue

        y = math.isqrt(y2)

        if y * y != y2:
            continue

        p = x - y
        q = x + y

        if verify_factor(n, p, q):
            return p, q, x

    return None, None, None


def recover_from_int_mask(
    n,
    start_x,
    mask,
):
    """
    Scan set bits of a Python integer bitset.

    This is particularly effective for the 60/66-bit experiments.
    """

    while mask:
        lowbit = mask & -mask

        idx = lowbit.bit_length() - 1

        x = start_x + idx

        y2 = x * x - n

        if y2 >= 0:
            y = math.isqrt(y2)

            if y * y == y2:
                p = x - y
                q = x + y

                if verify_factor(n, p, q):
                    return p, q, x

        mask ^= lowbit

    return None, None, None


def solve_with_subset(
    instance,
    masks,
    subset,
):
    representation = instance["representation"]

    mask = make_subset_mask(
        masks,
        subset,
        representation,
    )

    n = instance["n"]
    start_x = instance["start_x"]

    if representation == "bool":
        return recover_from_bool_mask(
            n,
            start_x,
            mask,
        )

    return recover_from_int_mask(
        n,
        start_x,
        mask,
    )


# ============================================================================
# PARETO FRONTIER
# ============================================================================

def pareto_frontier(candidates):
    """
    candidates:

        [(subset, survivors), ...]

    Objective dimensions:

        k = len(subset)
        survivors

    Candidate A dominates candidate B iff:

        k_A <= k_B
        survivors_A <= survivors_B

    and at least one inequality is strict.
    """

    ordered = sorted(
        candidates,
        key=lambda item: (
            len(item[0]),
            item[1],
        ),
    )

    frontier = []

    best_survivors = math.inf

    for subset, survivors in ordered:

        if survivors < best_survivors:
            frontier.append(
                (subset, survivors)
            )

            best_survivors = survivors

    return frontier


# ============================================================================
# BENCHMARKING
# ============================================================================

def benchmark_subset(
    instance,
    masks,
    subset,
    repeats=11,
):
    samples = []

    factor = None

    for _ in range(repeats):

        t0 = time.perf_counter()

        factor = solve_with_subset(
            instance,
            masks,
            subset,
        )

        dt = time.perf_counter() - t0

        samples.append(dt)

    return {
        "median": statistics.median(samples),
        "mean": statistics.mean(samples),
        "min": min(samples),
        "max": max(samples),
        "factor": factor,
    }


# ============================================================================
# MAIN EXPERIMENT
# ============================================================================

def run_instance(instance):

    name = instance["name"]
    base = tuple(sorted(instance["base"]))
    exp207 = tuple(sorted(instance["exp207"]))

    print()
    print("#" * 78)
    print(f"# {name} INSTANCE")
    print("#" * 78)

    masks = build_masks(instance)

    representation = instance["representation"]

    print()
    print(f"Representation = {representation}")

    # ------------------------------------------------------------------------
    # BASE
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(f"BASE INSTANCE - {name}")
    print("=" * 78)

    base_survivors = count_survivors(
        masks,
        base,
        representation,
    )

    print()
    print("Base subset:")
    print(base)

    print(
        f"Base size      = {len(base)}"
    )

    print(
        f"Base survivors = {base_survivors}"
    )

    # ------------------------------------------------------------------------
    # AVAILABLE PRIMES
    # ------------------------------------------------------------------------

    available = [
        p for p in PRIME_POOL
        if p not in base
    ]

    print()
    print("Available additional primes:")
    print(tuple(available))

    # ------------------------------------------------------------------------
    # EXHAUSTIVE SEARCH
    # ------------------------------------------------------------------------

    candidates = {}

    # Base itself must always be a candidate.
    candidates[base] = base_survivors

    # Explicitly preserve Exp 207 baseline.
    exp207_survivors = count_survivors(
        masks,
        exp207,
        representation,
    )

    candidates[exp207] = exp207_survivors

    for r in (1, 2, 3):

        print()
        print(
            f"SEARCHING {r}-PRIME EXTENSIONS"
        )

        t0 = time.perf_counter()

        evaluated = 0

        for added in itertools.combinations(
            available,
            r,
        ):

            subset = tuple(
                sorted(
                    base + tuple(added)
                )
            )

            survivors = count_survivors(
                masks,
                subset,
                representation,
            )

            candidates[subset] = survivors

            evaluated += 1

        dt = time.perf_counter() - t0

        print(
            f"  evaluated = {evaluated}"
        )

        print(
            f"  time      = {dt:.6f}s"
        )

    candidate_list = list(
        candidates.items()
    )

    # ------------------------------------------------------------------------
    # PARETO FRONTIER
    # ------------------------------------------------------------------------

    frontier = pareto_frontier(
        candidate_list
    )

    print()
    print("=" * 78)
    print("PARETO FRONTIER")
    print("=" * 78)

    for i, (subset, survivors) in enumerate(
        frontier,
        1,
    ):
        print(
            f"{i:3d}. "
            f"k={len(subset):2d} "
            f"survivors={survivors:<8d} "
            f"subset={subset}"
        )

    # ------------------------------------------------------------------------
    # CREATE BENCHMARK SET
    # ------------------------------------------------------------------------
    #
    # Benchmark:
    #
    #   - Every Pareto candidate.
    #   - Base.
    #   - Exp 207.
    #   - Best 5 survivor candidates at each subset size.
    #
    # This keeps the benchmark manageable while allowing the runtime optimum
    # to differ from the absolute survivor optimum.
    # ------------------------------------------------------------------------

    benchmark_subsets = {
        subset: survivors
        for subset, survivors in frontier
    }

    benchmark_subsets[base] = base_survivors
    benchmark_subsets[exp207] = exp207_survivors

    by_k = {}

    for subset, survivors in candidate_list:
        k = len(subset)

        by_k.setdefault(k, [])

        by_k[k].append(
            (subset, survivors)
        )

    for k, values in by_k.items():

        values.sort(
            key=lambda item: item[1]
        )

        for subset, survivors in values[:5]:
            benchmark_subsets[subset] = survivors

    benchmark_list = list(
        benchmark_subsets.items()
    )

    benchmark_list.sort(
        key=lambda item: (
            len(item[0]),
            item[1],
        )
    )

    print()
    print("=" * 78)
    print(
        "END-TO-END BENCHMARK"
    )
    print("=" * 78)

    results = []

    for i, (subset, survivors) in enumerate(
        benchmark_list,
        1,
    ):

        result = benchmark_subset(
            instance,
            masks,
            subset,
        )

        factor_p, factor_q, found_x = (
            result["factor"]
        )

        results.append(
            {
                "subset": subset,
                "survivors": survivors,
                "median": result["median"],
                "mean": result["mean"],
                "min": result["min"],
                "max": result["max"],
                "factor": result["factor"],
            }
        )

        print(
            f"{i:3d}. "
            f"runtime={result['median']:.9f}s "
            f"k={len(subset):2d} "
            f"survivors={survivors:<8d} "
            f"subset={subset}"
        )

        # Immediate correctness check.
        if factor_p is not None:

            assert (
                factor_p * factor_q
                == instance["n"]
            )

            assert (
                found_x
                == instance["true_x"]
            )

    # ------------------------------------------------------------------------
    # RUNTIME RANKING
    # ------------------------------------------------------------------------

    results.sort(
        key=lambda item: item["median"]
    )

    print()
    print("=" * 78)
    print("ACTUAL RUNTIME RANKING")
    print("=" * 78)

    for i, result in enumerate(
        results[:20],
        1,
    ):

        print(
            f"{i:3d}. "
            f"runtime={result['median']:.9f}s "
            f"k={len(result['subset']):2d} "
            f"survivors={result['survivors']:<8d} "
            f"subset={result['subset']}"
        )

    # ------------------------------------------------------------------------
    # FINAL WINNER
    # ------------------------------------------------------------------------

    winner = results[0]

    print()
    print("=" * 78)
    print("FINAL WINNER")
    print("=" * 78)

    print()
    print("Winner:")
    print(winner["subset"])

    print(
        f"size      = {len(winner['subset'])}"
    )

    print(
        f"survivors = {winner['survivors']}"
    )

    print(
        f"median    = {winner['median']:.9f}s"
    )

    print(
        f"mean      = {winner['mean']:.9f}s"
    )

    print(
        f"min       = {winner['min']:.9f}s"
    )

    print(
        f"max       = {winner['max']:.9f}s"
    )

    p, q, x = winner["factor"]

    print(
        f"factor    = ({p}, {q})"
    )

    print(
        f"found x   = {x}"
    )

    print(
        f"true x    = {instance['true_x']}"
    )

    # ------------------------------------------------------------------------
    # REPORT COMPARISON TO EXP 207
    # ------------------------------------------------------------------------

    exp207_result = None

    for result in results:
        if result["subset"] == exp207:
            exp207_result = result
            break

    if exp207_result is not None:

        speed_ratio = (
            exp207_result["median"]
            / winner["median"]
        )

        print()
        print(
            f"vs Exp207 = {speed_ratio:.3f}x"
        )

        if speed_ratio > 1:
            print(
                "Winner is faster than Exp207."
            )
        elif speed_ratio < 1:
            print(
                "Exp207 is faster than the new winner."
            )
        else:
            print(
                "Same measured runtime as Exp207."
            )

    print()
    print("-" * 78)

    return winner


# ============================================================================
# ENTRY POINT
# ============================================================================

def main():

    print("=" * 78)
    print("START EXPERIMENT 215")
    print("PARETO-FRONTIER RUNTIME SEARCH")
    print("=" * 78)

    all_winners = []

    for instance in INSTANCES:

        winner = run_instance(
            instance
        )

        all_winners.append(
            (
                instance["name"],
                winner,
            )
        )

    # ------------------------------------------------------------------------
    # GLOBAL SUMMARY
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("EXP 215 SUMMARY")
    print("=" * 78)

    for name, winner in all_winners:

        print()
        print(name)

        print(
            f"  winner    = {winner['subset']}"
        )

        print(
            f"  size      = {len(winner['subset'])}"
        )

        print(
            f"  survivors = {winner['survivors']}"
        )

        print(
            f"  median    = {winner['median']:.9f}s"
        )

        p, q, x = winner["factor"]

        print(
            f"  factor    = ({p}, {q})"
        )

        print(
            f"  x         = {x}"
        )

    print()
    print("=" * 78)
    print("FINISHED EXPERIMENT 215")
    print("=" * 78)


if __name__ == "__main__":
    main()