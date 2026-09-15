#!/usr/bin/env python3

import itertools
import math
import random
import statistics
import time

import numpy as np


# ============================================================================
# START EXPERIMENT 216
# RANDOM SEMIPRIME PARETO SEARCH
# ============================================================================
#
# Purpose:
#
#   Generate fresh random semiprimes p*q for each target bit size and test
#   whether the QR-prime optimization discovered on fixed benchmark numbers
#   generalizes to new numbers.
#
# Important implementation detail:
#
#   NEVER compute x*x directly with NumPy int64 here.
#
#   Instead:
#
#       x^2 - n (mod p)
#
#   is evaluated as:
#
#       ((x mod p)^2 - (n mod p)) mod p
#
#   Since p <= 199, (x mod p)^2 is tiny and cannot overflow.
#
# ============================================================================


# ============================================================================
# CONFIGURATION
# ============================================================================

TARGET_BITS = (
    48,
    54,
    60,
    66,
    70,
    74,
    78,
)

PRIME_LIMIT = 199

EXTENSION_DEPTHS = (
    1,
    2,
    3,
)

BENCHMARK_REPEATS = 11

TOP_SURVIVOR_CONTROLS_PER_K = 5

#
# This limits the Fermat search interval.
#
# The randomly generated instances are therefore not arbitrary semiprimes:
# they are random semiprimes with a controlled factor distance.
#
MAX_INTERVAL = 8_000_000

#
# None = fresh random seed every run.
#
# Put an integer here for reproducibility.
#
USER_SEED = None


# ============================================================================
# BASE SUBSETS
# ============================================================================
#
# These are the best known starting points from the earlier experiments.
#
# They are NOT assumed to be optimal for the random numbers.
#
# Experiment 216 will test whether extensions around these subsets still
# produce useful runtime behavior on fresh n.
#

BASE_SUBSETS = {

    48: (
        7, 11, 17, 29, 31, 37,
        41, 47, 53, 59, 71, 73,
    ),

    54: (
        3, 5, 7, 17, 19, 23, 29,
        37, 41, 43, 47, 61, 67,
    ),

    60: (
        5, 7, 11, 13, 17, 19,
        23, 41, 47, 53, 59, 73,
    ),

    66: (
        3, 5, 11, 17, 19, 31,
        41, 43, 53, 59, 61,
    ),

    70: (
            3, 5, 11, 17, 19, 31,
            41, 43, 53, 59, 61,
        ),

    74: (
            3, 5, 11, 17, 19, 31,
            41, 43, 53, 59, 61,
        ),

    78: (
            3, 5, 11, 17, 19, 31,
            41, 43, 53, 59, 61,
        ),
}


# ============================================================================
# PRIME GENERATION
# ============================================================================

def sieve_primes(limit):
    sieve = np.ones(
        limit + 1,
        dtype=np.bool_,
    )

    sieve[:2] = False

    for p in range(
        2,
        math.isqrt(limit) + 1,
    ):
        if sieve[p]:
            sieve[
                p * p :
                limit + 1 :
                p
            ] = False

    return tuple(
        int(x)
        for x in np.flatnonzero(sieve)
    )


ALL_PRIMES = sieve_primes(
    PRIME_LIMIT
)

#
# Prime 2 is excluded.
#
# For odd n:
#
#     x^2 - n == 0 (mod 2)
#
# for every integer x.
#
# Therefore it contributes no useful filtering information.
#

PRIME_POOL = tuple(
    p
    for p in ALL_PRIMES
    if p != 2
)


# ============================================================================
# MILLER-RABIN
# ============================================================================
#
# Deterministic for 64-bit integers.
#

MR_BASES_64 = (
    2,
    325,
    9375,
    28178,
    450775,
    9780504,
    1795265022,
)


def is_prime(n):
    if n < 2:
        return False

    small_primes = (
        2,
        3,
        5,
        7,
        11,
        13,
        17,
        19,
        23,
        29,
        31,
        37,
    )

    for p in small_primes:

        if n == p:
            return True

        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while (d & 1) == 0:
        d >>= 1
        s += 1

    for a in MR_BASES_64:

        if a % n == 0:
            continue

        x = pow(
            a,
            d,
            n,
        )

        if (
            x == 1
            or x == n - 1
        ):
            continue

        for _ in range(
            s - 1
        ):

            x = (
                x * x
            ) % n

            if x == n - 1:
                break

        else:
            return False

    return True


def next_prime(n):
    if n <= 2:
        return 2

    if n % 2 == 0:
        n += 1

    while not is_prime(n):
        n += 2

    return n


# ============================================================================
# RANDOM SEMIPRIME GENERATION
# ============================================================================

def generate_random_semiprime(
    bits,
    rng,
):
    """
    Generate random odd primes p < q such that:

        n = p*q

    has exactly `bits` bits.

    Factor distance is controlled so Fermat remains practical.
    """

    half = bits // 2

    p_min = (
        1 << (half - 1)
    )

    p_max = (
        1 << half
    )

    attempts = 0

    while True:

        attempts += 1

        p_candidate = rng.randrange(
            p_min,
            p_max,
        )

        p = next_prime(
            p_candidate
        )

        if p >= p_max:
            continue

        #
        # Controlled factor distance.
        #
        if bits <= 48:

            min_delta = max(
                10_000,
                p // 4,
            )

            max_delta = max(
                min_delta + 1,
                p // 2,
            )

        elif bits <= 54:

            min_delta = max(
                20_000,
                p // 20,
            )

            max_delta = max(
                min_delta + 1,
                p // 8,
            )

        elif bits <= 60:

            min_delta = max(
                20_000,
                p // 200,
            )

            max_delta = max(
                min_delta + 1,
                p // 50,
            )

        else:

            min_delta = max(
                20_000,
                p // 1000,
            )

            max_delta = max(
                min_delta + 1,
                p // 200,
            )

        delta = rng.randint(
            min_delta,
            max_delta,
        )

        q_candidate = (
            p + delta
        )

        q = next_prime(
            q_candidate
        )

        if q <= p:
            continue

        n = p * q

        if n.bit_length() != bits:
            continue

        start_x = math.isqrt(n)

        if (
            start_x * start_x
            < n
        ):
            start_x += 1

        true_x = (
            p + q
        ) // 2

        interval = (
            true_x
            - start_x
            + 1
        )

        if interval <= 0:
            continue

        if interval > MAX_INTERVAL:
            continue

        true_y = (
            q - p
        ) // 2

        #
        # Exact Fermat identity.
        #

        assert (
            true_x * true_x
            - n
            == true_y * true_y
        )

        assert (
            (
                true_x
                - true_y
            )
            *
            (
                true_x
                + true_y
            )
            == n
        )

        return {
            "bits": bits,
            "p": p,
            "q": q,
            "n": n,
            "start_x": start_x,
            "true_x": true_x,
            "true_y": true_y,
            "interval": interval,
            "attempts": attempts,
        }


# ============================================================================
# QR MASK GENERATION
# ============================================================================

def build_bool_qr_mask(
    n,
    start_x,
    length,
    p,
):
    """
    Build QR mask for one prime.

    True iff:

        x^2 - n

    is a quadratic residue modulo p.

    IMPORTANT:
        x*x is NOT computed directly.
        We first reduce x modulo p.
    """

    #
    # QR lookup table.
    #

    qr = np.zeros(
        p,
        dtype=np.bool_,
    )

    residues = np.arange(
        p,
        dtype=np.int64,
    )

    squares = (
        residues * residues
    ) % p

    qr[squares] = True

    #
    # n modulo p is tiny.
    #

    n_mod_p = n % p

    mask = np.empty(
        length,
        dtype=np.bool_,
    )

    chunk_size = 1_000_000

    for lo in range(
        0,
        length,
        chunk_size,
    ):

        hi = min(
            lo + chunk_size,
            length,
        )

        #
        # x itself may be ~2^33.
        # That's still safe in int64.
        #
        # We immediately reduce modulo p.
        #

        xs = np.arange(
            start_x + lo,
            start_x + hi,
            dtype=np.int64,
        )

        residues_x = (
            xs % p
        )

        #
        # residues_x < 199.
        #
        # Therefore this multiplication is tiny.
        #

        values = (
            (
                residues_x
                * residues_x
            )
            - n_mod_p
        ) % p

        mask[
            lo:hi
        ] = qr[
            values
        ]

    return mask


def bool_to_python_int(
    mask,
):
    packed = np.packbits(
        mask,
        bitorder="little",
    )

    return int.from_bytes(
        packed.tobytes(),
        byteorder="little",
        signed=False,
    )


def build_masks(
    instance,
):
    n = instance["n"]

    start_x = (
        instance["start_x"]
    )

    length = (
        instance["interval"]
    )

    bool_masks = {}
    int_masks = {}

    print()
    print(
        "Building QR masks for primes <= "
        f"{PRIME_LIMIT}"
    )

    total_t0 = (
        time.perf_counter()
    )

    for p in PRIME_POOL:

        t0 = (
            time.perf_counter()
        )

        mask = (
            build_bool_qr_mask(
                n=n,
                start_x=start_x,
                length=length,
                p=p,
            )
        )

        bool_masks[p] = mask

        int_masks[p] = (
            bool_to_python_int(
                mask
            )
        )

        dt = (
            time.perf_counter()
            - t0
        )

        print(
            f"  prime={p:<3d} "
            f"build={dt:.6f}s"
        )

    total_dt = (
        time.perf_counter()
        - total_t0
    )

    print(
        f"mask build total = "
        f"{total_dt:.6f}s"
    )

    return {
        "length": length,
        "bool": bool_masks,
        "int": int_masks,
    }


# ============================================================================
# MASK INTERSECTION
# ============================================================================

def make_subset_mask(
    masks,
    subset,
    representation,
):
    if representation == "bool":

        result = np.ones(
            masks["length"],
            dtype=np.bool_,
        )

        for p in subset:

            result &= (
                masks["bool"][p]
            )

            if not result.any():
                break

        return result

    if representation == "int":

        result = (
            1
            << masks["length"]
        ) - 1

        for p in subset:

            result &= (
                masks["int"][p]
            )

            if result == 0:
                break

        return result

    raise ValueError(
        f"Unknown representation: "
        f"{representation}"
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

    if representation == "bool":

        return int(
            np.count_nonzero(
                mask
            )
        )

    return mask.bit_count()


# ============================================================================
# FACTOR RECOVERY
# ============================================================================

def recover_from_bool_mask(
    instance,
    mask,
):
    n = instance["n"]

    start_x = (
        instance["start_x"]
    )

    indices = (
        np.flatnonzero(
            mask
        )
    )

    for idx in indices:

        x = (
            start_x
            + int(idx)
        )

        y2 = (
            x * x
            - n
        )

        if y2 < 0:
            continue

        y = math.isqrt(
            y2
        )

        if (
            y * y
            != y2
        ):
            continue

        p = x - y
        q = x + y

        if (
            p > 1
            and q > 1
            and p * q == n
        ):

            return (
                p,
                q,
                x,
            )

    return (
        None,
        None,
        None,
    )


def recover_from_int_mask(
    instance,
    mask,
):
    n = instance["n"]

    start_x = (
        instance["start_x"]
    )

    while mask:

        lowbit = (
            mask & -mask
        )

        idx = (
            lowbit.bit_length()
            - 1
        )

        x = (
            start_x
            + idx
        )

        y2 = (
            x * x
            - n
        )

        if y2 >= 0:

            y = math.isqrt(
                y2
            )

            if (
                y * y
                == y2
            ):

                p = x - y
                q = x + y

                if (
                    p > 1
                    and q > 1
                    and p * q == n
                ):

                    return (
                        p,
                        q,
                        x,
                    )

        mask ^= lowbit

    return (
        None,
        None,
        None,
    )


def solve_with_subset(
    instance,
    masks,
    subset,
):
    representation = (
        instance[
            "representation"
        ]
    )

    mask = make_subset_mask(
        masks,
        subset,
        representation,
    )

    if representation == "bool":

        return recover_from_bool_mask(
            instance,
            mask,
        )

    return recover_from_int_mask(
        instance,
        mask,
    )


# ============================================================================
# PARETO FRONTIER
# ============================================================================

def pareto_frontier(
    candidates,
):
    """
    Objective:

        minimize k
        minimize survivors

    Candidate A dominates candidate B when:

        k_A <= k_B
        survivors_A <= survivors_B

    and at least one is strict.
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

        if (
            survivors
            < best_survivors
        ):

            frontier.append(
                (
                    subset,
                    survivors,
                )
            )

            best_survivors = (
                survivors
            )

    return frontier


# ============================================================================
# BENCHMARK
# ============================================================================

def benchmark_subset(
    instance,
    masks,
    subset,
):
    samples = []

    factor = None

    for _ in range(
        BENCHMARK_REPEATS
    ):

        t0 = (
            time.perf_counter()
        )

        factor = (
            solve_with_subset(
                instance,
                masks,
                subset,
            )
        )

        dt = (
            time.perf_counter()
            - t0
        )

        samples.append(dt)

    return {
        "median": statistics.median(
            samples
        ),
        "mean": statistics.mean(
            samples
        ),
        "min": min(samples),
        "max": max(samples),
        "factor": factor,
    }


# ============================================================================
# RUN RANDOM INSTANCE
# ============================================================================

def run_instance(
    instance,
):
    bits = instance["bits"]

    print()
    print(
        "#" * 78
    )
    print(
        f"# RANDOM {bits}-BIT INSTANCE"
    )
    print(
        "#" * 78
    )

    # ------------------------------------------------------------------------
    # Random instance information
    # ------------------------------------------------------------------------

    print()
    print(
        "Generated random factors:"
    )

    print(
        f"p         = "
        f"{instance['p']}"
    )

    print(
        f"q         = "
        f"{instance['q']}"
    )

    print(
        f"n         = "
        f"{instance['n']}"
    )

    print(
        f"n bits    = "
        f"{instance['n'].bit_length()}"
    )

    print(
        f"start_x   = "
        f"{instance['start_x']}"
    )

    print(
        f"true_x    = "
        f"{instance['true_x']}"
    )

    print(
        f"true_y    = "
        f"{instance['true_y']}"
    )

    print(
        f"interval  = "
        f"{instance['interval']}"
    )

    print(
        f"generation attempts = "
        f"{instance['attempts']}"
    )

    representation = (
        instance[
            "representation"
        ]
    )

    print()
    print(
        f"Representation = "
        f"{representation}"
    )

    # ------------------------------------------------------------------------
    # Build masks
    # ------------------------------------------------------------------------

    masks = build_masks(
        instance
    )

    # ------------------------------------------------------------------------
    # Base subset
    # ------------------------------------------------------------------------

    base = tuple(
        sorted(
            BASE_SUBSETS[bits]
        )
    )

    print()
    print(
        "=" * 78
    )
    print(
        "BASE INSTANCE"
    )
    print(
        "=" * 78
    )

    print()
    print(
        "Base subset:"
    )

    print(base)

    base_survivors = (
        count_survivors(
            masks,
            base,
            representation,
        )
    )

    print(
        f"Base size      = "
        f"{len(base)}"
    )

    print(
        f"Base survivors = "
        f"{base_survivors}"
    )

    # ------------------------------------------------------------------------
    # Available primes
    # ------------------------------------------------------------------------

    available = [
        p
        for p in PRIME_POOL
        if p not in base
    ]

    print()
    print(
        "Available additional primes:"
    )

    print(
        tuple(available)
    )

    # ------------------------------------------------------------------------
    # Exhaustive 1/2/3 extensions
    # ------------------------------------------------------------------------

    candidates = {
        base: base_survivors
    }

    for depth in EXTENSION_DEPTHS:

        print()
        print(
            f"SEARCHING {depth}-PRIME "
            f"EXTENSIONS"
        )

        evaluated = 0

        t0 = (
            time.perf_counter()
        )

        for added in itertools.combinations(
            available,
            depth,
        ):

            subset = tuple(
                sorted(
                    base
                    + tuple(added)
                )
            )

            survivors = (
                count_survivors(
                    masks,
                    subset,
                    representation,
                )
            )

            candidates[
                subset
            ] = survivors

            evaluated += 1

        dt = (
            time.perf_counter()
            - t0
        )

        print(
            f"  evaluated = "
            f"{evaluated}"
        )

        print(
            f"  time      = "
            f"{dt:.6f}s"
        )

    candidate_list = list(
        candidates.items()
    )

    # ------------------------------------------------------------------------
    # Pareto frontier
    # ------------------------------------------------------------------------

    frontier = (
        pareto_frontier(
            candidate_list
        )
    )

    print()
    print(
        "=" * 78
    )
    print(
        "PARETO FRONTIER"
    )
    print(
        "=" * 78
    )

    for i, (
        subset,
        survivors,
    ) in enumerate(
        frontier,
        1,
    ):

        print(
            f"{i:3d}. "
            f"k={len(subset):2d} "
            f"survivors="
            f"{survivors:<8d} "
            f"subset={subset}"
        )

    # ------------------------------------------------------------------------
    # Benchmark set
    # ------------------------------------------------------------------------
    #
    # Every Pareto point plus the lowest-survivor few candidates for each k.
    #

    benchmark_candidates = {
        subset: survivors
        for subset, survivors in frontier
    }

    by_k = {}

    for subset, survivors in (
        candidate_list
    ):

        k = len(subset)

        by_k.setdefault(
            k,
            [],
        ).append(
            (
                subset,
                survivors,
            )
        )

    for k, values in by_k.items():

        values.sort(
            key=lambda item: item[1]
        )

        for subset, survivors in (
            values[
                :TOP_SURVIVOR_CONTROLS_PER_K
            ]
        ):

            benchmark_candidates[
                subset
            ] = survivors

    benchmark_list = list(
        benchmark_candidates.items()
    )

    benchmark_list.sort(
        key=lambda item: (
            len(item[0]),
            item[1],
        )
    )

    # ------------------------------------------------------------------------
    # Benchmark
    # ------------------------------------------------------------------------

    print()
    print(
        "=" * 78
    )
    print(
        "END-TO-END BENCHMARK"
    )
    print(
        "=" * 78
    )

    results = []

    for i, (
        subset,
        survivors,
    ) in enumerate(
        benchmark_list,
        1,
    ):

        result = benchmark_subset(
            instance,
            masks,
            subset,
        )

        p_found, q_found, x_found = (
            result["factor"]
        )

        #
        # Every benchmark must factor successfully.
        #

        if p_found is None:

            raise RuntimeError(
                "Factor recovery failed:"
                f" bits={bits}"
                f" subset={subset}"
            )

        assert (
            p_found
            * q_found
            == instance["n"]
        )

        assert (
            x_found
            == instance["true_x"]
        )

        results.append(
            {
                "subset": subset,
                "survivors": survivors,
                "median": result[
                    "median"
                ],
                "mean": result[
                    "mean"
                ],
                "min": result[
                    "min"
                ],
                "max": result[
                    "max"
                ],
                "factor": result[
                    "factor"
                ],
            }
        )

        print(
            f"{i:3d}. "
            f"runtime="
            f"{result['median']:.9f}s "
            f"k={len(subset):2d} "
            f"survivors="
            f"{survivors:<8d} "
            f"subset={subset}"
        )

    # ------------------------------------------------------------------------
    # Runtime ranking
    # ------------------------------------------------------------------------

    results.sort(
        key=lambda item:
            item["median"]
    )

    print()
    print(
        "=" * 78
    )
    print(
        "ACTUAL RUNTIME RANKING"
    )
    print(
        "=" * 78
    )

    for i, result in enumerate(
        results[:20],
        1,
    ):

        print(
            f"{i:3d}. "
            f"runtime="
            f"{result['median']:.9f}s "
            f"k={len(result['subset']):2d} "
            f"survivors="
            f"{result['survivors']:<8d} "
            f"subset="
            f"{result['subset']}"
        )

    # ------------------------------------------------------------------------
    # Winner
    # ------------------------------------------------------------------------

    winner = results[0]

    print()
    print(
        "=" * 78
    )
    print(
        "FINAL WINNER"
    )
    print(
        "=" * 78
    )

    print()
    print(
        "Winner:"
    )

    print(
        winner["subset"]
    )

    print(
        f"size      = "
        f"{len(winner['subset'])}"
    )

    print(
        f"survivors = "
        f"{winner['survivors']}"
    )

    print(
        f"median    = "
        f"{winner['median']:.9f}s"
    )

    print(
        f"mean      = "
        f"{winner['mean']:.9f}s"
    )

    print(
        f"min       = "
        f"{winner['min']:.9f}s"
    )

    print(
        f"max       = "
        f"{winner['max']:.9f}s"
    )

    p_found, q_found, x_found = (
        winner["factor"]
    )

    print(
        f"factor    = "
        f"({p_found}, {q_found})"
    )

    print(
        f"found x   = "
        f"{x_found}"
    )

    print(
        f"true x    = "
        f"{instance['true_x']}"
    )

    # ------------------------------------------------------------------------
    # Base comparison
    # ------------------------------------------------------------------------

    base_result = None

    for result in results:

        if result["subset"] == base:

            base_result = result
            break

    if base_result is not None:

        ratio = (
            base_result["median"]
            / winner["median"]
        )

        print()
        print(
            f"vs base = "
            f"{ratio:.3f}x"
        )

        if ratio > 1.0:

            print(
                "Winner is faster than "
                "the base subset."
            )

        elif ratio < 1.0:

            print(
                "Base subset is faster than "
                "the winner."
            )

        else:

            print(
                "Winner matches base runtime."
            )

    print()
    print(
        "-" * 78
    )

    return winner


# ============================================================================
# MAIN
# ============================================================================

def main():

    print(
        "=" * 78
    )

    print(
        "START EXPERIMENT 216"
    )

    print(
        "RANDOM SEMIPRIME PARETO SEARCH"
    )

    print(
        "=" * 78
    )

    # ------------------------------------------------------------------------
    # RNG
    # ------------------------------------------------------------------------

    if USER_SEED is None:

        seed = (
            random.SystemRandom()
            .randrange(
                0,
                2**63,
            )
        )

    else:

        seed = USER_SEED

    rng = random.Random(
        seed
    )

    print()
    print(
        f"RANDOM SEED = {seed}"
    )

    print()
    print(
        "Prime pool:"
    )

    print(
        PRIME_POOL
    )

    print(
        f"Prime count = "
        f"{len(PRIME_POOL)}"
    )

    # ------------------------------------------------------------------------
    # Generate all random instances first.
    #
    # This makes the complete experiment deterministic once the seed is known.
    # ------------------------------------------------------------------------

    instances = []

    for bits in TARGET_BITS:

        print()
        print(
            f"Generating random "
            f"{bits}-bit semiprime..."
        )

        instance = (
            generate_random_semiprime(
                bits,
                rng,
            )
        )

        #
        # Same representation choices as previous experiments.
        #

        if bits <= 54:

            instance[
                "representation"
            ] = "bool"

        else:

            instance[
                "representation"
            ] = "int"

        instances.append(
            instance
        )

        print(
            f"  p        = "
            f"{instance['p']}"
        )

        print(
            f"  q        = "
            f"{instance['q']}"
        )

        print(
            f"  n        = "
            f"{instance['n']}"
        )

        print(
            f"  interval = "
            f"{instance['interval']}"
        )

    # ------------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------------

    winners = []

    for instance in instances:

        winner = run_instance(
            instance
        )

        winners.append(
            (
                instance,
                winner,
            )
        )

    # ------------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------------

    print()
    print(
        "=" * 78
    )

    print(
        "EXP 216 SUMMARY"
    )

    print(
        "=" * 78
    )

    for instance, winner in winners:

        bits = instance["bits"]

        print()
        print(
            f"{bits}-BIT RANDOM INSTANCE"
        )

        print(
            f"  p         = "
            f"{instance['p']}"
        )

        print(
            f"  q         = "
            f"{instance['q']}"
        )

        print(
            f"  n         = "
            f"{instance['n']}"
        )

        print(
            f"  interval  = "
            f"{instance['interval']}"
        )

        print(
            f"  winner    = "
            f"{winner['subset']}"
        )

        print(
            f"  size      = "
            f"{len(winner['subset'])}"
        )

        print(
            f"  survivors = "
            f"{winner['survivors']}"
        )

        print(
            f"  median    = "
            f"{winner['median']:.9f}s"
        )

        p_found, q_found, x_found = (
            winner["factor"]
        )

        print(
            f"  factor    = "
            f"({p_found}, {q_found})"
        )

        print(
            f"  x         = "
            f"{x_found}"
        )

    print()
    print(
        "=" * 78
    )

    print(
        "FINISHED EXPERIMENT 216"
    )

    print(
        "=" * 78
    )


if __name__ == "__main__":
    main()