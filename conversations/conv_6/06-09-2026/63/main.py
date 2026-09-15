#!/usr/bin/env python3

# ============================================================
# START EXPERIMENT 219
#
# WIDE RANDOM SEMIPRIMES
#
# p and q are generated independently.
#
# IMPORTANT:
# We do NOT reject an instance because its Fermat interval is
# too large.  Wide random factors naturally produce enormous
# Fermat intervals.
#
# QR/Fermat benchmarking is only performed when the interval
# is below MAX_BENCHMARK_INTERVAL.
#
# ============================================================

import math
import random
import time
from statistics import median

import numpy as np


# ------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------

TARGET_BITS = (
    48,
    54,
    60,
    66,
    70,
    74,
    78,
    80,
    84,
    88,
)

# Require q to be at least this much larger than p.
#
# 0.20 means:
#
#       q >= 1.20 * p
#
MIN_GAP_RATIO = 0.20

# Only benchmark the actual x interval when it is reasonably
# small enough to materialize.
#
# IMPORTANT:
# This is NOT a generation restriction anymore.
MAX_BENCHMARK_INTERVAL = 50_000_000

# Number of benchmark repetitions.
BENCHMARK_REPEATS = 3

# Random seed.
RANDOM_SEED = 150314232428026902


# ------------------------------------------------------------
# QR PRIME POOL
# ------------------------------------------------------------

# 2 is excluded because for odd n:
#
#     x^2 - n == 0 (mod 2)
#
# so it provides no useful QR filtering.

PRIME_POOL = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113,
    127, 131, 137, 139, 149, 151, 157, 163,
    167, 173, 179, 181, 191, 193, 197, 199
]


# ------------------------------------------------------------
# BASE SUBSETS
# ------------------------------------------------------------

BASE_SUBSETS = {
    48: (
        7, 11, 17, 29, 31, 37,
        41, 47, 53, 59, 71, 73
    ),

    54: (
        3, 5, 7, 17, 19, 23,
        29, 37, 41, 43, 47, 61, 67
    ),

    60: (
        5, 7, 11, 13, 17, 19,
        23, 41, 47, 53, 59, 73
    ),

    66: (
        3, 5, 11, 17, 19, 31,
        41, 43, 53, 59, 61
    ),

    70: (
        3, 5, 11, 17, 19, 31,
        41, 43, 53, 59, 61
    ),

    74: (
        3, 5, 11, 17, 19, 31,
        41, 43, 53, 59, 61
    ),

    78: (
        3, 5, 11, 17, 19, 31,
        41, 43, 53, 59, 61
    ),

    80: (
        3, 5, 11, 17, 19, 31,
        41, 43, 53, 59, 61
    ),

    84: (
        3, 5, 11, 19, 31, 41,
        43, 53, 59, 61
    ),

    88: (
        3, 5, 11, 17, 19, 31,
        41, 43, 53, 59, 61
    ),
}


# ------------------------------------------------------------
# MILLER-RABIN
# ------------------------------------------------------------

def is_probable_prime(n: int, rounds: int = 20) -> bool:
    if n < 2:
        return False

    small_primes = (
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37, 41, 43
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

    for _ in range(rounds):
        a = random.randrange(2, n - 1)

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):
            x = (x * x) % n

            if x == n - 1:
                break

        else:
            return False

    return True


def random_prime(bits: int) -> int:
    while True:
        x = random.getrandbits(bits)

        # Force exact bit length.
        x |= (1 << (bits - 1))

        # Force odd.
        x |= 1

        if is_probable_prime(x):
            return x


# ------------------------------------------------------------
# RANDOM WIDE SEMIPRIME
# ------------------------------------------------------------

def generate_wide_semiprime(bits: int):
    """
    Generate two independent primes.

    The only constraint is:

        q/p >= 1 + MIN_GAP_RATIO

    There is deliberately NO Fermat-interval rejection here.
    """

    attempts = 0

    while True:
        attempts += 1

        p = random_prime(bits)
        q = random_prime(bits)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        gap = q - p
        gap_ratio = gap / p

        if gap_ratio < MIN_GAP_RATIO:
            continue

        n = p * q

        start_x = math.isqrt(n)

        if start_x * start_x < n:
            start_x += 1

        true_x = (p + q) // 2

        true_y = (q - p) // 2

        interval = true_x - start_x

        # Useful normalized measure.
        #
        # This tells us how enormous the Fermat search is
        # relative to sqrt(n).

        sqrt_n = math.isqrt(n)

        normalized_interval = interval / sqrt_n

        return {
            "bits": bits,
            "p": p,
            "q": q,
            "n": n,
            "gap": gap,
            "gap_ratio": gap_ratio,
            "start_x": start_x,
            "true_x": true_x,
            "true_y": true_y,
            "interval": interval,
            "normalized_interval": normalized_interval,
            "attempts": attempts,
        }


# ------------------------------------------------------------
# QR MASK
# ------------------------------------------------------------

def build_bool_qr_mask(xs, n, prime):
    """
    Return a boolean QR mask.

    Overflow is avoided by reducing x modulo prime before
    performing the square.
    """

    p = prime

    residues = np.zeros(p, dtype=np.bool_)

    for r in range(p):
        residues[(r * r) % p] = True

    x_mod_p = xs % p

    values = (
        x_mod_p * x_mod_p
        - (n % p)
    ) % p

    return residues[values]


# ------------------------------------------------------------
# QR BENCHMARK
# ------------------------------------------------------------

def benchmark_subset(instance, subset):

    interval = instance["interval"]

    if interval > MAX_BENCHMARK_INTERVAL:
        return {
            "skipped": True,
            "reason": (
                f"interval {interval:,} > "
                f"MAX_BENCHMARK_INTERVAL "
                f"{MAX_BENCHMARK_INTERVAL:,}"
            )
        }

    start_x = instance["start_x"]
    true_x = instance["true_x"]
    n = instance["n"]

    xs = np.arange(
        start_x,
        true_x + 1,
        dtype=np.int64
    )

    survivors = None
    timings = []

    for _ in range(BENCHMARK_REPEATS):

        t0 = time.perf_counter()

        mask = np.ones(
            len(xs),
            dtype=np.bool_
        )

        for prime in subset:
            mask &= build_bool_qr_mask(
                xs,
                n,
                prime
            )

        elapsed = time.perf_counter() - t0

        survivors = int(mask.sum())

        timings.append(elapsed)

    return {
        "skipped": False,
        "survivors": survivors,
        "median_time": median(timings),
        "min_time": min(timings),
        "max_time": max(timings),
    }


# ------------------------------------------------------------
# FERMAТ BASELINE
# ------------------------------------------------------------

def fermat_factor(
    n: int,
    start_x: int,
    max_iters=None
):
    """
    Ordinary Fermat search.

    Returns:
        (p, q, iterations)

    or None if max_iters is exceeded.
    """

    x = start_x

    # For odd n the valid x is odd.
    if (x & 1) == 0:
        x += 1

    iterations = 0

    while True:

        y2 = x * x - n

        if y2 >= 0:

            y = math.isqrt(y2)

            if y * y == y2:

                p = x - y
                q = x + y

                if (
                    p > 1
                    and q > 1
                    and p * q == n
                ):
                    return p, q, iterations

        x += 2
        iterations += 1

        if (
            max_iters is not None
            and iterations >= max_iters
        ):
            return None


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------

def main():

    if RANDOM_SEED is not None:
        random.seed(RANDOM_SEED)

    print("=" * 72)
    print("EXPERIMENT 219")
    print("WIDE RANDOM SEMIPRIME GENERATOR")
    print("=" * 72)

    print(
        f"random_seed                = "
        f"{RANDOM_SEED}"
    )

    print(
        f"min_gap_ratio              = "
        f"{MIN_GAP_RATIO}"
    )

    print(
        f"minimum q/p ratio          = "
        f"{1.0 + MIN_GAP_RATIO:.4f}"
    )

    print(
        f"max_benchmark_interval     = "
        f"{MAX_BENCHMARK_INTERVAL:,}"
    )

    print(
        f"target_bits                = "
        f"{TARGET_BITS}"
    )

    print()

    for bits in TARGET_BITS:

        print("-" * 72)
        print(f"[{bits}-BIT INSTANCE]")
        print("-" * 72)

        instance = generate_wide_semiprime(bits)

        p = instance["p"]
        q = instance["q"]
        n = instance["n"]

        print(f"p                         = {p}")
        print(f"q                         = {q}")
        print(f"n                         = {n}")

        print(
            f"q-p                       = "
            f"{instance['gap']:,}"
        )

        print(
            f"(q-p)/p                   = "
            f"{instance['gap_ratio']:.8f}"
        )

        print(
            f"q/p                       = "
            f"{q / p:.8f}"
        )

        print(
            f"start_x                   = "
            f"{instance['start_x']}"
        )

        print(
            f"true_x                    = "
            f"{instance['true_x']}"
        )

        print(
            f"true_y                    = "
            f"{instance['true_y']}"
        )

        print(
            f"Fermat interval           = "
            f"{instance['interval']:,}"
        )

        print(
            f"interval / sqrt(n)        = "
            f"{instance['normalized_interval']:.8f}"
        )

        print(
            f"generation attempts       = "
            f"{instance['attempts']}"
        )

        subset = BASE_SUBSETS.get(bits)

        print()

        if subset is None:
            print("No base QR subset configured.")
            print()
            continue

        print(
            f"base QR subset             = "
            f"{subset}"
        )

        print(
            f"subset size                = "
            f"{len(subset)}"
        )

        print()

        # ----------------------------------------------------
        # QR TEST
        # ----------------------------------------------------

        print("QR benchmark:")

        qr_result = benchmark_subset(
            instance,
            subset
        )

        if qr_result["skipped"]:

            print(
                f"  SKIPPED: "
                f"{qr_result['reason']}"
            )

        else:

            print(
                f"  survivors                = "
                f"{qr_result['survivors']}"
            )

            print(
                f"  median                   = "
                f"{qr_result['median_time']:.9f} s"
            )

            print(
                f"  minimum                  = "
                f"{qr_result['min_time']:.9f} s"
            )

            print(
                f"  maximum                  = "
                f"{qr_result['max_time']:.9f} s"
            )

        print()

        # ----------------------------------------------------
        # FERMAТ
        # ----------------------------------------------------
        #
        # Do NOT actually attempt huge Fermat intervals.
        #

        if instance["interval"] <= MAX_BENCHMARK_INTERVAL:

            print("Ordinary Fermat:")

            t0 = time.perf_counter()

            result = fermat_factor(
                n,
                instance["start_x"]
            )

            elapsed = time.perf_counter() - t0

            if result is None:

                print("  recovery                 = FAILED")

            else:

                fp, fq, iterations = result

                print(
                    f"  p                        = "
                    f"{fp}"
                )

                print(
                    f"  q                        = "
                    f"{fq}"
                )

                print(
                    f"  iterations               = "
                    f"{iterations:,}"
                )

                print(
                    f"  time                     = "
                    f"{elapsed:.9f} s"
                )

        else:

            print(
                "Ordinary Fermat            = SKIPPED "
                f"(interval too large)"
            )

        print()

    print("=" * 72)
    print("FINISHED EXPERIMENT 219")
    print("=" * 72)


if __name__ == "__main__":
    main()